!*****
!***** COMPRESSIBLE AERODYNAMICS & AEROACOUSTICS RESEARCH CODE (CANARD)
!*****

program canard
   use mo_kind,       ONLY : ieee64, ieee32, nr
   use mo_parameters, ONLY : zero, one, half, n45no, two, pi, gamm1, gam
   use mo_vars,       ONLY : times, m, tmax, timo, res, ss, ra0, ra1, nts,         &
                           & nsigi, tsam, nscrn, nrestart, nrecs, nrecd, nout,     &
                           & nn, nkrk, nlmx, nk, ndt, ndati, ndatafl, ndataav,     &
                           & ndata, nbsize, nbody, nbc, n, mbk, lxi, let, lze,     &
                           & lmx, ll, lis, lie, l, lim, fctr, dto, dts, dtko,      &
                           & dtk, dte, dtsum, dt, cinput, cfl, cdata, varr,        &
                           & vart, vmean, txx, tyy, tzz, txy, tyz, tzx, hxx,       &
                           & hyy, hzz, qo, qa, qb, de, mcd, ijk, xim, etm, zem,    &
                           & rr, umf, nnf, p, yaco, srefoo, srefp1dre, &
                           & lxim, letm, lzem, lpos
   use mo_vars,       ONLY : allocate_memory
   use mo_mpi,        ONLY : mpro, npro, myid, p_start, p_stop, p_barrier, p_sum,  &
                           & p_max
   use mo_io,         ONLY : read_inputo, allocate_io_memory, read_inputp,         &
                           & output_init, vminmax, read_restart_file,              &
                           & write_restart_file, write_output_file
   use mo_domdcomp,   ONLY : allocate_domdcomp, domdcomp_init, search_point,       &
                           & search_line, average_line, average_point
   use mo_grid,       ONLY : calc_grid, calc_grid_metrics
   use mo_gridgen,    ONLY : nthick
   use mo_sponge,     ONLY : spongeup, spongego
   use mo_gcbc,       ONLY : gcbc_init, gcbc_setup, gcbc_comm, gcbc_update,        &
                           & extracon, wall_condition_update, average_surface
   use mo_numerics,   ONLY : allocate_numerics, init_extracoeff_bounds,            &
                           & init_penta, mpigo, filte, read_input_numerics
   use mo_physics,    ONLY : init_physics, initialo, movef,                        &
                           & calc_viscous_shear_stress, calc_fluxes
   implicit none

!===== PREPARATION FOR PARALLEL COMPUTING

   CALL p_start

   allocate(lxim(0:mpro),letm(0:mpro),lzem(0:mpro),lpos(0:mpro))

   inquire(iolength=ll) real(1.0,kind=ieee32); nrecs=ll
   inquire(iolength=ll) real(1.0,kind=ieee64); nrecd=ll

!===== INPUT PARAMETERS

   call read_inputo
   call read_input_numerics
    
   cinput=cinput

   call init_physics

   call allocate_io_memory

   call allocate_domdcomp(mbk)

   call read_inputp

!===== DOMAIN DECOMPOSITION / BOUNDARY INFORMATION / SUBDOMAIN SIZES

   call domdcomp_init(mbk, nthick, nbody)

!===== WRITING START POSITIONS IN OUTPUT FILE

   call output_init

!===== ALLOCATION OF MAIN ARRAYS

   call allocate_memory
   call allocate_numerics(lim, nbsize)

!===== EXTRA COEFFICIENTS FOR DOMAIN BOUNDARIES

   call init_extracoeff_bounds

!===== PENTADIAGONAL MATRICES FOR DIFFERENCING & FILETERING

   call init_penta(lxi, let, lze, nbc)

!===== GRID INPUT & CALCULATION OF GRID METRICS

   call calc_grid(ss)
   call calc_grid_metrics(ss)

!===== EXTRA COEFFICIENTS FOR GCBC/GCIC

   call gcbc_init

!===== POINT JUNCTION SEARCH

   call search_point

!===== LINE JUNCTION SEARCH

   call search_line

!===== SETTING UP OUTPUT FILE & STORING GRID DATA

   open(0,file=cdata,status='unknown')
   close(0,status='delete') ! 'replace' not suitable as 'recl' may vary
   open(0,file=cdata,access='direct',form='unformatted',recl=nrecs*(lmx+1),status='new')
   do nn=1,3
      varr(:)=ss(:,nn) ! ss contains the grid data at this points from previous subroutines call
      write(0,rec=nn) varr(:)
      call vminmax(nn)
   end do
   close(0)

!===== SETTING UP SPONGE ZONE PARAMETERS

   call spongeup ! use ss which contains grid data

!===== INITIAL CONDITIONS

   if(nts==0) then
      n=0
      ndt=0
      dt=zero
      dts=zero
      dte=zero
      timo=zero
      call initialo ! use ss which contains grid data
   else
      call read_restart_file ! ss is not used
   end if
   qb(:,:)=zero

!============================================
!===== BEGINNING OF TIME MARCHING IN SOLUTION
!============================================

   if(myid==0) then
      open(1,file='signal.dat',access='direct',form='formatted',recl=16,status='replace')
      close(1)
   end if
   call p_barrier

   ndati=-1
   nsigi=-1
   dtsum=zero
   
   do while(timo<tmax.and.(dt/=zero.or.n<=2))

      if(myid**2+mod(n,nscrn)**2==0) then
         write(*,"(' n =',i8,'   time =',f12.5)") n,timo
      end if

!----- FILTERING & RE-INITIALISING

      do m=1,5
         call mpigo(qa(:,m), ijk, nbc, mcd, nbsize,1,n45no,3*(m-1)+1)
         call filte(qa(:,m), lxi, let, lze, ijk, nnf(1))
         call mpigo(qa(:,m), ijk, nbc, mcd, nbsize,1,n45no,3*(m-1)+2)
         call filte(qa(:,m), lxi, let, lze, ijk, nnf(2))
         call mpigo(qa(:,m), ijk, nbc, mcd, nbsize,1,n45no,3*(m-1)+3)
         call filte(qa(:,m), lxi, let, lze, ijk, nnf(3))
      end do
      qo(:,:)=qa(:,:)

!-------------------------------------
!----- BEGINNING OF RUNGE-KUTTA STAGES
!-------------------------------------
      do nk=1,nkrk

!----- MOVING FRAME VELOCITY & ACCELERATION BEFORE TIME ADVANCING

         dtko=dt*min(max(nk-2,0),1)/(nkrk-nk+3)
         dtk=dt*min(nk-1,1)/(nkrk-nk+2)
         call movef(dtko,dtk)

!----- TEMPORARY STORAGE OF PRIMITIVE VARIABLES & PRESSURE

         de(:,1)=one/qa(:,1)
         de(:,2)=qa(:,2)*de(:,1)
         de(:,3)=qa(:,3)*de(:,1)
         de(:,4)=qa(:,4)*de(:,1)

         p(:)=gamm1*(qa(:,5)-half*(qa(:,2)*de(:,2)+qa(:,3)*de(:,3)+qa(:,4)*de(:,4)))
         de(:,5)=gam*p(:)*de(:,1)
         ss(:,1)=srefp1dre*de(:,5)**1.5_nr/(de(:,5)+srefoo)

!----- DETERMINATION OF TIME STEP SIZE & OUTPUT TIME

         if(nk==1) then
            if(mod(n,10)==1) then
               ndt=n
               dts=dte
               if(dto<zero) then
                  rr(:,1)=xim(:,1)*xim(:,1)+xim(:,2)*xim(:,2)+xim(:,3)*xim(:,3)&
                         +etm(:,1)*etm(:,1)+etm(:,2)*etm(:,2)+etm(:,3)*etm(:,3)&
                         +zem(:,1)*zem(:,1)+zem(:,2)*zem(:,2)+zem(:,3)*zem(:,3)
                  rr(:,2)=abs(xim(:,1)*(de(:,2)+umf(1))+xim(:,2)*(de(:,3)+umf(2))+xim(:,3)*(de(:,4)+umf(3)))&
                         +abs(etm(:,1)*(de(:,2)+umf(1))+etm(:,2)*(de(:,3)+umf(2))+etm(:,3)*(de(:,4)+umf(3)))&
                         +abs(zem(:,1)*(de(:,2)+umf(1))+zem(:,2)*(de(:,3)+umf(2))+zem(:,3)*(de(:,4)+umf(3)))
                  ss(:,2)=abs(yaco(:))
                  res=maxval((sqrt(de(:,5)*rr(:,1))+rr(:,2))*ss(:,2))
                  call p_max(res, fctr)
                  ra0=cfl/fctr
                  ra1=ra0
#ifdef VISCOUS
                     res=maxval(de(:,1)*ss(:,1)*rr(:,1)*ss(:,2)*ss(:,2))
                     call p_max(res, fctr)
                     ra1=half/fctr
#endif
                  dte=min(ra0,ra1)
               else
                  dte=dto
               end if
            end if
            dt=dts+(dte-dts)*sin(0.05_nr*pi*(n-ndt))**two

            nout=0
            res=tsam+(ndati+1)*(tmax-tsam)/ndata
            if((timo-res)*(timo+dt-res)<=zero) then
               nout=1
               ndati=ndati+1
            end if
         end if

         call calc_viscous_shear_stress

         call calc_fluxes
!----- PREPARATION FOR GCBC & GCIC

         call gcbc_setup

!----- INTERNODE COMMNICATION FOR GCIC

         call gcbc_comm

!----- IMPLEMENTATION OF GCBC & GCIC

         call gcbc_update

!----- IMPLEMENTATION OF SPONGE CONDITION

         call spongego

!----- UPDATING CONSERVATIVE VARIABLES

         dtko=dt*min(nk-1,1)/(nkrk-nk+2)
         dtk=dt/(nkrk-nk+1)
         call movef(dtko,dtk)

         rr(:,1)=dtk*yaco(:)
         qa(:,1)=qo(:,1)-rr(:,1)*de(:,1)
         qa(:,2)=qo(:,2)-rr(:,1)*de(:,2)
         qa(:,3)=qo(:,3)-rr(:,1)*de(:,3)
         qa(:,4)=qo(:,4)-rr(:,1)*de(:,4)
         qa(:,5)=qo(:,5)-rr(:,1)*de(:,5)

         call extracon

!----- WALL TEMPERATURE/VELOCITY CONDITION
 
         call wall_condition_update

!----- POINT JUNCTION AVERAGING

         call average_point

!----- LINE JUNCTION AVERAGING

         call average_line

!----- INTERFACE SURFACE AVERAGING

         call average_surface

!-------------------------------
!----- END OF RUNGE-KUTTA STAGES
!-------------------------------

      end do

!---------------------
!----- ADVANCE IN TIME
!---------------------

      n=n+1
      timo=timo+dt

!----- RECORDING INTERMEDIATE RESULTS

      if(timo>tsam-(tmax-tsam)/ndata) then
         dtsum=dtsum+dt
         fctr=half*dt
         qb(:,:)=qb(:,:)+fctr*(qo(:,:)+qa(:,:))
         if(nout==1) then
            times(ndati)=timo-half*dtsum
            open(0,file=cdata,access='direct',form='unformatted',recl=nrecs*(lmx+1),status='old')
            if(n==1) then
               qb(:,:)=qo(:,:)
            else
               ra0=ndataav/dtsum
               ra1=one-ndataav
               qb(:,:)=ra0*qb(:,:)+ra1*qa(:,:)
            end if
            rr(:,1)=one/qb(:,1)
            do m=1,5
               select case(m)
               case(1)
                  varr(:)=qb(:,m)
               case(2:4)
                  varr(:)=rr(:,1)*qb(:,m)+umf(m-1)
               case(5)
                  varr(:)=gamm1*(qb(:,m)-half*rr(:,1)*(qb(:,2)*qb(:,2)+qb(:,3)*qb(:,3)+qb(:,4)*qb(:,4)))
               end select
               nn=3+5*ndati+m
               write(0,rec=nn) varr(:)
               call vminmax(nn)
            end do
            close(0)
            dtsum=zero
            qb(:,:)=zero
         end if
      end if

!==========================
!===== END OF TIME MARCHING
!==========================

   end do

   if(myid==0) then
      open(9,file='timeouts.dat',status='replace')
      write(9,'(es15.7)') times(:)
      close(9)
   end if

!===== GENERATING RESTART DATA FILE

   if(nrestart==1) then
      call write_restart_file
   end if

!===== POST-PROCESSING & GENERATING TECPLOT DATA FILE

   if(dt==zero) then
      if(myid==0) then
         write(*,*) "Overflow."
      end if
   else
      deallocate(qo,qa,qb,de,xim,etm,zem,rr,ss,p,yaco)
#ifdef VISCOUS
         deallocate(txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz)
#endif
      if(tmax>=tsam) then
         nlmx=(3+5*(ndata+1))*(lmx+1)-1
         ll=5*(lmx+1)-1
         allocate(vart(0:nlmx),vmean(0:ll))
         open(9,file=cdata,access='direct',form='unformatted',recl=nrecs*(nlmx+1),status='old')
         read(9,rec=1) vart(:)
         close(9,status='delete')

!----- CALCULATING UNSTEADY FLUCTUATIONS

         if(ndatafl==1) then
            fctr=half/(times(ndata)-times(0))
            vmean(:)=zero
            do n=0,ndata
               lis=(3+5*n)*(lmx+1)
               lie=lis+ll
               nn=n/ndata
               if(n*(n-ndata)==0) then
                  ra0=fctr*(times(n+1-nn)-times(n-nn))
               else
                  ra0=fctr*(times(n+1)-times(n-1))
               end if
               vmean(:)=vmean(:)+ra0*vart(lis:lie)
            end do
            do n=0,ndata
               lis=(3+5*n)*(lmx+1)
               lie=lis+ll
               vart(lis:lie)=vart(lis:lie)-vmean(:)
               do m=1,5
                  nn=3+5*n+m
                  l=lis+(m-1)*(lmx+1)
                  varr(:)=vart(l:l+lmx)
                  call vminmax(nn)
               end do
            end do
         end if

!----- COLLECTING DATA FROM SUBDOMAINS & BUILDING TECPLOT OUTPUT FILES

         call write_output_file

!-----
      end if
   end if

!===== END OF JOB

   if(myid==0) then
      write(*,*) "Finished."
   end if

   call p_stop

 end program canard

!*****
