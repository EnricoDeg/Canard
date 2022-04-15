!*****
!***** COMPRESSIBLE AERODYNAMICS & AEROACOUSTICS RESEARCH CODE (CANARD)
!*****

program main3d

   use mo_mpi, ONLY : mpro, npro, myid, p_start, p_stop, &
                      p_barrier, p_sum, p_max
   use mo_io
   use mo_domdcomp
   use mo_grid
   use mo_sponge
   use mo_gcbc
   use mo_numerics
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
   amachoo=sqrt(amach1*amach1+amach2*amach2+amach3*amach3)
   if(amachoo>sml) then
      reoo=reoo/amachoo
   end if
   srefoo=111/tempoo
   srefp1dre=(srefoo+one)/reoo
   sqrtrema=sqrt(reoo*amachoo)
   sqrtremai=one/max(sqrtrema,sml)
   uoo(:)=(/amach1,amach2,amach3/)

   call allocate_io_memory

   call allocate_domdcomp(mbk)

   call read_inputp

!===== DOMAIN DECOMPOSITION / BOUNDARY INFORMATION / SUBDOMAIN SIZES

   call domdcomp_init(mbk, nthick, nbody, nviscous)

!===== WRITING START POSITIONS IN OUTPUT FILE

   call output_init

!===== ALLOCATION OF MAIN ARRAYS

   call allocate_memory
   call allocate_numerics(lim)

!===== EXTRA COEFFICIENTS FOR DOMAIN BOUNDARIES

   call init_extracoeff_bounds

!===== PENTADIAGONAL MATRICES FOR DIFFERENCING & FILETERING

   call init_penta(lxi, let, lze, nbc)

!===== GRID INPUT & CALCULATION OF GRID METRICS

   call calc_grid
   call calc_grid_metrics

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
      varr(:)=ss(:,nn)
      write(0,rec=nn) varr(:)
      call vminmax(nn)
   end do
   close(0)

!===== SETTING UP SPONGE ZONE PARAMETERS

   call spongeup

!===== INITIAL CONDITIONS

   if(nts==0) then
      n=0
      ndt=0
      dt=zero
      dts=zero
      dte=zero
      timo=zero
      call initialo
   else
      call read_restart_file
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
         rr(:,1)=qa(:,m)
         call mpigo(1,nrone,n45no,3*(m-1)+1)
         call filte(nnf(1),1)
         call mpigo(1,nrone,n45no,3*(m-1)+2)
         call filte(nnf(2),1)
         call mpigo(1,nrone,n45no,3*(m-1)+3)
         call filte(nnf(3),1)
         qa(:,m)=rr(:,1)
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
                  if(nviscous==1) then
                     res=maxval(de(:,1)*ss(:,1)*rr(:,1)*ss(:,2)*ss(:,2))
                     call p_max(res, fctr)
                     ra1=half/fctr
                  end if
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

!----- VISCOUS SHEAR STRESSES & HEAT FLUXES

         if(nviscous==1) then
            de(:,1)=ss(:,1)

            rr(:,1)=de(:,2)
            m=2
            call mpigo(0,nrone,n45no,m)
            call deriv(3,1,m)
            call deriv(2,1,m)
            call deriv(1,1,m)
            txx(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
            hzz(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
            tzx(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

            rr(:,1)=de(:,3)
            m=3
            call mpigo(0,nrone,n45no,m)
            call deriv(3,1,m)
            call deriv(2,1,m)
            call deriv(1,1,m)
            txy(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
            tyy(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
            hxx(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

            rr(:,1)=de(:,4)
            m=4
            call mpigo(0,nrone,n45no,m)
            call deriv(3,1,m)
            call deriv(2,1,m)
            call deriv(1,1,m)
            hyy(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
            tyz(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
            tzz(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

            rr(:,1)=de(:,5)
            m=5
            call mpigo(0,nrone,n45no,m)
            call deriv(3,1,m)
            call deriv(2,1,m)
            call deriv(1,1,m)
            ss(:,1)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
            ss(:,2)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
            ss(:,3)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

            rr(:,1)=de(:,1)*yaco(:)
            rr(:,2)=gamm1prndtli*rr(:,1)
            de(:,5)=twothirds*(txx(:)+tyy(:)+tzz(:))

            txx(:)=rr(:,1)*(txx(:)+txx(:)-de(:,5))
            tyy(:)=rr(:,1)*(tyy(:)+tyy(:)-de(:,5))
            tzz(:)=rr(:,1)*(tzz(:)+tzz(:)-de(:,5))
            txy(:)=rr(:,1)*(txy(:)+hzz(:))
            tyz(:)=rr(:,1)*(tyz(:)+hxx(:))
            tzx(:)=rr(:,1)*(tzx(:)+hyy(:))
            hxx(:)=rr(:,2)*ss(:,1)+de(:,2)*txx(:)+de(:,3)*txy(:)+de(:,4)*tzx(:)
            hyy(:)=rr(:,2)*ss(:,2)+de(:,2)*txy(:)+de(:,3)*tyy(:)+de(:,4)*tyz(:)
            hzz(:)=rr(:,2)*ss(:,3)+de(:,2)*tzx(:)+de(:,3)*tyz(:)+de(:,4)*tzz(:)
         end if

!----- CALCULATION OF FLUX DERIVATIVES

         rr(:,1)=de(:,2)+umf(1)
         rr(:,2)=de(:,3)+umf(2)
         rr(:,3)=de(:,4)+umf(3)
         ss(:,1)=xim(:,1)*rr(:,1)+xim(:,2)*rr(:,2)+xim(:,3)*rr(:,3)
         ss(:,2)=etm(:,1)*rr(:,1)+etm(:,2)*rr(:,2)+etm(:,3)*rr(:,3)
         ss(:,3)=zem(:,1)*rr(:,1)+zem(:,2)*rr(:,2)+zem(:,3)*rr(:,3)

         rr(:,1)=qa(:,1)*ss(:,1)
         rr(:,2)=qa(:,1)*ss(:,2)
         rr(:,3)=qa(:,1)*ss(:,3)
         m=1
         call mpigo(0,nrall,n45no,m)
         call deriv(1,1,m)
         call deriv(2,2,m)
         call deriv(3,3,m)
         de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

         rr(:,1)=qa(:,2)*ss(:,1)+xim(:,1)*p(:)
         rr(:,2)=qa(:,2)*ss(:,2)+etm(:,1)*p(:)
         rr(:,3)=qa(:,2)*ss(:,3)+zem(:,1)*p(:)
         if(nviscous==1) then
            rr(:,1)=rr(:,1)-xim(:,1)*txx(:)-xim(:,2)*txy(:)-xim(:,3)*tzx(:)
            rr(:,2)=rr(:,2)-etm(:,1)*txx(:)-etm(:,2)*txy(:)-etm(:,3)*tzx(:)
            rr(:,3)=rr(:,3)-zem(:,1)*txx(:)-zem(:,2)*txy(:)-zem(:,3)*tzx(:)
         end if
         m=2
         call mpigo(0,nrall,n45no,m)
         call deriv(1,1,m)
         call deriv(2,2,m)
         call deriv(3,3,m)
         de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

         rr(:,1)=qa(:,3)*ss(:,1)+xim(:,2)*p(:)
         rr(:,2)=qa(:,3)*ss(:,2)+etm(:,2)*p(:)
         rr(:,3)=qa(:,3)*ss(:,3)+zem(:,2)*p(:)
         if(nviscous==1) then
            rr(:,1)=rr(:,1)-xim(:,1)*txy(:)-xim(:,2)*tyy(:)-xim(:,3)*tyz(:)
            rr(:,2)=rr(:,2)-etm(:,1)*txy(:)-etm(:,2)*tyy(:)-etm(:,3)*tyz(:)
            rr(:,3)=rr(:,3)-zem(:,1)*txy(:)-zem(:,2)*tyy(:)-zem(:,3)*tyz(:)
         end if
         m=3
         call mpigo(0,nrall,n45no,m)
         call deriv(1,1,m)
         call deriv(2,2,m)
         call deriv(3,3,m)
         de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

         rr(:,1)=qa(:,4)*ss(:,1)+xim(:,3)*p(:)
         rr(:,2)=qa(:,4)*ss(:,2)+etm(:,3)*p(:)
         rr(:,3)=qa(:,4)*ss(:,3)+zem(:,3)*p(:)
         if(nviscous==1) then
            rr(:,1)=rr(:,1)-xim(:,1)*tzx(:)-xim(:,2)*tyz(:)-xim(:,3)*tzz(:)
            rr(:,2)=rr(:,2)-etm(:,1)*tzx(:)-etm(:,2)*tyz(:)-etm(:,3)*tzz(:)
            rr(:,3)=rr(:,3)-zem(:,1)*tzx(:)-zem(:,2)*tyz(:)-zem(:,3)*tzz(:)
         end if
         m=4
         call mpigo(0,nrall,n45no,m)
         call deriv(1,1,m)
         call deriv(2,2,m)
         call deriv(3,3,m)
         de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

         de(:,5)=qa(:,5)+p(:)
         rr(:,1)=de(:,5)*ss(:,1)-p(:)*(umf(1)*xim(:,1)+umf(2)*xim(:,2)+umf(3)*xim(:,3))
         rr(:,2)=de(:,5)*ss(:,2)-p(:)*(umf(1)*etm(:,1)+umf(2)*etm(:,2)+umf(3)*etm(:,3))
         rr(:,3)=de(:,5)*ss(:,3)-p(:)*(umf(1)*zem(:,1)+umf(2)*zem(:,2)+umf(3)*zem(:,3))
         if(nviscous==1) then
            rr(:,1)=rr(:,1)-xim(:,1)*hxx(:)-xim(:,2)*hyy(:)-xim(:,3)*hzz(:)
            rr(:,2)=rr(:,2)-etm(:,1)*hxx(:)-etm(:,2)*hyy(:)-etm(:,3)*hzz(:)
            rr(:,3)=rr(:,3)-zem(:,1)*hxx(:)-zem(:,2)*hyy(:)-zem(:,3)*hzz(:)
         end if
         m=5
         call mpigo(0,nrall,n45no,m)
         call deriv(1,1,m)
         call deriv(2,2,m)
         call deriv(3,3,m)
         de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

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
      if(nviscous==1) then
         deallocate(txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz)
      end if
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

 end program main3d

!*****
