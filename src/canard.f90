!*****
!***** COMPRESSIBLE AERODYNAMICS & AEROACOUSTICS RESEARCH CODE (CANARD)
!*****

program canard
   use mo_kind,       ONLY : ieee64, ieee32, nr, ni, int64
   use mo_parameters, ONLY : zero, one, half, n45no, two, pi, gamm1, gam
   use mo_io,         ONLY : cdata
   use mo_mpi,        ONLY : p_get_n_processes, p_get_process_ID, p_start, p_stop, &
                           & p_barrier, p_sum, p_max
   use mo_io,         ONLY : read_input_main, allocate_io_memory,                  &
                           & output_init, vminmax, read_restart_file,              &
                           & write_restart_file, write_output_file,                &
                           & read_grid_parallel
   use mo_domdcomp,   ONLY : t_domdcomp
   use mo_grid,       ONLY : t_grid
   use mo_gridgen,    ONLY : t_grid_geom, read_input_gridgen, makegrid,            &
                           & get_grid_geometry
   use mo_sponge,     ONLY : spongeup, spongego, read_input_sponge
   use mo_gcbc,       ONLY : gcbc_init, gcbc_setup, gcbc_comm, gcbc_update,        &
                           & extracon, wall_condition_update, average_surface,     &
                           & read_input_gcbc
   use mo_numerics,   ONLY : t_numerics
   use mo_physics,    ONLY : t_physics
   use mo_timer,      ONLY : timer_init, timer_start, timer_stop, timer_print
   use mo_timer,      ONLY : timer_loop, timer_filter, timer_timestep,             &
                           & timer_VSstress, timer_fluxes, timer_GCBC,             &
                           & timer_averaging, timer_recording, timer_tot_output,   &
                           & timer_output, timer_total
   implicit none

   integer(kind=ni)    :: m, nn, ll, nout, lis, lie, l, ndati
   real(kind=nr)       :: res, ra0, ra1, fctr, dtko, dtk, dtsum
   integer(kind=int64) :: nlmx
   type(t_domdcomp)    :: p_domdcomp
   type(t_numerics)    :: p_numerics
   type(t_grid)        :: p_grid
   type(t_grid_geom)   :: p_grid_geom
   type(t_physics)     :: p_physics
   integer(kind=ni)    :: mbk
   integer(kind=ni)    :: nts, nscrn, ndata, ndatafl, ndataav
   integer(kind=ni)    :: nrestart
   real(kind=nr)       :: cfl, dto
   integer(kind=ni)    :: nbody
   real(kind=nr)       :: dts, dte
   real(kind=nr)       :: tsam
   real(kind=nr)       :: tmax
   integer(kind=ni)    :: nkrk
   real(kind=nr)       :: timo
   integer(kind=ni)    :: ndt
   integer(kind=ni)    :: nk
   integer(kind=ni)    :: lim
   integer(kind=ni)    :: n
   real(kind=nr)       :: dt
   integer(kind=ni)    :: j, k, kp, jp
   integer(kind=ni)    :: nrecs, myid, mpro
   logical             :: ltimer
   real(kind=nr), dimension(:), allocatable     :: times
   real(kind=nr), dimension(:), allocatable     :: p
   real(kind=nr), dimension(:,:), allocatable   :: qo
   real(kind=nr), dimension(:,:), allocatable   :: qb
   real(kind=nr), dimension(:,:), allocatable   :: qa
   real(kind=nr), dimension(:,:), allocatable   :: de
   real(kind=nr), dimension(:,:), allocatable   :: ss
   real(kind=nr), dimension(:,:), allocatable   :: rr
   real(kind=ieee32), dimension(:), allocatable :: vmean
   real(kind=ieee32), dimension(:), allocatable :: vart
   real(kind=ieee32), dimension(:), allocatable :: varr
   integer(kind=ni), dimension(:,:),   allocatable :: lio

!===== PREPARATION FOR PARALLEL COMPUTING

   CALL p_start

   myid = p_get_process_ID()
   mpro = p_get_n_processes() - 1

   inquire(iolength=ll) real(1.0,kind=ieee32); nrecs=ll

!===== INPUT PARAMETERS

   call read_input_main(mbk, nts, nscrn, ndata, ndatafl, ndataav, nrestart, cfl, &
                    dto, tsam, tmax, nkrk, nbody, ltimer)

   call p_numerics%read()

   call read_input_gcbc

   call p_physics%read()

   call read_input_gridgen
   
   call read_input_sponge

   call p_domdcomp%allocate(mbk,mpro)
   call p_domdcomp%read()
   call get_grid_geometry(p_grid_geom)

!===== DOMAIN DECOMPOSITION INITIALIZATION

   call p_domdcomp%init(mbk, p_grid_geom%nthick, nbody)
   lim=(p_domdcomp%lxi+1)+(p_domdcomp%let+1)+(p_domdcomp%lze+1)-1

!===== TIMERS INITIALIZATION
   if (ltimer) call timer_init()

!===== WRITING START POSITIONS IN OUTPUT FILE

   call allocate_io_memory(mbk, ndata)
   call output_init(p_domdcomp, mbk, ndata)

!===== ALLOCATION OF MAIN ARRAYS

   call p_physics%allocate(p_domdcomp%lmx)
   call p_numerics%allocate(lim, p_domdcomp%nbsize)

   ! main program local arrays
   allocate(times(0:ndata))
   allocate(qo(0:p_domdcomp%lmx,5))
   allocate(qb(0:p_domdcomp%lmx,5))
   allocate(qa(0:p_domdcomp%lmx,5))
   allocate(de(0:p_domdcomp%lmx,5))
   allocate(ss(0:p_domdcomp%lmx,3))
   allocate(rr(0:p_domdcomp%lmx,3))
   allocate(varr(0:p_domdcomp%lmx))
   allocate(p(0:p_domdcomp%lmx))

!===== EXTRA COEFFICIENTS FOR DOMAIN BOUNDARIES INITIALIZATION

   call p_numerics%init_extra

!===== PENTADIAGONAL MATRICES INITIALIZATION

   call p_numerics%init(p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, p_domdcomp%nbc, lim)

!===== GRID GENERATION & CALCULATION OF GRID METRICS
    
   allocate(lio(0:p_domdcomp%let,0:p_domdcomp%lze))
   call p_grid%allocate(p_domdcomp)

   do k=0,p_domdcomp%lze
      kp = k * ( p_domdcomp%leto + 1 ) * ( p_domdcomp%lxio + 1 )
      do j=0,p_domdcomp%let
         jp = j * ( p_domdcomp%lxio + 1 )
         lio(j,k) = jp + kp
      end do
   end do

   call makegrid(p_domdcomp%mb, p_domdcomp%lxio, p_domdcomp%leto, p_domdcomp%mo, mbk)
   call p_barrier
   call read_grid_parallel(p_domdcomp, ss, lio)
   call p_grid%grid_metrics(p_domdcomp, p_numerics, ss)

!===== EXTRA COEFFICIENTS FOR GCBC/GCIC INITIALIZATION

   call gcbc_init(p_domdcomp, p_grid%yaco)

!===== POINT JUNCTION SEARCH

   call p_domdcomp%search_point(mbk)

!===== LINE JUNCTION SEARCH

   call p_domdcomp%search_line(mbk)

!===== SETTING UP OUTPUT FILE & STORING GRID DATA

   open(0,file=cdata,status='unknown')
   close(0,status='delete') ! 'replace' not suitable as 'recl' may vary
   open(0,file=cdata,access='direct',form='unformatted',recl=nrecs*(p_domdcomp%lmx+1),status='new')
   do nn=1,3
      varr(:)=ss(:,nn) ! ss contains the grid data at this points from previous subroutines call
      write(0,rec=nn) varr(:)
      call vminmax(p_domdcomp, varr, nn)
   end do
   close(0)

!===== SETTING UP SPONGE ZONE PARAMETERS

   call spongeup(p_grid_geom, p_domdcomp%lmx, p_grid%yaco, de, ss) ! use ss which contains grid data

!===== INITIAL CONDITIONS

   if(nts==0) then
      n=0
      ndt=0
      dt=zero
      dts=zero
      dte=zero
      timo=zero
      call p_physics%init(p_domdcomp%lmx, qa, ss) ! use ss which contains grid data
   else
      call read_restart_file(p_domdcomp, qa, lio, dts, dte, timo, ndt, n, dt) ! ss is not used
   end if
   qb(:,:)=zero

!============================================
!===== BEGINNING OF TIME MARCHING IN SOLUTION
!============================================

   call p_barrier
   if (ltimer) call timer_start(timer_total)
   if (ltimer) call timer_start(timer_loop)

   ndati=-1
   dtsum=zero
   
   do while(timo<tmax.and.(dt/=zero.or.n<=2))

      if(myid**2+mod(n,nscrn)**2==0) then
         write(*,"(' n =',i8,'   time =',f12.5)") n,timo
      end if

!----- FILTERING & RE-INITIALISING
      if (ltimer) call timer_start(timer_filter)
      do m=1,5
         call p_numerics%mpigo(qa(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                             p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 3*(m-1)+1, &
                             p_domdcomp%lxi, p_domdcomp%let)
         call p_numerics%filte(qa(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                             p_domdcomp%lze, p_domdcomp%ijk, 1)
         call p_numerics%mpigo(qa(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                             p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 3*(m-1)+2, &
                             p_domdcomp%lxi, p_domdcomp%let)
         call p_numerics%filte(qa(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                             p_domdcomp%lze, p_domdcomp%ijk, 2)
         call p_numerics%mpigo(qa(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                             p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 3*(m-1)+3, &
                             p_domdcomp%lxi, p_domdcomp%let)
         call p_numerics%filte(qa(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                             p_domdcomp%lze, p_domdcomp%ijk, 3)
      end do
      if (ltimer) call timer_stop(timer_filter)
      qo(:,:)=qa(:,:)

!-------------------------------------
!----- BEGINNING OF RUNGE-KUTTA STAGES
!-------------------------------------
      do nk=1,nkrk

!----- MOVING FRAME VELOCITY & ACCELERATION BEFORE TIME ADVANCING

         dtko = dt * min( max( nk-2, 0 ),1 ) / ( nkrk - nk + 3 )
         dtk  = dt * min( nk-1, 1 ) / ( nkrk - nk + 2 )
         call p_physics%movef(dtko, dtk, timo)

!----- TEMPORARY STORAGE OF PRIMITIVE VARIABLES & PRESSURE

         de(:,1) = one / qa(:,1)
         de(:,2) = qa(:,2) * de(:,1)
         de(:,3) = qa(:,3) * de(:,1)
         de(:,4) = qa(:,4) * de(:,1)

         p(:)    = gamm1 * ( qa(:,5) - half * &
                           ( qa(:,2) * de(:,2) + qa(:,3) * de(:,3) + qa(:,4) * de(:,4) ) )
         de(:,5) = gam*p(:) * de(:,1)
         ss(:,1) = p_physics%srefp1dre * de(:,5)**1.5_nr / ( de(:,5) + p_physics%srefoo )

!----- DETERMINATION OF TIME STEP SIZE & OUTPUT TIME
         if(nk==1) then
            if (ltimer) call timer_start(timer_timestep)
            if(mod(n,10)==1) then
               ndt=n
               dts=dte
               if(dto<zero) then
                  call p_physics%calc_time_step(p_domdcomp%lmx, p_grid, de, ss(:,1), cfl, dte)
               else
                  dte = dto
               end if
            end if
            dt=dts+(dte-dts)*sin(0.05_nr*pi*(n-ndt))**two

            nout=0
            res=tsam+(ndati+1)*(tmax-tsam)/ndata
            if((timo-res)*(timo+dt-res)<=zero) then
               nout=1
               ndati=ndati+1
            end if
            if (ltimer) call timer_stop(timer_timestep)
         end if

         if (ltimer) call timer_start(timer_VSstress)
         call p_physics%calc_viscous_shear_stress(p_domdcomp, p_numerics, p_grid, de, ss(:,1))
         if (ltimer) call timer_stop(timer_VSstress)

         if (ltimer) call timer_start(timer_fluxes)
         call p_physics%calc_fluxes(p_domdcomp, p_numerics, p_grid, qa, p, de)
         if (ltimer) call timer_stop(timer_fluxes)

         if (ltimer) call timer_start(timer_GCBC)
!----- PREPARATION FOR GCBC & GCIC

         call gcbc_setup(p_domdcomp, p_numerics, p_grid, qa, p, de, p_physics%umf)

!----- INTERNODE COMMNICATION FOR GCIC

         call gcbc_comm(p_domdcomp, p_numerics)

!----- IMPLEMENTATION OF GCBC & GCIC

         call gcbc_update(p_domdcomp, p_numerics, p_grid, qa, p, de, nkrk, dt, p_physics%umf, p_physics%dudtmf)

         if (ltimer) call timer_stop(timer_GCBC)

!----- IMPLEMENTATION OF SPONGE CONDITION

         call spongego(p_domdcomp%lmx, qa, de)

!----- UPDATING CONSERVATIVE VARIABLES

         dtko=dt*min(nk-1,1)/(nkrk-nk+2)
         dtk=dt/(nkrk-nk+1)
         call p_physics%movef(dtko, dtk, timo)

         rr(:,1)=dtk*p_grid%yaco(:)
         qa(:,1)=qo(:,1)-rr(:,1)*de(:,1)
         qa(:,2)=qo(:,2)-rr(:,1)*de(:,2)
         qa(:,3)=qo(:,3)-rr(:,1)*de(:,3)
         qa(:,4)=qo(:,4)-rr(:,1)*de(:,4)
         qa(:,5)=qo(:,5)-rr(:,1)*de(:,5)

         call extracon(p_domdcomp, p_grid, p_physics, varr, qa, p, tmax, nkrk, timo, nk, dt)

!----- WALL TEMPERATURE/VELOCITY CONDITION
 
         call wall_condition_update(p_domdcomp, qa, p_physics%umf)

!----- POINT JUNCTION AVERAGING
         if (ltimer) call timer_start(timer_averaging)

         call p_domdcomp%average_point(qa)

!----- LINE JUNCTION AVERAGING

         call p_domdcomp%average_line(qa)

!----- INTERFACE SURFACE AVERAGING

         call average_surface(p_domdcomp, p_numerics, qa)

         if (ltimer) call timer_stop(timer_averaging)

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
         if (ltimer) call timer_start(timer_recording)
         dtsum=dtsum+dt
         fctr=half*dt
         qb(:,:)=qb(:,:)+fctr*(qo(:,:)+qa(:,:))
         if(nout==1) then
            times(ndati)=timo-half*dtsum
            open(0,file=cdata,access='direct',form='unformatted',recl=nrecs*(p_domdcomp%lmx+1),status='old')
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
                  varr(:)=rr(:,1)*qb(:,m)+p_physics%umf(m-1)
               case(5)
                  varr(:)=gamm1*(qb(:,m)-half*rr(:,1)*(qb(:,2)*qb(:,2)+qb(:,3)*qb(:,3)+qb(:,4)*qb(:,4)))
               end select
               nn=3+5*ndati+m
               write(0,rec=nn) varr(:)
               call vminmax(p_domdcomp, varr, nn)
            end do
            close(0)
            dtsum=zero
            qb(:,:)=zero
         end if
         if (ltimer) call timer_stop(timer_recording)
      end if

!==========================
!===== END OF TIME MARCHING
!==========================

   end do
   if (ltimer) call timer_stop(timer_loop)

!===== GENERATING RESTART DATA FILE

   if(nrestart==1) then
      call write_restart_file(p_domdcomp, qa, lio, dts, dte, timo, ndt, n, dt)
   end if

!===== POST-PROCESSING & GENERATING TECPLOT DATA FILE

   if(dt==zero) then
      if(myid==0) then
         write(*,*) "Overflow."
      end if
   else
      deallocate(qo,qa,qb,de,rr,ss,p)
      call p_grid%deallocate
      call p_physics%deallocate

      if(tmax>=tsam) then
         if (ltimer) call timer_start(timer_tot_output)
         nlmx=(3+5*(ndata+1))*(p_domdcomp%lmx+1)-1
         ll=5*(p_domdcomp%lmx+1)-1
         allocate(vart(0:nlmx),vmean(0:ll))
         open(9,file=cdata,access='direct',form='unformatted',recl=nrecs*(nlmx+1),status='old')
         read(9,rec=1) vart(:)
         close(9,status='delete')

!----- CALCULATING UNSTEADY FLUCTUATIONS

         if(ndatafl==1) then
            fctr=half/(times(ndata)-times(0))
            vmean(:)=zero
            do n=0,ndata
               lis=(3+5*n)*(p_domdcomp%lmx+1)
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
               lis=(3+5*n)*(p_domdcomp%lmx+1)
               lie=lis+ll
               vart(lis:lie)=vart(lis:lie)-vmean(:)
               do m=1,5
                  nn=3+5*n+m
                  l=lis+(m-1)*(p_domdcomp%lmx+1)
                  varr(:)=vart(l:l+p_domdcomp%lmx)
                  call vminmax(p_domdcomp, varr, nn)
               end do
            end do
         end if

!----- COLLECTING DATA FROM SUBDOMAINS & BUILDING TECPLOT OUTPUT FILES
         if (ltimer) call timer_start(timer_output)
         call write_output_file(p_domdcomp, mbk, ndata, times, nlmx, vart)
         if (ltimer) call timer_stop(timer_output)
!-----
         if (ltimer) call timer_stop(timer_tot_output)
      end if
   end if
   if (ltimer) call timer_stop(timer_total)
   
!===== TIMERS PRINT
   if (ltimer) call timer_print()

!===== END OF JOB

   if(myid==0) then
      write(*,*) "Finished."
   end if

   call p_stop

 end program canard

!*****