!*****
!***** COMPRESSIBLE AERODYNAMICS & AEROACOUSTICS RESEARCH CODE (CANARD)
!*****

module mo_canard_driver
   use mo_kind,       ONLY : ieee64, ieee32, nr, ni, int64
   use mo_parameters, ONLY : zero, one, half, n45no, two, pi, gamm1, gam, quarter
   use mo_mpi,        ONLY : p_get_n_processes, p_get_process_ID, p_barrier,       &
                           & p_sum, p_max, p_get_global_comm
   use mo_io,         ONLY : read_input_main, allocate_io_memory,                  &
                           & output_init, vminmax, read_restart_file,              &
                           & write_restart_file,                                   &
                           & read_grid_parallel, write_output_file
   use mo_io_server,  ONLY : io_server_stop, io_server_init,      &
                           & t_io_server_interface, io_server_write_output
   use mo_domdcomp,   ONLY : t_domdcomp
   use mo_grid,       ONLY : t_grid
   use mo_gridgen,    ONLY : t_grid_geom, read_input_gridgen, makegrid,            &
                           & get_grid_geometry
   use mo_sponge,     ONLY : spongeup, spongego, read_input_sponge
   use mo_gcbc,       ONLY : gcbc_init, gcbc_setup, gcbc_comm, gcbc_update,        &
                           & extracon, wall_condition_update, average_surface
   use mo_numerics,   ONLY : t_numerics
   use mo_physics,    ONLY : t_physics
   use mo_timer,      ONLY : timer_init, timer_start, timer_stop, timer_print
   use mo_timer,      ONLY : timer_loop, timer_filter, timer_timestep,             &
                           & timer_VSstress, timer_fluxes, timer_GCBC,             &
                           & timer_averaging, timer_recording, timer_tot_output,   &
                           & timer_output, timer_total
   implicit none
   PUBLIC
   contains

   subroutine canard_driver(laio, lmodel_role, mbk, ndata)
      logical, intent(in)          :: laio
      logical, intent(in)          :: lmodel_role
      integer(kind=ni), intent(in) :: mbk
      integer(kind=ni), intent(in) :: ndata
   
   integer(kind=ni),parameter    :: nkrk = 4

   integer(kind=ni)    :: m, nn, ll, nout, lis, lie, l, ndati, nnn
   real(kind=nr)       :: res, ra0, ra1, fctr, dtko, dtk, dtsum
   integer(kind=int64) :: nlmx
   type(t_domdcomp)    :: p_domdcomp
   type(t_numerics)    :: p_numerics
   type(t_grid)        :: p_grid
   type(t_grid_geom)   :: p_grid_geom
   type(t_physics)     :: p_physics
   type(t_io_server_interface) :: p_io_server_interface
   integer(kind=ni)    :: nts
   integer(kind=ni)    :: nrestart
   real(kind=nr)       :: cfl
   real(kind=nr)       :: dts, dte
   real(kind=nr)       :: tmax
   real(kind=nr)       :: timo
   integer(kind=ni)    :: ndt
   integer(kind=ni)    :: nk
   integer(kind=ni)    :: lim
   integer(kind=ni)    :: n
   integer(kind=ni)    :: comm_glob
   real(kind=nr)       :: dt
   integer(kind=ni)    :: j, k, kp, jp
   integer(kind=ni)    :: nrecs, myid, mpro, nvar
   logical             :: ltimer, loutput = .true.
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
   integer(kind=ni), dimension(:), allocatable   :: lpos_temp

!===== PREPARATION FOR PARALLEL COMPUTING

   myid = p_get_process_ID()
   mpro = p_get_n_processes() - 1
   comm_glob = p_get_global_comm()

   inquire(iolength=ll) real(1.0,kind=ieee32); nrecs=ll

!===== INPUT PARAMETERS
   if (myid==0) write(*,*) "Canard: read input parameters"

   call read_input_main(nts, nrestart, cfl, &
                        tmax, ltimer)
   if (ndata < 0) loutput = .false.

   call p_numerics%read()

   call p_physics%read()

   call read_input_gridgen
   
   call read_input_sponge

   call p_domdcomp%allocate(mbk,mpro)
   call p_domdcomp%read(mbk,lmodel_role)
   call get_grid_geometry(p_grid_geom)

!===== DOMAIN DECOMPOSITION INITIALIZATION
   if (myid==0) write(*,*) "Canard: initialize domain decomposition"

   call p_domdcomp%init(mbk)
   lim=(p_domdcomp%lxi+1)+(p_domdcomp%let+1)+(p_domdcomp%lze+1)-1

!===== TIMERS INITIALIZATION
   if (ltimer) call timer_init()

!===== WRITING START POSITIONS IN OUTPUT FILE
   if (myid==0) write(*,*) "Canard: initialize IO"
   allocate(lpos_temp(0:mpro))
   call allocate_io_memory(mbk, ndata)
   call output_init(p_domdcomp, mbk, ndata, lpos_temp=lpos_temp)

!===== IO SERVER INITIALIZATION
   if (laio .and. myid==0) write(*,*) "Canard: model intiialize io server"
   if (laio) call io_server_init(mbk, p_domdcomp, p_io_server_interface, mpro, lpos_temp, lmodel_role)
   deallocate(lpos_temp)

!===== ALLOCATION OF MAIN ARRAYS
   call p_physics%allocate(p_domdcomp%lmx)
   call p_numerics%allocate(lim, p_domdcomp%nbsize, p_domdcomp%lmx)

   ! main program local arrays
   allocate(times(0:ndata))
   times(:) = zero
   allocate(qo(0:p_domdcomp%lmx,5))
   allocate(qb(0:p_domdcomp%lmx,5))
   allocate(qa(0:p_domdcomp%lmx,5))
   allocate(de(0:p_domdcomp%lmx,5))
   allocate(ss(0:p_domdcomp%lmx,3))
   allocate(rr(0:p_domdcomp%lmx,3))
   allocate(varr(0:p_domdcomp%lmx))
   allocate(p(0:p_domdcomp%lmx))

!===== EXTRA COEFFICIENTS FOR DOMAIN BOUNDARIES INITIALIZATION
   if (myid==0) write(*,*) "Canard: initialize numerics"
   call p_numerics%init_extra

!===== PENTADIAGONAL MATRICES INITIALIZATION
   call p_numerics%init(p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, p_domdcomp%nbc, lim)

!===== GRID GENERATION & CALCULATION OF GRID METRICS
   if (myid==0) write(*,*) "Canard: grid generation"

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
   if (myid==0) write(*,*) "Canard: initialize GCBC"
   call gcbc_init(p_domdcomp, p_grid%yaco)

!===== POINT JUNCTION SEARCH
   call p_domdcomp%search_point(mbk)

!===== LINE JUNCTION SEARCH
   call p_domdcomp%search_line(mbk)

!===== SETTING UP OUTPUT FILE & STORING GRID DATA
   ndati=-1
   if (loutput) then
      nlmx=3*(p_domdcomp%lmx+1)-1      
      allocate(vart(0:nlmx))
      do nn=1,3
         ! ss contains the grid data at this points from previous subroutines call
         vart((nn-1)*(p_domdcomp%lmx+1):nn*(p_domdcomp%lmx+1)-1) = real(ss(0:p_domdcomp%lmx,nn), kind=ieee32)
      end do
      nvar = 3
      if (laio) then
         call io_server_write_output(vart, nvar, ndati, times(0), p_domdcomp, p_io_server_interface)
      else
         do nn=1,3
            call vminmax(p_domdcomp, vart((nn-1)*(p_domdcomp%lmx+1):nn*(p_domdcomp%lmx+1)-1), nn)
         end do
         call write_output_file(p_domdcomp, mbk, ndata, times, nlmx, vart, ndati, nvar)
      end if
      deallocate(vart)
   end if

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

   dtsum=zero
   
   do while(timo<tmax.and.(dt/=zero.or.n<=2))

      if(myid==0) then
         write(*,"(' n =',i8,'   time =',f12.5)") n,timo
      end if

!----- FILTERING & RE-INITIALISING
      if (ltimer) call timer_start(timer_filter)
      do m=1,5
         call p_numerics%mpigo(qa(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                             p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 3*(m-1)+1, &
                             p_domdcomp%lxi, p_domdcomp%let, m)
         call p_numerics%filte(qa(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                             p_domdcomp%lze, p_domdcomp%ijk, 1, m)
         call p_numerics%mpigo(qa(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                             p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 3*(m-1)+2, &
                             p_domdcomp%lxi, p_domdcomp%let, m)
         call p_numerics%filte(qa(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                             p_domdcomp%lze, p_domdcomp%ijk, 2, m)
         call p_numerics%mpigo(qa(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                             p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 3*(m-1)+3, &
                             p_domdcomp%lxi, p_domdcomp%let, m)
         call p_numerics%filte(qa(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                             p_domdcomp%lze, p_domdcomp%ijk, 3, m)
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
               call p_physics%calc_time_step(p_domdcomp%lmx, p_grid, de, ss(:,1), cfl, dte)
            end if
            dt=dts+(dte-dts)*sin(0.05_nr*pi*(n-ndt))**two

            nout=0
            res=(ndati+1)*(tmax)/ndata
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

         if ( nk == nkrk .and. ( timo + quarter - tmax )**two < ( half * dt )**two ) then
            call extracon(p_domdcomp, p_grid, p_physics, qa, p)
         end if

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
      if (loutput) then
         if(timo>-(tmax)/ndata) then
            if (ltimer) call timer_start(timer_recording)
            dtsum=dtsum+dt
            fctr=half*dt
            qb(:,:)=qb(:,:)+fctr*(qo(:,:)+qa(:,:))
            if(nout==1) then
               if (ltimer) call timer_start(timer_tot_output)
               times(ndati)=timo-half*dtsum
               if(n==1) then
                  qb(:,:)=qo(:,:)
               else
                  ra0=one/dtsum
                  qb(:,:)=ra0*qb(:,:)
               end if
               rr(:,1)=one/qb(:,1)               
               nlmx = 5*(p_domdcomp%lmx+1)-1
               allocate(vart(0:nlmx))
               do nn=1,5
                  select case(nn)
                  case(1)
                     vart((nn-1)*(p_domdcomp%lmx+1):nn*(p_domdcomp%lmx+1)-1) = &
                                real(qb(:,nn), kind=ieee32)
                  case(2:4)
                     vart((nn-1)*(p_domdcomp%lmx+1):nn*(p_domdcomp%lmx+1)-1) = &
                                real(rr(:,1)*qb(:,nn)+p_physics%umf(nn-1), kind=ieee32)
                  case(5)
                     vart((nn-1)*(p_domdcomp%lmx+1):nn*(p_domdcomp%lmx+1)-1) = &
                                real(gamm1*(qb(:,nn)-half*rr(:,1)*(qb(:,2)*qb(:,2)+qb(:,3)*qb(:,3)+qb(:,4)*qb(:,4))), kind=ieee32)
                  end select                  
               end do
               if (ltimer) call timer_start(timer_output)
               nvar = 5
               if (laio) then
                  call io_server_write_output(vart, nvar, ndati, times(ndati), p_domdcomp, p_io_server_interface)
               else
                  do nn=1,5
                     nnn=3+5*ndati+nn
                     call vminmax(p_domdcomp, vart((nn-1)*(p_domdcomp%lmx+1):nn*(p_domdcomp%lmx+1)-1), nnn)
                  end do
                  call write_output_file(p_domdcomp, mbk, ndata, times, nlmx, vart, ndati, nvar)
               end if
               if (ltimer) call timer_stop(timer_output)
               deallocate(vart)
               
               ! reset variables
               dtsum=zero
               qb(:,:)=zero
               if (ltimer) call timer_stop(timer_tot_output)
            end if
            if (ltimer) call timer_stop(timer_recording)
         end if
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

!===== FINAL CLEAN UP

   deallocate(qo,qa,qb,de,rr,ss,p)
   call p_grid%deallocate
   call p_physics%deallocate
   if (dt==zero .and. myid==0) write(*,*) "Canard: Overflow."

   if (laio) call io_server_stop

   if (ltimer) call timer_stop(timer_total)
   
!===== TIMERS PRINT
   if (ltimer) call timer_print()

!===== END OF JOB
   
   call p_barrier(comm_glob)
   if(myid==0) write(*,*) "Canard: Finished."

   end subroutine canard_driver

 end module mo_canard_driver

!*****
