program dcanard
   use mo_kind,     ONLY : ni, nr, ieee32
   use mo_mpi,      ONLY : p_start, p_stop, p_get_process_ID, p_get_n_processes, &
                         & p_get_global_comm
   use mo_io,       ONLY : read_input_driver, allocate_io_memory, output_init,   &
                         & read_restart_file
   use mo_domdcomp, ONLY : t_domdcomp
   implicit none
   
   integer(kind=ni)    :: nio
   integer(kind=ni)    :: mbk
   integer(kind=ni)    :: ndata
   logical             :: lmodel_role
   logical             :: laio
   integer(kind=ni)    :: ll, ix
   real(kind=nr)       :: dts, dte
   real(kind=nr)       :: timo
   integer(kind=ni)    :: ndt
   integer(kind=ni)    :: n
   real(kind=nr)       :: dt
   integer(kind=ni)    :: nrecs, myid, mpro
   integer(kind=ni)    :: comm_glob
   integer(kind=ni)    :: i, j, k, kp, jp
   type(t_domdcomp)    :: p_domdcomp
   character(len=400)  :: dirpath1, dirpath2, stol
   character(len=400)  :: outfile, filenum
   character(len=500)  :: arg
   integer(kind=ni)    :: fileunit
   integer(kind=ni), dimension(2) :: sa1, sa2
   real(kind=nr)       :: wtol = 1.0E-15_nr
   real(kind=nr), dimension(:,:), allocatable      :: qa1, qa2
   integer(kind=ni), dimension(:,:),   allocatable :: lio
   
   call p_start

   call read_input_driver(nio, mbk, ndata)
   
   myid = p_get_process_ID()
   mpro = p_get_n_processes() - 1
   comm_glob = p_get_global_comm()
   
   inquire(iolength=ll) real(1.0,kind=ieee32); nrecs=ll

   dirpath1 = "../run/validation/"
   dirpath2 = "../run/validation/"

   ! Parse command line arguments
   ix = 1
   do while (ix <= command_argument_count())
      call get_command_argument(ix, arg)
      select case (arg)
         case("-h","--help")
            if (myid == 0) call print_help()
            call p_stop
         case("-d","--directories")
            ix = ix + 1
            call get_command_argument(ix, dirpath1)
            ix = ix + 1
            call get_command_argument(ix, dirpath2)
            ix = ix + 1
         case("-t", "--tolerance")
            ix = ix + 1
            call get_command_argument(ix,stol)
            read(stol, *) wtol
            ix = ix + 1
         case default
            if (myid == 0) write(*,*) "wrong argument"
            call p_stop
      end select
   end do

   if (myid == 0) then
      print *, "Command line arguments:"
      print *, "directory 1: ", trim(dirpath1)
      print *, "directory 2: ", trim(dirpath2)
      print *, "tolerance  : ", wtol
   end if
   
   ! Initialize domain decomposition
   call p_domdcomp%allocate(mbk,mpro)
   call p_domdcomp%read(mbk,lmodel_role)
   call p_domdcomp%init(mbk)
   
   ! Initialize IO
   call allocate_io_memory(mbk, ndata)
   call output_init(p_domdcomp, mbk, ndata)
   
   allocate(qa1(0:p_domdcomp%lmx,5), qa2(0:p_domdcomp%lmx,5))
   allocate(lio(0:p_domdcomp%let,0:p_domdcomp%lze))
   
   do k=0,p_domdcomp%lze
      kp = k * ( p_domdcomp%leto + 1 ) * ( p_domdcomp%lxio + 1 )
      do j=0,p_domdcomp%let
         jp = j * ( p_domdcomp%lxio + 1 )
         lio(j,k) = jp + kp
      end do
   end do
   
   ! Read restart files
   call read_restart_file(p_domdcomp, qa1, lio, dts, dte, timo, ndt, n, dt, &
                          pathbase=dirpath1)
   
   call read_restart_file(p_domdcomp, qa2, lio, dts, dte, timo, ndt, n, dt, &
                          pathbase=dirpath2)
   

   ! Compare arrays
   write (filenum,'(I1)') myid
   outfile = "output/restart"//trim(filenum)//".dat"
   open(newunit = fileunit, &
        file    = trim(outfile), &
        status  = 'replace')
   sa1 = SHAPE(qa1)
   sa2 = SHAPE(qa2)
   if (all(sa1 == sa2)) then
      do j = 1,sa1(2)
         do i = 1,sa1(1)
            if (abs((REAL(qa1(i,j),KIND=NR)-REAL(qa2(i,j),KIND=NR))/(REAL(qa1(i,j),KIND=NR))) > wtol) then
               write(fileunit, '(a,i8,a,i8,a,e15.7,a,e15.7,a,e15.7,a,e15.7)')                        &
                     " i=",i," j=",j," v1=",qa1(i,j)," v2=",qa2(i,j),                      &
                     " abs=",abs(REAL(qa1(i,j),KIND=NR)-REAL(qa2(i,j),KIND=NR)),                       &
                     " rel=",abs((REAL(qa1(i,j),KIND=NR)-REAL(qa2(i,j),KIND=NR))/REAL(qa1(i,j),KIND=NR))
            end if
         end do
      end do
   end if
   close(fileunit)
   
   write(*,'(a,i8,a,i8,a,i8)') "Compared ", sa1(1), " elements for ", sa1(2), " variables in block ", p_domdcomp%mb
   
   call p_stop
   
   contains
   
   subroutine print_help()
      print '(a, /)', 'command-line options:'
      print '(a)',    '  -d, --directories     two directories with results to compare'
      print '(a)',    '  -t, --tolerance       relative difference tolerance'
      print '(a, /)', '  -h, --help            print usage information and exit'
   end subroutine print_help
   
end program dcanard