program main
   use mo_kind,   ONLY : ni
   use mo_mpi,    ONLY : p_start, p_stop, p_set_work_comm
   use mo_io,     ONLY : read_input_driver
   use mo_canard, ONLY : canard_driver
   use mo_aio,    ONLY : aio_driver
   implicit none

   integer(kind=ni)    :: nio
   integer(kind=ni)    :: mbk
   integer(kind=ni)    :: nbody
   integer(kind=ni)    :: ndata
   logical             :: lmodel_role
   logical             :: laio

   call p_start

   call read_input_driver(nio, mbk, nbody, ndata)

   if (nio > 0 .and. nio /= (mbk+1)) then
      write(*,*) "Currently only one IO server for each block supported"
      STOP 1
   end if
   
   call p_set_work_comm(nio, lmodel_role, laio)

   if (lmodel_role) then
      call canard_driver(laio, lmodel_role, mbk, nbody, ndata)
   else
      call aio_driver(lmodel_role, mbk, nbody, ndata)
   end if

   call p_stop

end program main