module mo_aio
   use mo_kind,       ONLY : ni
   use mo_mpi,        ONLY : p_get_n_processes, p_get_process_ID, p_barrier,       &
                           & p_get_global_comm
   use mo_domdcomp,   ONLY : t_domdcomp
   use mo_io,         ONLY : allocate_io_memory, output_init
   implicit none
   public
   contains

   subroutine aio_driver(lmodel_role, mbk, nbody, ndata)
      logical, intent(in)          :: lmodel_role
      integer(kind=ni), intent(in) :: mbk
      integer(kind=ni), intent(in) :: nbody
      integer(kind=ni), intent(in) :: ndata

      integer(kind=ni)    :: comm_glob, myid, mpro
      integer(kind=ni)    :: nthick ! need to be sent from canard driver
      type(t_domdcomp)    :: p_domdcomp

      myid      = p_get_process_ID()
      mpro      = p_get_n_processes() - 1
      comm_glob = p_get_global_comm()

!===== DOMAIN DECOMPOSITION INITIALIZATION

      call p_domdcomp%allocate(mbk,mpro)
      call p_domdcomp%read(lmodel_role)
      call p_domdcomp%init(mbk, nthick, nbody)

!===== WRITING START POSITIONS IN OUTPUT FILE

      call allocate_io_memory(mbk, ndata)
      call output_init(p_domdcomp, mbk, ndata)





      
      call p_barrier(comm_glob)

   end subroutine aio_driver


end module mo_aio