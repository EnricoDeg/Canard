module mo_aio_driver
   use mo_kind,                    ONLY : ni
   use mo_mpi,                     ONLY : p_get_n_processes, p_get_process_ID, p_barrier,       &
                                        & p_get_global_comm
   use mo_domdcomp,                ONLY : t_domdcomp
   use mo_io,                      ONLY : allocate_io_memory, output_init
   use mo_io_server,               ONLY : io_server_loop, t_model_interface
   implicit none
   public
   contains

   subroutine aio_driver(lmodel_role, mbk, ndata)
      logical, intent(in)          :: lmodel_role
      integer(kind=ni), intent(in) :: mbk
      integer(kind=ni), intent(in) :: ndata

      integer(kind=ni)    :: comm_glob, myid, mpro
      type(t_domdcomp)    :: p_domdcomp
      type(t_model_interface) :: p_model_interface

      myid      = p_get_process_ID()
      mpro      = p_get_n_processes() - 1
      comm_glob = p_get_global_comm()

!===== DOMAIN DECOMPOSITION INITIALIZATION
      if (myid==0) write(*,*) "IO server: initialize domain decomposition"

      call p_domdcomp%allocate(mbk,mpro)
      call p_domdcomp%read(mbk, lmodel_role)
      call p_domdcomp%init(mbk)

!===== WRITING START POSITIONS IN OUTPUT FILE
      if (myid==0) write(*,*) "IO server: initialize IO"

      call allocate_io_memory(mbk, ndata)
      call output_init(p_domdcomp, mbk, ndata)

!===== MAIN LOOP
      if (myid==0) write(*,*) "IO server: starting loop"

      call io_server_loop(mbk, ndata, p_domdcomp, p_model_interface, lmodel_role)

      call p_barrier(comm_glob)

   end subroutine aio_driver


end module mo_aio_driver