module mo_aio
   use mo_kind,       ONLY : ni
   use mo_mpi,        ONLY : p_get_n_processes, p_get_process_ID, p_barrier,       &
                           & p_get_global_comm
   implicit none
   public
   contains

   subroutine aio_driver
      integer(kind=ni)    :: comm_glob, myid, npro

      myid      = p_get_process_ID()
      npro      = p_get_n_processes()
      comm_glob = p_get_global_comm()

      write(*,*) "IO server number ", myid

      call p_barrier(comm_glob)

   end subroutine aio_driver


end module mo_aio