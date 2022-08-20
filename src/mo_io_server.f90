MODULE mo_io_server
   use mo_kind,                    ONLY : ni
   use mo_mpi,                     ONLY : p_get_process_ID, p_model2io
   use mo_model_ioserver_exchange, ONLY : get_cmd, send_cmd
   
   implicit none
   private
   integer(kind=ni), parameter :: CMD_IO_EXIT=1
   public :: io_server_loop, io_server_start, io_server_stop
   contains

   SUBROUTINE io_server_start(nthick, lmodel_role)
      integer(kind=ni), intent(inout) :: nthick
      LOGICAL, INTENT(IN)             :: lmodel_role

      call p_model2io(model=nthick, server=nthick, root=0, lmodel_role=lmodel_role)
      
   END SUBROUTINE io_server_start

   SUBROUTINE io_server_stop
      CALL send_cmd(CMD_IO_EXIT)
   END SUBROUTINE io_server_stop

   SUBROUTINE io_server_loop
      integer(kind=ni)    :: i_cmd
      integer(kind=ni)    :: myid

      myid      = p_get_process_ID()
      event_loop: DO
         i_cmd =  get_cmd()
         SELECT CASE(i_cmd)
         
         CASE(CMD_IO_EXIT)
            if (myid == 0) write(*,*) "(remote_stepon): CMD_IO_EXIT"
            EXIT event_loop

         CASE DEFAULT
            if (myid == 0) write(*,*) "(remote_stepon): UNSUPPORTED CMD"

         END SELECT
      END DO event_loop

   END SUBROUTINE io_server_loop

END MODULE mo_io_server