MODULE mo_io_server
   use mo_kind,                    ONLY : ni, nr
   use mo_mpi,                     ONLY : p_get_process_ID, p_model2io, p_get_n_processes
   use mo_model_ioserver_exchange, ONLY : get_cmd, send_cmd
   use mo_domdcomp,                ONLY : t_domdcomp
   
   implicit none
   private
   integer(kind=ni), parameter :: CMD_IO_EXIT=1

   type, public :: t_model_interface
      integer(kind=ni) :: mps
      integer(kind=ni) :: mpe
   end type t_model_interface

   type, public :: t_io_server_interface
      integer(kind=ni) :: mb
   end type t_io_server_interface

   INTERFACE io_server_init
      MODULE PROCEDURE io_server_init_model
      MODULE PROCEDURE io_server_init_io
   END INTERFACE io_server_init

   public :: io_server_loop, io_server_start, io_server_stop, io_server_init
   contains

   SUBROUTINE io_server_start(nthick, lmodel_role)
      integer(kind=ni), intent(inout) :: nthick
      logical, intent(in)             :: lmodel_role

      call p_model2io(model=nthick, server=nthick, root=0, lmodel_role=lmodel_role)

   END SUBROUTINE io_server_start

   SUBROUTINE io_server_init_model(mbk, p_domdcomp, p_io_server_interface, lmodel_role)
      integer(kind=ni), intent(in)               :: mbk
      type(t_domdcomp), intent(in)               :: p_domdcomp
      type(t_io_server_interface), intent(INOUT) :: p_io_server_interface
      LOGICAL, INTENT(IN)                        :: lmodel_role

      integer(kind=ni) :: myid, npro, i, j, temp

      myid = p_get_process_ID()
      npro = p_get_n_processes()
      p_io_server_interface%mb = npro + p_domdcomp%mb
      do i=0,mbk
         temp = p_domdcomp%mo(i)
         call p_model2io(model=temp, server=temp, root=0, lmodel_role=lmodel_role)
      end do
      do j=1,3
         do i=0,mbk
            temp = p_domdcomp%nbpc(i,j)
            call p_model2io(model=temp, server=temp, root=0, lmodel_role=lmodel_role)
         end do
      end do
      write(*,*) "model proc", myid, ": mb io server = ", p_io_server_interface%mb

   END SUBROUTINE io_server_init_model

   SUBROUTINE io_server_init_io(mbk, p_domdcomp, p_model_interface, lmodel_role)
      integer(kind=ni), intent(in)                  :: mbk
      type(t_domdcomp), intent(in)                  :: p_domdcomp
      type(t_model_interface), intent(inout)        :: p_model_interface
      LOGICAL, INTENT(IN)                           :: lmodel_role

      integer(kind=ni), dimension(:), allocatable   :: mo
      integer(kind=ni), dimension(:,:), allocatable :: nbpc
      integer(kind=ni) :: i, j, myid

      myid = p_get_process_ID()

      allocate(mo(0:mbk), nbpc(0:mbk,3))
      do i=0,mbk
         call p_model2io(model=mo(i), server=mo(i), root=0, lmodel_role=lmodel_role)         
      end do
      do j=1,3
         do i=0,mbk
            call p_model2io(model=nbpc(i,j), server=nbpc(i,j), root=0, lmodel_role=lmodel_role)            
         end do
      end do
      p_model_interface%mps = mo(p_domdcomp%mb)
      p_model_interface%mpe = p_model_interface%mps + nbpc(p_domdcomp%mb,1) * &
                                                      nbpc(p_domdcomp%mb,2) * &
                                                      nbpc(p_domdcomp%mb,3) - 1
      deallocate(mo, nbpc)

      write(*,*) "io server ", myid, ": mps = ", p_model_interface%mps, &
                                   " -- mpe = ", p_model_interface%mpe

   END SUBROUTINE io_server_init_io

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
