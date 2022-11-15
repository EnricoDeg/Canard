MODULE mo_io_server
   use mo_kind,                    ONLY : ni, nr, int64, ieee32
   use mo_parameters,              ONLY : zero
   use mo_mpi,                     ONLY : p_get_process_ID, p_model2io, p_get_n_processes, &
                                        & p_get_global_n_processes
   use mo_model_ioserver_exchange, ONLY : get_cmd, send_cmd, recv_model2io, send_model2io
   use mo_domdcomp,                ONLY : t_domdcomp
   use mo_io,                      ONLY : write_output_grid, vminmax
   
   implicit none
   private
   integer(kind=ni), parameter :: CMD_IO_EXIT         = 1
   integer(kind=ni), parameter :: CMD_IO_INIT         = 2
   integer(kind=ni), parameter :: CMD_WRITE_GRID_FILE = 3

   type, public :: t_model_interface
      integer(kind=ni) :: mps
      integer(kind=ni) :: mpe
      integer(kind=ni), dimension(:), allocatable :: lxim, letm, lzem
   end type t_model_interface

   type, public :: t_io_server_interface
      integer(kind=ni) :: mb
   end type t_io_server_interface

   INTERFACE io_server_init
      MODULE PROCEDURE io_server_init_model
      MODULE PROCEDURE io_server_init_io
   END INTERFACE io_server_init

   INTERFACE io_server_write_grid
      MODULE PROCEDURE io_server_write_grid_model
      MODULE PROCEDURE io_server_write_grid_io
   END INTERFACE io_server_write_grid

   public :: io_server_loop, io_server_stop, io_server_init
   public :: io_server_write_grid
   contains

   SUBROUTINE io_server_init_model(mbk, p_domdcomp, p_io_server_interface, lmodel_role)
      integer(kind=ni), intent(in)               :: mbk
      type(t_domdcomp), intent(in)               :: p_domdcomp
      type(t_io_server_interface), intent(INOUT) :: p_io_server_interface
      LOGICAL, INTENT(IN)                        :: lmodel_role

      integer(kind=ni) :: myid, npro, i, j, temp

      CALL send_cmd(CMD_IO_INIT)

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

      do i=0,npro-1
         temp = p_domdcomp%lxim(i)
         call p_model2io(model=temp, server=temp, root=0, lmodel_role=lmodel_role)
         temp = p_domdcomp%letm(i)
         call p_model2io(model=temp, server=temp, root=0, lmodel_role=lmodel_role)
         temp = p_domdcomp%lzem(i)
         call p_model2io(model=temp, server=temp, root=0, lmodel_role=lmodel_role)
      end do

   END SUBROUTINE io_server_init_model

   SUBROUTINE io_server_init_io(mbk, p_domdcomp, p_model_interface, lmodel_role)
      integer(kind=ni), intent(in)                  :: mbk
      type(t_domdcomp), intent(in)                  :: p_domdcomp
      type(t_model_interface), intent(inout)        :: p_model_interface
      LOGICAL, INTENT(IN)                           :: lmodel_role

      integer(kind=ni), dimension(:), allocatable   :: mo
      integer(kind=ni), dimension(:,:), allocatable :: nbpc
      integer(kind=ni) :: i, j, myid, npro, npro_gl, mpro_model

      myid = p_get_process_ID()
      npro = p_get_n_processes()
      npro_gl = p_get_global_n_processes()
      mpro_model = npro_gl - npro - 1
      allocate(p_model_interface%lxim(0:mpro_model), &
               p_model_interface%letm(0:mpro_model), &
               p_model_interface%lzem(0:mpro_model)  )

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

      do i=0,mpro_model
         call p_model2io(model=p_model_interface%lxim(i), server=p_model_interface%lxim(i), &
                         root=0, lmodel_role=lmodel_role)
         call p_model2io(model=p_model_interface%letm(i), server=p_model_interface%letm(i), &
                         root=0, lmodel_role=lmodel_role)
         call p_model2io(model=p_model_interface%lzem(i), server=p_model_interface%lzem(i), &
                         root=0, lmodel_role=lmodel_role)
      end do

      write(*,*) "io server ", myid, ": mps = ", p_model_interface%mps, &
                                   " -- mpe = ", p_model_interface%mpe

   END SUBROUTINE io_server_init_io

   SUBROUTINE io_server_write_grid_model(vara, p_domdcomp, p_io_server_interface)
      real(kind=ieee32), intent(in), dimension(:) :: vara
      type(t_domdcomp), intent(in)                :: p_domdcomp
      type(t_io_server_interface), intent(in)     :: p_io_server_interface

      integer(kind=ni) :: ljs, lje

      ljs = 0
      lje = 3 * ( p_domdcomp%lmx + 1 ) - 1
      call send_cmd(CMD_WRITE_GRID_FILE)
      call send_model2io(vara, p_io_server_interface%mb, ljs, lje)

   END SUBROUTINE io_server_write_grid_model

   SUBROUTINE io_server_write_grid_io(mbk, ndata, times, p_domdcomp, p_model_interface, lmodel_role)
      integer(kind=ni), intent(in)                  :: mbk
      integer(kind=ni), intent(in)                  :: ndata
      real(kind=nr), intent(in), dimension(0:ndata) :: times
      type(t_domdcomp), intent(in)                  :: p_domdcomp
      type(t_model_interface), intent(in)           :: p_model_interface
      logical, intent(in)                           :: lmodel_role

      integer(kind=ni)    :: myid, ltomb, llmb, mp, mps, mpe, nn
      integer(kind=int64) :: nlmx
      integer(kind=ni), allocatable, dimension(:) :: lis, lie
      real(kind=ieee32), allocatable, dimension(:)    :: vara

      ltomb = ( p_domdcomp%lxio + 1 ) * ( p_domdcomp%leto + 1 ) * &
              ( p_domdcomp%lzeo + 1 )
      llmb  = 3 * ltomb - 1
      nlmx  = int(llmb, kind=int64)
      mps = p_model_interface%mps
      mpe = p_model_interface%mpe
      allocate(vara(0:llmb))
      allocate(lis(0:mpe-mps), lie(0:mpe-mps))
      lis(0) = 0
      lie(0) = 3 * ( p_model_interface%lxim(mps) + 1 ) * &
                   ( p_model_interface%letm(mps) + 1 ) * &
                   ( p_model_interface%lzem(mps) + 1 ) - 1
      do mp=mps+1,mpe
         lis(mp-mps) = lie(mp-mps-1) + 1
         lie(mp-mps) = lis(mp-mps-1) + 3 * ( p_model_interface%lxim(mp) + 1 ) * &
                                           ( p_model_interface%letm(mp) + 1 ) * &
                                           ( p_model_interface%lzem(mp) + 1 ) - 1
      end do
      call recv_model2io(vara, mps, mpe, lis, lie)
      call vminmax(p_domdcomp, vara((nn-1)*(p_domdcomp%lmx+1):nn*(p_domdcomp%lmx+1)-1), nn)
      do nn=1,3
         call write_output_grid(p_domdcomp, mbk, ndata, times, nlmx, vara)
      end do
      deallocate(lis, lie)
      deallocate(vara)

   END SUBROUTINE io_server_write_grid_io

   SUBROUTINE io_server_stop
      CALL send_cmd(CMD_IO_EXIT)
   END SUBROUTINE io_server_stop

   SUBROUTINE io_server_loop(mbk, ndata, p_domdcomp, p_model_interface, lmodel_role)
      integer(kind=ni), intent(in)           :: mbk
      integer(kind=ni), intent(in)           :: ndata
      type(t_domdcomp), intent(in)           :: p_domdcomp
      type(t_model_interface), intent(inout) :: p_model_interface
      logical, intent(in)                    :: lmodel_role

      integer(kind=ni)    :: i_cmd
      integer(kind=ni)    :: myid
      real(kind=nr), dimension(:), allocatable :: times

      myid      = p_get_process_ID()
      allocate(times(0:ndata))
      times(:) = zero
      event_loop: DO
         i_cmd =  get_cmd()
         SELECT CASE(i_cmd)

         CASE(CMD_IO_INIT)
            if (myid == 0) write(*,*) "(remote_stepon): CMD_IO_INIT"
            call io_server_init(mbk, p_domdcomp, p_model_interface, lmodel_role)

         CASE(CMD_WRITE_GRID_FILE)
            if (myid == 0) write(*,*) "(remote_stepon): CMD_WRITE_GRID_FILE"
            call io_server_write_grid(mbk, ndata, times, p_domdcomp, &
                                      p_model_interface, lmodel_role)

         CASE(CMD_IO_EXIT)
            if (myid == 0) write(*,*) "(remote_stepon): CMD_IO_EXIT"
            EXIT event_loop

         CASE DEFAULT
            if (myid == 0) write(*,*) "(remote_stepon): UNSUPPORTED CMD"

         END SELECT
      END DO event_loop
      deallocate(times)

   END SUBROUTINE io_server_loop

END MODULE mo_io_server
