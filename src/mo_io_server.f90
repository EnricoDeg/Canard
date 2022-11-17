MODULE mo_io_server
   use mo_kind,                    ONLY : ni, nr, int64, ieee32
   use mo_parameters,              ONLY : zero
   use mo_mpi,                     ONLY : p_get_process_ID, p_model2io, p_get_n_processes, &
                                        & p_get_global_n_processes
   use mo_model_ioserver_exchange, ONLY : get_cmd, send_cmd, recv_model2io, send_model2io
   use mo_domdcomp,                ONLY : t_domdcomp
   use mo_io,                      ONLY : write_output_file_mb, vminmax
   
   implicit none
   private
   integer(kind=ni), parameter :: CMD_IO_EXIT           = 1
   integer(kind=ni), parameter :: CMD_IO_INIT           = 2
   integer(kind=ni), parameter :: CMD_WRITE_OUTPUT_FILE = 3

   type, public :: t_model_interface
      integer(kind=ni) :: mps
      integer(kind=ni) :: mpe
      integer(kind=ni), dimension(:), allocatable :: lxim, letm, lzem
      integer(kind=ni), dimension(:), allocatable :: lpos
   end type t_model_interface

   type, public :: t_io_server_interface
      integer(kind=ni) :: mb
   end type t_io_server_interface

   INTERFACE io_server_init
      MODULE PROCEDURE io_server_init_model
      MODULE PROCEDURE io_server_init_io
   END INTERFACE io_server_init

   INTERFACE io_server_write_output
      MODULE PROCEDURE io_server_write_output_model
      MODULE PROCEDURE io_server_write_output_io
   END INTERFACE io_server_write_output

   public :: io_server_loop, io_server_stop, io_server_init
   public :: io_server_write_output
   contains

   SUBROUTINE io_server_init_model(mbk, p_domdcomp, p_io_server_interface, mpro, lpos, lmodel_role)
      integer(kind=ni), intent(in)               :: mbk
      type(t_domdcomp), intent(in)               :: p_domdcomp
      type(t_io_server_interface), intent(INOUT) :: p_io_server_interface
      integer(kind=ni), intent(in)               :: mpro      
      integer(kind=ni), dimension(0:mpro), intent(in) :: lpos
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
         temp = lpos(i)
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
               p_model_interface%lzem(0:mpro_model), &
               p_model_interface%lpos(0:mpro_model)  )

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
         call p_model2io(model=p_model_interface%lpos(i), server=p_model_interface%lpos(i), &
                         root=0, lmodel_role=lmodel_role)
      end do      

   END SUBROUTINE io_server_init_io

   SUBROUTINE io_server_write_output_model(vara, mq, ndati, time, p_domdcomp, p_io_server_interface)
      real(kind=ieee32), intent(in), dimension(:) :: vara
      integer(kind=ni), intent(inout)             :: mq
      integer(kind=ni), intent(inout)             :: ndati
      real(kind=nr), intent(inout)                :: time
      type(t_domdcomp), intent(in)                :: p_domdcomp
      type(t_io_server_interface), intent(in)     :: p_io_server_interface

      integer(kind=ni) :: ljs, lje

      ljs = 0
      lje = mq * ( p_domdcomp%lmx + 1 ) - 1
      call send_cmd(CMD_WRITE_OUTPUT_FILE)
      call p_model2io(model=mq, server=mq, root=0, lmodel_role=.true.)
      call p_model2io(model=ndati, server=ndati, root=0, lmodel_role=.true.)
      if (mq == 5) then
         call p_model2io(model=time, server=time, root=0, lmodel_role=.true.)
      end if
      call send_model2io(vara, p_io_server_interface%mb, ljs, lje)

   END SUBROUTINE io_server_write_output_model

   SUBROUTINE io_server_write_output_io(mbk, ndata, times, p_domdcomp, p_model_interface)
      integer(kind=ni), intent(in)                  :: mbk
      integer(kind=ni), intent(in)                  :: ndata
      real(kind=nr), intent(inout), dimension(0:ndata) :: times
      type(t_domdcomp), intent(in)                  :: p_domdcomp
      type(t_model_interface), intent(in)           :: p_model_interface

      integer(kind=ni)    :: myid, ltomb, mp, mps, mpe, nn
      integer(kind=int64) :: llmb
      integer(kind=ni)    :: lisi, ljsi, j, k, m, ndati, nnn, mq
      integer(kind=ni),  allocatable, dimension(:) :: lis, lie
      real(kind=ieee32), allocatable, dimension(:) :: vara, varb
      
      call p_model2io(model=mq, server=mq, root=0, lmodel_role=.false.)
      call p_model2io(model=ndati, server=ndati, root=0, lmodel_role=.false.)
      if (mq == 5) then
         call p_model2io(model=times(ndati), server=times(ndati), root=0, lmodel_role=.false.)
      end if

      ltomb = ( p_domdcomp%lxio + 1 ) * ( p_domdcomp%leto + 1 ) * &
              ( p_domdcomp%lzeo + 1 )
      llmb  = mq * ltomb - 1
      mps = p_model_interface%mps
      mpe = p_model_interface%mpe
      allocate(vara(0:llmb), varb(0:llmb))
      allocate(lis(0:mpe-mps), lie(0:mpe-mps))

      ! receive data from IO client
      lis(0) = 0
      lie(0) = mq * (( p_model_interface%lxim(mps) + 1 ) * &
                   ( p_model_interface%letm(mps) + 1 ) * &
                   ( p_model_interface%lzem(mps) + 1 ) ) - 1
      do mp=mps+1,mpe
         lis(mp-mps) = lie(mp-mps-1) + 1
         lie(mp-mps) = lis(mp-mps  ) + mq * ( p_model_interface%lxim(mp) + 1 ) * &
                                           ( p_model_interface%letm(mp) + 1 ) * &
                                           ( p_model_interface%lzem(mp) + 1 ) - 1
      end do
      call recv_model2io(vara, mps, mpe, lis, lie)

      ! reorganize data
      lisi=0
      do mp=mps,mpe
         do m=1,mq
            do k=0,p_model_interface%lzem(mp)
               do j=0,p_model_interface%letm(mp)
                  ljsi = p_model_interface%lpos(mp) + ( m - 1 ) * ltomb + k * ( p_domdcomp%leto + 1 ) *     &
                                                                              ( p_domdcomp%lxio + 1 ) + j * &
                                                                              ( p_domdcomp%lxio + 1 )
                  varb(ljsi:ljsi+p_model_interface%lxim(mp)) = vara(lisi:lisi+p_model_interface%lxim(mp))
                  lisi = lisi + p_model_interface%lxim(mp) + 1
               end do
            end do
         end do
      end do
      do nn=1,mq
         if (mq == 5) then
            nnn=3+5*ndati+nn
         else
            nnn = nn
         end if
         call vminmax(p_domdcomp, varb((nn-1)*(ltomb):nn*(ltomb)-1), nnn)
       end do
      call write_output_file_mb(p_domdcomp, mbk, ndata, times, llmb, varb, mq, ndati)
      deallocate(lis, lie)
      deallocate(vara, varb)
      
   END SUBROUTINE io_server_write_output_io

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

         CASE(CMD_WRITE_OUTPUT_FILE)
            if (myid == 0) write(*,*) "(remote_stepon): CMD_WRITE_OUTPUT_FILE"
            call io_server_write_output(mbk, ndata, times, p_domdcomp, &
                                        p_model_interface)
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
