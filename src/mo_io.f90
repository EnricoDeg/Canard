!*****
!***** INPUT/OUTPUT MODULE
!*****

MODULE mo_io
   use mo_kind,       ONLY : ni, nr, int64, ieee32, int32, ieee64
   use mo_parameters, ONLY : zero
   use mo_domdcomp,   ONLY : t_domdcomp
   use mo_mpi,        ONLY : p_get_n_processes, p_get_process_ID,   &
                           & p_barrier, p_recv, p_send, p_null_req, &
                           & p_irecv, p_waitall
   use mo_utils,      ONLY : indx3

   IMPLICIT NONE
   PUBLIC

   integer(kind=int64), private, dimension(:),   allocatable :: lhmb
   real(kind=ieee32),   private, dimension(:,:), allocatable :: varm
   real(kind=ieee32),   private, dimension(:),   allocatable :: varmin, varmax
   real(kind=ieee32),   private, dimension(:),   allocatable :: vara, varb
   character(13),       private, dimension(:),   allocatable :: ctecplt, cthead
   character(4),        private, dimension(:),   allocatable :: cfilet

   character(1),        private, dimension(0:4)              :: cno
   integer(kind=ni),    private, dimension(0:4)              :: no
   real(kind=nr),       private, dimension(5)                :: cha, dha
   character(19),       private                              :: crestart
   character(4),        private, dimension(:), allocatable   :: czonet
   integer(kind=ni),    private, dimension(:), allocatable   :: lpos

   character(16), public :: cgrid

   CONTAINS

   SUBROUTINE read_input_main(nts, nrestart, &
                          cfl, tmax, ltimer)
      integer(kind=ni), intent(out) :: nts
      integer(kind=ni), intent(out) :: nrestart
      real(kind=nr),    intent(out) :: cfl
      real(kind=nr),    intent(out) :: tmax
      logical,          intent(out) :: ltimer
      integer(kind=ni) :: rc, fu
      character(16) :: cinput

      namelist /nml_canard/ nts,nrestart,cfl,tmax,ltimer

      open (action='read', file='input.canard', iostat=rc, newunit=fu)
      read (nml=nml_canard, iostat=rc, unit=fu)
      close(fu)
   
   END SUBROUTINE read_input_main

   SUBROUTINE read_input_driver(nio, mbk, ndata)
      integer(kind=ni), intent(out) :: nio
      integer(kind=ni), intent(out) :: mbk
      integer(kind=ni), intent(out) :: ndata
      character(16) :: cinput
      integer(kind=ni) :: rc, fu

      namelist /nml_driver/ nio,mbk,ndata

      open (action='read', file='input.canard', iostat=rc, newunit=fu)
      read (nml=nml_driver, iostat=rc, unit=fu)
      close(fu)

   END SUBROUTINE read_input_driver

   SUBROUTINE allocate_io_memory(mbk, ndata)
      integer(kind=ni), intent(in) :: mbk
      integer(kind=ni), intent(in) :: ndata
      integer(kind=ni) :: ll, mpro

      mpro = p_get_n_processes() -1
      
      ll=3+5*(ndata+1)
      allocate(cfilet(-1:ndata),     &
               ctecplt(-1:ndata),varm(0:1,0:mpro),  &
               varmin(ll),varmax(ll),cthead(0:mbk), &
               czonet(0:mbk),lhmb(0:mbk))
      allocate(lpos(0:mpro))
   END SUBROUTINE allocate_io_memory

!===== SUBROUTINE FOR INITIALIZING OUTPUT

   SUBROUTINE output_init(p_domdcomp, mbk, ndata, lpos_temp)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      integer(kind=ni), intent(in) :: mbk
      integer(kind=ni), intent(in) :: ndata
      integer(kind=ni), intent(inout), dimension(:), optional :: lpos_temp

      integer(kind=ni) :: n
      integer(kind=ni) :: mm, mp, kp, jp, i, j, k
      integer(kind=ni) :: myid

      myid = p_get_process_ID()
      cfilet(-1) = 'grid'
      do n=0,ndata
         no(2) = n / 100
         no(1) = mod(n,100) / 10
         no(0) = mod(n,10)
         cno   = achar(no+48)
         cfilet(n)  = 'n'//cno(2)//cno(1)//cno(0)
      end do
      do n=-1,ndata
         ctecplt(n) = 'data/'//cfilet(n)//'.plt'
      end do
      do mm=0,mbk
         no(2) = mm / 100
         no(1) = mod(mm,100) / 10
         no(0) = mod(mm,10)
         cno   = achar(no+48)
         czonet(mm) = 'z'//cno(2)//cno(1)//cno(0)
         cthead(mm) = 'data/'//czonet(mm)//'.plt'
      end do
      cgrid    = 'misc/grid'//czonet(p_domdcomp%mb)//'.dat'
      crestart = 'misc/restart'//czonet(p_domdcomp%mb)//'.dat'

      do mm=0,mbk
         lpos(p_domdcomp%mo(mm)) = 0
         do i=1,p_domdcomp%nbpc(mm,1)-1
            mp       = p_domdcomp%mo(mm) + i
            lpos(mp) = lpos(mp-1) + p_domdcomp%lxim(mp-1) + 1
         end do
         jp = p_domdcomp%nbpc(mm,1)
         do j=1,p_domdcomp%nbpc(mm,2)-1
            do i=0,p_domdcomp%nbpc(mm,1)-1
               mp       = p_domdcomp%mo(mm) + j * jp + i
               lpos(mp) = lpos(mp-jp) + ( p_domdcomp%lximb(mm) + 1 ) * &
                                        ( p_domdcomp%letm(mp-jp) + 1 )
            end do
         end do
         kp = p_domdcomp%nbpc(mm,1) * p_domdcomp%nbpc(mm,2)
         do k=1,p_domdcomp%nbpc(mm,3)-1
            do j=0,p_domdcomp%nbpc(mm,2)-1
               do i=0,p_domdcomp%nbpc(mm,1)-1
                  mp       = p_domdcomp%mo(mm) + k * kp + j * jp + i
                  lpos(mp) = lpos(mp-kp) + ( p_domdcomp%lximb(mm) + 1 ) * &
                                           ( p_domdcomp%letmb(mm) + 1 ) * &
                                           ( p_domdcomp%lzem(mp-kp) + 1 )
               end do
            end do
         end do
      end do

      if (present(lpos_temp)) lpos_temp(:) = lpos(:)

   END SUBROUTINE output_init

!===== SUBROUTINE FOR READING GRID 

   SUBROUTINE read_grid_parallel(p_domdcomp, ssk, lio)
      type(t_domdcomp), intent(IN)  :: p_domdcomp
      real(kind=nr),    intent(out) :: ssk(0:p_domdcomp%lmx,3)
      integer(kind=ni), dimension(0:p_domdcomp%let,0:p_domdcomp%lze), intent(in) :: lio
      integer(kind=ni) :: i, j, k, l, lq, lp
      integer(kind=ni) :: ll, nrecd
      integer(kind=ni) :: myid

      myid = p_get_process_ID()
      inquire(iolength=ll) real(1.0,kind=ieee64); nrecd=ll

      open(9,file=cgrid,access='direct',form='unformatted',recl=3*nrecd,status='old')
      lp = lpos(myid)
      do k=0,p_domdcomp%lze
         do j=0,p_domdcomp%let
            lq = lp + lio(j,k)
            do i=0,p_domdcomp%lxi
               l = indx3(i, j, k, 1, p_domdcomp%lxi, p_domdcomp%let)
               read(9,rec=lq+i+1) ssk(l,:)
            end do
         end do
      end do
      close(9)
      call p_barrier
      if ( myid == p_domdcomp%mo(p_domdcomp%mb) ) then
         open(9,file=cgrid,status='old')
         close(9,status='delete')
      end if

   END SUBROUTINE read_grid_parallel

!===== SUBROUTINE FOR READING RESTART FILE

   SUBROUTINE read_restart_file(p_domdcomp, qa, lio, dts, dte, timo, ndt, n, dt, pathbase)
      type(t_domdcomp), intent(IN)    :: p_domdcomp
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(inout) :: qa
      integer(kind=ni), dimension(0:p_domdcomp%let,0:p_domdcomp%lze), intent(in) :: lio
      real(kind=nr), intent(INOUT)    :: dts, dte
      real(kind=nr), intent(inout)    :: timo
      integer(kind=ni), intent(inout) :: ndt
      integer(kind=ni), intent(inout) :: n
      real(kind=nr), intent(inout)    :: dt
      character(len=*), intent(in), optional :: pathbase
      character(len=500) :: path_restart
      integer(kind=ni) :: lp, i, j, k, lq, l
      integer(kind=ni) :: ll, nrecd
      integer(kind=ni) :: myid

      if (present(pathbase)) then
         path_restart=trim(pathbase)//"/"//trim(crestart)
      else
         path_restart=trim(crestart)
      endif

      myid = p_get_process_ID()
      inquire(iolength=ll) real(1.0,kind=ieee64); nrecd=ll

      open(9,file=trim(path_restart),access='direct',form='unformatted',recl=5*nrecd,status='old')
      read(9,rec=1) cha(:)
      read(9,rec=2) dha(:)
      n    = cha(1)
      ndt  = cha(2)
      dt   = cha(3)
      dts  = cha(4)
      dte  = cha(5)
      timo = dha(1)
      lp   = lpos(myid)+2
      do k=0,p_domdcomp%lze
         do j=0,p_domdcomp%let
            lq = lp + lio(j,k)
            do i=0,p_domdcomp%lxi
               l = indx3(i, j, k, 1, p_domdcomp%lxi, p_domdcomp%let)
               read(9,rec=lq+i+1) qa(l,:)
            end do
         end do
      end do
      close(9)

   END SUBROUTINE read_restart_file

!===== SUBROUTINE FOR WRITING RESTART FILE

   SUBROUTINE write_restart_file(p_domdcomp, qa, lio, dts, dte, timo, ndt, n, dt)
      type(t_domdcomp), intent(IN)    :: p_domdcomp
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(inout) :: qa
      integer(kind=ni), dimension(0:p_domdcomp%let,0:p_domdcomp%lze), intent(in) :: lio
      real(kind=nr), intent(INOUT)    :: dts, dte
      real(kind=nr), intent(inout)    :: timo
      integer(kind=ni), intent(inout) :: ndt
      integer(kind=ni), intent(inout) :: n
      real(kind=nr), intent(inout)    :: dt
      integer(kind=ni) :: lp, i, j, k, lq, l
      integer(kind=ni) :: ll, nrecd
      integer(kind=ni) :: myid

      myid = p_get_process_ID()
      inquire(iolength=ll) real(1.0,kind=ieee64); nrecd=ll

      if ( myid == p_domdcomp%mo(p_domdcomp%mb) ) then
         open(9,file=crestart,access='direct',form='unformatted',recl=5*nrecd,status='replace')
         cha(:) = (/ real(n,kind=nr), real(ndt,kind=nr), dt, dts, dte /)
         dha(:) = (/ timo, zero, zero, zero, zero /)
         write(9,rec=1) cha(:)
         write(9,rec=2) dha(:)
         close(9)
      end if
      call p_barrier
      open(9,file=crestart,access='direct',form='unformatted',recl=5*nrecd,status='old')
      lp = lpos(myid) + 2
      do k=0,p_domdcomp%lze
         do j=0,p_domdcomp%let
            lq = lp + lio(j,k)
            do i=0,p_domdcomp%lxi
               l = indx3(i, j, k, 1, p_domdcomp%lxi, p_domdcomp%let)
               write(9,rec=lq+i+1) qa(l,:)
            end do
         end do
      end do
      close(9)

   END SUBROUTINE write_restart_file

!===== SUBROUTINE FOR WRITING OUTPUT FILE (IO SERVER)

   SUBROUTINE write_output_file_mb(p_domdcomp, mbk, ndata, times, llmb, vart, mq, n)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      integer(kind=ni), intent(in) :: mbk
      integer(kind=ni), intent(in) :: ndata
      real(kind=nr), dimension(0:ndata), intent(in) :: times
      integer(kind=int64), intent(in) :: llmb
      real(kind=ieee32), dimension(0:llmb), intent(inout) :: vart
      integer(kind=ni), intent(in) :: mq
      integer(kind=ni), intent(in) :: n

      real(kind=ieee32), dimension(:), allocatable :: varbm
      integer(kind=ni) :: mm, mp, j, k, mps, mpe, m, lmpi, lhf, itag, lh
      integer(kind=int64) :: llmo, lis, lie, ljs, lje
      integer(kind=ni) :: ltomb
      integer(kind=ni) :: ll, nrecs
      integer(kind=ni) :: myid, nn

      myid = p_get_process_ID()
      inquire(iolength=ll) real(1.0,kind=ieee32); nrecs=ll

      ltomb = ( p_domdcomp%lxio + 1 ) * &
              ( p_domdcomp%leto + 1 ) * &
              ( p_domdcomp%lzeo + 1 )
      lje  = -1
      ljs = lje + 1
      lje = ljs + mq * ( p_domdcomp%lmx + 1 ) - 1

      open(9,file=cthead(p_domdcomp%mb),access='stream',form='unformatted',status='replace')
      call techead(p_domdcomp, 9, n, p_domdcomp%mb, lh, mq, mbk, ndata, times)
      allocate(vara(0:lh+llmb))
      read(9,pos=1) vara(0:lh-1)
      close(9,status='delete')
      lhmb(p_domdcomp%mb) = lh + llmb + 1
      vara(lh:lh+llmb)    = vart(:)
      if(p_domdcomp%mb==0) then !---------------------
         do mm=1,mbk
            itag=2
            call p_recv(lhmb(mm), p_domdcomp%mo(mm), itag)
         end do
         llmo = sum(lhmb(:)) - 1
         allocate(varbm(0:llmo))
         lis = 0
         lie = lhmb(p_domdcomp%mb)-1
         varbm(lis:lie) = vara(:)
         do mm=1,mbk
            lis  = lie + 1
            lie  = lis + lhmb(mm) - 1
            lmpi = lie - lis + 1
            itag = 3
            call p_recv(varbm(lis:lie), p_domdcomp%mo(mm), itag, lmpi)
         end do
         open(0,file=ctecplt(n),status='unknown')
         close(0,status='delete') ! 'replace' not suitable as 'recl' may vary
         open(0,file=ctecplt(n),access='direct',form='unformatted',recl=nrecs*(llmo+1),status='new')
         write(0,rec=1) varbm(:)
         close(0)
         deallocate(varbm)
      else !-------------------
         itag = 2
         call p_send(lhmb(p_domdcomp%mb), p_domdcomp%mo(0), itag)
         lmpi = lhmb(p_domdcomp%mb)
         itag = 3
         call p_send(vara(:), p_domdcomp%mo(0), itag, lmpi)
      end if !-----------
      deallocate(vara)

   END SUBROUTINE write_output_file_mb

!===== SUBROUTINE FOR WRITING OUTPUT FILE

   SUBROUTINE write_output_file(p_domdcomp, mbk, ndata, times, nlmx, vart, n, mq)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      integer(kind=ni), intent(in) :: mbk
      integer(kind=ni), intent(in) :: ndata
      real(kind=nr), dimension(0:ndata), intent(in) :: times
      integer(kind=int64), intent(in) :: nlmx
      real(kind=ieee32), dimension(0:nlmx), intent(inout) :: vart
      integer(kind=ni), intent(in) :: n
      integer(kind=ni), intent(in) :: mq
      integer(kind=ni) :: mm, mp, j, k, mps, mpe, m, lmpi, lhf, itag, lh
      integer(kind=int64) :: llmo, llmb, lis, lie, ljs, lje
      integer(kind=ni) :: ltomb
      integer(kind=ni) :: ll, nrecs
      integer(kind=ni) :: myid

      myid = p_get_process_ID()
      inquire(iolength=ll) real(1.0,kind=ieee32); nrecs=ll

      ltomb = ( p_domdcomp%lxio + 1 ) * &
              ( p_domdcomp%leto + 1 ) * &
              ( p_domdcomp%lzeo + 1 )
      lje   = -1
         llmb = mq * ltomb - 1
         allocate(vara(0:llmb),varb(0:llmb))
         ljs = lje + 1
         lje = ljs + mq * ( p_domdcomp%lmx + 1 ) - 1
         if ( myid == p_domdcomp%mo(p_domdcomp%mb) ) then !===========================
            mps = p_domdcomp%mo(p_domdcomp%mb)
            mpe = mps + p_domdcomp%nbpc(p_domdcomp%mb,1) * &
                        p_domdcomp%nbpc(p_domdcomp%mb,2) * &
                        p_domdcomp%nbpc(p_domdcomp%mb,3)-1
            lis = 0
            lie = mq * ( p_domdcomp%lmx + 1 ) - 1
            vara(lis:lie) = vart(ljs:lje)
            do mp=mps+1,mpe
               lis = lie + 1
               lie = lis + mq * ( p_domdcomp%lxim(mp) + 1 ) * &
                                ( p_domdcomp%letm(mp) + 1 ) * &
                                ( p_domdcomp%lzem(mp) + 1 ) - 1
               lmpi = lie - lis + 1
               itag = 1
               call p_recv(vara(lis:lie), mp, itag, lmpi)
            end do
            lis=0
            do mp=mps,mpe
               do m=1,mq
                  do k=0,p_domdcomp%lzem(mp)
                     do j=0,p_domdcomp%letm(mp)
                        ljs = lpos(mp) + ( m - 1 ) * ltomb + k * ( p_domdcomp%leto + 1 ) *     &
                                                                 ( p_domdcomp%lxio + 1 ) + j * &
                                                                 ( p_domdcomp%lxio + 1 )
                        varb(ljs:ljs+p_domdcomp%lxim(mp)) = vara(lis:lis+p_domdcomp%lxim(mp))
                        lis = lis + p_domdcomp%lxim(mp) + 1
                     end do
                  end do
               end do
            end do
            deallocate(vara)
            call write_output_file_mb(p_domdcomp, mbk, ndata, times, llmb, varb, mq, n)
         else !===========================
            lmpi = lje - ljs + 1
            itag = 1
            call p_send(vart(ljs:lje), p_domdcomp%mo(p_domdcomp%mb), itag, lmpi)
            deallocate(vara)
         end if !=========================
         deallocate(varb)

   END SUBROUTINE write_output_file

!===== SUBROUTINE FOR GENERATING TECPLOT DATA FILE

   subroutine techead(p_domdcomp, nf, n, mb, lh, mq, mbk, ndata, times)
      type(t_domdcomp), intent(IN)    :: p_domdcomp
      integer(kind=ni), intent(in)    :: nf,n,mb
      integer(kind=ni), intent(inout) :: lh
      integer(kind=ni), intent(in)    :: mq
      integer(kind=ni), intent(in) :: mbk
      integer(kind=ni), intent(in) :: ndata
      real(kind=nr), dimension(0:ndata), intent(in) :: times
      integer(kind=ni)                :: mm, m, nn
      character(16) :: cinput

      lh=0
      if(mb==0) then
         write(nf,pos=4*lh+1) '#!TDV112'
         lh=lh+2
         write(nf,pos=4*lh+1) 1
         lh=lh+1 ! Header Section
         write(nf,pos=4*lh+1) int(min(n+2,2),kind=int32)
         lh=lh+1 ! File Type (0 = Full / 1 = Grid / 2 = Solution)
         cinput=cfilet(n)
         call strio(nf,lh,cinput) ! File Title
         write(nf,pos=4*lh+1) int(mq,kind=int32)
         lh=lh+1 ! Number of Variables
         if(n==-1) then
            cinput='x'
            call strio(nf,lh,cinput)
             cinput='y'
            call strio(nf,lh,cinput)
            cinput='z'
            call strio(nf,lh,cinput)
         else
            cinput='rho'
            call strio(nf,lh,cinput)
            cinput='u'
            call strio(nf,lh,cinput)
            cinput='v'
            call strio(nf,lh,cinput)
            cinput='w'
            call strio(nf,lh,cinput)
            cinput='p'
            call strio(nf,lh,cinput)
         end if
         do mm=0,mbk
            write(nf,pos=4*lh+1) 299.0
            lh=lh+1 ! Zone Marker
            cinput=czonet(mm)
            call strio(nf,lh,cinput) ! Zone Name
            write(nf,pos=4*lh+1) -1
            lh=lh+1 ! Parent Zone
            write(nf,pos=4*lh+1) int(mm+1,kind=int32)
            lh=lh+1 ! Strand ID
            write(nf,pos=4*lh+1) real(times(max(n,0)),kind=ieee64)
            lh=lh+2 ! Solution Time (Double)
            write(nf,pos=4*lh+1) -1
            lh=lh+1 ! (Not used. Set to -1.)
            write(nf,pos=4*lh+1) 0
            lh=lh+1 ! Zone Type
            write(nf,pos=4*lh+1) 0
            lh=lh+1 ! Specify Var Location
            write(nf,pos=4*lh+1) 0
            lh=lh+1 ! Raw Local 1-to-1 Face Neighbours Suppliled
            write(nf,pos=4*lh+1) 0
            lh=lh+1 ! Number of Miscellaneous Face Neighbour Connections
            write(nf,pos=4*lh+1) int(p_domdcomp%lximb(mm)+1,kind=int32)
            lh=lh+1 ! IMax
            write(nf,pos=4*lh+1) int(p_domdcomp%letmb(mm)+1,kind=int32)
            lh=lh+1 ! JMax
            write(nf,pos=4*lh+1) int(p_domdcomp%lzemb(mm)+1,kind=int32)
            lh=lh+1 ! KMax
            write(nf,pos=4*lh+1) 0
            lh=lh+1 ! No Auxillary Data Pairs
         end do
         write(nf,pos=4*lh+1) 357.0
         lh=lh+1 ! End of Header Marker
      end if
      write(nf,pos=4*lh+1) 299.0
      lh=lh+1 ! Zone Marker
      do m=1,mq
         write(nf,pos=4*lh+1) 1
         lh=lh+1 ! 1 = Float / 2 = Double
      end do
      write(nf,pos=4*lh+1) 0
      lh=lh+1 ! No Passive Variables
      write(nf,pos=4*lh+1) 0
      lh=lh+1 ! No Variable Sharing
      write(nf,pos=4*lh+1) -1
      lh=lh+1 ! Zero Based Zone Number to Share
      do m=1,mq
         nn=max(3+5*n,0)+m
         write(nf,pos=4*lh+1) real(varmin(nn),kind=ieee64)
         lh=lh+2 ! Minimum Value (Double) of Variables
         write(nf,pos=4*lh+1) real(varmax(nn),kind=ieee64)
         lh=lh+2 ! Maximum Value (Double) of Variables
      end do

   end subroutine techead

!===== SUBROUTINE FOR CHARACTER STRING CONVERSION

   subroutine strio(nfile,lh,cinput)

      integer(kind=ni),intent(in) :: nfile
      integer(kind=ni),intent(inout) :: lh
      character(16),intent(in) :: cinput
      integer(kind=ni) :: ll

      do ll=1,len_trim(cinput)
         write(nfile,pos=4*lh+1) ichar(cinput(ll:ll))
         lh=lh+1
      end do
      write(nfile,pos=4*lh+1) 0
      lh=lh+1

   end subroutine strio

!===== SUBROUTINE FOR FINDING VARIABLE MIN/MAX VALUES FOR TECPLOT DATA FILE

   subroutine vminmax(p_domdcomp, varr, nn)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      real(kind=ieee32), dimension(0:p_domdcomp%lmx), intent(in) :: varr
      integer(kind=ni),intent(in) :: nn
      integer(kind=ni) :: mp, mps, mpe, itag
      integer(kind=ni) :: myid

      myid = p_get_process_ID()
      varmin(nn)   = minval(varr)
      varmax(nn)   = maxval(varr)
      varm(:,myid) = (/varmin(nn),varmax(nn)/)

      call p_null_req
      itag = nn
      if ( myid == p_domdcomp%mo(p_domdcomp%mb) ) then
         mps = p_domdcomp%mo(p_domdcomp%mb)
         mpe = mps + p_domdcomp%nbpc(p_domdcomp%mb,1) * &
                     p_domdcomp%nbpc(p_domdcomp%mb,2) * &
                     p_domdcomp%nbpc(p_domdcomp%mb,3) - 1
         do mp=mps+1,mpe
            call p_irecv(varm(:,mp), mp, itag, 2)
         end do
         call p_waitall
         varmin(nn) = minval(varm(0,mps:mpe))
         varmax(nn) = maxval(varm(1,mps:mpe))
      else
         call p_send(varm(:,myid), p_domdcomp%mo(p_domdcomp%mb), itag, 2)
      end if

   end subroutine vminmax

END MODULE mo_io
