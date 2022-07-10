!*****
!***** INPUT/OUTPUT MODULE
!*****

MODULE mo_io
   use mo_kind,       ONLY : ni, nr, int64, ieee32, int32, ieee64
   use mo_parameters, ONLY : zero
   use mo_vars,       ONLY : lpos,             &
                           & mbk,   &
                           & n,                    &
                           & cnnode, cgrid, cdata,                    &
                           & nrecd, &
                           & dt, nrecs,              &
                           & varr, qa, vart
   use mo_grid,       ONLY : lio
   use mo_domdcomp,   ONLY : t_domdcomp
   use mo_mpi,        ONLY : mpro, myid, p_barrier, p_recv, p_send,         &
                             p_null_req, p_irecv, p_waitall
   use mo_utils,      ONLY : indx3
   use mo_gridgen,    ONLY : domh, doml0, doml1, let0, lxi0, lxi1,          &
                           & lxi2, lze0, nthick, skew, smg,         &
                           & smgvr, span, spx, szth0, szth1, wlea,          &
                           & wlew
   IMPLICIT NONE
   PUBLIC

   integer(kind=ni),    private, dimension(:),   allocatable :: idsgnl, lsgnl
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
   character(4),        private, dimension(:),allocatable    :: czonet

   CONTAINS

   SUBROUTINE read_inputo(nts, nscrn, ndata, ndatafl, ndataav, nrestart, &
                          cfl, dto, tsam, tmax, nkrk)
      integer(kind=ni), intent(out) :: nts, nscrn, ndata, ndatafl, ndataav
      integer(kind=ni), intent(out) :: nrestart
      real(kind=nr),    intent(out) :: cfl, dto
      real(kind=nr),    intent(out) :: tsam, tmax
      integer(kind=ni), intent(out) :: nkrk
      character(16) :: cinput

      open(9,file='inputo.dat',status='old')
      read(9,*) cinput,mbk
      read(9,*) cinput,nts
      read(9,*) cinput,nscrn
      read(9,*) cinput,ndata,ndatafl,ndataav
      read(9,*) cinput,nkrk
      read(9,*) cinput,nrestart
      read(9,*) cinput,cfl
      read(9,*) cinput,tmax,tsam
      read(9,*) cinput,dto
      close(9)

   END SUBROUTINE read_inputo

   SUBROUTINE read_inputp(nbody)
      integer(kind=ni), intent(out) :: nbody
      character(16) :: cinput

      open(9,file='inputp.dat',status='old')
      read(9,*) cinput,lxi0,lxi1,lxi2
      read(9,*) cinput,let0
      read(9,*) cinput,lze0
      read(9,*) cinput,nbody,nthick
      read(9,*) cinput,smg,smgvr
      read(9,*) cinput,doml0,doml1,domh
      read(9,*) cinput,span
      read(9,*) cinput,wlew,wlea
      read(9,*) cinput,szth0,szth1
      read(9,*) cinput,skew,spx
      close(9)

      allocate(idsgnl(0:lze0),lsgnl(0:lze0))

   END SUBROUTINE read_inputp

   SUBROUTINE allocate_io_memory(ndata)
      integer(kind=ni), intent(in) :: ndata
      integer(kind=ni) :: ll
      
      ll=3+5*(ndata+1)
      allocate(cfilet(-1:ndata),     &
               ctecplt(-1:ndata),varm(0:1,0:mpro),  &
               varmin(ll),varmax(ll),cthead(0:mbk), &
               czonet(0:mbk),lhmb(0:mbk))
   END SUBROUTINE allocate_io_memory

   SUBROUTINE output_init(p_domdcomp, ndata)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      integer(kind=ni), intent(in) :: ndata
      integer(kind=ni) :: mm, mp, kp, jp, i, j, k

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

      no(4) = myid / 10000
      no(3) = mod(myid,10000) / 1000
      no(2) = mod(myid,1000) / 100
      no(1) = mod(myid,100) / 10
      no(0) = mod(myid,10)
      cno   = achar(no+48)
      cnnode = cno(4)//cno(3)//cno(2)//cno(1)//cno(0)
      cdata  = 'misc/data'//cnnode//'.dat'

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

   END SUBROUTINE output_init

   SUBROUTINE read_restart_file(p_domdcomp, dts, dte, timo, ndt)
      type(t_domdcomp), intent(IN)    :: p_domdcomp
      real(kind=nr), intent(INOUT)    :: dts, dte
      real(kind=nr), intent(inout)    :: timo
      integer(kind=ni), intent(inout) :: ndt
      integer(kind=ni) :: lp, i, j, k, lq, l

      open(9,file=crestart,access='direct',form='unformatted',recl=5*nrecd,status='old')
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

   SUBROUTINE write_restart_file(p_domdcomp, dts, dte, timo, ndt)
      type(t_domdcomp), intent(IN)    :: p_domdcomp
      real(kind=nr), intent(INOUT)    :: dts, dte
      real(kind=nr), intent(inout)    :: timo
      integer(kind=ni), intent(inout) :: ndt
      integer(kind=ni) :: lp, i, j, k, lq, l

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

   SUBROUTINE write_output_file(p_domdcomp, ndata, times)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      integer(kind=ni), intent(in) :: ndata
      real(kind=nr), dimension(0:ndata), intent(in) :: times
      integer(kind=ni) :: mm, mp, j, k, mps, mpe, m, lmpi, lhf, itag, lh
      integer(kind=int64) :: llmo, llmb, lis, lie, ljs, lje
      integer(kind=ni) :: ltomb, mq

      ltomb = ( p_domdcomp%lxio + 1 ) * &
              ( p_domdcomp%leto + 1 ) * &
              ( p_domdcomp%lzeo + 1 )
      lje   = -1
      do n=-1,ndata
         mq   = 3 + 2 * min(n+1,1)
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
            open(9,file=cthead(p_domdcomp%mb),access='stream',form='unformatted',status='replace')
            call techead(p_domdcomp, 9, n, p_domdcomp%mb, lh, mq, ndata, times)
            deallocate(vara)
            allocate(vara(0:lh+llmb))
            read(9,pos=1) vara(0:lh-1)
            close(9,status='delete')
            lhmb(p_domdcomp%mb) = lh + llmb + 1
            vara(lh:lh+llmb)    = varb(:)
            if(p_domdcomp%mb==0) then !---------------------
               do mm=1,mbk
                  itag=2
                  call p_recv(lhmb(mm), p_domdcomp%mo(mm), itag)
               end do
               llmo = sum(lhmb(:)) - 1
               deallocate(varb)
               allocate(varb(0:llmo))
               lis = 0
               lie = lhmb(p_domdcomp%mb)-1
               varb(lis:lie) = vara(:)
               do mm=1,mbk
                  lis  = lie + 1
                  lie  = lis + lhmb(mm) - 1
                  lmpi = lie - lis + 1
                  itag = 3
                  call p_recv(varb(lis:lie), p_domdcomp%mo(mm), itag, lmpi)
               end do
               open(0,file=ctecplt(n),status='unknown')
               close(0,status='delete') ! 'replace' not suitable as 'recl' may vary
               open(0,file=ctecplt(n),access='direct',form='unformatted',recl=nrecs*(llmo+1),status='new')
               write(0,rec=1) varb(:)
               close(0)
            else !-------------------
               itag = 2
               call p_send(lhmb(p_domdcomp%mb), p_domdcomp%mo(0), itag)
               lmpi = lhmb(p_domdcomp%mb)
               itag = 3
               call p_send(vara(:), p_domdcomp%mo(0), itag, lmpi)
            end if !-----------
         else !===========================
            lmpi = lje - ljs + 1
            itag = 1
            call p_send(vart(ljs:lje), p_domdcomp%mo(p_domdcomp%mb), itag, lmpi)
         end if !=========================
         deallocate(vara,varb)
      end do
 
   END SUBROUTINE write_output_file

!===== SUBROUTINE FOR GENERATING TECPLOT DATA FILE

   subroutine techead(p_domdcomp, nf, n, mb, lh, mq, ndata, times)
      type(t_domdcomp), intent(IN)    :: p_domdcomp
      integer(kind=ni), intent(in)    :: nf,n,mb
      integer(kind=ni), intent(inout) :: lh
      integer(kind=ni), intent(in)    :: mq
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

   subroutine vminmax(p_domdcomp, nn)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      integer(kind=ni),intent(in) :: nn
      integer(kind=ni) :: mp, mps, mpe, itag

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