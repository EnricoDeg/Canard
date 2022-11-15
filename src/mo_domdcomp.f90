!*****
!***** DOMAIN DECOMPOSITION MODULE
!*****

MODULE mo_domdcomp
   use mo_kind,       ONLY : ni, nr
   use mo_parameters, ONLY : one
   use mo_mpi,        ONLY : p_null_req, p_irecv, p_isend, p_waitall,    &
                           & p_recv, p_send, p_bcast, p_get_n_processes, &
                           & p_get_process_ID
   use mo_utils,      ONLY : indx3
   IMPLICIT NONE
   PRIVATE

   integer(kind=ni), private, parameter                    :: ljpl=100

   type, public :: t_domdcomp
      ! private vars
      integer(kind=ni), private, dimension(3,0:1)             :: mmcd
      integer(kind=ni), private, dimension(3)                 :: ijkp
      integer(kind=ni), private, dimension(0:3,3,0:1)         :: ijp
      integer(kind=ni), private, dimension(0:3,3,0:1)         :: ijl
      integer(kind=ni), private, dimension(0:11,ljpl)         :: jlcd
      integer(kind=ni), private, dimension(0:7,ljpl)          :: jpcd
      integer(kind=ni), private, dimension(0:7)               :: njp
      integer(kind=ni), private, dimension(0:11)              :: njl

      integer(kind=ni), private, dimension(:,:,:), allocatable :: nbbc, mbcd
      integer(kind=ni), private, dimension(:,:), allocatable   :: jjp
      integer(kind=ni), private, dimension(:,:), allocatable   :: jjl
      integer(kind=ni), private, dimension(:), allocatable     :: imjp
      integer(kind=ni), private, dimension(:), allocatable     :: jptag
      integer(kind=ni), private, dimension(:), allocatable     :: imjl
      integer(kind=ni), private, dimension(:), allocatable     :: jltag

      ! public vars
      integer(kind=ni), public                              :: mb
      integer(kind=ni), public                              :: lmx
      integer(kind=ni), public                              :: lxio, leto, lzeo
      integer(kind=ni), public                              :: lxi, let, lze
      integer(kind=ni), public, dimension(3)                :: nbsize
      integer(kind=ni), public, dimension(3,3)              :: ijk
      integer(kind=ni), public, dimension(3,0:1)            :: nbc, mcd
      integer(kind=ni), public, dimension(:),   allocatable :: lxim, letm, lzem
      integer(kind=ni), public, dimension(:),   allocatable :: lximb, letmb, lzemb
      integer(kind=ni), public, dimension(:),   allocatable :: mo
      integer(kind=ni), public, dimension(:,:), allocatable :: nbpc

   contains

      ! private subroutines / functions
      procedure, private :: trimm
      procedure, private :: idsd3

      ! public subroutines / functions
      procedure, public :: allocate => allocate_domdcomp
      procedure, public :: init => domdcomp_init
      procedure, public :: read => domdcomp_read_input
      procedure, public :: search_point
      procedure, public :: search_line
      procedure, public :: average_point
      procedure, public :: average_line
   end type t_domdcomp

   CONTAINS

   subroutine allocate_domdcomp(this, nblocks, mpro)
      class(t_domdcomp), INTENT(INOUT) :: this
      integer(kind=ni), INTENT(IN)     :: nblocks
      integer(kind=ni), INTENT(IN)     :: mpro
      integer(kind=ni)                 :: lpp, lqq

      lpp=8*nblocks+7
      lqq=12*nblocks+11
      allocate(this%nbbc(0:nblocks,3,0:1), this%mbcd(0:nblocks,3,0:1))
      allocate(this%imjp(0:lpp), this%jptag(0:lpp), this%imjl(0:lqq), this%jltag(0:lqq))
      allocate(this%jjp(0:lpp,0:ljpl))
      allocate(this%jjl(0:lqq,0:ljpl))
      allocate(this%lximb(0:nblocks), this%letmb(0:nblocks), this%lzemb(0:nblocks))
      allocate(this%mo(0:nblocks), this%nbpc(0:nblocks,3))
      allocate(this%lxim(0:mpro), this%letm(0:mpro), this%lzem(0:mpro))

   end subroutine allocate_domdcomp

   subroutine domdcomp_init(this, nblocks)
      class(t_domdcomp), INTENT(INOUT) :: this
      integer(kind=ni), intent(in) :: nblocks
      integer(kind=ni) :: ipk, jpk, mmk, nnk, nstart, nend
      integer(kind=ni) :: llk, mpk, lpk, mak, lk, mm, mp, itag
      integer(kind=ni) :: myid, mpro

      myid = p_get_process_ID()
      mpro = p_get_n_processes() - 1

      this%mo(0) = 0
      do mmk = 1,nblocks
         this%mo(mmk) = this%mo(mmk-1) + this%nbpc(mmk-1,1) * &
                                         this%nbpc(mmk-1,2) * &
                                         this%nbpc(mmk-1,3)
      end do
      
      do mmk = 0,nblocks
         if ( myid >= this%mo(mmk) ) then
            this%mb = mmk
         end if
      end do

      this%lxio = this%lximb(this%mb)
      this%leto = this%letmb(this%mb)
      this%lzeo = this%lzemb(this%mb)

      ! domdcomp initialize
      this%ijkp(1) = mod( myid - this%mo(this%mb),    &
                                 this%nbpc(this%mb,1))
      this%ijkp(2) = mod((myid - this%mo(this%mb)) /  &
                                 this%nbpc(this%mb,1), this%nbpc(this%mb,2))
      this%ijkp(3) = mod((myid - this%mo(this%mb)) /  &
                                (this%nbpc(this%mb,1) * this%nbpc(this%mb,2)), &
                                 this%nbpc(this%mb,3))

      do nnk = 1,3
         nstart = mod(nnk,3)+1
         nend   = mod(nstart,3)+1
         do ipk = 0,1
            mm = this%mbcd(this%mb,nnk,ipk)
            if ( mm == -1 ) then
               this%mmcd(nnk,ipk) = -1
            else
               this%mmcd(nnk,ipk) = this%idsd3( (1 - ipk) * (this%nbpc(mm,nnk) - 1), &
                                                this%ijkp(nstart),                   &
                                                this%ijkp(nend),                     &
                                                mm,                             &
                                                nnk ) 
            end if
         end do
      end do

      do nnk = 1,3
         select case(nnk)
         case (1)
            llk = this%lxio
            mpk = 1
         case (2)
            llk = this%leto
            mpk = this%nbpc(this%mb,1)
         case (3)
            llk = this%lzeo
            mpk = this%nbpc(this%mb,1) * this%nbpc(this%mb,2)
         end select
         lpk = this%ijkp(nnk)
         mak = this%nbpc(this%mb,nnk)
         if( mak == 1 ) then
            lk = llk
            this%nbc(nnk,0) = this%nbbc(this%mb,nnk,0)
            this%nbc(nnk,1) = this%nbbc(this%mb,nnk,1)
            this%mcd(nnk,0) = this%mmcd(nnk,0)
            this%mcd(nnk,1) = this%mmcd(nnk,1)
         end if
         if( mak >= 2 ) then
            if( lpk == 0 ) then
               lk = llk - ( ( llk + 1 ) / mak ) * ( mak - 1 )
               this%nbc(nnk,0) = this%nbbc(this%mb,nnk,0)
               this%nbc(nnk,1) = 40
               this%mcd(nnk,0) = this%mmcd(nnk,0)
               this%mcd(nnk,1) = myid+mpk
            end if
            if ( lpk > 0 .and. lpk < mak-1 ) then
               lk = ( llk + 1 ) / mak - 1
               this%nbc(nnk,0) = 40
               this%nbc(nnk,1) = 40
               this%mcd(nnk,0) = myid-mpk
               this%mcd(nnk,1) = myid+mpk
            end if
            if ( lpk == mak - 1 ) then
               lk = ( llk + 1 ) / mak - 1
               this%nbc(nnk,0) = 40
               this%nbc(nnk,1) = this%nbbc(this%mb,nnk,1)
               this%mcd(nnk,0) = myid - mpk
               this%mcd(nnk,1) = this%mmcd(nnk,1)
            end if
         end if
         select case(nnk)
         case (1)
            this%lxi = lk
         case (2)
            this%let = lk
         case (3)
            this%lze = lk
         end select
      end do

      if ( myid == 0 ) then
         this%lxim(0) = this%lxi
         this%letm(0) = this%let
         this%lzem(0) = this%lze
         do mp = 1,mpro
            itag = 1
            call p_recv(this%lxim(mp), mp, itag)
            itag = 2
            call p_recv(this%letm(mp), mp, itag)
            itag = 3
            call p_recv(this%lzem(mp), mp, itag)
         end do
      else
         itag = 1
         call p_send(this%lxi, 0, itag)
         itag = 2
         call p_send(this%let, 0, itag)
         itag = 3
         call p_send(this%lze, 0, itag)
      end if
      call p_bcast(this%lxim(:), 0)
      call p_bcast(this%letm(:), 0)
      call p_bcast(this%lzem(:), 0)

      this%lmx = ( this%lxi + 1 ) * ( this%let + 1 ) * ( this%lze + 1 ) - 1

      this%ijk(1,1) = this%lxi
      this%ijk(2,1) = this%let
      this%ijk(3,1) = this%lze
      this%ijk(1,2) = this%let
      this%ijk(2,2) = this%lze
      this%ijk(3,2) = this%lxi
      this%ijk(1,3) = this%lze
      this%ijk(2,3) = this%lxi
      this%ijk(3,3) = this%let

      this%nbsize(:) = ( this%ijk(2,:) + 1 ) * ( this%ijk(3,:) + 1 )

   END SUBROUTINE domdcomp_init

   SUBROUTINE domdcomp_read_input(this, nblocks, lmodel_role)
      class(t_domdcomp), INTENT(INOUT) :: this
      integer, INTENT(IN)              :: nblocks
      logical, intent(in)              :: lmodel_role
      character(16) :: cinput
      integer :: mmk, nnk
      integer(kind=ni) :: rc, fu
      integer(kind=ni), dimension(0:nblocks,3) :: nbpc
      integer(kind=ni), dimension(0:nblocks)   :: lximb
      integer(kind=ni), dimension(0:nblocks)   :: letmb
      integer(kind=ni), dimension(0:nblocks)   :: lzemb
      integer(kind=ni), dimension(0:nblocks,3,2)   :: nbbc
      integer(kind=ni), dimension(0:nblocks,3,2)   :: mbcd

      namelist /nml_domdcomp/ nbpc, lximb, letmb, lzemb, nbbc, mbcd
      namelist /nml_aio/ nbpc, lximb, letmb, lzemb, nbbc, mbcd
      
      
      if (lmodel_role) then
         open (action='read', file='input.canard', iostat=rc, newunit=fu)
         read (nml=nml_domdcomp, iostat=rc, unit=fu)
         close(fu)
      else
         open (action='read', file='input.canard', iostat=rc, newunit=fu)
         read (nml=nml_aio, iostat=rc, unit=fu)
         close(fu)
      end if
      
      this%nbpc  = nbpc
      this%lximb = lximb
      this%letmb = letmb
      this%lzemb = lzemb
      this%nbbc  = nbbc
      this%mbcd  = mbcd
   
   END SUBROUTINE

   SUBROUTINE search_point(this, mbki)
      class(t_domdcomp), INTENT(INOUT) :: this
      integer(kind=ni), intent(IN)     :: mbki
      integer(kind=ni) :: is, ie, kk, jj, mm, mp, kp, nn
      integer(kind=ni) :: ll, l, jp, ip, ii, i, j, k
      integer(kind=ni) :: ks, ke, lp

      lp = 8 * mbki + 7
      this%imjp(:)  =  0
      this%jjp(:,:) = -1
      do ip = 0,1
         this%ijp(:,1,ip) = (/ 0, 2, 4, 6 /) + ip
         this%ijp(:,2,ip) = (/ 0, 4, 1, 5 /) + 2 * ip
         this%ijp(:,3,ip) = (/ 0, 1, 2, 3 /) + 4 * ip
      end do
      do i = 0,lp
         this%jjp(i,0) = i
      end do
      do mm = 0,mbki
         do nn = 1,3
            do ip = 0,1
               mp = this%mbcd(mm,nn,ip)
               if ( mp /= -1 ) then
                  do i = 0,3
                     is = 8 * mm + this%ijp(i,nn,ip)
                     ie = 8 * mp + this%ijp(i,nn,1-ip)
                     this%imjp(is) = this%imjp(is) + 1
                     this%jjp(is,this%imjp(is)) = this%jjp(ie,0)
                  end do
               end if
            end do
         end do
      end do
      do i = 0,lp
         do j = 0,lp
            ll = 0
            if ( (i-j)*this%imjp(i)*this%imjp(j) /= 0 ) then
               do k = 0,this%imjp(i)
                  do kk = 0,this%imjp(j)
                     if ( this%jjp(j,kk) - this%jjp(i,k) == 0 ) then
                        ll = ll + 1
                     end if
                  end do
               end do
               if ( ll /= 0 ) then
                  ks = this%imjp(i) + 1
                  ke = ks + this%imjp(j)
                  this%jjp(i,ks:ke) = this%jjp(j,0:ke-ks)
                  call this%trimm(this%jjp(i,:), this%imjp(i), ke)
               end if
            end if
         end do
      end do
      do i = 0,lp
         this%jptag(i) = minval(this%jjp(i,0:this%imjp(i)))
      end do
      this%njp(:)    = -1
      this%jpcd(:,:) = -1
      do i = 0,7
         ip = mod(i,2)
         jp = mod(i,4)/2
         kp = i / 4
         l = 8 * this%mb + i
         if ( this%ijkp(1) == ip * ( this%nbpc(this%mb,1) - 1 ) .and. &
              this%ijkp(2) == jp * ( this%nbpc(this%mb,2) - 1 ) .and. &
              this%ijkp(3) == kp * ( this%nbpc(this%mb,3) - 1 ) .and. &
              this%imjp(l) >= 3) then
            this%njp(i) = l
            do k = 1,this%imjp(l)
               mm = this%jjp(l,k) / 8
               j  = this%jjp(l,k) - 8 * mm
               ii = mod(j,2)
               jj = mod(j,4) / 2
               kk = j / 4
               this%jpcd(i,k) = this%idsd3( ii * ( this%nbpc(mm,1) - 1 ), &
                                       jj * ( this%nbpc(mm,2) - 1 ), &
                                       kk * ( this%nbpc(mm,3) - 1 ), mm, 1)
            end do
         end if
      end do

   END SUBROUTINE search_point

   SUBROUTINE search_line(this, mbki)
      class(t_domdcomp), INTENT(INOUT) :: this
      integer(kind=ni), intent(IN) :: mbki
      integer(kind=ni) :: ns, ne, is, ie, np, nq
      integer(kind=ni) :: kk, jj, mm, mp, nn, ll
      integer(kind=ni) :: l, jp, ip, ii, i, j, k
      integer(kind=ni) :: ks, ke, lp

      lp = 12 * mbki + 11
      this%imjl(:)  = 0
      this%jjl(:,:) = -1
      this%njl(0:3) = (/ 2, 2, 1, 1 /)
      do ip=0,1
         this%ijl(:,1,ip) = (/ 0, 1, 8, 10 /) + ip * this%njl(0:3)
         this%ijl(:,2,ip) = (/ 4, 5, 0, 2  /) + ip * this%njl(0:3)
         this%ijl(:,3,ip) = (/ 8, 9, 4, 6  /) + ip * this%njl(0:3)
      end do
      do i=0,lp
         this%jjl(i,0)=i
      end do
      do mm=0,mbki
         do nn=1,3
            do ip=0,1
               mp = this%mbcd(mm,nn,ip)
               if ( mp /= -1 ) then
                  do i=0,3
                     is = 12 * mm + this%ijl(i,nn,ip)
                     ie = 12 * mp + this%ijl(i,nn,1-ip)
                     this%imjl(is) = this%imjl(is) + 1
                     this%jjl(is,this%imjl(is)) = this%jjl(ie,0)
                  end do
               end if
            end do
         end do
      end do
      do i=0,lp
         do j=0,lp
            ll = 0
            if ( (i-j) * this%imjl(i) * this%imjl(j) /= 0 ) then
               do k=0,this%imjl(i)
                  do kk=0,this%imjl(j)
                     if ( this%jjl(j,kk) - this%jjl(i,k) == 0 ) then
                        ll = ll + 1
                     end if
                  end do
               end do
               if ( ll /= 0 ) then
                  ks = this%imjl(i) + 1
                  ke = ks + this%imjl(j)
                  this%jjl(i,ks:ke) = this%jjl(j,0:ke-ks)
                  call this%trimm(this%jjl(i,:), this%imjl(i), ke)
               end if
            end if
         end do
      end do
      do i=0,lp
         this%jltag(i) = minval(this%jjl(i,0:this%imjl(i)))
      end do
      this%njl(:)    = -1
      this%jlcd(:,:) = -1
      do i=0,11
         nn = i / 4 + 1
         ns = mod(nn,3) + 1
         ne = mod(ns,3) + 1
         ip = mod(i,4) / 2
         jp = mod(mod(i,4),2)
         l  = 12 * this%mb + i
         if( this%ijkp(nn) == ip * ( this%nbpc(this%mb,nn) - 1 ) .and. &
             this%ijkp(ns) == jp * ( this%nbpc(this%mb,ns) - 1 ) .and. &
             this%imjl(l)  >= 2 ) then
            this%njl(i) = l
            do k=1,this%imjl(l)
               mm = this%jjl(l,k) / 12
               j  = this%jjl(l,k) - 12 * mm
               np = j / 4 + 1
               nq = mod(np,3) + 1
               ii = mod(j,4) / 2
               jj = mod(mod(j,4),2)
               this%jlcd(i,k) = this%idsd3( ii * ( this%nbpc(mm,np) - 1 ), &
                                       jj * ( this%nbpc(mm,nq) - 1 ), &
                                       this%ijkp(ne), mm, np )
            end do
         end if
      end do

   END SUBROUTINE search_line

   SUBROUTINE average_point(this,qak)
      class(t_domdcomp), INTENT(INOUT) :: this
      real(kind=nr),    intent(inout), dimension(0:this%lmx,5) :: qak
      integer(kind=ni) :: kp, l, jp, ip, ii, i, j, lp, lq, itag
      real(kind=nr)    :: fctr
      real(kind=nr),dimension(:,:),allocatable :: rr

      allocate(rr(0:this%lmx,3))
      call p_null_req
      lp = -5
      lq = -5
      do i=0,7
         ii = this%njp(i)
         if ( ii /= -1 ) then
            itag = this%jptag(ii)
            ip = mod(i,2)
            jp = mod(i,4) / 2
            kp = i / 4
            lp = lp + 5
            l  = indx3( ip * this%lxi, jp * this%let, &
                        kp * this%lze, 1, this%lxi, this%let)
            rr(lp:lp+4,1)=qak(l,1:5)
            do j=1,this%imjp(ii)
               lq = lq + 5
               call p_isend(rr(lp:lp+4,1), this%jpcd(i,j), itag)
               call p_irecv(rr(lq:lq+4,2), this%jpcd(i,j), itag)
            end do
         end if
      end do
      call p_waitall
      lp = -5
      lq = -5
      do i=0,7
         ii = this%njp(i)
         if ( ii /= -1 ) then
            ip = mod(i,2)
            jp = mod(i,4) / 2
            kp = i / 4
            lp = lp + 5
            do j=1,this%imjp(ii)
               lq = lq + 5
               rr(lp:lp+4,1) = rr(lp:lp+4,1) + rr(lq:lq+4,2)
            end do
            fctr = one / ( this%imjp(ii) + 1 )
            l = indx3( ip * this%lxi, jp * this%let, &
                       kp * this%lze, 1, this%lxi, this%let)
            qak(l,1:5) = fctr * rr(lp:lp+4,1)
         end if
      end do
      deallocate(rr)

   END SUBROUTINE average_point


   SUBROUTINE average_line(this,qak)
      class(t_domdcomp), INTENT(INOUT) :: this
      real(kind=nr),    intent(inout), dimension(0:this%lmx,5) :: qak
      integer(kind=ni) :: nn, ll, l, jp, ip, ii, i, j, k
      integer(kind=ni) :: js, je, ks, ke, lp, lq, itag
      real(kind=nr)    :: fctr
      real(kind=nr),dimension(:,:),allocatable :: rr

      allocate(rr(0:this%lmx,3))
      call p_null_req
      lp = 0
      lq = 0
      do i=0,11
         ii = this%njl(i)
         if ( ii /= -1 ) then
            itag = this%jltag(ii)
            nn = i / 4 + 1
            ip = mod(i,4) / 2
            jp = mod(mod(i,4),2)
            ks = lp
            ke = ks + 5 * ( this%ijk(3,nn) + 1 ) - 1
            lp = lp + ke - ks + 1
            do k=0,this%ijk(3,nn)
               ll = ks + 5 * k
               l  = indx3( ip * this%ijk(1,nn), jp * this%ijk(2,nn), &
                           k, nn, this%lxi, this%let )
               rr(ll:ll+4,1) = qak(l,1:5)
            end do
            do j=1,this%imjl(ii)
               js = lq
               je = js + ke - ks
               lq = lq + je - js + 1
               call p_isend(rr(ks:ke,1), this%jlcd(i,j), itag)
               call p_irecv(rr(js:je,2), this%jlcd(i,j), itag)
            end do
         end if
      end do
      call p_waitall
      lp = 0
      lq = 0
      do i=0,11
         ii = this%njl(i)
         if ( ii /= -1 ) then
            nn = i / 4 + 1
            ip = mod(i,4) / 2
            jp = mod(mod(i,4),2)
            ks = lp
            ke = ks + 5 * ( this%ijk(3,nn) + 1 ) - 1
            lp = lp + ke - ks + 1
            do j=1,this%imjl(ii)
               js = lq
               je = js + ke - ks
               lq = lq + je - js + 1
               rr(ks:ke,1) = rr(ks:ke,1) + rr(js:je,2)
            end do
            fctr = one / ( this%imjl(ii) + 1 )
            do k=0,this%ijk(3,nn)
               ll = ks + 5 * k
               l  = indx3( ip * this%ijk(1,nn), jp * this%ijk(2,nn), &
                           k, nn, this%lxi, this%let)
               qak(l,1:5) = fctr * rr(ll:ll+4,1)
            end do
         end if
      end do
      deallocate(rr)

   END SUBROUTINE average_line

!===== SUBROUTINE FOR TRIMMING DOWN TO EFFECTIVE ARRAY

   subroutine trimm(this,itrim,lene,lens)
      class(t_domdcomp), INTENT(INOUT) :: this
      integer(kind=ni),dimension(0:ljpl),intent(inout) :: itrim
      integer(kind=ni),intent(inout) :: lene,lens
      integer(kind=ni) :: k,kk

      do k=0,lens-1
         do kk=k+1,lens
            if ( itrim(kk) - itrim(k) == 0 ) then
               itrim(kk) = -1
            end if
         end do
      end do
      lene = lens
      do k=0,lens
         if ( itrim(k) == -1 ) then
            lene = lene - 1
            itrim(k) = itrim(k+1)
            itrim(k+1) = -1
         end if
      end do

   end subroutine trimm

!===== FUNCTION FOR PROCESSOR-CORE INDEX TRANSFORMATION WITHIN THE BLOCK IN 3D

   function idsd3(this,i,j,k,mm,nn) result(lm)
      class(t_domdcomp), INTENT(INOUT) :: this
      integer(kind=ni),intent(in) :: i,j,k,mm,nn
      integer(kind=ni) :: lm

      select case(nn)
      case(1)
         lm = this%mo(mm) + ( k * this%nbpc(mm,2) + j ) * &
              this%nbpc(mm,1) + i
      case(2)
         lm = this%mo(mm) + ( j * this%nbpc(mm,2) + i ) * &
              this%nbpc(mm,1) + k
      case(3)
         lm = this%mo(mm) + ( i * this%nbpc(mm,2) + k ) * &
              this%nbpc(mm,1) + j
      end select

   end function idsd3

END MODULE mo_domdcomp
