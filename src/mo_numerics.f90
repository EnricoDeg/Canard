!*****
!***** 3D SOLVER MODULE
!*****

module mo_numerics
   use mo_parameters, ONLY : pi, half, zero, one, alpha01, alpha10,   &
                           & alpha, beta, beta02, beta13, alpha12,    &
                           & a01, a02, a03, a04, a10, a12, a13, a14,  &
                           & ab, aa, two, n45go,                      &
                           & BC_NON_REFLECTIVE, BC_WALL_INVISCID,     &
                           & BC_WALL_VISCOUS, BC_INTER_CURV,          &
                           & BC_INTER_STRAIGHT, BC_INTER_SUBDOMAINS,  &
                           & BC_PERIODIC
   use mo_kind,       ONLY : nr, ni
   use mo_mpi,        ONLY : p_null_req, p_isend, p_irecv, p_waitall, &
                             p_send
   use mo_utils,      ONLY : indx3, mtrxi
   implicit none
   private

   integer(kind=ni), private, parameter               :: lmd=11, lmf=11, lmp=max(lmd,lmf)
   integer(kind=ni), private, parameter, dimension(3) :: nnf = (/ 2,3,1 /)

   type, public :: t_numerics
      ! private vars
      real(kind=nr), private                             :: alphf, betf
      real(kind=nr), private, dimension(0:lmp,0:1,0:1)   :: pbci, pbco
      real(kind=nr), private, dimension(0:1,0:1)         :: pbcot
      real(kind=nr), private                             :: fa, fb, fc
      real(kind=nr), private, dimension(-2:2,0:2,0:1)    :: albef
      real(kind=nr), private, dimension(0:lmp)           :: sap
      real(kind=nr), private, dimension(0:4,0:2)         :: fbc
      integer(kind=ni), private, dimension(3,0:1,0:1)    :: ndf

      real(kind=nr), private, dimension(:,:), allocatable :: xu, yu
      real(kind=nr), private, dimension(:,:), allocatable :: xl, yl
      real(kind=nr), private :: fltk, fltrbc

      real(kind=nr), private, dimension(:,:,:), pointer :: send01, send02, send03
      real(kind=nr), private, dimension(:,:,:), pointer :: send11, send12, send13
      real(kind=nr), private, dimension(:,:,:,:), pointer :: recv01, recv02, recv03
      real(kind=nr), private, dimension(:,:,:,:), pointer :: recv11, recv12, recv13
      real(kind=nr), private, dimension(:,:,:), pointer :: send
      real(kind=nr), private, dimension(:,:,:,:), pointer :: recv
      integer(kind=ni), private, dimension(:),  allocatable :: li
      real(kind=nr), private, dimension(:),     allocatable :: sa, sb
      real(kind=nr), private, dimension(:),     allocatable :: sar, sbr

      ! public vars
      real(kind=nr), public, dimension(:,:,:), pointer :: drva1, drva2, drva3
      real(kind=nr), public, dimension(:,:,:), pointer :: drva

   contains

      ! private subroutines / functions
      procedure, private :: penta
      procedure, private :: fcbcm
      procedure, private :: fcint
      procedure, private :: sbcco

      ! public subroutines / functions
      procedure, public :: allocate => allocate_numerics
      procedure, public :: read => read_input_numerics
      procedure, public :: init_extra => init_extracoeff_bounds
      procedure, public :: init => init_penta
      procedure, public :: mpigo_1d
      procedure, public :: mpigo_2d
      generic,   public :: mpigo => mpigo_1d,mpigo_2d
      procedure, public :: deriv
      procedure, public :: filte

   end type t_numerics

   contains

   SUBROUTINE allocate_numerics(this, limk, nbsizek, lmxk)
      class(t_numerics), intent(inout)          :: this
      integer(kind=ni),intent(in)               :: limk
      integer(kind=ni),intent(in), dimension(3) :: nbsizek
      integer(kind=ni),intent(in)               :: lmxk
      integer(kind=ni) :: iik, jjk, kkk

      iik=nbsizek(1)-1
      jjk=nbsizek(2)-1
      kkk=nbsizek(3)-1

      allocate(this%xu(0:limk,3), this%yu(0:limk,3), this%xl(0:limk,2), this%yl(0:limk,2))
      allocate(this%li(0:limk),   this%sa(0:limk),   this%sb(0:limk))
      allocate(this%sar(0:lmxk),   this%sbr(0:lmxk))
      allocate(this%send01(0:iik,0:1,0:1), this%send02(0:jjk,0:1,0:1), this%send03(0:kkk,0:1,0:1))
      allocate(this%recv01(0:iik,0:1,0:1,1:5), this%recv02(0:jjk,0:1,0:1,1:5), this%recv03(0:kkk,0:1,0:1,1:5))
      allocate(this%send11(0:iik,0:2,0:1), this%send12(0:jjk,0:2,0:1), this%send13(0:kkk,0:2,0:1))
      allocate(this%recv11(0:iik,0:2,0:1,1:5), this%recv12(0:jjk,0:2,0:1,1:5), this%recv13(0:kkk,0:2,0:1,1:5))
      allocate(this%drva1(0:iik,5,0:1),    this%drva2(0:jjk,5,0:1),    this%drva3(0:kkk,5,0:1))

   END SUBROUTINE allocate_numerics

   SUBROUTINE read_input_numerics(this)
      class(t_numerics), intent(inout) :: this
      real(kind=nr) :: fltk, fltrbc      
      integer(kind=ni) :: fu, rc
      namelist /nml_numerics/ fltk, fltrbc

      open (action='read', file='input.canard', iostat=rc, newunit=fu)
      read (nml=nml_numerics, iostat=rc, unit=fu)
      close(fu)
      
      this%fltk = pi * fltk
      this%fltrbc = fltrbc

   END SUBROUTINE read_input_numerics

!===== EXTRA COEFFICIENTS FOR DOMAIN BOUNDARIES

   SUBROUTINE init_extracoeff_bounds(this)
      class(t_numerics), intent(inout) :: this
      integer(kind=ni) :: ntk, jjk, iik

      call this%fcbcm(this%fltk, this%fltrbc)
      call this%fcint(this%fltk, half, this%alphf, this%betf, this%fa, this%fb, this%fc)
      this%albef(:,0,1) = (/ zero,      zero,       one, this%alphf, this%betf /)
      this%albef(:,1,1) = (/ zero,      this%alphf, one, this%alphf, this%betf /)
      this%albef(:,2,1) = (/ this%betf, this%alphf, one, this%alphf, this%betf /)

      this%pbco(:,:,:) = zero
      this%pbci(:,:,:) = zero
      call this%sbcco
      do ntk = 0,1
         do jjk = 0,1
            iik = lmd + ntk * ( lmf - lmd )
            this%pbcot(jjk,ntk) = sum( this%pbco(0:iik,jjk,ntk) )
         end do
      end do

   END SUBROUTINE init_extracoeff_bounds

   SUBROUTINE init_penta(this, lxik, letk, lzek, nbck, lim)
      class(t_numerics), intent(inout) :: this
      integer(kind=ni), intent(in) :: lxik, letk, lzek
      integer(kind=ni), intent(in), dimension(3,0:1) :: nbck
      integer(kind=ni), intent(in) :: lim
      integer(kind=ni) :: nnk, ipk, istart, iend
      integer(kind=ni) :: nstart, nend, npk

      do nnk = 1,3
         select case(nnk)
         case(1)
            istart = 0
            iend = istart + lxik
         case(2)
            istart = lxik + 1
            iend = istart + letk
         case(3)
            istart = lxik + letk + 2
            iend = istart + lzek
         end select
         do ipk = 0,1
            npk=nbck(nnk,ipk)
            select case(npk)
            case(BC_NON_REFLECTIVE,BC_WALL_INVISCID,BC_WALL_VISCOUS,BC_INTER_CURV)
               this%ndf(nnk,ipk,0) = 0
               this%ndf(nnk,ipk,1) = 0
            case(BC_INTER_STRAIGHT,BC_INTER_SUBDOMAINS,BC_PERIODIC)
               this%ndf(nnk,ipk,0) = 1
               this%ndf(nnk,ipk,1) = 1
            end select
         end do
         nstart = this%ndf(nnk,0,0)
         nend   = this%ndf(nnk,1,0)
         call this%penta(this%xu(:,:), this%xl(:,:), istart, iend, nstart, nend, 0, lim)
         nstart = this%ndf(nnk,0,1)
         nend   = this%ndf(nnk,1,1)
         call this%penta(this%yu(:,:), this%yl(:,:), istart, iend, nstart, nend, 1, lim)
      end do
  
   END SUBROUTINE init_penta

!===== SUBROUTINE FOR CHOLESKY DECOMPOSITION OF PENTADIAGONAL MATRICES

   subroutine penta(this,xu,xl,is,ie,ns,ne,nt,lim)
      class(t_numerics), intent(inout) :: this
      integer(kind=ni),intent(in) :: is,ie,ns,ne,nt
      real(kind=nr),dimension(0:lim,3),intent(inout) :: xu
      real(kind=nr),dimension(0:lim,2),intent(inout) :: xl
      integer(kind=ni), intent(in) :: lim
      real(kind=nr),dimension(-2:2,0:2,0:1) :: albe
      real(kind=nr) :: alpho,beto
      integer(kind=ni) :: iik

      if(nt==0) then
         albe(:,0,0) = (/ zero, zero,    one, alpha01, beta02 /)
         albe(:,1,0) = (/ zero, alpha10, one, alpha12, beta13 /)
         albe(:,2,0) = (/ beta, alpha,   one, alpha,   beta /)
         alpho = alpha
         beto  = beta
         albe(:,0,1) = (/ zero, zero,  one, alpha, beta /)
         albe(:,1,1) = (/ zero, alpha, one, alpha, beta /)
         albe(:,2,1) = (/ beta, alpha, one, alpha, beta /)
      else
         albe  = this%albef
         alpho = this%alphf
         beto  = this%betf
      end if

      do iik = is,ie
         xl(iik,:) = one
         xu(iik,:) = one
      end do
      iik = is
      xu(iik,1) = one
      xu(iik,2) = albe(1,0,ns)
      xu(iik,3) = albe(2,0,ns)
      iik = is + 1
      xl(iik,2) = albe(-1,1,ns) * xu(iik-1,1)
      xu(iik,1) = one / ( one - xu(iik-1,2) * xl(iik,2) )
      xu(iik,2) = albe(1,1,ns) - xu(iik-1,3) * xl(iik,2)
      xu(iik,3) = albe(2,1,ns)
      iik = is + 2
      xl(iik,1) = albe(-2,2,ns) * xu(iik-2,1)
      xl(iik,2) = ( albe(-1,2,ns) - xu(iik-2,2) * xl(iik,1) ) * xu(iik-1,1)
      xu(iik,1) = one / ( one - xu(iik-2,3) * xl(iik,1) - xu(iik-1,2) * xl(iik,2))
      xu(iik,2) = albe(1,2,ns) - xu(iik-1,3) * xl(iik,2)
      xu(iik,3) = albe(2,2,ns)
      do iik = is+3,ie-3
         xl(iik,1) = beto * xu(iik-2,1)
         xl(iik,2) = ( alpho - xu(iik-2,2) * xl(iik,1) ) * xu(iik-1,1)
         xu(iik,1) = one / ( one - xu(iik-2,3) * xl(iik,1) - xu(iik-1,2) * xl(iik,2) )
         xu(iik,2) = alpho - xu(iik-1,3) * xl(iik,2)
         xu(iik,3) = beto
      end do
      iik = ie - 2
      xl(iik,1) = albe(2,2,ne) * xu(iik-2,1)
      xl(iik,2) = ( albe(1,2,ne) - xu(iik-2,2) * xl(iik,1) ) * xu(iik-1,1)
      xu(iik,1) = one / ( one - xu(iik-2,3) * xl(iik,1) - xu(iik-1,2) * xl(iik,2) )
      xu(iik,2) = albe(-1,2,ne) - xu(iik-1,3) * xl(iik,2)
      xu(iik,3) = albe(-2,2,ne)
      iik = ie - 1
      xl(iik,1) = albe(2,1,ne) * xu(iik-2,1)
      xl(iik,2) = ( albe(1,1,ne) - xu(iik-2,2) * xl(iik,1) ) * xu(iik-1,1)
      xu(iik,1) = one / ( one - xu(iik-2,3) * xl(iik,1) - xu(iik-1,2) * xl(iik,2) )
      xu(iik,2) = albe(-1,1,ne) - xu(iik-1,3) * xl(iik,2)
      iik = ie
      xl(iik,1) = albe(2,0,ne) * xu(iik-2,1)
      xl(iik,2) = ( albe(1,0,ne) - xu(iik-2,2) * xl(iik,1) ) * xu(iik-1,1)
      xu(iik,1) = one / ( one - xu(iik-2,3) * xl(iik,1) - xu(iik-1,2) * xl(iik,2) )
      do iik = is,ie
         xu(iik,2:3) = xu(iik,2:3) * xu(iik,1)
      end do

   end subroutine penta

!===== SUBROUTINE FOR BOUNDARY FILTER COEFFICIENTS

   subroutine fcbcm(this,fltk,fltrbc)
      class(t_numerics), intent(inout) :: this
      real(kind=nr),intent(in) :: fltk,fltrbc
      real(kind=nr) :: alphz,betz,za,zb,zc
      real(kind=nr) :: aok, fctrk, resk

      aok = log(fltrbc)
      call this%fcint(fltk, half, alphz, betz, za, zb, zc)
      fctrk = one / &
            ( one + alphz * fltrbc + betz * fltrbc**two)

      this%albef(:,0,0) = (/ zero,          &
                        zero,          &
                        one,           &
                        alphz * fctrk, &
                        betz * fctrk   /)
      resk = (fltrbc-1) * &
            (za + zc + (fltrbc+1) * (zb+fltrbc*zc)) / aok
      this%fbc(:,0) = (/ za - 5 * resk / 3,   &
                    zb + 10 * resk / 21, &
                    zc - 5 * resk / 42,  &
                    5 * resk / 252,      &
                    -resk / 630          /) * fctrk

      this%albef(:,1,0) = (/ zero,                  &
                        alphz + betz * fltrbc, &
                        one,                   &
                        alphz,                 &
                        betz                   /)
      resk = (fltrbc-1) * &
            (zb + zc * (fltrbc+1)) / aok
      this%fbc(:,1) = (/ za + zb + zc + 1627 * resk / 1260, &
                    za + 10 * resk / 21,               &
                    zb - 5 * resk / 42,                &
                    zc + 5 * resk / 252,               &
                    -resk / 630                        /)

      this%albef(:,2,0) = (/ betz,  &
                        alphz, &
                        one,   &
                        alphz, &
                        betz   /)
      resk = zc * (fltrbc-1) / aok
      this%fbc(:,2) = (/ zb + zc + 1627 * resk / 1260, &
                    za - 5 * resk / 3,            &
                    za - 5 * resk / 42,           &
                    zb + 5 * resk / 252,          &
                    zc - resk / 630               /)

   end subroutine fcbcm

!===== SUBROUTINE FOR INTERIOR FILTER COEFFICIENTS

   subroutine fcint(this,fltk,fltr,alphz,betz,za,zb,zc)
      class(t_numerics), intent(inout) :: this
      real(kind=nr),intent(in) :: fltk,fltr
      real(kind=nr),intent(inout) :: alphz,betz,za,zb,zc
      real(kind=nr),dimension(3) :: cosf
      real(kind=nr) :: fctrk

      cosf(1) = cos(fltk  )
      cosf(2) = cos(2*fltk)
      cosf(3) = cos(3*fltk)
      fctrk = 1 / ( 30 + 5 * (7-16*fltr) * cosf(1) + &
                   2 * (1+8*fltr) * cosf(2) - 3 * cosf(3))
      alphz = fctrk * (20 * (2*fltr-1) - 30 * cosf(1) + &
                      12 * (2 * fltr - 1) * cosf(2) - 2 * cosf(3))
      betz = half * fctrk * (2 * (13 - 8 * fltr) + &
             (33 - 48 * fltr) * cosf(1) + 6 * cosf(2) - cosf(3))
      za = 60 * (1-fltr) * cos(half*fltk)**4 * fctrk
      zb = -2 * za / 5
      zc = za / 15

   end subroutine fcint

!===== SUBROUTINE FOR SUBDOMAIN-BOUNDARY COEFFICIENTS

   subroutine sbcco(this)
      class(t_numerics), intent(inout) :: this
      real(kind=nr),dimension(:,:),allocatable :: ax,bx,rx,sx
      integer(kind=ni) :: ntk, llk, iik, istart, iend

      do ntk = 0,1

         if ( ntk == 0 ) then
            llk = lmd
            istart = 1
            iend = 2 * ( llk + 1 )
            allocate( ax(iend,iend), bx(iend,iend), &
                      rx(iend,iend), sx(iend,iend) )
            ax(:,:) = 0
            bx(:,:) = 0
            ax(istart,istart:istart+2) = (/ one, alpha01, beta02 /)
            bx(istart,istart:istart+4) = (/ -( a01 + a02 + a03 + a04 ), &
                                               a01,                     &
                                               a02,                     &
                                               a03,                     &
                                               a04                      /)
            ax(istart+1,istart:istart+3) = (/ alpha10, one, alpha12, beta13 /)
            bx(istart+1,istart:istart+4) = (/ a10,                      &
                                           -( a10 + a12 + a13 + a14 ),  &
                                              a12,                      &
                                              a13,                      &
                                              a14                       /)
            do iik = istart+2,iend-2
               ax(iik,iik-2:iik+2) = (/ beta, alpha, one, alpha, beta /)
               bx(iik,iik-2:iik+2) = (/ -ab, -aa, zero, aa, ab /)
            end do
            ax(iend-1,iend:iend-3:-1) =  ax(istart+1,istart:istart+3)
            bx(iend-1,iend:iend-4:-1) = -bx(istart+1,istart:istart+4)
            ax(iend,iend:iend-2:-1)   =  ax(istart,istart:istart+2)
            bx(iend,iend:iend-4:-1)   = -bx(istart,istart:istart+4)
         end if

         if ( ntk == 1 ) then
            llk = lmf
            istart = 1
            iend   = 2 * ( llk + 1 )
            allocate( ax(iend,iend), bx(iend,iend), &
                      rx(iend,iend), sx(iend,iend) )
            ax(:,:)                   = 0
            bx(:,:)                   = 0
            ax(istart,istart:istart+2)        = this%albef(0:2,0,0)
            bx(istart,istart+(/1,2,3,4,5/))   = this%fbc(:,0)
            bx(istart,istart)                 = -sum( this%fbc(:,0) )
            ax(istart+1,istart:istart+3)      = this%albef(-1:2,1,0)
            bx(istart+1,istart+(/0,2,3,4,5/)) = this%fbc(:,1)
            bx(istart+1,istart+1)             = -sum( this%fbc(:,1) )
            ax(istart+2,istart:istart+4)      = this%albef(-2:2,2,0)
            bx(istart+2,istart+(/0,1,3,4,5/)) = this%fbc(:,2)
            bx(istart+2,istart+2)             = -sum( this%fbc(:,2) )
            do iik = istart+3,iend-3
               ax(iik,iik-2:iik+2)          = (/ this%betf, this%alphf, one, this%alphf, this%betf /)
               bx(iik,iik-3:iik+3)          = (/ this%fc, this%fb, this%fa,            &
                                                 -2 * ( this%fa + this%fb + this%fc ), &
                                                 this%fa, this%fb, this%fc             /)
            end do
            ax(iend-2,iend:iend-4:-1) = ax(istart+2,istart:istart+4)
            bx(iend-2,iend:iend-5:-1) = bx(istart+2,istart:istart+5)
            ax(iend-1,iend:iend-3:-1) = ax(istart+1,istart:istart+3)
            bx(iend-1,iend:iend-5:-1) = bx(istart+1,istart:istart+5)
            ax(iend,iend:iend-2:-1)   = ax(istart,istart:istart+2)
            bx(iend,iend:iend-5:-1)   = bx(istart,istart:istart+5)
         end if

         call mtrxi( ax(:,:), sx(:,:), istart, iend )

         rx(:,:) = ax(:,:)

         iik = iend / 2 - 1
         rx(iik,iik+2)   = 0
         rx(iik+1,iik+2) = 0
         rx(iik+1,iik+3) = 0
         
         iik = iend / 2 + 2
         rx(iik,iik-2)   = 0
         rx(iik-1,iik-2) = 0
         rx(iik-1,iik-3) = 0
         
         ax(:,:) = matmul( rx(:,:), matmul( sx(:,:), bx(:,:) ) )
         
         iik = iend / 2 + 1
         this%pbco(llk:0:-1,0,ntk) = ax(iik,istart:istart+llk)
         this%pbci(0:llk,0,ntk)    = ax(iik,istart+llk+1:iend)
         
         iik = iend / 2 + 2
         this%pbco(llk:0:-1,1,ntk) = ax(iik,istart:istart+llk)
         this%pbci(0:llk,1,ntk)    = ax(iik,istart+llk+1:iend)
         
         deallocate( ax, bx, rx, sx )
           
      end do
   end subroutine sbcco

!===== SUBROUTINE FOR HALO EXCHANGE

   subroutine mpigo_1d(this, rfield, lmx, ijks, nbck, mcdk, nbsizek, nt, n45, itag, lxi, let, m)
      class(t_numerics), intent(inout) :: this
      real(kind=nr),    intent(inout), dimension(0:lmx) :: rfield
      integer(kind=ni), intent(in)                   :: lmx
      integer(kind=ni), intent(in), dimension(3,3)   :: ijks
      integer(kind=ni), intent(in), dimension(3,0:1) :: nbck
      integer(kind=ni), intent(in), dimension(3,0:1) :: mcdk
      integer(kind=ni), intent(in), dimension(3)     :: nbsizek
      integer(kind=ni), intent(in)                   :: nt,n45,itag
      integer(kind=ni), intent(in)                   :: lxi, let, m
      integer(kind=ni) :: mpk, nnk, nzk, ipk, iqk, istart, iend
      integer(kind=ni) :: iik, iii, jjj, kkk, kpp, jkk, lll
      real(kind=nr)    :: ra0k, resk

      select case(nt)
      case(0)
         mpk = lmd
      case(1)
         mpk = lmf
      end select

      call p_null_req
      do nnk = 1,3

         if ( nt == 0 ) then
            select case(nnk)
            case(1)
               this%send => this%send01
               this%recv => this%recv01
            case(2)
               this%send => this%send02
               this%recv => this%recv02
            case(3)
               this%send => this%send03
               this%recv => this%recv03
            end select
         else
            select case(nnk)
            case(1)
               this%send => this%send11
               this%recv => this%recv11
            case(2)
               this%send => this%send12
               this%recv => this%recv12
            case(3)
               this%send => this%send13
               this%recv => this%recv13
            end select
         end if
         
         do ipk = 0,1
            iqk    = 1 - ipk
            istart = ipk * ijks(1,nnk)
            iend   = 1 - 2 * ipk

            select case(nbck(nnk,ipk))
            case(BC_INTER_STRAIGHT)
               ra0k = zero
               iik  = 1
            case(BC_INTER_SUBDOMAINS)
               ra0k = zero
               iik  = 0
            case(BC_PERIODIC)
               ra0k = n45
               iik  = 1
            end select

            if ( this%ndf(nnk,ipk,nt) == 1 ) then
               do kkk = 0,ijks(3,nnk)
                  kpp = kkk * ( ijks(2,nnk) + 1 )
                  do jjj = 0,ijks(2,nnk)
                     jkk  = kpp + jjj
                     lll  = indx3(istart, jjj, kkk, nnk,lxi,let)
                     resk = ra0k * rfield(lll)
                     do iii = 0,mpk
                        lll = indx3(istart+iend*(iii+iik), jjj, kkk, nnk,lxi,let)
                        this%sap(iii) = rfield(lll)
                     end do
                     this%send(jkk,0,ipk)    = sum( this%pbco(0:mpk,0,nt) * this%sap(0:mpk) ) - resk * this%pbcot(0,nt)
                     this%send(jkk,1,ipk)    = sum( this%pbco(0:mpk,1,nt) * this%sap(0:mpk) ) - resk * this%pbcot(1,nt)
                     this%send(jkk,nt+1,ipk) = this%send(jkk,nt+1,ipk) + nt * ( this%sap(0) - resk - this%send(jkk,nt+1,ipk) )
                  end do
               end do
               if ( nt == 0 ) then
                  call p_isend(this%send(:,:,ipk), mcdk(nnk,ipk), itag+iqk, 2*nbsizek(nnk))
                  call p_irecv(this%recv(:,:,ipk,m), mcdk(nnk,ipk), itag+ipk, 2*nbsizek(nnk))
               else
                  call p_isend(this%send(:,:,ipk), mcdk(nnk,ipk), itag+iqk, 3*nbsizek(nnk))
                  call p_irecv(this%recv(:,:,ipk,m), mcdk(nnk,ipk), itag+ipk, 3*nbsizek(nnk))
               end if
            end if
         end do
      end do
      call p_waitall

      if ( n45 == n45go ) then
         do nnk = 1,3

            if ( nt == 0 ) then
               select case(nnk)
               case(1)
                  this%recv => this%recv01
               case(2)
                  this%recv => this%recv02
               case(3)
                  this%recv => this%recv03
               end select
            else
               select case(nnk)
               case(1)
                  this%recv => this%recv11
               case(2)
                  this%recv => this%recv12
               case(3)
                  this%recv => this%recv13
               end select
            end if
            
            do ipk = 0,1
               istart = ipk * ijks(1,nnk)
               if ( nbck(nnk,ipk) == BC_PERIODIC ) then
                  do kkk = 0,ijks(3,nnk)
                     kpp = kkk * ( ijks(2,nnk) + 1 )
                     do jjj = 0,ijks(2,nnk)
                        jkk = kpp + jjj
                        lll = indx3(istart, jjj, kkk, nnk,lxi,let)
                        this%recv(jkk,0,ipk,m)    = this%recv(jkk,0,ipk,m) + rfield(lll) * this%pbcot(0,nt)
                        this%recv(jkk,1,ipk,m)    = this%recv(jkk,1,ipk,m) + rfield(lll) * this%pbcot(1,nt)
                        this%recv(jkk,nt+1,ipk,m) = this%recv(jkk,nt+1,ipk,m) + nt * rfield(lll)
                     end do
                  end do
               end if
            end do
         end do
      end if

   end subroutine mpigo_1d

   subroutine mpigo_2d(this, rfield, lmx, ijks, nbck, mcdk, nbsizek, nt, nrt, n45, itag, lxi, let, m)
      class(t_numerics), intent(inout) :: this
      real(kind=nr),    intent(inout), dimension(0:lmx,3) :: rfield
      integer(kind=ni), intent(in)                   :: lmx
      integer(kind=ni), intent(in), dimension(3,3)   :: ijks
      integer(kind=ni), intent(in), dimension(3,0:1) :: nbck
      integer(kind=ni), intent(in), dimension(3,0:1) :: mcdk
      integer(kind=ni), intent(in), dimension(3)     :: nbsizek
      integer(kind=ni), intent(in)                   :: nt,nrt,n45,itag
      integer(kind=ni), intent(in)                   :: lxi, let, m
      integer(kind=ni) :: mpk, nnk, nzk, ipk, iqk, istart, iend
      integer(kind=ni) :: iik, iii, jjj, kkk, kpp, jkk, lll
      real(kind=nr)    :: ra0k, resk

      select case(nt)
      case(0)
         mpk = lmd
      case(1)
         mpk = lmf
      end select

      call p_null_req
      do nnk = 1,3
         nzk = ( 1 - nrt ) * ( nnk - 1 ) + 1

         if ( nt == 0 ) then
            select case(nnk)
            case(1)
               this%send => this%send01
               this%recv => this%recv01
            case(2)
               this%send => this%send02
               this%recv => this%recv02
            case(3)
               this%send => this%send03
               this%recv => this%recv03
            end select
         else
            select case(nnk)
            case(1)
               this%send => this%send11
               this%recv => this%recv11
            case(2)
               this%send => this%send12
               this%recv => this%recv12
            case(3)
               this%send => this%send13
               this%recv => this%recv13
            end select
         end if
         
         do ipk = 0,1
            iqk    = 1 - ipk
            istart = ipk * ijks(1,nnk)
            iend   = 1 - 2 * ipk

            select case(nbck(nnk,ipk))
            case(BC_INTER_STRAIGHT)
               ra0k = zero
               iik  = 1
            case(BC_INTER_SUBDOMAINS)
               ra0k = zero
               iik  = 0
            case(BC_PERIODIC)
               ra0k = n45
               iik  = 1
            end select

            if ( this%ndf(nnk,ipk,nt) == 1 ) then
               do kkk = 0,ijks(3,nnk)
                  kpp = kkk * ( ijks(2,nnk) + 1 )
                  do jjj = 0,ijks(2,nnk)
                     jkk  = kpp + jjj
                     lll  = indx3(istart, jjj, kkk, nnk,lxi,let)
                     resk = ra0k * rfield(lll,nzk)
                     do iii = 0,mpk
                        lll = indx3(istart+iend*(iii+iik), jjj, kkk, nnk,lxi,let)
                        this%sap(iii) = rfield(lll,nzk)
                     end do
                     this%send(jkk,0,ipk)    = sum( this%pbco(0:mpk,0,nt) * this%sap(0:mpk) ) - resk * this%pbcot(0,nt)
                     this%send(jkk,1,ipk)    = sum( this%pbco(0:mpk,1,nt) * this%sap(0:mpk) ) - resk * this%pbcot(1,nt)
                     this%send(jkk,nt+1,ipk) = this%send(jkk,nt+1,ipk) + nt * ( this%sap(0) - resk - this%send(jkk,nt+1,ipk) )
                  end do
               end do
               if ( nt == 0 ) then
                  call p_isend(this%send(:,:,ipk), mcdk(nnk,ipk), itag+iqk, 2*nbsizek(nnk))
                  call p_irecv(this%recv(:,:,ipk,m), mcdk(nnk,ipk), itag+ipk, 2*nbsizek(nnk))
               else
                  call p_isend(this%send(:,:,ipk), mcdk(nnk,ipk), itag+iqk, 3*nbsizek(nnk))
                  call p_irecv(this%recv(:,:,ipk,m), mcdk(nnk,ipk), itag+ipk, 3*nbsizek(nnk))
               end if
            end if
         end do
      end do
      call p_waitall

      if ( n45 == n45go ) then
         do nnk = 1,3
            nzk = ( 1 - nrt ) * ( nnk - 1 ) + 1

            if ( nt == 0 ) then
               select case(nnk)
               case(1)
                  this%recv => this%recv01
               case(2)
                  this%recv => this%recv02
               case(3)
                  this%recv => this%recv03
               end select
            else
               select case(nnk)
               case(1)
                  this%recv => this%recv11
               case(2)
                  this%recv => this%recv12
               case(3)
                  this%recv => this%recv13
               end select
            end if
            
            do ipk = 0,1
               istart = ipk * ijks(1,nnk)
               if ( nbck(nnk,ipk) == BC_PERIODIC ) then
                  do kkk = 0,ijks(3,nnk)
                     kpp = kkk * ( ijks(2,nnk) + 1 )
                     do jjj = 0,ijks(2,nnk)
                        jkk = kpp + jjj
                        lll = indx3(istart, jjj, kkk, nnk,lxi,let)
                        this%recv(jkk,0,ipk,m)    = this%recv(jkk,0,ipk,m) + rfield(lll,nzk) * this%pbcot(0,nt)
                        this%recv(jkk,1,ipk,m)    = this%recv(jkk,1,ipk,m) + rfield(lll,nzk) * this%pbcot(1,nt)
                        this%recv(jkk,nt+1,ipk,m) = this%recv(jkk,nt+1,ipk,m) + nt * rfield(lll,nzk)
                     end do
                  end do
               end if
            end do
         end do
      end if

   end subroutine mpigo_2d

!===== SUBROUTINE FOR COMPACT FINITE DIFFERENTIATING

   subroutine deriv(this, rfield, lmx, lxik, letk, lzek, ijks, nn, nz, m, luse_acc)
      class(t_numerics), intent(inout) :: this
      real(kind=nr),    intent(inout), dimension(0:lmx,3) :: rfield
      integer(kind=ni), intent(in)                  :: lmx
      integer(kind=ni), intent(in)                  :: lxik, letk, lzek
      integer(kind=ni), intent(in), dimension(3,3)  :: ijks
      integer(kind=ni), intent(in)                  :: nn, nz, m
      logical,          intent(in), optional        :: luse_acc
      logical          :: lacc
      integer(kind=ni) :: ntk, nstart, nend, istart, iend, ustart, uend
      integer(kind=ni) :: kkk, jjj, iii, kpp, jkk, lll, ijj, lm

      if (present(luse_acc)) then
        lacc = luse_acc
      else
        lacc = .false.
      end if

      ntk    = 0
      nstart = this%ndf(nn,0,0)
      nend   = this%ndf(nn,1,0)

      select case(nn)
      case(1)
         ustart =  0
         uend   =  ustart + lxik
         this%recv   => this%recv01
         this%drva   => this%drva1
      case(2)
         ustart =  lxik + 1
         uend   =  ustart + letk
         this%recv   => this%recv02
         this%drva   => this%drva2
      case(3)
         ustart =  lxik + letk + 2
         uend   =  ustart + lzek
         this%recv   => this%recv03
         this%drva   => this%drva3
      end select

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            do iii = 0,ijks(1,nn)
               lll = indx3(iii, jjj, kkk, nn, lxik, letk)
               ijj = iii + jjj * (ijks(1,nn)+1) + kkk * ((ijks(1,nn)+1)*(ijks(2,nn)+1))
               this%sar(ijj) = rfield(lll,nz)
            end do
         end do
      end do
      !$ACC END PARALLEL

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            istart = jjj * (ijks(1,nn)+1) + kkk * ((ijks(1,nn)+1)*(ijks(2,nn)+1)) 
            iend   = istart + ijks(1,nn)

            kpp = kkk * ( ijks(2,nn) + 1 )
            jkk = kpp + jjj

            select case(nstart)
            case(0)
               this%sbr(istart)   = sum( (/a01,a02,a03,a04/) * ( this%sar(istart+(/1,2,3,4/)) - this%sar(istart)   ) )
               this%sbr(istart+1) = sum( (/a10,a12,a13,a14/) * ( this%sar(istart+(/0,2,3,4/)) - this%sar(istart+1) ) )
            case(1)
               this%sbr(istart)   = sum( this%pbci(0:lmd,0,ntk) * this%sar(istart:istart+lmd) ) + this%recv(jkk,0,0,m)
               this%sbr(istart+1) = sum( this%pbci(0:lmd,1,ntk) * this%sar(istart:istart+lmd) ) + this%recv(jkk,1,0,m)
            end select
         end do
      end do
      !$ACC END PARALLEL

      !$ACC PARALLEL LOOP COLLAPSE(2) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            istart = jjj * (ijks(1,nn)+1) + kkk * ((ijks(1,nn)+1)*(ijks(2,nn)+1))
            iend   = istart + ijks(1,nn)
            !$ACC LOOP VECTOR
            do iii = istart+2,iend-2
               this%sbr(iii) = aa * ( this%sar(iii+1) - this%sar(iii-1) ) + ab * ( this%sar(iii+2) - this%sar(iii-2) )
            end do
         end do
      end do
      !$ACC END PARALLEL

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            istart = jjj * (ijks(1,nn)+1) + kkk * ((ijks(1,nn)+1)*(ijks(2,nn)+1))
            iend   = istart + ijks(1,nn)

            kpp = kkk * ( ijks(2,nn) + 1 )
            jkk = kpp + jjj

            select case(nend)
            case(0)
               this%sbr(iend)   = sum( (/a01,a02,a03,a04/) * ( this%sar(iend)   - this%sar(iend-(/1,2,3,4/)) ) )
               this%sbr(iend-1) = sum( (/a10,a12,a13,a14/) * ( this%sar(iend-1) - this%sar(iend-(/0,2,3,4/)) ) )
            case(1)
               this%sbr(iend)   = -sum( this%pbci(0:lmd,0,ntk) * this%sar(iend:iend-lmd:-1) ) - this%recv(jkk,0,1,m)
               this%sbr(iend-1) = -sum( this%pbci(0:lmd,1,ntk) * this%sar(iend:iend-lmd:-1) ) - this%recv(jkk,1,1,m)
            end select
         end do
      end do
      !$ACC END PARALLEL

      ! inner loop carry dependencies
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            istart = jjj * (ijks(1,nn)+1) + kkk * ((ijks(1,nn)+1)*(ijks(2,nn)+1))  
            iend   = istart + ijks(1,nn)

            this%sar(istart)   = this%sbr(istart)
            this%sar(istart+1) = this%sbr(istart+1) - this%xl(ustart+1,2) * this%sar(istart)
            !$ACC LOOP SEQ
            do iii = istart+2,iend
               this%sar(iii) = this%sbr(iii) - this%xl(iii-istart+ustart,1) * this%sar(iii-2) - &
                                               this%xl(iii-istart+ustart,2) * this%sar(iii-1)
            end do
         end do
      end do
      !$ACC END PARALLEL

      ! inner loop carry dependencies
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            istart = jjj * (ijks(1,nn)+1) + kkk * ((ijks(1,nn)+1)*(ijks(2,nn)+1))  
            iend   = istart + ijks(1,nn)

            kpp = kkk * ( ijks(2,nn) + 1 )
            jkk = kpp + jjj

            this%sbr(iend)   = this%xu(uend,1)   * this%sar(iend)
            this%sbr(iend-1) = this%xu(uend-1,1) * this%sar(iend-1) - this%xu(uend-1,2) * this%sbr(iend)
            !$ACC LOOP SEQ
            do iii = iend-2,istart,-1
               this%sbr(iii) = this%xu(iii-iend+uend,1) * this%sar(iii)   - & 
                               this%xu(iii-iend+uend,2) * this%sbr(iii+1) - &
                               this%xu(iii-iend+uend,3) * this%sbr(iii+2)
            end do
            this%drva(jkk,m,0) = this%sbr(istart)
            this%drva(jkk,m,1) = this%sbr(iend)
         end do
      end do
      !$ACC END PARALLEL

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            do iii = 0,ijks(1,nn)
               lll = indx3(iii, jjj, kkk, nn, lxik, letk)
               ijj = iii + jjj * (ijks(1,nn)+1) + kkk * ((ijks(1,nn)+1)*(ijks(2,nn)+1))
               rfield(lll,nn) = this%sbr(ijj)
            end do
         end do
      end do
      !$ACC END PARALLEL

   end subroutine deriv

!===== SUBROUTINE FOR COMPACT FILTERING

   subroutine filte(this, rfield, lmx, lxik, letk, lzek, ijks, inn, m, luse_acc)
      class(t_numerics), intent(inout) :: this
      real(kind=nr),    intent(inout), dimension(0:lmx) :: rfield
      integer(kind=ni), intent(in)                  :: lmx
      integer(kind=ni), intent(in)                  :: lxik, letk, lzek
      integer(kind=ni), intent(in), dimension(3,3)  :: ijks
      integer(kind=ni), intent(in)                  :: inn, m
      logical, intent(in), optional                 :: luse_acc
      logical          :: lacc
      integer(kind=ni) :: nstart, nend, istart, iend, ntk, nn, ustart, uend
      integer(kind=ni) :: kkk, jjj, iii, lll, kpp, jkk, ijj
      real(kind=nr)    :: resk, ra2k

      if (present(luse_acc)) then
        lacc = luse_acc
      else
        lacc = .false.
      end if

      nn = nnf(inn)

      ntk    = 1
      nstart = this%ndf(nn,0,1)
      nend   = this%ndf(nn,1,1)

      select case(nn)
      case(1)
         ustart =  0
         uend   =  ustart + lxik
         this%recv   => this%recv11
      case(2)
         ustart =  lxik + 1
         uend   =  ustart + letk
         this%recv   => this%recv12
      case(3)
         ustart =  lxik + letk + 2
         uend   =  ustart + lzek
         this%recv   => this%recv13
      end select

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            do iii = 0,ijks(1,nn)
               lll = indx3(iii, jjj, kkk, nn, lxik, letk)
               ijj = iii + jjj * (ijks(1,nn)+1) + kkk * ((ijks(1,nn)+1)*(ijks(2,nn)+1))
               this%sar(ijj) = rfield(lll)
            end do
         end do
      end do
      !$ACC END PARALLEL

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            istart = jjj * (ijks(1,nn)+1) + kkk * ((ijks(1,nn)+1)*(ijks(2,nn)+1)) 
            iend   = istart + ijks(1,nn)
            kpp = kkk * ( ijks(2,nn) + 1 )
            jkk = kpp + jjj

            select case(nstart)
            case(0)
               this%sbr(istart)   = sum( this%fbc(:,0) * ( this%sar(istart+(/1,2,3,4,5/)) - this%sar(istart)   ) )
               this%sbr(istart+1) = sum( this%fbc(:,1) * ( this%sar(istart+(/0,2,3,4,5/)) - this%sar(istart+1) ) )
               this%sbr(istart+2) = sum( this%fbc(:,2) * ( this%sar(istart+(/0,1,3,4,5/)) - this%sar(istart+2) ) )
            case(1)
               ra2k              = this%sar(istart+2) + this%sar(istart+2)
               this%sbr(istart)   = sum( this%pbci(0:lmf,0,ntk) * this%sar(istart:istart+lmf) ) + this%recv(jkk,0,0,m)
               this%sbr(istart+1) = sum( this%pbci(0:lmf,1,ntk) * this%sar(istart:istart+lmf) ) + this%recv(jkk,1,0,m)
               this%sbr(istart+2) = this%fa * ( this%sar(istart+1)  + this%sar(istart+3) - ra2k ) + &
                                    this%fb * ( this%sar(istart)    + this%sar(istart+4) - ra2k ) + &
                                    this%fc * ( this%recv(jkk,2,0,m) + this%sar(istart+5) - ra2k )
            end select
         end do
      end do
      !$ACC END PARALLEL

      !$ACC PARALLEL LOOP COLLAPSE(2) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            istart = jjj * (ijks(1,nn)+1) + kkk * ((ijks(1,nn)+1)*(ijks(2,nn)+1))
            iend   = istart + ijks(1,nn)
            !$ACC LOOP VECTOR
            do iii = istart+3,iend-3
               resk         = this%sar(iii) + this%sar(iii)
               this%sbr(iii) = this%fa * ( this%sar(iii-1) + this%sar(iii+1) - resk ) + &
                               this%fb * ( this%sar(iii-2) + this%sar(iii+2) - resk ) + &
                               this%fc * ( this%sar(iii-3) + this%sar(iii+3) - resk )
            end do
         end do
      end do
      !$ACC END PARALLEL

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            kpp = kkk * ( ijks(2,nn) + 1 )
            jkk = kpp + jjj
            
            select case(nend)
            case(0)
               this%sbr(iend)   = sum( this%fbc(:,0) * ( this%sar(iend-(/1,2,3,4,5/)) - this%sar(iend)   ) )
               this%sbr(iend-1) = sum( this%fbc(:,1) * ( this%sar(iend-(/0,2,3,4,5/)) - this%sar(iend-1) ) )
               this%sbr(iend-2) = sum( this%fbc(:,2) * ( this%sar(iend-(/0,1,3,4,5/)) - this%sar(iend-2) ) )
            case(1)
               ra2k             = this%sar(iend-2) + this%sar(iend-2)
               this%sbr(iend)   = sum( this%pbci(0:lmf,0,ntk) * this%sar(iend:iend-lmf:-1) ) + this%recv(jkk,0,1,m)
               this%sbr(iend-1) = sum( this%pbci(0:lmf,1,ntk) * this%sar(iend:iend-lmf:-1) ) + this%recv(jkk,1,1,m)
               this%sbr(iend-2) = this%fa * ( this%sar(iend-3) + this%sar(iend-1)    - ra2k ) + &
                                 this%fb * ( this%sar(iend-4) + this%sar(iend)      - ra2k ) + &
                                 this%fc * ( this%sar(iend-5) + this%recv(jkk,2,1,m) - ra2k )
            end select
         end do
      end do
      !$ACC END PARALLEL

      ! inner loop carry dependencies
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            istart = jjj * (ijks(1,nn)+1) + kkk * ((ijks(1,nn)+1)*(ijks(2,nn)+1))
            iend   = istart + ijks(1,nn)

            this%sar(istart)   = this%sbr(istart)
            this%sar(istart+1) = this%sbr(istart+1) - this%yl(ustart+1,2) * this%sar(istart)
            !$ACC LOOP SEQ
            do iii = istart+2,iend
               this%sar(iii) = this%sbr(iii) - this%yl(iii-istart+ustart,1) * this%sar(iii-2) - &
                                               this%yl(iii-istart+ustart,2) * this%sar(iii-1)
            end do
         end do
      end do
      !$ACC END PARALLEL

      ! inner loop carry dependencies
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            istart = jjj * (ijks(1,nn)+1) + kkk * ((ijks(1,nn)+1)*(ijks(2,nn)+1))
            iend   = istart + ijks(1,nn)

            this%sbr(iend)   = this%yu(uend,1)   * this%sar(iend)
            this%sbr(iend-1) = this%yu(uend-1,1) * this%sar(iend-1) - this%yu(uend-1,2) * this%sbr(iend)
            !$ACC LOOP SEQ
            do iii = iend-2,istart,-1
               this%sbr(iii) = this%yu(iii-iend+uend,1) * this%sar(iii)   - &
                               this%yu(iii-iend+uend,2) * this%sbr(iii+1) - &
                               this%yu(iii-iend+uend,3) * this%sbr(iii+2)
            end do
         end do
      end do
      !$ACC END PARALLEL

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) IF (lacc)
      do kkk = 0,ijks(3,nn)
         do jjj = 0,ijks(2,nn)
            do iii = 0,ijks(1,nn)
               lll = indx3(iii, jjj, kkk, nn, lxik, letk)
               ijj = iii + jjj * (ijks(1,nn)+1) + kkk * ((ijks(1,nn)+1)*(ijks(2,nn)+1))
               rfield(lll) = rfield(lll) + this%sbr(ijj)
            end do
         end do
      end do
      !$ACC END PARALLEL

   end subroutine filte

!=====

end module mo_numerics

!*****
