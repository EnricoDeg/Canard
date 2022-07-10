!*****
!***** 3D SOLVER MODULE
!*****

module mo_numerics
   use mo_parameters, ONLY : pi, half, zero, one, alpha01, alpha10,   &
                             alpha, beta, beta02, beta13, alpha12,    &
                             a01, a02, a03, a04, a10, a12, a13, a14, &
                             ab, aa, two, n45go
   use mo_kind,       ONLY : nr, ni
   use mo_mpi,        ONLY : p_null_req, p_isend, p_irecv, p_waitall, &
                             p_send, myid
   use mo_utils,      ONLY : indx3, mtrxi
   implicit none
   private

   integer(kind=ni), parameter :: lmd=11,lmf=11,lmp=max(lmd,lmf)
   real(kind=nr)               :: alphf, betf
   real(kind=nr), dimension(0:lmp,0:1,0:1) :: pbci,pbco
   real(kind=nr), dimension(0:1,0:1) :: pbcot
   real(kind=nr) :: fa,fb,fc
   real(kind=nr), dimension(-2:2,0:2,0:1) :: albef
   real(kind=nr), dimension(0:lmp) :: sap
   real(kind=nr), dimension(0:4,0:2) :: fbc
   integer(kind=ni), dimension(3,0:1,0:1) :: ndf
   integer(kind=ni),dimension(3) :: nnf

   real(kind=nr), dimension(:,:), allocatable :: xu,yu
   real(kind=nr), dimension(:,:), allocatable :: xl,yl
   character(16) :: ccinput
   real(kind=nr) :: fltk,fltrbc

   real(kind=nr), dimension(:,:,:), allocatable, target :: send01,send02,send03
   real(kind=nr), dimension(:,:,:), allocatable, target :: send11,send12,send13
   real(kind=nr), dimension(:,:,:), allocatable, target :: recv01,recv02,recv03
   real(kind=nr), dimension(:,:,:), allocatable, target :: recv11,recv12,recv13
   real(kind=nr), dimension(:,:,:), pointer :: send,recv
   integer(kind=ni), dimension(:),  allocatable :: li
   real(kind=nr), dimension(:),     allocatable :: sa,sb

   real(kind=nr), public, dimension(:,:,:), allocatable, target :: drva1,drva2,drva3
   real(kind=nr), public, dimension(:,:,:), pointer :: drva

   public :: allocate_numerics, read_input_numerics, init_extracoeff_bounds
   public :: init_penta, mpigo, deriv, filte

   INTERFACE mpigo
      MODULE PROCEDURE mpigo_1d
      MODULE PROCEDURE mpigo_2d
   END INTERFACE mpigo

   INTERFACE deriv
      MODULE PROCEDURE deriv_nooverwrite
      MODULE PROCEDURE deriv_overwrite
   END INTERFACE deriv

   contains

   SUBROUTINE allocate_numerics(limk, nbsizek)
      integer(kind=ni),intent(in) :: limk
      integer(kind=ni),intent(in), dimension(3) :: nbsizek
      integer(kind=ni) :: iik, jjk, kkk

      iik=nbsizek(1)-1
      jjk=nbsizek(2)-1
      kkk=nbsizek(3)-1

      allocate(xu(0:limk,3),yu(0:limk,3),xl(0:limk,2),yl(0:limk,2))
      allocate(li(0:limk),  sa(0:limk),  sb(0:limk))
      allocate(send01(0:iik,0:1,0:1),send02(0:jjk,0:1,0:1),send03(0:kkk,0:1,0:1))
      allocate(recv01(0:iik,0:1,0:1),recv02(0:jjk,0:1,0:1),recv03(0:kkk,0:1,0:1))
      allocate(send11(0:iik,0:2,0:1),send12(0:jjk,0:2,0:1),send13(0:kkk,0:2,0:1))
      allocate(recv11(0:iik,0:2,0:1),recv12(0:jjk,0:2,0:1),recv13(0:kkk,0:2,0:1))
      allocate(drva1(0:iik,5,0:1), drva2(0:jjk,5,0:1), drva3(0:kkk,5,0:1))

   END SUBROUTINE allocate_numerics

   SUBROUTINE read_input_numerics

      open(9,file='input.numerics',status='old')
      read(9,*) ccinput,fltk
      read(9,*) ccinput,fltrbc
      read(9,*) ccinput,nnf(:)
      close(9)
      fltk=pi*fltk

   END SUBROUTINE read_input_numerics

!===== EXTRA COEFFICIENTS FOR DOMAIN BOUNDARIES

   SUBROUTINE init_extracoeff_bounds
      integer(kind=ni) :: ntk, jjk, iik

      call fcbcm(fltk,fltrbc)
      call fcint(fltk,half,alphf,betf,fa,fb,fc)
      albef(:,0,1) = (/ zero, zero,  one, alphf, betf /)
      albef(:,1,1) = (/ zero, alphf, one, alphf, betf /)
      albef(:,2,1) = (/ betf, alphf, one, alphf, betf /)

      pbco(:,:,:) = zero
      pbci(:,:,:) = zero
      call sbcco
      do ntk = 0,1
         do jjk = 0,1
            iik = lmd + ntk * ( lmf - lmd )
            pbcot(jjk,ntk) = sum( pbco(0:iik,jjk,ntk) )
         end do
      end do

   END SUBROUTINE init_extracoeff_bounds

   SUBROUTINE init_penta(lxik, letk, lzek, nbck, lim)
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
            case(10,20,25,30)
               ndf(nnk,ipk,0) = 0
               ndf(nnk,ipk,1) = 0
            case(35,40,45)
               ndf(nnk,ipk,0) = 1
               ndf(nnk,ipk,1) = 1
            end select
         end do
         nstart = ndf(nnk,0,0)
         nend   = ndf(nnk,1,0)
         call penta(xu(:,:), xl(:,:), istart, iend, nstart, nend, 0, lim)
         nstart = ndf(nnk,0,1)
         nend   = ndf(nnk,1,1)
         call penta(yu(:,:), yl(:,:), istart, iend, nstart, nend, 1, lim)
      end do
  
   END SUBROUTINE init_penta

!===== SUBROUTINE FOR CHOLESKY DECOMPOSITION OF PENTADIAGONAL MATRICES

   subroutine penta(xu,xl,is,ie,ns,ne,nt,lim)

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
         albe  = albef
         alpho = alphf
         beto  = betf
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

   subroutine fcbcm(fltk,fltrbc)
 
      real(kind=nr),intent(in) :: fltk,fltrbc
      real(kind=nr) :: alphz,betz,za,zb,zc
      real(kind=nr) :: aok, fctrk, resk

      aok = log(fltrbc)
      call fcint(fltk, half, alphz, betz, za, zb, zc)
      fctrk = one / &
            ( one + alphz * fltrbc + betz * fltrbc**two)

      albef(:,0,0) = (/ zero,          &
                        zero,          &
                        one,           &
                        alphz * fctrk, &
                        betz * fctrk   /)
      resk = (fltrbc-1) * &
            (za + zc + (fltrbc+1) * (zb+fltrbc*zc)) / aok
      fbc(:,0) = (/ za - 5 * resk / 3,   &
                    zb + 10 * resk / 21, &
                    zc - 5 * resk / 42,  &
                    5 * resk / 252,      &
                    -resk / 630          /) * fctrk

      albef(:,1,0) = (/ zero,                  &
                        alphz + betz * fltrbc, &
                        one,                   &
                        alphz,                 &
                        betz                   /)
      resk = (fltrbc-1) * &
            (zb + zc * (fltrbc+1)) / aok
      fbc(:,1) = (/ za + zb + zc + 1627 * resk / 1260, &
                    za + 10 * resk / 21,               &
                    zb - 5 * resk / 42,                &
                    zc + 5 * resk / 252,               &
                    -resk / 630                        /)

      albef(:,2,0) = (/ betz,  &
                        alphz, &
                        one,   &
                        alphz, &
                        betz   /)
      resk = zc * (fltrbc-1) / aok
      fbc(:,2) = (/ zb + zc + 1627 * resk / 1260, &
                    za - 5 * resk / 3,            &
                    za - 5 * resk / 42,           &
                    zb + 5 * resk / 252,          &
                    zc - resk / 630               /)

   end subroutine fcbcm

!===== SUBROUTINE FOR INTERIOR FILTER COEFFICIENTS

   subroutine fcint(fltk,fltr,alphz,betz,za,zb,zc)
 
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

   subroutine sbcco

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
            ax(istart,istart:istart+2)        = albef(0:2,0,0)
            bx(istart,istart+(/1,2,3,4,5/))   = fbc(:,0)
            bx(istart,istart)                 = -sum( fbc(:,0) )
            ax(istart+1,istart:istart+3)      = albef(-1:2,1,0)
            bx(istart+1,istart+(/0,2,3,4,5/)) = fbc(:,1)
            bx(istart+1,istart+1)             = -sum( fbc(:,1) )
            ax(istart+2,istart:istart+4)      = albef(-2:2,2,0)
            bx(istart+2,istart+(/0,1,3,4,5/)) = fbc(:,2)
            bx(istart+2,istart+2)             = -sum( fbc(:,2) )
            do iik = istart+3,iend-3
               ax(iik,iik-2:iik+2)          = (/ betf, alphf, one, alphf, betf /)
               bx(iik,iik-3:iik+3)          = (/ fc, fb, fa,            &
                                           -2 * ( fa + fb + fc ), &
                                           fa, fb, fc             /)
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
         pbco(llk:0:-1,0,ntk) = ax(iik,istart:istart+llk)
         pbci(0:llk,0,ntk)    = ax(iik,istart+llk+1:iend)
         
         iik = iend / 2 + 2
         pbco(llk:0:-1,1,ntk) = ax(iik,istart:istart+llk)
         pbci(0:llk,1,ntk)    = ax(iik,istart+llk+1:iend)
         
         deallocate( ax, bx, rx, sx )
           
      end do
   end subroutine sbcco

!===== SUBROUTINE FOR HALO EXCHANGE

   subroutine mpigo_1d(rfield, lmx, ijks, nbck, mcdk, nbsizek, nt, n45, itag, lxi, let)
      real(kind=nr),    intent(inout), dimension(0:lmx) :: rfield
      integer(kind=ni), intent(in)                   :: lmx
      integer(kind=ni), intent(in), dimension(3,3)   :: ijks
      integer(kind=ni), intent(in), dimension(3,0:1) :: nbck
      integer(kind=ni), intent(in), dimension(3,0:1) :: mcdk
      integer(kind=ni), intent(in), dimension(3)     :: nbsizek
      integer(kind=ni), intent(in)                   :: nt,n45,itag
      integer(kind=ni), intent(in)                   :: lxi, let
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
               send => send01
               recv => recv01
            case(2)
               send => send02
               recv => recv02
            case(3)
               send => send03
               recv => recv03
            end select
         else
            select case(nnk)
            case(1)
               send => send11
               recv => recv11
            case(2)
               send => send12
               recv => recv12
            case(3)
               send => send13
               recv => recv13
            end select
         end if
         
         do ipk = 0,1
            iqk    = 1 - ipk
            istart = ipk * ijks(1,nnk)
            iend   = 1 - 2 * ipk

            select case(nbck(nnk,ipk))
            case(35)
               ra0k = zero
               iik  = 1
            case(40)
               ra0k = zero
               iik  = 0
            case(45)
               ra0k = n45
               iik  = 1
            end select

            if ( ndf(nnk,ipk,nt) == 1 ) then
               do kkk = 0,ijks(3,nnk)
                  kpp = kkk * ( ijks(2,nnk) + 1 )
                  do jjj = 0,ijks(2,nnk)
                     jkk  = kpp + jjj
                     lll  = indx3(istart, jjj, kkk, nnk,lxi,let)
                     resk = ra0k * rfield(lll)
                     do iii = 0,mpk
                        lll = indx3(istart+iend*(iii+iik), jjj, kkk, nnk,lxi,let)
                        sap(iii) = rfield(lll)
                     end do
                     send(jkk,0,ipk)    = sum( pbco(0:mpk,0,nt) * sap(0:mpk) ) - resk * pbcot(0,nt)
                     send(jkk,1,ipk)    = sum( pbco(0:mpk,1,nt) * sap(0:mpk) ) - resk * pbcot(1,nt)
                     send(jkk,nt+1,ipk) = send(jkk,nt+1,ipk) + nt * ( sap(0) - resk - send(jkk,nt+1,ipk) )
                  end do
               end do
               if ( nt == 0 ) then
                  call p_isend(send(:,:,ipk), mcdk(nnk,ipk), itag+iqk, 2*nbsizek(nnk))
                  call p_irecv(recv(:,:,ipk), mcdk(nnk,ipk), itag+ipk, 2*nbsizek(nnk))
               else
                  call p_isend(send(:,:,ipk), mcdk(nnk,ipk), itag+iqk, 3*nbsizek(nnk))
                  call p_irecv(recv(:,:,ipk), mcdk(nnk,ipk), itag+ipk, 3*nbsizek(nnk))
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
                  recv => recv01
               case(2)
                  recv => recv02
               case(3)
                  recv => recv03
               end select
            else
               select case(nnk)
               case(1)
                  recv => recv11
               case(2)
                  recv => recv12
               case(3)
                  recv => recv13
               end select
            end if
            
            do ipk = 0,1
               istart = ipk * ijks(1,nnk)
               if ( nbck(nnk,ipk) == 45 ) then
                  do kkk = 0,ijks(3,nnk)
                     kpp = kkk * ( ijks(2,nnk) + 1 )
                     do jjj = 0,ijks(2,nnk)
                        jkk = kpp + jjj
                        lll = indx3(istart, jjj, kkk, nnk,lxi,let)
                        recv(jkk,0,ipk)    = recv(jkk,0,ipk) + rfield(lll) * pbcot(0,nt)
                        recv(jkk,1,ipk)    = recv(jkk,1,ipk) + rfield(lll) * pbcot(1,nt)
                        recv(jkk,nt+1,ipk) = recv(jkk,nt+1,ipk) + nt * rfield(lll)
                     end do
                  end do
               end if
            end do
         end do
      end if

   end subroutine mpigo_1d

   subroutine mpigo_2d(rfield, lmx, ijks, nbck, mcdk, nbsizek, nt, nrt, n45, itag, lxi, let)
      real(kind=nr),    intent(inout), dimension(0:lmx,3) :: rfield
      integer(kind=ni), intent(in)                   :: lmx
      integer(kind=ni), intent(in), dimension(3,3)   :: ijks
      integer(kind=ni), intent(in), dimension(3,0:1) :: nbck
      integer(kind=ni), intent(in), dimension(3,0:1) :: mcdk
      integer(kind=ni), intent(in), dimension(3)     :: nbsizek
      integer(kind=ni), intent(in)                   :: nt,nrt,n45,itag
      integer(kind=ni), intent(in)                   :: lxi, let
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
               send => send01
               recv => recv01
            case(2)
               send => send02
               recv => recv02
            case(3)
               send => send03
               recv => recv03
            end select
         else
            select case(nnk)
            case(1)
               send => send11
               recv => recv11
            case(2)
               send => send12
               recv => recv12
            case(3)
               send => send13
               recv => recv13
            end select
         end if
         
         do ipk = 0,1
            iqk    = 1 - ipk
            istart = ipk * ijks(1,nnk)
            iend   = 1 - 2 * ipk

            select case(nbck(nnk,ipk))
            case(35)
               ra0k = zero
               iik  = 1
            case(40)
               ra0k = zero
               iik  = 0
            case(45)
               ra0k = n45
               iik  = 1
            end select

            if ( ndf(nnk,ipk,nt) == 1 ) then
               do kkk = 0,ijks(3,nnk)
                  kpp = kkk * ( ijks(2,nnk) + 1 )
                  do jjj = 0,ijks(2,nnk)
                     jkk  = kpp + jjj
                     lll  = indx3(istart, jjj, kkk, nnk,lxi,let)
                     resk = ra0k * rfield(lll,nzk)
                     do iii = 0,mpk
                        lll = indx3(istart+iend*(iii+iik), jjj, kkk, nnk,lxi,let)
                        sap(iii) = rfield(lll,nzk)
                     end do
                     send(jkk,0,ipk)    = sum( pbco(0:mpk,0,nt) * sap(0:mpk) ) - resk * pbcot(0,nt)
                     send(jkk,1,ipk)    = sum( pbco(0:mpk,1,nt) * sap(0:mpk) ) - resk * pbcot(1,nt)
                     send(jkk,nt+1,ipk) = send(jkk,nt+1,ipk) + nt * ( sap(0) - resk - send(jkk,nt+1,ipk) )
                  end do
               end do
               if ( nt == 0 ) then
                  call p_isend(send(:,:,ipk), mcdk(nnk,ipk), itag+iqk, 2*nbsizek(nnk))
                  call p_irecv(recv(:,:,ipk), mcdk(nnk,ipk), itag+ipk, 2*nbsizek(nnk))
               else
                  call p_isend(send(:,:,ipk), mcdk(nnk,ipk), itag+iqk, 3*nbsizek(nnk))
                  call p_irecv(recv(:,:,ipk), mcdk(nnk,ipk), itag+ipk, 3*nbsizek(nnk))
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
                  recv => recv01
               case(2)
                  recv => recv02
               case(3)
                  recv => recv03
               end select
            else
               select case(nnk)
               case(1)
                  recv => recv11
               case(2)
                  recv => recv12
               case(3)
                  recv => recv13
               end select
            end if
            
            do ipk = 0,1
               istart = ipk * ijks(1,nnk)
               if ( nbck(nnk,ipk) == 45 ) then
                  do kkk = 0,ijks(3,nnk)
                     kpp = kkk * ( ijks(2,nnk) + 1 )
                     do jjj = 0,ijks(2,nnk)
                        jkk = kpp + jjj
                        lll = indx3(istart, jjj, kkk, nnk,lxi,let)
                        recv(jkk,0,ipk)    = recv(jkk,0,ipk) + rfield(lll,nzk) * pbcot(0,nt)
                        recv(jkk,1,ipk)    = recv(jkk,1,ipk) + rfield(lll,nzk) * pbcot(1,nt)
                        recv(jkk,nt+1,ipk) = recv(jkk,nt+1,ipk) + nt * rfield(lll,nzk)
                     end do
                  end do
               end if
            end do
         end do
      end if

   end subroutine mpigo_2d

!===== SUBROUTINE FOR COMPACT FINITE DIFFERENTIATING

   subroutine deriv_nooverwrite(rfieldin, rfieldout, lmx, lxik, letk, lzek, ijks, nn, m)
      real(kind=nr),    intent(in),  dimension(0:lmx) :: rfieldin
      real(kind=nr),    intent(out), dimension(0:lmx) :: rfieldout
      integer(kind=ni), intent(in)                  :: lmx
      integer(kind=ni), intent(in)                  :: lxik, letk, lzek
      integer(kind=ni), intent(in), dimension(3,3)  :: ijks
      integer(kind=ni), intent(in)                  :: nn, m
      integer(kind=ni) :: ntk, nstart, nend, istart, iend
      integer(kind=ni) :: kkk, jjj, iii, kpp, jkk, lll

      ntk    = 0
      nstart = ndf(nn,0,0)
      nend   = ndf(nn,1,0)

      select case(nn)
      case(1)
         istart =  0
         iend   =  istart + lxik
         recv   => recv01
         drva   => drva1
      case(2)
         istart =  lxik + 1
         iend   =  istart + letk
         recv   => recv02
         drva   => drva2
      case(3)
         istart =  lxik + letk + 2
         iend   =  istart + lzek
         recv   => recv03
         drva   => drva3
      end select

      do kkk = 0,ijks(3,nn)
         kpp = kkk * ( ijks(2,nn) + 1 )
         do jjj=  0,ijks(2,nn)
            jkk = kpp + jjj
            do iii = istart,iend
               lll = indx3(iii-istart, jjj, kkk, nn, lxik, letk)
               li(iii) = lll
               sa(iii) = rfieldin(lll)
            end do

            select case(nstart)
            case(0)
               sb(istart)   = sum( (/a01,a02,a03,a04/) * ( sa(istart+(/1,2,3,4/)) - sa(istart)   ) )
               sb(istart+1) = sum( (/a10,a12,a13,a14/) * ( sa(istart+(/0,2,3,4/)) - sa(istart+1) ) )
            case(1)
               sb(istart)   = sum( pbci(0:lmd,0,ntk) * sa(istart:istart+lmd) ) + recv(jkk,0,0)
               sb(istart+1) = sum( pbci(0:lmd,1,ntk) * sa(istart:istart+lmd) ) + recv(jkk,1,0)
            end select

            do iii = istart+2,iend-2
               sb(iii) = aa * ( sa(iii+1) - sa(iii-1) ) + ab * ( sa(iii+2) - sa(iii-2) )
            end do

            select case(nend)
            case(0)
               sb(iend)   = sum( (/a01,a02,a03,a04/) * ( sa(iend)   - sa(iend-(/1,2,3,4/)) ) )
               sb(iend-1) = sum( (/a10,a12,a13,a14/) * ( sa(iend-1) - sa(iend-(/0,2,3,4/)) ) )
            case(1)
               sb(iend)   = -sum( pbci(0:lmd,0,ntk) * sa(iend:iend-lmd:-1) ) - recv(jkk,0,1)
               sb(iend-1) = -sum( pbci(0:lmd,1,ntk) * sa(iend:iend-lmd:-1) ) - recv(jkk,1,1)
            end select

            sa(istart)   = sb(istart)
            sa(istart+1) = sb(istart+1) - xl(istart+1,2) * sa(istart)
            do iii = istart+2,iend
               sa(iii) = sb(iii) - xl(iii,1) * sa(iii-2) - xl(iii,2) * sa(iii-1)
            end do

            sb(iend)   = xu(iend,1)   * sa(iend)
            sb(iend-1) = xu(iend-1,1) * sa(iend-1) - xu(iend-1,2) * sb(iend)
            do iii = iend-2,istart,-1
               sb(iii) = xu(iii,1) * sa(iii) - xu(iii,2) * sb(iii+1) - xu(iii,3) * sb(iii+2)
            end do

            do iii = istart,iend
               lll = li(iii)
               rfieldout(lll) = sb(iii)
            end do

            drva(jkk,m,0) = sb(istart)
            drva(jkk,m,1) = sb(iend)

         end do
      end do

   end subroutine deriv_nooverwrite

   subroutine deriv_overwrite(rfield, lmx, lxik, letk, lzek, ijks, nn, nz, m)
      real(kind=nr),    intent(inout), dimension(0:lmx,3) :: rfield
      integer(kind=ni), intent(in)                  :: lmx
      integer(kind=ni), intent(in)                  :: lxik, letk, lzek
      integer(kind=ni), intent(in), dimension(3,3)  :: ijks
      integer(kind=ni), intent(in)                  :: nn, nz, m
      integer(kind=ni) :: ntk, nstart, nend, istart, iend
      integer(kind=ni) :: kkk, jjj, iii, kpp, jkk, lll

      ntk    = 0
      nstart = ndf(nn,0,0)
      nend   = ndf(nn,1,0)

      select case(nn)
      case(1)
         istart =  0
         iend   =  istart + lxik
         recv   => recv01
         drva   => drva1
      case(2)
         istart =  lxik + 1
         iend   =  istart + letk
         recv   => recv02
         drva   => drva2
      case(3)
         istart =  lxik + letk + 2
         iend   =  istart + lzek
         recv   => recv03
         drva   => drva3
      end select

      do kkk = 0,ijks(3,nn)
         kpp = kkk * ( ijks(2,nn) + 1 )
         do jjj=  0,ijks(2,nn)
            jkk = kpp + jjj
            do iii = istart,iend
               lll = indx3(iii-istart, jjj, kkk, nn, lxik, letk)
               li(iii) = lll
               sa(iii) = rfield(lll,nz)
            end do

            select case(nstart)
            case(0)
               sb(istart)   = sum( (/a01,a02,a03,a04/) * ( sa(istart+(/1,2,3,4/)) - sa(istart)   ) )
               sb(istart+1) = sum( (/a10,a12,a13,a14/) * ( sa(istart+(/0,2,3,4/)) - sa(istart+1) ) )
            case(1)
               sb(istart)   = sum( pbci(0:lmd,0,ntk) * sa(istart:istart+lmd) ) + recv(jkk,0,0)
               sb(istart+1) = sum( pbci(0:lmd,1,ntk) * sa(istart:istart+lmd) ) + recv(jkk,1,0)
            end select

            do iii = istart+2,iend-2
               sb(iii) = aa * ( sa(iii+1) - sa(iii-1) ) + ab * ( sa(iii+2) - sa(iii-2) )
            end do

            select case(nend)
            case(0)
               sb(iend)   = sum( (/a01,a02,a03,a04/) * ( sa(iend)   - sa(iend-(/1,2,3,4/)) ) )
               sb(iend-1) = sum( (/a10,a12,a13,a14/) * ( sa(iend-1) - sa(iend-(/0,2,3,4/)) ) )
            case(1)
               sb(iend)   = -sum( pbci(0:lmd,0,ntk) * sa(iend:iend-lmd:-1) ) - recv(jkk,0,1)
               sb(iend-1) = -sum( pbci(0:lmd,1,ntk) * sa(iend:iend-lmd:-1) ) - recv(jkk,1,1)
            end select

            sa(istart)   = sb(istart)
            sa(istart+1) = sb(istart+1) - xl(istart+1,2) * sa(istart)
            do iii = istart+2,iend
               sa(iii) = sb(iii) - xl(iii,1) * sa(iii-2) - xl(iii,2) * sa(iii-1)
            end do

            sb(iend)   = xu(iend,1)   * sa(iend)
            sb(iend-1) = xu(iend-1,1) * sa(iend-1) - xu(iend-1,2) * sb(iend)
            do iii = iend-2,istart,-1
               sb(iii) = xu(iii,1) * sa(iii) - xu(iii,2) * sb(iii+1) - xu(iii,3) * sb(iii+2)
            end do

            do iii = istart,iend
               lll = li(iii)
               rfield(lll,nn) = sb(iii)
            end do

            drva(jkk,m,0) = sb(istart)
            drva(jkk,m,1) = sb(iend)

         end do
      end do

   end subroutine deriv_overwrite

!===== SUBROUTINE FOR COMPACT FILTERING

   subroutine filte(rfield, lmx, lxik, letk, lzek, ijks, inn)
      real(kind=nr),    intent(inout), dimension(0:lmx) :: rfield
      integer(kind=ni), intent(in)                  :: lmx
      integer(kind=ni), intent(in)                  :: lxik, letk, lzek
      integer(kind=ni), intent(in), dimension(3,3)  :: ijks
      integer(kind=ni), intent(in)                  :: inn
      integer(kind=ni) :: nstart, nend, istart, iend, ntk, nn
      integer(kind=ni) :: kkk, jjj, iii, lll, kpp, jkk
      real(kind=nr)    :: resk, ra2k

      nn = nnf(inn)

      ntk    = 1
      nstart = ndf(nn,0,1)
      nend   = ndf(nn,1,1)

      select case(nn)
      case(1)
         istart =  0
         iend   =  istart + lxik
         recv   => recv11
      case(2)
         istart =  lxik + 1
         iend   =  istart + letk
         recv   => recv12
      case(3)
         istart =  lxik + letk + 2
         iend   =  istart + lzek
         recv   => recv13
      end select

      do kkk = 0,ijks(3,nn)
         kpp = kkk * ( ijks(2,nn) + 1 )
         do jjj = 0,ijks(2,nn)
            jkk = kpp + jjj
            do iii = istart,iend
               lll     = indx3(iii-istart, jjj, kkk, nn, lxik, letk)
               li(iii) = lll
               sa(iii) = rfield(lll)
            end do

            select case(nstart)
            case(0)
               sb(istart)   = sum( fbc(:,0) * ( sa(istart+(/1,2,3,4,5/)) - sa(istart)   ) )
               sb(istart+1) = sum( fbc(:,1) * ( sa(istart+(/0,2,3,4,5/)) - sa(istart+1) ) )
               sb(istart+2) = sum( fbc(:,2) * ( sa(istart+(/0,1,3,4,5/)) - sa(istart+2) ) )
            case(1)
               ra2k         = sa(istart+2) + sa(istart+2)
               sb(istart)   = sum( pbci(0:lmf,0,ntk) * sa(istart:istart+lmf) ) + recv(jkk,0,0)
               sb(istart+1) = sum( pbci(0:lmf,1,ntk) * sa(istart:istart+lmf) ) + recv(jkk,1,0)
               sb(istart+2) = fa * ( sa(istart+1)  + sa(istart+3) - ra2k ) + &
                              fb * ( sa(istart)    + sa(istart+4) - ra2k ) + &
                              fc * ( recv(jkk,2,0) + sa(istart+5) - ra2k )
            end select

            do iii = istart+3,iend-3
               resk   = sa(iii) + sa(iii)
               sb(iii) = fa * ( sa(iii-1) + sa(iii+1) - resk ) + &
                         fb * ( sa(iii-2) + sa(iii+2) - resk ) + &
                         fc * ( sa(iii-3) + sa(iii+3) - resk )
            end do

            select case(nend)
            case(0)
               sb(iend)   = sum( fbc(:,0) * ( sa(iend-(/1,2,3,4,5/)) - sa(iend)   ) )
               sb(iend-1) = sum( fbc(:,1) * ( sa(iend-(/0,2,3,4,5/)) - sa(iend-1) ) )
               sb(iend-2) = sum( fbc(:,2) * ( sa(iend-(/0,1,3,4,5/)) - sa(iend-2) ) )
            case(1)
               ra2k       = sa(iend-2) + sa(iend-2)
               sb(iend)   = sum( pbci(0:lmf,0,ntk) * sa(iend:iend-lmf:-1) ) + recv(jkk,0,1)
               sb(iend-1) = sum( pbci(0:lmf,1,ntk) * sa(iend:iend-lmf:-1) ) + recv(jkk,1,1)
               sb(iend-2) = fa * ( sa(iend-3) + sa(iend-1)    - ra2k ) + &
                            fb * ( sa(iend-4) + sa(iend)      - ra2k ) + &
                            fc * ( sa(iend-5) + recv(jkk,2,1) - ra2k )
            end select

            sa(istart)   = sb(istart)
            sa(istart+1) = sb(istart+1) - yl(istart+1,2) * sa(istart)
            do iii = istart+2,iend
               sa(iii) = sb(iii) - yl(iii,1) * sa(iii-2) - yl(iii,2) * sa(iii-1)
            end do
            
            sb(iend)   = yu(iend,1)   * sa(iend)
            sb(iend-1) = yu(iend-1,1) * sa(iend-1) - yu(iend-1,2) * sb(iend)
            do iii = iend-2,istart,-1
               sb(iii) = yu(iii,1) * sa(iii) - yu(iii,2) * sb(iii+1) - yu(iii,3) * sb(iii+2)
            end do

            do iii = istart,iend
               lll         = li(iii)
               rfield(lll) = rfield(lll) + sb(iii)
            end do

         end do
      end do

   end subroutine filte

!=====

end module mo_numerics

!*****