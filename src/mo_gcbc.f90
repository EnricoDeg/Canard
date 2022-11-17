!*****
!***** GENERAL CHARACTERISTIC BOUNDARY CONDITIONS MODULE
!*****

MODULE mo_gcbc
   use mo_kind,       ONLY : ni, nr, ieee32
   use mo_parameters, ONLY : one, zero, sml, pi, half, beta13, beta02,    &
                           & beta, alpha12, alpha10, alpha, alpha01, two, &
                           & hamhamm1, gam, gamm1, hamm1
   use mo_physics,    ONLY : t_physics
   use mo_numerics,   ONLY : t_numerics
   use mo_grid,       ONLY : t_grid
   use mo_domdcomp,   ONLY : t_domdcomp
   use mo_mpi,        ONLY : p_null_req, p_irecv, p_isend, p_waitall
   use mo_utils,      ONLY : indx3, mtrxi
   IMPLICIT NONE
   PUBLIC

   integer(kind=ni), private, parameter                 :: mbci=4
   real(kind=nr),    private :: vn, vs
   real(kind=nr),    private, dimension(3) :: ve
   real(kind=nr),    private, dimension(mbci,mbci)      :: cbca,cbcs
   real(kind=nr),    private, dimension(mbci)           :: rbci,sbci
   real(kind=nr),    private, dimension(:), allocatable :: sbcc
   integer(kind=ni), private, dimension(:), allocatable :: nrr,npex
   real(kind=nr),    private, dimension(:,:,:), pointer :: drvb
   real(kind=nr),    private, dimension(:,:,:), allocatable, target :: drvb1,drvb2,drvb3
   real(kind=nr),    private, dimension(3) :: dm
   real(kind=nr),    private, dimension(5,5) :: xt
   real(kind=nr),    private, dimension(5) :: cha, dha
   real(kind=nr),    private :: aoi, hv2, ao

   private :: eleme, xtq2r, xtr2q

   CONTAINS

   SUBROUTINE gcbc_init(p_domdcomp, yaco)
      type(t_domdcomp),                              intent(IN) :: p_domdcomp
      real(kind=nr),    dimension(0:p_domdcomp%lmx), intent(in) :: yaco

      real(kind=nr), dimension(0:p_domdcomp%lmx,3) :: rr
      integer(kind=ni) :: ii, jj, kk, nn, np, lq, ll, l, ip, iq
      integer(kind=ni) :: i, j, k
      real(kind=nr)    :: res, fctr

      allocate(nrr(0:p_domdcomp%lmx), npex(0:p_domdcomp%lmx))

      ii = p_domdcomp%nbsize(1) - 1
      jj = p_domdcomp%nbsize(2) - 1
      kk = p_domdcomp%nbsize(3) - 1
      allocate(drvb1(0:ii,5,0:1), drvb2(0:jj,5,0:1), drvb3(0:kk,5,0:1))

      cbca(:,:)   = zero
      cbca(1,1:2) = (/ alpha01, beta02 /)
      cbca(2,1:3) = (/ one, alpha12, beta13 /)
      cbca(3,1:4) = (/ alpha, one, alpha, beta /)
      do i=4,mbci
         cbca(i,i-3:i)=(/beta,alpha,one,alpha/)
         if(i<mbci) then
            cbca(i,i+1)=beta
         end if
      end do

      rbci(:)   = zero
      rbci(1:3) = (/ one, alpha10, beta /)
      call mtrxi(cbca, cbcs, 1, mbci)
      sbci(:)   = -matmul(cbcs(:,:), rbci(:))
      fctr = pi / ( mbci + 1 )
      res  = zero
      do i=1,mbci
         res     = res + one
         sbci(i) = half * sbci(i) * ( one + cos( res * fctr ) )
      end do
      
      ll      = -1
      npex(:) = 0
      nrr(:)  = 0
      rr(:,1) = zero
      do nn=1,3
         do ip=0,1
            np = p_domdcomp%nbc(nn,ip)
            i  = ip * p_domdcomp%ijk(1,nn)
            iq = 1 - 2 * ip
            if ( ( np - 10 ) * ( np - 20 ) * ( np - 25 ) * ( np - 30 ) == 0 ) then
               do k=0,p_domdcomp%ijk(3,nn)
                  do j=0,p_domdcomp%ijk(2,nn)
                     l   = indx3(i, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                     ll  = ll + 1
                     res = one / yaco(l)
                     rr(l,1)  = rr(l,1) + one
                     rr(ll,2) = res
                     rr(ll,3) = l + sml
                     do ii=1,mbci
                        l  = indx3(i+iq*ii, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                        ll = ll + 1
                        rr(l,1)  = rr(l,1) + one
                        rr(ll,2) = res * sbci(ii)
                        rr(ll,3) = l + sml
                     end do
                  end do
               end do
            end if
            if ( ( np - 30 ) * ( np - 35 ) * ( np - 45 ) == 0 ) then
               do k=0,p_domdcomp%ijk(3,nn)
                  do j=0,p_domdcomp%ijk(2,nn)
                     l      = indx3(i, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                     nrr(l) = 1
                  end do
               end do
            end if
         end do
      end do
      do l=0,p_domdcomp%lmx
         nrr(l) = min(nrr(l)+npex(l),1)
      end do
      lq = ll
      allocate(sbcc(0:lq))
      do ll=0,lq
         l        = rr(ll,3)
         sbcc(ll) = rr(ll,2) / rr(l,1)
      end do

   END SUBROUTINE gcbc_init

   SUBROUTINE gcbc_setup(p_domdcomp, p_numerics, p_grid, qa, p, de, umf)
      type(t_domdcomp), intent(IN)    :: p_domdcomp
      type(t_numerics), intent(INOUT) :: p_numerics
      type(t_grid),     intent(IN)    :: p_grid
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(in) :: qa
      real(kind=nr), dimension(0:p_domdcomp%lmx), intent(in) :: p
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(in) :: de
      real(kind=nr), dimension(3), intent(in) :: umf
      integer(kind=ni) :: nn, np, l, ip, i, j, k, jk, kp
      real(kind=nr)    :: ra0, ra1
      real(kind=nr), dimension(:,:,:), pointer :: cm
      real(kind=nr), dimension(:,:,:), pointer :: drva

      do nn=1,3
         select case(nn)
         case(1)
            drva => p_numerics%drva1
            cm   => p_grid%cm1
         case(2)
            drva => p_numerics%drva2
            cm   => p_grid%cm2
         case(3)
            drva => p_numerics%drva3
            cm   => p_grid%cm3
         end select
         do ip=0,1
            np = p_domdcomp%nbc(nn,ip)
            i  = ip * p_domdcomp%ijk(1,nn)
            if ( ( np - 10 ) * ( np - 20 ) * ( np - 25 ) * ( np - 30 ) == 0 ) then
               ra0 = ( 20 - np ) * ( 25 - np ) * ( 30 - np ) / 3000
               ra1 = one - ra0
               do k=0,p_domdcomp%ijk(3,nn)
                  kp = k * ( p_domdcomp%ijk(2,nn) + 1 )
                  do j=0,p_domdcomp%ijk(2,nn)
                     jk = kp + j
                     l  = indx3(i, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                     call eleme(cm(jk,:,ip), qa(l,1), qa(l,2:4), p(l), umf)
                     call xtq2r(cm(jk,:,ip))
                     cha(:) = ra0 * drva(jk,:,ip) + ra1 * de(l,:)
                     drva(jk,:,ip) = matmul( xt(:,:), p_grid%yaco(l) * cha(:) )
                  end do
               end do
            end if
         end do
      end do

   END SUBROUTINE gcbc_setup

   SUBROUTINE gcbc_comm(p_domdcomp, p_numerics)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      type(t_numerics), intent(inout) :: p_numerics
      integer(kind=ni) :: nn, np, ip, iq, itag
      real(kind=nr), dimension(:,:,:), pointer :: drva

      call p_null_req
      itag = 30
      do nn=1,3
         select case(nn)
         case(1)
            drva => p_numerics%drva1
            drvb => drvb1
         case(2)
            drva => p_numerics%drva2
            drvb => drvb2
         case(3)
            drva => p_numerics%drva3
            drvb => drvb3
         end select
         do ip=0,1
            iq=1-ip
            np = p_domdcomp%nbc(nn,ip)
            if ( ( np - 30 ) * ( one + abs((np-20)*(np-25)) ) == 0 ) then
               call p_isend(drva(:,:,ip), p_domdcomp%mcd(nn,ip), itag+iq, &
                            5*p_domdcomp%nbsize(nn))
               call p_irecv(drvb(:,:,ip), p_domdcomp%mcd(nn,ip), itag+ip, &
                            5*p_domdcomp%nbsize(nn))
            end if
         end do
      end do
      call p_waitall

   END SUBROUTINE gcbc_comm

   SUBROUTINE gcbc_update(p_domdcomp, p_numerics, p_grid, qa, p, de, nkrk, dt, umf, dudtmf)
      type(t_domdcomp), intent(IN)    :: p_domdcomp
      type(t_numerics), intent(inout) :: p_numerics
      type(t_grid),     intent(IN)    :: p_grid
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(in) :: qa
      real(kind=nr), dimension(0:p_domdcomp%lmx), intent(in) :: p
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(inout) :: de
      integer(kind=ni), intent(in) :: nkrk
      real(kind=nr),    intent(in) :: dt
      real(kind=nr), dimension(3), intent(in) :: umf
      real(kind=nr), dimension(3), intent(in) :: dudtmf
      integer(kind=ni) :: ii, nn, np, ll, l, ip, iq, i, j, k
      integer(kind=ni) :: jk, kp
      real(kind=nr)    :: ra0, dtwi
      real(kind=nr), dimension(:,:,:), pointer :: cm
      real(kind=nr), dimension(:,:,:), pointer :: drva

      ll = -1
      do nn=1,3
         select case(nn)
         case(1)
            drva => p_numerics%drva1
            drvb => drvb1
            cm   => p_grid%cm1
         case(2)
            drva => p_numerics%drva2
            drvb => drvb2
            cm   => p_grid%cm2
         case(3)
            drva => p_numerics%drva3
            drvb => drvb3
            cm   => p_grid%cm3
         end select
         do ip=0,1
            np  = p_domdcomp%nbc(nn,ip)
            i   = ip * p_domdcomp%ijk(1,nn)
            iq  = 1 - 2 * ip
            ra0 = iq
            select case(np)
            case(10)
               do k=0,p_domdcomp%ijk(3,nn)
                  kp = k * ( p_domdcomp%ijk(2,nn) + 1 )
                  do j=0,p_domdcomp%ijk(2,nn)
                     jk = kp + j
                     l  = indx3(i, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                     call eleme(cm(jk,:,ip), qa(l,1), qa(l,2:4), p(l), umf)
                     cha(:) = drva(jk,:,ip)
                     dha(:) = drvb(jk,:,ip)
                     if ( ra0 * ( vn + vs + ao ) > zero ) then
                        cha(4) = zero
                     end if
                     if ( ra0 * ( vn + vs - ao ) > zero ) then
                        cha(5) = zero
                     end if
                     call xtr2q(cm(jk,:,ip))
                     dha(:) = matmul(xt(:,:),(cha(:)-drva(jk,:,ip)))
                     do ii=0,mbci
                        l  = indx3(i+iq*ii, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                        ll = ll + 1
                        de(l,:) = de(l,:) + sbcc(ll) * dha(:)
                     end do
                  end do
               end do
            case(20,25)
               dtwi = one / ( nkrk * dt + sml )
               do k=0,p_domdcomp%ijk(3,nn)
                  kp = k * ( p_domdcomp%ijk(2,nn) + 1 )
                  do j=0,p_domdcomp%ijk(2,nn)
                     jk = kp + j
                     l  = indx3(i, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                     call eleme(cm(jk,:,ip), qa(l,1), qa(l,2:4), p(l), umf)
                     cha(:) = drva(jk,:,ip)
                     dha(:) = drvb(jk,:,ip)
                     select case(npex(l))
                     case(0)
                        cha(4+ip) = cha(5-ip) + two * ra0 * aoi * qa(l,1) * &
                                   ( sum(cm(jk,:,ip)*dudtmf(:)) + dtwi * ( vn + vs ) )
                     end select
                     call xtr2q(cm(jk,:,ip))
                     dha(:) = matmul(xt(:,:),(cha(:)-drva(jk,:,ip)))
                     do ii=0,mbci
                        l  = indx3(i+iq*ii, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                        ll = ll + 1
                        de(l,:) = de(l,:) + sbcc(ll) * dha(:)
                     end do
                  end do
               end do
            case(30)
               do k=0,p_domdcomp%ijk(3,nn)
                  kp = k * ( p_domdcomp%ijk(2,nn) + 1 )
                  do j=0,p_domdcomp%ijk(2,nn)
                     jk = kp + j
                     l  = indx3(i, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                     call eleme(cm(jk,:,ip), qa(l,1), qa(l,2:4), p(l), umf)
                     cha(:) = drva(jk,:,ip)
                     dha(:) = drvb(jk,:,ip)
                     if ( ra0 * ( vn + vs ) > zero ) then
                        cha(1:3) = dha(1:3)
                     end if
                     if ( ra0 * ( vn + vs + ao ) > zero ) then
                        cha(4) = dha(4)
                     end if
                     if ( ra0 * ( vn + vs - ao ) > zero ) then
                        cha(5) = dha(5)
                     end if
                     call xtr2q(cm(jk,:,ip))
                     dha(:) = matmul(xt(:,:),(cha(:)-drva(jk,:,ip)))
                     do ii=0,mbci
                        l  = indx3(i+iq*ii, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                        ll = ll + 1
                        de(l,:) = de(l,:) + sbcc(ll) * dha(:)
                     end do
                  end do
               end do
            end select
         end do
      end do

   END SUBROUTINE gcbc_update

   SUBROUTINE average_surface(p_domdcomp, p_numerics, qa)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      type(t_numerics), intent(inout) :: p_numerics
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(inout) :: qa

      integer(kind=ni) :: nn, np, l, ip, iq, i, j, k, jk, kp
      integer(kind=ni) :: itag
      real(kind=nr), dimension(:,:,:), pointer :: drva
      real(kind=nr), dimension(0:p_domdcomp%lmx,3) :: rr

      call p_null_req
      itag = 30
      do nn=1,3
         select case(nn)
         case(1)
            drva => p_numerics%drva1
            drvb => drvb1
         case(2)
            drva => p_numerics%drva2
            drvb => drvb2
         case(3)
            drva => p_numerics%drva3
            drvb => drvb3
         end select
         do ip=0,1
            iq = 1 - ip
            np = p_domdcomp%nbc(nn,ip)
            i  = ip * p_domdcomp%ijk(1,nn)
            if ( ( np - 30 ) * ( np - 35 ) * ( np - 45 ) * &
                ( one + abs((np-20)*(np-25)) ) == 0 ) then
               do k=0,p_domdcomp%ijk(3,nn)
                  kp = k * ( p_domdcomp%ijk(2,nn) + 1 )
                  do j=0,p_domdcomp%ijk(2,nn)
                     l  = indx3(i, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                     jk = kp + j
                     drva(jk,:,ip) = qa(l,:)
                     rr(l,1) = one
                  end do
               end do
               call p_isend(drva(:,:,ip), p_domdcomp%mcd(nn,ip), itag+iq, &
                            5*p_domdcomp%nbsize(nn))
               call p_irecv(drvb(:,:,ip), p_domdcomp%mcd(nn,ip), itag+ip, &
                            5*p_domdcomp%nbsize(nn))
            end if
         end do
      end do
      call p_waitall
      do nn=1,3
         select case(nn)
         case(1)
            drvb => drvb1
         case(2)
            drvb => drvb2
         case(3)
            drvb => drvb3
         end select
         do ip=0,1
            np = p_domdcomp%nbc(nn,ip)
            i  = ip * p_domdcomp%ijk(1,nn)
            if( ( np - 30 ) * ( np - 35 ) * ( np - 45 ) * &
                ( one + abs((np-20)*(np-25)) ) == 0) then
               do k=0,p_domdcomp%ijk(3,nn)
                  kp = k * ( p_domdcomp%ijk(2,nn) + 1 )
                  do j=0,p_domdcomp%ijk(2,nn)
                     l  = indx3(i, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                     jk = kp + j
                     rr(l,1) = rr(l,1) + nrr(l)
                     rr(l,2) = nrr(l) / rr(l,1)
                     qa(l,:) = rr(l,2) * ( ( rr(l,1) - one ) * qa(l,:) + &
                               drvb(jk,:,ip) ) + ( 1 - nrr(l) ) * qa(l,:)
                  end do
               end do
            end if
         end do
      end do

   END SUBROUTINE average_surface

   SUBROUTINE wall_condition_update(p_domdcomp, qa, umf)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(inout) :: qa
      real(kind=nr), dimension(3), intent(in) :: umf
      integer(kind=ni) :: nn, np, l, ip, i, j, k
      real(kind=nr)    :: fctr, ra0

      do nn=1,3
         do ip=0,1
            np = p_domdcomp%nbc(nn,ip)
            i  = ip * p_domdcomp%ijk(1,nn)
            select case(np)
            case(20)
               do k=0,p_domdcomp%ijk(3,nn)
                  do j=0,p_domdcomp%ijk(2,nn)
                     l = indx3(i, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                     qa(l,5) = npex(l) * qa(l,5) + ( one - npex(l) ) * & 
                               ( hamhamm1 * qa(l,1)**gam + half * sum( qa(l,2:4) * qa(l,2:4) ) / qa(l,1) )
                  end do
               end do
            case(25)
               ra0 = hamhamm1
               do k=0,p_domdcomp%ijk(3,nn)
                  do j=0,p_domdcomp%ijk(2,nn)
                     l         = indx3(i, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                     fctr      = ( one - npex(l) ) * qa(l,1)
                     qa(l,2:4) = npex(l) * qa(l,2:4) - fctr * umf(:)
                     qa(l,5)   = npex(l) * qa(l,5) + ( one - npex(l) ) * &
                                 ( ra0 * qa(l,1) + half * sum( qa(l,2:4) * qa(l,2:4) ) / qa(l,1) )
                  end do
               end do
            end select
         end do
      end do

   END SUBROUTINE wall_condition_update

!===== EXTRA CONDITION

   subroutine extracon(p_domdcomp, p_grid, p_physics, qa, p)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      type(t_grid),     intent(IN) :: p_grid
      type(t_physics),     intent(IN) :: p_physics
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(in) :: qa
      real(kind=nr), dimension(0:p_domdcomp%lmx), intent(in) :: p
      
      real(kind=nr),dimension(3) :: vee
      integer(kind=ni) :: nn, l, ip, i, j, k, jk, kp
      real(kind=nr)    :: fctr, ra0, ra1, ra2, ra3
      real(kind=nr),dimension(3) :: rv
      
      nn = 2
      ip = 0
      i  = ip * p_domdcomp%ijk(1,nn)
      if ( p_domdcomp%nbc(nn,ip) == 25 ) then
         fctr = one / ( p_domdcomp%ijk(2,nn) + 1 )
         do k=0,p_domdcomp%ijk(3,nn)
            kp    = k * ( p_domdcomp%ijk(2,nn) + 1 )
            rv(:) = zero
            do j=0,p_domdcomp%ijk(2,nn)
               jk     = kp + j
               l      = indx3(i, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
               ra0    = two * acos(p_grid%cm2(jk,1,ip))
               ra1    = abs( half * sin(ra0) * ( p_physics%tyy(l) - p_physics%txx(l) ) + cos(ra0) * p_physics%txy(l) )
               ra2    = gam * p(l) / qa(l,1)
               ra3    = sqrt( ra1 * qa(l,1) ) * ( ra2 + p_physics%srefoo ) / ( p_physics%srefp1dre * ra2**1.5_nr )
               vee(1) = sqrt( ( p_grid%etm(l,2) * p_grid%zem(l,3) - p_grid%zem(l,2) * p_grid%etm(l,3) )**two + &
                              ( p_grid%etm(l,3) * p_grid%zem(l,1) - p_grid%zem(l,3) * p_grid%etm(l,1) )**two )
               vee(2) = p_grid%cm2(jk,1,ip) * ( p_grid%zem(l,2) * p_grid%xim(l,3)   - &
                                                p_grid%xim(l,2) * p_grid%zem(l,3) ) + &
                        p_grid%cm2(jk,2,ip) * ( p_grid%zem(l,3) * p_grid%xim(l,1)   - &
                                                p_grid%xim(l,3) * p_grid%zem(l,1) )
               vee(3) = p_grid%xim(l,1) * p_grid%etm(l,2) - p_grid%etm(l,1) * p_grid%xim(l,2)
               rv(:)  = rv(:) + ra3 * abs( vee(:) * p_grid%yaco(l) )
            end do
         end do
      end if
      
   end subroutine extracon

!===== SUBROUTINE FOR TRANSFORMATION FROM Q TO R IN GCBC/GCIC

   subroutine xtq2r(cm)

      real(kind=nr),dimension(3),intent(in) :: cm
      real(kind=nr),dimension(3) :: rv
      real(kind=nr) :: ho, bo, co
      
      ho    = gamm1 * aoi * aoi
      bo    = 1 - ho * hv2
      co    = aoi * vn
      dm(:) = aoi * cm(:)
      rv(:) = ho * ve(:)
      
      xt(1,1) = bo * cm(1) + dm(2) * ve(3) - dm(3) * ve(2)
      xt(1,2) = cm(1) * rv(1)
      xt(1,3) = cm(1) * rv(2) + dm(3)
      xt(1,4) = cm(1) * rv(3) - dm(2)
      xt(1,5) = -ho * cm(1)
      
      xt(2,1) = bo * cm(2) + dm(3) * ve(1) - dm(1) * ve(3)
      xt(2,2) = cm(2) * rv(1) - dm(3)
      xt(2,3) = cm(2) * rv(2)
      xt(2,4) = cm(2) * rv(3) + dm(1)
      xt(2,5) = -ho * cm(2)
      
      xt(3,1) = bo * cm(3) + dm(1) * ve(2) - dm(2) * ve(1)
      xt(3,2) = cm(3) * rv(1) + dm(2)
      xt(3,3) = cm(3) * rv(2) - dm(1)
      xt(3,4) = cm(3) * rv(3)
      xt(3,5) = -ho * cm(3)
      
      xt(4,1) = one - bo - co
      xt(4,2) = dm(1) - rv(1)
      xt(4,3) = dm(2) - rv(2)
      xt(4,4) = dm(3) - rv(3)
      xt(4,5) = ho
      
      xt(5,1) = one - bo + co
      xt(5,2) = -dm(1) - rv(1)
      xt(5,3) = -dm(2) - rv(2)
      xt(5,4) = -dm(3) - rv(3)
      xt(5,5) = ho
      
   end subroutine xtq2r

!===== SUBROUTINE FOR INVERSE TRANSFORMATION FROM R TO Q IN GCBC/GCIC

   subroutine xtr2q(cm)

      real(kind=nr),dimension(3),intent(in) :: cm
      real(kind=nr) :: bo, co

      bo    = hv2 + hamm1 * ao * ao
      co    = ao * vn
      dm(:) = ao * cm(:)

      xt(1,1) = cm(1)
      xt(1,2) = cm(2)
      xt(1,3) = cm(3)
      xt(1,4) = half
      xt(1,5) = half

      xt(2,1) = cm(1) * ve(1)
      xt(2,2) = cm(2) * ve(1) - dm(3)
      xt(2,3) = cm(3) * ve(1) + dm(2)
      xt(2,4) = half * ( ve(1) + dm(1) )
      xt(2,5) = xt(2,4) - dm(1)

      xt(3,1) = cm(1) * ve(2) + dm(3)
      xt(3,2) = cm(2) * ve(2)
      xt(3,3) = cm(3) * ve(2) - dm(1)
      xt(3,4) = half * ( ve(2) + dm(2) )
      xt(3,5) = xt(3,4) - dm(2)

      xt(4,1) = cm(1) * ve(3) - dm(2)
      xt(4,2) = cm(2) * ve(3) + dm(1)
      xt(4,3) = cm(3) * ve(3)
      xt(4,4) = half * ( ve(3) + dm(3) )
      xt(4,5) = xt(4,4) -  dm(3)

      xt(5,1) = hv2 * cm(1) + dm(3) * ve(2) - dm(2) * ve(3)
      xt(5,2) = hv2 * cm(2) + dm(1) * ve(3) - dm(3) * ve(1)
      xt(5,3) = hv2 * cm(3) + dm(2) * ve(1) - dm(1) * ve(2)
      xt(5,4) = half * ( bo + co )
      xt(5,5) = xt(5,4) - co

   end subroutine xtr2q

!===== SUBROUTINE FOR ELEMENTARY VARIABLES IN GCBC/GCIC

   subroutine eleme(cm, qa1, qa24, pp, umf)
      real(kind=nr), dimension(3), intent(in) :: cm
      real(kind=nr),               intent(in) :: qa1
      real(kind=nr), dimension(3), intent(in) :: qa24
      real(kind=nr),               intent(in) :: pp
      real(kind=nr), dimension(3), intent(in) :: umf
      real(kind=nr) :: rhoi

      rhoi  = one / qa1
      ao    = sqrt( gam * rhoi * pp )
      aoi   = one / ao
      ve(:) = rhoi * qa24(:)
      hv2   = half * ( ve(1) * ve(1) + ve(2) * ve(2) + ve(3) * ve(3) )
      vn    = cm(1) * ve(1) + cm(2) * ve(2) + cm(3) * ve(3)
      vs    = cm(1) * umf(1) + cm(2) * umf(2) + cm(3) * umf(3)

   end subroutine eleme

END MODULE mo_gcbc