!*****
!***** PHYSICS MODULE
!*****

MODULE mo_physics
   use mo_kind,       ONLY : nr, ni
   use mo_parameters, ONLY : sml, zero, one, pi, hamm1, hamhamm1, half, gam,      &
                           & gamm1, n45no, nrall, gamm1prndtli, nrone, twothirds
   use mo_mpi,        ONLY : p_max
   use mo_grid,       ONLY : t_grid
   use mo_domdcomp,   ONLY : t_domdcomp
   use mo_numerics,   ONLY : t_numerics
   implicit none
   public

   type, public :: t_physics
      ! private variables
      real(kind=nr),    private               :: reoo
      real(kind=nr),    private               :: tempoo
      real(kind=nr),    private               :: amach1, amach2, amach3
      real(kind=nr),    private               :: amachoo
      real(kind=nr),    private               :: timf
      integer(kind=ni), private               :: nsmf
      real(kind=nr),    private               :: sqrtrema, sqrtremai
      real(kind=nr),    private, dimension(3) :: uoo

      ! public variables
      real(kind=nr),    public                            :: srefoo, srefp1dre
      real(kind=nr),    public, dimension(3)              :: umf, dudtmf
      real(kind=nr),    public, dimension(:), allocatable :: txx, tyy, tzz, txy, tyz, tzx
      real(kind=nr),    public, dimension(:), allocatable :: hxx, hyy, hzz

   contains

      procedure, public :: allocate => allocate_physics_memory
      procedure, public :: deallocate => deallocate_physics_memory
      procedure, public :: read => read_input_physics
      procedure, public :: init => initialo
      procedure, public :: movef
      procedure, public :: calc_viscous_shear_stress
      procedure, public :: calc_fluxes
      procedure, public :: calc_time_step

   end type t_physics

   contains

   subroutine allocate_physics_memory(this, lmx)
      class(t_physics), INTENT(INOUT) :: this
      integer(kind=ni), intent(IN) :: lmx

#ifdef VISCOUS
      allocate(this%txx(0:lmx), this%tyy(0:lmx), this%tzz(0:lmx))
      allocate(this%txy(0:lmx), this%tyz(0:lmx), this%tzx(0:lmx))
      allocate(this%hxx(0:lmx), this%hyy(0:lmx), this%hzz(0:lmx))
#endif

   end subroutine allocate_physics_memory



   subroutine deallocate_physics_memory(this)
      class(t_physics), INTENT(INOUT) :: this

#ifdef VISCOUS
      deallocate(this%txx, this%tyy, this%tzz)
      deallocate(this%txy, this%tyz, this%tzx)
      deallocate(this%hxx, this%hyy, this%hzz)
#endif

   end subroutine deallocate_physics_memory

!===== INITIALIZE PHYSICS

   subroutine read_input_physics(this)
      class(t_physics), INTENT(INOUT) :: this
      character(16)    :: cinput
      real(kind=nr)    :: reoo
      real(kind=nr)    :: tempoo
      real(kind=nr)    :: amach1, amach2, amach3
      real(kind=nr)    :: amachoo
      real(kind=nr)    :: timf
      integer(kind=ni) :: nsmf
      integer(kind=ni) :: rc, fu

      namelist /nml_physics/ reoo, tempoo, amach1, amach2, amach3, &
                              timf, nsmf

      open (action='read', file='input.canard', iostat=rc, newunit=fu)
      read (nml=nml_physics, iostat=rc, unit=fu)
      close(fu)

      this%reoo   = reoo
      this%tempoo = tempoo
      this%amach1 = amach1
      this%amach2 = amach2
      this%amach3 = amach3
      this%timf   = timf
      this%nsmf   = nsmf
      
      this%amachoo = sqrt( this%amach1 * this%amach1 + &
                           this%amach2 * this%amach2 + &
                           this%amach3 * this%amach3 )
      if ( this%amachoo > sml ) then
         this%reoo = this%reoo / this%amachoo
      end if
      this%srefoo    = 111 / this%tempoo
      this%srefp1dre = ( this%srefoo + one ) / this%reoo
      this%sqrtrema  = sqrt( this%reoo * this%amachoo )
      this%sqrtremai = one / max( this%sqrtrema, sml )
      this%uoo(:)    = (/ this%amach1, this%amach2, this%amach3 /)

   end subroutine read_input_physics

!===== INITIAL CONDITIONS

   subroutine initialo(this, lmx, qa, ss)
      class(t_physics), INTENT(INOUT) :: this
      integer(kind=ni), intent(in) :: lmx
      real(kind=nr), dimension(0:lmx,5), intent(inout) :: qa
      real(kind=nr), dimension(0:lmx,3), intent(in)    :: ss
      real(kind=nr),dimension(3) :: vee
      real(kind=nr) :: radv, k1, k2, bo, hv2, ao
      integer(kind=ni) :: l

      radv = 1.0
      k1   = 12.5
      k2   = 1.0

      do l=0,lmx
         ao        = k2/2.0/pi * sqrt( exp( 1 - k1**2 * ( ss(l,1)**2 + ss(l,2)**2 ) / radv**2 ) )
         bo        = ( one - half * gamm1 * ao * ao )**hamm1
         qa(l,1)   = bo
         vee(:)    = (/ k1 * ss(l,2) * ao / radv, -k1 * ss(l,1) * ao / radv, zero/)
         hv2       = half * ( vee(1) * vee(1) + vee(2) * vee(2) + vee(3) * vee(3) )
         qa(l,2:4) = bo * vee(:)
         qa(l,5)   = hamhamm1 * bo**gam + hv2 * bo
      end do

   end subroutine initialo

!===== SUBROUTINE FOR MOVING FRAME VELOCITIES

   subroutine movef(this, dtko, dtk, timo)
      class(t_physics), INTENT(INOUT) :: this
      real(kind=nr), intent(in) :: dtko, dtk
      real(kind=nr), intent(in) :: timo
      real(kind=nr) :: ra0, ra1, ra2, dfdt, fctr, progmf

      if ( this%nsmf == 0 ) then
         ra0 = pi / this%timf
         ra1 = ra0 * min( timo, this%timf )
         ra2 = ra0 * min( timo+dtko, this%timf )

         fctr   = one - cos(ra1)
         dfdt   = ra0 * sin(ra2)
         progmf = half * ( fctr + dtk * dfdt )
         this%umf(:) = progmf * this%uoo(:)

         fctr      = sin(ra1)
         dfdt      = ra0 * cos(ra2)
         progmf    = half * ra0 * ( fctr + dtk * dfdt )
         this%dudtmf(:) = progmf * this%uoo(:)
      else
         this%umf(:)    = this%uoo(:)
         this%dudtmf(:) = zero
      end if

   end subroutine movef

   subroutine calc_time_step(this, lmx, p_grid, de, ssk, cfl, dte)
      class(t_physics), INTENT(INOUT) :: this
      integer(kind=ni), intent(in) :: lmx
      type(t_grid),     intent(in) :: p_grid
      real(kind=nr), dimension(0:lmx,5), intent(in) :: de
      real(kind=nr), dimension(0:lmx), intent(in) :: ssk
      real(kind=nr), intent(in) :: cfl
      real(kind=nr), intent(out) :: dte

      real(kind=nr) :: ra0, ra1, res, fctr
      real(kind=nr), dimension(0:lmx,3) :: rr
      real(kind=nr), dimension(0:lmx) :: ssi
      

      rr(:,1) = p_grid%xim(:,1) * p_grid%xim(:,1) + p_grid%xim(:,2) * p_grid%xim(:,2) + &
                p_grid%xim(:,3) * p_grid%xim(:,3) + p_grid%etm(:,1) * p_grid%etm(:,1) + &
                p_grid%etm(:,2) * p_grid%etm(:,2) + p_grid%etm(:,3) * p_grid%etm(:,3) + &
                p_grid%zem(:,1) * p_grid%zem(:,1) + p_grid%zem(:,2) * p_grid%zem(:,2) + &
                p_grid%zem(:,3) * p_grid%zem(:,3)
      rr(:,2) = abs( p_grid%xim(:,1) * ( de(:,2) + this%umf(1) )   + &
                     p_grid%xim(:,2) * ( de(:,3) + this%umf(2) )   + &
                     p_grid%xim(:,3) * ( de(:,4) + this%umf(3) ) ) + &
                abs( p_grid%etm(:,1) * ( de(:,2) + this%umf(1) )   + &
                     p_grid%etm(:,2) * ( de(:,3) + this%umf(2) )   + & 
                     p_grid%etm(:,3) * ( de(:,4) + this%umf(3) ) ) + &
                abs( p_grid%zem(:,1) * ( de(:,2) + this%umf(1) )   + &
                     p_grid%zem(:,2) * ( de(:,3) + this%umf(2) )   + &
                     p_grid%zem(:,3) * ( de(:,4) + this%umf(3) ) )
      ssi(:) = abs( p_grid%yaco(:) )
      res     = maxval( ( sqrt( de(:,5) * rr(:,1) ) + rr(:,2) ) * ssi(:) )
      call p_max(res, fctr)
      ra0 = cfl / fctr
      ra1 = ra0
#ifdef VISCOUS
      res = maxval( de(:,1) * ssk * rr(:,1) * ssi(:) * ssi(:) )
      call p_max(res, fctr)
      ra1 = half / fctr
#endif
      dte = min(ra0, ra1)
      
   end subroutine calc_time_step


!===== VISCOUS SHEAR STRESSES & HEAT FLUXES

   subroutine calc_viscous_shear_stress(this, p_domdcomp, p_numerics, p_grid, de, ssk, luse_acc)
      class(t_physics), intent(inout)                             :: this
      type(t_domdcomp), intent(inout)                             :: p_domdcomp
      type(t_numerics), intent(inout)                             :: p_numerics
      type(t_grid),     intent(inout)                             :: p_grid
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(inout) :: de
      real(kind=nr), dimension(0:p_domdcomp%lmx),   intent(in)    :: ssk
      logical,          intent(in), optional                      :: luse_acc

      real(kind=nr), dimension(0:p_domdcomp%lmx,3) :: ss
      real(kind=nr), dimension(0:p_domdcomp%lmx,3,2:5) :: rr
      integer(kind=ni) :: m, i
      logical          :: lacc

#ifdef VISCOUS
      if (present(luse_acc)) then
        lacc = luse_acc
      else
        lacc = .false.
      end if

      !$ACC PARALLEL LOOP GANG VECTOR IF (lacc)
      do i=0,p_domdcomp%lmx
         de(i,1) = ssk(i)
      end do
      !$ACC END PARALLEL

      ! Halo exchange
      do m=2,5
         !$ACC PARALLEL LOOP GANG VECTOR IF (lacc)
         do i=0,p_domdcomp%lmx
            rr(i,1,m) = de(i,m)
         end do
         !$ACC END PARALLEL
         call p_numerics%mpigo(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                               p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrone, n45no, m, &
                               p_domdcomp%lxi, p_domdcomp%let, m, luse_acc = .false.)
      end do

      m = 2
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 1, m, luse_acc = lacc)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 1, m, luse_acc = lacc)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m, luse_acc = lacc)
      !$ACC PARALLEL LOOP GANG VECTOR IF (lacc)
      do i=0,p_domdcomp%lmx
         this%txx(i) = p_grid%xim(i,1) * rr(i,1,m) + &
                       p_grid%etm(i,1) * rr(i,2,m) + &
                       p_grid%zem(i,1) * rr(i,3,m)
         this%hzz(i) = p_grid%xim(i,2) * rr(i,1,m) + &
                       p_grid%etm(i,2) * rr(i,2,m) + &
                       p_grid%zem(i,2) * rr(i,3,m)
         this%tzx(i) = p_grid%xim(i,3) * rr(i,1,m) + &
                       p_grid%etm(i,3) * rr(i,2,m) + &
                       p_grid%zem(i,3) * rr(i,3,m)
      end do
      !$ACC END PARALLEL

      m = 3
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 1, m, luse_acc = lacc)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, & 
                     p_domdcomp%ijk, 2, 1, m, luse_acc = lacc)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m, luse_acc = lacc)
      !$ACC PARALLEL LOOP GANG VECTOR IF (lacc)
      do i=0,p_domdcomp%lmx
         this%txy(i) = p_grid%xim(i,1) * rr(i,1,m) + &
                       p_grid%etm(i,1) * rr(i,2,m) + &
                       p_grid%zem(i,1) * rr(i,3,m)
         this%tyy(i) = p_grid%xim(i,2) * rr(i,1,m) + &
                       p_grid%etm(i,2) * rr(i,2,m) + &
                       p_grid%zem(i,2) * rr(i,3,m)
         this%hxx(i) = p_grid%xim(i,3) * rr(i,1,m) + & 
                       p_grid%etm(i,3) * rr(i,2,m) + &
                       p_grid%zem(i,3) * rr(i,3,m)
      end do
      !$ACC END PARALLEL

      m = 4
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 1, m, luse_acc = lacc)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 1, m, luse_acc = lacc)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m, luse_acc = lacc)
      !$ACC PARALLEL LOOP GANG VECTOR IF (lacc)
      do i=0,p_domdcomp%lmx
         this%hyy(i) = p_grid%xim(i,1) * rr(i,1,m) + &
                       p_grid%etm(i,1) * rr(i,2,m) + &
                       p_grid%zem(i,1) * rr(i,3,m)
         this%tyz(i) = p_grid%xim(i,2) * rr(i,1,m) + &
                       p_grid%etm(i,2) * rr(i,2,m) + &
                       p_grid%zem(i,2) * rr(i,3,m)
         this%tzz(i) = p_grid%xim(i,3) * rr(i,1,m) + &
                       p_grid%etm(i,3) * rr(i,2,m) + &
                       p_grid%zem(i,3) * rr(i,3,m)
      end do
      !$ACC END PARALLEL

      m = 5
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 1, m, luse_acc = lacc)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 1, m, luse_acc = lacc)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m, luse_acc = lacc)
      !$ACC PARALLEL LOOP GANG VECTOR IF (lacc)
      do i=0,p_domdcomp%lmx
         ss(i,1) = p_grid%xim(i,1) * rr(i,1,m) + &
                   p_grid%etm(i,1) * rr(i,2,m) + &
                   p_grid%zem(i,1) * rr(i,3,m)
         ss(i,2) = p_grid%xim(i,2) * rr(i,1,m) + &
                   p_grid%etm(i,2) * rr(i,2,m) + &
                   p_grid%zem(i,2) * rr(i,3,m)
         ss(i,3) = p_grid%xim(i,3) * rr(i,1,m) + &
                   p_grid%etm(i,3) * rr(i,2,m) + &
                   p_grid%zem(i,3) * rr(i,3,m)

         rr(i,1,m) = de(i,1) * p_grid%yaco(i)
         rr(i,2,m) = gamm1prndtli * rr(i,1,m)
         de(i,5)   = twothirds * ( this%txx(i) + this%tyy(i) + this%tzz(i) )

         this%txx(i) = rr(i,1,m) * ( this%txx(i) + this%txx(i) - de(i,5) )
         this%tyy(i) = rr(i,1,m) * ( this%tyy(i) + this%tyy(i) - de(i,5) )
         this%tzz(i) = rr(i,1,m) * ( this%tzz(i) + this%tzz(i) - de(i,5) )

         this%txy(i) = rr(i,1,m) * ( this%txy(i) + this%hzz(i) )
         this%tyz(i) = rr(i,1,m) * ( this%tyz(i) + this%hxx(i) )
         this%tzx(i) = rr(i,1,m) * ( this%tzx(i) + this%hyy(i) )

         this%hxx(i) = rr(i,2,m) * ss(i,1) + de(i,2) * this%txx(i) + &
                       de(i,3) * this%txy(i) + de(i,4) * this%tzx(i)
         this%hyy(i) = rr(i,2,m) * ss(i,2) + de(i,2) * this%txy(i) + &
                       de(i,3) * this%tyy(i) + de(i,4) * this%tyz(i)
         this%hzz(i) = rr(i,2,m) * ss(i,3) + de(i,2) * this%tzx(i) + &
                       de(i,3) * this%tyz(i) + de(i,4) * this%tzz(i)
      end do
      !$ACC END PARALLEL
#endif
   end subroutine calc_viscous_shear_stress


!===== CALCULATION OF FLUX DERIVATIVES

   subroutine calc_fluxes(this, p_domdcomp, p_numerics, p_grid, qa, p, de, luse_acc)
      class(t_physics), intent(inout)                             :: this
      type(t_domdcomp), intent(inout)                             :: p_domdcomp
      type(t_numerics), intent(inout)                             :: p_numerics
      type(t_grid),     intent(inout)                             :: p_grid
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(in)    :: qa
      real(kind=nr), dimension(0:p_domdcomp%lmx),   intent(in)    :: p
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(inout) :: de
      logical,          intent(in), optional                      :: luse_acc

      real(kind=nr), dimension(0:p_domdcomp%lmx,3) :: ss
      real(kind=nr), dimension(0:p_domdcomp%lmx,3,5) :: rr
      integer(kind=ni) :: m,i
      logical          :: lacc

      if (present(luse_acc)) then
        lacc = luse_acc
      else
        lacc = .false.
      end if

      !$ACC PARALLEL LOOP GANG VECTOR IF (lacc)
      do i=0,p_domdcomp%lmx
         rr(i,1,1) = de(i,2) + this%umf(1)
         rr(i,2,1) = de(i,3) + this%umf(2)
         rr(i,3,1) = de(i,4) + this%umf(3)

         ss(i,1) = p_grid%xim(i,1) * rr(i,1,1) + &
                   p_grid%xim(i,2) * rr(i,2,1) + &
                   p_grid%xim(i,3) * rr(i,3,1)
         ss(i,2) = p_grid%etm(i,1) * rr(i,1,1) + &
                   p_grid%etm(i,2) * rr(i,2,1) + &
                   p_grid%etm(i,3) * rr(i,3,1)
         ss(i,3) = p_grid%zem(i,1) * rr(i,1,1) + &
                   p_grid%zem(i,2) * rr(i,2,1) + &
                   p_grid%zem(i,3) * rr(i,3,1)

         rr(i,1,1) = qa(i,1) * ss(i,1)
         rr(i,2,1) = qa(i,1) * ss(i,2)
         rr(i,3,1) = qa(i,1) * ss(i,3)


         rr(i,1,2) = qa(i,2) * ss(i,1) + p_grid%xim(i,1) * p(i)
         rr(i,2,2) = qa(i,2) * ss(i,2) + p_grid%etm(i,1) * p(i)
         rr(i,3,2) = qa(i,2) * ss(i,3) + p_grid%zem(i,1) * p(i)
#ifdef VISCOUS
         rr(i,1,2) = rr(i,1,2) - p_grid%xim(i,1) * this%txx(i) - &
                                 p_grid%xim(i,2) * this%txy(i) - &
                                 p_grid%xim(i,3) * this%tzx(i)
         rr(i,2,2) = rr(i,2,2) - p_grid%etm(i,1) * this%txx(i) - &
                                 p_grid%etm(i,2) * this%txy(i) - &
                                 p_grid%etm(i,3) * this%tzx(i)
         rr(i,3,2) = rr(i,3,2) - p_grid%zem(i,1) * this%txx(i) - &
                                 p_grid%zem(i,2) * this%txy(i) - &
                                 p_grid%zem(i,3) * this%tzx(i)
#endif


         rr(i,1,3) = qa(i,3) * ss(i,1) + p_grid%xim(i,2) * p(i)
         rr(i,2,3) = qa(i,3) * ss(i,2) + p_grid%etm(i,2) * p(i)
         rr(i,3,3) = qa(i,3) * ss(i,3) + p_grid%zem(i,2) * p(i)
#ifdef VISCOUS
         rr(i,1,3) = rr(i,1,3) - p_grid%xim(i,1) * this%txy(i) - &
                                 p_grid%xim(i,2) * this%tyy(i) - &
                                 p_grid%xim(i,3) * this%tyz(i)
         rr(i,2,3) = rr(i,2,3) - p_grid%etm(i,1) * this%txy(i) - &
                                 p_grid%etm(i,2) * this%tyy(i) - &
                                 p_grid%etm(i,3) * this%tyz(i)
         rr(i,3,3) = rr(i,3,3) - p_grid%zem(i,1) * this%txy(i) - &
                                 p_grid%zem(i,2) * this%tyy(i) - &
                                 p_grid%zem(i,3) * this%tyz(i)
#endif


         rr(i,1,4) = qa(i,4) * ss(i,1) + p_grid%xim(i,3) * p(i)
         rr(i,2,4) = qa(i,4) * ss(i,2) + p_grid%etm(i,3) * p(i)
         rr(i,3,4) = qa(i,4) * ss(i,3) + p_grid%zem(i,3) * p(i)
#ifdef VISCOUS
         rr(i,1,4) = rr(i,1,4) - p_grid%xim(i,1) * this%tzx(i) - &
                                 p_grid%xim(i,2) * this%tyz(i) - &
                                 p_grid%xim(i,3) * this%tzz(i)
         rr(i,2,4) = rr(i,2,4) - p_grid%etm(i,1) * this%tzx(i) - &
                                 p_grid%etm(i,2) * this%tyz(i) - &
                                 p_grid%etm(i,3) * this%tzz(i)
         rr(i,3,4) = rr(i,3,4) - p_grid%zem(i,1) * this%tzx(i) - &
                                 p_grid%zem(i,2) * this%tyz(i) - & 
                                 p_grid%zem(i,3) * this%tzz(i)
#endif


         de(i,5)   = qa(i,5) + p(i)
         rr(i,1,5) = de(i,5) * ss(i,1) - p(i) *      &
                   ( this%umf(1) * p_grid%xim(i,1) + &
                     this%umf(2) * p_grid%xim(i,2) + &
                     this%umf(3) * p_grid%xim(i,3) )
         rr(i,2,5) = de(i,5) * ss(i,2) - p(i) *      &
                   ( this%umf(1) * p_grid%etm(i,1) + &
                     this%umf(2) * p_grid%etm(i,2) + &
                     this%umf(3) * p_grid%etm(i,3) )
         rr(i,3,5) = de(i,5) * ss(i,3) - p(i) *      &
                   ( this%umf(1) * p_grid%zem(i,1) + & 
                     this%umf(2) * p_grid%zem(i,2) + &
                     this%umf(3) * p_grid%zem(i,3) )
#ifdef VISCOUS
         rr(i,1,5) = rr(i,1,5) - p_grid%xim(i,1) * this%hxx(i) - &
                                 p_grid%xim(i,2) * this%hyy(i) - &
                                 p_grid%xim(i,3) * this%hzz(i)
         rr(i,2,5) = rr(i,2,5) - p_grid%etm(i,1) * this%hxx(i) - &
                                 p_grid%etm(i,2) * this%hyy(i) - &
                                 p_grid%etm(i,3) * this%hzz(i)
         rr(i,3,5) = rr(i,3,5) - p_grid%zem(i,1) * this%hxx(i) - & 
                                 p_grid%zem(i,2) * this%hyy(i) - &
                                 p_grid%zem(i,3) * this%hzz(i)
#endif
      end do
      !$ACC END PARALLEL

      ! Halo exchange
      do m=1,5
         call p_numerics%mpigo(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                        p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrall, n45no, m, &
                        p_domdcomp%lxi, p_domdcomp%let, m, luse_acc = .false.)
      end do

      do m=1,5
         call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                        p_domdcomp%ijk, 1, 1, m, luse_acc = lacc)
         call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                        p_domdcomp%ijk, 2, 2, m, luse_acc = lacc)
         call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                        p_domdcomp%ijk, 3, 3, m, luse_acc = lacc)
         !$ACC PARALLEL LOOP GANG VECTOR IF (lacc)
         do i=0,p_domdcomp%lmx
           de(i,m) = rr(i,1,m) + rr(i,2,m) + rr(i,3,m)
         end do
         !$ACC END PARALLEL
      end do

   end subroutine calc_fluxes

END MODULE mo_physics
