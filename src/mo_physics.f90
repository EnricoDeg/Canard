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

   subroutine calc_viscous_shear_stress(this, p_domdcomp, p_numerics, p_grid, de, ssk)
      class(t_physics), INTENT(INOUT) :: this
      type(t_domdcomp), intent(IN)    :: p_domdcomp
      type(t_numerics), intent(inout) :: p_numerics
      type(t_grid),     intent(in)    :: p_grid
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(inout) :: de
      real(kind=nr), dimension(0:p_domdcomp%lmx), intent(in)    :: ssk
      real(kind=nr), dimension(0:p_domdcomp%lmx,3) :: ss
      real(kind=nr), dimension(0:p_domdcomp%lmx,3,2:5) :: rr

      integer(kind=ni) :: m
#ifdef VISCOUS
      de(:,1) = ssk(:)

      ! Halo exchange
      do m=2,5
         rr(:,1,m) = de(:,m)
         call p_numerics%mpigo(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                               p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrone, n45no, m, &
                               p_domdcomp%lxi, p_domdcomp%let, m)
      end do

      m = 2
      call p_numerics%deriv(rr(:,1,m), rr(:,3,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                                   p_domdcomp%lze, p_domdcomp%ijk, 3, m)
      call p_numerics%deriv(rr(:,1,m), rr(:,2,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                                   p_domdcomp%lze, p_domdcomp%ijk, 2, m)
      call p_numerics%deriv(rr(:,1,m), rr(:,1,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                                   p_domdcomp%lze, p_domdcomp%ijk, 1, m)
      this%txx(:) = p_grid%xim(:,1) * rr(:,1,m) + p_grid%etm(:,1) * rr(:,2,m) + p_grid%zem(:,1) * rr(:,3,m)
      this%hzz(:) = p_grid%xim(:,2) * rr(:,1,m) + p_grid%etm(:,2) * rr(:,2,m) + p_grid%zem(:,2) * rr(:,3,m)
      this%tzx(:) = p_grid%xim(:,3) * rr(:,1,m) + p_grid%etm(:,3) * rr(:,2,m) + p_grid%zem(:,3) * rr(:,3,m)

      m = 3
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 1, m)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, & 
                     p_domdcomp%ijk, 2, 1, m)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m)
      this%txy(:) = p_grid%xim(:,1) * rr(:,1,m) + p_grid%etm(:,1) * rr(:,2,m) + p_grid%zem(:,1) * rr(:,3,m)
      this%tyy(:) = p_grid%xim(:,2) * rr(:,1,m) + p_grid%etm(:,2) * rr(:,2,m) + p_grid%zem(:,2) * rr(:,3,m)
      this%hxx(:) = p_grid%xim(:,3) * rr(:,1,m) + p_grid%etm(:,3) * rr(:,2,m) + p_grid%zem(:,3) * rr(:,3,m)

      m = 4
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 1, m)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 1, m)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m)
      this%hyy(:) = p_grid%xim(:,1) * rr(:,1,m) + p_grid%etm(:,1) * rr(:,2,m) + p_grid%zem(:,1) * rr(:,3,m)
      this%tyz(:) = p_grid%xim(:,2) * rr(:,1,m) + p_grid%etm(:,2) * rr(:,2,m) + p_grid%zem(:,2) * rr(:,3,m)
      this%tzz(:) = p_grid%xim(:,3) * rr(:,1,m) + p_grid%etm(:,3) * rr(:,2,m) + p_grid%zem(:,3) * rr(:,3,m)

      m = 5
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 1, m)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 1, m)
      call p_numerics%deriv(rr(:,:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m)
      ss(:,1) = p_grid%xim(:,1) * rr(:,1,m) + p_grid%etm(:,1) * rr(:,2,m) + p_grid%zem(:,1) * rr(:,3,m)
      ss(:,2) = p_grid%xim(:,2) * rr(:,1,m) + p_grid%etm(:,2) * rr(:,2,m) + p_grid%zem(:,2) * rr(:,3,m)
      ss(:,3) = p_grid%xim(:,3) * rr(:,1,m) + p_grid%etm(:,3) * rr(:,2,m) + p_grid%zem(:,3) * rr(:,3,m)

      rr(:,1,m) = de(:,1) * p_grid%yaco(:)
      rr(:,2,m) = gamm1prndtli * rr(:,1,m)
      de(:,5) = twothirds * ( this%txx(:) + this%tyy(:) + this%tzz(:) )

      this%txx(:) = rr(:,1,m) * ( this%txx(:) + this%txx(:) - de(:,5) )
      this%tyy(:) = rr(:,1,m) * ( this%tyy(:) + this%tyy(:) - de(:,5) )
      this%tzz(:) = rr(:,1,m) * ( this%tzz(:) + this%tzz(:) - de(:,5) )

      this%txy(:) = rr(:,1,m) * ( this%txy(:) + this%hzz(:) )
      this%tyz(:) = rr(:,1,m) * ( this%tyz(:) + this%hxx(:) )
      this%tzx(:) = rr(:,1,m) * ( this%tzx(:) + this%hyy(:) )

      this%hxx(:) = rr(:,2,m) * ss(:,1) + de(:,2) * this%txx(:) + &
                    de(:,3) * this%txy(:) + de(:,4) * this%tzx(:)
      this%hyy(:) = rr(:,2,m) * ss(:,2) + de(:,2) * this%txy(:) + &
                    de(:,3) * this%tyy(:) + de(:,4) * this%tyz(:)
      this%hzz(:) = rr(:,2,m) * ss(:,3) + de(:,2) * this%tzx(:) + &
                    de(:,3) * this%tyz(:) + de(:,4) * this%tzz(:)
#endif
   end subroutine calc_viscous_shear_stress


!===== CALCULATION OF FLUX DERIVATIVES

   subroutine calc_fluxes(this, p_domdcomp, p_numerics, p_grid, qa, p, de)
      class(t_physics), INTENT(INOUT) :: this
      type(t_domdcomp), intent(IN)    :: p_domdcomp
      type(t_numerics), intent(inout) :: p_numerics
      type(t_grid),     intent(in)    :: p_grid
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(in) :: qa
      real(kind=nr), dimension(0:p_domdcomp%lmx), intent(in) :: p
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(inout) :: de

      real(kind=nr), dimension(0:p_domdcomp%lmx,3) :: ss
      real(kind=nr), dimension(0:p_domdcomp%lmx,3) :: rr
      integer(kind=ni) :: m

      rr(:,1) = de(:,2) + this%umf(1)
      rr(:,2) = de(:,3) + this%umf(2)
      rr(:,3) = de(:,4) + this%umf(3)

      ss(:,1) = p_grid%xim(:,1) * rr(:,1) + p_grid%xim(:,2) * rr(:,2) + p_grid%xim(:,3) * rr(:,3)
      ss(:,2) = p_grid%etm(:,1) * rr(:,1) + p_grid%etm(:,2) * rr(:,2) + p_grid%etm(:,3) * rr(:,3)
      ss(:,3) = p_grid%zem(:,1) * rr(:,1) + p_grid%zem(:,2) * rr(:,2) + p_grid%zem(:,3) * rr(:,3)

      rr(:,1) = qa(:,1) * ss(:,1)
      rr(:,2) = qa(:,1) * ss(:,2)
      rr(:,3) = qa(:,1) * ss(:,3)
      m = 1
      call p_numerics%mpigo(rr, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                     p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrall, n45no, m, &
                     p_domdcomp%lxi, p_domdcomp%let, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 2, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 3, m)
      de(:,m) = rr(:,1) + rr(:,2) + rr(:,3)

      rr(:,1) = qa(:,2) * ss(:,1) + p_grid%xim(:,1) * p(:)
      rr(:,2) = qa(:,2) * ss(:,2) + p_grid%etm(:,1) * p(:)
      rr(:,3) = qa(:,2) * ss(:,3) + p_grid%zem(:,1) * p(:)
#ifdef VISCOUS
      rr(:,1) = rr(:,1) - p_grid%xim(:,1) * this%txx(:) - &
                          p_grid%xim(:,2) * this%txy(:) - p_grid%xim(:,3) * this%tzx(:)
      rr(:,2) = rr(:,2) - p_grid%etm(:,1) * this%txx(:) - &
                          p_grid%etm(:,2) * this%txy(:) - p_grid%etm(:,3) * this%tzx(:)
      rr(:,3) = rr(:,3) - p_grid%zem(:,1) * this%txx(:) - &
                          p_grid%zem(:,2) * this%txy(:) - p_grid%zem(:,3) * this%tzx(:)
#endif
      m = 2
      call p_numerics%mpigo(rr, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                     p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrall, n45no, m, &
                     p_domdcomp%lxi, p_domdcomp%let, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 2, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 3, m)
      de(:,m) = rr(:,1) + rr(:,2) + rr(:,3)

      rr(:,1) = qa(:,3) * ss(:,1) + p_grid%xim(:,2) * p(:)
      rr(:,2) = qa(:,3) * ss(:,2) + p_grid%etm(:,2) * p(:)
      rr(:,3) = qa(:,3) * ss(:,3) + p_grid%zem(:,2) * p(:)
#ifdef VISCOUS
      rr(:,1) = rr(:,1) - p_grid%xim(:,1) * this%txy(:) - &
                          p_grid%xim(:,2) * this%tyy(:) - p_grid%xim(:,3) * this%tyz(:)
      rr(:,2) = rr(:,2) - p_grid%etm(:,1) * this%txy(:) - &
                          p_grid%etm(:,2) * this%tyy(:) - p_grid%etm(:,3) * this%tyz(:)
      rr(:,3) = rr(:,3) - p_grid%zem(:,1) * this%txy(:) - &
                          p_grid%zem(:,2) * this%tyy(:) - p_grid%zem(:,3) * this%tyz(:)
#endif
      m = 3
      call p_numerics%mpigo(rr, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                     p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrall, n45no, m, &
                     p_domdcomp%lxi, p_domdcomp%let, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 2, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 3, m)
      de(:,m) = rr(:,1) + rr(:,2) + rr(:,3)

      rr(:,1) = qa(:,4) * ss(:,1) + p_grid%xim(:,3) * p(:)
      rr(:,2) = qa(:,4) * ss(:,2) + p_grid%etm(:,3) * p(:)
      rr(:,3) = qa(:,4) * ss(:,3) + p_grid%zem(:,3) * p(:)
#ifdef VISCOUS
      rr(:,1) = rr(:,1) - p_grid%xim(:,1) * this%tzx(:) - &
                          p_grid%xim(:,2) * this%tyz(:) - p_grid%xim(:,3) * this%tzz(:)
      rr(:,2) = rr(:,2) - p_grid%etm(:,1) * this%tzx(:) - &
                          p_grid%etm(:,2) * this%tyz(:) - p_grid%etm(:,3) * this%tzz(:)
      rr(:,3) = rr(:,3) - p_grid%zem(:,1) * this%tzx(:) - &
                          p_grid%zem(:,2) * this%tyz(:) - p_grid%zem(:,3) * this%tzz(:)
#endif
      m = 4
      call p_numerics%mpigo(rr, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                     p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrall, n45no, m, &
                     p_domdcomp%lxi, p_domdcomp%let, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 2, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 3, m)
      de(:,m) = rr(:,1) + rr(:,2) + rr(:,3)

      de(:,5) = qa(:,5) + p(:)
      rr(:,1) = de(:,5) * ss(:,1) - p(:) * &
              ( this%umf(1) * p_grid%xim(:,1) + this%umf(2) * p_grid%xim(:,2) + this%umf(3) * p_grid%xim(:,3) )
      rr(:,2) = de(:,5) * ss(:,2) - p(:) * &
              ( this%umf(1) * p_grid%etm(:,1) + this%umf(2) * p_grid%etm(:,2) + this%umf(3) * p_grid%etm(:,3) )
      rr(:,3) = de(:,5) * ss(:,3) - p(:) * &
              ( this%umf(1) * p_grid%zem(:,1) + this%umf(2) * p_grid%zem(:,2) + this%umf(3) * p_grid%zem(:,3) )
#ifdef VISCOUS
      rr(:,1) = rr(:,1) - p_grid%xim(:,1) * this%hxx(:) - p_grid%xim(:,2) * this%hyy(:) - p_grid%xim(:,3) * this%hzz(:)
      rr(:,2) = rr(:,2) - p_grid%etm(:,1) * this%hxx(:) - p_grid%etm(:,2) * this%hyy(:) - p_grid%etm(:,3) * this%hzz(:)
      rr(:,3) = rr(:,3) - p_grid%zem(:,1) * this%hxx(:) - p_grid%zem(:,2) * this%hyy(:) - p_grid%zem(:,3) * this%hzz(:)
#endif
      m = 5
      call p_numerics%mpigo(rr, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                     p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrall, n45no, m, &
                     p_domdcomp%lxi, p_domdcomp%let, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 2, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 3, m)
      de(:,m) = rr(:,1) + rr(:,2) + rr(:,3)

   end subroutine calc_fluxes

END MODULE mo_physics
