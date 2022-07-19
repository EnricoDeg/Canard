!*****
!***** PHYSICS MODULE
!*****

MODULE mo_physics
   use mo_kind,       ONLY : nr, ni
   use mo_parameters, ONLY : sml, zero, one, pi, hamm1, hamhamm1, half, gam,      &
                           & gamm1, n45no, nrall, gamm1prndtli, nrone, twothirds
   use mo_grid,       ONLY : t_grid
   use mo_domdcomp,   ONLY : t_domdcomp
   use mo_numerics,   ONLY : t_numerics
   implicit none
   public

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

   subroutine allocate_physics_memory(lmx)
      integer(kind=ni), intent(IN) :: lmx

#ifdef VISCOUS
      allocate(txx(0:lmx), tyy(0:lmx), tzz(0:lmx))
      allocate(txy(0:lmx), tyz(0:lmx), tzx(0:lmx))
      allocate(hxx(0:lmx), hyy(0:lmx), hzz(0:lmx))
#endif

   end subroutine allocate_physics_memory



   subroutine deallocate_physics_memory

#ifdef VISCOUS
      deallocate(txx, tyy, tzz, txy, tyz, tzx, hxx, hyy, hzz)
#endif

   end subroutine deallocate_physics_memory

!===== INITIALIZE PHYSICS

   subroutine read_input_physics
      character(16) :: cinput

      open(9,file='input.physics',status='old')
      read(9,*) cinput,reoo
      read(9,*) cinput,tempoo
      read(9,*) cinput,amach1
      read(9,*) cinput,amach2
      read(9,*) cinput,amach3
      read(9,*) cinput,timf
      read(9,*) cinput,nsmf
      close(9)

      amachoo = sqrt( amach1 * amach1 + amach2 * amach2 + amach3 * amach3 )
      if ( amachoo > sml ) then
         reoo = reoo / amachoo
      end if
      srefoo    = 111 / tempoo
      srefp1dre = ( srefoo + one ) / reoo
      sqrtrema  = sqrt( reoo * amachoo )
      sqrtremai = one / max( sqrtrema, sml )
      uoo(:)    = (/ amach1, amach2, amach3 /)

   end subroutine read_input_physics

!===== INITIAL CONDITIONS

   subroutine initialo(lmx, qa, ss)
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

   subroutine movef(dtko, dtk, timo)
      real(kind=nr), intent(in) :: dtko, dtk
      real(kind=nr), intent(in) :: timo
      real(kind=nr) :: ra0, ra1, ra2, dfdt, fctr, progmf

      if ( nsmf == 0 ) then
         ra0 = pi / timf
         ra1 = ra0 * min( timo, timf )
         ra2 = ra0 * min( timo+dtko, timf )

         fctr   = one - cos(ra1)
         dfdt   = ra0 * sin(ra2)
         progmf = half * ( fctr + dtk * dfdt )
         umf(:) = progmf * uoo(:)

         fctr      = sin(ra1)
         dfdt      = ra0 * cos(ra2)
         progmf    = half * ra0 * ( fctr + dtk * dfdt )
         dudtmf(:) = progmf * uoo(:)
      else
         umf(:)    = uoo(:)
         dudtmf(:) = zero
      end if

   end subroutine movef


!===== VISCOUS SHEAR STRESSES & HEAT FLUXES

   subroutine calc_viscous_shear_stress(p_domdcomp, p_numerics, p_grid, de, ssk)
      type(t_domdcomp), intent(IN)    :: p_domdcomp
      type(t_numerics), intent(inout) :: p_numerics
      type(t_grid),     intent(in)    :: p_grid
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(inout) :: de
      real(kind=nr), dimension(0:p_domdcomp%lmx), intent(in)    :: ssk
      real(kind=nr), dimension(0:p_domdcomp%lmx,3) :: ss
      real(kind=nr), dimension(0:p_domdcomp%lmx,3) :: rr

      integer(kind=ni) :: m
#ifdef VISCOUS
      de(:,1) = ssk(:)

      rr(:,1) = de(:,2)
      m = 2
      call p_numerics%mpigo(rr(:,1), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                          p_domdcomp%mcd, p_domdcomp%nbsize, 0, n45no, m, &
                          p_domdcomp%lxi, p_domdcomp%let)
      call p_numerics%deriv(rr(:,1), rr(:,3), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                                   p_domdcomp%lze, p_domdcomp%ijk, 3, m)
      call p_numerics%deriv(rr(:,1), rr(:,2), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                                   p_domdcomp%lze, p_domdcomp%ijk, 2, m)
      call p_numerics%deriv(rr(:,1), rr(:,1), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                                   p_domdcomp%lze, p_domdcomp%ijk, 1, m)
      txx(:) = p_grid%xim(:,1) * rr(:,1) + p_grid%etm(:,1) * rr(:,2) + p_grid%zem(:,1) * rr(:,3)
      hzz(:) = p_grid%xim(:,2) * rr(:,1) + p_grid%etm(:,2) * rr(:,2) + p_grid%zem(:,2) * rr(:,3)
      tzx(:) = p_grid%xim(:,3) * rr(:,1) + p_grid%etm(:,3) * rr(:,2) + p_grid%zem(:,3) * rr(:,3)

      rr(:,1) = de(:,3)
      m = 3
      call p_numerics%mpigo(rr, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                     p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrone, n45no, m, &
                     p_domdcomp%lxi, p_domdcomp%let)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 1, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, & 
                     p_domdcomp%ijk, 2, 1, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m)
      txy(:) = p_grid%xim(:,1) * rr(:,1) + p_grid%etm(:,1) * rr(:,2) + p_grid%zem(:,1) * rr(:,3)
      tyy(:) = p_grid%xim(:,2) * rr(:,1) + p_grid%etm(:,2) * rr(:,2) + p_grid%zem(:,2) * rr(:,3)
      hxx(:) = p_grid%xim(:,3) * rr(:,1) + p_grid%etm(:,3) * rr(:,2) + p_grid%zem(:,3) * rr(:,3)

      rr(:,1) = de(:,4)
      m = 4
      call p_numerics%mpigo(rr, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                     p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrone, n45no, m, &
                     p_domdcomp%lxi, p_domdcomp%let)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 1, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 1, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m)
      hyy(:) = p_grid%xim(:,1) * rr(:,1) + p_grid%etm(:,1) * rr(:,2) + p_grid%zem(:,1) * rr(:,3)
      tyz(:) = p_grid%xim(:,2) * rr(:,1) + p_grid%etm(:,2) * rr(:,2) + p_grid%zem(:,2) * rr(:,3)
      tzz(:) = p_grid%xim(:,3) * rr(:,1) + p_grid%etm(:,3) * rr(:,2) + p_grid%zem(:,3) * rr(:,3)

      rr(:,1) = de(:,5)
      m = 5
      call p_numerics%mpigo(rr, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                     p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrone, n45no, m, & 
                     p_domdcomp%lxi, p_domdcomp%let)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 1, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 1, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m)
      ss(:,1) = p_grid%xim(:,1) * rr(:,1) + p_grid%etm(:,1) * rr(:,2) + p_grid%zem(:,1) * rr(:,3)
      ss(:,2) = p_grid%xim(:,2) * rr(:,1) + p_grid%etm(:,2) * rr(:,2) + p_grid%zem(:,2) * rr(:,3)
      ss(:,3) = p_grid%xim(:,3) * rr(:,1) + p_grid%etm(:,3) * rr(:,2) + p_grid%zem(:,3) * rr(:,3)

      rr(:,1) = de(:,1) * p_grid%yaco(:)
      rr(:,2) = gamm1prndtli * rr(:,1)
      de(:,5) = twothirds * ( txx(:) + tyy(:) + tzz(:) )

      txx(:) = rr(:,1) * ( txx(:) + txx(:) - de(:,5) )
      tyy(:) = rr(:,1) * ( tyy(:) + tyy(:) - de(:,5) )
      tzz(:) = rr(:,1) * ( tzz(:) + tzz(:) - de(:,5) )
      txy(:) = rr(:,1) * ( txy(:) + hzz(:) )
      tyz(:) = rr(:,1) * ( tyz(:) + hxx(:) )
      tzx(:) = rr(:,1) * ( tzx(:) + hyy(:) )
      hxx(:) = rr(:,2) * ss(:,1) + de(:,2) * txx(:) + de(:,3) * txy(:) + de(:,4) * tzx(:)
      hyy(:) = rr(:,2) * ss(:,2) + de(:,2) * txy(:) + de(:,3) * tyy(:) + de(:,4) * tyz(:)
      hzz(:) = rr(:,2) * ss(:,3) + de(:,2) * tzx(:) + de(:,3) * tyz(:) + de(:,4) * tzz(:)
#endif
   end subroutine calc_viscous_shear_stress


!===== CALCULATION OF FLUX DERIVATIVES

   subroutine calc_fluxes(p_domdcomp, p_numerics, p_grid, qa, p, de)
      type(t_domdcomp), intent(IN)    :: p_domdcomp
      type(t_numerics), intent(inout) :: p_numerics
      type(t_grid),     intent(in)    :: p_grid
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(in) :: qa
      real(kind=nr), dimension(0:p_domdcomp%lmx), intent(in) :: p
      real(kind=nr), dimension(0:p_domdcomp%lmx,5), intent(inout) :: de

      real(kind=nr), dimension(0:p_domdcomp%lmx,3) :: ss
      real(kind=nr), dimension(0:p_domdcomp%lmx,3) :: rr
      integer(kind=ni) :: m

      rr(:,1) = de(:,2) + umf(1)
      rr(:,2) = de(:,3) + umf(2)
      rr(:,3) = de(:,4) + umf(3)
      ss(:,1) = p_grid%xim(:,1) * rr(:,1) + p_grid%xim(:,2) * rr(:,2) + p_grid%xim(:,3) * rr(:,3)
      ss(:,2) = p_grid%etm(:,1) * rr(:,1) + p_grid%etm(:,2) * rr(:,2) + p_grid%etm(:,3) * rr(:,3)
      ss(:,3) = p_grid%zem(:,1) * rr(:,1) + p_grid%zem(:,2) * rr(:,2) + p_grid%zem(:,3) * rr(:,3)

      rr(:,1) = qa(:,1) * ss(:,1)
      rr(:,2) = qa(:,1) * ss(:,2)
      rr(:,3) = qa(:,1) * ss(:,3)
      m = 1
      call p_numerics%mpigo(rr, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                     p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrall, n45no, m, &
                     p_domdcomp%lxi, p_domdcomp%let)
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
      rr(:,1) = rr(:,1) - p_grid%xim(:,1) * txx(:) - p_grid%xim(:,2) * txy(:) - p_grid%xim(:,3) * tzx(:)
      rr(:,2) = rr(:,2) - p_grid%etm(:,1) * txx(:) - p_grid%etm(:,2) * txy(:) - p_grid%etm(:,3) * tzx(:)
      rr(:,3) = rr(:,3) - p_grid%zem(:,1) * txx(:) - p_grid%zem(:,2) * txy(:) - p_grid%zem(:,3) * tzx(:)
#endif
      m = 2
      call p_numerics%mpigo(rr, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                     p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrall, n45no, m, &
                     p_domdcomp%lxi, p_domdcomp%let)
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
      rr(:,1) = rr(:,1) - p_grid%xim(:,1) * txy(:) - p_grid%xim(:,2) * tyy(:) - p_grid%xim(:,3) * tyz(:)
      rr(:,2) = rr(:,2) - p_grid%etm(:,1) * txy(:) - p_grid%etm(:,2) * tyy(:) - p_grid%etm(:,3) * tyz(:)
      rr(:,3) = rr(:,3) - p_grid%zem(:,1) * txy(:) - p_grid%zem(:,2) * tyy(:) - p_grid%zem(:,3) * tyz(:)
#endif
      m = 3
      call p_numerics%mpigo(rr, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                     p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrall, n45no, m, &
                     p_domdcomp%lxi, p_domdcomp%let)
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
      rr(:,1) = rr(:,1) - p_grid%xim(:,1) * tzx(:) - p_grid%xim(:,2) * tyz(:) - p_grid%xim(:,3) * tzz(:)
      rr(:,2) = rr(:,2) - p_grid%etm(:,1) * tzx(:) - p_grid%etm(:,2) * tyz(:) - p_grid%etm(:,3) * tzz(:)
      rr(:,3) = rr(:,3) - p_grid%zem(:,1) * tzx(:) - p_grid%zem(:,2) * tyz(:) - p_grid%zem(:,3) * tzz(:)
#endif
      m = 4
      call p_numerics%mpigo(rr, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                     p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrall, n45no, m, &
                     p_domdcomp%lxi, p_domdcomp%let)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 2, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 3, m)
      de(:,m) = rr(:,1) + rr(:,2) + rr(:,3)

      de(:,5) = qa(:,5) + p(:)
      rr(:,1) = de(:,5) * ss(:,1) - p(:) * &
              ( umf(1) * p_grid%xim(:,1) + umf(2) * p_grid%xim(:,2) + umf(3) * p_grid%xim(:,3) )
      rr(:,2) = de(:,5) * ss(:,2) - p(:) * &
              ( umf(1) * p_grid%etm(:,1) + umf(2) * p_grid%etm(:,2) + umf(3) * p_grid%etm(:,3) )
      rr(:,3) = de(:,5) * ss(:,3) - p(:) * &
              ( umf(1) * p_grid%zem(:,1) + umf(2) * p_grid%zem(:,2) + umf(3) * p_grid%zem(:,3) )
#ifdef VISCOUS
      rr(:,1) = rr(:,1) - p_grid%xim(:,1) * hxx(:) - p_grid%xim(:,2) * hyy(:) - p_grid%xim(:,3) * hzz(:)
      rr(:,2) = rr(:,2) - p_grid%etm(:,1) * hxx(:) - p_grid%etm(:,2) * hyy(:) - p_grid%etm(:,3) * hzz(:)
      rr(:,3) = rr(:,3) - p_grid%zem(:,1) * hxx(:) - p_grid%zem(:,2) * hyy(:) - p_grid%zem(:,3) * hzz(:)
#endif
      m = 5
      call p_numerics%mpigo(rr, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc,        &
                     p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrall, n45no, m, &
                     p_domdcomp%lxi, p_domdcomp%let)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 1, 1, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 2, 2, m)
      call p_numerics%deriv(rr, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, p_domdcomp%lze, &
                     p_domdcomp%ijk, 3, 3, m)
      de(:,m) = rr(:,1) + rr(:,2) + rr(:,3)

   end subroutine calc_fluxes

END MODULE mo_physics
