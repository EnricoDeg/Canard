!*****
!***** SPONGE IMPLEMENTATION
!*****

MODULE mo_sponge
   use mo_kind,       ONLY : ni, nr
   use mo_parameters, ONLY : two, zero, pi, sml, one, half, hamhamm1
   use mo_gridgen,    ONLY : t_grid_geom

   IMPLICIT NONE
   PUBLIC

   real(kind=nr),    private, dimension(:), allocatable :: asz, bsz
   integer(kind=ni), private, dimension(:), allocatable :: lcsz
   integer(kind=ni), private :: lsz
   real(kind=nr),    private :: szco

   CONTAINS

   subroutine read_input_sponge
      character(16) :: ccinput
      integer(kind=ni) :: fu, rc

      namelist /nml_sponge/ szco

      open (action='read', file='input.canard', iostat=rc, newunit=fu)
      read (nml=nml_sponge, iostat=rc, unit=fu)
      close(fu)

   end subroutine read_input_sponge

   !===== SETTING UP SPONGE ZONE PARAMETERS

   subroutine spongeup(p_grid_geom, lmx, yaco, de, ss)
      type(t_grid_geom), intent(in) :: p_grid_geom
      integer(kind=ni), INTENT(IN) :: lmx
      real(kind=nr), dimension(0:lmx),   intent(in)    :: yaco
      real(kind=nr), dimension(0:lmx,5), intent(inout) :: de
      real(kind=nr), dimension(0:lmx,3), intent(in)    :: ss
      integer(kind=ni) :: l, ll
      real(kind=nr) :: ra0, ra1, ra2, ra3
      real(kind=nr) :: tmpa, tmpb

      ll   = -1
      ra2  = zero
      tmpa = pi / p_grid_geom%szth0
      tmpb = pi / p_grid_geom%szth1
      do l=0,lmx
         ra3 = ra2 * ss(l,2)
         ra0 = tmpa * ( ss(l,1) - ( ra3 - p_grid_geom%doml0 + p_grid_geom%szth0 ) )
         ra1 = tmpb * ( ra3 + p_grid_geom%doml1 - p_grid_geom%szth1 - ss(l,1) )
         de(l,1) = szco * half * ( two + cos( max( min( ra0, pi ), zero ) ) + &
                                         cos( max( min( ra1, pi ), zero ) ) )
         de(l,2) = szco * half * ( one + cos( max( min( ra0, pi ), zero ) ) )
         if ( de(l,1) > sml ) then
            ll = ll + 1
            de(ll,5) = l + sml
         end if
      end do
      lsz = ll
      if ( lsz /= -1 ) then
         allocate(lcsz(0:lsz),asz(0:lsz),bsz(0:lsz))
         do ll=0,lsz
            l=de(ll,5)
            lcsz(ll)=l
            asz(ll)=de(l,1)/yaco(l)
            bsz(ll)=de(l,2)/yaco(l)
         end do
      end if

   end subroutine spongeup

 !===== SPONGE IMPLEMENTATION

   subroutine spongego(lmx, qa, de, luse_acc)
      integer(kind=ni), intent(in) :: lmx
      real(kind=nr), dimension(0:lmx,5), intent(in) :: qa
      real(kind=nr), dimension(0:lmx,5), intent(inout) :: de
      logical, intent(in), optional :: luse_acc
      integer(kind=ni) :: l, ll
      logical :: lacc

      IF (PRESENT(luse_acc)) THEN
         lacc = luse_acc
      ELSE
         lacc = .false.
      END IF

      !$ACC PARALLEL LOOP GANG VECTOR IF (lacc)
      do ll=0,lsz
         l=lcsz(ll)
         de(l,1)=de(l,1)+asz(ll)*(qa(l,1)-one)
         de(l,2:4)=de(l,2:4)+bsz(ll)*(qa(l,2:4)-zero)
         de(l,5)=de(l,5)+asz(ll)*(qa(l,5)-hamhamm1)
      end do
      !$ACC END PARALLEL

   end subroutine spongego

END MODULE mo_sponge
