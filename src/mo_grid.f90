!*****
!***** GRID MODULE
!*****

MODULE mo_grid
   use mo_kind,       ONLY : nr, ni
   use mo_parameters, ONLY : n45go, n45no, nrone, one, three
   use mo_domdcomp,   ONLY : t_domdcomp
   use mo_numerics,   ONLY : t_numerics
   use mo_utils,      ONLY : indx3
   
   IMPLICIT NONE
   PRIVATE

   type , public :: t_grid
      ! public variables
      real(kind=nr),    public, dimension(:),     allocatable :: yaco
      real(kind=nr),    public, dimension(:,:),   allocatable :: xim, etm, zem
      real(kind=nr),    public, dimension(:,:,:), pointer     :: cm1, cm2, cm3

   contains
      ! public subroutines
      procedure, public :: allocate     => allocate_grid_memory
      procedure, public :: deallocate   => deallocate_grid_memory
      procedure, public :: grid_metrics => calc_grid_metrics

   end type t_grid

   CONTAINS

   SUBROUTINE allocate_grid_memory(this, p_domdcomp)
      class(t_grid), INTENT(INOUT) :: this
      type(t_domdcomp), intent(IN) :: p_domdcomp
      integer(kind=ni) :: ii, jj, kk

      ii = p_domdcomp%nbsize(1)-1
      jj = p_domdcomp%nbsize(2)-1
      kk = p_domdcomp%nbsize(3)-1
      allocate(this%yaco(0:p_domdcomp%lmx))
      allocate(this%cm1(0:ii,3,0:1), this%cm2(0:jj,3,0:1), this%cm3(0:kk,3,0:1))
      allocate(this%xim(0:p_domdcomp%lmx,3), this%etm(0:p_domdcomp%lmx,3), this%zem(0:p_domdcomp%lmx,3))

   END SUBROUTINE allocate_grid_memory

   SUBROUTINE deallocate_grid_memory(this)
      class(t_grid), INTENT(INOUT) :: this

      deallocate(this%xim, this%etm, this%zem, this%yaco)

   END SUBROUTINE deallocate_grid_memory

   SUBROUTINE calc_grid_metrics(this, p_domdcomp, p_numerics, ssk)
      class(t_grid), INTENT(INOUT) :: this
      type(t_domdcomp), intent(IN) :: p_domdcomp
      type(t_numerics), intent(inout) :: p_numerics
      real(kind=nr), intent(in) :: ssk(0:p_domdcomp%lmx,3)
      integer(kind=ni) :: m, nn, ip, i, j, k, l, jk, kp
      real(kind=nr)    :: fctr
      real(kind=nr), dimension(:,:), allocatable :: dek, qok, qak, rrk
      real(kind=nr), dimension(3) :: rv

      allocate(dek(0:p_domdcomp%lmx,5), qok(0:p_domdcomp%lmx,5), &
               qak(0:p_domdcomp%lmx,5), rrk(0:p_domdcomp%lmx,3))

      rrk(:,1) = ssk(:,1)
      m = 1
      call p_numerics%mpigo(rrk, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                      p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrone, n45go, m, p_domdcomp%lxi, p_domdcomp%let,m)
      call p_numerics%deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 3, 1, m)
      call p_numerics%deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 2, 1, m)
      call p_numerics%deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 1, 1, m)
      qok(:,1) = rrk(:,1)
      qok(:,2) = rrk(:,2)
      qok(:,3) = rrk(:,3)

      rrk(:,1) = ssk(:,2)
      m = 2
      call p_numerics%mpigo(rrk, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                      p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrone, n45go, m, p_domdcomp%lxi, p_domdcomp%let,m)
      call p_numerics%deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 3, 1, m)
      call p_numerics%deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 2, 1, m)
      call p_numerics%deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 1, 1, m)
      qak(:,1) = rrk(:,1)
      qak(:,2) = rrk(:,2)
      qak(:,3) = rrk(:,3)

      rrk(:,1) = ssk(:,3)
      m = 3
      call p_numerics%mpigo(rrk, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                      p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrone, n45go, m, p_domdcomp%lxi, p_domdcomp%let,m)
      call p_numerics%deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 3, 1, m)
      call p_numerics%deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 2, 1, m)
      call p_numerics%deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 1, 1, m)
      dek(:,1) = rrk(:,1)
      dek(:,2) = rrk(:,2)
      dek(:,3) = rrk(:,3)

      this%xim(:,1) = qak(:,2) * dek(:,3) - dek(:,2) * qak(:,3)
      this%xim(:,2) = dek(:,2) * qok(:,3) - qok(:,2) * dek(:,3)
      this%xim(:,3) = qok(:,2) * qak(:,3) - qak(:,2) * qok(:,3)
      this%etm(:,1) = qak(:,3) * dek(:,1) - dek(:,3) * qak(:,1)
      this%etm(:,2) = dek(:,3) * qok(:,1) - qok(:,3) * dek(:,1)
      this%etm(:,3) = qok(:,3) * qak(:,1) - qak(:,3) * qok(:,1)
      this%zem(:,1) = qak(:,1) * dek(:,2) - dek(:,1) * qak(:,2)
      this%zem(:,2) = dek(:,1) * qok(:,2) - qok(:,1) * dek(:,2)
      this%zem(:,3) = qok(:,1) * qak(:,2) - qak(:,1) * qok(:,2)

      do m=1,3
         call p_numerics%mpigo(this%xim(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+1, p_domdcomp%lxi, p_domdcomp%let,m)
         call p_numerics%filte(this%xim(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, 1, m)
         call p_numerics%mpigo(this%xim(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+2, p_domdcomp%lxi, p_domdcomp%let,m)
         call p_numerics%filte(this%xim(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, 2, m)
         call p_numerics%mpigo(this%xim(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+3, p_domdcomp%lxi, p_domdcomp%let,m)
         call p_numerics%filte(this%xim(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, 3, m)

         call p_numerics%mpigo(this%etm(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+4, p_domdcomp%lxi, p_domdcomp%let,m)
         call p_numerics%filte(this%etm(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, 1, m)
         call p_numerics%mpigo(this%etm(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+5, p_domdcomp%lxi, p_domdcomp%let,m)
         call p_numerics%filte(this%etm(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, 2, m)
         call p_numerics%mpigo(this%etm(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+6, p_domdcomp%lxi, p_domdcomp%let,m)
         call p_numerics%filte(this%etm(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, 3, m)

         call p_numerics%mpigo(this%zem(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+7, p_domdcomp%lxi, p_domdcomp%let,m)
         call p_numerics%filte(this%zem(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, & 
                              p_domdcomp%lze, p_domdcomp%ijk, 1, m)
         call p_numerics%mpigo(this%zem(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+8, p_domdcomp%lxi, p_domdcomp%let,m)
         call p_numerics%filte(this%zem(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, 2, m)
         call p_numerics%mpigo(this%zem(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+9, p_domdcomp%lxi, p_domdcomp%let,m)
         call p_numerics%filte(this%zem(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, 3, m)
      end do
      this%yaco(:) = three / ( qok(:,1) * this%xim(:,1) + qok(:,2) * this%etm(:,1) + qok(:,3) * this%zem(:,1) &
                             + qak(:,1) * this%xim(:,2) + qak(:,2) * this%etm(:,2) + qak(:,3) * this%zem(:,2) &
                             + dek(:,1) * this%xim(:,3) + dek(:,2) * this%etm(:,3) + dek(:,3) * this%zem(:,3) )

      do nn=1,3
         do ip=0,1
            i = ip * p_domdcomp%ijk(1,nn)
            do k=0,p_domdcomp%ijk(3,nn)
               kp = k * ( p_domdcomp%ijk(2,nn) + 1 )
               do j=0,p_domdcomp%ijk(2,nn)
                  jk = kp + j
                  l = indx3(i, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                  select case(nn)
                  case(1)
                     rv(:) = this%yaco(l) * this%xim(l,:)
                     fctr  = one / sqrt(rv(1) * rv(1) + rv(2) * rv(2) + rv(3) * rv(3))
                     this%cm1(jk,:,ip) = fctr * rv(:)
                  case(2)
                     rv(:) = this%yaco(l) * this%etm(l,:)
                     fctr  = one / sqrt(rv(1) * rv(1) + rv(2) * rv(2) + rv(3) * rv(3))
                     this%cm2(jk,:,ip) = fctr * rv(:)
                  case(3)
                     rv(:) = this%yaco(l) * this%zem(l,:)
                     fctr  = one / sqrt(rv(1) * rv(1) + rv(2) * rv(2) + rv(3) * rv(3))
                     this%cm3(jk,:,ip) = fctr * rv(:)
                  end select
               end do
            end do
         end do
      end do
      deallocate(dek, qok, qak, rrk)

   END SUBROUTINE calc_grid_metrics

END MODULE mo_grid