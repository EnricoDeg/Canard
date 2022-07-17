!*****
!***** MAIN VARIABLES & DATA FOR 3D NAVIER-STOKES/EULER SOLVER
!*****

module mo_vars
   use mo_kind, ONLY : ni, nr, ieee32, int64
   implicit none
   PUBLIC

!===== ALLOCATABLE MAIN ARRAYS

   real(kind=nr),dimension(:,:),allocatable :: de

   real(kind=nr),dimension(:,:),allocatable :: rr, ss

!===== CONSTANT-SIZED MAIN VARIABLES
   integer(kind=ni) :: nrecs,nrecd
   integer(kind=ni) :: mbk

   real(kind=nr) :: srefoo,srefp1dre

   CONTAINS

   SUBROUTINE allocate_memory(lmx, nbsize)
      integer(kind=ni), intent(IN)               :: lmx
      integer(kind=ni), dimension(3), intent(IN) :: nbsize

      allocate(de(0:lmx,5))
      allocate(rr(0:lmx,3),ss(0:lmx,3))

   END SUBROUTINE allocate_memory

end module mo_vars
