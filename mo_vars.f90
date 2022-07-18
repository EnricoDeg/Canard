!*****
!***** MAIN VARIABLES & DATA FOR 3D NAVIER-STOKES/EULER SOLVER
!*****

module mo_vars
   use mo_kind, ONLY : ni, nr, ieee32, int64
   implicit none
   PUBLIC

!===== ALLOCATABLE MAIN ARRAYS
   real(kind=nr),dimension(:,:),allocatable :: rr

!===== CONSTANT-SIZED MAIN VARIABLES
   integer(kind=ni) :: nrecs,nrecd
   integer(kind=ni) :: mbk

   CONTAINS

   SUBROUTINE allocate_memory(lmx, nbsize)
      integer(kind=ni), intent(IN)               :: lmx
      integer(kind=ni), dimension(3), intent(IN) :: nbsize

      allocate(rr(0:lmx,3))

   END SUBROUTINE allocate_memory

end module mo_vars
