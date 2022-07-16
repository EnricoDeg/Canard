!*****
!***** MAIN VARIABLES & DATA FOR 3D NAVIER-STOKES/EULER SOLVER
!*****

module mo_vars
   use mo_kind, ONLY : ni, nr, ieee32, int64
   implicit none
   PUBLIC

!===== ALLOCATABLE MAIN ARRAYS

   real(kind=nr),dimension(:,:),allocatable :: qa, de

   real(kind=nr),dimension(:),allocatable :: p

   real(kind=nr),dimension(:,:),allocatable :: rr, ss

!===== CONSTANT-SIZED MAIN VARIABLES
   integer(kind=ni) :: nrecs,nrecd
   integer(kind=ni) :: mbk

   real(kind=nr) :: aoo,srefoo,srefp1dre
   real(kind=nr) :: hv2,ao,bo,aoi,sqrtrema,sqrtremai

   character(5) :: cnnode
   character(16) :: cgrid
   character(18) :: cdata

   CONTAINS

   SUBROUTINE allocate_memory(lmx, nbsize)
      integer(kind=ni), intent(IN)               :: lmx
      integer(kind=ni), dimension(3), intent(IN) :: nbsize

      allocate(qa(0:lmx,5),de(0:lmx,5))
      allocate(rr(0:lmx,3),ss(0:lmx,3))
      allocate(p(0:lmx))

   END SUBROUTINE allocate_memory

end module mo_vars
