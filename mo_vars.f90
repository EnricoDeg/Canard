!*****
!***** MAIN VARIABLES & DATA FOR 3D NAVIER-STOKES/EULER SOLVER
!*****

module mo_vars
   use mo_kind, ONLY : ni, nr, ieee32, int64
   implicit none
   PUBLIC

!===== CONSTANT-SIZED MAIN VARIABLES
   integer(kind=ni) :: nrecs,nrecd
   integer(kind=ni) :: mbk

end module mo_vars
