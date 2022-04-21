!*****
!***** PHYSICS MODULE
!*****

MODULE mo_physics
   use mo_vars
   implicit none

   contains

!===== INITIAL CONDITIONS

   subroutine initialo
      real(kind=nr),dimension(3) :: vee
      real(kind=nr) :: radv, k1, k2

      radv = 1.0
      k1 = 12.5
      k2 = 1.0

      do l=0,lmx
         ao = k2/2.0/pi * sqrt(exp(1 - k1**2 * (ss(l,1)**2 + ss(l,2)**2)/radv**2))
         bo=(one-half*gamm1*ao*ao)**hamm1
         qa(l,1)=bo;
         vee(:) = (/k1*ss(l,2)*ao/radv, -k1*ss(l,1)*ao/radv, zero/)
         hv2=half*(vee(1)*vee(1)+vee(2)*vee(2)+vee(3)*vee(3))
         qa(l,2:4)=bo*vee(:);
         qa(l,5)=hamhamm1*bo**gam+hv2*bo
      end do

   end subroutine initialo

!=====


END MODULE mo_physics