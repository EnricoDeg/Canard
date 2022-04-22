!*****
!***** PHYSICS MODULE
!*****

MODULE mo_physics
   use mo_vars
   implicit none

   real(kind=nr) :: reoo,tempoo,amach1,amach2,amach3
   real(kind=nr) :: amachoo
   real(kind=nr), dimension(3) :: uoo

   contains

   subroutine init_physics

      open(9,file='input.physics',status='old')
      read(9,*) cinput,reoo
      read(9,*) cinput,tempoo
      read(9,*) cinput,amach1
      read(9,*) cinput,amach2
      read(9,*) cinput,amach3
      close(9)

      amachoo=sqrt(amach1*amach1+amach2*amach2+amach3*amach3)
      if(amachoo>sml) then
         reoo=reoo/amachoo
      end if
      srefoo=111/tempoo
      srefp1dre=(srefoo+one)/reoo
      sqrtrema=sqrt(reoo*amachoo)
      sqrtremai=one/max(sqrtrema,sml)
      uoo(:)=(/amach1,amach2,amach3/)

   end subroutine init_physics

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

!===== SUBROUTINE FOR MOVING FRAME VELOCITIES

   subroutine movef(dtko,dtk)

      real(kind=nr),intent(in) :: dtko,dtk

      if(nsmf==0) then
         ra0=pi/timf
         ra1=ra0*min(timo,timf)
         ra2=ra0*min(timo+dtko,timf)

         fctr=one-cos(ra1)
         dfdt=ra0*sin(ra2)
         progmf=half*(fctr+dtk*dfdt)
         umf(:)=progmf*uoo(:)

         fctr=sin(ra1)
         dfdt=ra0*cos(ra2)
         progmf=half*ra0*(fctr+dtk*dfdt)
         dudtmf(:)=progmf*uoo(:)
      else
         umf(:)=uoo(:)
         dudtmf(:)=zero
      end if

   end subroutine movef

!=====


END MODULE mo_physics