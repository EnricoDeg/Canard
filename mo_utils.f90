!*****
!***** BASIC SUBROUTINES
!*****

module mo_utils

   use mainvar3d
   implicit none
   public

   contains

!===== SUBROUTINE FOR MATRIX INVERSION

   subroutine mtrxi(ax,sx,is,ie)

      integer(kind=ni),intent(in) :: is,ie
      real(kind=nr),dimension(is:ie,is:ie),intent(in) :: ax
      real(kind=nr),dimension(is:ie,is:ie),intent(inout) :: sx

      integer(kind=ni),dimension(1) :: imax
      integer(kind=ni),dimension(is:ie) :: ipvt
      real(kind=nr),dimension(is:ie,is:ie) :: rx
      real(kind=nr),dimension(is:ie) :: temp

      rx(:,:)=ax(:,:)
      ipvt(:)=(/(i,i=is,ie)/)
      do i=is,ie
         imax(:)=maxloc(abs(rx(i:ie,i)))
         m=i-1+imax(1)
         if(m/=i) then
            ipvt((/m,i/))=ipvt((/i,m/))
            rx((/m,i/),:)=rx((/i,m/),:)
         end if
         ra0=one/rx(i,i)
         temp(:)=rx(:,i)
         do j=is,ie
            ra1=ra0*rx(i,j)
            rx(:,j)=rx(:,j)-ra1*temp(:)
            rx(i,j)=ra1
         end do
         rx(:,i)=-ra0*temp(:)
         rx(i,i)=ra0
      end do
      sx(:,ipvt(:))=rx(:,:)

   end subroutine mtrxi

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

!===== FUNCTION FOR MAIN INDEX TRANSFORMATION IN 3D

   function indx3(i,j,k,nn) result(lm)

      integer(kind=ni),intent(in) :: i,j,k,nn
      integer(kind=ni) :: lm

      select case(nn)
      case(1)
         lm=(k*(let+1)+j)*(lxi+1)+i
      case(2)
         lm=(j*(let+1)+i)*(lxi+1)+k
      case(3)
         lm=(i*(let+1)+k)*(lxi+1)+j
      end select

   end function indx3

!=====

end module mo_utils

!*****