!*****
!***** BASIC SUBROUTINES
!*****

module mo_utils
   use mo_kind,       ONLY : ni, nr
   use mo_parameters, ONLY : one
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

      integer(kind=ni) :: iik, jjk, mk
      real(kind=nr)    :: ra0k, ra1k

      rx(:,:) = ax(:,:)
      ipvt(:) = (/(iik, iik = is,ie)/)
      do iik = is,ie
         imax(:) = maxloc( abs( rx(iik:ie,iik) ) )
         mk = iik - 1 + imax(1)
         if(mk /= iik) then
            ipvt((/mk,iik/)) = ipvt((/iik,mk/))
            rx((/mk,iik/),:) = rx((/iik,mk/),:)
         end if
         ra0k = one / rx(iik,iik)
         temp(:) = rx(:,iik)
         do jjk = is,ie
            ra1k = ra0k * rx(iik,jjk)
            rx(:,jjk) = rx(:,jjk) - ra1k * temp(:)
            rx(iik,jjk) = ra1k
         end do
         rx(:,iik) = -ra0k * temp(:)
         rx(iik,iik) = ra0k
      end do
      sx(:,ipvt(:)) = rx(:,:)

   end subroutine mtrxi

!===== FUNCTION FOR MAIN INDEX TRANSFORMATION IN 3D

   function indx3(i,j,k,nn, lxi, let) result(lm)
      !$ACC ROUTINE SEQ

      integer(kind=ni),intent(in) :: i,j,k,nn
      integer(kind=ni),intent(in) :: lxi, let
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
