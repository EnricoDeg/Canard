!*****
!***** BASIC SUBROUTINES
!*****

module subroutineso

   use mainvar3d
   implicit none
   public

   private :: fcbcm, fcint, sbcco
   real(kind=nr) :: alphf, betf

   contains

!===== EXTRA COEFFICIENTS FOR DOMAIN BOUNDARIES

   SUBROUTINE init_extracoeff_bounds

      call fcbcm(fltk,fltrbc)
      call fcint(fltk,half,alphf,betf,fa,fb,fc)
      albef(:,0,1)=(/zero,zero,one,alphf,betf/)
      albef(:,1,1)=(/zero,alphf,one,alphf,betf/)
      albef(:,2,1)=(/betf,alphf,one,alphf,betf/)

      pbco(:,:,:)=zero
      pbci(:,:,:)=zero
      call sbcco
      do nt=0,1
         do j=0,1
            ii=lmd+nt*(lmf-lmd)
            pbcot(j,nt)=sum(pbco(0:ii,j,nt))
         end do
      end do

   END SUBROUTINE init_extracoeff_bounds

!===== SUBROUTINE FOR CHOLESKY DECOMPOSITION OF PENTADIAGONAL MATRICES

   subroutine penta(xu,xl,is,ie,ns,ne,nt)

      integer(kind=ni),intent(in) :: is,ie,ns,ne,nt
      real(kind=nr),dimension(0:lim,3),intent(inout) :: xu
      real(kind=nr),dimension(0:lim,2),intent(inout) :: xl
      real(kind=nr),dimension(-2:2,0:2,0:1) :: albe
      real(kind=nr) :: alpho,beto

      if(nt==0) then
         albe(:,0,0)=(/zero,zero,one,alpha01,beta02/)
         albe(:,1,0)=(/zero,alpha10,one,alpha12,beta13/)
         albe(:,2,0)=(/beta,alpha,one,alpha,beta/)
         alpho=alpha; beto=beta
         albe(:,0,1)=(/zero,zero,one,alpha,beta/)
         albe(:,1,1)=(/zero,alpha,one,alpha,beta/)
         albe(:,2,1)=(/beta,alpha,one,alpha,beta/)
      else
         albe=albef; alpho=alphf; beto=betf
      end if

      do i=is,ie
         xl(i,:)=one; xu(i,:)=one
      end do
      i=is
      xu(i,1)=one
      xu(i,2)=albe(1,0,ns)
      xu(i,3)=albe(2,0,ns)
      i=is+1
      xl(i,2)=albe(-1,1,ns)*xu(i-1,1)
      xu(i,1)=one/(one-xu(i-1,2)*xl(i,2))
      xu(i,2)=albe(1,1,ns)-xu(i-1,3)*xl(i,2)
      xu(i,3)=albe(2,1,ns)
      i=is+2
      xl(i,1)=albe(-2,2,ns)*xu(i-2,1)
      xl(i,2)=(albe(-1,2,ns)-xu(i-2,2)*xl(i,1))*xu(i-1,1)
      xu(i,1)=one/(one-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
      xu(i,2)=albe(1,2,ns)-xu(i-1,3)*xl(i,2)
      xu(i,3)=albe(2,2,ns)
      do i=is+3,ie-3
         xl(i,1)=beto*xu(i-2,1)
         xl(i,2)=(alpho-xu(i-2,2)*xl(i,1))*xu(i-1,1)
         xu(i,1)=one/(one-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
         xu(i,2)=alpho-xu(i-1,3)*xl(i,2)
         xu(i,3)=beto
      end do
      i=ie-2
      xl(i,1)=albe(2,2,ne)*xu(i-2,1)
      xl(i,2)=(albe(1,2,ne)-xu(i-2,2)*xl(i,1))*xu(i-1,1)
      xu(i,1)=one/(one-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
      xu(i,2)=albe(-1,2,ne)-xu(i-1,3)*xl(i,2)
      xu(i,3)=albe(-2,2,ne)
      i=ie-1
      xl(i,1)=albe(2,1,ne)*xu(i-2,1)
      xl(i,2)=(albe(1,1,ne)-xu(i-2,2)*xl(i,1))*xu(i-1,1)
      xu(i,1)=one/(one-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
      xu(i,2)=albe(-1,1,ne)-xu(i-1,3)*xl(i,2)
      i=ie
      xl(i,1)=albe(2,0,ne)*xu(i-2,1)
      xl(i,2)=(albe(1,0,ne)-xu(i-2,2)*xl(i,1))*xu(i-1,1)
      xu(i,1)=one/(one-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
      do i=is,ie
         xu(i,2:3)=xu(i,2:3)*xu(i,1)
      end do

   end subroutine penta

!===== SUBROUTINE FOR BOUNDARY FILTER COEFFICIENTS

   subroutine fcbcm(fltk,fltrbc)
 
      real(kind=nr),intent(in) :: fltk,fltrbc
      real(kind=nr) :: alphz,betz,za,zb,zc

      ao=log(fltrbc)
      call fcint(fltk,half,alphz,betz,za,zb,zc)
      fctr=one/(one+alphz*fltrbc+betz*fltrbc**two)

      albef(:,0,0)=(/zero,zero,one,alphz*fctr,betz*fctr/)
      res=(fltrbc-1)*(za+zc+(fltrbc+1)*(zb+fltrbc*zc))/ao
      fbc(:,0)=(/za-5*res/3,zb+10*res/21,zc-5*res/42,5*res/252,-res/630/)*fctr

      albef(:,1,0)=(/zero,alphz+betz*fltrbc,one,alphz,betz/)
      res=(fltrbc-1)*(zb+zc*(fltrbc+1))/ao
      fbc(:,1)=(/za+zb+zc+1627*res/1260,za+10*res/21,zb-5*res/42,zc+5*res/252,-res/630/)

      albef(:,2,0)=(/betz,alphz,one,alphz,betz/)
      res=zc*(fltrbc-1)/ao
      fbc(:,2)=(/zb+zc+1627*res/1260,za-5*res/3,za-5*res/42,zb+5*res/252,zc-res/630/)

   end subroutine fcbcm

!===== SUBROUTINE FOR INTERIOR FILTER COEFFICIENTS

   subroutine fcint(fltk,fltr,alphz,betz,za,zb,zc)
 
      real(kind=nr),intent(in) :: fltk,fltr
      real(kind=nr),intent(inout) :: alphz,betz,za,zb,zc
      real(kind=nr),dimension(3) :: cosf

      cosf(1)=cos(fltk); cosf(2)=cos(2*fltk); cosf(3)=cos(3*fltk)
      fctr=1/(30+5*(7-16*fltr)*cosf(1)+2*(1+8*fltr)*cosf(2)-3*cosf(3))
      alphz=fctr*(20*(2*fltr-1)-30*cosf(1)+12*(2*fltr-1)*cosf(2)-2*cosf(3))
      betz=half*fctr*(2*(13-8*fltr)+(33-48*fltr)*cosf(1)+6*cosf(2)-cosf(3))
      za=60*(1-fltr)*cos(half*fltk)**4*fctr; zb=-2*za/5; zc=za/15

   end subroutine fcint

!===== SUBROUTINE FOR SUBDOMAIN-BOUNDARY COEFFICIENTS

   subroutine sbcco

      real(kind=nr),dimension(:,:),allocatable :: ax,bx,rx,sx

      do nt=0,1
         lp=2*nt-1
         if(nt==0) then
            ll=lmd
            is=1
            ie=2*(ll+1)
            allocate(ax(ie,ie),bx(ie,ie),rx(ie,ie),sx(ie,ie))
            ax(:,:)=0
            bx(:,:)=0
            ax(is,is:is+2)=(/one,alpha01,beta02/)
            bx(is,is:is+4)=(/-(a01+a02+a03+a04),a01,a02,a03,a04/)
            ax(is+1,is:is+3)=(/alpha10,one,alpha12,beta13/)
            bx(is+1,is:is+4)=(/a10,-(a10+a12+a13+a14),a12,a13,a14/)
            do i=is+2,ie-2
               ax(i,i-2:i+2)=(/beta,alpha,one,alpha,beta/)
               bx(i,i-2:i+2)=(/-ab,-aa,zero,aa,ab/)
            end do
            ax(ie-1,ie:ie-3:-1)=ax(is+1,is:is+3)
            bx(ie-1,ie:ie-4:-1)=-bx(is+1,is:is+4)
            ax(ie,ie:ie-2:-1)=ax(is,is:is+2)
            bx(ie,ie:ie-4:-1)=-bx(is,is:is+4)
         end if
         if(nt==1) then
            ll=lmf
            is=1
            ie=2*(ll+1)
            allocate(ax(ie,ie),bx(ie,ie),rx(ie,ie),sx(ie,ie))
            ax(:,:)=0
            bx(:,:)=0
            ax(is,is:is+2)=albef(0:2,0,0)
            bx(is,is+(/1,2,3,4,5/))=fbc(:,0)
            bx(is,is)=-sum(fbc(:,0))
            ax(is+1,is:is+3)=albef(-1:2,1,0)
            bx(is+1,is+(/0,2,3,4,5/))=fbc(:,1)
            bx(is+1,is+1)=-sum(fbc(:,1))
            ax(is+2,is:is+4)=albef(-2:2,2,0)
            bx(is+2,is+(/0,1,3,4,5/))=fbc(:,2)
            bx(is+2,is+2)=-sum(fbc(:,2))
            do i=is+3,ie-3
               ax(i,i-2:i+2)=(/betf,alphf,one,alphf,betf/)
               bx(i,i-3:i+3)=(/fc,fb,fa,-2*(fa+fb+fc),fa,fb,fc/)
            end do
            ax(ie-2,ie:ie-4:-1)=ax(is+2,is:is+4)
            bx(ie-2,ie:ie-5:-1)=bx(is+2,is:is+5)
            ax(ie-1,ie:ie-3:-1)=ax(is+1,is:is+3)
            bx(ie-1,ie:ie-5:-1)=bx(is+1,is:is+5)
            ax(ie,ie:ie-2:-1)=ax(is,is:is+2)
            bx(ie,ie:ie-5:-1)=bx(is,is:is+5)
         end if

         call mtrxi(ax(:,:),sx(:,:),is,ie)

         rx(:,:)=ax(:,:)
         i=ie/2-1
         rx(i,i+2)=0
         rx(i+1,i+2)=0
         rx(i+1,i+3)=0
         i=ie/2+2
         rx(i,i-2)=0
         rx(i-1,i-2)=0
         rx(i-1,i-3)=0
         ax(:,:)=matmul(rx(:,:),matmul(sx(:,:),bx(:,:)))
         i=ie/2+1
         pbco(ll:0:-1,0,nt)=ax(i,is:is+ll)
         pbci(0:ll,0,nt)=ax(i,is+ll+1:ie)
         i=ie/2+2
         pbco(ll:0:-1,1,nt)=ax(i,is:is+ll)
         pbci(0:ll,1,nt)=ax(i,is+ll+1:ie)
         deallocate(ax,bx,rx,sx)
      end do
   end subroutine sbcco

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

end module subroutineso

!*****