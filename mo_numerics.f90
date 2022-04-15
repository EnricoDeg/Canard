!*****
!***** 3D SOLVER MODULE
!*****

module mo_numerics

   use mo_mpi, ONLY : p_null_req, p_isend, p_irecv, p_waitall, &
                      p_send, myid
   use mainvar3d
   use mo_utils
   implicit none

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

   SUBROUTINE init_penta

      do nn=1,3
         select case(nn)
         case(1)
            is=0
            ie=is+lxi
         case(2)
            is=lxi+1
            ie=is+let
         case(3)
            is=lxi+let+2
            ie=is+lze
         end select
         do ip=0,1
            np=nbc(nn,ip)
            select case(np)
            case(10,20,25,30)
               ndf(nn,ip,0)=0
               ndf(nn,ip,1)=0
            case(35,40,45)
               ndf(nn,ip,0)=1
               ndf(nn,ip,1)=1
            end select
         end do
         ns=ndf(nn,0,0)
         ne=ndf(nn,1,0)
         call penta(xu(:,:),xl(:,:),is,ie,ns,ne,0)
         ns=ndf(nn,0,1)
         ne=ndf(nn,1,1)
         call penta(yu(:,:),yl(:,:),is,ie,ns,ne,1)
      end do
  
   END SUBROUTINE init_penta

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



!===== SUBROUTINE FOR MPI IMPLEMENTATION

   subroutine mpigo(nt,nrt,n45,itag)

      integer(kind=ni),intent(in) :: nt,nrt,n45,itag

      select case(nt)
      case(0)
         mp=lmd
      case(1)
         mp=lmf
      end select

      call p_null_req
      do nn=1,3
         nz=(1-nrt)*(nn-1)+1
         if(nt==0) then
            select case(nn)
            case(1)
               send=>send01
               recv=>recv01
            case(2)
               send=>send02
               recv=>recv02
            case(3)
               send=>send03
               recv=>recv03
            end select
         else
            select case(nn)
            case(1)
               send=>send11
               recv=>recv11
            case(2)
               send=>send12
               recv=>recv12
            case(3)
               send=>send13
               recv=>recv13
            end select
         end if
         do ip=0,1
            iq=1-ip
            is=ip*ijk(1,nn)
            ie=1-2*ip
            select case(nbc(nn,ip))
            case(35)
               ra0=zero
               ii=1
            case(40)
               ra0=zero
               ii=0
            case(45)
               ra0=n45
               ii=1
            end select
            if(ndf(nn,ip,nt)==1) then
               do k=0,ijk(3,nn)
                  kp=k*(ijk(2,nn)+1)
                  do j=0,ijk(2,nn)
                     jk=kp+j
                     l=indx3(is,j,k,nn)
                     res=ra0*rr(l,nz)
                     do i=0,mp
                        l=indx3(is+ie*(i+ii),j,k,nn)
                        sap(i)=rr(l,nz)
                     end do
                     send(jk,0,ip)=sum(pbco(0:mp,0,nt)*sap(0:mp))-res*pbcot(0,nt)
                     send(jk,1,ip)=sum(pbco(0:mp,1,nt)*sap(0:mp))-res*pbcot(1,nt)
                     send(jk,nt+1,ip)=send(jk,nt+1,ip)+nt*(sap(0)-res-send(jk,nt+1,ip))
                  end do
               end do
               if(nt==0) then
                  call p_isend(send(:,:,ip), mcd(nn,ip), itag+iq, 2*nbsize(nn))
                  call p_irecv(recv(:,:,ip), mcd(nn,ip), itag+ip, 2*nbsize(nn))
               else
                  call p_isend(send(:,:,ip), mcd(nn,ip), itag+iq, 3*nbsize(nn))
                  call p_irecv(recv(:,:,ip), mcd(nn,ip), itag+ip, 3*nbsize(nn))
               end if
            end if
         end do
      end do
      call p_waitall

      if(n45==n45go) then
         do nn=1,3
            nz=(1-nrt)*(nn-1)+1
            if(nt==0) then
               select case(nn)
               case(1)
                  recv=>recv01
               case(2)
                  recv=>recv02
               case(3)
                  recv=>recv03
               end select
            else
               select case(nn)
               case(1)
                  recv=>recv11
               case(2)
                  recv=>recv12
               case(3)
                  recv=>recv13
               end select
            end if
            do ip=0,1
               is=ip*ijk(1,nn)
               if(nbc(nn,ip)==45) then
                  do k=0,ijk(3,nn)
                     kp=k*(ijk(2,nn)+1)
                     do j=0,ijk(2,nn)
                        jk=kp+j
                        l=indx3(is,j,k,nn)
                        recv(jk,0,ip)=recv(jk,0,ip)+rr(l,nz)*pbcot(0,nt)
                        recv(jk,1,ip)=recv(jk,1,ip)+rr(l,nz)*pbcot(1,nt)
                        recv(jk,nt+1,ip)=recv(jk,nt+1,ip)+nt*rr(l,nz)
                     end do
                  end do
               end if
            end do
         end do
      end if

   end subroutine mpigo

!===== SUBROUTINE FOR COMPACT FINITE DIFFERENTIATING

   subroutine deriv(nn,nz,m)

      integer(kind=ni),intent(in) :: nn,nz,m

      nt=0
      ns=ndf(nn,0,0)
      ne=ndf(nn,1,0)

      select case(nn)
      case(1)
         is=0
         ie=is+lxi
         recv=>recv01
         drva=>drva1
      case(2)
         is=lxi+1
         ie=is+let
         recv=>recv02
         drva=>drva2
      case(3)
         is=lxi+let+2
         ie=is+lze
         recv=>recv03
         drva=>drva3
      end select

      do k=0,ijk(3,nn)
         kp=k*(ijk(2,nn)+1)
         do j=0,ijk(2,nn)
            jk=kp+j
            do i=is,ie
               l=indx3(i-is,j,k,nn)
               li(i)=l
               sa(i)=rr(l,nz)
            end do
            select case(ns)
            case(0)
               sb(is)=sum((/a01,a02,a03,a04/)*(sa(is+(/1,2,3,4/))-sa(is)))
               sb(is+1)=sum((/a10,a12,a13,a14/)*(sa(is+(/0,2,3,4/))-sa(is+1)))
            case(1)
               sb(is)=sum(pbci(0:lmd,0,nt)*sa(is:is+lmd))+recv(jk,0,0)
               sb(is+1)=sum(pbci(0:lmd,1,nt)*sa(is:is+lmd))+recv(jk,1,0)
            end select
            do i=is+2,ie-2
               sb(i)=aa*(sa(i+1)-sa(i-1))+ab*(sa(i+2)-sa(i-2))
            end do
            select case(ne)
            case(0)
               sb(ie)=sum((/a01,a02,a03,a04/)*(sa(ie)-sa(ie-(/1,2,3,4/))))
               sb(ie-1)=sum((/a10,a12,a13,a14/)*(sa(ie-1)-sa(ie-(/0,2,3,4/))))
            case(1)
               sb(ie)=-sum(pbci(0:lmd,0,nt)*sa(ie:ie-lmd:-1))-recv(jk,0,1)
               sb(ie-1)=-sum(pbci(0:lmd,1,nt)*sa(ie:ie-lmd:-1))-recv(jk,1,1)
            end select
            sa(is)=sb(is)
            sa(is+1)=sb(is+1)-xl(is+1,2)*sa(is)
            do i=is+2,ie
               sa(i)=sb(i)-xl(i,1)*sa(i-2)-xl(i,2)*sa(i-1)
            end do
            sb(ie)=xu(ie,1)*sa(ie)
            sb(ie-1)=xu(ie-1,1)*sa(ie-1)-xu(ie-1,2)*sb(ie)
            do i=ie-2,is,-1
               sb(i)=xu(i,1)*sa(i)-xu(i,2)*sb(i+1)-xu(i,3)*sb(i+2)
            end do
            do i=is,ie
               l=li(i)
               rr(l,nn)=sb(i)
            end do
            drva(jk,m,0)=sb(is)
            drva(jk,m,1)=sb(ie)
         end do
      end do

   end subroutine deriv

!===== SUBROUTINE FOR COMPACT FILTERING

   subroutine filte(nn,nz)

      integer(kind=ni),intent(in) :: nn,nz

      nt=1
      ns=ndf(nn,0,1)
      ne=ndf(nn,1,1)

      select case(nn)
      case(1)
         is=0
         ie=is+lxi
         recv=>recv11
      case(2)
         is=lxi+1
         ie=is+let
         recv=>recv12
      case(3)
         is=lxi+let+2
         ie=is+lze
         recv=>recv13
      end select

      do k=0,ijk(3,nn)
         kp=k*(ijk(2,nn)+1)
         do j=0,ijk(2,nn)
            jk=kp+j
            do i=is,ie
               l=indx3(i-is,j,k,nn)
               li(i)=l
               sa(i)=rr(l,nz)
            end do
            select case(ns)
            case(0)
               sb(is)=sum(fbc(:,0)*(sa(is+(/1,2,3,4,5/))-sa(is)))
               sb(is+1)=sum(fbc(:,1)*(sa(is+(/0,2,3,4,5/))-sa(is+1)))
               sb(is+2)=sum(fbc(:,2)*(sa(is+(/0,1,3,4,5/))-sa(is+2)))
            case(1)
               ra2=sa(is+2)+sa(is+2)
               sb(is)=sum(pbci(0:lmf,0,nt)*sa(is:is+lmf))+recv(jk,0,0)
               sb(is+1)=sum(pbci(0:lmf,1,nt)*sa(is:is+lmf))+recv(jk,1,0)
               sb(is+2)=fa*(sa(is+1)+sa(is+3)-ra2)+fb*(sa(is)+sa(is+4)-ra2)+fc*(recv(jk,2,0)+sa(is+5)-ra2)
            end select
            do i=is+3,ie-3
               res=sa(i)+sa(i)
               sb(i)=fa*(sa(i-1)+sa(i+1)-res)+fb*(sa(i-2)+sa(i+2)-res)+fc*(sa(i-3)+sa(i+3)-res)
            end do
            select case(ne)
            case(0)
               sb(ie)=sum(fbc(:,0)*(sa(ie-(/1,2,3,4,5/))-sa(ie)))
               sb(ie-1)=sum(fbc(:,1)*(sa(ie-(/0,2,3,4,5/))-sa(ie-1)))
               sb(ie-2)=sum(fbc(:,2)*(sa(ie-(/0,1,3,4,5/))-sa(ie-2)))
            case(1)
               ra2=sa(ie-2)+sa(ie-2)
               sb(ie)=sum(pbci(0:lmf,0,nt)*sa(ie:ie-lmf:-1))+recv(jk,0,1)
               sb(ie-1)=sum(pbci(0:lmf,1,nt)*sa(ie:ie-lmf:-1))+recv(jk,1,1)
               sb(ie-2)=fa*(sa(ie-3)+sa(ie-1)-ra2)+fb*(sa(ie-4)+sa(ie)-ra2)+fc*(sa(ie-5)+recv(jk,2,1)-ra2)
            end select
            sa(is)=sb(is)
            sa(is+1)=sb(is+1)-yl(is+1,2)*sa(is)
            do i=is+2,ie
               sa(i)=sb(i)-yl(i,1)*sa(i-2)-yl(i,2)*sa(i-1)
            end do
            sb(ie)=yu(ie,1)*sa(ie)
            sb(ie-1)=yu(ie-1,1)*sa(ie-1)-yu(ie-1,2)*sb(ie)
            do i=ie-2,is,-1
               sb(i)=yu(i,1)*sa(i)-yu(i,2)*sb(i+1)-yu(i,3)*sb(i+2)
            end do
            do i=is,ie
               l=li(i)
               rr(l,nz)=rr(l,nz)+sb(i)
            end do
         end do
      end do

   end subroutine filte

!=====

end module mo_numerics

!*****