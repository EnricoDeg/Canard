!*****
!***** 3D SOLVER MODULE
!*****

module subroutines3d

   use mo_mpi
   use mainvar3d
   use subroutineso
   implicit none

   contains

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

!===== SUBROUTINE FOR CALCULATING VORTICITY

   subroutine vorti

      ss(:,1)=one/qa(:,1); de(:,1:3)=zero

      rr(:,1)=ss(:,1)*qa(:,2)
      m=1
      call mpigo(0,nrone,n45no,m)
      call deriv(3,1,m)
      call deriv(2,1,m)
      call deriv(1,1,m)
      de(:,2)=de(:,2)+rr(:,1)*xim(:,3)+rr(:,2)*etm(:,3)+rr(:,3)*zem(:,3)
      de(:,3)=de(:,3)-rr(:,1)*xim(:,2)-rr(:,2)*etm(:,2)-rr(:,3)*zem(:,2)

      rr(:,1)=ss(:,1)*qa(:,3)
      m=2
      call mpigo(0,nrone,n45no,m)
      call deriv(3,1,m)
      call deriv(2,1,m)
      call deriv(1,1,m)
      de(:,3)=de(:,3)+rr(:,1)*xim(:,1)+rr(:,2)*etm(:,1)+rr(:,3)*zem(:,1)
      de(:,1)=de(:,1)-rr(:,1)*xim(:,3)-rr(:,2)*etm(:,3)-rr(:,3)*zem(:,3)

      rr(:,1)=ss(:,1)*qa(:,4)
      m=3
      call mpigo(0,nrone,n45no,m)
      call deriv(3,1,m)
      call deriv(2,1,m)
      call deriv(1,1,m)
      de(:,1)=de(:,1)+rr(:,1)*xim(:,2)+rr(:,2)*etm(:,2)+rr(:,3)*zem(:,2)
      de(:,2)=de(:,2)-rr(:,1)*xim(:,1)-rr(:,2)*etm(:,1)-rr(:,3)*zem(:,1)

      de(:,1)=de(:,1)*yaco(:); de(:,2)=de(:,2)*yaco(:); de(:,3)=de(:,3)*yaco(:)

   end subroutine vorti

!===== SUBROUTINE FOR FINDING VARIABLE MIN/MAX VALUES FOR TECPLOT DATA FILE

   subroutine vminmax(nn)

      integer(kind=ni),intent(in) :: nn

      varmin(nn)=minval(varr)
      varmax(nn)=maxval(varr)
      varm(:,myid)=(/varmin(nn),varmax(nn)/)

      call p_null_req
      itag=nn
      if(myid==mo(mb)) then
         mps=mo(mb)
         mpe=mps+nbpc(mb,1)*nbpc(mb,2)*nbpc(mb,3)-1
         do mp=mps+1,mpe
            call p_irecv(varm(:,mp), mp, itag, 2)
         end do
         call p_waitall
         varmin(nn)=minval(varm(0,mps:mpe))
         varmax(nn)=maxval(varm(1,mps:mpe))
      else
         call p_send(varm(:,myid), mo(mb), itag, 2)
      end if

   end subroutine vminmax

!=====

end module subroutines3d

!*****