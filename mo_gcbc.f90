!*****
!***** GENERAL CHARACTERISTIC BOUNDARY CONDITIONS MODULE
!*****

MODULE mo_gcbc
   use mainvar3d
   use mo_mpi
   use subroutineso
   use subroutines3d
   IMPLICIT NONE
   PUBLIC
   CONTAINS

   SUBROUTINE gcbc_init

      cbca(:,:)=zero
      cbca(1,1:2)=(/alpha01,beta02/)
      cbca(2,1:3)=(/one,alpha12,beta13/)
      cbca(3,1:4)=(/alpha,one,alpha,beta/)
      do i=4,mbci
         cbca(i,i-3:i)=(/beta,alpha,one,alpha/)
         if(i<mbci) then
            cbca(i,i+1)=beta
         end if
      end do
      rbci(:)=zero
      rbci(1:3)=(/one,alpha10,beta/)
      call mtrxi(cbca,cbcs,1,mbci)
      sbci(:)=-matmul(cbcs(:,:),rbci(:))
      fctr=pi/(mbci+1); res=zero
      do i=1,mbci
         res=res+one
         sbci(i)=half*sbci(i)*(one+cos(res*fctr))
      end do
      ll=-1
      npex(:)=0
      nrr(:)=0
      rr(:,1)=zero
      do nn=1,3
         do ip=0,1
            np=nbc(nn,ip)
            i=ip*ijk(1,nn)
            iq=1-2*ip
            if((np-10)*(np-20)*(np-25)*(np-30)==0) then
               do k=0,ijk(3,nn)
                  do j=0,ijk(2,nn)
                     l=indx3(i,j,k,nn)
                     ll=ll+1
                     res=one/yaco(l)
                     rr(l,1)=rr(l,1)+one
                     rr(ll,2)=res
                     rr(ll,3)=l+sml
                     do ii=1,mbci
                        l=indx3(i+iq*ii,j,k,nn)
                        ll=ll+1
                        rr(l,1)=rr(l,1)+one
                        rr(ll,2)=res*sbci(ii)
                        rr(ll,3)=l+sml
                     end do
                  end do
               end do
            end if
            if((np-30)*(np-35)*(np-45)==0) then
               do k=0,ijk(3,nn)
                  do j=0,ijk(2,nn)
                     l=indx3(i,j,k,nn)
                     nrr(l)=1
                  end do
               end do
            end if
         end do
      end do
      do l=0,lmx
         nrr(l)=min(nrr(l)+npex(l),1)
      end do
      lq=ll
      allocate(sbcc(0:lq))
      do ll=0,lq
         l=rr(ll,3)
         sbcc(ll)=rr(ll,2)/rr(l,1)
      end do

   END SUBROUTINE gcbc_init

   SUBROUTINE gcbc_setup

      do nn=1,3
         select case(nn)
         case(1)
            drva=>drva1
            cm=>cm1
         case(2)
            drva=>drva2
            cm=>cm2
         case(3)
            drva=>drva3
            cm=>cm3
         end select
         do ip=0,1
            np=nbc(nn,ip)
            i=ip*ijk(1,nn)
            if((np-10)*(np-20)*(np-25)*(np-30)==0) then
               ra0=(20-np)*(25-np)*(30-np)/3000
               ra1=one-ra0
               do k=0,ijk(3,nn)
                  kp=k*(ijk(2,nn)+1)
                  do j=0,ijk(2,nn)
                     jk=kp+j
                     l=indx3(i,j,k,nn)
                     call eleme(l,cm(jk,:,ip))
                     call xtq2r(cm(jk,:,ip))
                     cha(:)=ra0*drva(jk,:,ip)+ra1*de(l,:)
                     drva(jk,:,ip)=matmul(xt(:,:),yaco(l)*cha(:))
                  end do
               end do
            end if
         end do
      end do

   END SUBROUTINE gcbc_setup

   SUBROUTINE gcbc_comm

      call p_null_req
      itag=30
      do nn=1,3
         select case(nn)
         case(1)
            drva=>drva1
            drvb=>drvb1
         case(2)
            drva=>drva2
            drvb=>drvb2
         case(3)
            drva=>drva3
            drvb=>drvb3
         end select
         do ip=0,1
            iq=1-ip
            np=nbc(nn,ip)
            if((np-30)*(abs(nextgcic-1)+abs((np-20)*(np-25)))==0) then
               call p_isend(drva(:,:,ip), mcd(nn,ip), itag+iq, 5*nbsize(nn))
               call p_irecv(drvb(:,:,ip), mcd(nn,ip), itag+ip, 5*nbsize(nn))
            end if
         end do
      end do
      call p_waitall

   END SUBROUTINE gcbc_comm

   SUBROUTINE gcbc_update

      ll=-1
      do nn=1,3
         select case(nn)
         case(1)
            drva=>drva1
            drvb=>drvb1
            cm=>cm1
         case(2)
            drva=>drva2
            drvb=>drvb2
            cm=>cm2
         case(3)
            drva=>drva3
            drvb=>drvb3
            cm=>cm3
         end select
         do ip=0,1
            np=nbc(nn,ip)
            i=ip*ijk(1,nn)
            iq=1-2*ip
            ra0=iq
            select case(np)
            case(10)
               do k=0,ijk(3,nn)
                  kp=k*(ijk(2,nn)+1)
                  do j=0,ijk(2,nn)
                     jk=kp+j
                     l=indx3(i,j,k,nn)
                     call eleme(l,cm(jk,:,ip))
                     cha(:)=drva(jk,:,ip)
                     dha(:)=drvb(jk,:,ip)
                     if(ra0*(vn+vs+ao)>zero) then
                        cha(4)=zero
                     end if
                     if(ra0*(vn+vs-ao)>zero) then
                        cha(5)=zero
                     end if
                     call xtr2q(cm(jk,:,ip))
                     dha(:)=matmul(xt(:,:),(cha(:)-drva(jk,:,ip)))
                     do ii=0,mbci
                        l=indx3(i+iq*ii,j,k,nn)
                        ll=ll+1
                        de(l,:)=de(l,:)+sbcc(ll)*dha(:)
                     end do
                  end do
               end do
            case(20,25)
               dtwi=one/(nkrk*dt+sml)
               do k=0,ijk(3,nn)
                  kp=k*(ijk(2,nn)+1)
                  do j=0,ijk(2,nn)
                     jk=kp+j
                     l=indx3(i,j,k,nn)
                     call eleme(l,cm(jk,:,ip))
                     cha(:)=drva(jk,:,ip)
                     dha(:)=drvb(jk,:,ip)
                     select case(npex(l))
                     case(0)
                        cha(4+ip)=cha(5-ip)+two*ra0*aoi*qa(l,1)*(sum(cm(jk,:,ip)*dudtmf(:))+dtwi*(vn+vs))
                     end select
                     call xtr2q(cm(jk,:,ip))
                     dha(:)=matmul(xt(:,:),(cha(:)-drva(jk,:,ip)))
                     do ii=0,mbci
                        l=indx3(i+iq*ii,j,k,nn)
                        ll=ll+1
                        de(l,:)=de(l,:)+sbcc(ll)*dha(:)
                     end do
                  end do
               end do
            case(30)
               do k=0,ijk(3,nn)
                  kp=k*(ijk(2,nn)+1)
                  do j=0,ijk(2,nn)
                     jk=kp+j
                     l=indx3(i,j,k,nn)
                     call eleme(l,cm(jk,:,ip))
                     cha(:)=drva(jk,:,ip)
                     dha(:)=drvb(jk,:,ip)
                     if(ra0*(vn+vs)>zero) then
                        cha(1:3)=dha(1:3)
                     end if
                     if(ra0*(vn+vs+ao)>zero) then
                        cha(4)=dha(4)
                     end if
                     if(ra0*(vn+vs-ao)>zero) then
                        cha(5)=dha(5)
                     end if
                     call xtr2q(cm(jk,:,ip))
                     dha(:)=matmul(xt(:,:),(cha(:)-drva(jk,:,ip)))
                     do ii=0,mbci
                        l=indx3(i+iq*ii,j,k,nn)
                        ll=ll+1
                        de(l,:)=de(l,:)+sbcc(ll)*dha(:)
                     end do
                  end do
               end do
            end select
         end do
      end do

   END SUBROUTINE gcbc_update

   SUBROUTINE wall_condition_update

      do nn=1,3
         do ip=0,1
            np=nbc(nn,ip)
            i=ip*ijk(1,nn)
            select case(np)
            case(20)
               do k=0,ijk(3,nn)
                  do j=0,ijk(2,nn)
                     l=indx3(i,j,k,nn)
                     qa(l,5)=npex(l)*qa(l,5)+(one-npex(l))*(hamhamm1*qa(l,1)**gam+half*sum(qa(l,2:4)*qa(l,2:4))/qa(l,1))
                  end do
               end do
            case(25)
               ra0=hamhamm1*wtemp
               do k=0,ijk(3,nn)
                  do j=0,ijk(2,nn)
                     l=indx3(i,j,k,nn)
                     fctr=(one-npex(l))*qa(l,1)
                     qa(l,2:4)=npex(l)*qa(l,2:4)-fctr*umf(:)
                     qa(l,5)=npex(l)*qa(l,5)+(one-npex(l))*(ra0*qa(l,1)+half*sum(qa(l,2:4)*qa(l,2:4))/qa(l,1))
                  end do
               end do
            end select
         end do
      end do

   END SUBROUTINE wall_condition_update

!===== SUBROUTINE FOR TRANSFORMATION FROM Q TO R IN GCBC/GCIC

   subroutine xtq2r(cm)

      real(kind=nr),dimension(3),intent(in) :: cm

      ho=gamm1*aoi*aoi; bo=1-ho*hv2; co=aoi*vn; dm(:)=aoi*cm(:); rv(:)=ho*ve(:)

      xt(1,1)=bo*cm(1)+dm(2)*ve(3)-dm(3)*ve(2)
      xt(1,2)=cm(1)*rv(1)
      xt(1,3)=cm(1)*rv(2)+dm(3)
      xt(1,4)=cm(1)*rv(3)-dm(2)
      xt(1,5)=-ho*cm(1)

      xt(2,1)=bo*cm(2)+dm(3)*ve(1)-dm(1)*ve(3)
      xt(2,2)=cm(2)*rv(1)-dm(3)
      xt(2,3)=cm(2)*rv(2)
      xt(2,4)=cm(2)*rv(3)+dm(1)
      xt(2,5)=-ho*cm(2)

      xt(3,1)=bo*cm(3)+dm(1)*ve(2)-dm(2)*ve(1)
      xt(3,2)=cm(3)*rv(1)+dm(2)
      xt(3,3)=cm(3)*rv(2)-dm(1)
      xt(3,4)=cm(3)*rv(3)
      xt(3,5)=-ho*cm(3)

      xt(4,1)=one-bo-co
      xt(4,2)=dm(1)-rv(1)
      xt(4,3)=dm(2)-rv(2)
      xt(4,4)=dm(3)-rv(3)
      xt(4,5)=ho

      xt(5,1)=one-bo+co
      xt(5,2)=-dm(1)-rv(1)
      xt(5,3)=-dm(2)-rv(2)
      xt(5,4)=-dm(3)-rv(3)
      xt(5,5)=ho

   end subroutine xtq2r

!===== SUBROUTINE FOR INVERSE TRANSFORMATION FROM R TO Q IN GCBC/GCIC

   subroutine xtr2q(cm)

      real(kind=nr),dimension(3),intent(in) :: cm

      bo=hv2+hamm1*ao*ao; co=ao*vn; dm(:)=ao*cm(:)

      xt(1,1)=cm(1)
      xt(1,2)=cm(2)
      xt(1,3)=cm(3)
      xt(1,4)=half
      xt(1,5)=half

      xt(2,1)=cm(1)*ve(1)
      xt(2,2)=cm(2)*ve(1)-dm(3)
      xt(2,3)=cm(3)*ve(1)+dm(2)
      xt(2,4)=half*(ve(1)+dm(1))
      xt(2,5)=xt(2,4)-dm(1)

      xt(3,1)=cm(1)*ve(2)+dm(3)
      xt(3,2)=cm(2)*ve(2)
      xt(3,3)=cm(3)*ve(2)-dm(1)
      xt(3,4)=half*(ve(2)+dm(2))
      xt(3,5)=xt(3,4)-dm(2)

      xt(4,1)=cm(1)*ve(3)-dm(2)
      xt(4,2)=cm(2)*ve(3)+dm(1)
      xt(4,3)=cm(3)*ve(3)
      xt(4,4)=half*(ve(3)+dm(3))
      xt(4,5)=xt(4,4)-dm(3)

      xt(5,1)=hv2*cm(1)+dm(3)*ve(2)-dm(2)*ve(3)
      xt(5,2)=hv2*cm(2)+dm(1)*ve(3)-dm(3)*ve(1)
      xt(5,3)=hv2*cm(3)+dm(2)*ve(1)-dm(1)*ve(2)
      xt(5,4)=half*(bo+co)
      xt(5,5)=xt(5,4)-co

   end subroutine xtr2q

!===== SUBROUTINE FOR ELEMENTARY VARIABLES IN GCBC/GCIC

   subroutine eleme(l,cm)

      integer(kind=ni),intent(in) :: l
      real(kind=nr),dimension(3),intent(in) :: cm

      rhoi=one/qa(l,1)
      ao=sqrt(gam*rhoi*p(l))
      aoi=one/ao
      ve(:)=rhoi*qa(l,2:4)
      hv2=half*(ve(1)*ve(1)+ve(2)*ve(2)+ve(3)*ve(3))
      vn=cm(1)*ve(1)+cm(2)*ve(2)+cm(3)*ve(3)
      vs=cm(1)*umf(1)+cm(2)*umf(2)+cm(3)*umf(3)

   end subroutine eleme

END MODULE mo_gcbc