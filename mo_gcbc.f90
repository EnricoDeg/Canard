MODULE mo_gcbc
  use mainvar3d
  use mo_mpi
  use subroutineso
  use subroutines3d
  IMPLICIT NONE
  PUBLIC
  CONTAINS

  SUBROUTINE gcbc_init

    cbca(:,:)=zero; cbca(1,1:2)=(/alpha01,beta02/); cbca(2,1:3)=(/one,alpha12,beta13/); cbca(3,1:4)=(/alpha,one,alpha,beta/)
 do i=4,mbci
    cbca(i,i-3:i)=(/beta,alpha,one,alpha/); if(i<mbci) then; cbca(i,i+1)=beta; end if
 end do
    rbci(:)=zero; rbci(1:3)=(/one,alpha10,beta/)
    call mtrxi(cbca,cbcs,1,mbci); sbci(:)=-matmul(cbcs(:,:),rbci(:))
    fctr=pi/(mbci+1); res=zero
 do i=1,mbci; res=res+one
    sbci(i)=half*sbci(i)*(one+cos(res*fctr))
 end do
    ll=-1; npex(:)=0; nrr(:)=0; rr(:,1)=zero
 do nn=1,3; do ip=0,1; np=nbc(nn,ip); i=ip*ijk(1,nn); iq=1-2*ip
 if((np-10)*(np-20)*(np-25)*(np-30)==0) then
 do k=0,ijk(3,nn); do j=0,ijk(2,nn); l=indx3(i,j,k,nn)
    ll=ll+1; res=one/yaco(l); rr(l,1)=rr(l,1)+one; rr(ll,2)=res; rr(ll,3)=l+sml
 do ii=1,mbci; l=indx3(i+iq*ii,j,k,nn)
    ll=ll+1; rr(l,1)=rr(l,1)+one; rr(ll,2)=res*sbci(ii); rr(ll,3)=l+sml
 end do
 end do; end do
 end if
 if(abs(nextgcic-1)+abs((np-20)*(np-25))==0) then
 do k=0,ijk(3,nn); do j=0,ijk(2,nn); l=indx3(i,j,k,nn)
!    call extrabcc(npex(l))
 end do; end do
 end if
 if((np-30)*(np-35)*(np-45)==0) then
 do k=0,ijk(3,nn); do j=0,ijk(2,nn); l=indx3(i,j,k,nn)
    nrr(l)=1
 end do; end do
 end if
 end do; end do
 do l=0,lmx
    nrr(l)=min(nrr(l)+npex(l),1)
 end do
    lq=ll; allocate(sbcc(0:lq))
 do ll=0,lq; l=rr(ll,3)
    sbcc(ll)=rr(ll,2)/rr(l,1)
 end do

  END SUBROUTINE gcbc_init

  SUBROUTINE gcbc_setup

 do nn=1,3
 select case(nn)
 case(1); drva=>drva1; cm=>cm1; case(2); drva=>drva2; cm=>cm2; case(3); drva=>drva3; cm=>cm3
 end select
 do ip=0,1; np=nbc(nn,ip); i=ip*ijk(1,nn)
 if((np-10)*(np-20)*(np-25)*(np-30)==0) then
    ra0=(20-np)*(25-np)*(30-np)/3000; ra1=one-ra0
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    call eleme(l,cm(jk,:,ip)); call xtq2r(cm(jk,:,ip))
    cha(:)=ra0*drva(jk,:,ip)+ra1*de(l,:); drva(jk,:,ip)=matmul(xt(:,:),yaco(l)*cha(:))
 end do
 end do
 end if
 end do
 end do

  END SUBROUTINE gcbc_setup

  SUBROUTINE gcbc_comm

    ir=0; itag=30
 do nn=1,3
 select case(nn)
 case(1); drva=>drva1; drvb=>drvb1; case(2); drva=>drva2; drvb=>drvb2; case(3); drva=>drva3; drvb=>drvb3
 end select
 do ip=0,1; iq=1-ip; np=nbc(nn,ip)
 if((np-30)*(abs(nextgcic-1)+abs((np-20)*(np-25)))==0) then
    ir=ir+1; call MPI_ISEND(drva(:,:,ip),5*nbsize(nn),MPI_REAL8,mcd(nn,ip),itag+iq,icom,ireq(ir),ierr)
    ir=ir+1; call MPI_IRECV(drvb(:,:,ip),5*nbsize(nn),MPI_REAL8,mcd(nn,ip),itag+ip,icom,ireq(ir),ierr)
 end if
 end do
 end do
 if(ir/=0) then
    call MPI_WAITALL(ir,ireq,ista,ierr)
 end if

  END SUBROUTINE gcbc_comm

  SUBROUTINE gcbc_update

    ll=-1
 do nn=1,3
 select case(nn)
 case(1); drva=>drva1; drvb=>drvb1; cm=>cm1
 case(2); drva=>drva2; drvb=>drvb2; cm=>cm2
 case(3); drva=>drva3; drvb=>drvb3; cm=>cm3
 end select
 do ip=0,1; np=nbc(nn,ip); i=ip*ijk(1,nn); iq=1-2*ip; ra0=iq
 select case(np)
 case(10)
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    call eleme(l,cm(jk,:,ip)); cha(:)=drva(jk,:,ip); dha(:)=drvb(jk,:,ip)
    if(ra0*(vn+vs+ao)>zero) then; cha(4)=zero; end if
    if(ra0*(vn+vs-ao)>zero) then; cha(5)=zero; end if
    call xtr2q(cm(jk,:,ip)); dha(:)=matmul(xt(:,:),(cha(:)-drva(jk,:,ip)))
 do ii=0,mbci; l=indx3(i+iq*ii,j,k,nn)
    ll=ll+1; de(l,:)=de(l,:)+sbcc(ll)*dha(:)
 end do
 end do
 end do
 case(20,25); dtwi=one/(nkrk*dt+sml)
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    call eleme(l,cm(jk,:,ip)); cha(:)=drva(jk,:,ip); dha(:)=drvb(jk,:,ip)
 select case(npex(l))
 case(0); cha(4+ip)=cha(5-ip)+two*ra0*aoi*qa(l,1)*(sum(cm(jk,:,ip)*dudtmf(:))+dtwi*(vn+vs))
! case(1); call extrabcs
 end select
    call xtr2q(cm(jk,:,ip)); dha(:)=matmul(xt(:,:),(cha(:)-drva(jk,:,ip)))
 do ii=0,mbci; l=indx3(i+iq*ii,j,k,nn)
    ll=ll+1; de(l,:)=de(l,:)+sbcc(ll)*dha(:)
 end do
 end do
 end do
 case(30)
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    call eleme(l,cm(jk,:,ip)); cha(:)=drva(jk,:,ip); dha(:)=drvb(jk,:,ip)
    if(ra0*(vn+vs)>zero) then; cha(1:3)=dha(1:3); end if
    if(ra0*(vn+vs+ao)>zero) then; cha(4)=dha(4); end if
    if(ra0*(vn+vs-ao)>zero) then; cha(5)=dha(5); end if
    call xtr2q(cm(jk,:,ip)); dha(:)=matmul(xt(:,:),(cha(:)-drva(jk,:,ip)))
 do ii=0,mbci; l=indx3(i+iq*ii,j,k,nn)
    ll=ll+1; de(l,:)=de(l,:)+sbcc(ll)*dha(:)
 end do
 end do
 end do
 end select
 end do
 end do

  END SUBROUTINE gcbc_update

  SUBROUTINE wall_condition_update

 do nn=1,3; do ip=0,1; np=nbc(nn,ip); i=ip*ijk(1,nn)
 select case(np)
 case(20)
 do k=0,ijk(3,nn); do j=0,ijk(2,nn); l=indx3(i,j,k,nn)
    qa(l,5)=npex(l)*qa(l,5)+(one-npex(l))*(hamhamm1*qa(l,1)**gam+half*sum(qa(l,2:4)*qa(l,2:4))/qa(l,1))
 end do; end do
 case(25); ra0=hamhamm1*wtemp
 do k=0,ijk(3,nn); do j=0,ijk(2,nn); l=indx3(i,j,k,nn)
    fctr=(one-npex(l))*qa(l,1); qa(l,2:4)=npex(l)*qa(l,2:4)-fctr*umf(:)
    qa(l,5)=npex(l)*qa(l,5)+(one-npex(l))*(ra0*qa(l,1)+half*sum(qa(l,2:4)*qa(l,2:4))/qa(l,1))
 end do; end do
 end select
 end do; end do

  END SUBROUTINE wall_condition_update

END MODULE mo_gcbc