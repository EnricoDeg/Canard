!*****
!***** SUBROUTINES FOR 3D SOLVER
!*****

 module subroutines3d

 use mpi
 use mainvar3d
 use subroutineso
 implicit none

 contains

!===== SUBROUTINE FOR MPI IMPLEMENTATION

 subroutine mpigo(nt,nrt,n45,itag)

 integer(kind=ni),intent(in) :: nt,nrt,n45,itag

 select case(nt); case(0); mp=lmd; case(1); mp=lmf; case(2); mp=0; end select

    ir=0
 do nn=1,3; nz=(1-nrt)*(nn-1)+1
 select case(3*nt+nn)
 case(1); send=>send01; recv=>recv01; case(2); send=>send02; recv=>recv02; case(3); send=>send03; recv=>recv03
 case(4); send=>send11; recv=>recv11; case(5); send=>send12; recv=>recv12; case(6); send=>send13; recv=>recv13
 case(7); snde=>send21; rcve=>recv21; case(8); snde=>send22; rcve=>recv22; case(9); snde=>send23; rcve=>recv23
 end select
 do ip=0,1; iq=1-ip; is=ip*ijk(1,nn); ie=1-2*ip
 select case(nbc(nn,ip)); case(35); ra0=zero; ii=1; case(40); ra0=zero; ii=0; case(45); ra0=n45; ii=1; end select
 if(ndf(nn,ip,nt)==1) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(is,j,k,nn)
    res=ra0*rr(l,nz)
 do i=0,mp; l=indx3(is+ie*(i+ii),j,k,nn)
    sap(i)=rr(l,nz)
 end do
 select case(nt)
 case(0,1)
    send(jk,0,ip)=sum(pbco(0:mp,0,nt)*sap(0:mp))-res*pbcot(0,nt)
    send(jk,1,ip)=sum(pbco(0:mp,1,nt)*sap(0:mp))-res*pbcot(1,nt)
    send(jk,2,ip)=sap(0)-res
 if(nt==1) then
    send(jk,3,ip)=sap(1)-res
 end if
 case(2)
    snde(jk,ip)=sap(0)-res
 end select
 end do
 end do
 select case(nt)
 case(0,1)
    lmpi=(3+nt)*nbsize(nn)
    ir=ir+1; call MPI_ISEND(send(:,:,ip),lmpi,MPI_REAL8,mcd(nn,ip),itag+iq,icom,ireq(ir),ierr)
    ir=ir+1; call MPI_IRECV(recv(:,:,ip),lmpi,MPI_REAL8,mcd(nn,ip),itag+ip,icom,ireq(ir),ierr)
 case(2)
    ir=ir+1; call MPI_ISEND(snde(:,ip),nbsize(nn),MPI_REAL8,mcd(nn,ip),itag+iq,icom,ireq(ir),ierr)
    ir=ir+1; call MPI_IRECV(rcve(:,ip),nbsize(nn),MPI_REAL8,mcd(nn,ip),itag+ip,icom,ireq(ir),ierr)
 end select
 end if
 end do
 end do
 if(ir/=0) then
    call MPI_WAITALL(ir,ireq,ista,ierr)
 end if

 if(n45==n45go) then
 do nn=1,3; nz=(1-nrt)*(nn-1)+1
 select case(3*nt+nn)
 case(1); recv=>recv01; case(2); recv=>recv02; case(3); recv=>recv03
 case(4); recv=>recv11; case(5); recv=>recv12; case(6); recv=>recv13
 case(7); rcve=>recv21; case(8); rcve=>recv22; case(9); rcve=>recv23
 end select
 do ip=0,1; is=ip*ijk(1,nn)
 if(nbc(nn,ip)==45) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(is,j,k,nn)
 select case(nt)
 case(0,1)
    recv(jk,0,ip)=recv(jk,0,ip)+rr(l,nz)*pbcot(0,nt)
    recv(jk,1,ip)=recv(jk,1,ip)+rr(l,nz)*pbcot(1,nt)
    recv(jk,2,ip)=recv(jk,2,ip)+rr(l,nz)
 if(nt==1) then
    recv(jk,3,ip)=recv(jk,3,ip)+rr(l,nz)
 end if
 case(2)
    rcve(jk,ip)=rcve(jk,ip)+rr(l,nz)
 end select
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

    nt=0; ns=ndf(nn,0,nt); ne=ndf(nn,1,nt)

 select case(nn)
 case(1); is=0; ie=is+lxi; recv=>recv01; drva=>drva1
 case(2); is=lxi+1; ie=is+let; recv=>recv02; drva=>drva2
 case(3); is=lxi+let+2; ie=is+lze; recv=>recv03; drva=>drva3
 end select

 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j
 do i=is,ie; l=indx3(i-is,j,k,nn)
    li(i)=l; sa(i)=rr(l,nz)
 end do
 select case(ns)
 case(0)
    sb(is)=sum(dbc(0:5,0)*(sa(is+(/1,2,3,4,5,6/))-sa(is)))
    sb(is+1)=sum(dbc(0:5,1)*(sa(is+(/0,2,3,4,5,6/))-sa(is+1)))
!    sb(is+2)=aa2*(sa(is+3)-sa(is+1))+ab2*(sa(is+4)-sa(is))
    sb(is+2)=sum(dbc(0:5,2)*(sa(is+(/0,1,3,4,5,6/))-sa(is+2)))
 case(1)
    sb(is)=sum(pbci(0:lmd,0,nt)*sa(is:is+lmd))+recv(jk,0,0)
    sb(is+1)=sum(pbci(0:lmd,1,nt)*sa(is:is+lmd))+recv(jk,1,0)
    sb(is+2)=aa*(sa(is+3)-sa(is+1))+ab*(sa(is+4)-sa(is))+ac*(sa(is+5)-recv(jk,2,0))
 end select
 do i=is+3,ie-3
    sb(i)=aa*(sa(i+1)-sa(i-1))+ab*(sa(i+2)-sa(i-2))+ac*(sa(i+3)-sa(i-3))
 end do
 select case(ne)
 case(0)
    sb(ie)=sum(dbc(0:5,0)*(sa(ie)-sa(ie-(/1,2,3,4,5,6/))))
    sb(ie-1)=sum(dbc(0:5,1)*(sa(ie-1)-sa(ie-(/0,2,3,4,5,6/))))
!    sb(ie-2)=aa2*(sa(ie-1)-sa(ie-3))+ab2*(sa(ie)-sa(ie-4))
    sb(ie-2)=sum(dbc(0:5,2)*(sa(ie-2)-sa(ie-(/0,1,3,4,5,6/))))
 case(1)
    sb(ie)=-sum(pbci(0:lmd,0,nt)*sa(ie:ie-lmd:-1))-recv(jk,0,1)
    sb(ie-1)=-sum(pbci(0:lmd,1,nt)*sa(ie:ie-lmd:-1))-recv(jk,1,1)
    sb(ie-2)=aa*(sa(ie-1)-sa(ie-3))+ab*(sa(ie)-sa(ie-4))+ac*(recv(jk,2,1)-sa(ie-5))
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
    l=li(i); rr(l,nn)=sb(i)
 end do
    drva(jk,m,0)=sb(is); drva(jk,m,1)=sb(ie)
 end do
 end do

 end subroutine deriv

!===== SUBROUTINE FOR COMPACT FILTERING

 subroutine filte(nn)

 integer(kind=ni),intent(in) :: nn

    nt=1; ns=ndf(nn,0,nt); ne=ndf(nn,1,nt)

 select case(nn)
 case(1); is=0; ie=is+lxi; recv=>recv11
 case(2); is=lxi+1; ie=is+let; recv=>recv12
 case(3); is=lxi+let+2; ie=is+lze; recv=>recv13
 end select

    ao=quarter; bo=one/16.0_nr

 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j
 do i=is,ie; l=indx3(i-is,j,k,nn)
    li(i)=l; sa(i)=rr(l,1); sb(i)=sa(i)+sa(i)
 end do
 do i=is+1,ie-1
    sc(i)=sa(i-1)+sa(i+1)-sb(i)
 end do
 do i=is+2,ie-2
    sd(i)=sa(i-2)+sa(i+2)-sb(i)
 end do
 if(nshock==1) then
 select case(ns)
 case(0)
    se(is)=zero; se(is+1)=ao*sc(is+1)
 case(1)
    se(is)=ao*(recv(jk,2,0)+sa(is+1)-sb(is))+bo*(recv(jk,3,0)+sa(is+2)-sb(is))
    se(is+1)=ao*sc(is+1)+bo*(recv(jk,2,0)+sa(is+3)-sb(is+1))
 end select
 do i=is+2,ie-2
    se(i)=ao*sc(i)+bo*sd(i)
 end do
 select case(ne)
 case(0)
    se(ie)=zero; se(ie-1)=ao*sc(ie-1)
 case(1)
    se(ie)=ao*(sa(ie-1)+recv(jk,2,1)-sb(ie))+bo*(sa(ie-2)+recv(jk,3,1)-sb(ie))
    se(ie-1)=ao*sc(ie-1)+bo*(sa(ie-3)+recv(jk,2,1)-sb(ie-1))
 end select
 end if
 select case(ns)
 case(0)
    sb(is)=sum(fbc(0:4,0)*(sa(is+(/1,2,3,4,5/))-sa(is)))
    sb(is+1)=sum(fbc(0:4,1)*(sa(is+(/0,2,3,4,5/))-sa(is+1)))
    sb(is+2)=sum(fbc(0:4,2)*(sa(is+(/0,1,3,4,5/))-sa(is+2)))
 case(1)
    sb(is)=sum(pbci(0:lmf,0,nt)*sa(is:is+lmf))+recv(jk,0,0)
    sb(is+1)=sum(pbci(0:lmf,1,nt)*sa(is:is+lmf))+recv(jk,1,0)
    sb(is+2)=fa*sc(is+2)+fb*sd(is+2)+fc*(recv(jk,2,0)+sa(is+5)-sb(is+2))
 end select
 do i=is+3,ie-3
    sb(i)=fa*sc(i)+fb*sd(i)+fc*(sa(i-3)+sa(i+3)-sb(i))
 end do
 select case(ne)
 case(0)
    sb(ie)=sum(fbc(0:4,0)*(sa(ie-(/1,2,3,4,5/))-sa(ie)))
    sb(ie-1)=sum(fbc(0:4,1)*(sa(ie-(/0,2,3,4,5/))-sa(ie-1)))
    sb(ie-2)=sum(fbc(0:4,2)*(sa(ie-(/0,1,3,4,5/))-sa(ie-2)))
 case(1)
    sb(ie)=sum(pbci(0:lmf,0,nt)*sa(ie:ie-lmf:-1))+recv(jk,0,1)
    sb(ie-1)=sum(pbci(0:lmf,1,nt)*sa(ie:ie-lmf:-1))+recv(jk,1,1)
    sb(ie-2)=fa*sc(ie-2)+fb*sd(ie-2)+fc*(sa(ie-5)+recv(jk,2,1)-sb(ie-2))
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
    l=li(i); rr(l,1)=rr(l,1)+sb(i)
 end do
 if(nshock==1) then
 do i=is,ie
    l=li(i); rr(l,1)=rr(l,1)+(se(i)-sb(i))*de(l,nn)
 end do
 end if
 end do
 end do

 end subroutine filte

!===== SUBROUTINE FOR SHOCK DETECTOR

 subroutine shockdet(nn)

 integer(kind=ni),intent(in) :: nn

    nt=2; ns=ndf(nn,0,nt); ne=ndf(nn,1,nt)

 select case(nn)
 case(1); is=0; ie=is+lxi; rcve=>recv21
 case(2); is=lxi+1; ie=is+let; rcve=>recv22
 case(3); is=lxi+let+2; ie=is+lze; rcve=>recv23
 end select

 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j
 do i=is,ie; l=indx3(i-is,j,k,nn)
    li(i)=l; sa(i)=rr(l,1)
 end do
 select case(ns)
 case(0); sb(is)=zero; case(1); sb(is)=rcve(jk,0)+sa(is+1)-sa(is)-sa(is)
 end select
 do i=is+1,ie-1
    sb(i)=sa(i-1)+sa(i+1)-sa(is)-sa(is)
 end do
 select case(ne)
 case(0); sb(ie)=zero; case(1); sb(ie)=sa(ie-1)+rcve(jk,1)-sa(ie)-sa(ie)
 end select
 do i=is,ie
    l=li(i); de(l,nn)=one-one/(one+shockc*abs(sb(i)))
 end do
 end do
 end do

 end subroutine shockdet

!===== SUBROUTINE FOR CALCULATING VORTICITY

 subroutine vorti

    ss(:,1)=one/qa(:,1); de(:,1:3)=zero

    rr(:,1)=ss(:,1)*qa(:,2)
    m=1; call mpigo(0,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    de(:,2)=de(:,2)+rr(:,1)*xim(:,3)+rr(:,2)*etm(:,3)+rr(:,3)*zem(:,3)
    de(:,3)=de(:,3)-rr(:,1)*xim(:,2)-rr(:,2)*etm(:,2)-rr(:,3)*zem(:,2)

    rr(:,1)=ss(:,1)*qa(:,3)
    m=2; call mpigo(0,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    de(:,3)=de(:,3)+rr(:,1)*xim(:,1)+rr(:,2)*etm(:,1)+rr(:,3)*zem(:,1)
    de(:,1)=de(:,1)-rr(:,1)*xim(:,3)-rr(:,2)*etm(:,3)-rr(:,3)*zem(:,3)

    rr(:,1)=ss(:,1)*qa(:,4)
    m=3; call mpigo(0,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    de(:,1)=de(:,1)+rr(:,1)*xim(:,2)+rr(:,2)*etm(:,2)+rr(:,3)*zem(:,2)
    de(:,2)=de(:,2)-rr(:,1)*xim(:,1)-rr(:,2)*etm(:,1)-rr(:,3)*zem(:,1)

    de(:,1)=de(:,1)*yaco(:); de(:,2)=de(:,2)*yaco(:); de(:,3)=de(:,3)*yaco(:)

 end subroutine vorti

!===== SUBROUTINE FOR ELEMENTARY VARIABLES IN GCBC/GCIC

 subroutine eleme(l,cm)

 integer(kind=ni),intent(in) :: l
 real(kind=nr),dimension(3),intent(in) :: cm

    rhoi=one/qa(l,1); ao=sqrt(gam*rhoi*p(l)); aoi=one/ao
    ve(:)=rhoi*qa(l,2:4); hv2=half*(ve(1)*ve(1)+ve(2)*ve(2)+ve(3)*ve(3))
    vn=cm(1)*ve(1)+cm(2)*ve(2)+cm(3)*ve(3); vs=cm(1)*umf(1)+cm(2)*umf(2)+cm(3)*umf(3)

 end subroutine eleme

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

!===== SUBROUTINE FOR GENERATING TECPLOT DATA FILE

 subroutine techead(nf,n,mb,lh)

 integer(kind=ni),intent(in) :: nf,n,mb
 integer(kind=ni),intent(inout) :: lh

    lh=0
 if(mb==0) then
    write(nf,pos=4*lh+1) '#!TDV112'; lh=lh+2
    write(nf,pos=4*lh+1) 1; lh=lh+1 ! Header Section
    write(nf,pos=4*lh+1) int(min(n+2,2),kind=int32); lh=lh+1 ! File Type (0 = Full / 1 = Grid / 2 = Solution)
    cinput=cfilet(n); call strio(nf,lh,cinput) ! File Title
    write(nf,pos=4*lh+1) int(mq,kind=int32); lh=lh+1 ! Number of Variables
 if(n==-1) then
    cinput='x'; call strio(nf,lh,cinput)
    cinput='y'; call strio(nf,lh,cinput)
    cinput='z'; call strio(nf,lh,cinput)
 end if
 if(n>=0.and.n<=ndata+1) then
    cinput='rho'; call strio(nf,lh,cinput)
    cinput='u'; call strio(nf,lh,cinput)
    cinput='v'; call strio(nf,lh,cinput)
    cinput='w'; call strio(nf,lh,cinput)
    cinput='p'; call strio(nf,lh,cinput)
 end if
 if(n==ndata+2) then
    cinput='uu'; call strio(nf,lh,cinput)
    cinput='vv'; call strio(nf,lh,cinput)
    cinput='ww'; call strio(nf,lh,cinput)
    cinput='uv'; call strio(nf,lh,cinput)
    cinput='vw'; call strio(nf,lh,cinput)
    cinput='wu'; call strio(nf,lh,cinput)
 end if
 do mm=0,mbk
    write(nf,pos=4*lh+1) 299.0; lh=lh+1 ! Zone Marker
    cinput=czonet(mm); call strio(nf,lh,cinput) ! Zone Name
    write(nf,pos=4*lh+1) -1; lh=lh+1 ! Parent Zone
    write(nf,pos=4*lh+1) int(mm+1,kind=int32); lh=lh+1 ! Strand ID
    write(nf,pos=4*lh+1) real(times(min(max(n,0),ndata)),kind=ieee64); lh=lh+2 ! Solution Time (Double)
    write(nf,pos=4*lh+1) -1; lh=lh+1 ! (Not used. Set to -1.)
    write(nf,pos=4*lh+1) 0; lh=lh+1 ! Zone Type
    write(nf,pos=4*lh+1) 0; lh=lh+1 ! Specify Var Location
    write(nf,pos=4*lh+1) 0; lh=lh+1 ! Raw Local 1-to-1 Face Neighbours Suppliled
    write(nf,pos=4*lh+1) 0; lh=lh+1 ! Number of Miscellaneous Face Neighbour Connections
    write(nf,pos=4*lh+1) int(lximb(mm)+1,kind=int32); lh=lh+1 ! IMax
    write(nf,pos=4*lh+1) int(letmb(mm)+1,kind=int32); lh=lh+1 ! JMax
    write(nf,pos=4*lh+1) int(lzemb(mm)+1,kind=int32); lh=lh+1 ! KMax
    write(nf,pos=4*lh+1) 0; lh=lh+1 ! No Auxillary Data Pairs
 end do
    write(nf,pos=4*lh+1) 357.0; lh=lh+1 ! End of Header Marker
 end if
    write(nf,pos=4*lh+1) 299.0; lh=lh+1 ! Zone Marker
 do m=1,mq
    write(nf,pos=4*lh+1) 1; lh=lh+1 ! 1 = Float / 2 = Double
 end do
    write(nf,pos=4*lh+1) 0; lh=lh+1 ! No Passive Variables
    write(nf,pos=4*lh+1) 0; lh=lh+1 ! No Variable Sharing
    write(nf,pos=4*lh+1) -1; lh=lh+1 ! Zero Based Zone Number to Share
 do m=1,mq; nn=max(3+5*n,0)+m
    write(nf,pos=4*lh+1) real(varmin(nn),kind=ieee64); lh=lh+2 ! Minimum Value (Double) of Variables
    write(nf,pos=4*lh+1) real(varmax(nn),kind=ieee64); lh=lh+2 ! Maximum Value (Double) of Variables
 end do

 end subroutine techead

!===== SUBROUTINE FOR FINDING VARIABLE MIN/MAX VALUES FOR TECPLOT DATA FILE

 subroutine vminmax(nn)

 integer(kind=ni),intent(in) :: nn

    varmin(nn)=minval(varr); varmax(nn)=maxval(varr); varm(:,myid)=(/varmin(nn),varmax(nn)/)

    ir=0; itag=nn
 if(myid==mo(mb)) then
    mps=mo(mb); mpe=mps+nbpc(mb,1)*nbpc(mb,2)*nbpc(mb,3)-1
 do mp=mps+1,mpe
    ir=ir+1; call MPI_IRECV(varm(:,mp),2,MPI_REAL4,mp,itag,icom,ireq(ir),ierr)
 end do
 if(ir/=0) then
    call MPI_WAITALL(ir,ireq,ista,ierr)
 end if
    varmin(nn)=minval(varm(0,mps:mpe)); varmax(nn)=maxval(varm(1,mps:mpe))
 else
    call MPI_SEND(varm(:,myid),2,MPI_REAL4,mo(mb),itag,icom,ierr)
 end if

 end subroutine vminmax

!=====

 end module subroutines3d

!*****