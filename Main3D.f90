!*****
!***** 3D PARALLEL SOLVER
!*****

 program main3d

 use mpi
 use subroutineso
 use subroutines3d
 use problemcase
 implicit none

!===== PREPARATION FOR PARALLEL COMPUTING

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,npro,ierr)

    mpro=npro-1; icom=MPI_COMM_WORLD; info=MPI_INFO_NULL

    allocate(lxim(0:mpro),letm(0:mpro),lzem(0:mpro),lpos(0:mpro),vmpi(0:mpro))

    ll=max(npro,12); allocate(ista(MPI_STATUS_SIZE,ll),ireq(ll))

    inquire(iolength=ll) real(1.0,kind=ieee32); nrecs=ll
    inquire(iolength=ll) real(1.0,kind=ieee64); nrecd=ll

!===== INPUT PARAMETERS

    open(9,file='inputo.dat',status='old')
    read(9,*) cinput,mbk
    read(9,*) cinput,nts
    read(9,*) cinput,nscrn,nsgnl
    read(9,*) cinput,ndata,ndatafl,ndataav
    read(9,*) cinput,nkrk
    read(9,*) cinput,nviscous
    read(9,*) cinput,nsmf
    read(9,*) cinput,nrestart
    read(9,*) cinput,nextrabc,nextgcic
    read(9,*) cinput,nnf(:)
    read(9,*) cinput,reoo,tempoo
    read(9,*) cinput,amach1,amach2,amach3
    read(9,*) cinput,wtemp
    read(9,*) cinput,cfl
    read(9,*) cinput,tmax,timf,tsam
    read(9,*) cinput,fltk,fltrbc
    read(9,*) cinput,dto
    close(9)

    cinput=cinput; fltk=pi*fltk
    amachoo=sqrt(amach1*amach1+amach2*amach2+amach3*amach3)
 if(amachoo>sml) then
    reoo=reoo/amachoo
 end if
    srefoo=111/tempoo; srefp1dre=(srefoo+one)/reoo; sqrtrema=sqrt(reoo*amachoo); sqrtremai=one/max(sqrtrema,sml)
    uoo(:)=(/amach1,amach2,amach3/)

    ll=3+5*(ndata+1); lp=8*mbk+7; lq=12*mbk+11
    allocate(times(0:ndata),cfilet(-1:ndata),ctecplt(-1:ndata),varm(0:1,0:mpro),varmin(ll),varmax(ll))
    allocate(lximb(0:mbk),letmb(0:mbk),lzemb(0:mbk),lhmb(0:mbk),mo(0:mbk),nbpc(0:mbk,3),nbbc(0:mbk,3,0:1),mbcd(0:mbk,3,0:1))
    allocate(imjp(0:lp),jptag(0:lp),jjp(0:lp,0:ljpl))
    allocate(imjl(0:lq),jltag(0:lq),jjl(0:lq,0:ljpl))
    allocate(czonet(0:mbk),cthead(0:mbk))

    call inputext

!===== DOMAIN DECOMPOSITION & BOUNDARY INFORMATION

    mo(0)=0
 do mm=1,mbk
    mo(mm)=mo(mm-1)+nbpc(mm-1,1)*nbpc(mm-1,2)*nbpc(mm-1,3)
 end do
 do mm=0,mbk
    if(myid>=mo(mm)) then; mb=mm; end if
 end do
    lxio=lximb(mb); leto=letmb(mb); lzeo=lzemb(mb)

    cfilet(-1)='grid'
 do n=0,ndata
    no(2)=n/100; no(1)=mod(n,100)/10; no(0)=mod(n,10)
    cno=achar(no+48); cfilet(n)='n'//cno(2)//cno(1)//cno(0)
 end do
 do n=-1,ndata
    ctecplt(n)='data/'//cfilet(n)//'.plt'
 end do
 do mm=0,mbk
    no(2)=mm/100; no(1)=mod(mm,100)/10; no(0)=mod(mm,10)
    cno=achar(no+48); czonet(mm)='z'//cno(2)//cno(1)//cno(0)
    cthead(mm)='data/'//czonet(mm)//'.plt'
 end do
    cgrid='misc/grid'//czonet(mb)//'.dat'; crestart='misc/restart'//czonet(mb)//'.dat'

    no(4)=myid/10000; no(3)=mod(myid,10000)/1000; no(2)=mod(myid,1000)/100; no(1)=mod(myid,100)/10; no(0)=mod(myid,10)
    cno=achar(no+48); cnnode=cno(4)//cno(3)//cno(2)//cno(1)//cno(0)
    cdata='misc/data'//cnnode//'.dat'

    call domdcomp

    ijkp(1)=mod(myid-mo(mb),nbpc(mb,1))
    ijkp(2)=mod((myid-mo(mb))/nbpc(mb,1),nbpc(mb,2))
    ijkp(3)=mod((myid-mo(mb))/(nbpc(mb,1)*nbpc(mb,2)),nbpc(mb,3))

 do nn=1,3; ns=mod(nn,3)+1; ne=mod(ns,3)+1
 do ip=0,1; mm=mbcd(mb,nn,ip)
 if(mm==-1) then
    mmcd(nn,ip)=-1
 else
    mmcd(nn,ip)=idsd3((1-ip)*(nbpc(mm,nn)-1),ijkp(ns),ijkp(ne),mm,nn) 
 end if
 end do
 end do

 do nn=1,3
 select case(nn)
 case (1); ll=lxio; mp=1
 case (2); ll=leto; mp=nbpc(mb,1)
 case (3); ll=lzeo; mp=nbpc(mb,1)*nbpc(mb,2)
 end select
    lp=ijkp(nn); ma=nbpc(mb,nn)
 if(ma==1) then
    l=ll; nbc(nn,0)=nbbc(mb,nn,0); nbc(nn,1)=nbbc(mb,nn,1); mcd(nn,0)=mmcd(nn,0); mcd(nn,1)=mmcd(nn,1)
 end if
 if(ma>=2) then
 if(lp==0) then
    l=ll-((ll+1)/ma)*(ma-1); nbc(nn,0)=nbbc(mb,nn,0); nbc(nn,1)=40; mcd(nn,0)=mmcd(nn,0); mcd(nn,1)=myid+mp
 end if
 if(lp>0.and.lp<ma-1) then
    l=(ll+1)/ma-1; nbc(nn,0)=40; nbc(nn,1)=40; mcd(nn,0)=myid-mp; mcd(nn,1)=myid+mp
 end if
 if(lp==ma-1) then
    l=(ll+1)/ma-1; nbc(nn,0)=40; nbc(nn,1)=nbbc(mb,nn,1); mcd(nn,0)=myid-mp; mcd(nn,1)=mmcd(nn,1)
 end if
 end if
 select case(nn); case (1); lxi=l; case (2); let=l; case (3); lze=l; end select
 end do

!===== SUBDOMAIN SIZES & WRITING START POSITIONS IN OUTPUT FILE

 if(myid==0) then
    lxim(0)=lxi; letm(0)=let; lzem(0)=lze
 do mp=1,mpro
    itag=1; call MPI_RECV(lxim(mp),1,MPI_INTEGER4,mp,itag,icom,ista(:,mp),ierr)
    itag=2; call MPI_RECV(letm(mp),1,MPI_INTEGER4,mp,itag,icom,ista(:,mp),ierr)
    itag=3; call MPI_RECV(lzem(mp),1,MPI_INTEGER4,mp,itag,icom,ista(:,mp),ierr)
 end do
 else
    itag=1; call MPI_SEND(lxi,1,MPI_INTEGER4,0,itag,icom,ierr)
    itag=2; call MPI_SEND(let,1,MPI_INTEGER4,0,itag,icom,ierr)
    itag=3; call MPI_SEND(lze,1,MPI_INTEGER4,0,itag,icom,ierr)
 end if
    call MPI_BCAST(lxim(:),npro,MPI_INTEGER4,0,icom,ierr)
    call MPI_BCAST(letm(:),npro,MPI_INTEGER4,0,icom,ierr)
    call MPI_BCAST(lzem(:),npro,MPI_INTEGER4,0,icom,ierr)

    ltomb=(lxio+1)*(leto+1)*(lzeo+1)

    lmx=(lxi+1)*(let+1)*(lze+1)-1
    lim=(lxi+1)+(let+1)+(lze+1)-1

    ijk(1,1)=lxi; ijk(2,1)=let; ijk(3,1)=lze
    ijk(1,2)=let; ijk(2,2)=lze; ijk(3,2)=lxi
    ijk(1,3)=lze; ijk(2,3)=lxi; ijk(3,3)=let

    nbsize(:)=(ijk(2,:)+1)*(ijk(3,:)+1)

 do mm=0,mbk
    lpos(mo(mm))=0
 do i=1,nbpc(mm,1)-1
    mp=mo(mm)+i
    lpos(mp)=lpos(mp-1)+lxim(mp-1)+1
 end do
    jp=nbpc(mm,1)
 do j=1,nbpc(mm,2)-1; do i=0,nbpc(mm,1)-1
    mp=mo(mm)+j*jp+i
    lpos(mp)=lpos(mp-jp)+(lximb(mm)+1)*(letm(mp-jp)+1)
 end do; end do
    kp=nbpc(mm,1)*nbpc(mm,2)
 do k=1,nbpc(mm,3)-1; do j=0,nbpc(mm,2)-1; do i=0,nbpc(mm,1)-1
    mp=mo(mm)+k*kp+j*jp+i
    lpos(mp)=lpos(mp-kp)+(lximb(mm)+1)*(letmb(mm)+1)*(lzem(mp-kp)+1)
 end do; end do; end do
 end do

!===== ALLOCATION OF MAIN ARRAYS

    allocate(qo(0:lmx,5),qa(0:lmx,5),qb(0:lmx,5),de(0:lmx,5))
    allocate(xim(0:lmx,3),etm(0:lmx,3),zem(0:lmx,3),rr(0:lmx,3),ss(0:lmx,3))
    allocate(p(0:lmx),yaco(0:lmx),varr(0:lmx),nrr(0:lmx),npex(0:lmx))

 if(nviscous==1) then
    allocate(txx(0:lmx),tyy(0:lmx),tzz(0:lmx))
    allocate(txy(0:lmx),tyz(0:lmx),tzx(0:lmx))
    allocate(hxx(0:lmx),hyy(0:lmx),hzz(0:lmx))
 end if

    ii=nbsize(1)-1; jj=nbsize(2)-1; kk=nbsize(3)-1
    allocate(drva1(0:ii,5,0:1),drva2(0:jj,5,0:1),drva3(0:kk,5,0:1))
    allocate(drvb1(0:ii,5,0:1),drvb2(0:jj,5,0:1),drvb3(0:kk,5,0:1))
    allocate(send01(0:ii,0:1,0:1),send02(0:jj,0:1,0:1),send03(0:kk,0:1,0:1))
    allocate(recv01(0:ii,0:1,0:1),recv02(0:jj,0:1,0:1),recv03(0:kk,0:1,0:1))
    allocate(send11(0:ii,0:2,0:1),send12(0:jj,0:2,0:1),send13(0:kk,0:2,0:1))
    allocate(recv11(0:ii,0:2,0:1),recv12(0:jj,0:2,0:1),recv13(0:kk,0:2,0:1))
    allocate(cm1(0:ii,3,0:1),cm2(0:jj,3,0:1),cm3(0:kk,3,0:1))

    allocate(xu(0:lim,3),yu(0:lim,3),xl(0:lim,2),yl(0:lim,2),li(0:lim),sa(0:lim),sb(0:lim))

!===== EXTRA COEFFICIENTS FOR DOMAIN BOUNDARIES

    call fcbcm(fltk,fltrbc)
    call fcint(fltk,half,alphf,betf,fa,fb,fc)
    albef(:,0,1)=(/zero,zero,one,alphf,betf/)
    albef(:,1,1)=(/zero,alphf,one,alphf,betf/)
    albef(:,2,1)=(/betf,alphf,one,alphf,betf/)

    pbco(:,:,:)=zero; pbci(:,:,:)=zero; call sbcco
 do nt=0,1; do j=0,1; ii=lmd+nt*(lmf-lmd)
    pbcot(j,nt)=sum(pbco(0:ii,j,nt))
 end do; end do

!===== PENTADIAGONAL MATRICES FOR DIFFERENCING & FILETERING

 do nn=1,3
 select case(nn)
 case(1); is=0; ie=is+lxi; case(2); is=lxi+1; ie=is+let; case(3); is=lxi+let+2; ie=is+lze
 end select
 do ip=0,1; np=nbc(nn,ip)
 select case(np)
 case(10,20,25,30); ndf(nn,ip,0)=0; ndf(nn,ip,1)=0
 case(35,40,45); ndf(nn,ip,0)=1; ndf(nn,ip,1)=1
 end select
 end do
    ns=ndf(nn,0,0); ne=ndf(nn,1,0)
    call penta(xu(:,:),xl(:,:),is,ie,ns,ne,0)
    ns=ndf(nn,0,1); ne=ndf(nn,1,1)
    call penta(yu(:,:),yl(:,:),is,ie,ns,ne,1)
 end do

!===== GRID INPUT & CALCULATION OF GRID METRICS

    allocate(lio(0:let,0:lze))
 do k=0,lze; kp=k*(leto+1)*(lxio+1)
 do j=0,let; jp=j*(lxio+1)
    lio(j,k)=jp+kp
 end do
 end do
    call makegrid
    call MPI_BARRIER(icom,ierr)

    open(9,file=cgrid,access='direct',form='unformatted',recl=3*nrecd,status='old')
    lp=lpos(myid)
 do k=0,lze; do j=0,let; lq=lp+lio(j,k)
 do i=0,lxi; l=indx3(i,j,k,1)
    read(9,rec=lq+i+1) ss(l,:)
 end do
 end do; end do
    close(9)
    call MPI_BARRIER(icom,ierr)
 if(myid==mo(mb)) then
    open(9,file=cgrid,status='old'); close(9,status='delete')
 end if

    rr(:,1)=ss(:,1)
    m=1; call mpigo(0,nrone,n45go,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    qo(:,1)=rr(:,1); qo(:,2)=rr(:,2); qo(:,3)=rr(:,3)

    rr(:,1)=ss(:,2)
    m=2; call mpigo(0,nrone,n45go,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    qa(:,1)=rr(:,1); qa(:,2)=rr(:,2); qa(:,3)=rr(:,3)

    rr(:,1)=ss(:,3)
    m=3; call mpigo(0,nrone,n45go,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    de(:,1)=rr(:,1); de(:,2)=rr(:,2); de(:,3)=rr(:,3)

    xim(:,1)=qa(:,2)*de(:,3)-de(:,2)*qa(:,3)
    xim(:,2)=de(:,2)*qo(:,3)-qo(:,2)*de(:,3)
    xim(:,3)=qo(:,2)*qa(:,3)-qa(:,2)*qo(:,3)
    etm(:,1)=qa(:,3)*de(:,1)-de(:,3)*qa(:,1)
    etm(:,2)=de(:,3)*qo(:,1)-qo(:,3)*de(:,1)
    etm(:,3)=qo(:,3)*qa(:,1)-qa(:,3)*qo(:,1)
    zem(:,1)=qa(:,1)*de(:,2)-de(:,1)*qa(:,2)
    zem(:,2)=de(:,1)*qo(:,2)-qo(:,1)*de(:,2)
    zem(:,3)=qo(:,1)*qa(:,2)-qa(:,1)*qo(:,2)

!    rr(:,3)=qa(:,2)*ss(:,3); rr(:,2)=qa(:,3)*ss(:,3)
!    m=1; call mpigo(0,nrall,n45go,m); call deriv(3,3,m); call deriv(2,2,m); xim(:,m)=rr(:,3)-rr(:,2)
!    rr(:,3)=de(:,2)*ss(:,1); rr(:,2)=de(:,3)*ss(:,1)
!    m=2; call mpigo(0,nrall,n45go,m); call deriv(3,3,m); call deriv(2,2,m); xim(:,m)=rr(:,3)-rr(:,2)
!    rr(:,3)=qo(:,2)*ss(:,2); rr(:,2)=qo(:,3)*ss(:,2)
!    m=3; call mpigo(0,nrall,n45go,m); call deriv(3,3,m); call deriv(2,2,m); xim(:,m)=rr(:,3)-rr(:,2)
!
!    rr(:,1)=qa(:,3)*ss(:,3); rr(:,3)=qa(:,1)*ss(:,3)
!    m=1; call mpigo(0,nrall,n45go,m); call deriv(1,1,m); call deriv(3,3,m); etm(:,m)=rr(:,1)-rr(:,3)
!    rr(:,1)=de(:,3)*ss(:,1); rr(:,3)=de(:,1)*ss(:,1)
!    m=2; call mpigo(0,nrall,n45go,m); call deriv(1,1,m); call deriv(3,3,m); etm(:,m)=rr(:,1)-rr(:,3)
!    rr(:,1)=qo(:,3)*ss(:,2); rr(:,3)=qo(:,1)*ss(:,2)
!    m=3; call mpigo(0,nrall,n45go,m); call deriv(1,1,m); call deriv(3,3,m); etm(:,m)=rr(:,1)-rr(:,3)
!
!    rr(:,2)=qa(:,1)*ss(:,3); rr(:,1)=qa(:,2)*ss(:,3)
!    m=1; call mpigo(0,nrall,n45go,m); call deriv(2,2,m); call deriv(1,1,m); zem(:,m)=rr(:,2)-rr(:,1)
!    rr(:,2)=de(:,1)*ss(:,1); rr(:,1)=de(:,2)*ss(:,1)
!    m=2; call mpigo(0,nrall,n45go,m); call deriv(2,2,m); call deriv(1,1,m); zem(:,m)=rr(:,2)-rr(:,1)
!    rr(:,2)=qo(:,1)*ss(:,2); rr(:,1)=qo(:,2)*ss(:,2)
!    m=3; call mpigo(0,nrall,n45go,m); call deriv(2,2,m); call deriv(1,1,m); zem(:,m)=rr(:,2)-rr(:,1)

 do m=1,3
    rr(:,1)=xim(:,m)
    call mpigo(1,nrone,n45no,9*(m-1)+1); call filte(nnf(1),1)
    call mpigo(1,nrone,n45no,9*(m-1)+2); call filte(nnf(2),1)
    call mpigo(1,nrone,n45no,9*(m-1)+3); call filte(nnf(3),1)
    xim(:,m)=rr(:,1)
    rr(:,1)=etm(:,m)
    call mpigo(1,nrone,n45no,9*(m-1)+4); call filte(nnf(1),1)
    call mpigo(1,nrone,n45no,9*(m-1)+5); call filte(nnf(2),1)
    call mpigo(1,nrone,n45no,9*(m-1)+6); call filte(nnf(3),1)
    etm(:,m)=rr(:,1)
    rr(:,1)=zem(:,m)
    call mpigo(1,nrone,n45no,9*(m-1)+7); call filte(nnf(1),1)
    call mpigo(1,nrone,n45no,9*(m-1)+8); call filte(nnf(2),1)
    call mpigo(1,nrone,n45no,9*(m-1)+9); call filte(nnf(3),1)
    zem(:,m)=rr(:,1)
 end do
    yaco(:)=three/(qo(:,1)*xim(:,1)+qo(:,2)*etm(:,1)+qo(:,3)*zem(:,1)&
                  +qa(:,1)*xim(:,2)+qa(:,2)*etm(:,2)+qa(:,3)*zem(:,2)&
                  +de(:,1)*xim(:,3)+de(:,2)*etm(:,3)+de(:,3)*zem(:,3))

 do nn=1,3; do ip=0,1; i=ip*ijk(1,nn)
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
 select case(nn)
 case(1); rv(:)=yaco(l)*xim(l,:); fctr=one/sqrt(rv(1)*rv(1)+rv(2)*rv(2)+rv(3)*rv(3)); cm1(jk,:,ip)=fctr*rv(:)
 case(2); rv(:)=yaco(l)*etm(l,:); fctr=one/sqrt(rv(1)*rv(1)+rv(2)*rv(2)+rv(3)*rv(3)); cm2(jk,:,ip)=fctr*rv(:)
 case(3); rv(:)=yaco(l)*zem(l,:); fctr=one/sqrt(rv(1)*rv(1)+rv(2)*rv(2)+rv(3)*rv(3)); cm3(jk,:,ip)=fctr*rv(:)
 end select
 end do
 end do
 end do; end do

!===== EXTRA COEFFICIENTS FOR GCBC/GCIC

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
    call extrabcc(npex(l))
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

!===== POINT JUNCTION SEARCH

    lp=8*mbk+7; imjp(:)=0; jjp(:,:)=-1
 do ip=0,1
    ijp(:,1,ip)=(/0,2,4,6/)+ip; ijp(:,2,ip)=(/0,4,1,5/)+2*ip; ijp(:,3,ip)=(/0,1,2,3/)+4*ip
 end do
 do i=0,lp
    jjp(i,0)=i
 end do
 do mm=0,mbk; do nn=1,3; do ip=0,1; mp=mbcd(mm,nn,ip)
 if(mp/=-1) then
 do i=0,3; is=8*mm+ijp(i,nn,ip); ie=8*mp+ijp(i,nn,1-ip)
    imjp(is)=imjp(is)+1; jjp(is,imjp(is))=jjp(ie,0)
 end do
 end if
 end do; end do; end do
 do i=0,lp; do j=0,lp; ll=0
 if((i-j)*imjp(i)*imjp(j)/=0) then
 do k=0,imjp(i); do kk=0,imjp(j)
 if(jjp(j,kk)-jjp(i,k)==0) then; ll=ll+1; end if
 end do; end do
 if(ll/=0) then
    ks=imjp(i)+1; ke=ks+imjp(j); jjp(i,ks:ke)=jjp(j,0:ke-ks); call trimm(jjp(i,:),imjp(i),ke)
 end if
 end if
 end do; end do
 do i=0,lp
    jptag(i)=minval(jjp(i,0:imjp(i)))
 end do
    njp(:)=-1; jpcd(:,:)=-1
 do i=0,7
    ip=mod(i,2); jp=mod(i,4)/2; kp=i/4; l=8*mb+i
 if(ijkp(1)==ip*(nbpc(mb,1)-1).and.ijkp(2)==jp*(nbpc(mb,2)-1).and.ijkp(3)==kp*(nbpc(mb,3)-1).and.imjp(l)>=3) then; njp(i)=l
 do k=1,imjp(l)
    mm=jjp(l,k)/8; j=jjp(l,k)-8*mm; ii=mod(j,2); jj=mod(j,4)/2; kk=j/4
    jpcd(i,k)=idsd3(ii*(nbpc(mm,1)-1),jj*(nbpc(mm,2)-1),kk*(nbpc(mm,3)-1),mm,1)
 end do
 end if
 end do

!===== LINE JUNCTION SEARCH

    lp=12*mbk+11; imjl(:)=0; jjl(:,:)=-1; njl(0:3)=(/2,2,1,1/)
 do ip=0,1
    ijl(:,1,ip)=(/0,1,8,10/)+ip*njl(0:3); ijl(:,2,ip)=(/4,5,0,2/)+ip*njl(0:3); ijl(:,3,ip)=(/8,9,4,6/)+ip*njl(0:3)
 end do
 do i=0,lp
    jjl(i,0)=i
 end do
 do mm=0,mbk; do nn=1,3; do ip=0,1; mp=mbcd(mm,nn,ip)
 if(mp/=-1) then
 do i=0,3; is=12*mm+ijl(i,nn,ip); ie=12*mp+ijl(i,nn,1-ip)
    imjl(is)=imjl(is)+1; jjl(is,imjl(is))=jjl(ie,0)
 end do
 end if
 end do; end do; end do
 do i=0,lp; do j=0,lp; ll=0
 if((i-j)*imjl(i)*imjl(j)/=0) then
 do k=0,imjl(i); do kk=0,imjl(j)
 if(jjl(j,kk)-jjl(i,k)==0) then; ll=ll+1; end if
 end do; end do
 if(ll/=0) then
    ks=imjl(i)+1; ke=ks+imjl(j); jjl(i,ks:ke)=jjl(j,0:ke-ks); call trimm(jjl(i,:),imjl(i),ke)
 end if
 end if
 end do; end do
 do i=0,lp
    jltag(i)=minval(jjl(i,0:imjl(i)))
 end do
    njl(:)=-1; jlcd(:,:)=-1
 do i=0,11
    nn=i/4+1; ns=mod(nn,3)+1; ne=mod(ns,3)+1; ip=mod(i,4)/2; jp=mod(mod(i,4),2); l=12*mb+i
 if(ijkp(nn)==ip*(nbpc(mb,nn)-1).and.ijkp(ns)==jp*(nbpc(mb,ns)-1).and.imjl(l)>=2) then; njl(i)=l
 do k=1,imjl(l)
    mm=jjl(l,k)/12; j=jjl(l,k)-12*mm; np=j/4+1; nq=mod(np,3)+1; ii=mod(j,4)/2; jj=mod(mod(j,4),2)
    jlcd(i,k)=idsd3(ii*(nbpc(mm,np)-1),jj*(nbpc(mm,nq)-1),ijkp(ne),mm,np)
 end do
 end if
 end do

!===== SETTING UP OUTPUT FILE & STORING GRID DATA

    open(0,file=cdata,status='unknown'); close(0,status='delete') ! 'replace' not suitable as 'recl' may vary
    open(0,file=cdata,access='direct',form='unformatted',recl=nrecs*(lmx+1),status='new')
 do nn=1,3
    varr(:)=ss(:,nn); write(0,rec=nn) varr(:); call vminmax(nn)
 end do
    close(0)

!===== SETTING UP SPONGE ZONE PARAMETERS

    call spongeup

!===== INITIAL CONDITIONS

 if(nts==0) then
    n=0; ndt=0; dt=zero; dts=zero; dte=zero; timo=zero
    call initialo
 else
    open(9,file=crestart,access='direct',form='unformatted',recl=5*nrecd,status='old')
    read(9,rec=1) cha(:); read(9,rec=2) dha(:)
    n=cha(1); ndt=cha(2); dt=cha(3); dts=cha(4); dte=cha(5); timo=dha(1); lp=lpos(myid)+2
 do k=0,lze; do j=0,let; lq=lp+lio(j,k)
 do i=0,lxi; l=indx3(i,j,k,1)
    read(9,rec=lq+i+1) qa(l,:)
 end do
 end do; end do
    close(9)
 end if
    qb(:,:)=zero

!============================================
!===== BEGINNING OF TIME MARCHING IN SOLUTION
!============================================

 if(myid==0) then
    open(1,file='signal.dat',access='direct',form='formatted',recl=16,status='replace'); close(1)
 end if
    call MPI_BARRIER(icom,ierr)

    wts=MPI_WTIME()

    ndati=-1; nsigi=-1; dtsum=zero
 do while(timo<tmax.and.(dt/=zero.or.n<=2))

 if(myid**2+mod(n,nscrn)**2==0) then
    write(*,"(' n =',i8,'   time =',f12.5)") n,timo
 end if

!----- FILTERING & RE-INITIALISING

 do m=1,5
    rr(:,1)=qa(:,m)
    call mpigo(1,nrone,n45no,3*(m-1)+1); call filte(nnf(1),1)
    call mpigo(1,nrone,n45no,3*(m-1)+2); call filte(nnf(2),1)
    call mpigo(1,nrone,n45no,3*(m-1)+3); call filte(nnf(3),1)
    qa(:,m)=rr(:,1)
 end do
    qo(:,:)=qa(:,:)

!-------------------------------------
!----- BEGINNING OF RUNGE-KUTTA STAGES
!-------------------------------------

 do nk=1,nkrk

!----- MOVING FRAME VELOCITY & ACCELERATION BEFORE TIME ADVANCING

    dtko=dt*min(max(nk-2,0),1)/(nkrk-nk+3); dtk=dt*min(nk-1,1)/(nkrk-nk+2)
    call movef(dtko,dtk)

!----- TEMPORARY STORAGE OF PRIMITIVE VARIABLES & PRESSURE

    de(:,1)=one/qa(:,1)
    de(:,2)=qa(:,2)*de(:,1)
    de(:,3)=qa(:,3)*de(:,1)
    de(:,4)=qa(:,4)*de(:,1)

    p(:)=gamm1*(qa(:,5)-half*(qa(:,2)*de(:,2)+qa(:,3)*de(:,3)+qa(:,4)*de(:,4)))
    de(:,5)=gam*p(:)*de(:,1)
    ss(:,1)=srefp1dre*de(:,5)**1.5_nr/(de(:,5)+srefoo)

!----- DETERMINATION OF TIME STEP SIZE & OUTPUT TIME

 if(nk==1) then
 if(mod(n,10)==1) then; ndt=n; dts=dte
 if(dto<zero) then
    rr(:,1)=xim(:,1)*xim(:,1)+xim(:,2)*xim(:,2)+xim(:,3)*xim(:,3)&
           +etm(:,1)*etm(:,1)+etm(:,2)*etm(:,2)+etm(:,3)*etm(:,3)&
           +zem(:,1)*zem(:,1)+zem(:,2)*zem(:,2)+zem(:,3)*zem(:,3)
    rr(:,2)=abs(xim(:,1)*(de(:,2)+umf(1))+xim(:,2)*(de(:,3)+umf(2))+xim(:,3)*(de(:,4)+umf(3)))&
           +abs(etm(:,1)*(de(:,2)+umf(1))+etm(:,2)*(de(:,3)+umf(2))+etm(:,3)*(de(:,4)+umf(3)))&
           +abs(zem(:,1)*(de(:,2)+umf(1))+zem(:,2)*(de(:,3)+umf(2))+zem(:,3)*(de(:,4)+umf(3)))
    ss(:,2)=abs(yaco(:))
    res=maxval((sqrt(de(:,5)*rr(:,1))+rr(:,2))*ss(:,2))
    call MPI_ALLREDUCE(res,fctr,1,MPI_REAL8,MPI_MAX,icom,ierr)
    ra0=cfl/fctr; ra1=ra0
 if(nviscous==1) then
    res=maxval(de(:,1)*ss(:,1)*rr(:,1)*ss(:,2)*ss(:,2))
    call MPI_ALLREDUCE(res,fctr,1,MPI_REAL8,MPI_MAX,icom,ierr)
    ra1=half/fctr
 end if
    dte=min(ra0,ra1)
 else
    dte=dto
 end if
 end if
    dt=dts+(dte-dts)*sin(0.05_nr*pi*(n-ndt))**two

    nout=0; res=tsam+(ndati+1)*(tmax-tsam)/ndata
 if((timo-res)*(timo+dt-res)<=zero) then
    nout=1; ndati=ndati+1
 end if
 end if

!----- VISCOUS SHEAR STRESSES & HEAT FLUXES

 if(nviscous==1) then
    de(:,1)=ss(:,1)

    rr(:,1)=de(:,2)
    m=2; call mpigo(0,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    txx(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
    hzz(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
    tzx(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

    rr(:,1)=de(:,3)
    m=3; call mpigo(0,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    txy(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
    tyy(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
    hxx(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

    rr(:,1)=de(:,4)
    m=4; call mpigo(0,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    hyy(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
    tyz(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
    tzz(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

    rr(:,1)=de(:,5)
    m=5; call mpigo(0,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    ss(:,1)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
    ss(:,2)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
    ss(:,3)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

    rr(:,1)=de(:,1)*yaco(:)
    rr(:,2)=gamm1prndtli*rr(:,1)
    de(:,5)=twothirds*(txx(:)+tyy(:)+tzz(:))

    txx(:)=rr(:,1)*(txx(:)+txx(:)-de(:,5))
    tyy(:)=rr(:,1)*(tyy(:)+tyy(:)-de(:,5))
    tzz(:)=rr(:,1)*(tzz(:)+tzz(:)-de(:,5))
    txy(:)=rr(:,1)*(txy(:)+hzz(:))
    tyz(:)=rr(:,1)*(tyz(:)+hxx(:))
    tzx(:)=rr(:,1)*(tzx(:)+hyy(:))
    hxx(:)=rr(:,2)*ss(:,1)+de(:,2)*txx(:)+de(:,3)*txy(:)+de(:,4)*tzx(:)
    hyy(:)=rr(:,2)*ss(:,2)+de(:,2)*txy(:)+de(:,3)*tyy(:)+de(:,4)*tyz(:)
    hzz(:)=rr(:,2)*ss(:,3)+de(:,2)*tzx(:)+de(:,3)*tyz(:)+de(:,4)*tzz(:)
 end if

!----- CALCULATION OF FLUX DERIVATIVES

    rr(:,1)=de(:,2)+umf(1)
    rr(:,2)=de(:,3)+umf(2)
    rr(:,3)=de(:,4)+umf(3)
    ss(:,1)=xim(:,1)*rr(:,1)+xim(:,2)*rr(:,2)+xim(:,3)*rr(:,3)
    ss(:,2)=etm(:,1)*rr(:,1)+etm(:,2)*rr(:,2)+etm(:,3)*rr(:,3)
    ss(:,3)=zem(:,1)*rr(:,1)+zem(:,2)*rr(:,2)+zem(:,3)*rr(:,3)

    rr(:,1)=qa(:,1)*ss(:,1)
    rr(:,2)=qa(:,1)*ss(:,2)
    rr(:,3)=qa(:,1)*ss(:,3)
    m=1; call mpigo(0,nrall,n45no,m); call deriv(1,1,m); call deriv(2,2,m); call deriv(3,3,m)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    rr(:,1)=qa(:,2)*ss(:,1)+xim(:,1)*p(:)
    rr(:,2)=qa(:,2)*ss(:,2)+etm(:,1)*p(:)
    rr(:,3)=qa(:,2)*ss(:,3)+zem(:,1)*p(:)
 if(nviscous==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*txx(:)-xim(:,2)*txy(:)-xim(:,3)*tzx(:)
    rr(:,2)=rr(:,2)-etm(:,1)*txx(:)-etm(:,2)*txy(:)-etm(:,3)*tzx(:)
    rr(:,3)=rr(:,3)-zem(:,1)*txx(:)-zem(:,2)*txy(:)-zem(:,3)*tzx(:)
 end if
    m=2; call mpigo(0,nrall,n45no,m); call deriv(1,1,m); call deriv(2,2,m); call deriv(3,3,m)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    rr(:,1)=qa(:,3)*ss(:,1)+xim(:,2)*p(:)
    rr(:,2)=qa(:,3)*ss(:,2)+etm(:,2)*p(:)
    rr(:,3)=qa(:,3)*ss(:,3)+zem(:,2)*p(:)
 if(nviscous==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*txy(:)-xim(:,2)*tyy(:)-xim(:,3)*tyz(:)
    rr(:,2)=rr(:,2)-etm(:,1)*txy(:)-etm(:,2)*tyy(:)-etm(:,3)*tyz(:)
    rr(:,3)=rr(:,3)-zem(:,1)*txy(:)-zem(:,2)*tyy(:)-zem(:,3)*tyz(:)
 end if
    m=3; call mpigo(0,nrall,n45no,m); call deriv(1,1,m); call deriv(2,2,m); call deriv(3,3,m)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    rr(:,1)=qa(:,4)*ss(:,1)+xim(:,3)*p(:)
    rr(:,2)=qa(:,4)*ss(:,2)+etm(:,3)*p(:)
    rr(:,3)=qa(:,4)*ss(:,3)+zem(:,3)*p(:)
 if(nviscous==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*tzx(:)-xim(:,2)*tyz(:)-xim(:,3)*tzz(:)
    rr(:,2)=rr(:,2)-etm(:,1)*tzx(:)-etm(:,2)*tyz(:)-etm(:,3)*tzz(:)
    rr(:,3)=rr(:,3)-zem(:,1)*tzx(:)-zem(:,2)*tyz(:)-zem(:,3)*tzz(:)
 end if
    m=4; call mpigo(0,nrall,n45no,m); call deriv(1,1,m); call deriv(2,2,m); call deriv(3,3,m)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    de(:,5)=qa(:,5)+p(:)
    rr(:,1)=de(:,5)*ss(:,1)-p(:)*(umf(1)*xim(:,1)+umf(2)*xim(:,2)+umf(3)*xim(:,3))
    rr(:,2)=de(:,5)*ss(:,2)-p(:)*(umf(1)*etm(:,1)+umf(2)*etm(:,2)+umf(3)*etm(:,3))
    rr(:,3)=de(:,5)*ss(:,3)-p(:)*(umf(1)*zem(:,1)+umf(2)*zem(:,2)+umf(3)*zem(:,3))
 if(nviscous==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*hxx(:)-xim(:,2)*hyy(:)-xim(:,3)*hzz(:)
    rr(:,2)=rr(:,2)-etm(:,1)*hxx(:)-etm(:,2)*hyy(:)-etm(:,3)*hzz(:)
    rr(:,3)=rr(:,3)-zem(:,1)*hxx(:)-zem(:,2)*hyy(:)-zem(:,3)*hzz(:)
 end if
    m=5; call mpigo(0,nrall,n45no,m); call deriv(1,1,m); call deriv(2,2,m); call deriv(3,3,m)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

!----- PREPARATION FOR GCBC & GCIC

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

!----- INTERNODE COMMNICATION FOR GCIC

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

!----- IMPLEMENTATION OF GCBC & GCIC

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
 case(1); call extrabcs
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

!----- IMPLEMENTATION OF SPONGE CONDITION

    call spongego

!----- UPDATING CONSERVATIVE VARIABLES

    dtko=dt*min(nk-1,1)/(nkrk-nk+2); dtk=dt/(nkrk-nk+1)
    call movef(dtko,dtk)

    rr(:,1)=dtk*yaco(:)
    qa(:,1)=qo(:,1)-rr(:,1)*de(:,1)
    qa(:,2)=qo(:,2)-rr(:,1)*de(:,2)
    qa(:,3)=qo(:,3)-rr(:,1)*de(:,3)
    qa(:,4)=qo(:,4)-rr(:,1)*de(:,4)
    qa(:,5)=qo(:,5)-rr(:,1)*de(:,5)

    call extracon

!----- WALL TEMPERATURE/VELOCITY CONDITION

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

!----- POINT JUNCTION AVERAGING

    ir=0; lp=-5; lq=-5
 do i=0,7; ii=njp(i)
 if(ii/=-1) then; itag=jptag(ii); ip=mod(i,2); jp=mod(i,4)/2; kp=i/4; lp=lp+5
    l=indx3(ip*lxi,jp*let,kp*lze,1); rr(lp:lp+4,1)=qa(l,:)
 do j=1,imjp(ii); lq=lq+5
    ir=ir+1; call MPI_ISEND(rr(lp:lp+4,1),5,MPI_REAL8,jpcd(i,j),itag,icom,ireq(ir),ierr)
    ir=ir+1; call MPI_IRECV(rr(lq:lq+4,2),5,MPI_REAL8,jpcd(i,j),itag,icom,ireq(ir),ierr)
 end do
 end if
 end do
 if(ir/=0) then
    call MPI_WAITALL(ir,ireq,ista,ierr)
 end if
    lp=-5; lq=-5
 do i=0,7; ii=njp(i)
 if(ii/=-1) then; ip=mod(i,2); jp=mod(i,4)/2; kp=i/4; lp=lp+5
 do j=1,imjp(ii); lq=lq+5
    rr(lp:lp+4,1)=rr(lp:lp+4,1)+rr(lq:lq+4,2)
 end do
    fctr=one/(imjp(ii)+1); l=indx3(ip*lxi,jp*let,kp*lze,1); qa(l,:)=fctr*rr(lp:lp+4,1)
 end if
 end do

!----- LINE JUNCTION AVERAGING

    ir=0; lp=0; lq=0
 do i=0,11; ii=njl(i)
 if(ii/=-1) then; itag=jltag(ii)
    nn=i/4+1; ip=mod(i,4)/2; jp=mod(mod(i,4),2); ks=lp; ke=ks+5*(ijk(3,nn)+1)-1; lp=lp+ke-ks+1
 do k=0,ijk(3,nn); ll=ks+5*k
    l=indx3(ip*ijk(1,nn),jp*ijk(2,nn),k,nn); rr(ll:ll+4,1)=qa(l,:)
 end do
 do j=1,imjl(ii); js=lq; je=js+ke-ks; lq=lq+je-js+1
    ir=ir+1; call MPI_ISEND(rr(ks:ke,1),ke-ks+1,MPI_REAL8,jlcd(i,j),itag,icom,ireq(ir),ierr)
    ir=ir+1; call MPI_IRECV(rr(js:je,2),je-js+1,MPI_REAL8,jlcd(i,j),itag,icom,ireq(ir),ierr)
 end do
 end if
 end do
 if(ir/=0) then
    call MPI_WAITALL(ir,ireq,ista,ierr)
 end if
    lp=0; lq=0
 do i=0,11; ii=njl(i)
 if(ii/=-1) then
    nn=i/4+1; ip=mod(i,4)/2; jp=mod(mod(i,4),2); ks=lp; ke=ks+5*(ijk(3,nn)+1)-1; lp=lp+ke-ks+1
 do j=1,imjl(ii); js=lq; je=js+ke-ks; lq=lq+je-js+1
    rr(ks:ke,1)=rr(ks:ke,1)+rr(js:je,2)
 end do
    fctr=one/(imjl(ii)+1)
 do k=0,ijk(3,nn); ll=ks+5*k
    l=indx3(ip*ijk(1,nn),jp*ijk(2,nn),k,nn); qa(l,:)=fctr*rr(ll:ll+4,1)
 end do
 end if
 end do

!----- INTERFACE SURFACE AVERAGING

    ir=0; itag=30
 do nn=1,3
 select case(nn)
 case(1); drva=>drva1; drvb=>drvb1; case(2); drva=>drva2; drvb=>drvb2; case(3); drva=>drva3; drvb=>drvb3
 end select
 do ip=0,1; iq=1-ip; np=nbc(nn,ip); i=ip*ijk(1,nn)
 if((np-30)*(np-35)*(np-45)*(abs(nextgcic-1)+abs((np-20)*(np-25)))==0) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); l=indx3(i,j,k,nn); jk=kp+j
    drva(jk,:,ip)=qa(l,:); rr(l,1)=one
 end do
 end do
    ir=ir+1; call MPI_ISEND(drva(:,:,ip),5*nbsize(nn),MPI_REAL8,mcd(nn,ip),itag+iq,icom,ireq(ir),ierr)
    ir=ir+1; call MPI_IRECV(drvb(:,:,ip),5*nbsize(nn),MPI_REAL8,mcd(nn,ip),itag+ip,icom,ireq(ir),ierr)
 end if
 end do
 end do
 if(ir/=0) then
    call MPI_WAITALL(ir,ireq,ista,ierr)
 end if
 do nn=1,3
 select case(nn); case(1); drvb=>drvb1; case(2); drvb=>drvb2; case(3); drvb=>drvb3; end select
 do ip=0,1; np=nbc(nn,ip); i=ip*ijk(1,nn)
 if((np-30)*(np-35)*(np-45)*(abs(nextgcic-1)+abs((np-20)*(np-25)))==0) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); l=indx3(i,j,k,nn); jk=kp+j
    rr(l,1)=rr(l,1)+nrr(l); rr(l,2)=nrr(l)/rr(l,1)
    qa(l,:)=rr(l,2)*((rr(l,1)-one)*qa(l,:)+drvb(jk,:,ip))+(1-nrr(l))*qa(l,:)
 end do
 end do
 end if
 end do
 end do

!-------------------------------
!----- END OF RUNGE-KUTTA STAGES
!-------------------------------

 end do

!---------------------
!----- ADVANCE IN TIME
!---------------------

    n=n+1
    timo=timo+dt

!----- RECORDING INTERMEDIATE RESULTS

 if(timo>tsam-(tmax-tsam)/ndata) then
    dtsum=dtsum+dt; fctr=half*dt; qb(:,:)=qb(:,:)+fctr*(qo(:,:)+qa(:,:))
 if(nout==1) then
    times(ndati)=timo-half*dtsum
    open(0,file=cdata,access='direct',form='unformatted',recl=nrecs*(lmx+1),status='old')
 if(n==1) then
    qb(:,:)=qo(:,:)
 else
    ra0=ndataav/dtsum; ra1=one-ndataav; qb(:,:)=ra0*qb(:,:)+ra1*qa(:,:)
 end if
    rr(:,1)=one/qb(:,1)
 do m=1,5
 select case(m)
 case(1); varr(:)=qb(:,m); case(2:4); varr(:)=rr(:,1)*qb(:,m)+umf(m-1)
 case(5); varr(:)=gamm1*(qb(:,m)-half*rr(:,1)*(qb(:,2)*qb(:,2)+qb(:,3)*qb(:,3)+qb(:,4)*qb(:,4)))
 end select
    nn=3+5*ndati+m; write(0,rec=nn) varr(:); call vminmax(nn)
 end do
    close(0)
    dtsum=zero; qb(:,:)=zero
 end if
 end if

 if(timo>=tsam.and.mod(n,nsgnl)==0) then
    nsigi=nsigi+1; call signalgo
 end if

!==========================
!===== END OF TIME MARCHING
!==========================

 end do

    wte=MPI_WTIME(); res=wte-wts
    call MPI_ALLREDUCE(res,wtime,1,MPI_REAL8,MPI_SUM,icom,ierr)
 if(myid==0) then
    open(9,file='timeouts.dat',status='replace')
    write(9,'(es15.7)') times(:)
    close(9)
    open(9,file='walltime.dat',position='append',status='unknown')
    write(9,'(2es15.7)') real(npro,kind=nr),wtime/npro
    close(9)
 end if

!===== GENERATING RESTART DATA FILE

 if(nrestart==1) then
 if(myid==mo(mb)) then
    open(9,file=crestart,access='direct',form='unformatted',recl=5*nrecd,status='replace')
    cha(:)=(/real(n,kind=nr),real(ndt,kind=nr),dt,dts,dte/); dha(:)=(/timo,zero,zero,zero,zero/)
    write(9,rec=1) cha(:); write(9,rec=2) dha(:)
    close(9)
 end if
    call MPI_BARRIER(icom,ierr)
    open(9,file=crestart,access='direct',form='unformatted',recl=5*nrecd,status='old')
    lp=lpos(myid)+2
 do k=0,lze; do j=0,let; lq=lp+lio(j,k)
 do i=0,lxi; l=indx3(i,j,k,1)
    write(9,rec=lq+i+1) qa(l,:)
 end do
 end do; end do
    close(9)
 end if

!===== POST-PROCESSING & GENERATING TECPLOT DATA FILE

 if(dt==zero) then
 if(myid==0) then
    write(*,*) "Overflow."
 end if
 else
    deallocate(qo,qa,qb,de,xim,etm,zem,rr,ss,p,yaco)
 if(nviscous==1) then
    deallocate(txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz)
 end if
 if(tmax>=tsam) then
    nlmx=(3+5*(ndata+1))*(lmx+1)-1; ll=5*(lmx+1)-1; allocate(vart(0:nlmx),vmean(0:ll))
    open(9,file=cdata,access='direct',form='unformatted',recl=nrecs*(nlmx+1),status='old')
    read(9,rec=1) vart(:)
    close(9,status='delete')

!----- CALCULATING UNSTEADY FLUCTUATIONS

 if(ndatafl==1) then
    fctr=half/(times(ndata)-times(0)); vmean(:)=zero
 do n=0,ndata; lis=(3+5*n)*(lmx+1); lie=lis+ll; nn=n/ndata
 if(n*(n-ndata)==0) then
    ra0=fctr*(times(n+1-nn)-times(n-nn))
 else
    ra0=fctr*(times(n+1)-times(n-1))
 end if
    vmean(:)=vmean(:)+ra0*vart(lis:lie)
 end do
 do n=0,ndata; lis=(3+5*n)*(lmx+1); lie=lis+ll
    vart(lis:lie)=vart(lis:lie)-vmean(:)
 do m=1,5; nn=3+5*n+m; l=lis+(m-1)*(lmx+1)
    varr(:)=vart(l:l+lmx); call vminmax(nn)
 end do
 end do
 end if

!----- COLLECTING DATA FROM SUBDOMAINS & BUILDING TECPLOT OUTPUT FILES

    lje=-1
 do n=-1,ndata
    mq=3+2*min(n+1,1); llmb=mq*ltomb-1; allocate(vara(0:llmb),varb(0:llmb))
    ljs=lje+1; lje=ljs+mq*(lmx+1)-1
 if(myid==mo(mb)) then !===========================================================================
    mps=mo(mb); mpe=mps+nbpc(mb,1)*nbpc(mb,2)*nbpc(mb,3)-1
    lis=0; lie=mq*(lmx+1)-1; vara(lis:lie)=vart(ljs:lje)
 do mp=mps+1,mpe
    lis=lie+1; lie=lis+mq*(lxim(mp)+1)*(letm(mp)+1)*(lzem(mp)+1)-1; lmpi = lie-lis+1
    itag=1; call MPI_RECV(vara(lis:lie),lmpi,MPI_REAL4,mp,itag,icom,ista(:,mp),ierr)
 end do
    lis=0
 do mp=mps,mpe; do m=1,mq; do k=0,lzem(mp); do j=0,letm(mp)
    ljs=lpos(mp)+(m-1)*ltomb+k*(leto+1)*(lxio+1)+j*(lxio+1)
    varb(ljs:ljs+lxim(mp))=vara(lis:lis+lxim(mp)); lis=lis+lxim(mp)+1
 end do; end do; end do; end do
    open(9,file=cthead(mb),access='stream',form='unformatted',status='replace')
    call techead(9,n,mb,lh)
    deallocate(vara); allocate(vara(0:lh+llmb)); read(9,pos=1) vara(0:lh-1)
    close(9,status='delete')
    lhmb(mb)=lh+llmb+1; vara(lh:lh+llmb)=varb(:)
 if(mb==0) then !----------------------------------------------------------------------------------
 do mm=1,mbk
    itag=2; call MPI_RECV(lhmb(mm),1,MPI_INTEGER8,mo(mm),itag,icom,ista(:,mp),ierr)
 end do
    llmo=sum(lhmb(:))-1; deallocate(varb); allocate(varb(0:llmo))
    lis=0; lie=lhmb(mb)-1; varb(lis:lie)=vara(:)
 do mm=1,mbk
    lis=lie+1; lie=lis+lhmb(mm)-1; lmpi = lie-lis+1
    itag=3; call MPI_RECV(varb(lis:lie),lmpi,MPI_REAL4,mo(mm),itag,icom,ista(:,mp),ierr)
 end do
    open(0,file=ctecplt(n),status='unknown'); close(0,status='delete') ! 'replace' not suitable as 'recl' may vary
    open(0,file=ctecplt(n),access='direct',form='unformatted',recl=nrecs*(llmo+1),status='new')
    write(0,rec=1) varb(:)
    close(0)
 else !--------------------------------------------------------------------------------------------
    itag=2; call MPI_SEND(lhmb(mb),1,MPI_INTEGER8,mo(0),itag,icom,ierr)
    lmpi = lhmb(mb)
    itag=3; call MPI_SEND(vara(:),lmpi,MPI_REAL4,mo(0),itag,icom,ierr)
 end if !------------------------------------------------------------------------------------------
 else !============================================================================================
    lmpi = lje-ljs+1
    itag=1; call MPI_SEND(vart(ljs:lje),lmpi,MPI_REAL4,mo(mb),itag,icom,ierr)
 end if !==========================================================================================
    deallocate(vara,varb)
 end do

!-----

 end if
 end if

!===== END OF JOB

 if(myid==0) then
    write(*,*) "Finished."
 end if

    call MPI_FINALIZE(ierr)

 end program main3d

!*****
