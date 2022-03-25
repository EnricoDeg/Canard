!*****
!***** COMPRESSIBLE AERODYNAMICS & AEROACOUSTICS RESEARCH CODE (CANARD)
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
    read(9,*) cinput,reuoo,tempoo
    read(9,*) cinput,amachoo,aoa1,aoa2
    read(9,*) cinput,wtemp
    read(9,*) cinput,cfl
    read(9,*) cinput,tmax,timf,tsam
    read(9,*) cinput,fltk,fltke,fltexr
    read(9,*) cinput,shockr
    read(9,*) cinput,dto
    close(9)

    nnf(:,0)=(/1,2,3/); nnf(:,1)=(/2,3,1/); nnf(:,2)=(/3,1,2/); nnf(:,3)=(/3,2,1/); nnf(:,4)=(/2,1,3/); nnf(:,5)=(/1,3,2/)

    cinput=cinput; fltk=pi*fltk; fltke=pi*fltke; aoa1=pi*aoa1/180.0_nr; aoa2=pi*aoa2/180.0_nr
    amach1=amachoo*cos(aoa1)*cos(aoa2); amach2=amachoo*sin(aoa1)*cos(aoa2); amach3=amachoo*sin(aoa2)
 if(amachoo>sml) then
    reaoo=reuoo/amachoo
 end if
    srefoo=111/tempoo; srefp1dre=(srefoo+one)/reaoo; sqrtrema=sqrt(reaoo*amachoo); sqrtremai=one/max(sqrtrema,sml)
    uoo(:)=(/amach1,amach2,amach3/)
 if(shockr<=zero) then; nshock=0; else; nshock=1; end if

    ll=3+5*(ndata+2)+6; lp=8*mbk+7; lq=12*mbk+11
    allocate(times(0:ndata),cfilet(-1:ndata+2),ctecplt(-1:ndata+2),varm(0:1,0:mpro),varmin(ll),varmax(ll))
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
    cfilet(ndata+1)='mean'
    cfilet(ndata+2)='rstt'
 do n=-1,ndata+2
    ctecplt(n)='data/'//cfilet(n)//'.plt'
 end do
 do mm=0,mbk
    no(2)=mm/100; no(1)=mod(mm,100)/10; no(0)=mod(mm,10)
    cno=achar(no+48); czonet(mm)='z'//cno(2)//cno(1)//cno(0)
    cthead(mm)='data/'//czonet(mm)//'.plt'
 end do
    crestart='data/restart'//czonet(mb)//'.dat'
    cgrid='misc/grid'//czonet(mb)//'.dat'

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
    lh=(ll+1)/ma-1; ip=mod(ll+1,ma)
 if(ip==0) then; ii=0; else; ii=1; end if
 if(lp==0) then
    l=lh; nbc(nn,0)=nbbc(mb,nn,0); nbc(nn,1)=40; mcd(nn,0)=mmcd(nn,0); mcd(nn,1)=myid+mp
 end if
 if(lp>0.and.lp<=ma-ip-1) then
    l=lh; nbc(nn,0)=40; nbc(nn,1)=40; mcd(nn,0)=myid-mp; mcd(nn,1)=myid+mp
 end if
 if(lp>=ma-ip.and.lp<ma-1) then
    l=lh+ii; nbc(nn,0)=40; nbc(nn,1)=40; mcd(nn,0)=myid-mp; mcd(nn,1)=myid+mp
 end if
 if(lp==ma-1) then
    l=lh+ii; nbc(nn,0)=40; nbc(nn,1)=nbbc(mb,nn,1); mcd(nn,0)=myid-mp; mcd(nn,1)=mmcd(nn,1)
 end if
! if(lp==0) then
!    l=ll-((ll+1)/ma)*(ma-1); nbc(nn,0)=nbbc(mb,nn,0); nbc(nn,1)=40; mcd(nn,0)=mmcd(nn,0); mcd(nn,1)=myid+mp
! end if
! if(lp>0.and.lp<ma-1) then
!    l=(ll+1)/ma-1; nbc(nn,0)=40; nbc(nn,1)=40; mcd(nn,0)=myid-mp; mcd(nn,1)=myid+mp
! end if
! if(lp==ma-1) then
!    l=(ll+1)/ma-1; nbc(nn,0)=40; nbc(nn,1)=nbbc(mb,nn,1); mcd(nn,0)=myid-mp; mcd(nn,1)=mmcd(nn,1)
! end if
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

    ltomb=(lxio+1)*(leto+1)*(lzeo+1)-1

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

    allocate(qo(0:lmx,5),qa(0:lmx,5),de(0:lmx,5))
    allocate(xim(0:lmx,3),etm(0:lmx,3),zem(0:lmx,3),rr(0:lmx,3),ss(0:lmx,3))
    allocate(p(0:lmx),yaco(0:lmx),ayaco(0:lmx),amet(0:lmx),varr(0:lmx),nrr(0:lmx),npex(0:lmx))

 if(nviscous==1) then
    allocate(txx(0:lmx),tyy(0:lmx),tzz(0:lmx))
    allocate(txy(0:lmx),tyz(0:lmx),tzx(0:lmx))
    allocate(hxx(0:lmx),hyy(0:lmx),hzz(0:lmx))
 end if

    ii=nbsize(1)-1; jj=nbsize(2)-1; kk=nbsize(3)-1
    allocate(drva1(0:ii,5,0:1),drva2(0:jj,5,0:1),drva3(0:kk,5,0:1))
    allocate(drvb1(0:ii,5,0:1),drvb2(0:jj,5,0:1),drvb3(0:kk,5,0:1))
    allocate(send01(0:ii,0:2,0:1),send02(0:jj,0:2,0:1),send03(0:kk,0:2,0:1))
    allocate(recv01(0:ii,0:2,0:1),recv02(0:jj,0:2,0:1),recv03(0:kk,0:2,0:1))
    allocate(send11(0:ii,0:3,0:1),send12(0:jj,0:3,0:1),send13(0:kk,0:3,0:1))
    allocate(recv11(0:ii,0:3,0:1),recv12(0:jj,0:3,0:1),recv13(0:kk,0:3,0:1))
    allocate(send21(0:ii,0:1),send22(0:jj,0:1),send23(0:kk,0:1))
    allocate(recv21(0:ii,0:1),recv22(0:jj,0:1),recv23(0:kk,0:1))
    allocate(cm1(0:ii,3,0:1),cm2(0:jj,3,0:1),cm3(0:kk,3,0:1))

    allocate(xu(0:lim,3),yu(0:lim,3),xl(0:lim,2),yl(0:lim,2),li(0:lim),sa(0:lim),sb(0:lim),sc(0:lim),sd(0:lim),se(0:lim))

    allocate(vsamp(0:lmx,5),vmean(0:lmx,5),v2mean(0:lmx,6))

!===== EXTRA COEFFICIENTS FOR DOMAIN BOUNDARIES

    albed(:,0,0)=(/zero,zero,one,alpha01,beta02/)
    albed(:,1,0)=(/zero,alpha10,one,alpha12,beta13/)
!    albed(:,2,0)=(/beta2,alpha2,one,alpha2,beta2/)
    albed(:,2,0)=(/beta20,alpha21,one,alpha23,beta24/)
    dbc(:,0)=(/a01,a02,a03,a04,a05,a06/)
    dbc(:,1)=(/a10,a12,a13,a14,a15,a16/)
!    dbc(:,2)=(/-ab2,-aa2,aa2,ab2,zero,zero/)
    dbc(:,2)=(/a20,a21,a23,a24,a25,a26/)
    albed(:,0,1)=(/zero,zero,one,alpha,beta/)
    albed(:,1,1)=(/zero,alpha,one,alpha,beta/)
    albed(:,2,1)=(/beta,alpha,one,alpha,beta/)

    call fcbcm(fltk,fltke,fltexr)
    call fcint(fltk,half,alphf,betf,fa,fb,fc)
    albef(:,0,1)=(/zero,zero,one,alphf,betf/)
    albef(:,1,1)=(/zero,alphf,one,alphf,betf/)
    albef(:,2,1)=(/betf,alphf,one,alphf,betf/)

    pbco(:,:,:)=zero; pbci(:,:,:)=zero
 do nt=0,1; ii=lmd+nt*(lmf-lmd)
    call sbcco(nt)
 do j=0,1
    pbcot(j,nt)=sum(pbco(0:ii,j,nt))
 end do
 end do

!===== PENTADIAGONAL MATRICES FOR DIFFERENCING & FILETERING

 do nn=1,3
 select case(nn)
 case(1); is=0; ie=is+lxi; case(2); is=lxi+1; ie=is+let; case(3); is=lxi+let+2; ie=is+lze
 end select
 do ip=0,1; np=nbc(nn,ip)
 select case(np)
 case(10,20,25,30); ndf(nn,ip,0)=0; ndf(nn,ip,1:2)=0
 case(35,40,45); ndf(nn,ip,0)=1; ndf(nn,ip,1:2)=1
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

    yaco(:)=three/(qo(:,1)*xim(:,1)+qo(:,2)*etm(:,1)+qo(:,3)*zem(:,1)&
                  +qa(:,1)*xim(:,2)+qa(:,2)*etm(:,2)+qa(:,3)*zem(:,2)&
                  +de(:,1)*xim(:,3)+de(:,2)*etm(:,3)+de(:,3)*zem(:,3))

    ayaco(:)=abs(yaco(:))
    amet(:)=xim(:,1)*xim(:,1)+xim(:,2)*xim(:,2)+xim(:,3)*xim(:,3)&
           +etm(:,1)*etm(:,1)+etm(:,2)*etm(:,2)+etm(:,3)*etm(:,3)&
           +zem(:,1)*zem(:,1)+zem(:,2)*zem(:,2)+zem(:,3)*zem(:,3)

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

!    cbca(:,:)=zero; cbca(1,1:2)=(/alpha01,beta02/); cbca(2,1:3)=(/one,alpha12,beta13/); cbca(3,1:4)=(/alpha2,one,alpha2,beta2/)
    cbca(:,:)=zero; cbca(1,1:2)=(/alpha01,beta02/); cbca(2,1:3)=(/one,alpha12,beta13/); cbca(3,1:4)=(/alpha21,one,alpha23,beta24/)
 do i=4,mbci
    cbca(i,i-3:i)=(/beta,alpha,one,alpha/); if(i<mbci) then; cbca(i,i+1)=beta; end if
 end do
!    rbci(:)=zero; rbci(1:3)=(/one,alpha10,beta2/)
    rbci(:)=zero; rbci(1:3)=(/one,alpha10,beta20/)
    call mtrxi(cbca(:,:),cbcs(:,:),1,mbci); sbci(:)=-matmul(cbcs(:,:),rbci(:))
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
    n=0; ndt=0; dt=zero; dts=zero; dte=zero; timo=zero; umf(:)=nsmf*uoo(:); shockc=zero
    call initialo
 else
    open(9,file=crestart,access='direct',form='unformatted',recl=5*nrecd,status='old')
    read(9,rec=1) cha(:); read(9,rec=2) dha(:)
    n=cha(1)+sml; ndt=cha(2)+sml; dt=cha(3); dts=cha(4); dte=cha(5); timo=dha(1); umf(:)=dha(2:4); shockc=dha(5)
    lp=lpos(myid)+2
 do k=0,lze; do j=0,let; lq=lp+lio(j,k)
 do i=0,lxi; l=indx3(i,j,k,1)
    read(9,rec=lq+i+1) qa(l,:)
 end do
 end do; end do
    close(9)
    tsam=timo
 end if
    dtsum=zero; vsamp(:,:)=zero; vmean(:,:)=zero; v2mean(:,:)=zero

!============================================
!===== BEGINNING OF TIME MARCHING IN SOLUTION
!============================================

 if(myid==0) then
    open(1,file='signal.dat',access='direct',form='formatted',recl=16,status='replace'); close(1)
 end if
    call MPI_BARRIER(icom,ierr)

    wts=MPI_WTIME()

    ndati=-1; nsigi=-1
 do while(timo<tmax.and.(dt/=zero.or.n<=2))

!----- FILTERING & RE-INITIALISING

 if(nshock==1) then
    rr(:,3)=gamgamm1*(qa(:,1)*qa(:,5)-half*(qa(:,2)*qa(:,2)+qa(:,3)*qa(:,3)+qa(:,4)*qa(:,4)))
 if(mod(n,10)==0) then
    de(:,1)=qa(:,2)-umf(1)*qa(:,1); de(:,2)=qa(:,3)-umf(2)*qa(:,1); de(:,3)=qa(:,4)-umf(3)*qa(:,1)
    rr(:,2)=(de(:,1)*de(:,1)+de(:,2)*de(:,2)+de(:,3)*de(:,3))/rr(:,3)
    res=maxval(rr(:,2))
    call MPI_ALLREDUCE(res,fctr,1,MPI_REAL8,MPI_MAX,icom,ierr)
    shockc=ham*shockr*max(sqrt(fctr)-onethird,zero)
 end if
    rr(:,1)=ham*rr(:,3)/qa(:,1)
    call mpigo(2,nrone,n45no,99); call shockdet(1); call shockdet(2); call shockdet(3)
 end if

    fctr=sqrt(umf(1)*umf(1)+umf(2)*umf(2)+umf(3)*umf(3))/max(amachoo,sml)
    rv(:)=(/fltk+fctr*(0.5_nr*pi-fltk),fltke+fctr*(0.75_nr*pi-fltke),fltexr/)
    call fcbcm(rv(1),rv(2),rv(3))
    call fcint(rv(1),half,alphf,betf,fa,fb,fc)
    albef(:,0,1)=(/zero,zero,one,alphf,betf/)
    albef(:,1,1)=(/zero,alphf,one,alphf,betf/)
    albef(:,2,1)=(/betf,alphf,one,alphf,betf/)

    nt=1; pbco(:,:,nt)=zero; pbci(:,:,nt)=zero; call sbcco(nt)
 do j=0,1
    pbcot(j,nt)=sum(pbco(0:lmf,j,nt))
 end do
 do nn=1,3
 select case(nn)
 case(1); is=0; ie=is+lxi; case(2); is=lxi+1; ie=is+let; case(3); is=lxi+let+2; ie=is+lze
 end select
    ns=ndf(nn,0,nt); ne=ndf(nn,1,nt)
    call penta(yu(:,:),yl(:,:),is,ie,ns,ne,nt)
 end do

    lp=mod(n,6)
 do m=1,5
    rr(:,1)=qa(:,m)
    call mpigo(1,nrone,n45no,3*(m-1)+1); call filte(nnf(1,lp))
    call mpigo(1,nrone,n45no,3*(m-1)+2); call filte(nnf(2,lp))
    call mpigo(1,nrone,n45no,3*(m-1)+3); call filte(nnf(3,lp))
    qa(:,m)=rr(:,1)
 end do
    qo(:,:)=qa(:,:)

 if(mod(n,nscrn)==0) then
 if(nshock==1) then
    res=sqrt(maxval(de(:,1)*de(:,1)+de(:,2)*de(:,2)+de(:,3)*de(:,3)))
    call MPI_ALLREDUCE(res,fctr,1,MPI_REAL8,MPI_MAX,icom,ierr)
 else
    fctr=zero
 end if
 if(myid==0) then
    write(*,"(' n =',i8,'   time =',es12.4,'   shockd =',es12.4)") n,timo,fctr
 end if
 end if

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
    rr(:,1)=abs(xim(:,1)*(de(:,2)-umf(1))+xim(:,2)*(de(:,3)-umf(2))+xim(:,3)*(de(:,4)-umf(3)))&
           +abs(etm(:,1)*(de(:,2)-umf(1))+etm(:,2)*(de(:,3)-umf(2))+etm(:,3)*(de(:,4)-umf(3)))&
           +abs(zem(:,1)*(de(:,2)-umf(1))+zem(:,2)*(de(:,3)-umf(2))+zem(:,3)*(de(:,4)-umf(3)))
    res=maxval(ayaco(:)*(sqrt(amet(:)*de(:,5))+rr(:,1)))
    call MPI_ALLREDUCE(res,fctr,1,MPI_REAL8,MPI_MAX,icom,ierr)
    ra0=cfl/fctr; ra1=ra0
 if(nviscous==1) then
    res=maxval(amet(:)*ayaco(:)*ayaco(:)*de(:,1)*ss(:,1))
    call MPI_ALLREDUCE(res,fctr,1,MPI_REAL8,MPI_MAX,icom,ierr)
    ra1=half/fctr
 end if
    dte=min(ra0,ra1)
 else
    dte=dto
 end if
 end if
    dt=dts+(dte-dts)*sin(0.05_nr*pi*(n-ndt))**two
 end if

!----- VISCOUS SHEAR STRESSES & HEAT FLUXES

 if(timo>=timf) then
    nvisc=nviscous
 else
    nvisc=nviscous*(nk/nkrk)
 end if
 if(nvisc==1) then
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

    rr(:,1)=de(:,2)-umf(1)
    rr(:,2)=de(:,3)-umf(2)
    rr(:,3)=de(:,4)-umf(3)
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
 if(nvisc==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*txx(:)-xim(:,2)*txy(:)-xim(:,3)*tzx(:)
    rr(:,2)=rr(:,2)-etm(:,1)*txx(:)-etm(:,2)*txy(:)-etm(:,3)*tzx(:)
    rr(:,3)=rr(:,3)-zem(:,1)*txx(:)-zem(:,2)*txy(:)-zem(:,3)*tzx(:)
 end if
    m=2; call mpigo(0,nrall,n45no,m); call deriv(1,1,m); call deriv(2,2,m); call deriv(3,3,m)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    rr(:,1)=qa(:,3)*ss(:,1)+xim(:,2)*p(:)
    rr(:,2)=qa(:,3)*ss(:,2)+etm(:,2)*p(:)
    rr(:,3)=qa(:,3)*ss(:,3)+zem(:,2)*p(:)
 if(nvisc==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*txy(:)-xim(:,2)*tyy(:)-xim(:,3)*tyz(:)
    rr(:,2)=rr(:,2)-etm(:,1)*txy(:)-etm(:,2)*tyy(:)-etm(:,3)*tyz(:)
    rr(:,3)=rr(:,3)-zem(:,1)*txy(:)-zem(:,2)*tyy(:)-zem(:,3)*tyz(:)
 end if
    m=3; call mpigo(0,nrall,n45no,m); call deriv(1,1,m); call deriv(2,2,m); call deriv(3,3,m)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    rr(:,1)=qa(:,4)*ss(:,1)+xim(:,3)*p(:)
    rr(:,2)=qa(:,4)*ss(:,2)+etm(:,3)*p(:)
    rr(:,3)=qa(:,4)*ss(:,3)+zem(:,3)*p(:)
 if(nvisc==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*tzx(:)-xim(:,2)*tyz(:)-xim(:,3)*tzz(:)
    rr(:,2)=rr(:,2)-etm(:,1)*tzx(:)-etm(:,2)*tyz(:)-etm(:,3)*tzz(:)
    rr(:,3)=rr(:,3)-zem(:,1)*tzx(:)-zem(:,2)*tyz(:)-zem(:,3)*tzz(:)
 end if
    m=4; call mpigo(0,nrall,n45no,m); call deriv(1,1,m); call deriv(2,2,m); call deriv(3,3,m)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    de(:,5)=qa(:,5)+p(:)
    rr(:,1)=de(:,5)*ss(:,1)+p(:)*(umf(1)*xim(:,1)+umf(2)*xim(:,2)+umf(3)*xim(:,3))
    rr(:,2)=de(:,5)*ss(:,2)+p(:)*(umf(1)*etm(:,1)+umf(2)*etm(:,2)+umf(3)*etm(:,3))
    rr(:,3)=de(:,5)*ss(:,3)+p(:)*(umf(1)*zem(:,1)+umf(2)*zem(:,2)+umf(3)*zem(:,3))
 if(nvisc==1) then
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
    lmpi=5*nbsize(nn)
    ir=ir+1; call MPI_ISEND(drva(:,:,ip),lmpi,MPI_REAL8,mcd(nn,ip),itag+iq,icom,ireq(ir),ierr)
    ir=ir+1; call MPI_IRECV(drvb(:,:,ip),lmpi,MPI_REAL8,mcd(nn,ip),itag+ip,icom,ireq(ir),ierr)
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
    if(ra0*(vn-vs+ao)>zero) then; cha(4)=zero; end if
    if(ra0*(vn-vs-ao)>zero) then; cha(5)=zero; end if
    call xtr2q(cm(jk,:,ip)); dha(:)=matmul(xt(:,:),(cha(:)-drva(jk,:,ip)))
 do ii=0,mbci; l=indx3(i+iq*ii,j,k,nn)
    ll=ll+1; de(l,:)=de(l,:)+sbcc(ll)*dha(:)
 end do
 end do
 end do
 case(20,25); dtwi=two*ra0/(nkrk*dt+sml)
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    call eleme(l,cm(jk,:,ip)); cha(:)=drva(jk,:,ip); dha(:)=drvb(jk,:,ip)
 select case(npex(l))
 case(0); cha(4+ip)=cha(5-ip)+dtwi*(vn-vs)*aoi*qa(l,1)
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
    if(ra0*(vn-vs)>zero) then; cha(1:3)=dha(1:3); end if
    if(ra0*(vn-vs+ao)>zero) then; cha(4)=dha(4); end if
    if(ra0*(vn-vs-ao)>zero) then; cha(5)=dha(5); end if
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
 if(np==25) then
 do k=0,ijk(3,nn); do j=0,ijk(2,nn); l=indx3(i,j,k,nn)
    fctr=(1-npex(l))*qa(l,1); qa(l,2:4)=npex(l)*qa(l,2:4)+fctr*umf(:)
 end do; end do
 if(wtemp>zero) then; ra0=hamhamm1*wtemp
 do k=0,ijk(3,nn); do j=0,ijk(2,nn); l=indx3(i,j,k,nn)
    qa(l,5)=npex(l)*qa(l,5)+(1-npex(l))*(ra0*qa(l,1)+half*sum(qa(l,2:4)*qa(l,2:4))/qa(l,1))
 end do; end do
 end if
 end if
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
 do j=1,imjl(ii); js=lq; je=js+ke-ks; lq=lq+je-js+1; kk=ke-ks+1; jj=je-js+1
    ir=ir+1; call MPI_ISEND(rr(ks:ke,1),kk,MPI_REAL8,jlcd(i,j),itag,icom,ireq(ir),ierr)
    ir=ir+1; call MPI_IRECV(rr(js:je,2),jj,MPI_REAL8,jlcd(i,j),itag,icom,ireq(ir),ierr)
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
    lmpi=5*nbsize(nn)
    ir=ir+1; call MPI_ISEND(drva(:,:,ip),lmpi,MPI_REAL8,mcd(nn,ip),itag+iq,icom,ireq(ir),ierr)
    ir=ir+1; call MPI_IRECV(drvb(:,:,ip),lmpi,MPI_REAL8,mcd(nn,ip),itag+ip,icom,ireq(ir),ierr)
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

 if(timo-dt>=tsam) then
    fctr=half*gamm1
    de(:,1)=half*(qo(:,1)+qa(:,1)); rr(:,1)=half/de(:,1)    
    de(:,2)=rr(:,1)*(qo(:,2)+qa(:,2)); de(:,3)=rr(:,1)*(qo(:,3)+qa(:,3)); de(:,4)=rr(:,1)*(qo(:,4)+qa(:,4))
    de(:,5)=fctr*(qo(:,5)+qa(:,5)-de(:,1)*(de(:,2)*de(:,2)+de(:,3)*de(:,3)+de(:,4)*de(:,4)))
 do m=1,6; mps=mod(m-1,3)+2; mpe=mod(m,3)+2
 if(m<=5) then
    vsamp(:,m)=vsamp(:,m)+dt*de(:,m)
 end if
 select case(m)
 case(1:3); v2mean(:,m)=v2mean(:,m)+dt*de(:,m+1)*de(:,m+1)
 case(4:6); v2mean(:,m)=v2mean(:,m)+dt*de(:,mps)*de(:,mpe)
 end select
 end do
    dtsum=dtsum+dt; ra0=(tmax-tsam)/ndata; ra1=min(tsam+ra0*(half+(ndati+1)),tmax)
 if((timo-dt-ra1)*(timo-ra1)<zero) then
    ndati=ndati+1
 select case(ndati/ndata+(ndata+ndati)/(ndata+1))
 case(0); times(ndati)=timo-dtsum; case(1); times(ndati)=timo-half*dtsum; case(2); times(ndati)=timo
 end select
    vmean(:,:)=vmean(:,:)+vsamp(:,:)
    open(0,file=cdata,access='direct',form='unformatted',recl=nrecs*(lmx+1),status='old')
 if(n==1) then
    res=zero; fctr=one
 else
    res=ndataav/dtsum; fctr=one-ndataav
 end if
 do m=1,5
    varr(:)=res*vsamp(:,m)+fctr*de(:,m)
    nn=3+5*ndati+m; write(0,rec=nn) varr(:)
 end do
    close(0)
    dtsum=zero; vsamp(:,:)=zero
 if(ndati==ndata) then
    res=one/(times(ndata)-times(0)); vmean(:,:)=res*vmean(:,:); v2mean(:,:)=res*v2mean(:,:)
 end if
 end if
 end if

    call signalgo

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

 if(dt==zero) then
 if(myid==0) then
    write(*,*) "Overflow."
 end if
 else
 if(nrestart==1) then
 if(myid==mo(mb)) then
    allocate(vrss(0:ltomb,5),vrst(0:ltomb,5))
    lis=0; lie=lmx; vrss(lis:lie,:)=qa(:,:)
    mps=mo(mb); mpe=mps+nbpc(mb,1)*nbpc(mb,2)*nbpc(mb,3)-1
 do mp=mps+1,mpe
    lis=lie+1; lie=lis+(lxim(mp)+1)*(letm(mp)+1)*(lzem(mp)+1)-1; lmpi=lie-lis+1
 do m=1,5
    itag=m; call MPI_RECV(vrss(lis:lie,m),lmpi,MPI_REAL8,mp,itag,icom,ista(:,mp),ierr)
 end do
 end do
    lie=-1
 do mp=mps,mpe; do k=0,lzem(mp); do j=0,letm(mp)
    lis=lie+1; lie=lis+lxim(mp); ljs=lpos(mp)+k*(leto+1)*(lxio+1)+j*(lxio+1); lje=ljs+(lie-lis)
    vrst(ljs:lje,:)=vrss(lis:lie,:)
 end do; end do; end do
    cha(:)=(/real(n,kind=nr),real(ndt,kind=nr),dt,dts,dte/); dha(:)=(/timo,umf(1),umf(2),umf(3),shockc/)
    open(9,file=crestart,status='unknown'); close(9,status='delete') ! 'replace' not suitable as 'recl' may vary
    open(9,file=crestart,access='direct',form='unformatted',recl=5*nrecd,status='new')
    write(9,rec=1) cha(:); write(9,rec=2) dha(:)
 do l=0,ltomb
    write(9,rec=l+3) vrst(l,:)
 end do
    close(9)
    deallocate(vrss,vrst)
 else
    lmpi=lmx+1
 do m=1,5
    itag=m; call MPI_SEND(qa(:,m),lmpi,MPI_REAL8,mo(mb),itag,icom,ierr)
 end do
 end if
 end if

!===== POST-PROCESSING & GENERATING TECPLOT DATA FILE

    deallocate(qo,qa,de,xim,etm,zem,rr,ss,p,yaco)
 if(nviscous==1) then
    deallocate(txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz)
 end if
 if(tmax>=tsam) then
    nlmx=(3+5*(ndata+2)+6)*(lmx+1)-1; allocate(vart(0:nlmx))
    open(9,file=cdata,access='stream',form='unformatted',status='old')
    read(9,pos=1) vart(0:nlmx-(5+6)*(lmx+1))
    close(9,status='delete')
 do n=0,ndata+2; lis=(3+5*n)*(lmx+1)
 if(n<=ndata) then; mq=5
 if(ndatafl==1) then
 do m=1,5; l=lis+(m-1)*(lmx+1)
    vart(l:l+lmx)=vart(l:l+lmx)-vmean(:,m)
 end do
 end if
 end if
 if(n==ndata+1) then; mq=5
 do m=1,5; l=lis+(m-1)*(lmx+1)
    vart(l:l+lmx)=vmean(:,m)
 end do
 end if
 if(n==ndata+2) then; mq=6
 do m=1,6; l=lis+(m-1)*(lmx+1); mps=mod(m-1,3)+2; mpe=mod(m,3)+2
 select case(m)
 case(1:3); vart(l:l+lmx)=v2mean(:,m)!-vmean(:,m+1)*vmean(:,m+1)
 case(4:6); vart(l:l+lmx)=v2mean(:,m)!-vmean(:,mps)*vmean(:,mpe)
 end select
 end do
 end if
 do m=1,mq; l=lis+(m-1)*(lmx+1); nn=3+5*n+m
    varr(:)=vart(l:l+lmx); call vminmax(nn)
 end do
 end do

!----- COLLECTING DATA FROM SUBDOMAINS & BUILDING TECPLOT OUTPUT FILES

    lje=-1
 do n=-1,ndata+2
 select case((n+1)/(ndata+3)+((ndata+3)+(n+1))/(ndata+4))
 case(0); mq=3; case(1); mq=5; case(2); mq=6
 end select
    llmb=mq*(ltomb+1)-1; allocate(vara(0:llmb),varb(0:llmb))
    ljs=lje+1; lje=ljs+mq*(lmx+1)-1
 if(myid==mo(mb)) then !===========================================================================
    mps=mo(mb); mpe=mps+nbpc(mb,1)*nbpc(mb,2)*nbpc(mb,3)-1
    lis=0; lie=mq*(lmx+1)-1; vara(lis:lie)=vart(ljs:lje)
 do mp=mps+1,mpe
    lis=lie+1; lie=lis+mq*(lxim(mp)+1)*(letm(mp)+1)*(lzem(mp)+1)-1; lmpi=lie-lis+1
    itag=1; call MPI_RECV(vara(lis:lie),lmpi,MPI_REAL4,mp,itag,icom,ista(:,mp),ierr)
 end do
    lie=-1
 do mp=mps,mpe; do m=1,mq; do k=0,lzem(mp); do j=0,letm(mp)
    lis=lie+1; lie=lis+lxim(mp); llis=lpos(mp)+(m-1)*(ltomb+1)+k*(leto+1)*(lxio+1)+j*(lxio+1); llie=llis+lie-lis
    varb(llis:llie)=vara(lis:lie)
 end do; end do; end do; end do
    open(9,file=cthead(mb),access='stream',form='unformatted',status='replace')
    call techead(9,n,mb,lh)
    deallocate(vara); allocate(vara(0:lh+llmb)); read(9,pos=1) vara(0:lh-1)
    close(9,status='delete')
    lhmb(mb)=lh+llmb+1; vara(lh:lh+llmb)=varb(:)
 if(mb==0) then !----------------------------------------------------------------------------------
 do mm=1,mbk
    itag=2; call MPI_RECV(lhmb(mm),1,MPI_INTEGER8,mo(mm),itag,icom,ista(:,mo(mm)),ierr)
 end do
    llmo=sum(lhmb(:))-1; deallocate(varb); allocate(varb(0:llmo))
    lis=0; lie=lhmb(mb)-1; varb(lis:lie)=vara(:)
 do mm=1,mbk
    lis=lie+1; lie=lis+lhmb(mm)-1; lmpi=lie-lis+1
    itag=3; call MPI_RECV(varb(lis:lie),lmpi,MPI_REAL4,mo(mm),itag,icom,ista(:,mo(mm)),ierr)
 end do
    open(0,file=ctecplt(n),access='stream',form='unformatted',status='replace')
    write(0,pos=1) varb(:)
    close(0)
 else !--------------------------------------------------------------------------------------------
    itag=2; call MPI_SEND(lhmb(mb),1,MPI_INTEGER8,mo(0),itag,icom,ierr)
    lmpi=lhmb(mb)
    itag=3; call MPI_SEND(vara(:),lmpi,MPI_REAL4,mo(0),itag,icom,ierr)
 end if !------------------------------------------------------------------------------------------
 else !============================================================================================
    lmpi=lje-ljs+1
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
