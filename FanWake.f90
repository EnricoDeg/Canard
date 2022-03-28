!*****
!***** 3D AEROFOIL-TURBULENCE INTERACTION
!*****

 module problemcase

 use mpi
 use subroutineso
 use gridgen
 implicit none

 integer(kind=ni) :: nbody,nthick,ngridv
 integer(kind=ni),dimension(:),allocatable :: idsgnl,lsgnl
 real(kind=nr) :: smg,smgvr,doml0,doml1,domh,span,wlew,wlea,szth0,szth1,szco,skew,spx
 real(kind=nr) :: vtxr,vtxs, radv, k1, k2

 contains

!===== SUBINPUT PARAMETERS

 subroutine inputext

    open(9,file='inputp.dat',status='old')
    read(9,*) cinput,lxi0,lxi1,lxi2
    read(9,*) cinput,let0
    read(9,*) cinput,lze0
    read(9,*) cinput,nbody,nthick
    read(9,*) cinput,ngridv
    read(9,*) cinput,nbpc(:,1)
    read(9,*) cinput,nbpc(:,2)
    read(9,*) cinput,nbpc(:,3)
    read(9,*) cinput,smg,smgvr
    read(9,*) cinput,doml0,doml1,domh
    read(9,*) cinput,span
    read(9,*) cinput,wlew,wlea
    read(9,*) cinput,szth0,szth1,szco
    read(9,*) cinput,skew,spx
    close(9)

    lximb(:)=(/lxi0,lxi1,lxi2,lxi0,lxi1,lxi2/)
    letmb(:)=(/let0,let0,let0,let0,let0,let0/)
    lzemb(:)=(/lze0,lze0,lze0,lze0,lze0,lze0/)

    allocate(idsgnl(0:lze0),lsgnl(0:lze0))

 end subroutine inputext

!===== DOMAIN DECOMPOSITION & BOUNDARY INFORMATION

 subroutine domdcomp

    ip=30*nthick+35*(1-nthick); jp=35*(1-nbody)+nbody*(20+5*nviscous)
 do mm=0,mbk
 select case(mm)
 case(0,3); nbbc(mm,1,:)=(/10,ip/); mbcd(mm,1,:)=(/-1,mm+1/)
 case(1,4); nbbc(mm,1,:)=(/ip,ip/); mbcd(mm,1,:)=(/mm-1,mm+1/)
 case(2,5); nbbc(mm,1,:)=(/ip,10/); mbcd(mm,1,:)=(/mm-1,-1/)
 end select
 select case(mm)
 case(0,2); nbbc(mm,2,:)=(/45,ip/); mbcd(mm,2,:)=(/mm+3,mm+3/)
 case(1); nbbc(mm,2,:)=(/45,ip/); mbcd(mm,2,:)=(/mm+3,mm+3/)
 case(3,5); nbbc(mm,2,:)=(/ip,45/); mbcd(mm,2,:)=(/mm-3,mm-3/)
 case(4); nbbc(mm,2,:)=(/ip,45/); mbcd(mm,2,:)=(/mm-3,mm-3/)
 end select
    nbbc(mm,3,:)=(/45,45/); mbcd(mm,3,:)=(/mm,mm/)
 end do

 end subroutine domdcomp

!===== GRID GENERATION

 subroutine makegrid

    call gridaerofoil(ngridv,nthick,smg,smgvr,doml0,doml1,domh,span,wlew,wlea,szth0,szth1,skew,spx)

 end subroutine makegrid

!===== INITIAL CONDITIONS

 subroutine initialo

!vtxr=0.2_nr*span; 
!vtxs=-0.02_nr; 
!fctr=vtxs/twopi; 
!ra0=-half-span; 
!ra1=-span*tan(11*pi/180); 
!ra2=one/vtxr; ra3=one/vtxr**2
radv = 1.0
k1 = 12.5
k2 = 1.0

do l=0,lmx
!    rv(1)=ra3*((ss(l,1)-ra0)**two+(ss(l,2)-ra1)**two+ss(l,3)**two);
!    rv(2)=ra2*(1.5_nr/(rv(1)+1.0e-32)-one)
!    ao=fctr*exp(half*(one-rv(1)))*rv(1)**0.75_nr;
    ao = k2/2.0/pi * sqrt(exp(1 - k1**2 * (ss(l,1)**2 + ss(l,2)**2)/radv**2))
    bo=(one-half*gamm1*ao*ao)**hamm1
    qa(l,1)=bo;
!    ve(:)=rv(2)*ao*(/-(ss(l,2)-ra1),(ss(l,1)-ra0),zero/);
    ve(:) = (/k1*ss(l,2)*ao/radv, -k1*ss(l,1)*ao/radv, zero/)
    hv2=half*(ve(1)*ve(1)+ve(2)*ve(2)+ve(3)*ve(3))
    qa(l,2:4)=bo*ve(:);
    qa(l,5)=hamhamm1*bo**gam+hv2*bo
 end do

 end subroutine initialo

!===== SETTING UP SPONGE ZONE PARAMETERS

 subroutine spongeup

    ll=-1; ra2=skew/domh; tmpa=pi/szth0; tmpb=pi/szth1
 do l=0,lmx
    ra3=ra2*ss(l,2)
    ra0=tmpa*(ss(l,1)-(ra3-doml0+szth0)); ra1=tmpb*(ra3+doml1-szth1-ss(l,1))
    de(l,1)=szco*half*(two+cos(max(min(ra0,pi),zero))+cos(max(min(ra1,pi),zero)))
    de(l,2)=szco*half*(one+cos(max(min(ra0,pi),zero)))
 if(de(l,1)>sml) then
    ll=ll+1; de(ll,5)=l+sml
 end if
 end do
    lsz=ll
 if(lsz/=-1) then
    allocate(lcsz(0:lsz),asz(0:lsz),bsz(0:lsz))
 do ll=0,lsz; l=de(ll,5); lcsz(ll)=l
    asz(ll)=de(l,1)/yaco(l); bsz(ll)=de(l,2)/yaco(l)
 end do
 end if

 end subroutine spongeup

!===== SPONGE IMPLEMENTATION

 subroutine spongego

 do ll=0,lsz; l=lcsz(ll)
    de(l,1)=de(l,1)+asz(ll)*(qa(l,1)-one)
    de(l,2:4)=de(l,2:4)+bsz(ll)*(qa(l,2:4)-zero)
    de(l,5)=de(l,5)+asz(ll)*(qa(l,5)-hamhamm1)
 end do

 end subroutine spongego

!===== EXTRA CONDITION

 subroutine extracon

 if(nk==nkrk.and.(timo+quarter-tmax)**two<(half*dt)**two) then
    nn=2; ip=0; i=ip*ijk(1,nn)
 if(nbc(nn,ip)==25) then
    open(2,file=cdata,access='direct',form='unformatted',recl=nrecs*(lmx+1),status='old')
    read(2,rec=1) varr(:); ss(:,1)=varr(:)
    close(2)
    open(2,file='misc/wall'//cnnode//'.dat',status='replace')
    write(2,*) 'variables=x,v1,v2,v3'
    write(2,*) 'zone'
    fctr=one/(ijk(2,nn)+1)
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1); rv(:)=zero
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    ra0=two*acos(cm2(jk,1,ip)); ra1=abs(half*sin(ra0)*(tyy(l)-txx(l))+cos(ra0)*txy(l))
    ra2=gam*p(l)/qa(l,1); ra3=sqrt(ra1*qa(l,1))*(ra2+srefoo)/(srefp1dre*ra2**1.5_nr)
    ve(1)=sqrt((etm(l,2)*zem(l,3)-zem(l,2)*etm(l,3))**two+(etm(l,3)*zem(l,1)-zem(l,3)*etm(l,1))**two)
    ve(2)=cm2(jk,1,ip)*(zem(l,2)*xim(l,3)-xim(l,2)*zem(l,3))+cm2(jk,2,ip)*(zem(l,3)*xim(l,1)-xim(l,3)*zem(l,1))
    ve(3)=xim(l,1)*etm(l,2)-etm(l,1)*xim(l,2)
    rv(:)=rv(:)+ra3*abs(ve(:)*yaco(l))
 end do
    write(2,'(2e16.8)') ss(l,1),fctr*rv(:)
 end do
    close(2)
 end if
 end if

 end subroutine extracon

!===== EXTRA BOUNDARY CONDITION LOCATION (NULL)

 subroutine extrabcc(nflag)

 integer(kind=ni),intent(inout) :: nflag

 end subroutine extrabcc

!===== EXTRA BOUNDARY CONDITION IMPLEMENTATION (NULL)

 subroutine extrabcs

 end subroutine extrabcs

!===== SIGNAL RECORDING

 subroutine signalgo

!    open(1,file='signal.dat',access='direct',form='formatted',recl=16,status='old')
!    close(1)

 end subroutine signalgo

!=====

 end module problemcase

!*****
