!*****
!***** 3D AEROFOIL-TURBULENCE INTERACTION
!*****
 module problemcase
 use mpi
 use subroutineso
 use gridgen
 implicit none

 integer(kind=ni) :: nbody,nthick,ngridv
 real(kind=nr) :: smg,smgvr,doml0,doml1,domh,span,wlew,wlea,szth0,szth1,szco,skew,spx
 real(kind=nr) :: vtxr,vtxs, radv, k1, k2

 contains

!===== GRID GENERATION

 subroutine makegrid

    call gridaerofoil(ngridv,nthick,smg,smgvr,doml0,doml1,domh,span,wlew,wlea,szth0,szth1,skew,spx)

 end subroutine makegrid

!===== INITIAL CONDITIONS

 subroutine initialo

radv = 1.0
k1 = 12.5
k2 = 1.0

do l=0,lmx
    ao = k2/2.0/pi * sqrt(exp(1 - k1**2 * (ss(l,1)**2 + ss(l,2)**2)/radv**2))
    bo=(one-half*gamm1*ao*ao)**hamm1
    qa(l,1)=bo;
    ve(:) = (/k1*ss(l,2)*ao/radv, -k1*ss(l,1)*ao/radv, zero/)
    hv2=half*(ve(1)*ve(1)+ve(2)*ve(2)+ve(3)*ve(3))
    qa(l,2:4)=bo*ve(:);
    qa(l,5)=hamhamm1*bo**gam+hv2*bo
 end do

 end subroutine initialo

!
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

! subroutine extrabcc(nflag)

! integer(kind=ni),intent(inout) :: nflag

! end subroutine extrabcc

!===== EXTRA BOUNDARY CONDITION IMPLEMENTATION (NULL)

! subroutine extrabcs
!
! end subroutine extrabcs

!===== SIGNAL RECORDING

 subroutine signalgo

!    open(1,file='signal.dat',access='direct',form='formatted',recl=16,status='old')
!    close(1)

 end subroutine signalgo

!=====

 end module problemcase

!*****
