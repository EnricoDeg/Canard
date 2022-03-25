!*****
!***** 3D FLAT-PLATE GRID GENERATION
!*****

 module gridgen

 use subroutineso
 implicit none

 integer(kind=ni) :: lxi0,lxi1,lxi2,let0,lze0,let1,lete1,lets2,lettr
 integer(kind=ni) :: lxit,lett,lxie0,lxis1,lxie1,lxis2,lete0,lets1,lxisz,im,jm

 integer(kind=ni),dimension(-1:2) :: ilag

 real(kind=nr),dimension(4,2) :: xcij,ycij
 real(kind=nr),dimension(4) :: xco,yco
 real(kind=nr),dimension(-1:2) :: alag,blag

 real(kind=nr),dimension(:,:),allocatable :: xx,yy,zz
 real(kind=nr),dimension(:),allocatable :: xyzmb,zs
 real(kind=nr),dimension(:,:),allocatable :: xp,yp,xq,yq
 real(kind=nr),dimension(:),allocatable :: pxi,qet
 real(kind=nr),dimension(:,:),allocatable :: qeto

 real(kind=nr) :: ts,te,shs,she,shswle
 real(kind=nr) :: xa,xb,xc,xd,xe,xo,ya,yb,yc,yd,yo,sho,pp,qq
 real(kind=nr) :: am,tmp,tmpa,tmpb,gf
 
 character(len=8) :: nfle, mbchar
 contains

!===== GRID GENERATION

 subroutine gridaerofoil(ngridv,nthick,smg,smgvr,doml0,doml1,domh,span,wlew,wlea,szth0,szth1,skew,spx)

 integer(kind=ni),intent(in) :: ngridv,nthick
 real(kind=nr),intent(in) :: smg,smgvr,doml0,doml1,domh,span,wlew,wlea,szth0,szth1,skew,spx

    lxit=lxi0+lxi1+lxi2+2; lett=2*let0+1
    lxie0=lxi0; lxis1=lxie0+1; lxie1=lxis1+lxi1; lxis2=lxie1+1
    lete0=let0; lets1=lete0+1
    lete1=lets1+let1; lets2=lete1+1; lettr=2*let0+let1+2;

    shs=smg; she=shs

    np=3*(lxio+1)*(leto+1)-1
    allocate(xx(0:lxit,0:lett),yy(0:lxit,0:lett),zz(0:lxit,0:lett))
    allocate(xp(0:lxit,0:3),yp(0:lxit,0:3),xq(0:lett,0:3),yq(0:lett,0:3))
    allocate(xyzmb(0:np),zs(0:lze0),pxi(0:lxit),qet(0:lett),qeto(0:lett,2))
 
 
 if(myid==mo(mb)) then
 open(9,file=cgrid,status='unknown'); close(9,status='delete') ! 'replace' not suitable as 'recl' may vary
 open(9,file=cgrid,access='direct',form='unformatted',recl=nrecd*(np+1),status='new')
 
write(mbchar, "(I0)") mb
 open(97, file='grid/v2c_sweep25_'//trim(mbchar)//'.dat')

 do k=0,lze0
 
 zs(k) = span*(real(lze0-k,kind=nr)/lze0-half)
 
! write(nfle,"(I0)") k
! write(mbchar, "(I0)") mb
! open(97, file='grid/oat15a_out_'//trim(mbchar)//'.dat')
 select case(mb);
        case(0); js=0; je=lete0; is=0; ie=lxie0;
        case(1); js=0; je=lete0; is=lxis1; ie=lxie1;
        case(2); js=0; je=lete0; is=lxis2; ie=lxit;
        case(3); js=lets1; je=lett; is=0; ie=lxie0;
        case(4); js=lets1; je=lett; is=lxis1; ie=lxie1;
        case(5); js=lets2; je=lettr; is=lxis2; ie=lxit;
        case(6); js=lets1; je=lete1; is=lxis2; ie=lxit;
 end select
 
 do j=0,je-js
   do i=0,ie-is
     read(97,*) xx(i,j), yy(i,j)
   end do
 end do
 
 l=-3
 do j=0,je-js; do i=0,ie-is; l=l+3
    xyzmb(l:l+2)=(/xx(i,j),yy(i,j),zs(k)/)
 end do; end do
 write(9,rec=k+1) xyzmb(:)
 end do
  
 close(9)
 
 !----- GRID CHECKING
 if(ngridv==1) then
    open(8,file='misc/gridview'//czonet(mb)//'.dat')
    write(8,*) 'variables=v1,v2'; write(8,"('zone i=',i4,' j=',i4)") ie-is+1,je-js+1
 do j=0,je-js; do i=0,ie-is
    write(8,'(2es15.7)') xx(i,j),yy(i,j)
 end do; end do
    close(8)
 else
    open(8,file='misc/gridview'//czonet(mb)//'.dat'); close(8,status='delete')
 end if

 end if
    deallocate(xx,yy,zz,xyzmb,xp,yp,xq,yq,pxi,qet)
  
   end subroutine gridaerofoil
   
   
 end module gridgen

!*****
