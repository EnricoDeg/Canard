!*****
!***** BASIC SUBROUTINES
!*****

 module subroutineso

! use mainvar2d
 use mainvar3d
 implicit none

 contains

!===== SUBROUTINE FOR CHOLESKY DECOMPOSITION OF PENTADIAGONAL MATRICES

 subroutine penta(xu,xl,is,ie,ns,ne,nt)

 integer(kind=ni),intent(in) :: is,ie,ns,ne,nt
 real(kind=nr),dimension(0:lim,3),intent(inout) :: xu
 real(kind=nr),dimension(0:lim,2),intent(inout) :: xl
 real(kind=nr),dimension(-2:2,0:2,0:1) :: albe
 real(kind=nr) :: alpho,beto

 if(nt==0) then
    albe=albed; alpho=alpha; beto=beta
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

 subroutine fcbcm(fltk,fltke,fltexr)
 
 real(kind=nr),intent(in) :: fltk,fltke,fltexr
 real(kind=nr),dimension(0:2) :: fltkbc
 real(kind=nr) :: alphz,betz,za,zb,zc

    ao=log(fltexr); fltkbc(:)=(fltke-fltk)*(one+cos(pi*(/3,4,5/)/6))+fltk

    call fcint(fltkbc(0),half,alphz,betz,za,zb,zc); fctr=one/(one+alphz*fltexr+betz*fltexr**two)
    albef(:,0,0)=(/zero,zero,one,alphz*fctr,betz*fctr/); res=(fltexr-1)*(za+zc+(fltexr+1)*(zb+fltexr*zc))/ao
    fbc(:,0)=(/za-5*res/3,zb+10*res/21,zc-5*res/42,5*res/252,-res/630,zero/)*fctr

    call fcint(fltkbc(1),half,alphz,betz,za,zb,zc)
    albef(:,1,0)=(/zero,alphz+betz*fltexr,one,alphz,betz/); res=(fltexr-1)*(zb+zc*(fltexr+1))/ao
    fbc(:,1)=(/za+zb+zc+1627*res/1260,za+10*res/21,zb-5*res/42,zc+5*res/252,-res/630,zero/)

    call fcint(fltkbc(2),half,alphz,betz,za,zb,zc)
    albef(:,2,0)=(/betz,alphz,one,alphz,betz/); res=zc*(fltexr-1)/ao
    fbc(:,2)=(/zb+zc+1627*res/1260,za-5*res/3,za-5*res/42,zb+5*res/252,zc-res/630,zero/)

!    fltkbc(:)=(fltke-fltk)*(one+cos(pi*(/3,4,5/)/6))+fltk
!    
!    call fcint(fltkbc(0),half,alphz,betz,za,zb,zc); fctr=one+(alphz+betz)*(one-fltexr); res=one+fltexr
!    albef(:,0,0)=(/zero,zero,fctr,alphz*res,betz*res/)/fctr
!    fbc(:,0)=(/za*res,zb*res,zc*res,zero,zero/)/fctr
!
!    call fcint(fltkbc(1),half,alphz,betz,za,zb,zc); fctr=one+betz*fltexr
!    albef(:,1,0)=(/zero,alphz+betz*(one-fltexr),fctr,alphz,betz/)/fctr
!    fbc(:,1)=(/za+zb+zc-(zb+zc)*fltexr,za+zc*fltexr,zb,zc,zero/)/fctr
!
!    call fcint(fltkbc(2),half,alphz,betz,za,zb,zc)
!    albef(:,2,0)=(/betz,alphz,one,alphz,betz/)
!    fbc(:,2)=(/zb+zc*(one-fltexr),za+zc*fltexr,za,zb,zc/)

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

 subroutine sbcco(nt)

 integer(kind=ni),intent(in) :: nt
 real(kind=nr),dimension(:,:),allocatable :: ax,bx,rx,sx
 real(kind=nr),dimension(-2:2,0:2,0:1) :: albe
 real(kind=nr),dimension(0:5,0:2) :: abc

    lp=2*nt-1
 select case(nt)
 case(0); ll=lmd; albe=albed; abc=dbc
 case(1); ll=lmf; albe=albef; abc=fbc
 end select
    is=1; ie=2*(ll+1)
    allocate(ax(is:ie,is:ie),bx(is:ie,is:ie),rx(is:ie,is:ie),sx(is:ie,is:ie)); ax(:,:)=zero; bx(:,:)=zero

    ax(is,is:is+2)=albe(0:2,0,0); bx(is,is+(/1,2,3,4,5,6/))=abc(0:5,0); bx(is,is)=-sum(abc(0:5,0))
    ax(is+1,is:is+3)=albe(-1:2,1,0); bx(is+1,is+(/0,2,3,4,5,6/))=abc(0:5,1); bx(is+1,is+1)=-sum(abc(0:5,1))
    ax(is+2,is:is+4)=albe(-2:2,2,0); bx(is+2,is+(/0,1,3,4,5,6/))=abc(0:5,2); bx(is+2,is+2)=-sum(abc(0:5,2))
 do i=is+3,ie-3
    ax(i,i-2:i+2)=(1-nt)*(/beta,alpha,one,alpha,beta/)+nt*(/betf,alphf,one,alphf,betf/)
    bx(i,i-3:i+3)=(1-nt)*(/-ac,-ab,-aa,zero,aa,ab,ac/)+nt*(/fc,fb,fa,-two*(fa+fb+fc),fa,fb,fc/)
 end do
    ax(ie-2,ie:is:-1)=ax(is+2,is:ie); bx(ie-2,ie:is:-1)=lp*bx(is+2,is:ie)
    ax(ie-1,ie:is:-1)=ax(is+1,is:ie); bx(ie-1,ie:is:-1)=lp*bx(is+1,is:ie)
    ax(ie,ie:is:-1)=ax(is,is:ie); bx(ie,ie:is:-1)=lp*bx(is,is:ie)

    call mtrxi(ax(:,:),sx(:,:),is,ie)

    rx(:,:)=ax(:,:)
    i=ie/2-1; rx(i,i+2)=zero; rx(i+1,i+2)=zero; rx(i+1,i+3)=zero
    i=ie/2+2; rx(i,i-2)=zero; rx(i-1,i-2)=zero; rx(i-1,i-3)=zero
    ax(:,:)=matmul(rx(:,:),matmul(sx(:,:),bx(:,:)))
    i=ie/2+1; pbco(ll:0:-1,0,nt)=ax(i,is:is+ll); pbci(0:ll,0,nt)=ax(i,is+ll+1:ie)
    i=ie/2+2; pbco(ll:0:-1,1,nt)=ax(i,is:is+ll); pbci(0:ll,1,nt)=ax(i,is+ll+1:ie)

    deallocate(ax,bx,rx,sx)

 end subroutine sbcco

!===== SUBROUTINE FOR MATRIX INVERSION

 subroutine mtrxi(ax,sx,is,ie)

 integer(kind=ni),intent(in) :: is,ie
 real(kind=nr),dimension(is:ie,is:ie),intent(in) :: ax
 real(kind=nr),dimension(is:ie,is:ie),intent(inout) :: sx

 integer(kind=ni),dimension(is:ie) :: ipvt
 integer(kind=ni),dimension(1) :: imax
 integer(kind=ni) :: i,j,m
 real(kind=nr),dimension(is:ie,is:ie) :: rx
 real(kind=nr),dimension(is:ie) :: temp

    rx(:,:)=ax(:,:); ipvt(:)=(/(i,i=is,ie)/)
 do i=is,ie
    imax(:)=maxloc(abs(rx(i:ie,i))); m=i-1+imax(1)
 if(m/=i) then
    ipvt((/m,i/))=ipvt((/i,m/)); rx((/m,i/),:)=rx((/i,m/),:)
 end if
    ra0=one/rx(i,i); temp(:)=rx(:,i)
 do j=is,ie
    ra1=ra0*rx(i,j); rx(:,j)=rx(:,j)-ra1*temp(:); rx(i,j)=ra1
 end do
    rx(:,i)=-ra0*temp(:); rx(i,i)=ra0
 end do
    sx(:,ipvt(:))=rx(:,:)

 end subroutine mtrxi

!===== SUBROUTINE FOR MOVING FRAME VELOCITIES

 subroutine movef(dtko,dtk)

 real(kind=nr),intent(in) :: dtko,dtk

 if(nsmf==1) then
    ra0=one/timf; ra1=ra0*min(timo,timf); ra2=ra0*min(timo+dtko,timf)

    fctr=(one-ra1)**three
    dfdt=three*ra0*(one-ra2)**two
    progmf=fctr+dtk*dfdt
    umf(:)=progmf*uoo(:)
 else
    umf(:)=zero
 end if

 end subroutine movef

!===== SUBROUTINE FOR BLASIUS LAMINAR BOUNDARY LAYER

 subroutine lambl(x,y,blu,blv,blm,lbl)

 integer(kind=ni),intent(in) :: lbl
 real(kind=nr),dimension(0:lbl),intent(in) :: x,y
 real(kind=nr),dimension(0:lbl),intent(inout) :: blu,blv,blm

 real(kind=nr) :: blas,spr,eta,etb

    blas=0.3320573362151963_nr; spr=sqrt(prndtl)
    ra0=half*pi; ra1=3*blas*blas/560; ra2=11*blas/420
 do i=0,lbl
    fctr=one/sqrt(x(i))
    eta=fctr*sqrtrema*y(i); etb=spr*eta
    res=0.000001_nr*eta**four*exp(eta**1.3625_nr)
    blu(i)=(blas*eta+ra1*eta**four+res)/(one+ra2*eta**three+res)
    blv(i)=0.8604_nr*fctr*sqrtremai*(sin(ra0*tanh(exp(0.177_nr*eta)-one)))**1.96_nr
    res=0.000001_nr*etb**four*exp(etb**1.3625_nr)
    blm(i)=(blas*etb+ra1*etb**four+res)/(one+ra2*etb**three+res)
 end do

 end subroutine lambl

!===== SUBROUTINE FOR GRID LINE GENERATION

 subroutine gridf(x,xxi,xo,xn,dxo,dxn,lxi,mxin,ip)

 integer(kind=ni),intent(in) :: lxi,mxin,ip
 real(kind=nr),dimension(0:lxi),intent(inout) :: x,xxi
 real(kind=nr),intent(in) :: xo,xn,dxo,dxn

 integer(kind=ni) :: i,ii
 real(kind=nr) :: dxoo,dxnn,aa,bb,cc,dd,xi,fctr

    dxoo=dxo; dxnn=dxn; xi=mxin
 if(dxo==free) then
    dxoo=two*(xn-xo)/xi-dxnn
 end if
 if(dxn==free) then
    dxnn=two*(xn-xo)/xi-dxoo
 end if
    aa=6.0_nr*(xn-xo)-3.0_nr*xi*(dxoo+dxnn)
    bb=15.0_nr*(xo-xn)+xi*(8.0_nr*dxoo+7.0_nr*dxnn)
    cc=10.0_nr*(xn-xo)-xi*(6.0_nr*dxoo+4.0_nr*dxnn)
    dd=xi*dxoo; fctr=one/xi
 do i=0,mxin
    ii=i+ip; xi=i*fctr
    x(ii)=aa*xi**five+bb*xi**four+cc*xi**three+dd*xi+xo
    xxi(ii)=fctr*(five*aa*xi**four+four*bb*xi**three+three*cc*xi**two+dd)
 end do

 end subroutine gridf

!===== SUBROUTINE FOR CHARACTER STRING CONVERSION

 subroutine strio(nfile,lh,cinput)

 integer(kind=ni),intent(in) :: nfile
 integer(kind=ni),intent(inout) :: lh
 character(16),intent(in) :: cinput
 integer(kind=ni) :: ll

 do ll=1,len_trim(cinput)
    write(nfile,pos=4*lh+1) ichar(cinput(ll:ll)); lh=lh+1
 end do
    write(nfile,pos=4*lh+1) 0; lh=lh+1

 end subroutine strio

!===== SUBROUTINE FOR TRIMMING DOWN TO EFFECTIVE ARRAY

 subroutine trimm(itrim,lene,lens)

 integer(kind=ni),dimension(0:ljpl),intent(inout) :: itrim
 integer(kind=ni),intent(inout) :: lene,lens
 integer(kind=ni) :: k,kk

 do k=0,lens-1; do kk=k+1,lens
 if(itrim(kk)-itrim(k)==0) then; itrim(kk)=-1; end if
 end do; end do
    lene=lens
 do k=0,lens
 if(itrim(k)==-1) then; lene=lene-1; itrim(k)=itrim(k+1); itrim(k+1)=-1; end if
 end do

 end subroutine trimm

!===== FUNCTION FOR PROCESSOR-CORE INDEX TRANSFORMATION WITHIN THE BLOCK IN 3D

 function idsd3(i,j,k,mm,nn) result(lm)

 integer(kind=ni),intent(in) :: i,j,k,mm,nn
 integer(kind=ni) :: lm

 select case(nn)
 case(1); lm=mo(mm)+(k*nbpc(mm,2)+j)*nbpc(mm,1)+i
 case(2); lm=mo(mm)+(j*nbpc(mm,2)+i)*nbpc(mm,1)+k
 case(3); lm=mo(mm)+(i*nbpc(mm,2)+k)*nbpc(mm,1)+j
 end select

 end function idsd3

!===== FUNCTION FOR MAIN INDEX TRANSFORMATION IN 3D

 function indx3(i,j,k,nn) result(lm)

 integer(kind=ni),intent(in) :: i,j,k,nn
 integer(kind=ni) :: lm

 select case(nn)
 case(1); lm=(k*(let+1)+j)*(lxi+1)+i
 case(2); lm=(j*(let+1)+i)*(lxi+1)+k
 case(3); lm=(i*(let+1)+k)*(lxi+1)+j
 end select

 end function indx3

!===== FUNCTION FOR MAIN INDEX TRANSFORMATION IN 2D

 function indx2(i,j,nn) result(lm)

 integer(kind=ni),intent(in) :: i,j,nn
 integer(kind=ni) :: lm

 select case(nn)
 case(1); lm=j*(lxi+1)+i
 case(2); lm=i*(lxi+1)+j
 end select

 end function indx2

!=====

 end module subroutineso

!*****