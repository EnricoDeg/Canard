!*****
!***** GRID GENERATION MODULE
!*****

module mo_gridgen
   use mo_kind,       ONLY : ni, nr, ieee64
   use mo_parameters, ONLY : zero, twopi, two, three, onethird, halfpi,     &
                           & free, half, one, four, five, sml, pi, hamhamm1
   use mo_io,         ONLY : cgrid
   use mo_mpi,        ONLY : p_get_process_ID
   implicit none
   public

   type, public :: t_grid_geom
      real(kind=nr)    :: szth0
      real(kind=nr)    :: szth1
      real(kind=nr)    :: doml0
      real(kind=nr)    :: doml1
      real(kind=nr)    :: domh
   end type t_grid_geom 

   integer(kind=ni), private, parameter :: lnaca=90
   integer(kind=ni) :: lxi0,lxi1,lxi2,let0,lze0 ! input vars
   integer(kind=ni) :: gridtype ! 0:toy, 1:aerofoil

   real(kind=nr), private, dimension(0:lnaca,2) :: xnaca,ynaca,xf,yf
   real(kind=nr), private, dimension(4,2) :: xcij,ycij
   real(kind=nr), private, dimension(4) :: xco,yco
   real(kind=nr), private, dimension(-1:2) :: alag,blag

   real(kind=nr), private, dimension(:,:), allocatable :: xx,yy,zz
   real(kind=nr), private, dimension(:),   allocatable :: xyzmb,zs
   real(kind=nr), private, dimension(:,:), allocatable :: xp,yp,xq,yq
   real(kind=nr), private, dimension(:),   allocatable :: pxi,qet
   real(kind=nr), private, dimension(:,:), allocatable :: qeto

   real(kind=nr), private :: ts,te,shs,she,shswle
   real(kind=nr), private :: xa,xb,xc,xd,xe,xo,ya,yb,yc,yd,yo,sho,pp,qq
   real(kind=nr), private :: am,tmp,tmpa,tmpb,gf

   integer(kind=ni), private :: nthick ! input vars
   real(kind=nr),    private :: smg,smgvr,doml0,doml1,domh,span,szth0,szth1,spx ! input vars

   real(kind=nr), private, dimension(5,5) :: xt


   CONTAINS

   SUBROUTINE read_input_gridgen
      
      character(16) :: cinput
      integer(kind=ni) :: fu, rc

      namelist /nml_gridgen/ gridtype,lxi0,let0,lze0,let0,lze0,doml0, &
                             doml1,domh,span,szth0,szth1,nthick, &
                             smg,smgvr,lxi1,lxi2

      open (action='read', file='input.canard', iostat=rc, newunit=fu)
      read (nml=nml_gridgen, iostat=rc, unit=fu)
      close(fu)

   END SUBROUTINE read_input_gridgen

   SUBROUTINE get_grid_geometry(p_grid_geom)
      type(t_grid_geom), intent(out) :: p_grid_geom
   
      p_grid_geom%szth0  = szth0
      p_grid_geom%szth1  = szth1
      p_grid_geom%doml0  = doml0
      p_grid_geom%doml1  = doml1
      p_grid_geom%domh   = domh

   END SUBROUTINE get_grid_geometry

!===== GRID GENERATION

   subroutine makegrid(mb, lxio, leto, mo, nblocks)
      integer(kind=ni), intent(in) :: mb, lxio, leto, nblocks
      integer(kind=ni), intent(in), dimension(0:nblocks) :: mo
      
      if (gridtype == 0) then
         call gridtoy(mb, lxio, leto, mo, nblocks, doml0, doml1, domh, span)
      else if (gridtype == 1) then
         call gridtoy_multi(mb, lxio, leto, mo, nblocks, doml0, doml1, domh, span)
      else if (gridtype == 2) then
         call gridaerofoil(mb, lxio, leto, mo, nblocks, &
                           nthick, smg, smgvr, doml0, doml1, domh, span, szth0, szth1, spx)
      else
         write(*,*) 'Unknown grid type'
         STOP 1
      end if

   end subroutine makegrid

!===== GRID GENERATION

   subroutine gridaerofoil(mb, lxio, leto, mo, nblocks, &
                           nthick, smg, smgvr, doml0, doml1, domh, span, szth0, szth1, spx)
      integer(kind=ni), intent(in) :: mb, lxio, leto, nblocks
      integer(kind=ni), intent(in), dimension(0:nblocks) :: mo
      integer(kind=ni), intent(in) :: nthick
      real(kind=nr),    intent(in) :: smg,smgvr,doml0,doml1,domh,span,szth0,szth1,spx
      integer(kind=ni)             :: m, n, np, ll, lh, k, js, je, jp, jj
      integer(kind=ni)             :: l, j, i, is, ie, ip, ii, im, jm, lxisz
      integer(kind=ni)             :: lxit, lett, lxie0, lxis1, lxie1, lxis2, lete0, lets1
      real(kind=nr)                :: res, ra0, ra1, ra2, ra3, fctr
      real(kind=nr) ,dimension(5)  :: cha, dha
      integer(kind=ni) :: nrecd, myid

      inquire(iolength=ll) real(1.0,kind=ieee64); nrecd=ll
      myid = p_get_process_ID()

      lxit=lxi0+lxi1+lxi2+2
      lett=2*let0+1
      lxie0=lxi0
      lxis1=lxie0+1
      lxie1=lxis1+lxi1
      lxis2=lxie1+1
      lete0=let0
      lets1=lete0+1

      shs=smg
      she=shs

      np=3*(lxio+1)*(leto+1)-1
      allocate(xx(0:lxit,0:lett),yy(0:lxit,0:lett),zz(0:lxit,0:lett))
      allocate(xp(0:lxit,0:3),yp(0:lxit,0:3),xq(0:lett,0:3),yq(0:lett,0:3))
      allocate(xyzmb(0:np),zs(0:lze0),pxi(0:lxit),qet(0:lett),qeto(0:lett,2))

!----- AEROFOIL SECTION

      if(myid==mo(mb)) then
         open(9,file=cgrid,status='unknown')
         close(9,status='delete')
         open(9,file=cgrid,access='direct',form='unformatted',recl=nrecd*(np+1),status='new')
         if(nthick==1) then
            open(8,file='aerofoil.dat')
            do n=1,2
               do i=0,lnaca
                  read(8,*) xnaca(i,n),ynaca(i,n) 
               end do
            end do
            close(8)
            is=0
            ie=lnaca
            do m=0,3
               do i=is,ie,ie-is
                  res=half*(xnaca(i,1)+xnaca(i,2))
                  xnaca(i,1)=res
                  xnaca(i,2)=res
                  res=half*(ynaca(i,1)+ynaca(i,2))
                  ynaca(i,1)=res
                  ynaca(i,2)=res
               end do
               cha(:)=(/-one,four,10.0_nr,four,-one/)/16.0_nr
               xf(:,:)=xnaca(:,:)
               yf(:,:)=ynaca(:,:)
               do n=1,2
                  do i=is+2,ie-2
                     xf(i,n)=sum(cha(:)*xnaca(i-2:i+2,n))
                     yf(i,n)=sum(cha(:)*ynaca(i-2:i+2,n))
                  end do
                  xf(is,n)=sum(cha(:)*(/xnaca(is+2:is+1:-1,3-n),xnaca(is:is+2,n)/))
                  yf(is,n)=sum(cha(:)*(/ynaca(is+2:is+1:-1,3-n),ynaca(is:is+2,n)/))
                  xf(is+1,n)=sum(cha(:)*(/xnaca(is+1,3-n),xnaca(is:is+3,n)/))
                  yf(is+1,n)=sum(cha(:)*(/ynaca(is+1,3-n),ynaca(is:is+3,n)/))
                  xf(ie-1,n)=sum(cha(:)*(/xnaca(ie-3:ie,n),xnaca(ie-1,3-n)/))
                  yf(ie-1,n)=sum(cha(:)*(/ynaca(ie-3:ie,n),ynaca(ie-1,3-n)/))
                  xf(ie,n)=sum(cha(:)*(/xnaca(ie-2:ie,n),xnaca(ie-1:ie-2:-1,3-n)/))
                  yf(ie,n)=sum(cha(:)*(/ynaca(ie-2:ie,n),ynaca(ie-1:ie-2:-1,3-n)/))
               end do
               do n=1,2
                  do i=is,ie
                     xnaca(i,n)=xf(i,n)
                     ynaca(i,n)=yf(i,n)
                  end do
               end do
            end do
            xb=half*(xnaca(is,1)+xnaca(is,2))
            xc=half*(xnaca(ie,1)+xnaca(ie,2))
            yb=half*(ynaca(is,1)+ynaca(is,2))
            yc=half*(ynaca(ie,1)+ynaca(ie,2))
         end if

         do n=1,2
            xp(lxis1,n)=xb
            yp(lxis1,n)=yb
            xp(lxie1,n)=xc
            yp(lxie1,n)=yc
            if(nthick==0) then
               ll=0
               tmpa=shs
               tmpb=she
            else
               ll=12 ! "ll" must be equal to or larger than 4.
               do i=lxis1+1,lxis1+ll
                  fctr=1.0_nr**(i-lxis1-1)
                  if(i<=lxis1+2.or.abs(yp(i-1,n))<abs(ynaca(is+2,n))) then
                     call aerofint(i,i-1,n,2,fctr*shs)
                  else
                     call aerofint(i,i-1,n,1,fctr*shs)
                  end if
               end do
               do i=lxie1-1,lxie1-ll,-1
                  fctr=1.0_nr**(lxie1-1-i)
                  if(i>=lxie1-2.or.abs(yp(i+1,n))<abs(ynaca(ie-2,n))) then
                     call aerofint(i,i+1,n,2,fctr*she)
                  else
                     call aerofint(i,i+1,n,1,fctr*she)
                  end if
               end do
               tmpa=sum(xp(lxis1+ll-4:lxis1+ll,n)*(/3.0_nr,-16.0_nr,36.0_nr,-48.0_nr,25.0_nr/))/12.0_nr
               tmpb=sum(xp(lxie1-ll+4:lxie1-ll:-1,n)*(/3.0_nr,-16.0_nr,36.0_nr,-48.0_nr,25.0_nr/))/12.0_nr
            end if
            ip=lxis1+ll
            im=lxi1-2*ll
            call gridf(xp(:,n),pxi,xp(lxis1+ll,n),xp(lxie1-ll,n),tmpa,-tmpb,lxit,im,ip)
            do i=lxis1+ll+1,lxie1-ll-1
               yp(i,n)=xylagran(i,n,1)
            end do
         end do

!----- SPANWISE EXTENSION

         do k=0,lze0
            zs(k)=span*(real(lze0-k,kind=nr)/lze0-half)

            xa=-doml0
            xd=doml1-szth1
            xe=doml1
            ya=-domh
            yd=domh

            shswle=shs

!----- WAVY LEADING EDGE

            tmp=one
            if(nthick==0) then
               xb=half-tmp
               yb=zero
               xc=half
               yc=zero
            else
               do n=1,2
                  do i=lxis1,lxie1
                     xp(i,n)=tmp*(xp(i,n)-xc)+xc
                     yp(i,n)=tmp*(yp(i,n)-yc)+yc
                  end do
               end do
               xb=half*(xp(lxis1,1)+xp(lxis1,2))
               yb=half*(yp(lxis1,1)+yp(lxis1,2))
            end if

!----- HORIZONTAL INTERFACE POINTS IN XI-DIRECTION

            n=1
            sho=szth0/16
            ip=0
            im=lxi0
            call gridf(xp(:,n),pxi,xa,xb,sho,shswle,lxit,im,ip)
            if(k==0) then
               lh=minloc(abs(xa+szth0-xp(0:lxi0,n)),1)-1
               lxisz=lh
            end if
            ip=lxis2
            im=lxi2-lxisz
            call gridf(xp(:,n),pxi,xc,xd,she,free,lxit,im,ip)
            ip=ip+im
            im=lxisz
            call gridf(xp(:,n),pxi,xd,xe,pxi(ip),free,lxit,im,ip)
            tmpa=pxi(lxisz)
            tmpb=pxi(lxit-lxisz)
            do m=1,2
               select case(m)
               case(1)
                  is=0
                  ie=lxie0
                  ii=is
                  xo=xb-xa
                  yo=yb
               case(2)
                  is=lxis2
                  ie=lxit
                  ii=ie
                  xo=xd-xc
                  yo=yc
               end select
               am=yo/sin(halfpi)**two
               do i=is,ie
                  gf=(-1)**(m+1)*(xp(i,n)-xp(ii,n))/xo
                  yp(i,n)=am*sin(halfpi*gf)**two
               end do
               do i=is,ie
                  xp(i,n+1)=xp(i,n)
                  yp(i,n+1)=yp(i,n)
               end do
            end do

!----- TOP & BOTTOM BOUNDARY POINTS IN XI-DIRECTION

            n=0
            fctr=three
            sho=szth0/16
            if(nthick==0) then
               do i=0,lxit
                  xp(i,n)=xp(i,n+1)
               end do
            else
               ip=0
               im=lxi0
               call gridf(xp(:,n),pxi,xa,xb-spx,sho,fctr*shs,lxit,im,ip)
               ip=lxis1
               im=lxi1
               call gridf(xp(:,n),pxi,xb-spx,xc+spx,fctr*shs,fctr*she,lxit,im,ip)
               ip=lxis2
               im=lxi2-lxisz
               call gridf(xp(:,n),pxi,xc+spx,xd,fctr*she,free,lxit,im,ip)
               ip=ip+im
               im=lxisz
               call gridf(xp(:,n),pxi,xd,xe,pxi(ip),free,lxit,im,ip)
            end if
            do i=0,lxit
               yp(i,n)=ya
               xp(i,n+3)=xp(i,n)
               yp(i,n+3)=yd
            end do

!----- VERTICAL INTERFACE POINTS IN ETA-DIRECTION

            if(nthick==0) then
               fctr=smgvr
            else
               fctr=smgvr*sqrt(half)
            end if
            do m=1,2
               select case(m)
               case(1)
                  yo=yb
                  sho=fctr*shs
               case(2)
                  yo=yc
                  sho=fctr*she*sqrt(two)
               end select
               jj=3*let0/4
               ra0=onethird*domh
               jp=let0-jj
               jm=jj
               call gridf(yq(:,m),qet,yo-ra0,yo,free,sho,lett,jm,jp)
               jp=0
               jm=let0-jj
               call gridf(yq(:,m),qet,ya,yo-ra0,free,qet(jp+jm),lett,jm,jp)
               jp=lets1
               jm=jj
               call gridf(yq(:,m),qet,yo,yo+ra0,sho,free,lett,jm,jp)
               jp=jp+jm
               jm=let0-jj
               call gridf(yq(:,m),qet,yo+ra0,yd,qet(jp),free,lett,jm,jp)
               qeto(:,m)=qet(:)
            end do
            if(nthick==0) then
               fctr=zero
               do j=0,lett
                  xq(j,1)=xb+fctr*(yq(j,1)-ya)
                  xq(j,2)=xc+fctr*(yq(j,2)-ya)
               end do
            else
               do n=0,2,2
                  select case(n)
                  case(0)
                     js=0
                     je=lete0
                  case(2)
                     js=lets1
                     je=lett
                  end select
                  do m=1,2
                     res=zero
                     select case(m)
                     case(1)
                        i=lxis1
                        pp=(1-n/2)*res-(n/2)*one
                        qq=(1-n/2)*one+(n/2)*res
                     case(2)
                        i=lxie1
                        pp=(1-n/2)*res+(n/2)*one
                        qq=(n/2-1)*one+(n/2)*res
                     end select
                     tmpa=xp(i,n+1)-xp(i,n)
                     tmpb=yp(i,n+1)-yp(i,n)
                     ra1=pp
                     ra2=(three*tmpa-(two*pp+qq)*tmpb)/tmpb**two
                     ra3=((pp+qq)*tmpb-two*tmpa)/tmpb**three
                     do j=js,je
                        tmp=yq(j,m)-yp(i,n)
                        xq(j,m)=xp(i,n)+ra1*tmp+ra2*tmp**two+ra3*tmp**three
                     end do
                  end do
               end do
            end if

!----- LEFT & RIGHT BOUNDARY POINTS IN ETA-DIRECTION

            if(nthick==0) then
               fctr=two*smgvr
            else
               fctr=two*smgvr*sqrt(half)
            end if
            do m=0,3,3
               select case(m)
               case(0)
                  yo=zero
                  sho=fctr*shs
               case(3)
                  yo=zero
                  sho=fctr*she*sqrt(two)
               end select
               jj=3*let0/4
               ra0=onethird*domh
               jp=let0-jj
               jm=jj
               call gridf(yq(:,m),qet,yo-ra0,yo,free,sho,lett,jm,jp)
               jp=0
               jm=let0-jj
               call gridf(yq(:,m),qet,ya,yo-ra0,free,qet(jp+jm),lett,jm,jp)
               jp=lets1
               jm=jj
               call gridf(yq(:,m),qet,yo,yo+ra0,sho,free,lett,jm,jp)
               jp=jp+jm
               jm=let0-jj
               call gridf(yq(:,m),qet,yo+ra0,yd,qet(jp),free,lett,jm,jp)
            end do
            fctr=zero
            do j=0,lett
               xq(j,0)=xa+fctr*(yq(j,0)-ya)
               xq(j,3)=xe+fctr*(yq(j,3)-ya)
            end do

!----- GRID OUTPUT

            do n=0,2,2
               select case(n)
               case(0)
                  js=0
                  je=lete0
               case(2)
                  js=lets1
                  je=lett
               end select
               do m=0,2
                  select case(m)
                  case(0)
                     is=0
                     ie=lxie0
                  case(1)
                     is=lxis1
                     ie=lxie1
                  case(2)
                     is=lxis2
                     ie=lxit
                  end select
                  do j=js,je
                     xx((/is,ie/),j)=xq(j,(/m,m+1/))
                     yy((/is,ie/),j)=yq(j,(/m,m+1/))
                  end do
                  do i=is,ie
                     xx(i,(/js,je/))=xp(i,(/n,n+1/))
                     yy(i,(/js,je/))=yp(i,(/n,n+1/))
                  end do
                  xco(:)=(/xx(is,js),xx(ie,js),xx(ie,je),xx(is,je)/)
                  yco(:)=(/yy(is,js),yy(ie,js),yy(ie,je),yy(is,je)/)
                  do j=js+1,je-1
                     do i=is+1,ie-1
                        xcij(:,1)=(/xp(i,n),xp(i,n),xp(i,n+1),xp(i,n+1)/)-xco(:)
                        xcij(:,2)=(/xq(j,m),xq(j,m+1),xq(j,m+1),xq(j,m)/)-xco(:)
                        ycij(:,1)=(/yp(i,n),yp(i,n),yp(i,n+1),yp(i,n+1)/)-yco(:)
                        ycij(:,2)=(/yq(j,m),yq(j,m+1),yq(j,m+1),yq(j,m)/)-yco(:)
                        xt(1:4,1)=(/(i-is)*(j-js),(ie-i)*(j-js),(ie-i)*(je-j),(i-is)*(je-j)/)**two
                        xt(1,2)=xt(2,1)*xt(3,1)*xt(4,1)
                        xt(2,2)=xt(3,1)*xt(4,1)*xt(1,1)
                        xt(3,2)=xt(4,1)*xt(1,1)*xt(2,1)
                        xt(4,2)=xt(1,1)*xt(2,1)*xt(3,1)
                        res=one/sum(xt(1:4,2))
                        cha(1:4)=xcij(:,1)+xcij(:,2)+xco(:)
                        dha(1:4)=ycij(:,1)+ycij(:,2)+yco(:)
                        xx(i,j)=res*sum(cha(1:4)*xt(1:4,2))
                        yy(i,j)=res*sum(dha(1:4)*xt(1:4,2))
                     end do
                  end do
               end do
            end do
 
            l=-3
            do j=js,je
               do i=is,ie
                  l=l+3
                  xyzmb(l:l+2)=(/xx(i,j),yy(i,j),zs(k)/)
               end do
            end do
            write(9,rec=k+1) xyzmb(:)
         end do
         close(9)

      end if
      deallocate(xx,yy,zz,xyzmb,xp,yp,xq,yq,pxi,qet)

   end subroutine gridaerofoil


   subroutine gridtoy(mb, lxio, leto, mo, nblocks, doml0, doml1, domh, span)
      integer(kind=ni), intent(in) :: mb, lxio, leto, nblocks
      integer(kind=ni), intent(in), dimension(0:nblocks) :: mo
      real(kind=nr),    intent(in) :: doml0, doml1, domh, span
      integer(kind=ni)             :: np, ll, k, js, je
      integer(kind=ni)             :: l, j, i, is, ie
      integer(kind=ni)             :: lxit, lett, lxie0, lxis1, lxie1, lxis2, lete0, lets1
      real(kind=nr)                :: ra0, ra1, ra2, ra3
      integer(kind=ni)             :: nrecd, myid

      inquire(iolength=ll) real(1.0,kind=ieee64); nrecd=ll
      myid = p_get_process_ID()

      lxit=lxi0+1
      lett=let0+1
      lxie0=lxi0
      lxis1=lxie0+1
      lxie1=lxis1+lxi1
      lxis2=lxie1+1
      lete0=let0
      lets1=lete0+1

      np=3*(lxio+1)*(leto+1)-1
      allocate(xx(0:lxit,0:lett),yy(0:lxit,0:lett))
      allocate(xyzmb(0:np), zs(0:lze0))

      if(myid==mo(mb)) then
         open(9,file=cgrid,status='unknown')
         close(9,status='delete')
         open(9,file=cgrid,access='direct',form='unformatted',recl=nrecd*(np+1),status='new')

         do k=0,lze0
            zs(k)=span*(real(lze0-k,kind=nr)/lze0-half)

         select case(mb)
            case(0)
               js=0
               je=lete0
               ra3=0
            end select
            select case(mb)
            case(0)
               is=0
               ie=lxie0
               ra2=0
            end select
 
            ra0 = (doml1+doml0)/(lxi0)
            ra1 = 2.0_nr*domh/(let0)
            l=-3
            do j=js,je
               do i=is,ie
                  l=l+3
                  xx(i,j) = -doml0 + (i-ra2)*ra0
                  yy(i,j) = -domh + (j-ra3)*ra1
                  xyzmb(l:l+2)=(/xx(i,j),yy(i,j),zs(k)/)
               end do
            end do
            write(9,rec=k+1) xyzmb(:)
         end do
         close(9)

      end if
      deallocate(xx, yy, xyzmb, zs)

   end subroutine gridtoy

   subroutine gridtoy_multi(mb, lxio, leto, mo, nblocks, doml0, doml1, domh, span)
      integer(kind=ni), intent(in) :: mb, lxio, leto, nblocks
      integer(kind=ni), intent(in), dimension(0:nblocks) :: mo
      real(kind=nr),    intent(in) :: doml0, doml1, domh, span
      integer(kind=ni)             :: np, ll, k, js, je
      integer(kind=ni)             :: l, j, i, is, ie
      integer(kind=ni)             :: lxit, lett, lxie0, lxis1, lxie1, lxis2, lete0, lets1
      real(kind=nr)                :: ra0, ra1, ra2, ra3
      integer(kind=ni)             :: nrecd, myid

      inquire(iolength=ll) real(1.0,kind=ieee64); nrecd=ll
      myid = p_get_process_ID()

      lxit=lxi0+lxi1+lxi2+2
      lett=2*let0+1
      lxie0=lxi0
      lxis1=lxie0+1
      lxie1=lxis1+lxi1
      lxis2=lxie1+1
      lete0=let0
      lets1=lete0+1

      np=3*(lxio+1)*(leto+1)-1
      allocate(xx(0:lxit,0:lett),yy(0:lxit,0:lett))
      allocate(xyzmb(0:np), zs(0:lze0))

      if(myid==mo(mb)) then
         open(9,file=cgrid,status='unknown')
         close(9,status='delete')
         open(9,file=cgrid,access='direct',form='unformatted',recl=nrecd*(np+1),status='new')

         do k=0,lze0
            zs(k)=span*(real(lze0-k,kind=nr)/lze0-half)

         select case(mb)
            case(0,1,2)
               js=0
               je=lete0
               ra3=0
            case(3,4,5)
               js=lets1
               je=lett
               ra3=1
            end select
            select case(mb)
            case(0,3)
               is=0
               ie=lxie0
               ra2=0
            case(1,4)
               is=lxis1
               ie=lxie1
               ra2=1
            case(2,5)
               is=lxis2
               ie=lxit
               ra2=2;
            end select

            ra0 = (doml1+doml0)/(lxi0+lxi1+lxi2)
            ra1 = domh/(2.0*let0)*2.0
            l=-3
            do j=js,je
               do i=is,ie
                  l=l+3
                  xx(i,j) = -doml0 + (i-ra2)*ra0
                  yy(i,j) = -domh + (j-ra3)*ra1
                  xyzmb(l:l+2)=(/xx(i,j),yy(i,j),zs(k)/)
               end do
            end do
            write(9,rec=k+1) xyzmb(:)
         end do
         close(9)

      end if
      deallocate(xx, yy, xyzmb, zs)

   end subroutine gridtoy_multi

!===== AEROFOIL INTERPOLATION

   subroutine aerofint(i,ii,n,nn,shg)

      integer(kind=ni), intent(in) :: i,ii,n,nn
      real(kind=nr),    intent(in) :: shg
      real(kind=nr)                :: err,fshg0,fshg1
      real(kind=nr)                :: res, ra0, ra1, ra2, ra3

      err=1.0e-4_nr
      fshg0=0.8_nr
      fshg1=(two-fshg0)/fshg0
      select case(nn)
      case(1)
         ra0=fshg0*(xp(ii,n)-xp(2*ii-i,n))
         xp(i,n)=xp(ii,n)+ra0
         yp(i,n)=xylagran(i,n,1)
      case(2)
         ra0=(-1)**n*shg*fshg0
         yp(i,n)=yp(ii,n)+ra0
         xp(i,n)=xylagran(i,n,2)
      end select
      ra1=fshg1*ra0
      ra2=errfn(i,ii,n,shg)
      do while(abs(ra2)>err)
         select case(nn)
         case(1)
            xp(i,n)=xp(ii,n)+ra1
            yp(i,n)=xylagran(i,n,1)
         case(2)
            yp(i,n)=yp(ii,n)+ra1
            xp(i,n)=xylagran(i,n,2)
         end select
         ra3=errfn(i,ii,n,shg)
         res=ra1-ra3*(ra1-ra0)/(ra3-ra2)
         ra0=ra1
         ra1=res
         ra2=ra3
      end do

   end subroutine aerofint

!===== ERROR FUNCTION

   function errfn(i,ii,n,shg) result(err)

      integer(kind=ni),intent(in) :: i,ii,n
      real(kind=nr),intent(in) :: shg
      real(kind=nr) :: err

      err=sqrt((xp(i,n)-xp(ii,n))**two+(yp(i,n)-yp(ii,n))**two)/shg-one

   end function errfn

!===== LAGRANGIAN INTERPOLATION

   function xylagran(i,n,nn) result(xy)

      integer(kind=ni),intent(in) :: i,n,nn
      integer(kind=ni) :: is,ie,ii,ip,jj
      real(kind=nr),dimension(-1:lnaca+2) :: xlag,ylag
      real(kind=nr),dimension(-1:2) :: alag,blag
      real(kind=nr) :: ao,bo,xy

      is=0
      ie=lnaca
      select case(nn)
      case(1)
         ii=minloc((xnaca(is:ie-1,n)-xp(i,n))*(xnaca(is+1:ie,n)-xp(i,n)),1)-1
         ip=min(max(ii,is+1),ie-2)
         xy=zero
         alag(:)=xp(i,n)-xnaca(ip-1:ip+2,n)
         do jj=-1,2
            blag(:)=xnaca(ip+jj,n)-xnaca(ip-1:ip+2,n)
            ao=one
            bo=one
            do ii=-1,2
               if(ii/=jj) then
                  ao=ao*alag(ii)
                  bo=bo*blag(ii)
               end if
            end do
            xy=xy+ao*ynaca(ip+jj,n)/bo
         end do
      case(2)
         ii=3
         ao=abs(xp(i,n)-xnaca(is,n))
         bo=abs(xp(i,n)-xnaca(ie,n))
         if(ao<bo) then
            ip=minloc((ynaca(is:is+ii-1,n)-yp(i,n))*(ynaca(is+1:is+ii,n)-yp(i,n)),1)-1
         else
            ip=minloc((ynaca(ie-ii+1:ie,n)-yp(i,n))*(ynaca(ie-ii:ie-1,n)-yp(i,n)),1)-1+ie-ii
         end if
         xlag(is:ie)=xnaca(:,n)
         xlag(is-1)=xnaca(is+1,3-n)
         xlag((/ie+1,ie+2/))=xnaca((/ie-1,ie-2/),3-n)
         ylag(is:ie)=ynaca(:,n)
         ylag(is-1)=ynaca(is+1,3-n)
         ylag((/ie+1,ie+2/))=ynaca((/ie-1,ie-2/),3-n)
         xy=zero
         alag(:)=yp(i,n)-ylag(ip-1:ip+2)
         do jj=-1,2
            blag(:)=ylag(ip+jj)-ylag(ip-1:ip+2)
            ao=one
            bo=one
            do ii=-1,2
               if(ii/=jj) then
                  ao=ao*alag(ii)
                  bo=bo*blag(ii)
               end if
            end do
            xy=xy+ao*xlag(ip+jj)/bo
         end do
      end select

   end function xylagran

!===== SUBROUTINE FOR GRID LINE GENERATION

   subroutine gridf(x,xxi,xo,xn,dxo,dxn,lxi,mxin,ip)

      integer(kind=ni),intent(in) :: lxi,mxin,ip
      real(kind=nr),dimension(0:lxi),intent(inout) :: x,xxi
      real(kind=nr),intent(in) :: xo,xn,dxo,dxn

      integer(kind=ni) :: i,ii
      real(kind=nr) :: dxoo,dxnn,aa,bb,cc,dd,xi,fctr

      dxoo=dxo
      dxnn=dxn
      xi=mxin
      if(int((dxo+sml)/free,kind=ni)==1) then
         dxoo=two*(xn-xo)/xi-dxnn
      end if
      if(int((dxn+sml)/free,kind=ni)==1) then
         dxnn=two*(xn-xo)/xi-dxoo
      end if
      aa=6.0_nr*(xn-xo)-3.0_nr*xi*(dxoo+dxnn)
      bb=15.0_nr*(xo-xn)+xi*(8.0_nr*dxoo+7.0_nr*dxnn)
      cc=10.0_nr*(xn-xo)-xi*(6.0_nr*dxoo+4.0_nr*dxnn)
      dd=xi*dxoo
      fctr=one/xi
      do i=0,mxin
         ii=i+ip
         xi=i*fctr
         x(ii)=aa*xi**five+bb*xi**four+cc*xi**three+dd*xi+xo
         xxi(ii)=fctr*(five*aa*xi**four+four*bb*xi**three+three*cc*xi**two+dd)
      end do

   end subroutine gridf
   
 end module mo_gridgen