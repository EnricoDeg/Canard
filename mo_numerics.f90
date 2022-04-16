!*****
!***** 3D SOLVER MODULE
!*****

module mo_numerics

   use mo_mpi, ONLY : p_null_req, p_isend, p_irecv, p_waitall, &
                      p_send, myid
   use mainvar3d
   use mo_utils
   implicit none

   private :: fcbcm, fcint, sbcco

   integer(kind=ni), private, parameter :: lmd=11,lmf=11,lmp=max(lmd,lmf)
   real(kind=nr),    private :: alphf, betf
   real(kind=nr),    private, dimension(0:lmp,0:1,0:1) :: pbci,pbco
   real(kind=nr),    private, dimension(0:1,0:1) :: pbcot
   real(kind=nr),    private :: fa,fb,fc
   real(kind=nr),    private, dimension(-2:2,0:2,0:1) :: albef
   real(kind=nr),    private, dimension(0:lmp) :: sap
   real(kind=nr),    private, dimension(0:4,0:2) :: fbc
   integer(kind=ni), private, dimension(3,0:1,0:1) :: ndf

   real(kind=nr),    private, dimension(:,:), allocatable :: xu,yu
   real(kind=nr),    private, dimension(:,:), allocatable :: xl,yl
   character(16),    private :: ccinput
   real(kind=nr),    private :: fltk,fltrbc
   

   contains

   SUBROUTINE allocate_numerics(limk)
      integer(kind=ni),intent(in) :: limk

      allocate(xu(0:limk,3),yu(0:limk,3),xl(0:limk,2),yl(0:limk,2))

   END SUBROUTINE allocate_numerics

   SUBROUTINE read_input_numerics

      open(9,file='input.numerics',status='old')
      read(9,*) ccinput,fltk
      read(9,*) ccinput,fltrbc
      close(9)
      fltk=pi*fltk

   END SUBROUTINE read_input_numerics

!===== EXTRA COEFFICIENTS FOR DOMAIN BOUNDARIES

   SUBROUTINE init_extracoeff_bounds
      integer(kind=ni) :: ntk, jjk, iik

      call fcbcm(fltk,fltrbc)
      call fcint(fltk,half,alphf,betf,fa,fb,fc)
      albef(:,0,1) = (/ zero,zero, one,alphf,betf /)
      albef(:,1,1) = (/ zero,alphf,one,alphf,betf /)
      albef(:,2,1) = (/ betf,alphf,one,alphf,betf /)

      pbco(:,:,:) = zero
      pbci(:,:,:) = zero
      call sbcco
      do ntk = 0,1
         do jjk = 0,1
            iik = lmd + ntk * ( lmf - lmd )
            pbcot(jjk,ntk) = sum( pbco(0:iik,jjk,ntk) )
         end do
      end do

   END SUBROUTINE init_extracoeff_bounds

   SUBROUTINE init_penta(lxik, letk, lzek, nbck)
      integer(kind=ni), intent(in) :: lxik, letk, lzek
      integer(kind=ni), intent(in), dimension(3,0:1) :: nbck
      integer(kind=ni) :: nnk, ipk, istart, iend
      integer(kind=ni) :: nstart, nend, npk

      do nnk = 1,3
         select case(nnk)
         case(1)
            istart = 0
            iend = istart + lxik
         case(2)
            istart = lxik + 1
            iend = istart + letk
         case(3)
            istart = lxik + letk + 2
            iend = istart + lzek
         end select
         do ipk = 0,1
            npk=nbck(nnk,ipk)
            select case(npk)
            case(10,20,25,30)
               ndf(nnk,ipk,0) = 0
               ndf(nnk,ipk,1) = 0
            case(35,40,45)
               ndf(nnk,ipk,0) = 1
               ndf(nnk,ipk,1) = 1
            end select
         end do
         nstart = ndf(nnk,0,0)
         nend   = ndf(nnk,1,0)
         call penta(xu(:,:), xl(:,:), istart, iend, nstart, nend, 0)
         nstart = ndf(nnk,0,1)
         nend   = ndf(nnk,1,1)
         call penta(yu(:,:), yl(:,:), istart, iend, nstart, nend, 1)
      end do
  
   END SUBROUTINE init_penta

!===== SUBROUTINE FOR CHOLESKY DECOMPOSITION OF PENTADIAGONAL MATRICES

   subroutine penta(xu,xl,is,ie,ns,ne,nt)

      integer(kind=ni),intent(in) :: is,ie,ns,ne,nt
      real(kind=nr),dimension(0:lim,3),intent(inout) :: xu
      real(kind=nr),dimension(0:lim,2),intent(inout) :: xl
      real(kind=nr),dimension(-2:2,0:2,0:1) :: albe
      real(kind=nr) :: alpho,beto
      integer(kind=ni) :: iik

      if(nt==0) then
         albe(:,0,0) = (/ zero, zero,    one, alpha01, beta02 /)
         albe(:,1,0) = (/ zero, alpha10, one, alpha12, beta13 /)
         albe(:,2,0) = (/ beta, alpha,   one, alpha,   beta /)
         alpho = alpha
         beto  = beta
         albe(:,0,1) = (/ zero, zero,  one, alpha, beta /)
         albe(:,1,1) = (/ zero, alpha, one, alpha, beta /)
         albe(:,2,1) = (/ beta, alpha, one, alpha, beta /)
      else
         albe  = albef
         alpho = alphf
         beto  = betf
      end if

      do iik = is,ie
         xl(iik,:) = one
         xu(iik,:) = one
      end do
      iik = is
      xu(iik,1) = one
      xu(iik,2) = albe(1,0,ns)
      xu(iik,3) = albe(2,0,ns)
      iik = is + 1
      xl(iik,2) = albe(-1,1,ns) * xu(iik-1,1)
      xu(iik,1) = one / ( one - xu(iik-1,2) * xl(iik,2) )
      xu(iik,2) = albe(1,1,ns) - xu(iik-1,3) * xl(iik,2)
      xu(iik,3) = albe(2,1,ns)
      iik = is + 2
      xl(iik,1) = albe(-2,2,ns) * xu(iik-2,1)
      xl(iik,2) = ( albe(-1,2,ns) - xu(iik-2,2) * xl(iik,1) ) * xu(iik-1,1)
      xu(iik,1) = one / ( one - xu(iik-2,3) * xl(iik,1) - xu(iik-1,2) * xl(iik,2))
      xu(iik,2) = albe(1,2,ns) - xu(iik-1,3) * xl(iik,2)
      xu(iik,3) = albe(2,2,ns)
      do iik = is+3,ie-3
         xl(iik,1) = beto * xu(iik-2,1)
         xl(iik,2) = ( alpho - xu(iik-2,2) * xl(iik,1) ) * xu(iik-1,1)
         xu(iik,1) = one / ( one - xu(iik-2,3) * xl(iik,1) - xu(iik-1,2) * xl(iik,2) )
         xu(iik,2) = alpho - xu(iik-1,3) * xl(iik,2)
         xu(iik,3) = beto
      end do
      iik = ie - 2
      xl(iik,1) = albe(2,2,ne) * xu(iik-2,1)
      xl(iik,2) = ( albe(1,2,ne) - xu(iik-2,2) * xl(iik,1) ) * xu(iik-1,1)
      xu(iik,1) = one / ( one - xu(iik-2,3) * xl(iik,1) - xu(iik-1,2) * xl(iik,2) )
      xu(iik,2) = albe(-1,2,ne) - xu(iik-1,3) * xl(iik,2)
      xu(iik,3) = albe(-2,2,ne)
      iik = ie - 1
      xl(iik,1) = albe(2,1,ne) * xu(iik-2,1)
      xl(iik,2) = ( albe(1,1,ne) - xu(iik-2,2) * xl(iik,1) ) * xu(iik-1,1)
      xu(iik,1) = one / ( one - xu(iik-2,3) * xl(iik,1) - xu(iik-1,2) * xl(iik,2) )
      xu(iik,2) = albe(-1,1,ne) - xu(iik-1,3) * xl(iik,2)
      iik = ie
      xl(iik,1) = albe(2,0,ne) * xu(iik-2,1)
      xl(iik,2) = ( albe(1,0,ne) - xu(iik-2,2) * xl(iik,1) ) * xu(iik-1,1)
      xu(iik,1) = one / ( one - xu(iik-2,3) * xl(iik,1) - xu(iik-1,2) * xl(iik,2) )
      do iik = is,ie
         xu(iik,2:3) = xu(iik,2:3) * xu(iik,1)
      end do

   end subroutine penta

!===== SUBROUTINE FOR BOUNDARY FILTER COEFFICIENTS

   subroutine fcbcm(fltk,fltrbc)
 
      real(kind=nr),intent(in) :: fltk,fltrbc
      real(kind=nr) :: alphz,betz,za,zb,zc
      real(kind=nr) :: aok, fctrk, resk

      aok = log(fltrbc)
      call fcint(fltk, half, alphz, betz, za, zb, zc)
      fctrk = one / &
            ( one + alphz * fltrbc + betz * fltrbc**two)

      albef(:,0,0) = (/ zero,          &
                        zero,          &
                        one,           &
                        alphz * fctrk, &
                        betz * fctrk   /)
      resk = (fltrbc-1) * &
            (za + zc + (fltrbc+1) * (zb+fltrbc*zc)) / aok
      fbc(:,0) = (/ za - 5 * resk / 3,   &
                    zb + 10 * resk / 21, &
                    zc - 5 * resk / 42,  &
                    5 * resk / 252,      &
                    -resk / 630          /) * fctrk

      albef(:,1,0) = (/ zero,                  &
                        alphz + betz * fltrbc, &
                        one,                   &
                        alphz,                 &
                        betz                   /)
      resk = (fltrbc-1) * &
            (zb + zc * (fltrbc+1)) / aok
      fbc(:,1) = (/ za + zb + zc + 1627 * resk / 1260, &
                    za + 10 * resk / 21,               &
                    zb - 5 * resk / 42,                &
                    zc + 5 * resk / 252,               &
                    -resk / 630                        /)

      albef(:,2,0) = (/ betz,  &
                        alphz, &
                        one,   &
                        alphz, &
                        betz   /)
      resk = zc * (fltrbc-1) / aok
      fbc(:,2) = (/ zb + zc + 1627 * resk / 1260, &
                    za - 5 * resk / 3,            &
                    za - 5 * resk / 42,           &
                    zb + 5 * resk / 252,          &
                    zc - resk / 630               /)

   end subroutine fcbcm

!===== SUBROUTINE FOR INTERIOR FILTER COEFFICIENTS

   subroutine fcint(fltk,fltr,alphz,betz,za,zb,zc)
 
      real(kind=nr),intent(in) :: fltk,fltr
      real(kind=nr),intent(inout) :: alphz,betz,za,zb,zc
      real(kind=nr),dimension(3) :: cosf
      real(kind=nr) :: fctrk

      cosf(1) = cos(fltk  )
      cosf(2) = cos(2*fltk)
      cosf(3) = cos(3*fltk)
      fctrk = 1 / ( 30 + 5 * (7-16*fltr) * cosf(1) + &
                   2 * (1+8*fltr) * cosf(2) - 3 * cosf(3))
      alphz = fctrk * (20 * (2*fltr-1) - 30 * cosf(1) + &
                      12 * (2 * fltr - 1) * cosf(2) - 2 * cosf(3))
      betz = half * fctrk * (2 * (13 - 8 * fltr) + &
             (33 - 48 * fltr) * cosf(1) + 6 * cosf(2) - cosf(3))
      za = 60 * (1-fltr) * cos(half*fltk)**4 * fctrk
      zb = -2 * za / 5
      zc = za / 15

   end subroutine fcint

!===== SUBROUTINE FOR SUBDOMAIN-BOUNDARY COEFFICIENTS

   subroutine sbcco

      real(kind=nr),dimension(:,:),allocatable :: ax,bx,rx,sx
      integer(kind=ni) :: ntk, llk, iik, istart, iend

      do ntk = 0,1

         if ( ntk == 0 ) then
            llk = lmd
            istart = 1
            iend = 2 * ( llk + 1 )
            allocate( ax(iend,iend), bx(iend,iend), &
                      rx(iend,iend), sx(iend,iend) )
            ax(:,:) = 0
            bx(:,:) = 0
            ax(istart,istart:istart+2) = (/ one, alpha01, beta02 /)
            bx(istart,istart:istart+4) = (/ -( a01 + a02 + a03 + a04 ), &
                                               a01,                     &
                                               a02,                     &
                                               a03,                     &
                                               a04                      /)
            ax(istart+1,istart:istart+3) = (/ alpha10, one, alpha12, beta13 /)
            bx(istart+1,istart:istart+4) = (/ a10,                      &
                                           -( a10 + a12 + a13 + a14 ),  &
                                              a12,                      &
                                              a13,                      &
                                              a14                       /)
            do iik = istart+2,iend-2
               ax(iik,iik-2:iik+2) = (/ beta, alpha, one, alpha, beta /)
               bx(iik,iik-2:iik+2) = (/ -ab, -aa, zero, aa, ab /)
            end do
            ax(iend-1,iend:iend-3:-1) =  ax(istart+1,istart:istart+3)
            bx(iend-1,iend:iend-4:-1) = -bx(istart+1,istart:istart+4)
            ax(iend,iend:iend-2:-1)   =  ax(istart,istart:istart+2)
            bx(iend,iend:iend-4:-1)   = -bx(istart,istart:istart+4)
         end if

         if ( ntk == 1 ) then
            llk = lmf
            istart = 1
            iend   = 2 * ( llk + 1 )
            allocate( ax(iend,iend), bx(iend,iend), &
                      rx(iend,iend), sx(iend,iend) )
            ax(:,:)                   = 0
            bx(:,:)                   = 0
            ax(istart,istart:istart+2)        = albef(0:2,0,0)
            bx(istart,istart+(/1,2,3,4,5/))   = fbc(:,0)
            bx(istart,istart)                 = -sum( fbc(:,0) )
            ax(istart+1,istart:istart+3)      = albef(-1:2,1,0)
            bx(istart+1,istart+(/0,2,3,4,5/)) = fbc(:,1)
            bx(istart+1,istart+1)             = -sum( fbc(:,1) )
            ax(istart+2,istart:istart+4)      = albef(-2:2,2,0)
            bx(istart+2,istart+(/0,1,3,4,5/)) = fbc(:,2)
            bx(istart+2,istart+2)             = -sum( fbc(:,2) )
            do iik = istart+3,iend-3
               ax(iik,iik-2:iik+2)          = (/ betf, alphf, one, alphf, betf /)
               bx(iik,iik-3:iik+3)          = (/ fc, fb, fa,            &
                                           -2 * ( fa + fb + fc ), &
                                           fa, fb, fc             /)
            end do
            ax(iend-2,iend:iend-4:-1) = ax(istart+2,istart:istart+4)
            bx(iend-2,iend:iend-5:-1) = bx(istart+2,istart:istart+5)
            ax(iend-1,iend:iend-3:-1) = ax(istart+1,istart:istart+3)
            bx(iend-1,iend:iend-5:-1) = bx(istart+1,istart:istart+5)
            ax(iend,iend:iend-2:-1)   = ax(istart,istart:istart+2)
            bx(iend,iend:iend-5:-1)   = bx(istart,istart:istart+5)
         end if

         call mtrxi( ax(:,:), sx(:,:), istart, iend )

         rx(:,:) = ax(:,:)

         iik = iend / 2 - 1
         rx(iik,iik+2)   = 0
         rx(iik+1,iik+2) = 0
         rx(iik+1,iik+3) = 0
         
         iik = iend / 2 + 2
         rx(iik,iik-2)   = 0
         rx(iik-1,iik-2) = 0
         rx(iik-1,iik-3) = 0
         
         ax(:,:) = matmul( rx(:,:), matmul( sx(:,:), bx(:,:) ) )
         
         iik = iend / 2 + 1
         pbco(llk:0:-1,0,ntk) = ax(iik,istart:istart+llk)
         pbci(0:llk,0,ntk)    = ax(iik,istart+llk+1:iend)
         
         iik = iend / 2 + 2
         pbco(llk:0:-1,1,ntk) = ax(iik,istart:istart+llk)
         pbci(0:llk,1,ntk)    = ax(iik,istart+llk+1:iend)
         
         deallocate( ax, bx, rx, sx )
           
      end do
   end subroutine sbcco

!===== SUBROUTINE FOR MPI IMPLEMENTATION

   subroutine mpigo(nt,nrt,n45,itag)

      integer(kind=ni),intent(in) :: nt,nrt,n45,itag

      select case(nt)
      case(0)
         mp=lmd
      case(1)
         mp=lmf
      end select

      call p_null_req
      do nn=1,3
         nz=(1-nrt)*(nn-1)+1
         if(nt==0) then
            select case(nn)
            case(1)
               send=>send01
               recv=>recv01
            case(2)
               send=>send02
               recv=>recv02
            case(3)
               send=>send03
               recv=>recv03
            end select
         else
            select case(nn)
            case(1)
               send=>send11
               recv=>recv11
            case(2)
               send=>send12
               recv=>recv12
            case(3)
               send=>send13
               recv=>recv13
            end select
         end if
         do ip=0,1
            iq=1-ip
            is=ip*ijk(1,nn)
            ie=1-2*ip
            select case(nbc(nn,ip))
            case(35)
               ra0=zero
               ii=1
            case(40)
               ra0=zero
               ii=0
            case(45)
               ra0=n45
               ii=1
            end select
            if(ndf(nn,ip,nt)==1) then
               do k=0,ijk(3,nn)
                  kp=k*(ijk(2,nn)+1)
                  do j=0,ijk(2,nn)
                     jk=kp+j
                     l=indx3(is,j,k,nn)
                     res=ra0*rr(l,nz)
                     do i=0,mp
                        l=indx3(is+ie*(i+ii),j,k,nn)
                        sap(i)=rr(l,nz)
                     end do
                     send(jk,0,ip)=sum(pbco(0:mp,0,nt)*sap(0:mp))-res*pbcot(0,nt)
                     send(jk,1,ip)=sum(pbco(0:mp,1,nt)*sap(0:mp))-res*pbcot(1,nt)
                     send(jk,nt+1,ip)=send(jk,nt+1,ip)+nt*(sap(0)-res-send(jk,nt+1,ip))
                  end do
               end do
               if(nt==0) then
                  call p_isend(send(:,:,ip), mcd(nn,ip), itag+iq, 2*nbsize(nn))
                  call p_irecv(recv(:,:,ip), mcd(nn,ip), itag+ip, 2*nbsize(nn))
               else
                  call p_isend(send(:,:,ip), mcd(nn,ip), itag+iq, 3*nbsize(nn))
                  call p_irecv(recv(:,:,ip), mcd(nn,ip), itag+ip, 3*nbsize(nn))
               end if
            end if
         end do
      end do
      call p_waitall

      if(n45==n45go) then
         do nn=1,3
            nz=(1-nrt)*(nn-1)+1
            if(nt==0) then
               select case(nn)
               case(1)
                  recv=>recv01
               case(2)
                  recv=>recv02
               case(3)
                  recv=>recv03
               end select
            else
               select case(nn)
               case(1)
                  recv=>recv11
               case(2)
                  recv=>recv12
               case(3)
                  recv=>recv13
               end select
            end if
            do ip=0,1
               is=ip*ijk(1,nn)
               if(nbc(nn,ip)==45) then
                  do k=0,ijk(3,nn)
                     kp=k*(ijk(2,nn)+1)
                     do j=0,ijk(2,nn)
                        jk=kp+j
                        l=indx3(is,j,k,nn)
                        recv(jk,0,ip)=recv(jk,0,ip)+rr(l,nz)*pbcot(0,nt)
                        recv(jk,1,ip)=recv(jk,1,ip)+rr(l,nz)*pbcot(1,nt)
                        recv(jk,nt+1,ip)=recv(jk,nt+1,ip)+nt*rr(l,nz)
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

      nt=0
      ns=ndf(nn,0,0)
      ne=ndf(nn,1,0)

      select case(nn)
      case(1)
         is=0
         ie=is+lxi
         recv=>recv01
         drva=>drva1
      case(2)
         is=lxi+1
         ie=is+let
         recv=>recv02
         drva=>drva2
      case(3)
         is=lxi+let+2
         ie=is+lze
         recv=>recv03
         drva=>drva3
      end select

      do k=0,ijk(3,nn)
         kp=k*(ijk(2,nn)+1)
         do j=0,ijk(2,nn)
            jk=kp+j
            do i=is,ie
               l=indx3(i-is,j,k,nn)
               li(i)=l
               sa(i)=rr(l,nz)
            end do
            select case(ns)
            case(0)
               sb(is)=sum((/a01,a02,a03,a04/)*(sa(is+(/1,2,3,4/))-sa(is)))
               sb(is+1)=sum((/a10,a12,a13,a14/)*(sa(is+(/0,2,3,4/))-sa(is+1)))
            case(1)
               sb(is)=sum(pbci(0:lmd,0,nt)*sa(is:is+lmd))+recv(jk,0,0)
               sb(is+1)=sum(pbci(0:lmd,1,nt)*sa(is:is+lmd))+recv(jk,1,0)
            end select
            do i=is+2,ie-2
               sb(i)=aa*(sa(i+1)-sa(i-1))+ab*(sa(i+2)-sa(i-2))
            end do
            select case(ne)
            case(0)
               sb(ie)=sum((/a01,a02,a03,a04/)*(sa(ie)-sa(ie-(/1,2,3,4/))))
               sb(ie-1)=sum((/a10,a12,a13,a14/)*(sa(ie-1)-sa(ie-(/0,2,3,4/))))
            case(1)
               sb(ie)=-sum(pbci(0:lmd,0,nt)*sa(ie:ie-lmd:-1))-recv(jk,0,1)
               sb(ie-1)=-sum(pbci(0:lmd,1,nt)*sa(ie:ie-lmd:-1))-recv(jk,1,1)
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
               l=li(i)
               rr(l,nn)=sb(i)
            end do
            drva(jk,m,0)=sb(is)
            drva(jk,m,1)=sb(ie)
         end do
      end do

   end subroutine deriv

!===== SUBROUTINE FOR COMPACT FILTERING

   subroutine filte(nn,nz)

      integer(kind=ni),intent(in) :: nn,nz

      nt=1
      ns=ndf(nn,0,1)
      ne=ndf(nn,1,1)

      select case(nn)
      case(1)
         is=0
         ie=is+lxi
         recv=>recv11
      case(2)
         is=lxi+1
         ie=is+let
         recv=>recv12
      case(3)
         is=lxi+let+2
         ie=is+lze
         recv=>recv13
      end select

      do k=0,ijk(3,nn)
         kp=k*(ijk(2,nn)+1)
         do j=0,ijk(2,nn)
            jk=kp+j
            do i=is,ie
               l=indx3(i-is,j,k,nn)
               li(i)=l
               sa(i)=rr(l,nz)
            end do
            select case(ns)
            case(0)
               sb(is)=sum(fbc(:,0)*(sa(is+(/1,2,3,4,5/))-sa(is)))
               sb(is+1)=sum(fbc(:,1)*(sa(is+(/0,2,3,4,5/))-sa(is+1)))
               sb(is+2)=sum(fbc(:,2)*(sa(is+(/0,1,3,4,5/))-sa(is+2)))
            case(1)
               ra2=sa(is+2)+sa(is+2)
               sb(is)=sum(pbci(0:lmf,0,nt)*sa(is:is+lmf))+recv(jk,0,0)
               sb(is+1)=sum(pbci(0:lmf,1,nt)*sa(is:is+lmf))+recv(jk,1,0)
               sb(is+2)=fa*(sa(is+1)+sa(is+3)-ra2)+fb*(sa(is)+sa(is+4)-ra2)+fc*(recv(jk,2,0)+sa(is+5)-ra2)
            end select
            do i=is+3,ie-3
               res=sa(i)+sa(i)
               sb(i)=fa*(sa(i-1)+sa(i+1)-res)+fb*(sa(i-2)+sa(i+2)-res)+fc*(sa(i-3)+sa(i+3)-res)
            end do
            select case(ne)
            case(0)
               sb(ie)=sum(fbc(:,0)*(sa(ie-(/1,2,3,4,5/))-sa(ie)))
               sb(ie-1)=sum(fbc(:,1)*(sa(ie-(/0,2,3,4,5/))-sa(ie-1)))
               sb(ie-2)=sum(fbc(:,2)*(sa(ie-(/0,1,3,4,5/))-sa(ie-2)))
            case(1)
               ra2=sa(ie-2)+sa(ie-2)
               sb(ie)=sum(pbci(0:lmf,0,nt)*sa(ie:ie-lmf:-1))+recv(jk,0,1)
               sb(ie-1)=sum(pbci(0:lmf,1,nt)*sa(ie:ie-lmf:-1))+recv(jk,1,1)
               sb(ie-2)=fa*(sa(ie-3)+sa(ie-1)-ra2)+fb*(sa(ie-4)+sa(ie)-ra2)+fc*(sa(ie-5)+recv(jk,2,1)-ra2)
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
               l=li(i)
               rr(l,nz)=rr(l,nz)+sb(i)
            end do
         end do
      end do

   end subroutine filte

!=====

end module mo_numerics

!*****