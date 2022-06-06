!*****
!***** DOMAIN DECOMPOSITION MODULE
!*****

MODULE mo_domdcomp
   use mo_kind, ONLY : ni, nr
   use mo_parameters, ONLY : one
   use mo_vars, ONLY : ijk, lxim, letm, lzem, lximb, letmb, lzemb, &
                     & mo, nbsize, nbpc, mcd, nbc,   &
                     & lxio, leto, lzeo, lxi, let, lze, mb, ltomb, &
                     & lmx, lim, mbk
   use mo_mpi, ONLY : myid
   use mo_mpi, ONLY : p_null_req, p_irecv, p_isend, p_waitall, &
                      p_recv, p_send, p_bcast, mpro
   use mo_utils, ONLY : indx3
   IMPLICIT NONE
   PUBLIC

   integer(kind=ni), private, parameter                    :: ljpl=100

   integer(kind=ni), private, dimension(3,0:1)             :: mmcd
   integer(kind=ni), private, dimension(3)                 :: ijkp
   integer(kind=ni), private, dimension(0:3,3,0:1)         :: ijp
   integer(kind=ni), private, dimension(0:3,3,0:1)         :: ijl
   integer(kind=ni), private, dimension(0:11,ljpl)         :: jlcd
   integer(kind=ni), private, dimension(0:7,ljpl)          :: jpcd
   integer(kind=ni), private, dimension(0:7)               :: njp
   integer(kind=ni), private, dimension(0:11)              :: njl


   integer(kind=ni), private, dimension(:,:,:), allocatable :: nbbc, mbcd
   integer(kind=ni), private, dimension(:,:), allocatable   :: jjp
   integer(kind=ni), private, dimension(:,:), allocatable   :: jjl
   integer(kind=ni), private, dimension(:), allocatable     :: imjp
   integer(kind=ni), private, dimension(:), allocatable     :: jptag
   integer(kind=ni), private, dimension(:), allocatable     :: imjl
   integer(kind=ni), private, dimension(:), allocatable     :: jltag


   CONTAINS

   subroutine allocate_domdcomp(nblocks)
      integer(kind=ni), INTENT(IN) :: nblocks
      integer(kind=ni)             :: lpp, lqq

      lpp=8*nblocks+7
      lqq=12*nblocks+11
      allocate(nbbc(0:nblocks,3,0:1), mbcd(0:nblocks,3,0:1))
      allocate(imjp(0:lpp), jptag(0:lpp), imjl(0:lqq), jltag(0:lqq))
      allocate(jjp(0:lpp,0:ljpl))
      allocate(jjl(0:lqq,0:ljpl))
      allocate(lximb(0:nblocks), letmb(0:nblocks), lzemb(0:nblocks))
      allocate(mo(0:nblocks), nbpc(0:nblocks,3))

   end subroutine allocate_domdcomp

   subroutine domdcomp_init(nblocks, nkthick, nkbody)
      integer(kind=ni), intent(in) :: nblocks, nkthick
      integer(kind=ni), intent(in) :: nkbody
      integer(kind=ni) :: ipk, jpk, mmk, nnk, nstart, nend
      integer(kind=ni) :: llk, mpk, lpk, mak, lk, mm, mp, itag

      mo(0) = 0
      do mmk = 1,nblocks
         mo(mmk) = mo(mmk-1) + nbpc(mmk-1,1) * &
                               nbpc(mmk-1,2) * &
                               nbpc(mmk-1,3)
      end do
      
      do mmk = 0,nblocks
         if(myid >= mo(mmk)) then
            mb = mmk
         end if
      end do

      lxio = lximb(mb)
      leto = letmb(mb)
      lzeo = lzemb(mb)

      ! multiblock info
      ipk = 30 * nkthick + 35 * (1 - nkthick)
#ifdef VISCOUS
      jpk = 35 * (1 - nkbody) + nkbody * (20 + 5)
#else
      jpk = 35 * (1 - nkbody) + nkbody * (20)
#endif
      do mmk = 0,nblocks
         select case(mmk)
         case(0,3)
            nbbc(mmk,1,:) = (/ 10, ipk /)
            mbcd(mmk,1,:) = (/ -1, mmk+1 /)
         case(1,4)
            nbbc(mmk,1,:) = (/ ipk, ipk  /)
            mbcd(mmk,1,:) = (/ mmk-1, mmk+1 /)
         case(2,5)
            nbbc(mmk,1,:) = (/ipk,10/)
            mbcd(mmk,1,:) = (/mmk-1,-1/)
         end select

         select case(mmk)
         case(0,2)
            nbbc(mmk,2,:) = (/45,ipk/)
            mbcd(mmk,2,:) = (/mmk+3,mmk+3/)
         case(1)
            nbbc(mmk,2,:) = (/45,ipk/)
            mbcd(mmk,2,:) = (/mmk+3,mmk+3/)
         case(3,5)
            nbbc(mmk,2,:) = (/ipk,45/)
            mbcd(mmk,2,:) = (/mmk-3,mmk-3/)
         case(4)
            nbbc(mmk,2,:) = (/ipk,45/)
            mbcd(mmk,2,:) = (/mmk-3,mmk-3/)
         end select

         nbbc(mmk,3,:) = (/45,45/)
         mbcd(mmk,3,:) = (/mmk,mmk/)
      end do

      ! domdcomp initialize
      ijkp(1) = mod(myid-mo(mb),nbpc(mb,1))
      ijkp(2) = mod((myid-mo(mb))/nbpc(mb,1),nbpc(mb,2))
      ijkp(3) = mod((myid-mo(mb))/(nbpc(mb,1)*nbpc(mb,2)),nbpc(mb,3))

      do nnk = 1,3
         nstart = mod(nnk,3)+1
         nend   = mod(nstart,3)+1
         do ipk = 0,1
            mm = mbcd(mb,nnk,ipk)
            if(mm == -1) then
               mmcd(nnk,ipk) = -1
            else
               mmcd(nnk,ipk) = idsd3( (1 - ipk) * (nbpc(mm,nnk) - 1), &
                                      ijkp(nstart),                   &
                                      ijkp(nend),                     &
                                      mm,                             &
                                      nnk ) 
            end if
         end do
      end do

      do nnk = 1,3
         select case(nnk)
         case (1)
            llk = lxio
            mpk = 1
         case (2)
            llk = leto
            mpk = nbpc(mb,1)
         case (3)
            llk = lzeo
            mpk = nbpc(mb,1) * nbpc(mb,2)
         end select
         lpk = ijkp(nnk)
         mak = nbpc(mb,nnk)
         if(mak == 1) then
            lk = llk
            nbc(nnk,0) = nbbc(mb,nnk,0)
            nbc(nnk,1) = nbbc(mb,nnk,1)
            mcd(nnk,0) = mmcd(nnk,0)
            mcd(nnk,1) = mmcd(nnk,1)
         end if
         if(mak >= 2) then
            if(lpk == 0) then
               lk = llk - ( ( llk + 1 ) / mak ) * ( mak - 1 )
               nbc(nnk,0) = nbbc(mb,nnk,0)
               nbc(nnk,1) = 40
               mcd(nnk,0) = mmcd(nnk,0)
               mcd(nnk,1) = myid+mpk
            end if
            if(lpk > 0 .and. lpk < mak-1) then
               lk = ( llk + 1 ) / mak - 1
               nbc(nnk,0) = 40
               nbc(nnk,1) = 40
               mcd(nnk,0) = myid-mpk
               mcd(nnk,1) = myid+mpk
            end if
            if(lpk == mak-1) then
               lk = ( llk + 1 ) / mak - 1
               nbc(nnk,0) = 40
               nbc(nnk,1) = nbbc(mb,nnk,1)
               mcd(nnk,0) = myid - mpk
               mcd(nnk,1) = mmcd(nnk,1)
            end if
         end if
         select case(nnk)
         case (1)
            lxi=lk
         case (2)
            let=lk
         case (3)
            lze=lk
         end select
      end do

      if(myid==0) then
         lxim(0)=lxi
         letm(0)=let
         lzem(0)=lze
         do mp=1,mpro
            itag=1
            call p_recv(lxim(mp), mp, itag)
            itag=2
            call p_recv(letm(mp), mp, itag)
            itag=3
            call p_recv(lzem(mp), mp, itag)
         end do
      else
         itag=1
         call p_send(lxi, 0, itag)
         itag=2
         call p_send(let, 0, itag)
         itag=3
         call p_send(lze, 0, itag)
      end if
      call p_bcast(lxim(:), 0)
      call p_bcast(letm(:), 0)
      call p_bcast(lzem(:), 0)

      ltomb=(lxio+1)*(leto+1)*(lzeo+1)

      lmx=(lxi+1)*(let+1)*(lze+1)-1
      lim=(lxi+1)+(let+1)+(lze+1)-1

      ijk(1,1)=lxi
      ijk(2,1)=let
      ijk(3,1)=lze
      ijk(1,2)=let
      ijk(2,2)=lze
      ijk(3,2)=lxi
      ijk(1,3)=lze
      ijk(2,3)=lxi
      ijk(3,3)=let

      nbsize(:)=(ijk(2,:)+1)*(ijk(3,:)+1)

   END SUBROUTINE domdcomp_init

   SUBROUTINE search_point
      integer(kind=ni) :: is, ie, kk, jj, mm, mp, kp, nn
      integer(kind=ni) :: ll, l, jp, ip, ii, i, j, k
      integer(kind=ni) :: ks, ke, lp

      lp=8*mbk+7
      imjp(:)=0
      jjp(:,:)=-1
      do ip=0,1
         ijp(:,1,ip)=(/0,2,4,6/)+ip
         ijp(:,2,ip)=(/0,4,1,5/)+2*ip
         ijp(:,3,ip)=(/0,1,2,3/)+4*ip
      end do
      do i=0,lp
         jjp(i,0)=i
      end do
      do mm=0,mbk
         do nn=1,3
            do ip=0,1
               mp=mbcd(mm,nn,ip)
               if(mp/=-1) then
                  do i=0,3
                     is=8*mm+ijp(i,nn,ip)
                     ie=8*mp+ijp(i,nn,1-ip)
                     imjp(is)=imjp(is)+1
                     jjp(is,imjp(is))=jjp(ie,0)
                  end do
               end if
            end do
         end do
      end do
      do i=0,lp
         do j=0,lp
            ll=0
            if((i-j)*imjp(i)*imjp(j)/=0) then
               do k=0,imjp(i)
                  do kk=0,imjp(j)
                     if(jjp(j,kk)-jjp(i,k)==0) then
                        ll=ll+1
                     end if
                  end do
               end do
               if(ll/=0) then
                  ks=imjp(i)+1
                  ke=ks+imjp(j)
                  jjp(i,ks:ke)=jjp(j,0:ke-ks)
                  call trimm(jjp(i,:),imjp(i),ke)
               end if
            end if
         end do
      end do
      do i=0,lp
         jptag(i)=minval(jjp(i,0:imjp(i)))
      end do
      njp(:)=-1
      jpcd(:,:)=-1
      do i=0,7
         ip=mod(i,2)
         jp=mod(i,4)/2
         kp=i/4
         l=8*mb+i
         if( ijkp(1)==ip*(nbpc(mb,1)-1) .and. &
             ijkp(2)==jp*(nbpc(mb,2)-1) .and. &
             ijkp(3)==kp*(nbpc(mb,3)-1) .and. &
             imjp(l)>=3) then
            njp(i)=l
            do k=1,imjp(l)
               mm=jjp(l,k)/8
               j=jjp(l,k)-8*mm
               ii=mod(j,2)
               jj=mod(j,4)/2
               kk=j/4
               jpcd(i,k)=idsd3(ii*(nbpc(mm,1)-1),jj*(nbpc(mm,2)-1),kk*(nbpc(mm,3)-1),mm,1)
            end do
         end if
      end do

   END SUBROUTINE search_point

   SUBROUTINE search_line
      integer(kind=ni) :: ns, ne, is, ie, np, nq
      integer(kind=ni) :: kk, jj, mm, mp, nn, ll
      integer(kind=ni) :: l, jp, ip, ii, i, j, k
      integer(kind=ni) :: ks, ke, lp

      lp=12*mbk+11
      imjl(:)=0
      jjl(:,:)=-1
      njl(0:3)=(/2,2,1,1/)
      do ip=0,1
         ijl(:,1,ip)=(/0,1,8,10/)+ip*njl(0:3)
         ijl(:,2,ip)=(/4,5,0,2/)+ip*njl(0:3)
         ijl(:,3,ip)=(/8,9,4,6/)+ip*njl(0:3)
      end do
      do i=0,lp
         jjl(i,0)=i
      end do
      do mm=0,mbk
         do nn=1,3
            do ip=0,1
               mp=mbcd(mm,nn,ip)
               if(mp/=-1) then
                  do i=0,3
                     is=12*mm+ijl(i,nn,ip)
                     ie=12*mp+ijl(i,nn,1-ip)
                     imjl(is)=imjl(is)+1
                     jjl(is,imjl(is))=jjl(ie,0)
                  end do
               end if
            end do
         end do
      end do
      do i=0,lp
         do j=0,lp
            ll=0
            if((i-j)*imjl(i)*imjl(j)/=0) then
               do k=0,imjl(i)
                  do kk=0,imjl(j)
                     if(jjl(j,kk)-jjl(i,k)==0) then
                        ll=ll+1
                     end if
                  end do
               end do
               if(ll/=0) then
                  ks=imjl(i)+1
                  ke=ks+imjl(j)
                  jjl(i,ks:ke)=jjl(j,0:ke-ks)
                  call trimm(jjl(i,:),imjl(i),ke)
               end if
            end if
         end do
      end do
      do i=0,lp
         jltag(i)=minval(jjl(i,0:imjl(i)))
      end do
      njl(:)=-1
      jlcd(:,:)=-1
      do i=0,11
         nn=i/4+1
         ns=mod(nn,3)+1
         ne=mod(ns,3)+1
         ip=mod(i,4)/2
         jp=mod(mod(i,4),2)
         l=12*mb+i
         if( ijkp(nn)==ip*(nbpc(mb,nn)-1) .and. &
             ijkp(ns)==jp*(nbpc(mb,ns)-1) .and. &
             imjl(l)>=2) then
            njl(i)=l
            do k=1,imjl(l)
               mm=jjl(l,k)/12
               j=jjl(l,k)-12*mm
               np=j/4+1
               nq=mod(np,3)+1
               ii=mod(j,4)/2
               jj=mod(mod(j,4),2)
               jlcd(i,k)=idsd3(ii*(nbpc(mm,np)-1),jj*(nbpc(mm,nq)-1),ijkp(ne),mm,np)
            end do
         end if
      end do

   END SUBROUTINE search_line

   SUBROUTINE average_point(qak)
      real(kind=nr),    intent(inout), dimension(0:lmx,5) :: qak
      integer(kind=ni) :: kp, l, jp, ip, ii, i, j, lp, lq, itag
      real(kind=nr)    :: fctr
      real(kind=nr),dimension(:,:),allocatable :: rr

      allocate(rr(0:lmx,3))
      call p_null_req
      lp=-5
      lq=-5
      do i=0,7
         ii=njp(i)
         if(ii/=-1) then
            itag=jptag(ii)
            ip=mod(i,2)
            jp=mod(i,4)/2
            kp=i/4
            lp=lp+5
            l=indx3(ip*lxi,jp*let,kp*lze,1)
            rr(lp:lp+4,1)=qak(l,:)
            do j=1,imjp(ii)
               lq=lq+5
               call p_isend(rr(lp:lp+4,1), jpcd(i,j), itag)
               call p_irecv(rr(lq:lq+4,2), jpcd(i,j), itag)
            end do
         end if
      end do
      call p_waitall
      lp=-5
      lq=-5
      do i=0,7
         ii=njp(i)
         if(ii/=-1) then
            ip=mod(i,2)
            jp=mod(i,4)/2
            kp=i/4
            lp=lp+5
            do j=1,imjp(ii)
               lq=lq+5
               rr(lp:lp+4,1)=rr(lp:lp+4,1)+rr(lq:lq+4,2)
            end do
            fctr=one/(imjp(ii)+1)
            l=indx3(ip*lxi,jp*let,kp*lze,1)
            qak(l,:)=fctr*rr(lp:lp+4,1)
         end if
      end do
      deallocate(rr)

   END SUBROUTINE average_point


   SUBROUTINE average_line(qak)
      real(kind=nr),    intent(inout), dimension(0:lmx,5) :: qak
      integer(kind=ni) :: nn, ll, l, jp, ip, ii, i, j, k
      integer(kind=ni) :: js, je, ks, ke, lp, lq, itag
      real(kind=nr)    :: fctr
      real(kind=nr),dimension(:,:),allocatable :: rr

      allocate(rr(0:lmx,3))
      call p_null_req
      lp=0
      lq=0
      do i=0,11
         ii=njl(i)
         if(ii/=-1) then
            itag=jltag(ii)
            nn=i/4+1
            ip=mod(i,4)/2
            jp=mod(mod(i,4),2)
            ks=lp
            ke=ks+5*(ijk(3,nn)+1)-1
            lp=lp+ke-ks+1
            do k=0,ijk(3,nn)
               ll=ks+5*k
               l=indx3(ip*ijk(1,nn),jp*ijk(2,nn),k,nn)
               rr(ll:ll+4,1)=qak(l,:)
            end do
            do j=1,imjl(ii)
               js=lq
               je=js+ke-ks
               lq=lq+je-js+1
               call p_isend(rr(ks:ke,1), jlcd(i,j), itag)
               call p_irecv(rr(js:je,2), jlcd(i,j), itag)
            end do
         end if
      end do
      call p_waitall
      lp=0
      lq=0
      do i=0,11
         ii=njl(i)
         if(ii/=-1) then
            nn=i/4+1
            ip=mod(i,4)/2
            jp=mod(mod(i,4),2)
            ks=lp
            ke=ks+5*(ijk(3,nn)+1)-1
            lp=lp+ke-ks+1
            do j=1,imjl(ii)
               js=lq
               je=js+ke-ks
               lq=lq+je-js+1
               rr(ks:ke,1)=rr(ks:ke,1)+rr(js:je,2)
            end do
            fctr=one/(imjl(ii)+1)
            do k=0,ijk(3,nn)
               ll=ks+5*k
               l=indx3(ip*ijk(1,nn),jp*ijk(2,nn),k,nn)
               qak(l,:)=fctr*rr(ll:ll+4,1)
            end do
         end if
      end do
      deallocate(rr)

   END SUBROUTINE average_line

!===== SUBROUTINE FOR TRIMMING DOWN TO EFFECTIVE ARRAY

   subroutine trimm(itrim,lene,lens)

      integer(kind=ni),dimension(0:ljpl),intent(inout) :: itrim
      integer(kind=ni),intent(inout) :: lene,lens
      integer(kind=ni) :: k,kk

      do k=0,lens-1
         do kk=k+1,lens
            if(itrim(kk)-itrim(k)==0) then
               itrim(kk)=-1
            end if
         end do
      end do
      lene=lens
      do k=0,lens
         if(itrim(k)==-1) then
            lene=lene-1
            itrim(k)=itrim(k+1)
            itrim(k+1)=-1
         end if
      end do

   end subroutine trimm

!===== FUNCTION FOR PROCESSOR-CORE INDEX TRANSFORMATION WITHIN THE BLOCK IN 3D

   function idsd3(i,j,k,mm,nn) result(lm)

      integer(kind=ni),intent(in) :: i,j,k,mm,nn
      integer(kind=ni) :: lm

      select case(nn)
      case(1)
         lm=mo(mm)+(k*nbpc(mm,2)+j)*nbpc(mm,1)+i
      case(2)
         lm=mo(mm)+(j*nbpc(mm,2)+i)*nbpc(mm,1)+k
      case(3)
         lm=mo(mm)+(i*nbpc(mm,2)+k)*nbpc(mm,1)+j
      end select

   end function idsd3

END MODULE mo_domdcomp