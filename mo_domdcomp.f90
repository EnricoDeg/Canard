!*****
!***** DOMAIN DECOMPOSITION MODULE
!*****

MODULE mo_domdcomp
   use mainvar3d
   use mo_mpi
   use problemcase
   IMPLICIT NONE
   PUBLIC

   CONTAINS

   subroutine domdcomp_info

      ip=30*nthick+35*(1-nthick)
      jp=35*(1-nbody)+nbody*(20+5*nviscous)
      do mm=0,mbk
         select case(mm)
         case(0,3)
            nbbc(mm,1,:)=(/10,ip/)
            mbcd(mm,1,:)=(/-1,mm+1/)
         case(1,4)
            nbbc(mm,1,:)=(/ip,ip/)
            mbcd(mm,1,:)=(/mm-1,mm+1/)
         case(2,5)
            nbbc(mm,1,:)=(/ip,10/)
            mbcd(mm,1,:)=(/mm-1,-1/)
         end select

         select case(mm)
         case(0,2)
            nbbc(mm,2,:)=(/45,ip/)
            mbcd(mm,2,:)=(/mm+3,mm+3/)
         case(1)
            nbbc(mm,2,:)=(/45,ip/)
            mbcd(mm,2,:)=(/mm+3,mm+3/)
         case(3,5)
            nbbc(mm,2,:)=(/ip,45/)
            mbcd(mm,2,:)=(/mm-3,mm-3/)
         case(4)
            nbbc(mm,2,:)=(/ip,45/)
            mbcd(mm,2,:)=(/mm-3,mm-3/)
         end select

         nbbc(mm,3,:)=(/45,45/)
         mbcd(mm,3,:)=(/mm,mm/)
      end do

   end subroutine domdcomp_info

   SUBROUTINE domdcomp_init

      ijkp(1)=mod(myid-mo(mb),nbpc(mb,1))
      ijkp(2)=mod((myid-mo(mb))/nbpc(mb,1),nbpc(mb,2))
      ijkp(3)=mod((myid-mo(mb))/(nbpc(mb,1)*nbpc(mb,2)),nbpc(mb,3))

      do nn=1,3
         ns=mod(nn,3)+1
         ne=mod(ns,3)+1
         do ip=0,1
            mm=mbcd(mb,nn,ip)
            if(mm==-1) then
               mmcd(nn,ip)=-1
            else
               mmcd(nn,ip)=idsd3((1-ip)*(nbpc(mm,nn)-1),ijkp(ns),ijkp(ne),mm,nn) 
            end if
         end do
      end do

      do nn=1,3
         select case(nn)
         case (1)
            ll=lxio
            mp=1
         case (2)
            ll=leto
            mp=nbpc(mb,1)
         case (3)
            ll=lzeo
            mp=nbpc(mb,1)*nbpc(mb,2)
         end select
         lp=ijkp(nn)
         ma=nbpc(mb,nn)
         if(ma==1) then
            l=ll
            nbc(nn,0)=nbbc(mb,nn,0)
            nbc(nn,1)=nbbc(mb,nn,1)
            mcd(nn,0)=mmcd(nn,0)
            mcd(nn,1)=mmcd(nn,1)
         end if
         if(ma>=2) then
            if(lp==0) then
               l=ll-((ll+1)/ma)*(ma-1)
               nbc(nn,0)=nbbc(mb,nn,0)
               nbc(nn,1)=40
               mcd(nn,0)=mmcd(nn,0)
               mcd(nn,1)=myid+mp
            end if
            if(lp>0.and.lp<ma-1) then
               l=(ll+1)/ma-1
               nbc(nn,0)=40
               nbc(nn,1)=40
               mcd(nn,0)=myid-mp
               mcd(nn,1)=myid+mp
            end if
            if(lp==ma-1) then
               l=(ll+1)/ma-1
               nbc(nn,0)=40
               nbc(nn,1)=nbbc(mb,nn,1)
               mcd(nn,0)=myid-mp
               mcd(nn,1)=mmcd(nn,1)
            end if
         end if
         select case(nn)
         case (1)
            lxi=l
         case (2)
            let=l
         case (3)
            lze=l
         end select
      end do

   END SUBROUTINE domdcomp_init

   SUBROUTINE subdomain_init

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

   END SUBROUTINE subdomain_init

   SUBROUTINE search_point

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

   SUBROUTINE average_point

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
            rr(lp:lp+4,1)=qa(l,:)
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
            qa(l,:)=fctr*rr(lp:lp+4,1)
         end if
      end do

   END SUBROUTINE average_point


   SUBROUTINE average_line

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
               rr(ll:ll+4,1)=qa(l,:)
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
               qa(l,:)=fctr*rr(ll:ll+4,1)
            end do
         end if
      end do

   END SUBROUTINE average_line

   SUBROUTINE average_surface

      call p_null_req
      itag=30
      do nn=1,3
         select case(nn)
         case(1)
            drva=>drva1
            drvb=>drvb1
         case(2)
            drva=>drva2
            drvb=>drvb2
         case(3)
            drva=>drva3
            drvb=>drvb3
         end select
         do ip=0,1
            iq=1-ip
            np=nbc(nn,ip)
            i=ip*ijk(1,nn)
            if( (np-30) * (np-35) * (np-45) * &
                (abs(nextgcic-1)+abs((np-20)*(np-25))) == 0) then
               do k=0,ijk(3,nn)
                  kp=k*(ijk(2,nn)+1)
                  do j=0,ijk(2,nn)
                     l=indx3(i,j,k,nn)
                     jk=kp+j
                     drva(jk,:,ip)=qa(l,:)
                     rr(l,1)=one
                  end do
               end do
               call p_isend(drva(:,:,ip), mcd(nn,ip), itag+iq, 5*nbsize(nn))
               call p_irecv(drvb(:,:,ip), mcd(nn,ip), itag+ip, 5*nbsize(nn))
            end if
         end do
      end do
      call p_waitall
      do nn=1,3
         select case(nn)
         case(1)
            drvb=>drvb1
         case(2)
            drvb=>drvb2
         case(3)
            drvb=>drvb3
         end select
         do ip=0,1
            np=nbc(nn,ip)
            i=ip*ijk(1,nn)
            if( (np-30) * (np-35) * (np-45) * &
                (abs(nextgcic-1)+abs((np-20)*(np-25))) == 0) then
               do k=0,ijk(3,nn)
                  kp=k*(ijk(2,nn)+1)
                  do j=0,ijk(2,nn)
                     l=indx3(i,j,k,nn)
                     jk=kp+j
                     rr(l,1)=rr(l,1)+nrr(l)
                     rr(l,2)=nrr(l)/rr(l,1)
                     qa(l,:)=rr(l,2)*((rr(l,1)-one)*qa(l,:)+drvb(jk,:,ip))+(1-nrr(l))*qa(l,:)
                  end do
               end do
            end if
         end do
      end do

   END SUBROUTINE average_surface

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