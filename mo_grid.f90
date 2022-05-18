!*****
!***** GRID MODULE
!*****

MODULE mo_grid
   use mo_mpi, ONLY : myid, p_barrier
   use mo_vars
   use mo_numerics
   use mo_utils
   use mo_gridgen
   
   IMPLICIT NONE
   PUBLIC
   CONTAINS

   SUBROUTINE calc_grid

      allocate(lio(0:let,0:lze))
      do k=0,lze
         kp=k*(leto+1)*(lxio+1)
         do j=0,let
            jp=j*(lxio+1)
            lio(j,k)=jp+kp
         end do
      end do
      call makegrid
      call p_barrier

      open(9,file=cgrid,access='direct',form='unformatted',recl=3*nrecd,status='old')
      lp=lpos(myid)
      do k=0,lze
         do j=0,let
            lq=lp+lio(j,k)
            do i=0,lxi
               l=indx3(i,j,k,1)
               read(9,rec=lq+i+1) ss(l,:)
            end do
         end do
      end do
      close(9)
      call p_barrier
      if(myid==mo(mb)) then
         open(9,file=cgrid,status='old')
         close(9,status='delete')
      end if

   END SUBROUTINE calc_grid

   SUBROUTINE calc_grid_metrics

      real(kind=nr),dimension(:,:),allocatable :: dek, qok, qak, rrk
      allocate(dek(0:lmx,5), qok(0:lmx,5), qak(0:lmx,5), rrk(0:lmx,3))

      rrk(:,1)=ss(:,1)
      m=1
      call mpigo(rrk, ijk, nbc, mcd, nbsize,0,nrone,n45go,m)
      call deriv(rrk, lxi, let, lze, ijk, 3, 1, m)
      call deriv(rrk, lxi, let, lze, ijk, 2, 1, m)
      call deriv(rrk, lxi, let, lze, ijk, 1, 1, m)
      qok(:,1)=rrk(:,1)
      qok(:,2)=rrk(:,2)
      qok(:,3)=rrk(:,3)

      rrk(:,1)=ss(:,2)
      m=2
      call mpigo(rrk, ijk, nbc, mcd, nbsize,0,nrone,n45go,m)
      call deriv(rrk, lxi, let, lze, ijk, 3, 1, m)
      call deriv(rrk, lxi, let, lze, ijk, 2, 1, m)
      call deriv(rrk, lxi, let, lze, ijk, 1, 1, m)
      qak(:,1)=rrk(:,1)
      qak(:,2)=rrk(:,2)
      qak(:,3)=rrk(:,3)

      rrk(:,1)=ss(:,3)
      m=3
      call mpigo(rrk, ijk, nbc, mcd, nbsize,0,nrone,n45go,m)
      call deriv(rrk, lxi, let, lze, ijk, 3, 1, m)
      call deriv(rrk, lxi, let, lze, ijk, 2, 1, m)
      call deriv(rrk, lxi, let, lze, ijk, 1, 1, m)
      dek(:,1)=rrk(:,1)
      dek(:,2)=rrk(:,2)
      dek(:,3)=rrk(:,3)

      xim(:,1)=qak(:,2)*dek(:,3)-dek(:,2)*qak(:,3)
      xim(:,2)=dek(:,2)*qok(:,3)-qok(:,2)*dek(:,3)
      xim(:,3)=qok(:,2)*qak(:,3)-qak(:,2)*qok(:,3)
      etm(:,1)=qak(:,3)*dek(:,1)-dek(:,3)*qak(:,1)
      etm(:,2)=dek(:,3)*qok(:,1)-qok(:,3)*dek(:,1)
      etm(:,3)=qok(:,3)*qak(:,1)-qak(:,3)*qok(:,1)
      zem(:,1)=qak(:,1)*dek(:,2)-dek(:,1)*qak(:,2)
      zem(:,2)=dek(:,1)*qok(:,2)-qok(:,1)*dek(:,2)
      zem(:,3)=qok(:,1)*qak(:,2)-qak(:,1)*qok(:,2)

      do m=1,3
         call mpigo(xim(:,m), ijk, nbc, mcd, nbsize,1,n45no,9*(m-1)+1)
         call filte(xim(:,m), lxi, let, lze, ijk, nnf(1))
         call mpigo(xim(:,m), ijk, nbc, mcd, nbsize,1,n45no,9*(m-1)+2)
         call filte(xim(:,m), lxi, let, lze, ijk, nnf(2))
         call mpigo(xim(:,m), ijk, nbc, mcd, nbsize,1,n45no,9*(m-1)+3)
         call filte(xim(:,m), lxi, let, lze, ijk, nnf(3))

         call mpigo(etm(:,m), ijk, nbc, mcd, nbsize,1,n45no,9*(m-1)+4)
         call filte(etm(:,m), lxi, let, lze, ijk, nnf(1))
         call mpigo(etm(:,m), ijk, nbc, mcd, nbsize,1,n45no,9*(m-1)+5)
         call filte(etm(:,m), lxi, let, lze, ijk, nnf(2))
         call mpigo(etm(:,m), ijk, nbc, mcd, nbsize,1,n45no,9*(m-1)+6)
         call filte(etm(:,m), lxi, let, lze, ijk, nnf(3))

         call mpigo(zem(:,m), ijk, nbc, mcd, nbsize,1,n45no,9*(m-1)+7)
         call filte(zem(:,m), lxi, let, lze, ijk, nnf(1))
         call mpigo(zem(:,m), ijk, nbc, mcd, nbsize,1,n45no,9*(m-1)+8)
         call filte(zem(:,m), lxi, let, lze, ijk, nnf(2))
         call mpigo(zem(:,m), ijk, nbc, mcd, nbsize,1,n45no,9*(m-1)+9)
         call filte(zem(:,m), lxi, let, lze, ijk, nnf(3))
      end do
      yaco(:)=three/(qok(:,1)*xim(:,1)+qok(:,2)*etm(:,1)+qok(:,3)*zem(:,1)&
                    +qak(:,1)*xim(:,2)+qak(:,2)*etm(:,2)+qak(:,3)*zem(:,2)&
                    +dek(:,1)*xim(:,3)+dek(:,2)*etm(:,3)+dek(:,3)*zem(:,3))

      do nn=1,3
         do ip=0,1
            i=ip*ijk(1,nn)
            do k=0,ijk(3,nn)
               kp=k*(ijk(2,nn)+1)
               do j=0,ijk(2,nn)
                  jk=kp+j
                  l=indx3(i,j,k,nn)
                  select case(nn)
                  case(1)
                     rv(:)=yaco(l)*xim(l,:)
                     fctr=one/sqrt(rv(1)*rv(1)+rv(2)*rv(2)+rv(3)*rv(3))
                     cm1(jk,:,ip)=fctr*rv(:)
                  case(2)
                     rv(:)=yaco(l)*etm(l,:)
                     fctr=one/sqrt(rv(1)*rv(1)+rv(2)*rv(2)+rv(3)*rv(3))
                     cm2(jk,:,ip)=fctr*rv(:)
                  case(3)
                     rv(:)=yaco(l)*zem(l,:)
                     fctr=one/sqrt(rv(1)*rv(1)+rv(2)*rv(2)+rv(3)*rv(3))
                     cm3(jk,:,ip)=fctr*rv(:)
                  end select
               end do
            end do
         end do
      end do
      deallocate(dek, qok, qak, rrk)

   END SUBROUTINE calc_grid_metrics

END MODULE mo_grid