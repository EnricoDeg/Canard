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

      rr(:,1)=ss(:,1)
      m=1
      call mpigo(rr, ijk, nbc, mcd, nbsize,0,nrone,n45go,m)
      call deriv(rr, lxi, let, lze, ijk, 3, 1, m)
      call deriv(rr, lxi, let, lze, ijk, 2, 1, m)
      call deriv(rr, lxi, let, lze, ijk, 1, 1, m)
      qo(:,1)=rr(:,1)
      qo(:,2)=rr(:,2)
      qo(:,3)=rr(:,3)

      rr(:,1)=ss(:,2)
      m=2
      call mpigo(rr, ijk, nbc, mcd, nbsize,0,nrone,n45go,m)
      call deriv(rr, lxi, let, lze, ijk, 3, 1, m)
      call deriv(rr, lxi, let, lze, ijk, 2, 1, m)
      call deriv(rr, lxi, let, lze, ijk, 1, 1, m)
      qa(:,1)=rr(:,1)
      qa(:,2)=rr(:,2)
      qa(:,3)=rr(:,3)

      rr(:,1)=ss(:,3)
      m=3
      call mpigo(rr, ijk, nbc, mcd, nbsize,0,nrone,n45go,m)
      call deriv(rr, lxi, let, lze, ijk, 3, 1, m)
      call deriv(rr, lxi, let, lze, ijk, 2, 1, m)
      call deriv(rr, lxi, let, lze, ijk, 1, 1, m)
      de(:,1)=rr(:,1)
      de(:,2)=rr(:,2)
      de(:,3)=rr(:,3)

      xim(:,1)=qa(:,2)*de(:,3)-de(:,2)*qa(:,3)
      xim(:,2)=de(:,2)*qo(:,3)-qo(:,2)*de(:,3)
      xim(:,3)=qo(:,2)*qa(:,3)-qa(:,2)*qo(:,3)
      etm(:,1)=qa(:,3)*de(:,1)-de(:,3)*qa(:,1)
      etm(:,2)=de(:,3)*qo(:,1)-qo(:,3)*de(:,1)
      etm(:,3)=qo(:,3)*qa(:,1)-qa(:,3)*qo(:,1)
      zem(:,1)=qa(:,1)*de(:,2)-de(:,1)*qa(:,2)
      zem(:,2)=de(:,1)*qo(:,2)-qo(:,1)*de(:,2)
      zem(:,3)=qo(:,1)*qa(:,2)-qa(:,1)*qo(:,2)

      do m=1,3
         rr(:,1)=xim(:,m)
         call mpigo(rr, ijk, nbc, mcd, nbsize,1,nrone,n45no,9*(m-1)+1)
         call filte(rr(:,1), lxi, let, lze, ijk, nnf(1),1)
         call mpigo(rr, ijk, nbc, mcd, nbsize,1,nrone,n45no,9*(m-1)+2)
         call filte(rr(:,1), lxi, let, lze, ijk, nnf(2),1)
         call mpigo(rr, ijk, nbc, mcd, nbsize,1,nrone,n45no,9*(m-1)+3)
         call filte(rr(:,1), lxi, let, lze, ijk, nnf(3),1)
         xim(:,m)=rr(:,1)
         rr(:,1)=etm(:,m)
         call mpigo(rr, ijk, nbc, mcd, nbsize,1,nrone,n45no,9*(m-1)+4)
         call filte(rr(:,1), lxi, let, lze, ijk, nnf(1),1)
         call mpigo(rr, ijk, nbc, mcd, nbsize,1,nrone,n45no,9*(m-1)+5)
         call filte(rr(:,1), lxi, let, lze, ijk, nnf(2),1)
         call mpigo(rr, ijk, nbc, mcd, nbsize,1,nrone,n45no,9*(m-1)+6)
         call filte(rr(:,1), lxi, let, lze, ijk, nnf(3),1)
         etm(:,m)=rr(:,1)
         rr(:,1)=zem(:,m)
         call mpigo(rr, ijk, nbc, mcd, nbsize,1,nrone,n45no,9*(m-1)+7)
         call filte(rr(:,1), lxi, let, lze, ijk, nnf(1),1)
         call mpigo(rr, ijk, nbc, mcd, nbsize,1,nrone,n45no,9*(m-1)+8)
         call filte(rr(:,1), lxi, let, lze, ijk, nnf(2),1)
         call mpigo(rr, ijk, nbc, mcd, nbsize,1,nrone,n45no,9*(m-1)+9)
         call filte(rr(:,1), lxi, let, lze, ijk, nnf(3),1)
         zem(:,m)=rr(:,1)
      end do
      yaco(:)=three/(qo(:,1)*xim(:,1)+qo(:,2)*etm(:,1)+qo(:,3)*zem(:,1)&
                    +qa(:,1)*xim(:,2)+qa(:,2)*etm(:,2)+qa(:,3)*zem(:,2)&
                    +de(:,1)*xim(:,3)+de(:,2)*etm(:,3)+de(:,3)*zem(:,3))

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

   END SUBROUTINE calc_grid_metrics

END MODULE mo_grid