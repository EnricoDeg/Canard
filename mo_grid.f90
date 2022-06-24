!*****
!***** GRID MODULE
!*****

MODULE mo_grid
   use mo_mpi,        ONLY : myid, p_barrier
   use mo_kind,       ONLY : nr, ni
   use mo_parameters, ONLY : n45go, n45no, nrone, one, three
   use mo_vars,       ONLY : cgrid, mbk, nrecd, nnf, lpos
   use mo_domdcomp,   ONLY : t_domdcomp
   use mo_numerics,   ONLY : mpigo, deriv, filte
   use mo_utils,      ONLY : indx3
   use mo_gridgen,    ONLY : makegrid
   
   IMPLICIT NONE
   PUBLIC

   integer(kind=ni), public, dimension(:,:),   allocatable          :: lio
   real(kind=nr),    public, dimension(:),     allocatable          :: yaco
   real(kind=nr),    public, dimension(:,:),   allocatable          :: xim, etm, zem
   real(kind=nr),    public, dimension(:,:,:), allocatable, target  :: cm1, cm2, cm3

   CONTAINS

   SUBROUTINE allocate_grid(p_domdcomp)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      integer(kind=ni) :: ii, jj, kk

      ii = p_domdcomp%nbsize(1)-1
      jj = p_domdcomp%nbsize(2)-1
      kk = p_domdcomp%nbsize(3)-1
      allocate(lio(0:p_domdcomp%let,0:p_domdcomp%lze), yaco(0:p_domdcomp%lmx))
      allocate(cm1(0:ii,3,0:1), cm2(0:jj,3,0:1), cm3(0:kk,3,0:1))
      allocate(xim(0:p_domdcomp%lmx,3), etm(0:p_domdcomp%lmx,3), zem(0:p_domdcomp%lmx,3))

   END SUBROUTINE allocate_grid

   SUBROUTINE calc_grid(p_domdcomp, ssk)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      real(kind=nr), intent(out) :: ssk(0:p_domdcomp%lmx,3)
      integer(kind=ni) :: i, j, k, l, jk, kp, lq, jp, lp

      do k=0,p_domdcomp%lze
         kp = k * ( p_domdcomp%leto + 1 ) * ( p_domdcomp%lxio + 1 )
         do j=0,p_domdcomp%let
            jp = j * ( p_domdcomp%lxio + 1 )
            lio(j,k) = jp + kp
         end do
      end do
      call makegrid(p_domdcomp%mb, p_domdcomp%lxio, p_domdcomp%leto, p_domdcomp%mo, mbk)
      call p_barrier

      open(9,file=cgrid,access='direct',form='unformatted',recl=3*nrecd,status='old')
      lp = lpos(myid)
      do k=0,p_domdcomp%lze
         do j=0,p_domdcomp%let
            lq = lp + lio(j,k)
            do i=0,p_domdcomp%lxi
               l = indx3(i, j, k, 1, p_domdcomp%lxi, p_domdcomp%let)
               read(9,rec=lq+i+1) ssk(l,:)
            end do
         end do
      end do
      close(9)
      call p_barrier
      if ( myid == p_domdcomp%mo(p_domdcomp%mb) ) then
         open(9,file=cgrid,status='old')
         close(9,status='delete')
      end if

   END SUBROUTINE calc_grid

   SUBROUTINE calc_grid_metrics(p_domdcomp, ssk)
      type(t_domdcomp), intent(IN) :: p_domdcomp
      real(kind=nr), intent(in) :: ssk(0:p_domdcomp%lmx,3)
      integer(kind=ni) :: m, nn, ip, i, j, k, l, jk, kp
      real(kind=nr)    :: fctr
      real(kind=nr), dimension(:,:), allocatable :: dek, qok, qak, rrk
      real(kind=nr), dimension(3) :: rv

      allocate(dek(0:p_domdcomp%lmx,5), qok(0:p_domdcomp%lmx,5), &
               qak(0:p_domdcomp%lmx,5), rrk(0:p_domdcomp%lmx,3))

      rrk(:,1)=ssk(:,1)
      m=1
      call mpigo(rrk, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                      p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrone, n45go, m, p_domdcomp%lxi, p_domdcomp%let)
      call deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 3, 1, m)
      call deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 2, 1, m)
      call deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 1, 1, m)
      qok(:,1)=rrk(:,1)
      qok(:,2)=rrk(:,2)
      qok(:,3)=rrk(:,3)

      rrk(:,1)=ssk(:,2)
      m=2
      call mpigo(rrk, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                      p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrone, n45go, m, p_domdcomp%lxi, p_domdcomp%let)
      call deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 3, 1, m)
      call deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 2, 1, m)
      call deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 1, 1, m)
      qak(:,1)=rrk(:,1)
      qak(:,2)=rrk(:,2)
      qak(:,3)=rrk(:,3)

      rrk(:,1)=ssk(:,3)
      m=3
      call mpigo(rrk, p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                      p_domdcomp%mcd, p_domdcomp%nbsize, 0, nrone, n45go, m, p_domdcomp%lxi, p_domdcomp%let)
      call deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 3, 1, m)
      call deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 2, 1, m)
      call deriv(rrk, p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                      p_domdcomp%lze, p_domdcomp%ijk, 1, 1, m)
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
         call mpigo(xim(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+1, p_domdcomp%lxi, p_domdcomp%let)
         call filte(xim(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, nnf(1))
         call mpigo(xim(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+2, p_domdcomp%lxi, p_domdcomp%let)
         call filte(xim(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, nnf(2))
         call mpigo(xim(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+3, p_domdcomp%lxi, p_domdcomp%let)
         call filte(xim(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, nnf(3))

         call mpigo(etm(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+4, p_domdcomp%lxi, p_domdcomp%let)
         call filte(etm(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, nnf(1))
         call mpigo(etm(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+5, p_domdcomp%lxi, p_domdcomp%let)
         call filte(etm(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, nnf(2))
         call mpigo(etm(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+6, p_domdcomp%lxi, p_domdcomp%let)
         call filte(etm(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, nnf(3))

         call mpigo(zem(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+7, p_domdcomp%lxi, p_domdcomp%let)
         call filte(zem(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, & 
                              p_domdcomp%lze, p_domdcomp%ijk, nnf(1))
         call mpigo(zem(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+8, p_domdcomp%lxi, p_domdcomp%let)
         call filte(zem(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, nnf(2))
         call mpigo(zem(:,m), p_domdcomp%lmx, p_domdcomp%ijk, p_domdcomp%nbc, &
                              p_domdcomp%mcd, p_domdcomp%nbsize, 1, n45no, 9*(m-1)+9, p_domdcomp%lxi, p_domdcomp%let)
         call filte(zem(:,m), p_domdcomp%lmx, p_domdcomp%lxi, p_domdcomp%let, &
                              p_domdcomp%lze, p_domdcomp%ijk, nnf(3))
      end do
      yaco(:)=three/(qok(:,1)*xim(:,1)+qok(:,2)*etm(:,1)+qok(:,3)*zem(:,1)&
                    +qak(:,1)*xim(:,2)+qak(:,2)*etm(:,2)+qak(:,3)*zem(:,2)&
                    +dek(:,1)*xim(:,3)+dek(:,2)*etm(:,3)+dek(:,3)*zem(:,3))

      do nn=1,3
         do ip=0,1
            i = ip * p_domdcomp%ijk(1,nn)
            do k=0,p_domdcomp%ijk(3,nn)
               kp = k * ( p_domdcomp%ijk(2,nn) + 1 )
               do j=0,p_domdcomp%ijk(2,nn)
                  jk = kp + j
                  l = indx3(i, j, k, nn, p_domdcomp%lxi, p_domdcomp%let)
                  select case(nn)
                  case(1)
                     rv(:) = yaco(l) * xim(l,:)
                     fctr  = one / sqrt(rv(1) * rv(1) + rv(2) * rv(2) + rv(3) * rv(3))
                     cm1(jk,:,ip) = fctr * rv(:)
                  case(2)
                     rv(:) = yaco(l) * etm(l,:)
                     fctr  = one / sqrt(rv(1) * rv(1) + rv(2) * rv(2) + rv(3) * rv(3))
                     cm2(jk,:,ip) = fctr * rv(:)
                  case(3)
                     rv(:) = yaco(l) * zem(l,:)
                     fctr  = one / sqrt(rv(1) * rv(1) + rv(2) * rv(2) + rv(3) * rv(3))
                     cm3(jk,:,ip) = fctr * rv(:)
                  end select
               end do
            end do
         end do
      end do
      deallocate(dek, qok, qak, rrk)

   END SUBROUTINE calc_grid_metrics

END MODULE mo_grid