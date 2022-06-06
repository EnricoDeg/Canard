!*****
!***** MAIN VARIABLES & DATA FOR 3D NAVIER-STOKES/EULER SOLVER
!*****

module mo_vars
   use mo_kind, ONLY : ni, nr, ieee32, int64
   implicit none
   PUBLIC

!===== ALLOCATABLE MAIN ARRAYS
   integer(kind=ni),dimension(:,:),allocatable :: lio
   integer(kind=ni),dimension(:),allocatable :: lpos

   real(kind=nr),dimension(:,:),allocatable :: qo,qa,qb,de
   real(kind=nr),dimension(:),allocatable :: txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz

   real(kind=nr),dimension(:,:),allocatable :: xim,etm,zem
   real(kind=nr),dimension(:),allocatable :: p,yaco

   real(kind=nr),dimension(:,:),allocatable :: rr,ss

   real(kind=nr),dimension(:),allocatable :: times

   real(kind=nr),dimension(:,:,:),pointer :: drva,cm

   real(kind=nr),dimension(:,:,:),allocatable,target :: drva1,drva2,drva3
   real(kind=nr),dimension(:,:,:),allocatable,target :: cm1,cm2,cm3

   real(kind=ieee32),dimension(:),allocatable :: varr,vart,vmean

!===== CONSTANT-SIZED MAIN VARIABLES
   integer(kind=ni),dimension(3) :: nnf
   integer(kind=ni) :: lim,nrecs,nrecd
   integer(kind=ni) :: mbk, n, ndt, nk, mq
   integer(kind=ni) :: nts,nscrn,nsgnl,ndata,ndatafl,ndataav,nkrk,nsmf,nrestart,nextrabc,nextgcic

   real(kind=nr),dimension(5,5) :: xt
   real(kind=nr),dimension(5) :: cha,dha
   real(kind=nr),dimension(3) :: umf,dudtmf
   real(kind=nr) :: wtemp,cfl,tmax,timf,dto
   real(kind=nr) :: aoo,srefoo,srefp1dre
   real(kind=nr) :: dt,dts,dte,timo,tsam
   real(kind=nr) :: hv2,ao,bo,aoi,sqrtrema,sqrtremai

   real(kind=nr) :: szco
   integer(kind=ni) :: nbody

   character(5) :: cnnode
   character(16) :: cinput
   character(16) :: cgrid
   character(18) :: cdata

   CONTAINS

   SUBROUTINE allocate_memory(lmx, nbsize)
      integer(kind=ni), intent(IN)               :: lmx
      integer(kind=ni), dimension(3), intent(IN) :: nbsize
      integer(kind=ni) :: ii, jj, kk

      allocate(qo(0:lmx,5),qa(0:lmx,5),qb(0:lmx,5),de(0:lmx,5))
      allocate(xim(0:lmx,3),etm(0:lmx,3),zem(0:lmx,3),rr(0:lmx,3),ss(0:lmx,3))
      allocate(p(0:lmx),yaco(0:lmx),varr(0:lmx))

#ifdef VISCOUS
         allocate(txx(0:lmx), tyy(0:lmx), tzz(0:lmx))
         allocate(txy(0:lmx), tyz(0:lmx), tzx(0:lmx))
         allocate(hxx(0:lmx), hyy(0:lmx), hzz(0:lmx))
#endif

      ii=nbsize(1)-1
      jj=nbsize(2)-1
      kk=nbsize(3)-1
      allocate(drva1(0:ii,5,0:1), drva2(0:jj,5,0:1), drva3(0:kk,5,0:1))
      allocate(cm1(0:ii,3,0:1)  , cm2(0:jj,3,0:1)  , cm3(0:kk,3,0:1))

   END SUBROUTINE allocate_memory

end module mo_vars
