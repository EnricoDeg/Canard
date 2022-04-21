!*****
!***** MAIN VARIABLES & DATA FOR 3D NAVIER-STOKES/EULER SOLVER
!*****

module mainvar3d
   use mo_kind
   use mo_parameters
   implicit none
   PUBLIC

!===== ALLOCATABLE MAIN ARRAYS

   integer(kind=ni),dimension(:,:),allocatable :: nbpc,lio
   integer(kind=ni),dimension(:),allocatable :: lcsz
   integer(kind=ni),dimension(:),allocatable :: lxim,letm,lzem,lpos
   integer(kind=ni),dimension(:),allocatable :: lximb,letmb,lzemb,mo

   real(kind=nr),dimension(:,:),allocatable :: qo,qa,qb,de
   real(kind=nr),dimension(:),allocatable :: txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz

   real(kind=nr),dimension(:,:),allocatable :: xim,etm,zem
   real(kind=nr),dimension(:),allocatable :: p,yaco

   real(kind=nr),dimension(:,:),allocatable :: rr,ss

   real(kind=nr),dimension(:),allocatable :: asz,bsz
   real(kind=nr),dimension(:),allocatable :: times

   real(kind=nr),dimension(:,:,:),pointer :: drva,cm

   real(kind=nr),dimension(:,:,:),allocatable,target :: drva1,drva2,drva3
   real(kind=nr),dimension(:,:,:),allocatable,target :: cm1,cm2,cm3

   real(kind=ieee32),dimension(:,:),allocatable :: varm
   real(kind=ieee32),dimension(:),allocatable :: varr,vart,vara,varb,vmean,varmin,varmax

   character(13),dimension(:),allocatable :: ctecplt,cthead
   character(4),dimension(:),allocatable :: cfilet
   character(4),dimension(:),allocatable :: czonet

!===== CONSTANT-SIZED MAIN VARIABLES

   integer(kind=ni),dimension(3,3) :: ijk
   integer(kind=ni),dimension(3,0:1) :: nbc,mcd,nsz
   integer(kind=ni),dimension(0:11) :: njl
   integer(kind=ni),dimension(0:7) :: njp
   integer(kind=ni),dimension(0:4) :: no
   integer(kind=ni),dimension(3) :: ms,me,nbsize,nnf
   integer(kind=ni) :: lxio,leto,lzeo,lxi,let,lze,lmx,lim,lsz,nrecs,nrecd
   integer(kind=ni) :: i,ii,is,ie,ip,iq,j,jj,js,je,jp,jq,jk,k,kk,ks,ke,kp,l,lh,ll,lp,lq,ltomb
   integer(kind=ni) :: m,ma,mb,mm,mp,mq,mbk,mps,mpe,n,ndt,nn,nk,ns,ne,np,nq,nt,nz,ndati,nsigi,nout,nfile
   integer(kind=ni) :: nts,nscrn,nsgnl,ndata,ndatafl,ndataav,nviscous,nkrk,nsmf,nrestart,nextrabc,nextgcic
   integer(kind=ni) :: itag, lmpi
   integer(kind=int64) :: nlmx,llmb,llmo,lis,lie,ljs,lje

   real(kind=nr),dimension(5,5) :: xt
   real(kind=nr),dimension(5) :: cha,dha
   real(kind=nr),dimension(3) :: dm,rv,uoo,umf,dudtmf
   real(kind=nr) :: ra0,ra1,ra2,ra3,res,fctr,dfdt
   real(kind=nr) :: reoo,tempoo,amach1,amach2,amach3,wtemp,cfl,tmax,timf,dto
   real(kind=nr) :: aoo,amachoo,srefoo,srefp1dre
   real(kind=nr) :: dt,dts,dte,dtk,dtko,dtsum,dtwi,timo,tsam,wts,wte,wtime
   real(kind=nr) :: hv2,ao,bo,co,ho,aoi,rhoi,progmf,sqrtrema,sqrtremai

   character(1),dimension(0:4) :: cno
   character(3) :: cndata
   character(5) :: cnnode
   character(13) :: coutput
   character(16) :: cinput
   character(16) :: cgrid
   character(18) :: cdata
   character(19) :: crestart

   CONTAINS

   SUBROUTINE allocate_memory

      allocate(qo(0:lmx,5),qa(0:lmx,5),qb(0:lmx,5),de(0:lmx,5))
      allocate(xim(0:lmx,3),etm(0:lmx,3),zem(0:lmx,3),rr(0:lmx,3),ss(0:lmx,3))
      allocate(p(0:lmx),yaco(0:lmx),varr(0:lmx))

      if(nviscous==1) then
         allocate(txx(0:lmx), tyy(0:lmx), tzz(0:lmx))
         allocate(txy(0:lmx), tyz(0:lmx), tzx(0:lmx))
         allocate(hxx(0:lmx), hyy(0:lmx), hzz(0:lmx))
      end if

      ii=nbsize(1)-1
      jj=nbsize(2)-1
      kk=nbsize(3)-1
      allocate(drva1(0:ii,5,0:1), drva2(0:jj,5,0:1), drva3(0:kk,5,0:1))
      allocate(cm1(0:ii,3,0:1)  , cm2(0:jj,3,0:1)  , cm3(0:kk,3,0:1))

   END SUBROUTINE allocate_memory

end module mainvar3d

!*****