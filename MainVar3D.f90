!*****
!***** MAIN VARIABLES & DATA FOR 3D NAVIER-STOKES/EULER SOLVER
!*****

module mainvar3d

   implicit none
   PUBLIC

!===== CONSTANT PARAMETERS

   integer,parameter :: icray=1

   integer,parameter :: int32=selected_int_kind(9),int64=selected_int_kind(18)
   integer,parameter :: ieee32=selected_real_kind(6,37),ieee64=selected_real_kind(15,307)
   integer,parameter :: ni=int32,nr=ieee64

   integer(kind=ni),parameter :: nrall=0,nrone=1,n45no=0,n45go=1
   integer(kind=ni),parameter :: lmd=11,lmf=11,lmp=max(lmd,lmf),mbci=4
   integer(kind=ni),parameter :: liofs=16,liofl=24,ljpl=100

   character(len=*),parameter :: fmts='es15.8',fmtl='es23.16',fmtsa=fmts//',a',fmtla=fmtl//',a'

   real(kind=nr),parameter :: zero=0,quarter=0.25_nr,half=0.5_nr,one=1,two=2,three=3,four=4,five=5
   real(kind=nr),parameter :: onethird=one/three,twothirds=two/three
   real(kind=nr),parameter :: pi=acos(-one),halfpi=pi/two,twopi=two*pi,sqrt2=sqrt(two),sqrt2i=one/sqrt2
   real(kind=nr),parameter :: sml=1.0e-8_nr,free=1.0e+6_nr
   real(kind=nr),parameter :: gam=1.4_nr,gamm1=gam-one,ham=one/gam,hamm1=one/gamm1,hamhamm1=ham*hamm1
   real(kind=nr),parameter :: prndtl=0.71_nr,gamm1prndtli=one/(gamm1*prndtl)

   real(kind=nr),parameter :: alpha=4/9.0_nr
   real(kind=nr),parameter :: beta=1/36.0_nr
   real(kind=nr),parameter :: aa=20/27.0_nr
   real(kind=nr),parameter :: ab=25/216.0_nr

   real(kind=nr),parameter :: alpha10=0.09166666665902128_nr
   real(kind=nr),parameter :: alpha12=1.3500000000438646_nr
   real(kind=nr),parameter :: beta13=0.2666666666794667_nr
   real(kind=nr),parameter :: a10=-0.3506944444280639_nr
   real(kind=nr),parameter :: a12=0.8249999999974599_nr
   real(kind=nr),parameter :: a13=0.7444444444773477_nr
   real(kind=nr),parameter :: a14=0.014583333334139971_nr

   real(kind=nr),parameter :: alpha01=8.000000000449464_nr
   real(kind=nr),parameter :: beta02=6.000000000542_nr
   real(kind=nr),parameter :: a01=-6.666666667794732_nr
   real(kind=nr),parameter :: a02=9.000000000052879_nr
   real(kind=nr),parameter :: a03=1.333333333178555_nr
   real(kind=nr),parameter :: a04=-0.0833333333455545_nr

!===== ALLOCATABLE MAIN ARRAYS

   integer(kind=ni),dimension(:,:,:),allocatable :: nbbc,mbcd
   integer(kind=ni),dimension(:,:),allocatable :: nbpc,lio,jjp,jjl
   integer(kind=ni),dimension(:),allocatable :: imjp,imjl,jptag,jltag
   integer(kind=ni),dimension(:),allocatable :: li,lcsz
   integer(kind=ni),dimension(:),allocatable :: lxim,letm,lzem,lpos
   integer(kind=ni),dimension(:),allocatable :: lximb,letmb,lzemb,mo,nrr,npex
   integer(kind=int64),dimension(:),allocatable :: lhmb

   real(kind=nr),dimension(:,:),allocatable :: qo,qa,qb,de
   real(kind=nr),dimension(:),allocatable :: txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz

   real(kind=nr),dimension(:,:),allocatable :: xim,etm,zem
   real(kind=nr),dimension(:),allocatable :: p,yaco
   real(kind=nr),dimension(:),allocatable :: sbcc

   real(kind=nr),dimension(:,:),allocatable :: rr,ss

   real(kind=nr),dimension(:,:),allocatable :: xu,yu
   real(kind=nr),dimension(:,:),allocatable :: xl,yl
   real(kind=nr),dimension(:),allocatable :: sa,sb

   real(kind=nr),dimension(:),allocatable :: asz,bsz
   real(kind=nr),dimension(:),allocatable :: times,vmpi

   real(kind=nr),dimension(:,:,:),pointer :: drva,drvb,send,recv,cm

   real(kind=nr),dimension(:,:,:),allocatable,target :: drva1,drva2,drva3
   real(kind=nr),dimension(:,:,:),allocatable,target :: drvb1,drvb2,drvb3
   real(kind=nr),dimension(:,:,:),allocatable,target :: send01,send02,send03
   real(kind=nr),dimension(:,:,:),allocatable,target :: send11,send12,send13
   real(kind=nr),dimension(:,:,:),allocatable,target :: recv01,recv02,recv03
   real(kind=nr),dimension(:,:,:),allocatable,target :: recv11,recv12,recv13
   real(kind=nr),dimension(:,:,:),allocatable,target :: cm1,cm2,cm3

   real(kind=ieee32),dimension(:,:),allocatable :: varm
   real(kind=ieee32),dimension(:),allocatable :: varr,vart,vara,varb,vmean,varmin,varmax

   character(13),dimension(:),allocatable :: ctecplt,cthead
   character(4),dimension(:),allocatable :: cfilet
   character(4),dimension(:),allocatable :: czonet

!===== CONSTANT-SIZED MAIN VARIABLES

   integer(kind=ni),dimension(0:3,3,0:1) :: ijl,ijp
   integer(kind=ni),dimension(3,0:1,0:1) :: ndf
   integer(kind=ni),dimension(0:11,ljpl) :: jlcd
   integer(kind=ni),dimension(0:7,ljpl) :: jpcd
   integer(kind=ni),dimension(3,3) :: ijk
   integer(kind=ni),dimension(3,0:1) :: nbc,mcd,mmcd,nsz
   integer(kind=ni),dimension(0:11) :: njl
   integer(kind=ni),dimension(0:7) :: njp
   integer(kind=ni),dimension(0:4) :: no
   integer(kind=ni),dimension(3) :: ms,me,nbsize,ijkp,nnf
   integer(kind=ni) :: lxio,leto,lzeo,lxi,let,lze,lmx,lim,lsz,nrecs,nrecd
   integer(kind=ni) :: i,ii,is,ie,ip,iq,j,jj,js,je,jp,jq,jk,k,kk,ks,ke,kp,l,lh,ll,lp,lq,ltomb
   integer(kind=ni) :: m,ma,mb,mm,mp,mq,mbk,mps,mpe,n,ndt,nn,nk,ns,ne,np,nq,nt,nz,ndati,nsigi,nout,nfile
   integer(kind=ni) :: nts,nscrn,nsgnl,ndata,ndatafl,ndataav,nviscous,nkrk,nsmf,nrestart,nextrabc,nextgcic
   integer(kind=int64) :: nlmx,llmb,llmo,lis,lie,ljs,lje

   real(kind=nr),dimension(0:lmp,0:1,0:1) :: pbci,pbco
   real(kind=nr),dimension(-2:2,0:2,0:1) :: albef
   real(kind=nr),dimension(mbci,mbci) :: cbca,cbcs
   real(kind=nr),dimension(5,5) :: xt
   real(kind=nr),dimension(0:4,0:2) :: fbc
   real(kind=nr),dimension(0:1,0:1) :: pbcot
   real(kind=nr),dimension(0:lmp) :: sap
   real(kind=nr),dimension(mbci) :: rbci,sbci
   real(kind=nr),dimension(5) :: cha,dha
   real(kind=nr),dimension(3) :: ve,dm,rv,uoo,umf,dudtmf
   real(kind=nr) :: alphf,betf,fa,fb,fc
   real(kind=nr) :: ra0,ra1,ra2,ra3,res,fctr,dfdt
   real(kind=nr) :: reoo,tempoo,amach1,amach2,amach3,wtemp,cfl,tmax,timf,fltk,fltrbc,dto
   real(kind=nr) :: aoo,amachoo,srefoo,srefp1dre
   real(kind=nr) :: dt,dts,dte,dtk,dtko,dtsum,dtwi,timo,tsam,wts,wte,wtime
   real(kind=nr) :: vn,vs,hv2,ao,bo,co,ho,aoi,rhoi,progmf,sqrtrema,sqrtremai

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
      allocate(p(0:lmx),yaco(0:lmx),varr(0:lmx),nrr(0:lmx),npex(0:lmx))

      if(nviscous==1) then
         allocate(txx(0:lmx),tyy(0:lmx),tzz(0:lmx))
         allocate(txy(0:lmx),tyz(0:lmx),tzx(0:lmx))
         allocate(hxx(0:lmx),hyy(0:lmx),hzz(0:lmx))
      end if

      ii=nbsize(1)-1
      jj=nbsize(2)-1
      kk=nbsize(3)-1
      allocate(drva1(0:ii,5,0:1),drva2(0:jj,5,0:1),drva3(0:kk,5,0:1))
      allocate(drvb1(0:ii,5,0:1),drvb2(0:jj,5,0:1),drvb3(0:kk,5,0:1))
      allocate(send01(0:ii,0:1,0:1),send02(0:jj,0:1,0:1),send03(0:kk,0:1,0:1))
      allocate(recv01(0:ii,0:1,0:1),recv02(0:jj,0:1,0:1),recv03(0:kk,0:1,0:1))
      allocate(send11(0:ii,0:2,0:1),send12(0:jj,0:2,0:1),send13(0:kk,0:2,0:1))
      allocate(recv11(0:ii,0:2,0:1),recv12(0:jj,0:2,0:1),recv13(0:kk,0:2,0:1))
      allocate(cm1(0:ii,3,0:1),cm2(0:jj,3,0:1),cm3(0:kk,3,0:1))
      allocate(xu(0:lim,3),yu(0:lim,3),xl(0:lim,2),yl(0:lim,2),li(0:lim),sa(0:lim),sb(0:lim))

   END SUBROUTINE allocate_memory

end module mainvar3d

!*****