!*****
!***** MAIN VARIABLES & DATA FOR 3D NAVIER-STOKES/EULER SOLVER
!*****

 module mainvar3d

 implicit none

!===== CONSTANT PARAMETERS

 integer,parameter :: icray=1

 integer,parameter :: int32=selected_int_kind(9),int64=selected_int_kind(18)
 integer,parameter :: ieee32=selected_real_kind(6,37),ieee64=selected_real_kind(15,307)
 integer,parameter :: ni=int32,nr=ieee64

 integer(kind=ni),parameter :: nrall=0,nrone=1,n45no=0,n45go=1
 integer(kind=ni),parameter :: lmd=11,lmf=11,lmp=max(lmd,lmf),mbci=4
 integer(kind=ni),parameter :: liofs=16,liofl=24,ljpl=100,nfilters=1000

 character(len=*),parameter :: fmts='es15.8',fmtl='es23.16',fmtsa=fmts//',a',fmtla=fmtl//',a'

 real(kind=nr),parameter :: zero=0,quarter=0.25_nr,half=0.5_nr,one=1,two=2,three=3,four=4,five=5
 real(kind=nr),parameter :: onethird=one/three,twothirds=two/three
 real(kind=nr),parameter :: pi=acos(-one),halfpi=pi/two,twopi=two*pi,sqrt2=sqrt(two),sqrt2i=one/sqrt2
 real(kind=nr),parameter :: sml=1.0e-8_nr,free=65536.0_nr
 real(kind=nr),parameter :: gam=1.4_nr,gamm1=gam-one,ham=one/gam,hamm1=one/gamm1
 real(kind=nr),parameter :: gamgamm1=gam*gamm1,hamhamm1=ham*hamm1,gamhamm1=gam*hamm1,hamgamm1=ham*gamm1
 real(kind=nr),parameter :: prndtl=0.71_nr,gamm1prndtli=one/(gamm1*prndtl)

 real(nr),parameter :: alpha=0.5862704032801503_nr
 real(nr),parameter :: beta=0.09549533555017055_nr
 real(nr),parameter :: aa=0.6431406736919156_nr
 real(nr),parameter :: ab=0.2586011023495066_nr
 real(nr),parameter :: ac=0.007140953479797375_nr

 real(nr),parameter :: beta20=3.558088039657372e-2_nr
 real(nr),parameter :: alpha21=4.553293665405087e-1_nr
 real(nr),parameter :: alpha23=5.485834660256362e-1_nr
 real(nr),parameter :: beta24=6.309453322148130e-2_nr
 real(nr),parameter :: a20=-1.370437811855902e-1_nr
 real(nr),parameter :: a21=-7.011151014992384e-1_nr
 real(nr),parameter :: a23=7.172015214343915e-1_nr
 real(nr),parameter :: a24=2.007061784274521e-1_nr
 real(nr),parameter :: a25=3.086671576887690e-3_nr
 real(nr),parameter :: a26=-1.220776765442191e-4_nr

 real(nr),parameter :: alpha10=9.048950052297307e-2_nr
 real(nr),parameter :: alpha12=1.626495004100577e0_nr
 real(nr),parameter :: beta13=4.991726552839925e-1_nr
 real(nr),parameter :: a10=-3.412394328068927e-1_nr
 real(nr),parameter :: a12=5.140577541790093e-1_nr
 real(nr),parameter :: a13=1.085757714299012e0_nr
 real(nr),parameter :: a14=6.848088581095074e-2_nr
 real(nr),parameter :: a15=-4.253760742663719e-3_nr
 real(nr),parameter :: a16=1.833859722828081e-4_nr

 real(nr),parameter :: alpha01=6.274536975899427e0_nr
 real(nr),parameter :: beta02=3.720061050024955e0_nr
 real(nr),parameter :: a01=-3.712509096607685e0_nr
 real(nr),parameter :: a02=6.421147876250940e0_nr
 real(nr),parameter :: a03=6.634843901364587e-1_nr
 real(nr),parameter :: a04=-2.753768568242051e-2_nr
 real(nr),parameter :: a05=-3.885109816678388e-3_nr
 real(nr),parameter :: a06=6.557485725845055e-4_nr

! real(nr),parameter :: beta20=0.03250008295108466_nr
! real(nr),parameter :: alpha21=0.3998040493524358_nr
! real(nr),parameter :: alpha23=0.7719261277615860_nr
! real(nr),parameter :: beta24=0.1626635931256900_nr
! real(nr),parameter :: a20=-0.1219006056449124_nr
! real(nr),parameter :: a21=-0.6301651351188667_nr
! real(nr),parameter :: a23=0.6521195063966084_nr
! real(nr),parameter :: a24=0.3938843551210350_nr
! real(nr),parameter :: a25=0.01904944407973912_nr
! real(nr),parameter :: a26=-0.001027260523947668_nr
!
! real(nr),parameter :: alpha10=0.08360703307833438_nr
! real(nr),parameter :: alpha12=2.058102869495757_nr
! real(nr),parameter :: beta13=0.9704052014790193_nr
! real(nr),parameter :: a10=-0.3177447290722621_nr
! real(nr),parameter :: a12=-0.02807631929593225_nr
! real(nr),parameter :: a13=1.593461635747659_nr
! real(nr),parameter :: a14=0.2533027046976367_nr
! real(nr),parameter :: a15=-0.03619652460174756_nr
! real(nr),parameter :: a16=0.004080281419108407_nr
!
! real(nr),parameter :: alpha01=5.912678614078549_nr
! real(nr),parameter :: beta02=3.775623951744012_nr
! real(nr),parameter :: a01=-3.456878182643609_nr
! real(nr),parameter :: a02=5.839043358834730_nr
! real(nr),parameter :: a03=1.015886726041007_nr
! real(nr),parameter :: a04=-0.2246526470654333_nr
! real(nr),parameter :: a05=0.08564940889936562_nr
! real(nr),parameter :: a06=-0.01836710059356763_nr

!===== ALLOCATABLE MAIN ARRAYS

 integer(kind=ni),dimension(:,:,:),allocatable :: nbbc,mbcd
 integer(kind=ni),dimension(:,:),allocatable :: nbpc,lio,jjp,jjl
 integer(kind=ni),dimension(:),allocatable :: imjp,imjl,jptag,jltag
 integer(kind=ni),dimension(:),allocatable :: li,lcsz
 integer(kind=ni),dimension(:),allocatable :: lxim,letm,lzem,lpos
 integer(kind=ni),dimension(:),allocatable :: lximb,letmb,lzemb,mo,nrr,npex
 integer(kind=int64),dimension(:),allocatable :: lhmb

 real(kind=nr),dimension(:,:),allocatable :: qo,qa,de
 real(kind=nr),dimension(:),allocatable :: txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz

 real(kind=nr),dimension(:,:),allocatable :: xim,etm,zem
 real(kind=nr),dimension(:),allocatable :: p,yaco,ayaco,amet
 real(kind=nr),dimension(:),allocatable :: sbcc

 real(kind=nr),dimension(:,:),allocatable :: rr,ss,vrss,vrst

 real(kind=nr),dimension(:,:),allocatable :: xu,yu
 real(kind=nr),dimension(:,:),allocatable :: xl,yl
 real(kind=nr),dimension(:),allocatable :: sa,sb,sc,sd,se

 real(kind=nr),dimension(:,:),allocatable :: asz
 real(kind=nr),dimension(:),allocatable :: times,vmpi

 real(kind=nr),dimension(:,:,:),pointer :: drva,drvb,send,recv,cm
 real(kind=nr),dimension(:,:),pointer :: snde,rcve

 real(kind=nr),dimension(:,:,:),allocatable,target :: drva1,drva2,drva3
 real(kind=nr),dimension(:,:,:),allocatable,target :: drvb1,drvb2,drvb3
 real(kind=nr),dimension(:,:,:),allocatable,target :: send01,send02,send03
 real(kind=nr),dimension(:,:,:),allocatable,target :: recv01,recv02,recv03
 real(kind=nr),dimension(:,:,:),allocatable,target :: send11,send12,send13
 real(kind=nr),dimension(:,:,:),allocatable,target :: recv11,recv12,recv13
 real(kind=nr),dimension(:,:,:),allocatable,target :: cm1,cm2,cm3
 real(kind=nr),dimension(:,:),allocatable,target :: send21,send22,send23
 real(kind=nr),dimension(:,:),allocatable,target :: recv21,recv22,recv23

 real(kind=ieee32),dimension(:,:),allocatable :: varm,vsamp,vmean,v2mean
 real(kind=ieee32),dimension(:),allocatable :: varr,vart,vara,varb,varmin,varmax

 character(13),dimension(:),allocatable :: ctecplt,cthead
 character(4),dimension(:),allocatable :: cfilet
 character(4),dimension(:),allocatable :: czonet

!===== CONSTANT-SIZED MAIN VARIABLES

 integer(kind=ni),dimension(0:3,3,0:1) :: ijl,ijp
 integer(kind=ni),dimension(3,0:1,0:2) :: ndf
 integer(kind=ni),dimension(0:11,ljpl) :: jlcd
 integer(kind=ni),dimension(0:7,ljpl) :: jpcd
 integer(kind=ni),dimension(3,0:5) :: nnf
 integer(kind=ni),dimension(3,3) :: ijk
 integer(kind=ni),dimension(3,0:1) :: nbc,mcd,mmcd,nsz
 integer(kind=ni),dimension(0:11) :: njl
 integer(kind=ni),dimension(0:7) :: njp
 integer(kind=ni),dimension(0:4) :: no
 integer(kind=ni),dimension(3) :: ms,me,nbsize,ijkp
 integer(kind=ni) :: lxio,leto,lzeo,lxi,let,lze,lmx,lim,lsz,nrecs,nrecd
 integer(kind=ni) :: i,ii,is,ie,ip,iq,j,jj,js,je,jp,jq,jk,k,kk,ks,ke,kp,l,lh,ll,lmpi,lp,lq,ltomb
 integer(kind=ni) :: m,ma,mb,mm,mp,mq,mbk,mps,mpe,n,ndt,nn,nk,ns,ne,np,nq,nt,nz,ndati,nsigi,nfile,nshock
 integer(kind=ni) :: nts,nscrn,nsgnl,ndata,ndatafl,ndataav,nviscous,nvisc,nkrk,nsmf,nrestart,nextrabc,nextgcic
 integer(kind=int64) :: nlmx,llmb,llmo,lis,lie,ljs,lje,llis,llie

 real(kind=nr),dimension(0:lmp,0:1,0:1) :: pbci,pbco
 real(kind=nr),dimension(-2:2,0:2,0:1) :: albed,albef
 real(kind=nr),dimension(mbci,mbci) :: cbca,cbcs
 real(kind=nr),dimension(5,5) :: xt
 real(kind=nr),dimension(0:5,0:2) :: dbc,fbc
 real(kind=nr),dimension(0:1,0:1) :: pbcot
 real(kind=nr),dimension(0:lmp) :: sap
 real(kind=nr),dimension(mbci) :: rbci,sbci
 real(kind=nr),dimension(5) :: cha,dha
 real(kind=nr),dimension(3) :: ve,dm,rv,uoo,umf
 real(kind=nr) :: alphf,betf,fa,fb,fc
 real(kind=nr) :: ra0,ra1,ra2,ra3,res,fctr,dfdt,shockc
 real(kind=nr) :: reuoo,reaoo,tempoo,amachoo,aoa1,aoa2,wtemp,cfl,tmax,timf,fltk,fltke,fltexr,shockr,dto
 real(kind=nr) :: aoo,amach1,amach2,amach3,srefoo,srefp1dre
 real(kind=nr) :: dt,dts,dte,dtk,dtko,dtsum,dtwi,timo,tsam,wts,wte,wtime
 real(kind=nr) :: vn,vs,hv2,ao,bo,co,ho,aoi,rhoi,progmf,sqrtrema,sqrtremai

 character(1),dimension(0:4) :: cno
 character(3) :: cndata
 character(5) :: cnnode
 character(13) :: coutput
 character(16) :: cinput
 character(17) :: cgrid
 character(18) :: cdata
 character(20) :: crestart

!===== INTEGER VARIABLES FOR MPI COMMANDS

 integer(kind=ni),dimension(:,:),allocatable :: ista
 integer(kind=ni),dimension(:),allocatable :: ireq
 integer(kind=ni) :: ir,mpro,npro,myid,itag,info,icom,ierr

!=====

 end module mainvar3d

!*****