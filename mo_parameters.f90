MODULE mo_parameters
   use mo_kind
   implicit none
   PUBLIC

   !===== CONSTANT PARAMETERS

   integer(kind=ni),parameter :: nrall=0,nrone=1,n45no=0,n45go=1
   integer(kind=ni),parameter :: liofs=16,liofl=24

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

END MODULE mo_parameters