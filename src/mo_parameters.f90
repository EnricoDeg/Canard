MODULE mo_parameters
   use mo_kind, ONLY : ni, nr
   implicit none
   PUBLIC

   !===== CONSTANT PARAMETERS

   !----- halo exchange
   integer(kind=ni), parameter :: nrall=0
   integer(kind=ni), parameter :: nrone=1
   integer(kind=ni), parameter :: n45no=0
   integer(kind=ni), parameter :: n45go=1

   !----- numbers
   real(kind=nr),    parameter :: zero      = 0
   real(kind=nr),    parameter :: quarter   = 0.25_nr
   real(kind=nr),    parameter :: half      = 0.5_nr
   real(kind=nr),    parameter :: one       = 1
   real(kind=nr),    parameter :: two       = 2
   real(kind=nr),    parameter :: three     = 3
   real(kind=nr),    parameter :: four      = 4
   real(kind=nr),    parameter :: five      = 5
   real(kind=nr),    parameter :: onethird  = one / three
   real(kind=nr),    parameter :: twothirds = two / three
   real(kind=nr),    parameter :: pi        = acos(-one)
   real(kind=nr),    parameter :: halfpi    = pi / two
   real(kind=nr),    parameter :: twopi     = two * pi
   real(kind=nr),    parameter :: sqrt2     = sqrt(two)
   real(kind=nr),    parameter :: sqrt2i    = one / sqrt2
   real(kind=nr),    parameter :: sml       = 1.0e-8_nr
   real(kind=nr),    parameter :: free      = 1.0e+6_nr

   !----- boundary contidions
   integer(kind=ni), parameter :: BC_NON_REFLECTIVE   = 10
   integer(kind=ni), parameter :: BC_WALL_INVISCID    = 20
   integer(kind=ni), parameter :: BC_WALL_VISCOUS     = 25
   integer(kind=ni), parameter :: BC_INTER_CURV       = 30
   integer(kind=ni), parameter :: BC_INTER_STRAIGHT   = 35
   integer(kind=ni), parameter :: BC_INTER_SUBDOMAINS = 40
   integer(kind=ni), parameter :: BC_PERIODIC         = 45

   !----- physics
   real(kind=nr),    parameter :: gam          = 1.4_nr
   real(kind=nr),    parameter :: gamm1        = gam - one
   real(kind=nr),    parameter :: ham          = one / gam
   real(kind=nr),    parameter :: hamm1        = one / gamm1
   real(kind=nr),    parameter :: hamhamm1     = ham * hamm1
   real(kind=nr),    parameter :: prndtl       = 0.71_nr
   real(kind=nr),    parameter :: gamm1prndtli = one / ( gamm1 * prndtl )

   !----- numerics
   real(kind=nr),    parameter :: alpha = 4 / 9.0_nr
   real(kind=nr),    parameter :: beta  = 1 / 36.0_nr
   real(kind=nr),    parameter :: aa    = 20 / 27.0_nr
   real(kind=nr),    parameter :: ab    = 25 / 216.0_nr

   real(kind=nr),    parameter :: alpha10 = 0.09166666665902128_nr
   real(kind=nr),    parameter :: alpha12 = 1.3500000000438646_nr
   real(kind=nr),    parameter :: beta13  = 0.2666666666794667_nr
   real(kind=nr),    parameter :: a10     = -0.3506944444280639_nr
   real(kind=nr),    parameter :: a12     = 0.8249999999974599_nr
   real(kind=nr),    parameter :: a13     = 0.7444444444773477_nr
   real(kind=nr),    parameter :: a14     = 0.014583333334139971_nr

   real(kind=nr),    parameter :: alpha01 = 8.000000000449464_nr
   real(kind=nr),    parameter :: beta02  = 6.000000000542_nr
   real(kind=nr),    parameter :: a01     = -6.666666667794732_nr
   real(kind=nr),    parameter :: a02     = 9.000000000052879_nr
   real(kind=nr),    parameter :: a03     = 1.333333333178555_nr
   real(kind=nr),    parameter :: a04     = -0.0833333333455545_nr

END MODULE mo_parameters