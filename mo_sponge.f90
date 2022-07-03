!*****
!***** SPONGE IMPLEMENTATION
!*****

MODULE mo_sponge
   use mo_kind,       ONLY : ni, nr
   use mo_parameters, ONLY : two, zero, pi, sml, one, half, hamhamm1
   use mo_vars,       ONLY : de, ss, qa
   use mo_grid,       ONLY : yaco
   use mo_gridgen,    ONLY : szth0, szth1, skew, doml0, doml1, domh

   IMPLICIT NONE
   PUBLIC

   real(kind=nr),    private, dimension(:), allocatable :: asz, bsz
   integer(kind=ni), private, dimension(:), allocatable :: lcsz
   integer(kind=ni), private :: lsz
   real(kind=nr),    private :: szco

   CONTAINS

   subroutine read_input_sponge
      character(16) :: ccinput

      open(9,file='input.sponge',status='old')
      read(9,*) ccinput,szco
      close(9)

   end subroutine read_input_sponge

   !===== SETTING UP SPONGE ZONE PARAMETERS

   subroutine spongeup(lmx)
      integer(kind=ni), INTENT(IN) :: lmx
      integer(kind=ni) :: l, ll
      real(kind=nr) :: ra0, ra1, ra2, ra3
      real(kind=nr) :: tmpa, tmpb

      ll=-1
      ra2=skew/domh
      tmpa=pi/szth0
      tmpb=pi/szth1
      do l=0,lmx
         ra3=ra2*ss(l,2)
         ra0=tmpa*(ss(l,1)-(ra3-doml0+szth0))
         ra1=tmpb*(ra3+doml1-szth1-ss(l,1))
         de(l,1)=szco*half*(two+cos(max(min(ra0,pi),zero))+cos(max(min(ra1,pi),zero)))
         de(l,2)=szco*half*(one+cos(max(min(ra0,pi),zero)))
         if(de(l,1)>sml) then
            ll=ll+1
            de(ll,5)=l+sml
         end if
      end do
      lsz=ll
      if(lsz/=-1) then
         allocate(lcsz(0:lsz),asz(0:lsz),bsz(0:lsz))
         do ll=0,lsz
            l=de(ll,5)
            lcsz(ll)=l
            asz(ll)=de(l,1)/yaco(l)
            bsz(ll)=de(l,2)/yaco(l)
         end do
      end if

   end subroutine spongeup

 !===== SPONGE IMPLEMENTATION

   subroutine spongego
      integer(kind=ni) :: l, ll

      do ll=0,lsz
         l=lcsz(ll)
         de(l,1)=de(l,1)+asz(ll)*(qa(l,1)-one)
         de(l,2:4)=de(l,2:4)+bsz(ll)*(qa(l,2:4)-zero)
         de(l,5)=de(l,5)+asz(ll)*(qa(l,5)-hamhamm1)
      end do

   end subroutine spongego

END MODULE mo_sponge