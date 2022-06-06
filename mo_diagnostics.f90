MODULE mo_diagnostics
   use mo_kind,       ONLY : ni
   use mo_parameters, ONLY : one, nrone, n45no
   use mo_vars,       ONLY : rr, ss, de, qa, yaco, xim, etm, zem
   use mo_domdcomp,   ONLY : ijk, nbsize, nbc, mcd
   use mo_numerics,   ONLY : mpigo, deriv
   public

   contains

   !===== SUBROUTINE FOR CALCULATING VORTICITY

   subroutine vorti
      integer(kind=ni) :: m

      ss(:,1)=one/qa(:,1); de(:,1:3)=zero

      rr(:,1)=ss(:,1)*qa(:,2)
      m=1
      call mpigo(rr, ijk, nbc, mcd, nbsize,0,nrone,n45no,m)
      call deriv(rr, lxi, let, lze, ijk, 3, 1, m)
      call deriv(rr, lxi, let, lze, ijk, 2, 1, m)
      call deriv(rr, lxi, let, lze, ijk, 1, 1, m)
      de(:,2)=de(:,2)+rr(:,1)*xim(:,3)+rr(:,2)*etm(:,3)+rr(:,3)*zem(:,3)
      de(:,3)=de(:,3)-rr(:,1)*xim(:,2)-rr(:,2)*etm(:,2)-rr(:,3)*zem(:,2)

      rr(:,1)=ss(:,1)*qa(:,3)
      m=2
      call mpigo(rr, ijk, nbc, mcd, nbsize,0,nrone,n45no,m)
      call deriv(rr, lxi, let, lze, ijk, 3, 1, m)
      call deriv(rr, lxi, let, lze, ijk, 2, 1, m)
      call deriv(rr, lxi, let, lze, ijk, 1, 1, m)
      de(:,3)=de(:,3)+rr(:,1)*xim(:,1)+rr(:,2)*etm(:,1)+rr(:,3)*zem(:,1)
      de(:,1)=de(:,1)-rr(:,1)*xim(:,3)-rr(:,2)*etm(:,3)-rr(:,3)*zem(:,3)

      rr(:,1)=ss(:,1)*qa(:,4)
      m=3
      call mpigo(rr, ijk, nbc, mcd, nbsize,0,nrone,n45no,m)
      call deriv(rr, lxi, let, lze, ijk, 3, 1, m)
      call deriv(rr, lxi, let, lze, ijk, 2, 1, m)
      call deriv(rr, lxi, let, lze, ijk, 1, 1, m)
      de(:,1)=de(:,1)+rr(:,1)*xim(:,2)+rr(:,2)*etm(:,2)+rr(:,3)*zem(:,2)
      de(:,2)=de(:,2)-rr(:,1)*xim(:,1)-rr(:,2)*etm(:,1)-rr(:,3)*zem(:,1)

      de(:,1)=de(:,1)*yaco(:); de(:,2)=de(:,2)*yaco(:); de(:,3)=de(:,3)*yaco(:)

   end subroutine vorti

END MODULE mo_diagnostics