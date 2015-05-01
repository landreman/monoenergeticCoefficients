!#include <finclude/petscvecdef.h>

subroutine populateRHS()

  !use petscvec
  use indices
  use stel_kinds
  use variables, only: G, I, masterProc, B, dBdtheta, dBdzeta, Ntheta, Nzeta, Nxi, brhs
  use cyclic_red, only: ns0, nsn

  implicit none

  integer :: itheta, izeta, L, index
  real(rprec) :: valueToInsert

  if (masterProc) then
     print *,"Entering populateRHS"
  end if

  brhs = 0.0d+0

  L = 0
  if ((L+1 >= ns0) .and. (L+1 <= nsn)) then
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           index = getBlockIndex(itheta,izeta)
           valueToInsert = (1.0d+0)/B(itheta,izeta)*(G*dBdtheta(itheta,izeta)-I*dBdzeta(itheta,izeta))
           brhs(index,L+1) = valueToInsert*4/(3.0d+0)
        end do
     end do
  end if

  L = 2
  if ((L+1 >= ns0) .and. (L+1 <= nsn)) then
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           index = getBlockIndex(itheta,izeta)
           valueToInsert = (1.0d+0)/B(itheta,izeta)*(G*dBdtheta(itheta,izeta)-I*dBdzeta(itheta,izeta))
           brhs(index,L+1) = valueToInsert*2/(3.0d+0)
        end do
     end do
  end if

end subroutine populateRHS
