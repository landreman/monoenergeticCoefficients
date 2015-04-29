!#include <finclude/petscvecdef.h>

subroutine populateRHS()

  !use petscvec
  use indices
  use stel_kinds
  use variables, only: G, I, masterProc, B, dBdtheta, dBdzeta, Ntheta, Nzeta, Nxi, brhs, ns0, nsn

  implicit none

  integer :: itheta, izeta, L, index
  real(rprec) :: valueToInsert

  if (masterProc) then
     print *,"Entering populateRHS"
  end if

  brhs = 0.0d+0

  do itheta = 1,Ntheta
     do izeta = 1,Nzeta
        index = getBlockIndex(itheta,izeta)
        if ((index >= ns0) .and. (index <= nsn)) then
           valueToInsert = (1.0d+0)/B(itheta,izeta)*(G*dBdtheta(itheta,izeta)-I*dBdzeta(itheta,izeta))

           L = 0
           brhs(index,L+1) = valueToInsert*4/(3.0d+0)

           L = 2
           brhs(index,L+1) = valueToInsert*2/(3.0d+0)
        end if
     end do
  end do

end subroutine populateRHS
