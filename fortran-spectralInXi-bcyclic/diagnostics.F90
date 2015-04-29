!#include <finclude/petsckspdef.h>

subroutine diagnostics(solution)

  !use petscksp
  use stel_kinds
  use variables
  use indices

  implicit none

  integer :: itheta, izeta, L, index
  real(rprec) :: flux, flow, VPrime, spatialPart
  double precision :: sendBuffer, recvBuffer
  integer :: ierr

  if (masterProc) then
     print *,"Entering diagnostics"
  end if

  flux = 0
  flow = 0

  do itheta = 1,Ntheta
     do izeta = 1,Nzeta
        index = getBlockIndex(itheta,izeta)
        if ((index >= ns0) .and. (index <= nsn)) then
           spatialPart = (G * dBdtheta(itheta,izeta) - I * dBdzeta(itheta,izeta))/(B(itheta,izeta) ** 3) &
                * thetaWeights(itheta)*zetaWeights(izeta)

           L = 0
           flux = flux + brhs(index,L+1) * spatialPart * 8/3

           L = 2
           flux = flux + brhs(index,L+1) * spatialPart * 4/15

           L = 1
           flow = flow + brhs(index,L+1) / B(itheta,izeta) * thetaWeights(itheta)*zetaWeights(izeta)
        end if
     end do
  end do
          
  VPrime = 0
  do itheta = 1,Ntheta
     do izeta = 1,Nzeta
        VPrime = VPrime + thetaWeights(itheta)*zetaWeights(izeta)/(B(itheta,izeta) ** 2)
     end do
  end do
     
  flow = flow * 4 / (3*sqrtpi*G*VPrime)
  flux = -2 / (sqrtpi*G*G*VPrime)*flux

  ! Sum results from all procs:                                                                    
  sendBuffer = flux
  call MPI_Reduce(sendBuffer, recvBuffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  flux = recvBuffer(1)

  sendBuffer = flow
  call MPI_Reduce(sendBuffer, recvBuffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  flow = recvBuffer(1)

  if (masterProc) then   
     print *,"Results: VPrime = ",VPrime,", flux = ",flux,", flow = ",flow
  end if

end subroutine diagnostics
