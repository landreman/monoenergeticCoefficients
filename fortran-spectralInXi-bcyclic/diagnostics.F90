!#include <finclude/petsckspdef.h>

subroutine diagnostics()

  !use petscksp
  use cyclic_red, only: nsn, ns0
  use indices
  use stel_kinds
  use variables

  implicit none

  INCLUDE 'mpif.h' 
  integer :: itheta, izeta, L, index
  real(rprec) :: flux, flow, VPrime, spatialPart
  double precision :: sendBuffer(1), recvBuffer(1)
  integer :: ierr

  integer :: clockStop
  real :: elapsedTime

  character(len=6) :: filename="output"
  integer :: fileUnit=11, didFileAccessWork

  if (masterProc) then
     print *,"Entering diagnostics."
  end if

!!$  print *,"Here comes solution."
!!$  do L=0,(Nxi-1)
!!$     print *,"L=",L
!!$     do itheta = 1,Ntheta
!!$        print *,"  itheta=",itheta
!!$        do izeta = 1,Nzeta
!!$           index = getBlockIndex(itheta,izeta)
!!$           print *,brhs(index,L+1)
!!$        end do
!!$     end do
!!$  end do


  flux = 0
  flow = 0

  L = 0
  if ((L+1 >= ns0) .and. (L+1 <= nsn)) then
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           index = getBlockIndex(itheta,izeta)
           spatialPart = (G * dBdtheta(itheta,izeta) - I * dBdzeta(itheta,izeta))/(B(itheta,izeta) ** 3) &
                * thetaWeights(itheta)*zetaWeights(izeta)
           flux = flux + brhs(index,L+1) * spatialPart * 8/3
        end do
     end do
  end if

  L = 2
  if ((L+1 >= ns0) .and. (L+1 <= nsn)) then
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           index = getBlockIndex(itheta,izeta)
           spatialPart = (G * dBdtheta(itheta,izeta) - I * dBdzeta(itheta,izeta))/(B(itheta,izeta) ** 3) &
                * thetaWeights(itheta)*zetaWeights(izeta)
           flux = flux + brhs(index,L+1) * spatialPart * 4/15
        end do
     end do
  end if

  L = 1
  if ((L+1 >= ns0) .and. (L+1 <= nsn)) then
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           index = getBlockIndex(itheta,izeta)
           flow = flow + brhs(index,L+1) / B(itheta,izeta) * thetaWeights(itheta)*zetaWeights(izeta)
        end do
     end do
  end if
          
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

  call system_clock(clockStop)
  elapsedTime = real(clockStop-clockStart)/clockRate

  if (masterProc) then   
     print *,"Results:"
     print *,"  Flux = ",flux
     print *,"  Flow = ",flow
     print *,"  Time for solve (seconds) = ",elapsedTime


     open(unit=fileUnit, file=filename, action="write", iostat = didFileAccessWork)
     if (didFileAccessWork /= 0) then
        print *,"Error opening ", trim(filename)
        stop
     else
        write (unit=fileUnit,fmt="(a,es22.15)") "Flux = ", flux
        write (unit=fileUnit,fmt="(a,es22.15)") "Flow = ", flow
        write (unit=fileUnit,fmt="(a,es22.15)") "Time = ", elapsedTime
        close(unit=fileUnit)
     end if
  end if

end subroutine diagnostics
