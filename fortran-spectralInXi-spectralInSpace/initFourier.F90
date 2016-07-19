!#include <petsc/finclude/petscsysdef.h>

subroutine initFourier()

  !use petscsys
  use variables

  implicit none

  integer :: j

  print *,"Entering initFourier"

  ! Set up real-space grids for theta and zeta:
  Ntheta = 2*mmax+3
  Nzeta  = 2*nmax+3
  allocate(theta(Ntheta))
  allocate(zeta(Nzeta))
  theta = [( 2*pi*j/Ntheta, j=0,Ntheta-1 )]
  zeta  = [( 2*pi*j/Nzeta,  j=0,Nzeta-1  )] / Nperiods
  thetaWeight = 2*pi/Ntheta
  zetaWeight  = 2*pi/Nzeta

  ! Evaluate |B| and its derivatives on this grid:
  call computeB()

  call chooseFourierModes()

  allocate(ddtheta(NFourier2,NFourier2))
  allocate(ddzeta (NFourier2,NFourier2))
  call FourierDifferentiationMatrices(NFourier, xm, xn, ddtheta, ddzeta)

  matrixSize = NFourier2*Nxi+1

end subroutine initFourier


