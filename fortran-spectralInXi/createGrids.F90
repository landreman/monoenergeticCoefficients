#include <petsc/finclude/petscsysdef.h>

subroutine createGrids()

  use petscsys

  use variables, only: Nperiods, pi, matrixSize, &
       B, dBdtheta, dBdzeta, &
       Ntheta, Nzeta, Nxi, theta, zeta, &
       thetaWeights, zetaWeights, &
       ddtheta, ddzeta, &
       thetaGridScheme, zetaGridScheme

  implicit none

  PetscScalar, dimension(:,:), allocatable :: d2dtheta2
  PetscScalar, dimension(:,:), allocatable :: d2dzeta2
  PetscScalar :: zetaMax

  matrixSize = Ntheta*Nzeta*Nxi

  ! *****************************************************
  ! Set up theta grid
  ! *****************************************************

  allocate(theta(Ntheta))
  allocate(thetaWeights(Ntheta))
  allocate(ddtheta(Ntheta,Ntheta))
  allocate(d2dtheta2(Ntheta,Ntheta))

  call uniformDiffMatrices(Ntheta, 0.0d+0, 2*pi, thetaGridScheme, theta, thetaWeights, &
       ddtheta, d2dtheta2)

  deallocate(d2dtheta2)

  ! *****************************************************
  ! Set up zeta grid
  ! *****************************************************

  allocate(zeta(Nzeta))
  allocate(zetaWeights(Nzeta))
  allocate(ddzeta(Nzeta,Nzeta))
  allocate(d2dzeta2(Nzeta,Nzeta))

  zetaMax = 2*pi/Nperiods

  call uniformDiffMatrices(Nzeta, 0.0d+0, zetaMax, zetaGridScheme, zeta, zetaWeights, &
       ddzeta, d2dzeta2)

  deallocate(d2dzeta2)

  ! *****************************************************
  ! Compute B and its derivatives
  ! *****************************************************

  allocate(B(Ntheta,Nzeta))
  allocate(dBdtheta(Ntheta,Nzeta))
  allocate(dBdzeta(Ntheta,Nzeta))
  call computeB()

end subroutine createGrids


