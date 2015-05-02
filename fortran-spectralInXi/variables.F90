
module variables

  implicit none

!#include <finclude/petscsys.h>
#include <finclude/petscsysdef.h>

  PetscInt :: Ntheta, Nzeta, Nxi, Nperiods, helicity_l, matrixSize
  PetscReal :: nu, epsilon_t, epsilon_h, iota, G, I
  PetscOffset :: offset
  PetscInt :: myRank, numProcs
  logical :: masterProc

  PetscScalar, dimension(:), allocatable :: theta, zeta
  PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights
  PetscScalar, dimension(:,:), allocatable :: ddtheta
  PetscScalar, dimension(:,:), allocatable :: ddzeta
  PetscScalar, dimension(:,:), allocatable :: B, dBdtheta, dBdzeta
  PetscScalar :: diagonalShift

  PetscScalar, parameter :: pi = 3.14159265358979d+0
  PetscScalar, parameter :: sqrtpi = 1.77245385090552d+0

  PetscInt :: thetaGridScheme, zetaGridScheme

  integer :: clockStart, clockRate

end module variables
