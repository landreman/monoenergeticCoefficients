
module variables

  implicit none

!#include <finclude/petscsys.h>
#include <petsc/finclude/petscsysdef.h>

  PetscInt :: NFourier, NFourier2, Nxi, Nperiods, helicity_l, matrixSize
  PetscReal :: nu, epsilon_t, epsilon_h, iota, G, I
  PetscOffset :: offset
  PetscInt :: myRank, numProcs
  logical :: masterProc

  integer :: mmax, nmax, Ntheta, Nzeta
  integer, dimension(:), allocatable :: xm, xn
  PetscScalar, dimension(:), allocatable :: theta, zeta
  PetscScalar :: thetaWeight, zetaWeight
  PetscScalar, dimension(:,:), allocatable :: ddtheta, ddzeta
  PetscScalar, dimension(:,:), allocatable :: B, dBdtheta, dBdzeta

  PetscScalar, parameter :: pi = 3.14159265358979d+0
  PetscScalar, parameter :: sqrtpi = 1.77245385090552d+0

  integer :: clockStart, clockRate

end module variables
