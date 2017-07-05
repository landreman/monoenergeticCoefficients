
module variables

  implicit none

!#include <finclude/petscsys.h>
#include <petsc/finclude/petscmatdef.h>

  PetscInt :: Ntheta, Nzeta, Nxi, Nperiods, helicity_l, matrixSize
  PetscReal :: nu, nu_hat, E, epsilon_t, epsilon_h, iota, G, I
  PetscOffset :: offset
  PetscInt :: myRank, numProcs
  logical :: masterProc

  PetscScalar, dimension(:), allocatable :: theta, zeta
  PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights
  PetscScalar, dimension(:,:), allocatable :: ddtheta
  PetscScalar, dimension(:,:), allocatable :: ddzeta
  PetscScalar, dimension(:,:), allocatable :: B, dBdtheta, dBdzeta, sqrt_g
  PetscScalar :: diagonalShift, thetaWeight, zetaWeight, VPrime, FSAB2

  PetscScalar, parameter :: pi = 3.14159265358979d+0
  PetscScalar, parameter :: sqrtpi = 1.77245385090552d+0
  PetscScalar, parameter :: one = 1.0d+0
  PetscScalar, parameter :: zero = 0.0d+0

  PetscInt :: thetaGridScheme, zetaGridScheme

  integer :: clockStart, clockRate

  integer :: L_scaling_option
  PetscScalar, dimension(:), allocatable :: L_scaling
  Vec :: weights_vec, CHat_inverse
  Mat :: V_matrix, pick_out_p0_matrix, injection_matrix
  integer :: fieldsplit_option, constraint_option

end module variables
