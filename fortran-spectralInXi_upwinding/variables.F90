
module variables

  implicit none

!#include <finclude/petscsys.h>
#include <petsc/finclude/petscsysdef.h>

  PetscInt :: Ntheta, Nzeta, Nxi, Nperiods, helicity_l, matrixSize
  PetscReal :: nu, epsilon_t, epsilon_h, iota, G, I
  PetscOffset :: offset
  PetscInt :: myRank, numProcs
  logical :: masterProc

  PetscScalar, dimension(:), allocatable :: theta, zeta
  PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights
  PetscScalar, dimension(:,:), allocatable :: ddtheta_sum, ddtheta_difference
  PetscScalar, dimension(:,:), allocatable :: ddzeta_sum, ddzeta_difference
  PetscScalar, dimension(:,:), allocatable :: ddtheta_sum_preconditioner, ddtheta_difference_preconditioner
  PetscScalar, dimension(:,:), allocatable :: ddzeta_sum_preconditioner, ddzeta_difference_preconditioner
  PetscScalar, dimension(:,:), allocatable :: B, dBdtheta, dBdzeta

  PetscScalar, parameter :: pi = 3.14159265358979d+0
  PetscScalar, parameter :: sqrtpi = 1.77245385090552d+0
  PetscScalar, parameter :: zero = 0.0d+0
  PetscScalar, parameter :: one = 1.0d+0
  PetscScalar, parameter :: two = 2.0d+0

  PetscInt :: theta_derivative_option, preconditioner_theta_derivative_option
  PetscInt :: zeta_derivative_option, preconditioner_zeta_derivative_option

  integer :: clockStart, clockRate
  integer :: constraint_option

end module variables
