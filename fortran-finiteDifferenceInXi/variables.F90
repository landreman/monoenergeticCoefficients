#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

module variables

  implicit none

!!#include <finclude/petscsys.h>
!#include <finclude/petscsysdef.h>

  PetscInt :: Ntheta, Nzeta, Nxi, Nperiods, helicity_l, matrixSize
  PetscReal :: nu, E, epsilon_t, epsilon_h, iota, G, I
  PetscOffset :: offset
  PetscInt :: myRank, numProcs
  PetscInt :: ithetaMin, ithetaMax, localNtheta
  PetscInt :: izetaMin, izetaMax, localNzeta
  logical :: masterProc

  PetscScalar, dimension(:), allocatable :: theta, zeta, xi
  PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights, xiWeights

  PetscScalar, dimension(:,:), allocatable, target :: ddtheta_plus
  PetscScalar, dimension(:,:), allocatable, target :: ddtheta_minus
  PetscScalar, dimension(:,:), allocatable, target :: ddtheta_plus_preconditioner
  PetscScalar, dimension(:,:), allocatable, target :: ddtheta_minus_preconditioner

  PetscScalar, dimension(:,:), allocatable, target :: ddzeta_plus
  PetscScalar, dimension(:,:), allocatable, target :: ddzeta_minus
  PetscScalar, dimension(:,:), allocatable, target :: ddzeta_plus_preconditioner
  PetscScalar, dimension(:,:), allocatable, target :: ddzeta_minus_preconditioner

  PetscScalar, dimension(:,:), allocatable, target :: ddxi_plus
  PetscScalar, dimension(:,:), allocatable, target :: ddxi_minus
  PetscScalar, dimension(:,:), allocatable, target :: ddxi_plus_preconditioner
  PetscScalar, dimension(:,:), allocatable, target :: ddxi_minus_preconditioner

  PetscScalar, dimension(:,:), allocatable, target :: pitch_angle_scattering_operator
  PetscScalar, dimension(:,:), allocatable, target :: pitch_angle_scattering_operator_preconditioner

  PetscScalar, dimension(:,:), allocatable :: B, dBdtheta, dBdzeta
  PetscReal :: VPrime, FSAB2

  PetscScalar, parameter :: pi = 3.14159265358979d+0
  PetscScalar, parameter :: sqrtpi = 1.77245385090552d+0
  PetscScalar, parameter :: zero = 0.0d+0, one = 1.0d+0, two = 2.0d+0

  PetscInt :: theta_derivative_option, preconditioner_theta_derivative_option
  PetscInt :: zeta_derivative_option, preconditioner_zeta_derivative_option
  PetscInt :: xi_derivative_option, preconditioner_xi_derivative_option
  PetscInt :: pitch_angle_scattering_option, preconditioner_pitch_angle_scattering_option
  PetscInt :: xi_quadrature_option, constraint_option

end module variables
