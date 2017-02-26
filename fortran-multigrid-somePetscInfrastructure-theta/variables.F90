#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

module variables

  implicit none

!!#include <finclude/petscsys.h>
!#include <finclude/petscsysdef.h>

  KSP :: main_ksp
  PC :: preconditioner_context

  PetscInt :: Ntheta, Nzeta, Nxi, Nperiods, helicity_l, matrixSize
  PetscReal :: nu, E, epsilon_t, epsilon_h, iota, G, I
  PetscOffset :: offset
  PetscInt :: myRank, numProcs
  logical :: masterProc
  PetscScalar :: zetaMax, VPrime

  PetscScalar, parameter :: pi = 3.14159265358979d+0
  PetscScalar, parameter :: sqrtpi = 1.77245385090552d+0
  PetscScalar, parameter :: zero = 0.0d+0, one = 1.0d+0, two = 2.0d+0

  PetscInt :: theta_derivative_option, preconditioner_theta_derivative_option
  PetscInt :: zeta_derivative_option, preconditioner_zeta_derivative_option
  PetscInt :: xi_derivative_option, preconditioner_xi_derivative_option
  PetscInt :: pitch_angle_scattering_option, preconditioner_pitch_angle_scattering_option
  PetscInt :: xi_quadrature_option, constraint_option

  PetscBool :: coarsen_theta, coarsen_zeta, coarsen_xi
  PetscInt :: N_levels
  PetscInt, allocatable, dimension(:) :: Ntheta_levels, Nzeta_levels, Nxi_levels
  PetscInt :: Ntheta_min, Nzeta_min, Nxi_min

  PetscInt :: smoothing_option, restriction_option, coarsen_option=1, defect_option

  type :: multigrid_level
     PetscInt :: Ntheta, Nzeta, Nxi, matrixSize

     PetscScalar, allocatable, dimension(:) :: theta, zeta, xi
     PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights, xiWeights

     PetscScalar, dimension(:,:), allocatable :: ddtheta_plus
     PetscScalar, dimension(:,:), allocatable :: ddtheta_minus
     PetscScalar, dimension(:,:), allocatable :: ddtheta_plus_preconditioner
     PetscScalar, dimension(:,:), allocatable :: ddtheta_minus_preconditioner
     
     PetscScalar, dimension(:,:), allocatable :: ddzeta_plus
     PetscScalar, dimension(:,:), allocatable :: ddzeta_minus
     PetscScalar, dimension(:,:), allocatable :: ddzeta_plus_preconditioner
     PetscScalar, dimension(:,:), allocatable :: ddzeta_minus_preconditioner

     PetscScalar, dimension(:,:), allocatable :: ddxi_plus
     PetscScalar, dimension(:,:), allocatable :: ddxi_minus
     PetscScalar, dimension(:,:), allocatable :: ddxi_plus_preconditioner
     PetscScalar, dimension(:,:), allocatable :: ddxi_minus_preconditioner

     PetscScalar, dimension(:,:), allocatable :: pitch_angle_scattering_operator
     PetscScalar, dimension(:,:), allocatable :: pitch_angle_scattering_operator_preconditioner

     PetscScalar, dimension(:,:), allocatable :: B, dBdtheta, dBdzeta

     PetscInt :: ithetaMin, ithetaMax
     PetscInt :: izetaMin, izetaMax

     Mat :: low_order_matrix, high_order_matrix, mixed_order_matrix

     Mat :: smoothing_matrix

  end type multigrid_level

  type (multigrid_level), allocatable, dimension(:), target :: levels

  Mat, allocatable, dimension(:) :: multigrid_prolongation_matrices, multigrid_restriction_matrices

  PetscReal :: flux, flow
  PetscReal :: zeta_diffusion, theta_diffusion

end module variables
