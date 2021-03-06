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

  KSP :: outer_ksp, inner_ksp
  PC :: outer_preconditioner, inner_preconditioner

  PetscInt :: Ntheta, Nzeta, Nxi, Nperiods, helicity_l, matrixSize
  PetscReal :: nu, E, epsilon_t, epsilon_h, iota, G, I
  PetscOffset :: offset
  PetscInt :: myRank, numProcs
  logical :: masterProc
  PetscScalar :: zetaMax, VPrime, FSAB2

  PetscScalar, parameter :: pi = 3.14159265358979d+0
  PetscScalar, parameter :: sqrtpi = 1.77245385090552d+0
  PetscScalar, parameter :: zero = 0.0d+0, one = 1.0d+0, two = 2.0d+0

  PetscInt :: theta_derivative_option, preconditioner_theta_derivative_option
  PetscInt :: zeta_derivative_option, preconditioner_zeta_derivative_option
  PetscInt :: constraint_option

  PetscBool :: coarsen_theta, coarsen_zeta, coarsen_xi
  PetscInt :: N_levels
  PetscInt, allocatable, dimension(:) :: Ntheta_levels, Nzeta_levels, Nxi_levels
  PetscInt :: Ntheta_min, Nzeta_min, Nxi_min

  PetscInt :: smoothing_option, restriction_option, coarsen_option=1, defect_option

  type :: multigrid_level
     PetscInt :: Ntheta, Nzeta, Nxi, matrixSize

     PetscScalar, allocatable, dimension(:) :: theta, zeta
     PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights

     PetscScalar, dimension(:,:), allocatable :: ddtheta_plus
     PetscScalar, dimension(:,:), allocatable :: ddtheta_minus
     PetscScalar, dimension(:,:), allocatable :: ddtheta_plus_preconditioner
     PetscScalar, dimension(:,:), allocatable :: ddtheta_minus_preconditioner
     
     PetscScalar, dimension(:,:), allocatable :: ddtheta_sum
     PetscScalar, dimension(:,:), allocatable :: ddtheta_difference
     PetscScalar, dimension(:,:), allocatable :: ddtheta_sum_preconditioner
     PetscScalar, dimension(:,:), allocatable :: ddtheta_difference_preconditioner
     
     PetscScalar, dimension(:,:), allocatable :: ddzeta_plus
     PetscScalar, dimension(:,:), allocatable :: ddzeta_minus
     PetscScalar, dimension(:,:), allocatable :: ddzeta_plus_preconditioner
     PetscScalar, dimension(:,:), allocatable :: ddzeta_minus_preconditioner

     PetscScalar, dimension(:,:), allocatable :: ddzeta_sum
     PetscScalar, dimension(:,:), allocatable :: ddzeta_difference
     PetscScalar, dimension(:,:), allocatable :: ddzeta_sum_preconditioner
     PetscScalar, dimension(:,:), allocatable :: ddzeta_difference_preconditioner

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
  logical :: shift_L0_in_smoother

  integer :: f_scaling_option
  integer :: L_scaling_option
  PetscScalar, dimension(:), allocatable :: f_scaling, L_scaling
  PetscScalar :: omega, upwinding_scale_factor, preconditioner_upwinding_scale_factor

  integer :: preconditioning_option ! 1=just apply inner_KSP, 2=also explicitly handle the sources and constraints.

end module variables
