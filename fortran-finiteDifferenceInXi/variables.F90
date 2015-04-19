
module variables

  implicit none

!#include <finclude/petscsys.h>
#include <finclude/petscsysdef.h>

  PetscInt :: Ntheta, Nzeta, Nxi, Nperiods, helicity_l, matrixSize
  PetscReal :: nu, epsilon_t, epsilon_h, iota, G, I
  PetscOffset :: offset
  PetscInt :: myRank, numProcs
  logical :: masterProc

  PetscScalar, dimension(:), allocatable :: theta, zeta, xi
  PetscScalar, dimension(:), allocatable :: thetaWeights, zetaWeights, xiWeights
  PetscScalar, dimension(:,:), allocatable :: ddtheta_positiveXi
  PetscScalar, dimension(:,:), allocatable :: ddtheta_negativeXi
  PetscScalar, dimension(:,:), allocatable :: ddtheta_preconditioner_positiveXi
  PetscScalar, dimension(:,:), allocatable :: ddtheta_preconditioner_negativeXi
  PetscScalar, dimension(:,:), allocatable :: ddzeta_positiveXi
  PetscScalar, dimension(:,:), allocatable :: ddzeta_negativeXi
  PetscScalar, dimension(:,:), allocatable :: ddzeta_preconditioner_positiveXi
  PetscScalar, dimension(:,:), allocatable :: ddzeta_preconditioner_negativeXi
  PetscScalar, dimension(:,:), allocatable :: ddxi, LorentzOperator
  PetscScalar, dimension(:,:), allocatable :: B, dBdtheta, dBdzeta

  PetscScalar, parameter :: pi = 3.14159265358979d+0
  PetscScalar, parameter :: sqrtpi = 1.77245385090552d+0

  PetscInt :: thetaGridScheme, zetaGridScheme, thetaGridScheme_pc, zetaGridScheme_pc

end module variables
