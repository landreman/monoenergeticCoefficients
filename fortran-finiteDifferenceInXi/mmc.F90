! Main program

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

program mmc

  use petscksp

  use variables

  implicit none

!#include <finclude/petscsys.h>
!#include <finclude/petscvec.h>
!#include <finclude/petscmat.h>
!#include <finclude/petscksp.h>

  PetscErrorCode ierr
  PetscBool :: wasSet
  KSP :: ksp
  Vec :: rhs, solution
  PetscInt :: userContext ! Not used
  Mat :: matrix, pcMatrix
#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR > 6))
  PetscViewerAndFormat vf
#endif

  !external populateMatrix, populateRHS

  call PETSCInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_COMM_SIZE(PETSC_COMM_WORLD, numProcs, ierr)
  call MPI_COMM_RANK(PETSC_COMM_WORLD, myRank, ierr)

  print *,"I am proc",myRank," of ",numProcs
  masterProc = (myRank==0)

  ! Set defaults:
  nu = 0.1d+0
  epsilon_t = -0.07053d+0
  epsilon_h = 0.05067d+0
  iota = 0.4542d+0
  G = 3.7481d+0
  I = 0d+0
  Nperiods = 10
  helicity_l = 2
  Ntheta = 13
  Nzeta = 15
  Nxi = 16

  theta_derivative_option = 8
  preconditioner_theta_derivative_option = 4
  zeta_derivative_option = 8
  preconditioner_zeta_derivative_option = 4
  xi_derivative_option = 3
  preconditioner_xi_derivative_option = 2
  pitch_angle_scattering_option = 3
  preconditioner_pitch_angle_scattering_option = 2

  xi_quadrature_option = 3
  constraint_option = 2

  print *,"AAA"
#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR > 6))
#define new_argument PETSC_NULL_OBJECT,
#else
#define new_argument
#endif

  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Ntheta', Ntheta, wasSet, ierr)
  print *,"BBB"
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Nzeta', Nzeta, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Nxi', Nxi, wasSet, ierr)
  call PetscOptionsGetReal(new_argument PETSC_NULL_CHARACTER, '-nu', nu, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-theta_derivative_option', theta_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-preconditioner_theta_derivative_option', preconditioner_theta_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-zeta_derivative_option', zeta_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-preconditioner_zeta_derivative_option', preconditioner_zeta_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-xi_derivative_option', xi_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-preconditioner_xi_derivative_option', preconditioner_xi_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-pitch_angle_scattering_option', pitch_angle_scattering_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-preconditioner_pitch_angle_scattering_option', preconditioner_pitch_angle_scattering_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-xi_quadrature_option', xi_quadrature_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-constraint_option', constraint_option, wasSet, ierr)

  if (masterProc) then
     print *,"Ntheta = ",Ntheta
     print *,"Nzeta = ",Nzeta
     print *,"Nxi = ",Nxi
     print *,"nu = ",nu
     print *,"theta_derivative_option = ",theta_derivative_option
     print *,"preconditioner_theta_derivative_option = ",preconditioner_theta_derivative_option
     print *,"zeta_derivative_option = ",zeta_derivative_option
     print *,"preconditioner_zeta_derivative_option = ",preconditioner_zeta_derivative_option
     print *,"xi_derivative_option = ",xi_derivative_option
     print *,"preconditioner_xi_derivative_option = ",preconditioner_xi_derivative_option
     print *,"pitch_angle_scattering_option = ",pitch_angle_scattering_option
     print *,"preconditioner_pitch_angle_scattering_option = ",preconditioner_pitch_angle_scattering_option
     print *,"xi_quadrature_option = ",xi_quadrature_option
     print *,"constraint_option = ",constraint_option
  end if

  if (constraint_option<0 .or. constraint_option>2) stop "Invalid constraint_option"

  call createGrids()

  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)

  call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrixSize, rhs, ierr)
  call VecDuplicate(rhs, solution, ierr)

  call populateRHS(rhs)

  call preallocateMatrix(matrix,1)
  call preallocateMatrix(pcMatrix,0)

  call populateMatrix(matrix,1)
  call populateMatrix(pcMatrix,0)

  call KSPSetOperators(ksp, matrix, pcMatrix, ierr)

!  call KSPSetComputeRHS(ksp,populateRHS,userContext,ierr)
!  call KSPSetComputeOperators(ksp,populateMatrix,userContext,ierr)
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 7))
  call KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
#else
  call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, vf, ierr) 
  call KSPMonitorSet(ksp, KSPMonitorDefault, vf, PetscViewerAndFormatDestroy, ierr)
#endif
  call KSPSetFromOptions(ksp,ierr)

  if (masterProc) then
     print *,"Beginning solve..."
  end if
  call KSPSolve(ksp, rhs, solution, ierr)
  if (masterProc) then
     print *,"Done!"
  end if

  call diagnostics(solution)

  call KSPDestroy(ksp,ierr)

  call PETScFinalize(ierr)

end program
