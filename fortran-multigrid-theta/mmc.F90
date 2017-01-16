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
  PC :: preconditioner_context
  Vec :: rhs, solution
  PetscInt :: userContext ! Not used
  Mat :: high_order_matrix_on_finest_level
  PetscInt :: VecLocalSize
#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR > 6))
  PetscViewerAndFormat vf
#endif

  !external populateMatrix, populateRHS
  external apply_multigrid_cycle

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

  refine_theta = PETSC_FALSE
  refine_zeta = PETSC_TRUE
  refine_xi = PETSC_FALSE
  Ntheta_min = 5
  Nzeta_min = 7
  Nxi_min = 9

  smoothing_option = 1
  restriction_option = 1

#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR > 6))
#define new_argument PETSC_NULL_OBJECT,
#else
#define new_argument
#endif

  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Ntheta', Ntheta, wasSet, ierr)
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
  call PetscOptionsGetBool(new_argument PETSC_NULL_CHARACTER, '-refine_theta', refine_theta, wasSet, ierr)
  call PetscOptionsGetBool(new_argument PETSC_NULL_CHARACTER, '-refine_zeta', refine_zeta, wasSet, ierr)
  call PetscOptionsGetBool(new_argument PETSC_NULL_CHARACTER, '-refine_xi', refine_xi, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Ntheta_min', Ntheta_min, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Nzeta_min', Nzeta_min, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Nxi_min', Nxi_min, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-smoothing_option', smoothing_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-restriction_option', restriction_option, wasSet, ierr)

  ! Make sure Ntheta and Nzeta are odd:
  if (mod(Ntheta, 2) == 0) then
     Ntheta = Ntheta + 1
  end if
  if (mod(Nzeta, 2) == 0) then
     Nzeta = Nzeta + 1
  end if

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
     print *,"refine_theta = ",refine_theta
     print *,"refine_zeta = ",refine_zeta
     print *,"refine_xi = ",refine_xi
     print *,"Ntheta_min = ",Ntheta_min
     print *,"Nzeta_min = ",Nzeta_min
     print *,"Nxi_min = ",Nxi_min
     print *,"smoothing_option = ",smoothing_option
     print *,"restriction_option = ",restriction_option
  end if

  if (constraint_option<0 .or. constraint_option>2) stop "Invalid constraint_option"

  call setup_multigrid()

  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)

  call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrixSize, rhs, ierr)
  call VecDuplicate(rhs, solution, ierr)

  call populateRHS(rhs)

  ! We need the high-order 'real' matrix on the finest level:
  call preallocateMatrix(high_order_matrix_on_finest_level,1,1)
  call populateMatrix(high_order_matrix_on_finest_level,1,1)
  ! In the next line, the 3rd argument (the Mat object for the KSP preconditioner)
  ! is not used for anything because we set a custom preconditioner.
  call KSPSetOperators(ksp, high_order_matrix_on_finest_level, high_order_matrix_on_finest_level, ierr)

  ! Create the object for the multigrid preconditioner:
  call KSPGetPC(ksp,preconditioner_context,ierr)
  call PCSetType(preconditioner_context,PCSHELL,ierr)
  call PCShellSetApply(preconditioner_context,apply_multigrid_cycle,ierr)


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
