! Main program

#include <petsc/finclude/petsckspdef.h>

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
  PC :: pc
  Vec :: rhs, solution
  PetscInt :: userContext ! Not used
  Mat :: matrix, preconditioner_matrix
  PetscViewerAndFormat :: vf 

  external populateMatrix, populateRHS

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
  theta_derivative_option = 3
  zeta_derivative_option  = 3
  preconditioner_theta_derivative_option = 3
  preconditioner_zeta_derivative_option  = 3
  constraint_option = 1

  call readInput()

  ! Command-line arguments will override input.namelist:

  call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER, '-Ntheta', Ntheta, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER, '-Nzeta', Nzeta, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER, '-Nxi', Nxi, wasSet, ierr)
  call PetscOptionsGetReal(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER, '-nu', nu, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER, '-theta_derivative_option', theta_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER, '-zeta_derivative_option', zeta_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER, '-preconditioner_theta_derivative_option', preconditioner_theta_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER, '-preconditioner_zeta_derivative_option', preconditioner_zeta_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER, '-constraint_option', constraint_option, wasSet, ierr)

  if (masterProc) then
     print *,"Ntheta = ",Ntheta
     print *,"Nzeta = ",Nzeta
     print *,"Nxi = ",Nxi
     print *,"nu = ",nu
     print *,"theta_derivative_option    = ",theta_derivative_option
     print *,"zeta_derivative_option     = ",zeta_derivative_option
     print *,"preconditioner_theta_derivative_option    = ",preconditioner_theta_derivative_option
     print *,"preconditioner_zeta_derivative_option     = ",preconditioner_zeta_derivative_option
     print *,"constraint_option = ",constraint_option
  end if

  call createGrids()

  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)

  call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrixSize, rhs, ierr)
  call VecDuplicate(rhs, solution, ierr)

  call populateRHS(rhs)

  call preallocateMatrix(preconditioner_matrix,0)
  call populateMatrix(preconditioner_matrix,0)
  call preallocateMatrix(matrix,1)
  call populateMatrix(matrix,1)

  call KSPSetOperators(ksp, matrix, preconditioner_matrix, ierr)
  call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, vf, ierr) 
  !call KSPMonitorSet(ksp, KSPMonitorDefault, vf, PetscViewerAndFormatDestroy, ierr)
  call KSPMonitorSet(ksp, KSPMonitorTrueResidualNorm, vf, PetscViewerAndFormatDestroy, ierr)
  call KSPGetPC(ksp, pc, ierr)
  call PCSetType(pc, PCLU, ierr)
  call KSPSetFromOptions(ksp,ierr)

  if (masterProc) then
     print *,"Beginning solve..."
  end if
  call system_clock(clockStart, clockRate)
  call KSPSolve(ksp, rhs, solution, ierr)
  if (masterProc) then
     print *,"Done!"
  end if

  call diagnostics(solution)

  call KSPDestroy(ksp,ierr)

  call PETScFinalize(ierr)

end program
