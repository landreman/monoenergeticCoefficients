! Main program

#include <finclude/petsckspdef.h>

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
  Mat :: matrix

  external populateMatrix, populateRHS

  call PETSCInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_COMM_SIZE(PETSC_COMM_WORLD, numProcs, ierr)
  call MPI_COMM_RANK(PETSC_COMM_WORLD, myRank, ierr)

  print *,"I am proc",myRank," of ",numProcs
  masterProc = (myRank==0)

  ! Set defaults:
  nu = 0.1d+0
  diagonalShift = nu * 1.0d-5
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
  thetaGridScheme = 10
  zetaGridScheme = 10
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-Ntheta', Ntheta, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-Nzeta', Nzeta, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-Nxi', Nxi, wasSet, ierr)
  call PetscOptionsGetReal(PETSC_NULL_CHARACTER, '-nu', nu, wasSet, ierr)
  call PetscOptionsGetReal(PETSC_NULL_CHARACTER, '-diagonalShift', diagonalShift, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-thetaGridScheme', thetaGridScheme, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-zetaGridScheme', zetaGridScheme, wasSet, ierr)

  if (masterProc) then
     print *,"Ntheta = ",Ntheta
     print *,"Nzeta = ",Nzeta
     print *,"Nxi = ",Nxi
     print *,"nu = ",nu
     print *,"diagonalShift = ",diagonalShift
     print *,"thetaGridScheme    = ",thetaGridScheme
     print *,"zetaGridScheme     = ",zetaGridScheme
  end if

  call createGrids()

  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)

  call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrixSize, rhs, ierr)
  call VecDuplicate(rhs, solution, ierr)

  call populateRHS(rhs)

  call preallocateMatrix(matrix)
  call populateMatrix(matrix)

  call KSPSetOperators(ksp, matrix, matrix, ierr)

!  call KSPSetComputeRHS(ksp,populateRHS,userContext,ierr)
!  call KSPSetComputeOperators(ksp,populateMatrix,userContext,ierr)
  call KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
  call KSPGetPC(ksp, pc, ierr)
  call PCSetType(pc, PCLU, ierr)
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
