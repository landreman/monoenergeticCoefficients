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
  KSP :: outer_ksp
  PC :: outer_preconditioner, inner_preconditioner
  Vec :: rhs, solution
  PetscInt :: userContext ! Not used
  Mat :: matrix, preconditionerMatrix
  PetscViewerAndFormat :: vf
  PetscViewer :: viewer
  integer, dimension(:), allocatable :: is_array
  IS :: IS_f0F1, IS_source_constraint
  integer :: j

  external apply_preconditioner

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
  L_scaling_option = 2
  fieldsplit_option = 0

  print *,"AAA"
  call readInput()
  print *,"BBB"
  ! Command-line arguments will override input.namelist:

  call PetscOptionsInsertString(PETSC_NULL_OBJECT, "-options_left", ierr) ! This helps for spotting misspelled options!
  call PetscOptionsInsertString(PETSC_NULL_OBJECT, "-outer_ksp_view", ierr)
  print *,"CCC"

  call PetscOptionsGetInt(PETSC_NULL_OBJECT, PETSC_NULL_CHARACTER, '-Ntheta', Ntheta, wasSet, ierr)
  print *,"DDD"
  call PetscOptionsGetInt(PETSC_NULL_OBJECT, PETSC_NULL_CHARACTER, '-Nzeta', Nzeta, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OBJECT, PETSC_NULL_CHARACTER, '-Nxi', Nxi, wasSet, ierr)
  call PetscOptionsGetReal(PETSC_NULL_OBJECT, PETSC_NULL_CHARACTER, '-nu', nu, wasSet, ierr)
  call PetscOptionsGetReal(PETSC_NULL_OBJECT, PETSC_NULL_CHARACTER, '-diagonalShift', diagonalShift, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OBJECT, PETSC_NULL_CHARACTER, '-thetaGridScheme', thetaGridScheme, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OBJECT, PETSC_NULL_CHARACTER, '-zetaGridScheme', zetaGridScheme, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OBJECT, PETSC_NULL_CHARACTER, '-L_scaling_option', L_scaling_option, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OBJECT, PETSC_NULL_CHARACTER, '-fieldsplit_option', fieldsplit_option, wasSet, ierr)

  if (masterProc) then
     print *,"Ntheta = ",Ntheta
     print *,"Nzeta = ",Nzeta
     print *,"Nxi = ",Nxi
     print *,"nu = ",nu
     print *,"diagonalShift = ",diagonalShift
     print *,"thetaGridScheme    = ",thetaGridScheme
     print *,"zetaGridScheme     = ",zetaGridScheme
     print *,"L_scaling_option = ",L_scaling_option
  end if

  call createGrids()

  call KSPCreate(PETSC_COMM_WORLD,outer_ksp,ierr)
  call KSPAppendOptionsPrefix(outer_ksp, 'outer_', ierr)
  call KSPGetPC(outer_ksp, outer_preconditioner, ierr)
  call PCSetType(outer_preconditioner, PCSHELL, ierr)
  call PCShellSetApply(outer_preconditioner, apply_preconditioner, ierr)

  call KSPCreate(PETSC_COMM_WORLD,inner_ksp,ierr)
  call KSPAppendOptionsPrefix(inner_ksp, 'inner_', ierr)
  call KSPGetPC(inner_ksp, inner_preconditioner, ierr)
  call KSPSetType(inner_ksp, KSPPREONLY, ierr)

  call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrixSize, rhs, ierr)
  call VecDuplicate(rhs, solution, ierr)

  call populateMatrix(matrix,1)
  call populateMatrix(preconditionerMatrix,0)
  call populateRHS(rhs)

  print *,"mmc MMMM"
  call KSPSetOperators(outer_ksp, matrix, preconditionerMatrix, ierr)
  call KSPSetOperators(inner_ksp, matrix, preconditionerMatrix, ierr)
  print *,"mmc NNNN"
!  call KSPSetComputeRHS(ksp,populateRHS,userContext,ierr)
!  call KSPSetComputeOperators(ksp,populateMatrix,userContext,ierr)

!  call KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr) !Old syntax
  call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, vf, ierr)
  !call KSPMonitorSet(ksp, KSPMonitorDefault, vf, PetscViewerAndFormatDestroy, ierr)
  call KSPMonitorSet(outer_ksp, KSPMonitorTrueResidualNorm, vf, PetscViewerAndFormatDestroy, ierr)

  if (fieldsplit_option == 0) then
     ! No fieldsplit
     call PCSetType(inner_preconditioner, PCLU, ierr)
     call PCFactorSetMatSolverPackage(inner_preconditioner, MATSOLVERMUMPS, ierr)
     print *,"NOT using a fieldsplit to separate off the sources/constraints."
  else
     ! Set up a fieldsplit to separate off the sources/constraints
     print *,"Setting up a fieldsplit to separate off the sources/constraints."

     call PCSetType(inner_preconditioner, PCFIELDSPLIT, ierr)
     call PCFIELDSPLITSetType(inner_preconditioner, PC_COMPOSITE_ADDITIVE, ierr) ! This line sets block Jacobi.

     ! Create an index set 'IS_all' which represents all indices of the big matrix and vectors, which each processor
     ! owning the indices it usually owns.
     !allocate(IS_array(matrixSize))
     !call VecGetOwnershipRange(solution, first_row_this_proc_owns, last_row_this_proc_owns, ierr)
     !IS_array(1:last_row_this_proc_owns-first_row_this_proc_owns) = [( j, j = first_row_this_proc_owns, last_row_this_proc_owns-1 )]
     !call ISCreateGeneral(PETSC_COMM_WORLD,last_row_this_proc_owns-first_row_this_proc_owns,IS_array,PETSC_COPY_VALUES,IS_all,ierr)
     
     ! The next 2 lines only work in serial.
     IS_array = [( j, j=0,matrixSize-3 )]
     call ISCreateGeneral(PETSC_COMM_WORLD,matrixSize-2,IS_array,PETSC_COPY_VALUES,IS_f0F1,ierr)
     
     ! Create an index set 'IS_source_constraint' that represents the indices for the sources and constraints.
     ! In this index set, the master processor owns everything, unlike the global matrix and vectors in which
     ! one or more processors at the end of the communicator own the sources/constraints.
     IS_array(1) = matrixSize-2 ! Remember -1 since PETSc uses 0-based indices
     IS_array(2) = matrixSize-1 ! Remember -1 since PETSc uses 0-based indices
     call ISCreateGeneral(PETSC_COMM_WORLD,2,IS_array,PETSC_COPY_VALUES,IS_source_constraint,ierr)

     call PCFieldSplitSetIS(inner_preconditioner,'f0F1',IS_f0F1,ierr)
     call PCFieldSplitSetIS(inner_preconditioner,'constraints',IS_source_constraint,ierr)

  end if
  call KSPSetFromOptions(outer_ksp,ierr)
  call KSPSetFromOptions(inner_ksp,ierr)

  call KSPSetUp(inner_ksp, ierr)
  call KSPView(inner_ksp, PETSC_VIEWER_STDOUT_WORLD,ierr)

  if (masterProc) then
     print *,"Beginning solve..."
  end if
  call system_clock(clockStart, clockRate)
  call KSPSolve(outer_ksp, rhs, solution, ierr)
  if (masterProc) then
     print *,"Done!"
  end if

  !call VecView(solution, PETSC_VIEWER_STDOUT_WORLD,ierr)
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mmc_solution.dat", FILE_MODE_WRITE, viewer, ierr)
  call VecView(solution, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)

  call diagnostics(solution)

  call KSPDestroy(outer_ksp,ierr)

  call PETScFinalize(ierr)

end program
