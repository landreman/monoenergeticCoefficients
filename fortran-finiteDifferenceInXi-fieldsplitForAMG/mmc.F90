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
  PC :: pc
  IS :: IS_main, IS_source_constraint
  KSP :: sub_ksps(2)
  Mat :: sub_Amat, sub_Pmat
  integer, dimension(:), allocatable :: is_array
  integer :: j, num_fieldsplits
  MatNullSpace :: nullspace
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
  E = 0
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

  theta_derivative_option = 10
  preconditioner_theta_derivative_option = 4
  zeta_derivative_option = 8
  preconditioner_zeta_derivative_option = 4
  xi_derivative_option = 3
  preconditioner_xi_derivative_option = 2
  pitch_angle_scattering_option = 3
  preconditioner_pitch_angle_scattering_option = 2
  theta_upwinding_factor = 0.2
  zeta_upwinding_factor = 0.0

  xi_quadrature_option = 3
  constraint_option = 1
  do_fieldsplit = .false.
  geometry_option = 1

  call read_input()
  ! Command-line arguments will override input.namelist

#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR > 6))
#define new_argument PETSC_NULL_OBJECT,
#else
#define new_argument
#endif

  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Ntheta', Ntheta, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Nzeta', Nzeta, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Nxi', Nxi, wasSet, ierr)
  call PetscOptionsGetReal(new_argument PETSC_NULL_CHARACTER, '-nu', nu, wasSet, ierr)
  call PetscOptionsGetReal(new_argument PETSC_NULL_CHARACTER, '-E', E, wasSet, ierr)
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
  call PetscOptionsGetBool(new_argument PETSC_NULL_CHARACTER, '-do_fieldsplit', do_fieldsplit, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-geometry_option', geometry_option, wasSet, ierr)
  if (geometry_option==2) Nperiods=5

  if (masterProc) then
     print *,"Ntheta = ",Ntheta
     print *,"Nzeta = ",Nzeta
     print *,"Nxi = ",Nxi
     print *,"nu = ",nu
     print *,"E = ",E
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
     print *,"do_fieldsplit = ",do_fieldsplit
     print *,"geometry_option = ",geometry_option
  end if

  if (do_fieldsplit .and. constraint_option .ne. 1) stop "do_fieldsplit requires constraint_option=1."

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

  if (do_fieldsplit) then
     call KSPGetPC(ksp, pc, ierr)
     call PCSetType(pc, PCFIELDSPLIT, ierr)

     allocate(is_array(matrixSize-1))
     do j=1,matrixSize-1
        is_array(j)=j-1
     end do
     call ISCreateGeneral(PETSC_COMM_WORLD,matrixSize-1,is_array,PETSC_COPY_VALUES,IS_main,ierr)
     !call ISView(IS_main,PETSC_VIEWER_STDOUT_WORLD,ierr)
     call PCFieldSplitSetIS(pc,PETSC_NULL_CHARACTER,IS_main,ierr)
     
     is_array(1) = matrixSize-1
     call ISCreateGeneral(PETSC_COMM_WORLD,1,is_array,PETSC_COPY_VALUES,IS_source_constraint,ierr)
     !call ISView(IS_source_constraint,PETSC_VIEWER_STDOUT_WORLD,ierr)
     call PCFieldSplitSetIS(pc,PETSC_NULL_CHARACTER,IS_source_constraint,ierr)
  end if

!  call KSPSetComputeRHS(ksp,populateRHS,userContext,ierr)
!  call KSPSetComputeOperators(ksp,populateMatrix,userContext,ierr)
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 7))
  call KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
#else
  call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, vf, ierr) 
  !call KSPMonitorSet(ksp, KSPMonitorDefault, vf, PetscViewerAndFormatDestroy, ierr)
  call KSPMonitorSet(ksp, KSPMonitorTrueResidualNorm, vf, PetscViewerAndFormatDestroy, ierr)
#endif
  call KSPSetFromOptions(ksp,ierr)

  if (do_fieldsplit) then
     ! Add a null space to the main block
     call KSPSetUp(ksp,ierr)
     call PCFieldSplitGetSubKSP(pc,num_fieldsplits, sub_ksps, ierr)
     call KSPGetOperators(sub_ksps(1), sub_Amat, sub_Pmat, ierr)
     print *,"Does sub_Amat==sub_Pmat?",sub_Amat==sub_Pmat
     call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
     call MatSetNullSpace(sub_Pmat,nullspace,ierr)
     !call KSPSetNullSpace(sub_ksps(1),nullspace,ierr)
     call MatNullSpaceDestroy(nullspace,ierr)
  end if

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
