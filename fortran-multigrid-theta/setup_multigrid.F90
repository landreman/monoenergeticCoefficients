#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

subroutine setup_multigrid()

  use petscmat
  use petscksp
  use variables

  implicit none

  integer :: j, k, level
  PC :: preconditioner_context
  Vec :: temp_vec
  PetscErrorCode :: ierr
  character(len=100) :: filename
  PetscViewer :: viewer

  if (masterProc) print *,"Entering setup_multigrid"

  call set_grid_resolutions()

  do j = 1,N_levels
     call create_grids(j)
  end do

  call computeB()

  ! Always build the high-order matrix at the finest level:
  call preallocateMatrix(levels(1)%high_order_matrix,1,1)
  call populateMatrix(levels(1)%high_order_matrix,1,1)

  if (defect_option==4) then
     call preallocateMatrix(levels(N_levels)%mixed_order_matrix,10,N_levels)
     call populateMatrix(levels(N_levels)%mixed_order_matrix,10,N_levels)
  end if

  do j = 1,N_levels
    
     if (defect_option > 1 .and. j>1) then
        ! Build the high-order matrix at this level:
        call preallocateMatrix(levels(j)%high_order_matrix,1,j)
        call populateMatrix(levels(j)%high_order_matrix,1,j)
     end if

     ! Build the low-order matrix at this level:
     call preallocateMatrix(levels(j)%low_order_matrix,0,j)
     call populateMatrix(levels(j)%low_order_matrix,0,j)

     ! Initialize some vectors on this level for the multigrid iterations:
     call MatCreateVecs(levels(j)%low_order_matrix, PETSC_NULL_OBJECT, levels(j)%residual_vec, ierr)
     call VecDuplicate(levels(j)%residual_vec, levels(j)%solution_vec, ierr)
     call VecDuplicate(levels(j)%residual_vec, levels(j)%temp_vec, ierr)
     call VecDuplicate(levels(j)%residual_vec, levels(j)%rhs_vec, ierr)
     call VecDuplicate(levels(j)%residual_vec, levels(j)%smoother_shift, ierr)
  end do

  matrixSize = levels(1)%matrixSize

  allocate(multigrid_prolongation_matrices(N_levels-1))
  allocate(multigrid_restriction_matrices(N_levels-1))

  do j=1,N_levels-1
     call restriction_prolongation_matrices(j)
  end do

  ! *****************************************************
  ! Build matrices and vectors needed for smoothing
  ! *****************************************************

  select case (smoothing_option)
  case (1)
     do level = 1,N_levels-1  ! We don't need to smooth on the coarsest level, where we do a direct solve.

        ! Make a matrix which is just like the low-order matrix, but without the sources and constraints:
        call preallocateMatrix(levels(level)%Jacobi_iteration_matrix,4,level)
        call populateMatrix(levels(level)%Jacobi_iteration_matrix,4,level)

        ! Extract the diagonal:        
        call MatCreateVecs(levels(level)%Jacobi_iteration_matrix, levels(level)%omega_times_inverse_diagonal,PETSC_NULL_OBJECT,ierr)
        call MatGetDiagonal(levels(level)%Jacobi_iteration_matrix, levels(level)%omega_times_inverse_diagonal,ierr)
        call VecReciprocal(levels(level)%omega_times_inverse_diagonal,ierr)
        call VecScale(levels(level)%omega_times_inverse_diagonal, Jacobi_omega, ierr)

        call VecDuplicate(levels(level)%omega_times_inverse_diagonal, temp_vec, ierr)
        call VecCopy(levels(level)%omega_times_inverse_diagonal, temp_vec, ierr)
        call VecScale(temp_vec, -one, ierr) ! So now temp_vec = -omega * (D^{-1})
        call MatDiagonalScale(levels(level)%Jacobi_iteration_matrix, temp_vec, PETSC_NULL_OBJECT, ierr)

        ! At this point Jacobi_iteration_matrix is -omega * (D^{-1}) * (L+U+D). We need to change the diagonal elements to 1-omega:
        call VecSet(temp_vec, one-Jacobi_omega, ierr)
        !call VecSetValue(temp_vec, levels(level)%matrixSize-1, one, INSERT_VALUES, ierr) ! Sets the last element to 1.0
        call VecSetValue(temp_vec, levels(level)%matrixSize-1, zero, INSERT_VALUES, ierr) ! Sets the last element to 0.0
        call VecAssemblyBegin(temp_vec, ierr)
        call VecAssemblyEnd(temp_vec, ierr)
        call MatDiagonalSet(levels(level)%Jacobi_iteration_matrix, temp_vec, INSERT_VALUES, ierr)
        call VecDestroy(temp_vec)

        ! Set last element of omega_times_inverse_diagonal to be 0, for use in the smoothing iterations:
        call VecSetValue(levels(level)%omega_times_inverse_diagonal, levels(level)%matrixSize-1, zero, INSERT_VALUES, ierr) ! Sets the last element to 0.0
        call VecAssemblyBegin(levels(level)%omega_times_inverse_diagonal, ierr)
        call VecAssemblyEnd(levels(level)%omega_times_inverse_diagonal, ierr)

        write (filename,fmt="(a,i1,a)") "mmc_Jacobi_iteration_matrix_level_",level,".dat"
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), FILE_MODE_WRITE, viewer, ierr)
        call MatView(levels(level)%Jacobi_iteration_matrix, viewer, ierr)
        call PetscViewerDestroy(viewer, ierr)
     end do

  case (2)
     do level = 1,N_levels-1  ! We don't need to smooth on the coarsest level, where we do a direct solve.

        ! [(1-omega)I + omega D] u = omega b = [(1-omega)I - omega(L+U)] u 

        ! Make matrices which are just like the low-order matrix, but without the sources and constraints:
        ! Diagonal part, and off-diagonal in zeta:
        call preallocateMatrix(levels(level)%smoothing_diagonal_matrix,5,level)
        call populateMatrix(levels(level)%smoothing_diagonal_matrix,5,level)
        ! Off-diagonal terms in theta and xi:
        call preallocateMatrix(levels(level)%smoothing_off_diagonal_matrix,6,level)
        call populateMatrix(levels(level)%smoothing_off_diagonal_matrix,6,level)

        print *,"11111"
        ! Scale and shift diagonals:
        call MatScale(levels(level)%smoothing_diagonal_matrix, Jacobi_omega, ierr)
        !call MatCreateVecs(levels(level)%smoothing_diagonal_matrix, temp_vec, PETSC_NULL_OBJECT,ierr)
        print *,"@@@@@@@"
        !call VecSet(temp_vec, one-Jacobi_omega, ierr)
        print *,"2222222"
        !call MatDiagonalSet(levels(level)%smoothing_diagonal_matrix, temp_vec, ADD_VALUES, ierr)
        call MatShift(levels(level)%smoothing_diagonal_matrix, one-Jacobi_omega, ierr)
        print *,"333333"

        call MatScale(levels(level)%smoothing_off_diagonal_matrix, -Jacobi_omega, ierr)
        print *,"3.5  3.5   3.5"
        !call MatDiagonalSet(levels(level)%smoothing_off_diagonal_matrix, temp_vec, ADD_VALUES, ierr)
        call MatShift(levels(level)%smoothing_off_diagonal_matrix, one-Jacobi_omega, ierr)
        print *,"3.7 3.7 3.7"
        !call VecDestroy(temp_vec)
        print *,"44444"
        call KSPCreate(PETSC_COMM_WORLD,levels(level)%smoothing_ksp,ierr)
        print *,"4.5 4.5 4.5"
        call KSPSetOperators(levels(level)%smoothing_ksp, levels(level)%smoothing_diagonal_matrix, levels(level)%smoothing_diagonal_matrix, ierr)
        print *,"55555"
        call KSPGetPC(levels(level)%smoothing_ksp,preconditioner_context,ierr)
        print *,"5.5 5.5 5.5"
        call PCSetType(preconditioner_context,PCLU,ierr)
        call KSPSetType(levels(level)%smoothing_ksp, KSPPREONLY, ierr)
        call PCFactorSetMatSolverPackage(preconditioner_context, MATSOLVERMUMPS, ierr)

        print *,"666666"
        write (filename,fmt="(a,i1,a)") "mmc_smoothing_diagonal_matrix_level_",level,".dat"
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), FILE_MODE_WRITE, viewer, ierr)
        call MatView(levels(level)%smoothing_diagonal_matrix, viewer, ierr)
        call PetscViewerDestroy(viewer, ierr)

        write (filename,fmt="(a,i1,a)") "mmc_smoothing_off_diagonal_matrix_level_",level,".dat"
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), FILE_MODE_WRITE, viewer, ierr)
        call MatView(levels(level)%smoothing_off_diagonal_matrix, viewer, ierr)
        call PetscViewerDestroy(viewer, ierr)
     end do

  case (3)
     ! Same as PETSC's native multigrid: KSPRICHARDSON with PCSOR

     do level = 1,N_levels-1  ! We don't need to smooth on the coarsest level, where we do a direct solve.

        ! Make a matrix which is just like the low-order matrix, but without the sources and constraints:
        call preallocateMatrix(levels(level)%smoothing_off_diagonal_matrix,4,level)
        call populateMatrix(levels(level)%smoothing_off_diagonal_matrix,4,level)

        call KSPCreate(PETSC_COMM_WORLD,levels(level)%smoothing_ksp,ierr)
        call KSPSetOperators(levels(level)%smoothing_ksp, levels(level)%smoothing_off_diagonal_matrix, levels(level)%smoothing_off_diagonal_matrix, ierr)
        call KSPGetPC(levels(level)%smoothing_ksp,preconditioner_context,ierr)
        call PCSetType(preconditioner_context,PCSOR,ierr)
        call KSPSetType(levels(level)%smoothing_ksp, KSPRICHARDSON, ierr)
        call KSPSetTolerances(levels(level)%smoothing_ksp, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, 2, ierr) ! Is the correct number of iterations 2 or 1?
        call KSPSetNormType(levels(level)%smoothing_ksp, KSP_NORM_NONE, ierr)
        call KSPSetInitialGuessNonzero(levels(level)%smoothing_ksp, PETSC_TRUE, ierr)
     end do
  case default
     print *,"Invalid smoothing_option:",smoothing_option
  end select


!!$  ! We need the low-order preconditioner matrix on the coarsest level:
!!$  call preallocateMatrix(preconditioner_matrix_on_coarsest_level,0,N_levels)
!!$  call populateMatrix(preconditioner_matrix_on_coarsest_level,0,N_levels)

  ! Set up the direct solver for the coarsest level
  call KSPCreate(PETSC_COMM_WORLD,ksp_on_coarsest_level,ierr)
  if (defect_option==2) then
     if (masterProc) print *,"Using HIGH order matrix for the direct solve on the coarest level."
     call KSPSetOperators(ksp_on_coarsest_level, levels(N_levels)%high_order_matrix, levels(N_levels)%high_order_matrix, ierr)
  elseif (defect_option==4) then
     if (masterProc) print *,"Using MIXED discretization matrix for the direct solve on the coarest level."
     call KSPSetOperators(ksp_on_coarsest_level, levels(N_levels)%mixed_order_matrix, levels(N_levels)%mixed_order_matrix, ierr)
  else
     if (masterProc) print *,"Using LOW order matrix for the direct solve on the coarest level."
     call KSPSetOperators(ksp_on_coarsest_level, levels(N_levels)%low_order_matrix, levels(N_levels)%low_order_matrix, ierr)
  end if
  call KSPGetPC(ksp_on_coarsest_level,preconditioner_context,ierr)
  call PCSetType(preconditioner_context,PCLU,ierr)
  call KSPSetType(ksp_on_coarsest_level, KSPPREONLY, ierr)
  call PCFactorSetMatSolverPackage(preconditioner_context, MATSOLVERMUMPS, ierr)

!  stop

end subroutine setup_multigrid


