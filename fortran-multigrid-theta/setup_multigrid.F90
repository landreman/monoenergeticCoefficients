#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

subroutine setup_multigrid()

  use petscmat
  use variables

  implicit none

  integer :: j, k, level, itheta, izeta
  PC :: preconditioner_context
  Vec :: temp_vec
  PetscErrorCode :: ierr

  call set_grid_resolutions()

  do j = 1,N_levels
     call create_grids(j, levels(j)%Ntheta, levels(j)%Nzeta, levels(j)%Nxi, levels(j)%matrixSize, levels(j)%theta, levels(j)%zeta, levels(j)%xi, &
          levels(j)%thetaWeights, levels(j)%zetaWeights, levels(j)%xiWeights, levels(j)%ithetaMin, levels(j)%ithetaMax, levels(j)%localNtheta, &
          levels(j)%izetaMin, levels(j)%izetaMax, levels(j)%localNzeta, &
          levels(j)%ddtheta_plus, levels(j)%ddtheta_minus, levels(j)%ddtheta_plus_preconditioner, levels(j)%ddtheta_minus_preconditioner, &
          levels(j)%ddzeta_plus, levels(j)%ddzeta_minus, levels(j)%ddzeta_plus_preconditioner, levels(j)%ddzeta_minus_preconditioner, &
          levels(j)%ddxi_plus, levels(j)%ddxi_minus, levels(j)%ddxi_plus_preconditioner, levels(j)%ddxi_minus_preconditioner, &
          levels(j)%pitch_angle_scattering_operator, levels(j)%pitch_angle_scattering_operator_preconditioner)

     
     ! Compute B and its derivatives     
     allocate(levels(j)%B(levels(j)%Ntheta,levels(j)%Nzeta))
     allocate(levels(j)%dBdtheta(levels(j)%Ntheta,levels(j)%Nzeta))
     allocate(levels(j)%dBdzeta(levels(j)%Ntheta,levels(j)%Nzeta))
     call computeB(levels(j)%Ntheta, levels(j)%theta, levels(j)%Nzeta, levels(j)%zeta, levels(j)%B, levels(j)%dBdtheta, levels(j)%dBdzeta)

     if (j==1) then
        VPrime = 0
        do itheta = 1,Ntheta
           do izeta = 1,Nzeta
              VPrime = VPrime + levels(1)%thetaWeights(itheta) * levels(1)%zetaWeights(izeta) / (levels(1)%B(itheta,izeta) ** 2)
           end do
        end do
     end if

     ! Build the low-order matrix at this level:
     call preallocateMatrix(levels(j)%low_order_matrix,0,j)
     call populateMatrix(levels(j)%low_order_matrix,0,j)

     ! Initialize some vectors on this level for the multigrid iterations:
     call MatCreateVecs(levels(j)%low_order_matrix, PETSC_NULL_OBJECT, levels(j)%residual_vec, ierr)
     call VecDuplicate(levels(j)%residual_vec, levels(j)%solution_vec, ierr)
     call VecDuplicate(levels(j)%residual_vec, levels(j)%temp_vec, ierr)
     call VecDuplicate(levels(j)%residual_vec, levels(j)%smoother_shift, ierr)
  end do

!!$  allocate(multigrid_interpolation_matrices(N_levels-1))
!!$  allocate(multigrid_restriction_matrices(N_levels-1))
!!$
!!$  do j=1,N_levels-1
!!$     call interpolation_restriction(multigrid_interpolation_matrices(j), multigrid_restriction_matrices(j), &
!!$          levels(j)%Ntheta, levels(j+1)%Ntheta, levels(j)%theta, levels(j+1)%theta, &
!!$          levels(j)%Nzeta, levels(j+1)%Nzeta, levels(j)%zeta, levels(j+1)%zeta, &
!!$          levels(j)%Nxi, levels(j+1)%Nxi, levels(j)%xi, levels(j+1)%xi)
!!$  end do

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
        call MatDiagonalSet(levels(level)%Jacobi_iteration_matrix, temp_vec, INSERT_VALUES, ierr)
        call VecDestroy(temp_vec)

        ! Set last element of omega_times_inverse_diagonal to be 0?
     end do
  case default
     print *,"Invalid smoothing_option:",smoothing_option
  end select


!!$  ! We need the low-order preconditioner matrix on the coarsest level:
!!$  call preallocateMatrix(preconditioner_matrix_on_coarsest_level,0,N_levels)
!!$  call populateMatrix(preconditioner_matrix_on_coarsest_level,0,N_levels)

  ! Set up the direct solver for the coarsest level
  call KSPCreate(PETSC_COMM_WORLD,ksp_on_coarsest_level,ierr)
  call KSPSetOperators(ksp_on_coarsest_level, levels(N_levels)%low_order_matrix, levels(N_levels)%low_order_matrix, ierr)
  call KSPGetPC(ksp_on_coarsest_level,preconditioner_context,ierr)
  call PCSetType(preconditioner_context,PCLU,ierr)
  call KSPSetType(ksp_on_coarsest_level, KSPPREONLY, ierr)
  call PCFactorSetMatSolverPackage(preconditioner_context, MATSOLVERMUMPS, ierr)

!  stop

end subroutine setup_multigrid


