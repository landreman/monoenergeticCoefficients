#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

subroutine apply_multigrid_cycle(preconditioner_context, input_vec, output_vec, ierr)

  use variables, only: levels, masterProc, N_levels, N_smoothing

  implicit none

  PC :: preconditioner_context
  Vec :: input_vec, output_vec
  PetscErrorCode :: ierr
  integer :: level, reason

  if (masterProc) print *,"Beginning apply_multigrid_cycle"

  call VecCopy(input_vec, levels(1)%rhs_vec, ierr)

  ! Proceed from the finest level to the coarsest level.
  do level = 1, N_levels-1
     ! smoother_shift = omega * (D^{-1}) * rhs_vec
     call VecPointwiseMult(levels(level)%smoother_shift, levels(level)%omega_times_inverse_diagonal, levels(level)%rhs_vec, ierr)
     ! Make sure last element is correct?

     ! Apply pre-smoothing. Note that since the initial guess is 0, then 1 iteration of smoothing is equivalent 
     ! to setting the solution vector to the smoother_shift.
     call VecCopy(levels(level)%smoother_shift, levels(level)%solution_vec, ierr)
     do j_smoothing = 2,N_smoothing
        ! temp_vec = Jacobi_iteration_matrix * solution_vec + smoother_shift
        call MatMultAdd(levels(level)%Jacobi_iteration_matrix, levels(level)%solution_vec, levels(level)%smoother_shift, levels(level)%temp_vec, ierr)
        call VecCopy(levels(level)%temp_vec, levels(level)%solution_vec, ierr)
     end do

     ! Construct the residual = rhs_vec - matrix * solution_vec
     call VecCopy(levels(level)%rhs_vec, levels(level)%residual_vec, ierr)
     call MatMult(levels(level)%temp_vec, levels(level)%low_order_matrix, levels(level)%solution_vec, ierr)
     call VecAXPY(levels(level)%residual_vec, -one, levels(level)%temp_vec, ierr)

     ! Restrict residual to the next level:
     call MatMult(multigrid_restriction_matrices(level), levels(level)%residual_vec, levels(level+1)%rhs_vec, ierr)
  end do

  ! Solve the system directly on the coarsest level.
  if (masterProc) print *,"About to apply the direct solver on the coarsest level."
  call KSPSolve(ksp_on_coarsest_level, levels(N_levels)%rhs_vec, levels(N_levels)%solution_vec, ierr)
  call KSPGetConvergedReason(ksp_on_coarsest_level, reason, ierr)
  if (reason<=0) then
     if (masterProc) print *,"Error! KSP on coarsest level did not converge. KSPConvergedReason=",reason
  end if

  ! Proceed from the coarsest level to the finest level.
  do level = 1,N_levels-1
     ! Prolong solution from the previous level, and add result to the solution there.
     call MatMultAdd(multigrid_prolongation_matrices(level), levels(level+1)%solution_vec, levels(level)%solution_vec, levels(level)%temp_vec, ierr)
     call VecCopy(levels(level)%temp_vec, levels(level)%solution_vec, ierr)

     ! Apply post-smoothing
     do j_smoothing = 1,N_smoothing
        ! temp_vec = Jacobi_iteration_matrix * solution_vec + smoother_shift
        call MatMultAdd(levels(level)%Jacobi_iteration_matrix, levels(level)%solution_vec, levels(level)%smoother_shift, levels(level)%temp_vec, ierr)
        call VecCopy(levels(level)%temp_vec, levels(level)%solution_vec, ierr)
     end do
  end do

  call VecCopy(levels(1)%solution_vec, output_vec, ierr)

end subroutine apply_multigrid_cycle


