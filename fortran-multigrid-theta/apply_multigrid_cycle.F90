#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

subroutine apply_multigrid_cycle(preconditioner_context, input_vec, output_vec, ierr)

  use petscsys
  use variables, only: levels, masterProc, N_levels, N_pre_smoothing, N_post_smoothing, one, ksp_on_coarsest_level, &
       multigrid_restriction_matrices, multigrid_prolongation_matrices, Jacobi_omega, smoothing_option, defect_option

  implicit none

  PC :: preconditioner_context
  Vec :: input_vec, output_vec
  PetscErrorCode :: ierr
  integer :: level, reason, j_smoothing, vec_size

  if (masterProc) print *,"  Beginning apply_multigrid_cycle"

  !call VecView(input_vec, PETSC_VIEWER_STDOUT_WORLD, ierr)

  call VecCopy(input_vec, levels(1)%rhs_vec, ierr)

  ! Proceed from the finest level to the coarsest level.
  do level = 1, N_levels-1
     if (masterProc) print "(a,i4)","     Applying pre-smoothing to level",level

!!$     print *,"111"
!!$     call VecGetSize(levels(level)%smoother_shift, vec_size, ierr)
!!$     print *,"1Size of levels(level)%smoother_shift:",vec_size
!!$     call VecGetSize(levels(level)%omega_times_inverse_diagonal, vec_size, ierr)
!!$     print *,"1Size of levels(level)%omega_times_inverse_diagonal:",vec_size
!!$     call VecGetSize(levels(level)%residual_vec, vec_size, ierr)
!!$     print *,"1Size of levels(level)%residual_vec:",vec_size
!!$     call VecGetSize(levels(level)%rhs_vec, vec_size, ierr)
!!$     print *,"1Size of levels(level)%rhs_vec:",vec_size
!!$     call VecGetSize(levels(level)%solution_vec, vec_size, ierr)
!!$     print *,"1Size of levels(level)%solution_vec:",vec_size
!!$     call VecGetSize(levels(level)%temp_vec, vec_size, ierr)
!!$     print *,"1Size of levels(level)%temp_vec:",vec_size

     select case (smoothing_option)
     case (1)
        ! smoother_shift = omega * (D^{-1}) * rhs_vec
        call VecPointwiseMult(levels(level)%smoother_shift, levels(level)%omega_times_inverse_diagonal, levels(level)%rhs_vec, ierr)
        ! Make sure last element is correct?

        ! Apply pre-smoothing. Note that since the initial guess is 0, then 1 iteration of smoothing is equivalent 
        ! to setting the solution vector to the smoother_shift.
        call VecCopy(levels(level)%smoother_shift, levels(level)%solution_vec, ierr)
        do j_smoothing = 2,N_pre_smoothing
           ! temp_vec = Jacobi_iteration_matrix * solution_vec + smoother_shift
           call MatMultAdd(levels(level)%Jacobi_iteration_matrix, levels(level)%solution_vec, levels(level)%smoother_shift, levels(level)%temp_vec, ierr)
           call VecCopy(levels(level)%temp_vec, levels(level)%solution_vec, ierr)
        end do

     case (2)
        ! [(1-omega)I + omega D] u = omega b = [(1-omega)I - omega(L+U)] u

        ! Set smoother_shift = omega * rhs_vec:
        call VecCopy(levels(level)%rhs_vec, levels(level)%smoother_shift, ierr)
        call VecScale(levels(level)%smoother_shift, Jacobi_omega, ierr)

        do j_smoothing = 1,N_pre_smoothing
           ! Form temp_vec = smoother_shift + smoothing_off_diagonal_matrix * solution_vec
           if (j_smoothing==1) then
              ! For first iteration, where the initial guess for the solution is 0, we can skip the matrix multiplication:
              call KSPSolve(levels(level)%smoothing_ksp, levels(level)%smoother_shift, levels(level)%solution_vec, ierr)
           else
              call MatMultAdd(levels(level)%smoothing_off_diagonal_matrix, levels(level)%solution_vec, levels(level)%smoother_shift, levels(level)%temp_vec, ierr)
              call KSPSolve(levels(level)%smoothing_ksp, levels(level)%temp_vec, levels(level)%solution_vec, ierr)
           end if
           call KSPGetConvergedReason(levels(level)%smoothing_ksp, reason, ierr)
           if (reason <= 0) then
              if (masterProc) print *,"KSP failed in post-smoothing! j_smoothing=",j_smoothing,", reason=",reason
              stop
           end if
        end do

     case (3)
        call VecSet(levels(level)%solution_vec, 0.0d+0, ierr)
        call KSPSolve(levels(level)%smoothing_ksp, levels(level)%rhs_vec, levels(level)%solution_vec, ierr)
        call KSPGetConvergedReason(levels(level)%smoothing_ksp, reason, ierr)
        if (reason <= 0) then
           if (masterProc) print *,"KSP failed in post-smoothing! j_smoothing=",j_smoothing,", reason=",reason
           stop
        end if

     case default
        stop "Invalid smoothing_option in apply_multigrid_cycle 1"
     end select

     ! Construct the residual = rhs_vec - matrix * solution_vec
     call VecCopy(levels(level)%rhs_vec, levels(level)%residual_vec, ierr)
     select case (defect_option)
     case (1)
        if (masterProc) print *,"    Computing residual using the LOW order matrix."
        call MatMult(levels(level)%low_order_matrix, levels(level)%solution_vec, levels(level)%temp_vec, ierr)
     case (2,3,4)
        if (masterProc) print *,"    Computing residual using the HIGH order matrix."
        call MatMult(levels(level)%high_order_matrix, levels(level)%solution_vec, levels(level)%temp_vec, ierr)
     case default
        print *,"Invalid defect_option:",defect_option
        stop
     end select
     call VecAXPY(levels(level)%residual_vec, -one, levels(level)%temp_vec, ierr)

     ! Restrict residual to the next level:
     call MatMult(multigrid_restriction_matrices(level), levels(level)%residual_vec, levels(level+1)%rhs_vec, ierr)
  end do

  ! Solve the system directly on the coarsest level.
  if (masterProc) print *,"    About to apply the direct solver on the coarsest level."
  call KSPSolve(ksp_on_coarsest_level, levels(N_levels)%rhs_vec, levels(N_levels)%solution_vec, ierr)
  call KSPGetConvergedReason(ksp_on_coarsest_level, reason, ierr)
  if (reason<=0) then
     if (masterProc) print *,"Error! KSP on coarsest level did not converge. KSPConvergedReason=",reason
  end if

  ! Proceed from the coarsest level to the finest level.
  do level = N_levels-1,1,-1
     if (masterProc) print "(a,i4)","     Applying post-smoothing to level",level
     ! Prolong solution from the previous level, and add result to the solution here.
     call MatMultAdd(multigrid_prolongation_matrices(level), levels(level+1)%solution_vec, levels(level)%solution_vec, levels(level)%temp_vec, ierr)
     call VecCopy(levels(level)%temp_vec, levels(level)%solution_vec, ierr)

     ! Apply post-smoothing
     select case (smoothing_option)
     case (1)
        do j_smoothing = 1,N_post_smoothing
           ! temp_vec = Jacobi_iteration_matrix * solution_vec + smoother_shift
           call MatMultAdd(levels(level)%Jacobi_iteration_matrix, levels(level)%solution_vec, levels(level)%smoother_shift, levels(level)%temp_vec, ierr)
           call VecCopy(levels(level)%temp_vec, levels(level)%solution_vec, ierr)
        end do
     case (2)
        do j_smoothing = 1,N_post_smoothing
           ! Form temp_vec = smoother_shift + smoothing_off_diagonal_matrix * solution_vec
           call MatMultAdd(levels(level)%smoothing_off_diagonal_matrix, levels(level)%solution_vec, levels(level)%smoother_shift, levels(level)%temp_vec, ierr)
           call KSPSolve(levels(level)%smoothing_ksp, levels(level)%temp_vec, levels(level)%solution_vec, ierr)
           call KSPGetConvergedReason(levels(level)%smoothing_ksp, reason, ierr)
           if (reason <= 0) then
              if (masterProc) print *,"KSP failed in post-smoothing! reason=",reason
              stop
           end if
        end do
     case (3)
        call KSPSolve(levels(level)%smoothing_ksp, levels(level)%temp_vec, levels(level)%solution_vec, ierr) ! Could also use solution_vec as both arguments here I think.
        call KSPGetConvergedReason(levels(level)%smoothing_ksp, reason, ierr)
        if (reason <= 0) then
           if (masterProc) print *,"KSP failed in post-smoothing! j_smoothing=",j_smoothing,", reason=",reason
           stop
        end if
     case default
        stop "Invalid smoothing_option in apply_multigrid_cycle 2"
     end select
  end do

  call VecCopy(levels(1)%solution_vec, output_vec, ierr)

end subroutine apply_multigrid_cycle


