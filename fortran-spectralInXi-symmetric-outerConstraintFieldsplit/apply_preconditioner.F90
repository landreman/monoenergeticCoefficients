  ! This subroutine could probably be sped up by re-using Vecs instead of creating/destroying them over and over again?
#include <petsc/finclude/petsckspdef.h>

  subroutine apply_preconditioner(outer_preconditioner, input_Vec, output_Vec, ierr)

    use petscksp
    use variables, only: masterProc, numProcs, myRank, inner_KSP, matrixSize, one, fieldsplit_option, Ntheta, Nzeta, Nxi, sqrt_g
    use indices

    implicit none

    PetscErrorCode :: ierr
    Mat :: matrix, preconditioner_matrix, temp_Mat
    PC :: outer_preconditioner
    Vec :: input_Vec, output_Vec
    KSPConvergedReason :: reason
    integer, dimension(:), allocatable :: IS_array
    integer :: IS_array_index, j, L, index, itheta, izeta
    ! For these next variables, the 'save' attribute means these variables can be initialized in the first pass
    ! through this subroutine, and used again on subsequent calls to this subroutine.
    KSP, save :: constraints_times_sources_KSP
    Mat, save :: sources_Mat, constraints_Mat, constraints_times_sources_Mat, pb_Mat, cq_Mat
    IS, save :: IS_all, IS_source_constraint
    logical, save :: first_call = .true.
    Vec :: temp_Vec_1, temp_Vec_2, ainv_times_stuff_Vec, cq_times_r_Vec, y_Vec, s_Vec, pq_Vec
    PC :: constraints_times_sources_PC
!    logical :: verbose = .true.
    logical :: verbose = .false.
    integer :: first_row_this_proc_owns, last_row_this_proc_owns
    PetscLogEvent :: event
    PetscViewer :: viewer

    call PetscLogEventRegister("apply_preconditi", 0, event, ierr)
    call PetscLogEventBegin(event,ierr)

    select case (fieldsplit_option)
       case (0)
          if (masterProc) print *,"apply_preconditioner called. Only applying inner_KSP."
          !print *,"Here comes KSPView on inner_KSP:"
          !call KSPView(inner_KSP, PETSC_VIEWER_STDOUT_WORLD,ierr)
          if (verbose) then
             print *,"000 Here comes input_Vec for the inner_KSP solve:"
             call VecView(input_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
             print *,"0.3 0.3 0.3 Here comes output_Vec before solve"
             call VecView(output_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
             print *,"0.7 0.7 0.7 about to call KSPSolve."
          end if
          call KSPSolve(inner_KSP, input_Vec, output_Vec, ierr)
          if (verbose) then
             print *,"0.9 0.9 0.9 Here comes output_Vec after solve"
             call VecView(output_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
             print *,"111"
          end if
          call KSPGetConvergedReason(inner_KSP, reason, ierr)
          print *,"222"
          if (reason <= 0) then
             print *,"WARNING: inner KSP failed with reason",reason
          end if

       case (1,2)
          if (masterProc) print *,"apply_preconditioner called. Handling sources/constraints."
          
          ! We will solve the block system
          !
          ! | a  b | |x|   |r|
          ! |      | | | = | |
          ! | c  0 | |y|   |s|
          !
          ! using the properties c q a = 0 and a p b = 0, where q and p are diagonal. 
          ! The solution is found in my note 20170703-01. The solution turns out to be
          ! x = (I - p b (c p b)^{-1} c) a^{-1} (I - b (c q b)^{-1} c q) r + p b (c p b)^{-1} s,
          ! y = (c q b)^{-1} c q r.
          ! For our problem, p=q has the action of applying 1/sqrt_g.

          ! First, if we haven't already done so, extract the sub-matrices b and c.
          if (first_call) then
             if (masterProc) print *,"Extracting sources and constraints for the shell preconditioner."

             ! Create an index set 'IS_all' which represents all indices of the big matrix and vectors, which each processor
             ! owning the indices it usually owns.
             allocate(IS_array(matrixSize))
             call VecGetOwnershipRange(input_Vec, first_row_this_proc_owns, last_row_this_proc_owns, ierr)
             IS_array(1:last_row_this_proc_owns-first_row_this_proc_owns) = [( j, j = first_row_this_proc_owns, last_row_this_proc_owns-1 )]
             call ISCreateGeneral(PETSC_COMM_WORLD,last_row_this_proc_owns-first_row_this_proc_owns,IS_array,PETSC_COPY_VALUES,IS_all,ierr)
             !print *,"Here comes IS_all:"
             !call ISView(IS_all, PETSC_VIEWER_STDOUT_WORLD,ierr)
             ! The next 2 lines only work in serial.
             !IS_array = [( j, j=0,matrixSize-1 )]
             !call ISCreateGeneral(PETSC_COMM_WORLD,matrixSize,IS_array,PETSC_COPY_VALUES,IS_all,ierr)
          
             ! Create an index set 'IS_source_constraint' that represents the indices for the sources and constraints.
             ! In this index set, the master processor owns everything, unlike the global matrix and vectors in which
             ! one or more processors at the end of the communicator own the sources/constraints.
             IS_array(1) = matrixSize-2 ! Remember -1 since PETSc uses 0-based indices
             IS_array(2) = matrixSize-1 ! Remember -1 since PETSc uses 0-based indices
             call ISCreateGeneral(PETSC_COMM_WORLD,2,IS_array,PETSC_COPY_VALUES,IS_source_constraint,ierr)
             !print *,"Here comes IS_source_constraint:"
             !call ISView(IS_source_constraint, PETSC_VIEWER_STDOUT_WORLD,ierr)

             ! Now 'slice' the big matrix to get the smaller non-square matrices that represent the sources and constraints:
             ! All rows, last columns -> sources_Mat
             call PCGetOperators(outer_preconditioner, matrix, preconditioner_matrix, ierr)
             call MatGetSubMatrix(preconditioner_matrix, IS_all, IS_source_constraint, MAT_INITIAL_MATRIX, sources_Mat, ierr)
             ! Last rows, all columns -> constraints_Mat
             call MatGetSubMatrix(preconditioner_matrix, IS_source_constraint, IS_all, MAT_INITIAL_MATRIX, constraints_Mat, ierr)

             ! Create the diagonal matrix p=q used for scaling.
             call VecDuplicate(input_vec, pq_Vec, ierr)
             call VecSet(pq_Vec, one, ierr)
             L=0
             do itheta = 1,Ntheta
                do izeta = 1,Nzeta
                   index = getIndex(itheta, izeta, L+1)
                   call VecSetValue(pq_Vec, index, one/sqrt_g(itheta,izeta), INSERT_VALUES, ierr)
                   call VecSetValue(pq_Vec, index + Ntheta*Nzeta*Nxi, one/sqrt_g(itheta,izeta), INSERT_VALUES, ierr)
                end do
             end do
             call VecAssemblyBegin(pq_Vec, ierr)
             call VecAssemblyEnd(pq_Vec, ierr)

             ! Set up the matrices p*b and c*q, which will be used later:
             call MatConvert(sources_Mat, MATSAME, MAT_INITIAL_MATRIX, pb_Mat, ierr) ! Creates a new matrix and copies an existing matrix into it.
             call MatConvert(constraints_Mat, MATSAME, MAT_INITIAL_MATRIX, cq_Mat, ierr) ! Creates a new matrix and copies an existing matrix into it.
             call MatDiagonalScale(pb_Mat, pq_Vec, PETSC_NULL_OBJECT, ierr) ! Left-multiply by p=q
             call MatDiagonalScale(cq_Mat, PETSC_NULL_OBJECT, pq_Vec, ierr) ! Right-multiply by p=q
             call VecDestroy(pq_Vec, ierr)

             ! Verify that c*q*a=0:
             call MatMatMult(cq_Mat, preconditioner_matrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, temp_Mat, ierr)
             print *,"Here comes c*q*matrix:"
             !call MatView(temp_Mat, PETSC_VIEWER_STDOUT_WORLD,ierr)
             call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "c_times_q_times_wholeMatrix.dat", FILE_MODE_WRITE, viewer, ierr)
             call MatView(temp_Mat, viewer, ierr)
             call PetscViewerDestroy(viewer, ierr)
             call MatDestroy(temp_Mat, ierr)

             ! Verify that a*p*b=0:
             call MatMatMult(preconditioner_matrix, pb_Mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, temp_Mat, ierr)
             print *,"Here comes matrix*p*b:"
             !call MatView(temp_Mat, PETSC_VIEWER_STDOUT_WORLD,ierr)
             call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "wholeMatrix_times_p_times_b.dat", FILE_MODE_WRITE, viewer, ierr)
             call MatView(temp_Mat, viewer, ierr)
             call PetscViewerDestroy(viewer, ierr)
             call MatDestroy(temp_Mat, ierr)

             ! The matrix given by the product of the constraints * p * sources will be used later. Compute it now.
             !call MatMatMult(constraints_Mat, sources_Mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, constraints_times_sources_Mat, ierr)
             call MatMatMult(constraints_Mat, pb_Mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, constraints_times_sources_Mat, ierr)
             if (verbose) then
                if (masterProc) print *,"Here comes sources_Mat:"
                call MatView(sources_Mat,PETSC_VIEWER_STDOUT_WORLD,ierr)
                if (masterProc) print *,"Here comes constraints_Mat:"
                call MatView(constraints_Mat,PETSC_VIEWER_STDOUT_WORLD,ierr)
                if (masterProc) print *,"Here comes constraints_times_sources_Mat:"
                call MatView(constraints_times_sources_Mat,PETSC_VIEWER_STDOUT_WORLD,ierr)
             end if
                if (masterProc) print *,"Here comes constraints_times_sources_Mat:"
                call MatView(constraints_times_sources_Mat,PETSC_VIEWER_STDOUT_WORLD,ierr)

             ! Set up the KSP for applying (c p b)^{-1}:
             call KSPCreate(PETSC_COMM_WORLD, constraints_times_sources_KSP, ierr)
             call KSPAppendOptionsPrefix(constraints_times_sources_KSP, 'constraints_times_sources_', ierr)
             call KSPSetType(constraints_times_sources_KSP, KSPPREONLY, ierr)
             call KSPGetPC(constraints_times_sources_KSP, constraints_times_sources_PC, ierr)
             !call PCSetType(constraints_times_sources_PC, PCLU, ierr)
             call PCSetType(constraints_times_sources_PC, PCREDUNDANT, ierr)
             call KSPSetOperators(constraints_times_sources_KSP, constraints_times_sources_Mat, constraints_times_sources_Mat, ierr)
             call KSPSetFromOptions(constraints_times_sources_KSP, ierr)
          end if

          if (verbose) then
             if (masterProc) print *,"Here is input_Vec:"
             call VecView(input_Vec, PETSC_VIEWER_STDOUT_WORLD, ierr)
          end if

          ! Compute the part of the solution in the DKE rows (x) arising from s:
          ! x = p b (c p b)^{-1} s
          call VecGetSubVector(input_Vec, IS_source_constraint, s_Vec, ierr) ! Extract s from the global input_Vec.
          if (verbose) then
             if (masterProc) print *,"Here is s_Vec:"
             call VecView(s_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
          end if
          call MatCreateVecs(constraints_times_sources_Mat, PETSC_NULL_OBJECT, temp_Vec_1, ierr) ! Create temp_Vec_1 of the appropriate size.
          
          if (masterProc) print *,"Here is temp_Vec_1:"
          call VecView(temp_Vec_1, PETSC_VIEWER_STDOUT_WORLD,ierr)

          call KSPSolve(constraints_times_sources_KSP, s_Vec, temp_Vec_1, ierr) ! Apply (c p b)^{-1}
          call KSPGetConvergedReason(constraints_times_sources_KSP, reason, ierr)
          if (reason <= 0 .and. masterProc) print *,"WARNING: constraints_times_sources_KSP failed (call 1) with reason",reason
          if (verbose) then
             if (masterProc) print *,"Here is (c p b)^{-1} s:"
             call VecView(temp_Vec_1, PETSC_VIEWER_STDOUT_WORLD,ierr)
          end if
          if (first_call) then
             if (masterProc) print *,"Here comes constraints_times_sources_KSP:"
             call KSPView(constraints_times_sources_KSP, PETSC_VIEWER_STDOUT_WORLD,ierr)
          end if
          call VecRestoreSubVector(input_Vec, IS_source_constraint,s_Vec, ierr)
          call MatMult(pb_Mat, temp_Vec_1, output_Vec, ierr) ! Left-multiply by p b, and store the result in output_Vec.
          call VecDestroy(temp_Vec_1, ierr)
          if (verbose) then
             if (masterProc) print *,"Here is output_Vec after the first term has been added:"
             call VecView(output_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
          end if

          ! Add the part of the solution in the DKE rows (x) arising from r:
          ! x += (I - p b (c p b)^{-1} c) a^{-1} (I - b (c q b)^{-1} c q) r
          call MatCreateVecs(constraints_times_sources_Mat, PETSC_NULL_OBJECT, cq_times_r_Vec, ierr) ! Create cq_times_r_Vec of the appropriate size.
          call MatMult(cq_Mat, input_Vec, cq_times_r_Vec, ierr) ! Form c*q*r.
          call MatCreateVecs(constraints_times_sources_Mat, PETSC_NULL_OBJECT, temp_Vec_1, ierr) ! Create temp_Vec_1 of the appropriate size.
          call KSPSolve(constraints_times_sources_KSP, cq_times_r_Vec, temp_Vec_1, ierr) ! Apply (c p b)^{-1}
          call KSPGetConvergedReason(constraints_times_sources_KSP, reason, ierr)
          if (reason <= 0 .and. masterProc) print *,"WARNING: constraints_times_sources_KSP failed (call 2) with reason",reason
          if (verbose) then
             if (masterProc) print *,"Here is (c q b)^{-1} c q r:"
             call VecView(temp_Vec_1, PETSC_VIEWER_STDOUT_WORLD, ierr)
          end if
          print *,myRank,'AAAA'
          call VecDuplicate(input_Vec, temp_Vec_2, ierr) ! Create temp_Vec_2 of the appropriate size.
          call MatMult(sources_Mat, temp_Vec_1, temp_Vec_2, ierr) ! Multiply by b to get b (c q b)^{-1} c q r. Result is now in temp_Vec_2.
          call VecDestroy(temp_Vec_1, ierr)
          call VecDuplicate(input_Vec, temp_Vec_1, ierr) ! Create temp_Vec_1 of the appropriate size.
          call VecWAXPY(temp_Vec_1, -one, temp_Vec_2, input_Vec, ierr) ! Now temp_Vec_1 holds (I - b (c q b)^{-1} c q) r.
          call VecDestroy(temp_Vec_2, ierr)
          call VecDuplicate(input_Vec, ainv_times_stuff_Vec, ierr) ! Create ainv_times_stuff_Vec of the appropriate size. 
          print *,myRank,'MMMMM'
          ! VVVVVVV  Next comes the big solve, involving multigrid. VVVVVVVVV
          call KSPSolve(inner_KSP, temp_Vec_1, ainv_times_stuff_Vec, ierr)
          ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          print *,myRank,'NNNN'
          call KSPGetConvergedReason(inner_KSP, reason, ierr)
          if (reason <= 0 .and. masterProc) print *,"WARNING: inner KSP failed with reason",reason
          call MatCreateVecs(constraints_times_sources_Mat, PETSC_NULL_OBJECT, temp_Vec_2, ierr) ! Create temp_Vec_2 of the appropriate size.
          call MatMult(constraints_Mat, ainv_times_stuff_Vec, temp_Vec_2, ierr) ! Form c*a^{-1}*....
          call MatCreateVecs(constraints_times_sources_Mat, PETSC_NULL_OBJECT, temp_Vec_1, ierr) ! Create temp_Vec_1 of the appropriate size.
          call KSPSolve(constraints_times_sources_KSP, temp_Vec_2, temp_Vec_1, ierr) ! Apply (c p b)^{-1}
          call KSPGetConvergedReason(constraints_times_sources_KSP, reason, ierr)
          if (reason <= 0 .and. masterProc) print *,"WARNING: constraints_times_sources_KSP failed (call 3) with reason",reason
          call VecDestroy(temp_Vec_2, ierr)
          call VecDuplicate(input_Vec, temp_Vec_2, ierr) ! Create temp_Vec_2 of the appropriate size.
          call MatMult(pb_Mat, temp_Vec_1, temp_Vec_2, ierr) ! Multiply by p b to get p b (c p b)^{-1} c a^{-1} ... Result is now in temp_Vec_2.
          call VecDestroy(temp_Vec_1, ierr)
          call VecDuplicate(input_Vec, temp_Vec_1, ierr) ! Create temp_Vec_1 of the appropriate size.
          call VecWAXPY(temp_Vec_1, -one, temp_Vec_2, ainv_times_stuff_Vec, ierr) ! Now temp_Vec_1 holds (I - p b (c p b)^{-1} c) a^{-1} ....
          call VecDestroy(temp_Vec_2, ierr)
          call VecAXPY(output_Vec, one, temp_Vec_1, ierr) ! Add the result to the term in output_Vec we computed previously.
          call VecDestroy(temp_Vec_1, ierr)

          if (verbose) then
             if (masterProc) print *,"Here is output_Vec after the DKE rows have been finished:"
             call VecView(output_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
          end if

          ! Compute the part of the solution in the constraint rows (y):
          ! y = (c q b)^{-1} c q r
          ! I think this over-writes what was previously in the constraint rows of output_Vec, but I'm not certain.
          call VecGetSubVector(output_Vec, IS_source_constraint, y_Vec, ierr) ! Extract y from the global output_Vec.
          call KSPSolve(constraints_times_sources_KSP, cq_times_r_Vec, y_Vec, ierr) ! Apply (c q b)^{-1}
          call KSPGetConvergedReason(constraints_times_sources_KSP, reason, ierr)
          if (reason <= 0 .and. masterProc) print *,"WARNING: constraints_times_sources_KSP failed (call 4) with reason",reason
          call VecRestoreSubVector(output_Vec, IS_source_constraint, y_Vec, ierr)
          call VecDestroy(cq_times_r_Vec, ierr)

          if (verbose) then
             if (masterProc) print *,"Here is the final output_Vec:"
             call VecView(output_Vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
          end if

       case default
          print *,"Error! Invalid fieldsplit_option:",fieldsplit_option
          stop
    end select

    if (first_call) then
       first_call = .false.
       if (masterProc) print *,"Here comes inner_KSP:"
       call KSPView(inner_KSP, PETSC_VIEWER_STDOUT_WORLD,ierr)
    end if

    call PetscLogEventEnd(event,ierr)

  end subroutine apply_preconditioner
