#include <petsc/finclude/petscvecdef.h>

subroutine populateRHS(vec)

  use petscvec

  use indices
  use variables

  implicit none


  Vec :: vec, rhs_of_original_dke, R_with_weights, temp_vec, temp_vec2, little_vec
  PetscErrorCode :: ierr

  PetscInt :: itheta, izeta, L, index
  PetscReal :: valueToInsert
  PetscViewer :: viewer

  !call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrixSize, rhs_of_original_dke, ierr)
  call VecDuplicate(vec, rhs_of_original_dke, ierr)

  call VecSet(rhs_of_original_dke, 0.0d+0, ierr)

  if (masterProc) then
     ! For simplicity, do everything on proc 1.

     print *,"Entering populateRHS"
  
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           valueToInsert = (1.0d+0)/(B(itheta,izeta)*B(itheta,izeta)*B(itheta,izeta)*sqrt_g(itheta,izeta)) &
                * (G*dBdtheta(itheta,izeta)-I*dBdzeta(itheta,izeta))

           L=0
           index = getIndex(itheta,izeta,L+1)
           call VecSetValue(rhs_of_original_dke, index, valueToInsert*4/(3.0d+0) / L_scaling(L+1), INSERT_VALUES, ierr)

           L=2
           index = getIndex(itheta,izeta,L+1)
           call VecSetValue(rhs_of_original_dke, index, valueToInsert*2/(3.0d+0) / L_scaling(L+1), INSERT_VALUES, ierr)

        end do
     end do
  end if

  call VecAssemblyBegin(rhs_of_original_dke, ierr)
  call VecAssemblyEnd(rhs_of_original_dke, ierr)

  call VecDuplicate(rhs_of_original_dke, R_with_weights, ierr)
  call VecPointwiseMult(R_with_weights, rhs_of_original_dke, weights_vec, ierr)

  ! Next, assemble temp_vec2 = R_with_weights + (-V') * weights_matrix * CHat_inv * rhs_of_original_dke

  call VecDuplicate(rhs_of_original_dke, temp_vec, ierr)
  call VecPointwiseMult(temp_vec, Chat_inverse, R_with_weights, ierr)
  call VecScale(temp_vec, -one, ierr)

  call VecDuplicate(rhs_of_original_dke, temp_vec2, ierr)
  print *,"RHS MMM"
  call MatMultTransposeAdd(V_matrix, temp_vec, R_with_weights, temp_vec2, ierr)
  print *,"RHS NNN"

  ! Finally, assemble rhs = temp_vec2 + (R0 for the \mathcal{F}(1) equation)
  call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, Ntheta*Nzeta, little_vec, ierr)
  call MatMultTranspose(pick_out_p0_matrix, R_with_weights, little_vec, ierr) ! This line picks out R_0 from R
  print *,"RHS PPP"
  call MatMultTransposeAdd(injection_matrix, little_vec, temp_vec2, vec, ierr)
  print *,"RHS QQQ"

  !call VecView(vec, PETSC_VIEWER_STDOUT_WORLD,ierr)
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mmc_rhs.dat", FILE_MODE_WRITE, viewer, ierr)
  call VecView(vec, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)

  if (masterProc) print *,"Leaving populateRHS."

end subroutine populateRHS
