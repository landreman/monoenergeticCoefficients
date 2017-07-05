#include <petsc/finclude/petscmatdef.h>

subroutine populateMatrix(matrix)

  use petscmat

  use indices
  use variables
  use sparsify

  implicit none

  !#include <finclude/petscsys.h>
  !#include <finclude/petscvec.h>
  !#include <finclude/petscmat.h>
  !#include <finclude/petscdm.h>
  !#include <finclude/petscdmda.h>

  Mat :: matrix

  PetscErrorCode :: ierr
  MatNullSpace :: nullspace
  character(len=100) :: filename
  PetscViewer :: viewer

  PetscInt :: itheta, ithetaRow, ithetaCol
  PetscInt :: izeta, izetaRow, izetaCol
  PetscInt :: L, ell, N
  PetscInt :: rowIndex, colIndex, index
  PetscScalar, dimension(:,:), allocatable :: ddthetaToUse, ddzetaToUse
  PetscScalar :: temp, factor
  Vec :: real_collision_operator, temp_vec, null_vecs(2)
  Mat :: temp_mat, temp_mat2
  MatNullSpace :: null_space
  PetscBool :: null_space_test_result

  if (masterProc) then
     print *,"Entering populateMatrix."
  end if

  ! For simplicity, do all work on the master proc.

  nu_hat = nu / (G + iota * I)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now build a diagonal matrix that has the effect of applying dtheta dzeta sqrt_g \int dxi p_L (...) 
  ! to a vector of amplitudes of p_ell polynomials.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrixSize, weights_vec, ierr)
  if (masterProc) then
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           do L = 0,(Nxi-1)
              index = getIndex(itheta, izeta, L+1)
              call VecSetValue(weights_vec, index, L_scaling(L+1) * L_scaling(L+1) * (2 / (2 * L + one)) &
                   * (thetaWeight * zetaWeight / VPrime) * sqrt_g(itheta,izeta), INSERT_VALUES, ierr)
           end do
        end do
     end do
  end if
  call VecAssemblyBegin(weights_vec, ierr)
  call VecAssemblyEnd(weights_vec, ierr)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now build the 'real' collision operator, as a Vec
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Preallocate for only 1 element on the diagonal:
  !call MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
  !     1, PETSC_NULL_INTEGER, 0, PETSC_NULL_INTEGER, real_collision_operator, ierr)
  call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrixSize, real_collision_operator, ierr)
  if (masterProc) then
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           do L = 1,(Nxi-1) ! We can skip 0
              index = getIndex(itheta, izeta, L+1)
              call VecSetValue(real_collision_operator, index, -(nu_hat/2)*L*(L+1), INSERT_VALUES, ierr)
              !call MatSetValue(real_collision_operator, index, index, -(nu_hat/2)*L*(L+1), INSERT_VALUES, ierr)
           end do
        end do
     end do
  end if
  call VecAssemblyBegin(real_collision_operator, ierr)
  call VecAssemblyEnd(real_collision_operator, ierr)
  !call MatAssemblyBegin(real_collision_operator, MAT_FINAL_ASSEMBLY, ierr)
  !call MatAssemblyEnd(real_collision_operator, MAT_FINAL_ASSEMBLY, ierr)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now build \hat{C}^{-1}, the regularized inverse of the collision operator
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Preallocate for only 1 element on the diagonal:
  !call MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
  !     1, PETSC_NULL_INTEGER, 0, PETSC_NULL_INTEGER, CHat_inverse, ierr)
  call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrixSize, CHat_inverse, ierr)
  if (masterProc) then
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           do L = 1,(Nxi-1) ! We skip 0 to avoid divide-by-0
              index = getIndex(itheta, izeta, L+1)
              call VecSetValue(CHat_inverse, index, -2/(nu_hat*L*(L+1)), INSERT_VALUES, ierr)
              !call MatSetValue(CHat_inverse, index, index, -2/(nu_hat*L*(L+1)), INSERT_VALUES, ierr)
           end do
        end do
     end do
  end if
  call VecAssemblyBegin(CHat_inverse, ierr)
  call VecAssemblyEnd(CHat_inverse, ierr)
  !call MatAssemblyBegin(CHat_inverse, MAT_FINAL_ASSEMBLY, ierr)
  !call MatAssemblyEnd(CHat_inverse, MAT_FINAL_ASSEMBLY, ierr)

  print *,"MMMM"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now build the V operator
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  N = 3*(5+5-1) ! tridiagonal in xi. For each of these bands in xi, we are pentadiagonal in theta and zeta.
  call MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
       N, PETSC_NULL_INTEGER, N, PETSC_NULL_INTEGER, V_matrix, ierr)

  if (masterProc) then
     ! Add d/dtheta parallel streaming term
     itheta = -1
     izetaRow = -1
     izetaCol = -1
     do ithetaRow = 1,Ntheta
        do izeta = 1,Nzeta
           do L = 0,(Nxi-1)
              rowIndex = getIndex(ithetaRow, izeta, L+1)
              do ithetaCol = 1,Ntheta

                 ! Diagonal term (ExB)
                 ell = L
                 colIndex = getIndex(ithetaCol, izeta, ell+1)
                 call MatSetValueSparse(V_matrix, rowIndex, colIndex, &
                      iota*E/(FSAB2*sqrt_g(ithetaRow,izeta)) * ddtheta(ithetaRow,ithetaCol), ADD_VALUES, ierr)

                 ! Super-diagonal term
                 if (L < Nxi-1) then
                    ell = L + 1
                    colIndex = getIndex(ithetaCol, izeta, ell+1)
                    call MatSetValueSparse(V_matrix, rowIndex, colIndex, &
                      iota/(B(ithetaRow,izeta)*sqrt_g(ithetaRow,izeta)) * ddtheta(ithetaRow,ithetaCol) * L_scaling(ell+1)/L_scaling(L+1)*ell/(2*ell+1), &
                      ADD_VALUES, ierr)
                 end if

                 ! Sub-diagonal term
                 if (L > 0) then
                    ell = L - 1
                    colIndex = getIndex(ithetaCol, izeta, ell+1)
                    call MatSetValueSparse(V_matrix, rowIndex, colIndex, &
                      iota/(B(ithetaRow,izeta)*sqrt_g(ithetaRow,izeta)) * ddtheta(ithetaRow,ithetaCol) * L_scaling(ell+1)/L_scaling(L+1)*(ell+1)/(2*ell+1), &
                      ADD_VALUES, ierr)
                 end if

              end do
           end do
        end do
     end do

     ! Add d/dzeta parallel streaming term
     izeta = -1
     ithetaRow = -1
     ithetaCol = -1
     do izetaRow = 1,Nzeta
        do itheta = 1,Ntheta
           do L = 0,(Nxi-1)
              rowIndex = getIndex(itheta, izetaRow, L+1)
              do izetaCol = 1,Nzeta

                 ! Diagonal term (ExB):
                 ell = L
                 colIndex = getIndex(itheta, izetaCol, ell+1)
                 call MatSetValueSparse(V_matrix, rowIndex, colIndex, (-I/G)*iota*E/(FSAB2*sqrt_g(itheta,izetaRow)) * ddzeta(izetaRow,izetaCol), ADD_VALUES, ierr)

                 ! Super-diagonal term
                 if (L < Nxi-1) then
                    ell = L + 1
                    colIndex = getIndex(itheta, izetaCol, ell+1)
                    call MatSetValueSparse(V_matrix, rowIndex, colIndex, &
                         1/(B(itheta,izetaRow)*sqrt_g(itheta,izetaRow)) * ddzeta(izetaRow,izetaCol) * L_scaling(ell+1)/L_scaling(L+1)*ell/(2*ell+one), &
                         ADD_VALUES, ierr)
                 end if

                 ! Sub-diagonal term
                 if (L > 0) then
                    ell = L - 1
                    colIndex = getIndex(itheta, izetaCol, ell+1)
                    call MatSetValueSparse(V_matrix, rowIndex, colIndex, &
                         1/(B(itheta,izetaRow)*sqrt_g(itheta,izetaRow)) * ddzeta(izetaRow,izetaCol) * L_scaling(ell+1)/L_scaling(L+1)*(ell+1)/(2*ell+one), &
                         ADD_VALUES, ierr)
                 end if

              end do
           end do
        end do
     end do

     ! Add mirror term
     izetaRow = -1
     izetaCol = -1
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           temp = -(dBdtheta(itheta,izeta)*iota + dBdzeta(itheta,izeta)) / (2*B(itheta,izeta)*B(itheta,izeta)*sqrt_g(itheta,izeta))
           do L = 0,(Nxi-1)
              rowIndex = getIndex(itheta,izeta,L+1)

              ! Super-diagonal term
              if (L < Nxi-1) then
                 ell = L + 1
                 colIndex = getIndex(itheta,izeta,ell+1)
                 call MatSetValueSparse(V_matrix, rowIndex, colIndex, &
                      L_scaling(ell+1)/L_scaling(L+1)*(ell+1)*ell/(2*ell+1)*temp, ADD_VALUES, ierr)
              end if

              ! Super-diagonal term
              if (L > 0) then
                 ell = L - 1
                 colIndex = getIndex(itheta,izeta,ell+1)
                 call MatSetValueSparse(V_matrix, rowIndex, colIndex, &
                      -L_scaling(ell+1)/L_scaling(L+1)*(ell+1)*ell/(2*ell+1)*temp, ADD_VALUES, ierr)
              end if

           end do
        end do
     end do

  end if
  call MatAssemblyBegin(V_matrix, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(V_matrix, MAT_FINAL_ASSEMBLY, ierr)

  print *,"PPPP"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Add sources and constraints
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call preallocateMatrix(matrix)
  if (masterProc) then
     if (constraint_option==1) then
        print *,"Adding constraint rows and source columns"
        L = 0
        do itheta = 1,Ntheta
           do izeta = 1,Nzeta
              index = getIndex(itheta,izeta,L+1)
              factor = (thetaWeight * zetaWeight / VPrime) * sqrt_g(itheta,izeta)
              call MatSetValue(matrix, matrixSize-2, index, factor, ADD_VALUES, ierr)
              call MatSetValue(matrix, index, matrixSize-2, factor, ADD_VALUES, ierr)
              call MatSetValue(matrix, matrixSize-1, index + Ntheta*Nzeta*Nxi, factor, ADD_VALUES, ierr)
              call MatSetValue(matrix, index + Ntheta*Nzeta*Nxi, matrixSize-1, factor, ADD_VALUES, ierr)
           end do
        end do
     else
        print *,"NOT adding constraint rows and source columns"
     end if

     ! Add 0's on the diagonal, to avoid errors when we try to shift the diagonal later
     do index = 0,(matrixSize-1) ! Remember PETSc indices are 0-based
        call MatSetValue(matrix, index, index, zero, ADD_VALUES, ierr)
     end do
  end if

  

  call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)

  print *,"QQQ"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Construct the main block of the matrix:
  ! matrix = weights_matrix * (-real_collision_operator) + (-V') * weights_matrix * CHat_inv * V
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !call MatDuplicate(real_collision_operator, MAT_COPY_VALUES, temp_mat, ierr)
  !call MatScale(temp_mat, -one, ierr)
  !call MatDiagonalScale(temp_mat, weights_vec, PETSC_NULL_VEC, ierr)
  ! Now temp_mat holds -weights_vec * real_collision_operator

  call VecDuplicate(weights_vec, temp_vec, ierr)
  call VecPointwiseMult(temp_vec, weights_vec, real_collision_operator, ierr)
  call VecScale(temp_vec, -one, ierr)
  ! Now temp_vec holds -weights_vec * real_collision_operator 
  call MatDiagonalSet(matrix, temp_vec, INSERT_VALUES, ierr)
  print *,"SSSS"

  call VecPointwiseMult(temp_vec, weights_vec, CHat_inverse, ierr)
  call VecScale(temp_vec, -one, ierr)
  call MatDuplicate(V_matrix, MAT_COPY_VALUES, temp_mat, ierr)
  call MatDiagonalScale(temp_mat, temp_vec, PETSC_NULL_OBJECT, ierr)
  ! Now temp_mat holds -weights_matrix * CHat_inv * V. Next, left-multiply by V^T:
  call MatTransposeMatMult(V_matrix, temp_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, temp_mat2, ierr)
  ! Add result to the main matrix:
  call MatAXPY(matrix, one, temp_mat2, DIFFERENT_NONZERO_PATTERN, ierr)

  print *,"UUUU"

  call MatDestroy(temp_mat, ierr)
  call MatDestroy(temp_mat2, ierr)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now handle the part of the matrix involving \mathcal{F}(1)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Create a non-square (tall, skinny) matrix containing 1s on the diagonal, to pick out the part of V that acts on p_0:
  call MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, matrixSize, Ntheta*Nzeta, &
       1, PETSC_NULL_INTEGER, 1, PETSC_NULL_INTEGER, pick_out_p0_matrix, ierr)
  if (masterProc) then
     ! We assume here that xi is the most significant coordinate!!!
     index = -1
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           index = index + 1
           call MatSetValue(pick_out_p0_matrix, index, index, one, INSERT_VALUES, ierr)
        end do
     end do
  end if
  call MatAssemblyBegin(pick_out_p0_matrix, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(pick_out_p0_matrix, MAT_FINAL_ASSEMBLY, ierr) 

  print *,"VVVV"

  call MatMatMult(V_matrix, pick_out_p0_matrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, temp_mat, ierr)
  print *,"WWWW"
  call VecCopy(weights_vec, temp_vec, ierr)
  call VecScale(temp_vec, -one, ierr)
  call MatDiagonalScale(temp_mat, temp_vec, PETSC_NULL_OBJECT, ierr)
  ! Now temp_mat corresponds to -weights_vec * V(:,p0)
  print *,"XXX"
  ! Next, create a non-square (short, wide) matrix with 1s on the diagonal of the off-diagonal block
  call MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, Ntheta*Nzeta, matrixSize, &
       1, PETSC_NULL_INTEGER, 1, PETSC_NULL_INTEGER, injection_matrix, ierr)
  if (masterProc) then
     ! We assume here that xi is the most significant coordinate!!!
     rowIndex = -1
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           rowIndex = rowIndex + 1
           colIndex = rowIndex + Ntheta*Nzeta*Nxi
           call MatSetValue(injection_matrix, rowIndex, colIndex, one, INSERT_VALUES, ierr)
        end do
     end do
  end if
  call MatAssemblyBegin(injection_matrix, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(injection_matrix, MAT_FINAL_ASSEMBLY, ierr) 
  print *,"YYY"

  ! Now expand temp_mat into a square (matrixSize x matrixSize) matrix
  call MatMatMult(temp_mat, injection_matrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, temp_mat2, ierr)
  call MatDestroy(temp_mat, ierr)
  call MatAXPY(matrix, one, temp_mat2, DIFFERENT_NONZERO_PATTERN, ierr)
  ! Also add the transpose
  call MatTranspose(temp_mat2, MAT_INITIAL_MATRIX, temp_mat, ierr)
  call MatAXPY(matrix, one, temp_mat, DIFFERENT_NONZERO_PATTERN, ierr)

  if (constraint_option==2) then
     if (masterProc) print *,"Attaching a null space to the matrix."
     call VecDuplicate(weights_vec, null_vecs(1), ierr)
     call VecDuplicate(weights_vec, null_vecs(2), ierr)
     call VecSet(null_vecs(1), zero, ierr)
     call VecSet(null_vecs(2), zero, ierr)
     if (masterProc) then
        L = 0
        do itheta = 1,Ntheta
           do izeta = 1,Nzeta
              index = getIndex(itheta,izeta,L+1)
              call VecSetValue(null_vecs(1), index, one, INSERT_VALUES, ierr)
              call VecSetValue(null_vecs(2), index + Ntheta*Nzeta*Nxi, one, INSERT_VALUES, ierr)
           end do
        end do
     end if
     call VecAssemblyBegin(null_vecs(1), ierr)
     call VecAssemblyBegin(null_vecs(2), ierr)
     call VecAssemblyEnd(null_vecs(1), ierr)
     call VecAssemblyEnd(null_vecs(2), ierr)

     !print *,"Here comes null_vecs(1):"
     !call VecView(null_vecs(1), PETSC_VIEWER_STDOUT_WORLD,ierr)
     call PetscViewerASCIIOpen(PETSC_COMM_WORLD, "null_vecs_1.dat", viewer, ierr)
     call PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_INDEX, ierr)
     call VecView(null_vecs(1), viewer, ierr)
     call PetscViewerDestroy(viewer, ierr)

     !print *,"Here comes null_vecs(2):"
     !call VecView(null_vecs(2), PETSC_VIEWER_STDOUT_WORLD,ierr)
     call PetscViewerASCIIOpen(PETSC_COMM_WORLD, "null_vecs_2.dat", viewer, ierr)
     call PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_INDEX, ierr)
     call VecView(null_vecs(2), viewer, ierr)
     call PetscViewerDestroy(viewer, ierr)

     call MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, 2, null_vecs, null_space, ierr)
     call MatNullSpaceTest(null_space, matrix, null_space_test_result, ierr)
     if (masterProc) print *,"Result of null space test:",null_space_test_result,null_space_test_result.eqv.PETSC_TRUE
     call MatSetNullSpace(matrix, null_space, ierr)
     call MatSetTransposeNullSpace(matrix, null_space, ierr)
  else
     if (masterProc) print *,"NOT attaching a null space to the matrix."
  end if

!!$  ! The matrix has a 1D null space, with the null vector corresponding to a constant:
!!$  call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
!!$  call MatSetNullSpace(matrix,nullspace,ierr)
!!$  call MatNullSpaceDestroy(nullspace,ierr)

  !write (filename,fmt="(a,i1,a)") "mmc_matrix_",whichMatrix,".dat"
  filename = "mmc_matrix.dat"
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), FILE_MODE_WRITE, viewer, ierr)
  call MatView(matrix, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)

  filename = "mmc_V.dat"
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), FILE_MODE_WRITE, viewer, ierr)
  call MatView(V_matrix, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)


end subroutine populateMatrix


