#include <petsc/finclude/petscmatdef.h>

subroutine populateMatrix(matrix)

  use petscmat
  use FourierConvolutionMatrixMod
  use FourierTransformMod
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
  !MatNullSpace :: nullspace
  !character(len=100) :: filename
  !PetscViewer :: viewer

  PetscInt :: imn_row, imn_col, imn
  PetscInt :: L, ell
  PetscInt :: rowIndex, colIndex, index
  PetscScalar, dimension(:), allocatable :: tempFourierVector
  PetscScalar, dimension(:,:), allocatable :: tempFourierMatrix
  PetscScalar :: temp

  if (masterProc) then
     print *,"Entering populateMatrix."
     print *,"matrixSize=",matrixSize
  end if

  if (masterProc) then
     ! For simplicity, do all work on the master proc.

     allocate(tempFourierVector(NFourier2))
     allocate(tempFourierMatrix(NFourier2,NFourier2))

     ! Add parallel streaming term
     call FourierTransform(B,tempFourierVector)
     call FourierConvolutionMatrix(tempFourierVector,tempFourierMatrix)
     tempFourierMatrix = matmul(tempFourierMatrix, iota*ddtheta+ddzeta)
     do L = 0,(Nxi-1)
        do imn_row = 1,NFourier2
           rowIndex = getIndex(imn_row, L+1)
           do imn_col = 1,NFourier2
              ! Super-diagonal term
              if (L < Nxi-1) then
                 ell = L + 1
                 colIndex = getIndex(imn_col, ell+1)
                 call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      tempFourierMatrix(imn_row,imn_col)*(L+1)/(2*L+3), ADD_VALUES, ierr)
              end if

              ! Sub-diagonal term
              if (L > 0) then
                 ell = L - 1
                 colIndex = getIndex(imn_col, ell+1)
                 call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      tempFourierMatrix(imn_row,imn_col)*L/(2*L-1), ADD_VALUES, ierr)
              end if

           end do
        end do
     end do

     ! Add mirror term
     call FourierTransform(iota*dBdtheta+dBdzeta,tempFourierVector)
     call FourierConvolutionMatrix(tempFourierVector,tempFourierMatrix)
     do imn_row = 1,NFourier2
        do imn_col = 1,NFourier2
           temp = -(0.5d+0)*tempFourierMatrix(imn_row,imn_col)
           do L = 0,(Nxi-1)
              rowIndex = getIndex(imn_row,L+1)

              ! Super-diagonal term
              if (L < Nxi-1) then
                 ell = L + 1
                 colIndex = getIndex(imn_col,ell+1)
                 call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      temp*(L+1)*(L+2)/(2*L+3), ADD_VALUES, ierr)
              end if

              ! Super-diagonal term
              if (L > 0) then
                 ell = L - 1
                 colIndex = getIndex(imn_col,ell+1)
                 call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      temp*(-L)*(L-1)/(2*L-1), ADD_VALUES, ierr)
              end if

           end do
        end do
     end do

     ! Add collision operator
     do imn = 1,NFourier2
        do L = 0,(Nxi-1)
           index = getIndex(imn,L+1)
           temp = nu/2*L*(L+1)           
           call MatSetValueSparse(matrix, index, index, temp, ADD_VALUES, ierr)
        end do
     end do

     ! Add source
     rowIndex=0 ! 0-based indices!
     colIndex = matrixSize-1 ! 0-based indices!
     temp = 1.0d+0
     call MatSetValue(matrix, rowIndex, colIndex, temp, ADD_VALUES, ierr)

     ! Add constraint
     call FourierTransform(1/(B*B),tempFourierVector)
     call FourierConvolutionMatrix(tempFourierVector,tempFourierMatrix)
     rowIndex = matrixSize-1 ! 0-based indices!
     L=0
     do imn_col = 1,NFourier2
        colIndex = getIndex(imn_col, L+1)
        call MatSetValueSparse(matrix,rowIndex,colIndex,tempFourierMatrix(1,imn_col), ADD_VALUES, ierr)
     end do

     ! Add some 0's on the diagonal, to keep PETSC happy:
     rowIndex = matrixSize-1  ! 0-based indices!
     colIndex = matrixSize-1  ! 0-based indices!
     temp = 0.0d+0
     call MatSetValue(matrix, rowIndex, colIndex, temp, ADD_VALUES, ierr)
     L=0
     do imn = 1,NFourier2
        index = getIndex(imn,L+1)
        call MatSetValue(matrix, index, index, temp, ADD_VALUES, ierr)
     end do

     deallocate(tempFourierVector,tempFourierMatrix)
  end if

  call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)

!!$  ! The matrix has a 1D null space, with the null vector corresponding to a constant:
!!$  call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
!!$  call MatSetNullSpace(matrix,nullspace,ierr)
!!$  call MatNullSpaceDestroy(nullspace,ierr)

!!$  !write (filename,fmt="(a,i1,a)") "mmc_matrix_",whichMatrix,".dat"
!!$  filename = "mmc_matrix.dat"
!!$  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), FILE_MODE_WRITE, viewer, ierr)
!!$  call MatView(matrix, viewer, ierr)
!!$  call PetscViewerDestroy(viewer, ierr)


end subroutine populateMatrix


