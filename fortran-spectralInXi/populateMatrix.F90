#include <finclude/petscmatdef.h>

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
  PetscInt :: L, ell
  PetscInt :: rowIndex, colIndex, index
  PetscScalar, dimension(:,:), allocatable :: ddthetaToUse, ddzetaToUse
  PetscScalar :: temp

  if (masterProc) then
     print *,"Entering populateMatrix."
  end if

  if (masterProc) then
     ! For simplicity, do all work on the master proc.

     ! Add d/dtheta parallel streaming term
     itheta = -1
     izetaRow = -1
     izetaCol = -1
     do ithetaRow = 1,Ntheta
        do izeta = 1,Nzeta
           do L = 0,(Nxi-1)
              rowIndex = getIndex(ithetaRow, izeta, L+1)
              do ithetaCol = 1,Ntheta

                 ! Super-diagonal term
                 if (L < Nxi-1) then
                    ell = L + 1
                    colIndex = getIndex(ithetaCol, izeta, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      iota*B(ithetaRow,izeta)*ddtheta(ithetaRow,ithetaCol)*(L+1)/(2*L+3), ADD_VALUES, ierr)
                 end if

                 ! Sub-diagonal term
                 if (L > 0) then
                    ell = L - 1
                    colIndex = getIndex(ithetaCol, izeta, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      iota*B(ithetaRow,izeta)*ddtheta(ithetaRow,ithetaCol)*L/(2*L-1), ADD_VALUES, ierr)
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

                 ! Super-diagonal term
                 if (L < Nxi-1) then
                    ell = L + 1
                    colIndex = getIndex(itheta, izetaCol, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         B(itheta,izetaRow)*ddzeta(izetaRow,izetaCol)*(L+1)/(2*L+3), ADD_VALUES, ierr)
                 end if

                 ! Sub-diagonal term
                 if (L > 0) then
                    ell = L - 1
                    colIndex = getIndex(itheta, izetaCol, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         B(itheta,izetaRow)*ddzeta(izetaRow,izetaCol)*L/(2*L-1), ADD_VALUES, ierr)
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
           temp = -(0.5d+0)*(dBdtheta(itheta,izeta)*iota + dBdzeta(itheta,izeta))
           do L = 0,(Nxi-1)
              rowIndex = getIndex(itheta,izeta,L+1)

              ! Super-diagonal term
              if (L < Nxi-1) then
                 ell = L + 1
                 colIndex = getIndex(itheta,izeta,ell+1)
                 call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      temp*(L+1)*(L+2)/(2*L+3), ADD_VALUES, ierr)
              end if

              ! Super-diagonal term
              if (L > 0) then
                 ell = L - 1
                 colIndex = getIndex(itheta,izeta,ell+1)
                 call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      temp*(-L)*(L-1)/(2*L-1), ADD_VALUES, ierr)
              end if

           end do
        end do
     end do

     ! Add collision operator
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           do L = 0,(Nxi-1)
              index = getIndex(itheta,izeta,L+1)
              temp = nu/2*L*(L+1)
              if (L==0) then
                 temp = diagonalShift
              end if

              call MatSetValueSparse(matrix, index, index, temp, ADD_VALUES, ierr)
           end do
        end do
     end do

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


