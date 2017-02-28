#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

subroutine populateMatrix(matrix, whichMatrix, level)

  use petscmat

  use indices
  use variables, Ntheta_fine => Ntheta, Nzeta_fine => Nzeta, Nxi_fine => Nxi, matrixSize_fine => matrixSize
  use sparsify

  implicit none

  !#include <finclude/petscsys.h>
  !#include <finclude/petscvec.h>
  !#include <finclude/petscmat.h>
  !#include <finclude/petscdm.h>
  !#include <finclude/petscdmda.h>

  Mat :: matrix
  PetscInt, intent(in) :: whichMatrix, level

  PetscErrorCode :: ierr
  MatNullSpace :: nullspace
  character(len=100) :: filename
  PetscViewer :: viewer

  PetscInt :: itheta, ithetaRow, ithetaCol
  PetscInt :: izeta, izetaRow, izetaCol
  PetscInt :: ixi, ixiRow, ixiCol
  PetscInt :: rowIndex, colIndex, j, irow
  PetscScalar, dimension(:,:), pointer :: derivative_matrix_to_use, pitch_angle_scattering_operator_to_use
  PetscScalar :: temp, factor
  PetscInt, allocatable, dimension(:) :: row_indices, col_indices
  PetscScalar, allocatable, dimension(:,:) :: values

  ! For convenience, use some short variable names to refer to quantities on this level:
  integer :: Ntheta, Nzeta, Nxi, matrixSize
  integer :: ithetaMin, ithetaMax, izetaMin, izetaMax
  PetscScalar, dimension(:), pointer :: xi, thetaWeights, zetaWeights, xiWeights
  PetscScalar, dimension(:,:), pointer :: B, dBdtheta, dBdzeta
  
  ! Values for whichMatrix:
  ! 0 = low-order preconditioner matrix, used only on the coarsest grid.
  ! 1 = high-order 'real' matrix, used only on the finest grid.
  ! 4 = Like the low-order preconditioner matrix, except the sources and constraints are not included, and there is a 1 on the corresponding diagonal
  ! 5 = Like 4, but only the diagonal and the off-diagonal-in-zeta terms.
  ! 6 = Like 4, but only the off-diagonal-in theta, xi, and x terms.
  ! 10 = Like 1, except low-order derivatives are used for the d/dzeta terms.

  if (masterProc) then
     print "(a,i4,a,i4)"," Entering populateMatrix for level ",level,", whichMatrix = ",whichMatrix
  end if

  ! For convenience, use some short variable names to refer to quantities on this level:
  Ntheta = levels(level)%Ntheta
  Nzeta  = levels(level)%Nzeta
  Nxi    = levels(level)%Nxi
  matrixSize = levels(level)%matrixSize
  ithetaMin = levels(level)%ithetaMin
  ithetaMax = levels(level)%ithetaMax
  izetaMin  = levels(level)%izetaMin
  izetaMax  = levels(level)%izetaMax
  xi => levels(level)%xi
  thetaWeights => levels(level)%thetaWeights
  zetaWeights  => levels(level)%zetaWeights
  xiWeights    => levels(level)%xiWeights
  B => levels(level)%B
  dBdtheta => levels(level)%dBdtheta
  dBdzeta => levels(level)%dBdzeta

  ! Ensure all diagonal entries are set
  if (masterProc .and. whichMatrix==6) then
  !if (masterProc) then
     do rowIndex = 0,matrixSize-1
        call MatSetValue(matrix,rowIndex,rowIndex,zero,ADD_VALUES,ierr)
     end do
  end if

  ! Add d/dtheta parallel streaming term
  itheta = -1
  izetaRow = -1
  izetaCol = -1
  do ixi = 1,Nxi
     do ithetaRow = ithetaMin,ithetaMax
        do izeta = izetaMin,izetaMax
           factor = iota*B(ithetaRow,izeta)*xi(ixi) + E*iota*B(ithetaRow,izeta)*B(ithetaRow,izeta)/FSAB2
           if (factor>0) then
              if (whichMatrix == 1 .or. whichMatrix == 10) then
                 derivative_matrix_to_use => levels(level)%ddtheta_plus
              else
                 derivative_matrix_to_use => levels(level)%ddtheta_plus_preconditioner
              end if
           else
              if (whichMatrix == 1 .or. whichMatrix == 10) then
                 derivative_matrix_to_use => levels(level)%ddtheta_minus
              else
                 derivative_matrix_to_use => levels(level)%ddtheta_minus_preconditioner
              end if
           end if
           rowIndex = getIndex(level,ithetaRow, izeta, ixi)
           do ithetaCol = 1,Ntheta
              if (whichMatrix==5 .and. ithetaCol .ne. ithetaRow) cycle ! For whichMatrix==5, add only the diagonal.
              if (whichMatrix==6 .and. ithetaCol == ithetaRow) cycle   ! For whichMatrix==6, add only the off-diagonal elements.
              colIndex = getIndex(level,ithetaCol, izeta, ixi)
              call MatSetValueSparse(matrix, rowIndex, colIndex, &
                   factor*derivative_matrix_to_use(ithetaRow,ithetaCol), ADD_VALUES, ierr)
           end do
        end do
     end do
  end do
  
  ! Add d/dzeta parallel streaming term
  if (whichMatrix .ne. 6) then
     izeta = -1
     ithetaRow = -1
     ithetaCol = -1
     do ixi = 1,Nxi
        do izetaRow = izetaMin,izetaMax
           do itheta = ithetaMin,ithetaMax
              factor = B(itheta,izetaRow)*xi(ixi) - E*iota*I*B(itheta,izetaRow)*B(itheta,izetaRow)/(G*FSAB2)
              if (factor>0) then
                 if (whichMatrix == 1) then
                    derivative_matrix_to_use => levels(level)%ddzeta_plus
                 else
                    derivative_matrix_to_use => levels(level)%ddzeta_plus_preconditioner
                 end if
              else
                 if (whichMatrix == 1) then
                    derivative_matrix_to_use => levels(level)%ddzeta_minus
                 else
                    derivative_matrix_to_use => levels(level)%ddzeta_minus_preconditioner
                 end if
              end if
              rowIndex = getIndex(level,itheta, izetaRow, ixi)
              do izetaCol = 1,Nzeta
                 colIndex = getIndex(level,itheta, izetaCol, ixi)
                 call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      factor*derivative_matrix_to_use(izetaRow,izetaCol), ADD_VALUES, ierr)
              end do
           end do
        end do
     end do
  end if
  
  ! Add mirror term and collision operator
  izetaRow = -1
  izetaCol = -1
  ixi = -1
  if (whichMatrix==1 .or. whichMatrix == 10) then
     pitch_angle_scattering_operator_to_use => levels(level)%pitch_angle_scattering_operator
  else
     pitch_angle_scattering_operator_to_use => levels(level)%pitch_angle_scattering_operator_preconditioner
  end if
  do itheta = ithetaMin,ithetaMax
     do izeta = izetaMin,izetaMax
        do ixiRow = 1,Nxi
           rowIndex = getIndex(level,itheta,izeta,ixiRow)
           temp = -(0.5d+0)*(1-xi(ixiRow)*xi(ixiRow))*(dBdtheta(itheta,izeta)*iota + dBdzeta(itheta,izeta))
           if (temp>0) then
              if (whichMatrix == 1 .or. whichMatrix == 10) then
                 derivative_matrix_to_use => levels(level)%ddxi_plus
              else
                 derivative_matrix_to_use => levels(level)%ddxi_plus_preconditioner
              end if
           else
              if (whichMatrix == 1 .or. whichMatrix == 10) then
                 derivative_matrix_to_use => levels(level)%ddxi_minus
              else
                 derivative_matrix_to_use => levels(level)%ddxi_minus_preconditioner
              end if
           end if
           do ixiCol = 1,Nxi
              if (whichMatrix==5 .and. ixiCol .ne. ixiRow) cycle ! For whichMatrix==5, add only the diagonal.
              if (whichMatrix==6 .and. ixiCol  ==  ixiRow) cycle ! For whichMatrix==6, add only the off-diagonal elements.
              colIndex = getIndex(level,itheta,izeta,ixiCol)
              call MatSetValueSparse(matrix, rowIndex, colIndex, &
                   temp * derivative_matrix_to_use(ixiRow,ixiCol) &  ! Mirror term
                   - nu * pitch_angle_scattering_operator_to_use(ixiRow,ixiCol), &    ! Collision operator
                   ADD_VALUES, ierr)
           end do
        end do
     end do
  end do

  ! Add source term
  if (constraint_option==1 .and. masterProc .and. (whichMatrix==0 .or. whichMatrix==1 .or. whichMatrix == 10)) then
     do irow = 1,matrixSize-1
        call MatSetValue(matrix,irow-1,matrixSize-1,one,ADD_VALUES,ierr) ! -1 because PETSc uses 0-based indices.
     end do
  end if

  ! Add constraint row
  if (constraint_option==1 .and. masterProc .and. (whichMatrix==0 .or. whichMatrix==1 .or. whichMatrix == 10)) then
     allocate(row_indices(1))
     allocate(col_indices(matrixSize-1))
     allocate(values(matrixSize-1,1))
     row_indices = matrixSize - 1  ! -1 because PETSc uses 0-based indexing
     col_indices = [(j, j = 0,matrixSize-2 )]
     do izeta = 1,Nzeta
        do itheta = 1,Ntheta
           do ixi = 1,Nxi
              values(getIndex(level,itheta,izeta,ixi)+1,1) = (xiWeights(ixi)*thetaWeights(itheta)*zetaWeights(izeta))/(VPrime*B(itheta,izeta)*B(itheta,izeta))
           end do
        end do
     end do
     call MatSetValues(matrix, 1, row_indices, matrixSize-1, col_indices, values, ADD_VALUES, ierr)
     deallocate(row_indices, col_indices, values)
  end if

  ! For smoothing, we may need a 1 on the diagonal element corresponding to the sources/constraints:
  if (constraint_option==1 .and. masterProc .and. (whichMatrix==4 .or. whichMatrix==5)) then
     call MatSetValue(matrix,matrixSize-1,matrixSize-1,one,ADD_VALUES,ierr) ! -1 because PETSc uses 0-based indices.
  end if

  ! Done with adding matrix elements.
  
  call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)

  if (constraint_option==2 .and. whichMatrix==1) then
     ! The matrix has a 1D null space, with the null vector corresponding to a constant:
     call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
     call MatSetNullSpace(matrix,nullspace,ierr)
     call MatNullSpaceDestroy(nullspace,ierr)
  end if

  write (filename,fmt="(a,i1,a,i1,a)") "mmc_matrix_level_",level,'_whichMatrix_',whichMatrix,".dat"
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), FILE_MODE_WRITE, viewer, ierr)
  call MatView(matrix, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)




end subroutine populateMatrix


