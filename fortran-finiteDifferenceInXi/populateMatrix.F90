#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

subroutine populateMatrix(matrix, whichMatrix)

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
  PetscInt, intent(in) :: whichMatrix

  PetscErrorCode :: ierr
  MatNullSpace :: nullspace
  character(len=100) :: filename
  PetscViewer :: viewer

  PetscInt :: itheta, ithetaRow, ithetaCol
  PetscInt :: izeta, izetaRow, izetaCol
  PetscInt :: ixi, ixiRow, ixiCol
  PetscInt :: rowIndex, colIndex, irow, j
  PetscScalar, dimension(:,:), pointer :: derivative_matrix_to_use, pitch_angle_scattering_operator_to_use
  PetscScalar :: temp, factor
  PetscInt, allocatable, dimension(:) :: row_indices, col_indices
  PetscScalar, allocatable, dimension(:,:) :: values
  Vec :: null_Vecs(1)

  if (masterProc) then
     print *,"Entering populateMatrix for whichMatrix = ",whichMatrix
  end if

  ! Add d/dtheta parallel streaming term
  itheta = -1
  izetaRow = -1
  izetaCol = -1
  do ixi = 1,Nxi
     do ithetaRow = ithetaMin,ithetaMax
        do izeta = izetaMin,izetaMax
           factor = iota*B(ithetaRow,izeta)*xi(ixi) + E*iota*B(ithetaRow,izeta)*B(ithetaRow,izeta)/FSAB2
           if (factor > 0) then
              if (whichMatrix == 1) then
                 derivative_matrix_to_use => ddtheta_plus
              else
                 derivative_matrix_to_use => ddtheta_plus_preconditioner
              end if
           else
              if (whichMatrix == 1) then
                 derivative_matrix_to_use => ddtheta_minus
              else
                 derivative_matrix_to_use => ddtheta_minus_preconditioner
              end if
           end if
           rowIndex = getIndex(ithetaRow, izeta, ixi)
           do ithetaCol = 1,Ntheta
              colIndex = getIndex(ithetaCol, izeta, ixi)
              call MatSetValueSparse(matrix, rowIndex, colIndex, &
                   factor*derivative_matrix_to_use(ithetaRow,ithetaCol), ADD_VALUES, ierr)
           end do
        end do
     end do
  end do
  
  ! Add d/dzeta parallel streaming term
  izeta = -1
  ithetaRow = -1
  ithetaCol = -1
  do ixi = 1,Nxi
     do izetaRow = izetaMin,izetaMax
        do itheta = ithetaMin,ithetaMax
           factor = B(itheta,izetaRow)*xi(ixi) - E*iota*I*B(itheta,izetaRow)*B(itheta,izetaRow)/(G*FSAB2)
           if (factor>0) then
              if (whichMatrix == 1) then
                 derivative_matrix_to_use => ddzeta_plus
              else
                 derivative_matrix_to_use => ddzeta_plus_preconditioner
              end if
           else
              if (whichMatrix == 1) then
                 derivative_matrix_to_use => ddzeta_minus
              else
                 derivative_matrix_to_use => ddzeta_minus_preconditioner
              end if
           end if
           rowIndex = getIndex(itheta, izetaRow, ixi)
           do izetaCol = 1,Nzeta
              colIndex = getIndex(itheta, izetaCol, ixi)
              call MatSetValueSparse(matrix, rowIndex, colIndex, &
                   factor*derivative_matrix_to_use(izetaRow,izetaCol), ADD_VALUES, ierr)
           end do
        end do
     end do
  end do
  
  ! Add mirror term and collision operator
  izetaRow = -1
  izetaCol = -1
  ixi = -1
  if (whichMatrix==1) then
     pitch_angle_scattering_operator_to_use => pitch_angle_scattering_operator
  else
     pitch_angle_scattering_operator_to_use => pitch_angle_scattering_operator_preconditioner
  end if
  do itheta = ithetaMin,ithetaMax
     do izeta = izetaMin,izetaMax
        do ixiRow = 1,Nxi
           rowIndex = getIndex(itheta,izeta,ixiRow)
           temp = -(0.5d+0)*(1-xi(ixiRow)*xi(ixiRow))*(dBdtheta(itheta,izeta)*iota + dBdzeta(itheta,izeta))
           if (temp>0) then
              if (whichMatrix == 1) then
                 derivative_matrix_to_use => ddxi_plus
              else
                 derivative_matrix_to_use => ddxi_plus_preconditioner
              end if
           else
              if (whichMatrix == 1) then
                 derivative_matrix_to_use => ddxi_minus
              else
                 derivative_matrix_to_use => ddxi_minus_preconditioner
              end if
           end if
           do ixiCol = 1,Nxi
              colIndex = getIndex(itheta,izeta,ixiCol)
              call MatSetValueSparse(matrix, rowIndex, colIndex, &
                   temp * derivative_matrix_to_use(ixiRow,ixiCol) &  ! Mirror term
                   - nu * pitch_angle_scattering_operator_to_use(ixiRow,ixiCol), &    ! Collision operator
                   ADD_VALUES, ierr)
           end do
        end do
     end do
  end do
  
  ! Add source term
  if (constraint_option==1 .and. masterProc) then
     do irow = 1,matrixSize-1
        call MatSetValue(matrix,irow-1,matrixSize-1,one,ADD_VALUES,ierr) ! -1 because PETSc uses 0-based indices.
     end do
  end if


  ! Add constraint row
  if (constraint_option==1 .and. masterProc) then
     allocate(row_indices(1))
     allocate(col_indices(matrixSize-1))
     allocate(values(matrixSize-1,1))
     row_indices = matrixSize - 1  ! -1 because PETSc uses 0-based indexing
     col_indices = [(j, j = 0,matrixSize-2 )]
     do izeta = 1,Nzeta
        do itheta = 1,Ntheta
           do ixi = 1,Nxi
              values(getIndex(itheta,izeta,ixi)+1,1) = (xiWeights(ixi)*thetaWeights(itheta)*zetaWeights(izeta))/(VPrime*B(itheta,izeta)*B(itheta,izeta))
           end do
        end do
     end do
     call MatSetValues(matrix, 1, row_indices, matrixSize-1, col_indices, values, ADD_VALUES, ierr)
     deallocate(row_indices, col_indices, values)
  end if

  ! Done with adding matrix elements.
  call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)

  if (constraint_option==2 .or. constraint_option==3) then
     ! The matrix has a 1D null space, with the null vector corresponding to a constant:
     call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
     call MatSetNullSpace(matrix,nullspace,ierr)
     call MatNullSpaceDestroy(nullspace,ierr)
     if (masterProc) print *,"Adding null space."
  end if

  if (constraint_option==3) then
     call MatCreateVecs(matrix, PETSC_NULL_OBJECT, null_Vecs(1), ierr)
     if (masterProc) then
        print *,"Adding transpose null space."
        do izeta = 1,Nzeta
           do itheta = 1,Ntheta
              do ixi = 1,Nxi
                 call VecSetValue(null_Vecs(1), getIndex(itheta,izeta,ixi), (xiWeights(ixi)*thetaWeights(itheta)*zetaWeights(izeta))/(VPrime*B(itheta,izeta)*B(itheta,izeta)), INSERT_VALUES, ierr)
              end do
           end do
        end do
     end if
     call VecAssemblyBegin(null_Vecs(1), ierr)
     call VecAssemblyEnd(null_Vecs(1), ierr)
     call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_FALSE,1,null_Vecs,nullspace,ierr)
     call MatSetTransposeNullSpace(matrix,nullspace,ierr)
     call MatNullSpaceDestroy(nullspace,ierr)
     call VecDestroy(null_Vecs(1),ierr)

     !call MatView(matrix,PETSC_VIEWER_STDOUT_WORLD,ierr)
  end if

  write (filename,fmt="(a,i1,a)") "mmc_matrix_",whichMatrix,".dat"
  !call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), FILE_MODE_WRITE, viewer, ierr)
  !call MatView(matrix, viewer, ierr)
  !call PetscViewerDestroy(viewer, ierr)




end subroutine populateMatrix


