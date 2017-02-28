#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

subroutine preallocateMatrix(matrix, whichMatrix, level)

  use petscmat

  use variables, only: levels, numProcs, masterProc, constraint_option

  implicit none

  integer, intent(in) :: whichMatrix, level
  Mat :: matrix
  integer :: predictedNNZForEachRowOfPreconditioner, predictedNNZForEachRowOfTotalMatrix
  integer, dimension(:), allocatable :: predictedNNZsForEachRow, predictedNNZsForEachRowDiagonal
  PetscErrorCode :: ierr
  integer :: tempInt1, i, itheta, izeta, ispecies, index
  integer :: firstRowThisProcOwns, lastRowThisProcOwns, numLocalRows
  integer :: Ntheta, Nzeta, Nxi, matrixSize

  if (masterProc) then
     print "(a,i4,a,i4)"," Beginning preallocation for level ",level,", whichMatrix = ",whichMatrix
  end if

  Ntheta = levels(level)%Ntheta
  Nzeta  = levels(level)%Nzeta
  Nxi    = levels(level)%Nxi
  matrixSize = levels(level)%matrixSize

  allocate(predictedNNZsForEachRow(matrixSize))
  allocate(predictedNNZsForEachRowDiagonal(matrixSize))
  ! Set tempInt1 to the expected number of nonzeros in a row of the kinetic equation block:
  tempInt1 = (7+7)*5 + 1 ! The ddtheta and ddzeta terms are penta-diagonal in L, with up to 7 nonzeros per row in each L. Add 1 for source.
  if (tempInt1 > matrixSize) then
     tempInt1 = matrixSize
  end if
  predictedNNZsForEachRow = tempInt1

  predictedNNZsForEachRowDiagonal = predictedNNZsForEachRow


  select case (constraint_option)
  case (0,2)
  case (1)
     if (whichMatrix==4) then
        predictedNNZsForEachRow(matrixSize) = 1
     else
        predictedNNZsForEachRow(matrixSize) = Ntheta*Nzeta + 1  ! +1 for diagonal
     end if
  case default
     print *,"Invalid constraint_option:",constraint_option
     stop
  end select

  predictedNNZsForEachRowDiagonal = predictedNNZsForEachRow

  ! 1 = New method with lower, more precise estimated number-of-nonzeros.
  ! This method is more complicated, but it should use much less memory.

  call MatCreate(PETSC_COMM_WORLD, matrix, ierr)
  !call MatSetType(matrix, MATMPIAIJ, ierr)
  call MatSetType(matrix, MATAIJ, ierr)

  numLocalRows = PETSC_DECIDE
  call PetscSplitOwnership(PETSC_COMM_WORLD, numLocalRows, matrixSize, ierr)

  call MatSetSizes(matrix, numLocalRows, numLocalRows, PETSC_DETERMINE, PETSC_DETERMINE, ierr)

  ! We first pre-allocate assuming number-of-nonzeros = 0, because due to a quirk in PETSc,
  ! MatGetOwnershipRange only works after MatXXXSetPreallocation is called:
  if (numProcs == 1) then
     call MatSeqAIJSetPreallocation(matrix, 0, PETSC_NULL_INTEGER, ierr)
  else
     call MatMPIAIJSetPreallocation(matrix, 0, PETSC_NULL_INTEGER, 0, PETSC_NULL_INTEGER, ierr)
  end if

  call MatGetOwnershipRange(matrix, firstRowThisProcOwns, lastRowThisProcOwns, ierr)
  !print *,"I am proc ",myRank," and I own rows ",firstRowThisProcOwns," to ",lastRowThisProcOwns-1

  ! To avoid a PETSc error message, the predicted nnz for each row of the diagonal blocks must be no greater than the # of columns this proc owns:
  ! But we must not lower the predicted nnz for the off-diagonal blocks, because then the total predicted nnz for the row
  ! would be too low.
  tempInt1 = lastRowThisProcOwns - firstRowThisProcOwns
  do i=firstRowThisProcOwns+1,lastRowThisProcOwns
     if (predictedNNZsForEachRowDiagonal(i) > tempInt1) then
        predictedNNZsForEachRowDiagonal(i) = tempInt1
     end if
  end do

  ! Now, set the real estimated number-of-nonzeros:
  if (numProcs == 1) then
     call MatSeqAIJSetPreallocation(matrix, 0, predictedNNZsForEachRow(firstRowThisProcOwns+1:lastRowThisProcOwns), ierr)
  else
     call MatMPIAIJSetPreallocation(matrix, &
          0, predictedNNZsForEachRowDiagonal(firstRowThisProcOwns+1:lastRowThisProcOwns), &
          0, predictedNNZsForEachRow(firstRowThisProcOwns+1:lastRowThisProcOwns), ierr)
  end if


  ! If any mallocs are required during matrix assembly, do not generate an error:
  !call MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)


end subroutine preallocateMatrix
