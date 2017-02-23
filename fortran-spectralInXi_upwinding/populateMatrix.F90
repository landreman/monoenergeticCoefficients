#include <petsc/finclude/petscmatdef.h>

subroutine populateMatrix(matrix, whichMatrix)

  use petscmat

  use indices
  use variables
  use sparsify

  implicit none

  Mat :: matrix

  PetscErrorCode :: ierr
  MatNullSpace :: nullspace
  character(len=100) :: filename
  PetscViewer :: viewer
  integer :: whichMatrix

  PetscInt :: itheta, ithetaRow, ithetaCol
  PetscInt :: izeta, izetaRow, izetaCol
  PetscInt :: L, ell
  PetscInt :: rowIndex, colIndex, index
  PetscScalar, dimension(:,:), allocatable :: ddtheta_sum_to_use, ddtheta_difference_to_use
  PetscScalar, dimension(:,:), allocatable :: ddzeta_sum_to_use, ddzeta_difference_to_use
  PetscScalar :: temp, factor, L_factor
  Vec :: null_vec(1)

  if (masterProc) then
     print *,"Entering populateMatrix with whichMatrix =",whichMatrix
  end if

  if (masterProc) then
     ! For simplicity, do all work on the master proc.

     ! Add d/dtheta parallel streaming term
     allocate(ddtheta_sum_to_use(Ntheta,Ntheta))
     allocate(ddtheta_difference_to_use(Ntheta,Ntheta))
     do L=0,(Nxi-1)

        if (whichMatrix>0) then
           ddtheta_sum_to_use = ddtheta_sum
           ddtheta_difference_to_use = ddtheta_difference
        else
           ddtheta_sum_to_use = ddtheta_sum_preconditioner
           ddtheta_difference_to_use = ddtheta_difference_preconditioner
        end if
        if (iota < 0) ddtheta_difference_to_use = -ddtheta_difference_to_use
        
        do ithetaRow = 1,Ntheta
           do izeta = 1,Nzeta
              factor = iota * B(ithetaRow,izeta)
              rowIndex = getIndex(ithetaRow, izeta, L+1)
              do ithetaCol = 1,Ntheta
                    
                 ! Diagonal-in-L term:
                 ell = L
                 L_factor = (2*L*L+2*L-one)/((2*L+3)*(2*L-one))
                 colIndex = getIndex(ithetaCol, izeta, ell+1)
                 call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      L_factor*factor*ddtheta_difference_to_use(ithetaRow,ithetaCol), ADD_VALUES, ierr)
                 
                 ! Super-diagonal-in-L term
                 if (L < Nxi-1) then
                    ell = L + 1
                    L_factor = (L+1)/(two*L+3)
                    colIndex = getIndex(ithetaCol, izeta, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         L_factor*factor*ddtheta_sum_to_use(ithetaRow,ithetaCol), ADD_VALUES, ierr)
                 end if
                 
                 ! Sub-diagonal-in-L term
                 if (L > 0) then
                    ell = L - 1
                    L_factor = L/(2*L-one)
                    colIndex = getIndex(ithetaCol, izeta, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         L_factor*factor*ddtheta_sum_to_use(ithetaRow,ithetaCol), ADD_VALUES, ierr)
                 end if
                 
                 if (whichMatrix>0) then
                    ! Super-super-diagonal-in-L term
                    if (L < Nxi-2) then
                       ell = L + 2
                       L_factor = (L+two)*(L+one)/((two*L+5)*(two*L+3))
                       colIndex = getIndex(ithetaCol, izeta, ell+1)
                       call MatSetValueSparse(matrix, rowIndex, colIndex, &
                            L_factor*factor*ddtheta_difference_to_use(ithetaRow,ithetaCol), ADD_VALUES, ierr)
                    end if
                    
                    ! Sub-sub-diagonal-in-L term
                    if (L > 1) then
                       ell = L - 2
                       L_factor = (L-one)*L/((two*L-3)*(two*L-one))
                       colIndex = getIndex(ithetaCol, izeta, ell+1)
                       call MatSetValueSparse(matrix, rowIndex, colIndex, &
                            L_factor*factor*ddtheta_difference_to_use(ithetaRow,ithetaCol), ADD_VALUES, ierr)
                    end if
                 end if
              end do
           end do
        end do
     end do
     deallocate(ddtheta_sum_to_use, ddtheta_difference_to_use)


!!$     itheta = -1
!!$     izetaRow = -1
!!$     izetaCol = -1
!!$     do ithetaRow = 1,Ntheta
!!$        do izeta = 1,Nzeta
!!$           do L = 0,(Nxi-1)
!!$              rowIndex = getIndex(ithetaRow, izeta, L+1)
!!$              do ithetaCol = 1,Ntheta
!!$
!!$                 ! Super-diagonal term
!!$                 if (L < Nxi-1) then
!!$                    ell = L + 1
!!$                    colIndex = getIndex(ithetaCol, izeta, ell+1)
!!$                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
!!$                      iota*B(ithetaRow,izeta)*ddtheta(ithetaRow,ithetaCol)*(L+1)/(2*L+3), ADD_VALUES, ierr)
!!$                 end if
!!$
!!$                 ! Sub-diagonal term
!!$                 if (L > 0) then
!!$                    ell = L - 1
!!$                    colIndex = getIndex(ithetaCol, izeta, ell+1)
!!$                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
!!$                      iota*B(ithetaRow,izeta)*ddtheta(ithetaRow,ithetaCol)*L/(2*L-1), ADD_VALUES, ierr)
!!$                 end if
!!$
!!$              end do
!!$           end do
!!$        end do
!!$     end do

     ! Add d/dzeta parallel streaming term
     allocate(ddzeta_sum_to_use(Nzeta,Nzeta))
     allocate(ddzeta_difference_to_use(Nzeta,Nzeta))
     do L=0,(Nxi-1)

        if (whichMatrix>0) then
           ddzeta_sum_to_use = ddzeta_sum
           ddzeta_difference_to_use = ddzeta_difference
        else
           ddzeta_sum_to_use = ddzeta_sum_preconditioner
           ddzeta_difference_to_use = ddzeta_difference_preconditioner
        end if
        
        do izetaRow = 1,Nzeta
           do itheta = 1,Ntheta
              factor = B(itheta,izetaRow)
              rowIndex = getIndex(itheta, izetaRow, L+1)
              do izetaCol = 1,Nzeta
                    
                 ! Diagonal-in-L term:
                 ell = L
                 L_factor = (2*L*L+2*L-one)/((2*L+3)*(2*L-one))
                 colIndex = getIndex(itheta, izetaCol, ell+1)
                 call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      L_factor*factor*ddzeta_difference_to_use(izetaRow,izetaCol), ADD_VALUES, ierr)
                 
                 ! Super-diagonal-in-L term
                 if (L < Nxi-1) then
                    ell = L + 1
                    L_factor = (L+1)/(two*L+3)
                    colIndex = getIndex(itheta, izetaCol, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         L_factor*factor*ddzeta_sum_to_use(izetaRow,izetaCol), ADD_VALUES, ierr)
                 end if
                 
                 ! Sub-diagonal-in-L term
                 if (L > 0) then
                    ell = L - 1
                    L_factor = L/(2*L-one)
                    colIndex = getIndex(itheta, izetaCol, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         L_factor*factor*ddzeta_sum_to_use(izetaRow,izetaCol), ADD_VALUES, ierr)
                 end if
                 
                 if (whichMatrix>0) then
                    ! Super-super-diagonal-in-L term
                    if (L < Nxi-2) then
                       ell = L + 2
                       L_factor = (L+two)*(L+one)/((two*L+5)*(two*L+3))
                       colIndex = getIndex(itheta, izetaCol, ell+1)
                       call MatSetValueSparse(matrix, rowIndex, colIndex, &
                            L_factor*factor*ddzeta_difference_to_use(izetaRow,izetaCol), ADD_VALUES, ierr)
                    end if
                    
                    ! Sub-sub-diagonal-in-L term
                    if (L > 1) then
                       ell = L - 2
                       L_factor = (L-one)*L/((two*L-3)*(two*L-one))
                       colIndex = getIndex(itheta, izetaCol, ell+1)
                       call MatSetValueSparse(matrix, rowIndex, colIndex, &
                            L_factor*factor*ddzeta_difference_to_use(izetaRow,izetaCol), ADD_VALUES, ierr)
                    end if
                 end if
              end do
           end do
        end do
     end do
     deallocate(ddzeta_sum_to_use, ddzeta_difference_to_use)

!!$     izeta = -1
!!$     ithetaRow = -1
!!$     ithetaCol = -1
!!$     do izetaRow = 1,Nzeta
!!$        do itheta = 1,Ntheta
!!$           do L = 0,(Nxi-1)
!!$              rowIndex = getIndex(itheta, izetaRow, L+1)
!!$              do izetaCol = 1,Nzeta
!!$
!!$                 ! Super-diagonal term
!!$                 if (L < Nxi-1) then
!!$                    ell = L + 1
!!$                    colIndex = getIndex(itheta, izetaCol, ell+1)
!!$                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
!!$                         B(itheta,izetaRow)*ddzeta(izetaRow,izetaCol)*(L+1)/(2*L+3), ADD_VALUES, ierr)
!!$                 end if
!!$
!!$                 ! Sub-diagonal term
!!$                 if (L > 0) then
!!$                    ell = L - 1
!!$                    colIndex = getIndex(itheta, izetaCol, ell+1)
!!$                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
!!$                         B(itheta,izetaRow)*ddzeta(izetaRow,izetaCol)*L/(2*L-1), ADD_VALUES, ierr)
!!$                 end if
!!$
!!$              end do
!!$           end do
!!$        end do
!!$     end do

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
              ! Let's use MatSetValue instead of MatSetValueSparse here to ensure the diagonal is set even if it is 0.
              call MatSetValue(matrix, index, index, temp, ADD_VALUES, ierr)
           end do
        end do
     end do

     if (constraint_option==1) then
        ! Add source
        do itheta = 1,Ntheta
           do izeta = 1,Nzeta
              index = getIndex(itheta,izeta,1)
              call MatSetValue(matrix, index, matrixSize-1, one, ADD_VALUES, ierr)
           end do
        end do

        ! Add constraint
        do itheta = 1,Ntheta
           do izeta = 1,Nzeta
              index = getIndex(itheta,izeta,1)
              call MatSetValue(matrix, matrixSize-1, index, one, ADD_VALUES, ierr)
           end do
        end do

        ! Add 0 on diagonal to keep PETSc happy
        call MatSetValue(matrix, matrixSize-1, matrixSize-1, zero, ADD_VALUES, ierr)
     end if

  end if

  call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)

  if (constraint_option==2) then
     ! The matrix has a 1D null space, with the null vector corresponding to a constant in the L=0 component:
     call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrixSize, null_vec(1), ierr)
     call VecSet(null_vec(1), zero, ierr)
     if (masterProc) then
        do itheta = 1,Ntheta
           do izeta = 1,Nzeta
              index = getIndex(itheta,izeta,1)
              call VecSetValue(null_vec(1), index, one, INSERT_VALUES, ierr)
           end do
        end do        
     end if
     call VecAssemblyBegin(null_vec(1), ierr)
     call VecAssemblyEnd(null_vec(1), ierr)

     call MatNullSpaceCreate(MPI_COMM_WORLD,PETSC_FALSE,1,null_vec,nullspace, ierr)
!!$  call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr) ! This line creates a constant null space including all xi
     call MatSetNullSpace(matrix,nullspace,ierr)
     call MatNullSpaceDestroy(nullspace,ierr)
  end if

!!$  !write (filename,fmt="(a,i1,a)") "mmc_matrix_",whichMatrix,".dat"
!!$  filename = "mmc_matrix.dat"
!!$  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), FILE_MODE_WRITE, viewer, ierr)
!!$  call MatView(matrix, viewer, ierr)
!!$  call PetscViewerDestroy(viewer, ierr)


end subroutine populateMatrix


