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
  PetscInt :: L, ell
  PetscInt :: rowIndex, colIndex, j, irow, index
  PetscScalar, dimension(:,:), pointer :: derivative_matrix_to_use
  PetscScalar :: temp, factor, L_factor, upwinding_scale_factor_to_use

  ! For convenience, use some short variable names to refer to quantities on this level:
  integer :: Ntheta, Nzeta, Nxi, matrixSize
  integer :: ithetaMin, ithetaMax, izetaMin, izetaMax
  PetscScalar, dimension(:), pointer :: thetaWeights, zetaWeights
  PetscScalar, dimension(:,:), pointer :: B, dBdtheta, dBdzeta
  PetscScalar, dimension(:,:), pointer :: ddtheta_sum_to_use, ddtheta_difference_to_use
  PetscScalar, dimension(:,:), pointer :: ddzeta_sum_to_use, ddzeta_difference_to_use
  
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
  thetaWeights => levels(level)%thetaWeights
  zetaWeights  => levels(level)%zetaWeights
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
  if (whichMatrix == 1) then
     ddtheta_sum_to_use => levels(level)%ddtheta_sum
     ddtheta_difference_to_use => levels(level)%ddtheta_difference
     upwinding_scale_factor_to_use = upwinding_scale_factor
  else
     ddtheta_sum_to_use => levels(level)%ddtheta_sum_preconditioner
     ddtheta_difference_to_use => levels(level)%ddtheta_difference_preconditioner
     upwinding_scale_factor_to_use = preconditioner_upwinding_scale_factor
  end if
  do L=0,(Nxi-1)
     do ithetaRow = ithetaMin,ithetaMax
        do izeta = izetaMin,izetaMax
           factor = iota * B(ithetaRow,izeta)
           rowIndex = getIndex(level, ithetaRow, izeta, L+1)
           do ithetaCol = 1,Ntheta
              
              ! Diagonal-in-L term:
              ell = L
              L_factor = (2*L*L+2*L-one)/((2*L+3)*(2*L-one)) * L_scaling(L+1) * f_scaling(ell+1)
              colIndex = getIndex(level, ithetaCol, izeta, ell+1)
              call MatSetValueSparse(matrix, rowIndex, colIndex, &
                   upwinding_scale_factor_to_use*L_factor*abs(factor)*ddtheta_difference_to_use(ithetaRow,ithetaCol), ADD_VALUES, ierr)
              ! Need the abs() around factor for the ddtheta_difference terms to handle the case iota<0.

              !if (whichMatrix .ne. 4) then
              if (.true.) then
                 ! Super-diagonal-in-L term
                 if (L < Nxi-1) then
                    ell = L + 1
                    L_factor = (L+1)/(two*L+3) * L_scaling(L+1) * f_scaling(ell+1)
                    colIndex = getIndex(level, ithetaCol, izeta, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         L_factor*factor*ddtheta_sum_to_use(ithetaRow,ithetaCol), ADD_VALUES, ierr)
                 end if
              
                 ! Sub-diagonal-in-L term
                 if (L > 0) then
                    ell = L - 1
                    L_factor = L/(2*L-one) * L_scaling(L+1) * f_scaling(ell+1)
                    colIndex = getIndex(level, ithetaCol, izeta, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         L_factor*factor*ddtheta_sum_to_use(ithetaRow,ithetaCol), ADD_VALUES, ierr)
                 end if
              end if
              
              !if (whichMatrix .ne. 4) then
              !if (whichMatrix>0) then
              if (.true.) then
                 ! Super-super-diagonal-in-L term
                 if (L < Nxi-2) then
                    ell = L + 2
                    L_factor = (L+two)*(L+one)/((two*L+5)*(two*L+3)) * L_scaling(L+1) * f_scaling(ell+1)
                    colIndex = getIndex(level, ithetaCol, izeta, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         upwinding_scale_factor_to_use*L_factor*abs(factor)*ddtheta_difference_to_use(ithetaRow,ithetaCol), ADD_VALUES, ierr)
                    ! Need the abs() around factor for the ddtheta_difference terms to handle the case iota<0.
                 end if
                 
                 ! Sub-sub-diagonal-in-L term
                 if (L > 1) then
                    ell = L - 2
                    L_factor = (L-one)*L/((two*L-3)*(two*L-one)) * L_scaling(L+1) * f_scaling(ell+1)
                    colIndex = getIndex(level, ithetaCol, izeta, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         upwinding_scale_factor_to_use*L_factor*abs(factor)*ddtheta_difference_to_use(ithetaRow,ithetaCol), ADD_VALUES, ierr)
                    ! Need the abs() around factor for the ddtheta_difference terms to handle the case iota<0.
                 end if
              end if
           end do
        end do
     end do
  end do

!!$  ! Add d/dtheta parallel streaming term
!!$  itheta = -1
!!$  izetaRow = -1
!!$  izetaCol = -1
!!$  do ixi = 1,Nxi
!!$     do ithetaRow = ithetaMin,ithetaMax
!!$        do izeta = izetaMin,izetaMax
!!$           factor = iota*B(ithetaRow,izeta)*xi(ixi) + E*iota*B(ithetaRow,izeta)*B(ithetaRow,izeta)/FSAB2
!!$           if (factor>0) then
!!$              if (whichMatrix == 1 .or. whichMatrix == 10) then
!!$                 derivative_matrix_to_use => levels(level)%ddtheta_plus
!!$              else
!!$                 derivative_matrix_to_use => levels(level)%ddtheta_plus_preconditioner
!!$              end if
!!$           else
!!$              if (whichMatrix == 1 .or. whichMatrix == 10) then
!!$                 derivative_matrix_to_use => levels(level)%ddtheta_minus
!!$              else
!!$                 derivative_matrix_to_use => levels(level)%ddtheta_minus_preconditioner
!!$              end if
!!$           end if
!!$           rowIndex = getIndex(level,ithetaRow, izeta, ixi)
!!$           do ithetaCol = 1,Ntheta
!!$              if (whichMatrix==5 .and. ithetaCol .ne. ithetaRow) cycle ! For whichMatrix==5, add only the diagonal.
!!$              if (whichMatrix==6 .and. ithetaCol == ithetaRow) cycle   ! For whichMatrix==6, add only the off-diagonal elements.
!!$              colIndex = getIndex(level,ithetaCol, izeta, ixi)
!!$              call MatSetValueSparse(matrix, rowIndex, colIndex, &
!!$                   factor*derivative_matrix_to_use(ithetaRow,ithetaCol), ADD_VALUES, ierr)
!!$           end do
!!$        end do
!!$     end do
!!$  end do
  
  
  ! Add d/dzeta parallel streaming term
  if (whichMatrix == 1) then
     ddzeta_sum_to_use => levels(level)%ddzeta_sum
     ddzeta_difference_to_use => levels(level)%ddzeta_difference
  else
     ddzeta_sum_to_use => levels(level)%ddzeta_sum_preconditioner
     ddzeta_difference_to_use => levels(level)%ddzeta_difference_preconditioner
  end if
  do L=0,(Nxi-1)     
     do izetaRow = izetaMin,izetaMax
        do itheta = ithetaMin,ithetaMax
           factor = B(itheta,izetaRow)
           rowIndex = getIndex(level, itheta, izetaRow, L+1)
           do izetaCol = 1,Nzeta
              
              ! Diagonal-in-L term:
              ell = L
              L_factor = (2*L*L+2*L-one)/((2*L+3)*(2*L-one)) * L_scaling(L+1) * f_scaling(ell+1)
              colIndex = getIndex(level, itheta, izetaCol, ell+1)
              call MatSetValueSparse(matrix, rowIndex, colIndex, &
                   upwinding_scale_factor_to_use*L_factor*factor*ddzeta_difference_to_use(izetaRow,izetaCol), ADD_VALUES, ierr)
              
              !if (whichMatrix .ne. 4) then
              if (.true.) then
                 ! Super-diagonal-in-L term
                 if (L < Nxi-1) then
                    ell = L + 1
                    L_factor = (L+1)/(two*L+3) * L_scaling(L+1) * f_scaling(ell+1)
                    colIndex = getIndex(level, itheta, izetaCol, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         L_factor*factor*ddzeta_sum_to_use(izetaRow,izetaCol), ADD_VALUES, ierr)
                 end if
                 
                 ! Sub-diagonal-in-L term
                 if (L > 0) then
                    ell = L - 1
                    L_factor = L/(2*L-one) * L_scaling(L+1) * f_scaling(ell+1)
                    colIndex = getIndex(level, itheta, izetaCol, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         L_factor*factor*ddzeta_sum_to_use(izetaRow,izetaCol), ADD_VALUES, ierr)
                 end if
              end if
              
              !if (whichMatrix .ne. 4) then
              !if (whichMatrix>0) then
              if (.true.) then
                 ! Super-super-diagonal-in-L term
                 if (L < Nxi-2) then
                    ell = L + 2
                    L_factor = (L+two)*(L+one)/((two*L+5)*(two*L+3)) * L_scaling(L+1) * f_scaling(ell+1)
                    colIndex = getIndex(level, itheta, izetaCol, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         upwinding_scale_factor_to_use*L_factor*factor*ddzeta_difference_to_use(izetaRow,izetaCol), ADD_VALUES, ierr)
                 end if
                 
                 ! Sub-sub-diagonal-in-L term
                 if (L > 1) then
                    ell = L - 2
                    L_factor = (L-one)*L/((two*L-3)*(two*L-one)) * L_scaling(L+1) * f_scaling(ell+1)
                    colIndex = getIndex(level, itheta, izetaCol, ell+1)
                    call MatSetValueSparse(matrix, rowIndex, colIndex, &
                         upwinding_scale_factor_to_use*L_factor*factor*ddzeta_difference_to_use(izetaRow,izetaCol), ADD_VALUES, ierr)
                 end if
              end if
           end do
        end do
     end do
  end do

!!$  ! Add d/dzeta parallel streaming term
!!$  if (whichMatrix .ne. 6) then
!!$     izeta = -1
!!$     ithetaRow = -1
!!$     ithetaCol = -1
!!$     do ixi = 1,Nxi
!!$        do izetaRow = izetaMin,izetaMax
!!$           do itheta = ithetaMin,ithetaMax
!!$              factor = B(itheta,izetaRow)*xi(ixi) - E*iota*I*B(itheta,izetaRow)*B(itheta,izetaRow)/(G*FSAB2)
!!$              if (factor>0) then
!!$                 if (whichMatrix == 1) then
!!$                    derivative_matrix_to_use => levels(level)%ddzeta_plus
!!$                 else
!!$                    derivative_matrix_to_use => levels(level)%ddzeta_plus_preconditioner
!!$                 end if
!!$              else
!!$                 if (whichMatrix == 1) then
!!$                    derivative_matrix_to_use => levels(level)%ddzeta_minus
!!$                 else
!!$                    derivative_matrix_to_use => levels(level)%ddzeta_minus_preconditioner
!!$                 end if
!!$              end if
!!$              rowIndex = getIndex(level,itheta, izetaRow, ixi)
!!$              do izetaCol = 1,Nzeta
!!$                 colIndex = getIndex(level,itheta, izetaCol, ixi)
!!$                 call MatSetValueSparse(matrix, rowIndex, colIndex, &
!!$                      factor*derivative_matrix_to_use(izetaRow,izetaCol), ADD_VALUES, ierr)
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$  end if


!!$  ! Add mirror term
!!$  izetaRow = -1
!!$  izetaCol = -1
!!$  ixi = -1
!!$  do itheta = ithetaMin,ithetaMax
!!$     do izeta = izetaMin,izetaMax
!!$        do ixiRow = 1,Nxi
!!$           rowIndex = getIndex(level,itheta,izeta,ixiRow)
!!$           temp = -(0.5d+0)*(1-xi(ixiRow)*xi(ixiRow))*(dBdtheta(itheta,izeta)*iota + dBdzeta(itheta,izeta))
!!$           if (temp>0) then
!!$              if (whichMatrix == 1 .or. whichMatrix == 10) then
!!$                 derivative_matrix_to_use => levels(level)%ddxi_plus
!!$              else
!!$                 derivative_matrix_to_use => levels(level)%ddxi_plus_preconditioner
!!$              end if
!!$           else
!!$              if (whichMatrix == 1 .or. whichMatrix == 10) then
!!$                 derivative_matrix_to_use => levels(level)%ddxi_minus
!!$              else
!!$                 derivative_matrix_to_use => levels(level)%ddxi_minus_preconditioner
!!$              end if
!!$           end if
!!$           do ixiCol = 1,Nxi
!!$              if (whichMatrix==5 .and. ixiCol .ne. ixiRow) cycle ! For whichMatrix==5, add only the diagonal.
!!$              if (whichMatrix==6 .and. ixiCol  ==  ixiRow) cycle ! For whichMatrix==6, add only the off-diagonal elements.
!!$              colIndex = getIndex(level,itheta,izeta,ixiCol)
!!$              call MatSetValueSparse(matrix, rowIndex, colIndex, &
!!$                   temp * derivative_matrix_to_use(ixiRow,ixiCol) &  ! Mirror term
!!$                   - nu * pitch_angle_scattering_operator_to_use(ixiRow,ixiCol), &    ! Collision operator
!!$                   ADD_VALUES, ierr)
!!$           end do
!!$        end do
!!$     end do
!!$  end do

  ! Add mirror term                                                                                                                                                        
  izetaRow = -1
  izetaCol = -1
  do itheta = ithetaMin,ithetaMax
     do izeta = izetaMin,izetaMax
        temp = -(0.5d+0)*(dBdtheta(itheta,izeta)*iota + dBdzeta(itheta,izeta))
        do L = 0,(Nxi-1)
           rowIndex = getIndex(level,itheta,izeta,L+1)
           
           ! Super-diagonal term                                                                                                                                      
           if (L < Nxi-1) then
              ell = L + 1
              colIndex = getIndex(level,itheta,izeta,ell+1)
              call MatSetValueSparse(matrix, rowIndex, colIndex, &
                   temp*(L+1)*(L+2)/(2*L+3) * L_scaling(L+1) * f_scaling(ell+1), ADD_VALUES, ierr)
           end if
           
           ! Sub-diagonal term                                                                                                                           
           if (L > 0) then
              ell = L - 1
              colIndex = getIndex(level,itheta,izeta,ell+1)
              call MatSetValueSparse(matrix, rowIndex, colIndex, &
                   temp*(-L)*(L-1)/(2*L-1) * L_scaling(L+1) * f_scaling(ell+1), ADD_VALUES, ierr)
           end if
           
        end do
     end do
  end do
  
  ! Add collision operator                                                                                                                                  
  do itheta = ithetaMin,ithetaMax
     do izeta = izetaMin,izetaMax
        do L = 0,(Nxi-1)
           index = getIndex(level,itheta,izeta,L+1)
           ell = L
           temp = nu/2*L*(L+1) * L_scaling(L+1) * f_scaling(ell+1)
           if (whichMatrix==4 .and. L==0 .and. shift_L0_in_smoother) temp = nu * L_scaling(L+1) * f_scaling(ell+1)  ! These diagonal elements are otherwise small, so shift them a little for the smoother matrix.
           ! Let's use MatSetValue instead of MatSetValueSparse here to ensure the diagonal is set even if it is 0. 
           call MatSetValue(matrix, index, index, temp, ADD_VALUES, ierr)
        end do
     end do
  end do
  
  if (constraint_option==1 .and. masterProc .and. (whichMatrix==0 .or. whichMatrix==1 .or. whichMatrix == 10)) then
     ! Add source                                                                                                                                                          
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           index = getIndex(level,itheta,izeta,1)
           call MatSetValue(matrix, index, matrixSize-1, one, ADD_VALUES, ierr)
        end do
     end do

     ! Add constraint                                                                                                                                                      
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           index = getIndex(level,itheta,izeta,1)
           call MatSetValue(matrix, matrixSize-1, index, one, ADD_VALUES, ierr)
        end do
     end do

     ! Add 0 on diagonal to keep PETSc happy                                                                                                                               
     call MatSetValue(matrix, matrixSize-1, matrixSize-1, zero, ADD_VALUES, ierr)
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

  write (filename,fmt="(a,i1,a,i1,a)") "mmc_matrix_level_",level,'_whichMatrix_',whichMatrix,".txt"
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD, trim(filename), viewer, ierr)
  call MatView(matrix, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)




end subroutine populateMatrix


