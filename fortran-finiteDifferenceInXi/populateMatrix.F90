#include <finclude/petscmatdef.h>

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
  PetscInt :: rowIndex, colIndex
  PetscScalar, dimension(:,:), allocatable :: ddthetaToUse, ddzetaToUse
  PetscScalar :: temp

  if (masterProc) then
     print *,"Entering populateMatrix for whichMatrix = ",whichMatrix
  end if

  allocate(ddthetaToUse(Ntheta,Ntheta))
  allocate(ddzetaToUse(Nzeta,Nzeta))

  if (masterProc) then
     ! For simplicity, do all work on the master proc.

     !if (whichMatrix==0) then
     if (.false.) then
        print *,"theta:"
        print *,theta
        print *,"zeta:"
        print *,zeta
        print *,"xi:"
        print *,xi
        print *,"ddtheta_positiveXi:"
        do itheta=1,Ntheta
           print *,ddtheta_positiveXi(itheta,:)
        end do
        print *,"ddtheta_preconditioner_positiveXi:"
        do itheta=1,Ntheta
           print *,ddtheta_preconditioner_positiveXi(itheta,:)
        end do
        print *,"ddtheta_preconditioner_negativeXi:"
        do itheta=1,Ntheta
           print *,ddtheta_preconditioner_negativeXi(itheta,:)
        end do
        print *,"ddzeta_negativeXi:"
        do izeta=1,Nzeta
           print *,ddzeta_negativeXi(izeta,:)
        end do
        print *,"ddzeta_preconditioner_positiveXi:"
        do izeta=1,Nzeta
           print *,ddzeta_preconditioner_positiveXi(izeta,:)
        end do
        print *,"ddzeta_preconditioner_negativeXi:"
        do izeta=1,Nzeta
           print *,ddzeta_preconditioner_negativeXi(izeta,:)
        end do
        print *,"B:"
        do itheta=1,Ntheta
           print *,B(itheta,:)
        end do
     end if

     ! Add d/dtheta parallel streaming term
     itheta = -1
     izetaRow = -1
     izetaCol = -1
     do ixi = 1,Nxi
        if (whichMatrix == 1) then
           if (iota*xi(ixi)>0) then
              ddthetaToUse = ddtheta_positiveXi
           else
              ddthetaToUse = ddtheta_negativeXi
           end if
        else
           if (iota*xi(ixi)>0) then
              ddthetaToUse = ddtheta_preconditioner_positiveXi
           else
              ddthetaToUse = ddtheta_preconditioner_negativeXi
           end if
        end if
        do ithetaRow = 1,Ntheta
           do izeta = 1,Nzeta
              rowIndex = getIndex(ithetaRow, izeta, ixi)
              do ithetaCol = 1,Ntheta
                 colIndex = getIndex(ithetaCol, izeta, ixi)
                 call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      iota*B(ithetaRow,izeta)*xi(ixi)*ddthetaToUse(ithetaRow,ithetaCol), ADD_VALUES, ierr)
              end do
           end do
        end do
     end do

     ! Add d/dzeta parallel streaming term
     izeta = -1
     ithetaRow = -1
     ithetaCol = -1
     do ixi = 1,Nxi
        if (whichMatrix == 1) then
           if (xi(ixi)>0) then
              ddzetaToUse = ddzeta_positiveXi
           else
              ddzetaToUse = ddzeta_negativeXi
           end if
        else
           if (xi(ixi)>0) then
              ddzetaToUse = ddzeta_preconditioner_positiveXi
           else
              ddzetaToUse = ddzeta_preconditioner_negativeXi
           end if
        end if
        do izetaRow = 1,Nzeta
           do itheta = 1,Ntheta
              rowIndex = getIndex(itheta, izetaRow, ixi)
              do izetaCol = 1,Nzeta
                 colIndex = getIndex(itheta, izetaCol, ixi)
                 call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      B(itheta,izetaRow)*xi(ixi)*ddzetaToUse(izetaRow,izetaCol), ADD_VALUES, ierr)
              end do
           end do
        end do
     end do

     ! Add mirror term and collision operator
     izetaRow = -1
     izetaCol = -1
     ixi = -1
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           do ixiRow = 1,Nxi
              rowIndex = getIndex(itheta,izeta,ixiRow)
              temp = -(0.5d+0)*(1-xi(ixiRow)*xi(ixiRow))*(dBdtheta(itheta,izeta)*iota + dBdzeta(itheta,izeta))
              do ixiCol = 1,Nxi
                 colIndex = getIndex(itheta,izeta,ixiCol)
                 call MatSetValueSparse(matrix, rowIndex, colIndex, &
                      temp*ddxi(ixiRow,ixiCol) &  ! Mirror term
                      - nu * LorentzOperator(ixiRow,ixiCol), &    ! Collision operator
                      ADD_VALUES, ierr)
              end do
           end do
        end do
     end do

  end if

  call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)

  ! The matrix has a 1D null space, with the null vector corresponding to a constant:
!  call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
!  call MatSetNullSpace(matrix,nullspace,ierr)
!  call MatNullSpaceDestroy(nullspace,ierr)

  write (filename,fmt="(a,i1,a)") "mmc_matrix_",whichMatrix,".dat"
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), FILE_MODE_WRITE, viewer, ierr)
  call MatView(matrix, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)

  deallocate(ddthetaToUse, ddzetaToUse)




end subroutine populateMatrix


