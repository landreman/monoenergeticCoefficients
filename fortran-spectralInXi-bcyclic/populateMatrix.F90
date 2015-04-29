!#include <finclude/petscmatdef.h>

subroutine populateMatrix(matrix)

  !use petscmat

  use stel_kinds
  use indices
  use variables

  implicit none

  integer :: itheta, ithetaRow, ithetaCol
  integer :: izeta, izetaRow, izetaCol
  integer :: L, ell
  integer :: rowIndex, colIndex, index
  real(rprec) :: temp

  if (masterProc) then
     print *,"Entering populateMatrix."
  end if
  
  ! Add d/dtheta parallel streaming term
  itheta = -1
  izetaRow = -1
  izetaCol = -1
  do ithetaRow = 1,Ntheta
     do izeta = 1,Nzeta
        rowIndex = getBlockIndex(ithetaRow, izeta)
        if ((rowIndex >= ns0) .and. (rowIndex <= nsn)) then
           do L = 0,(Nxi-1)
              do ithetaCol = 1,Ntheta

                 ! Super-diagonal term
                 if (L < Nxi-1) then
                    ell = L + 1
                    colIndex = getBlockIndex(ithetaCol, izeta)
                    ublk(rowIndex, colIndex, L+1) = &
                      iota*B(ithetaRow,izeta)*ddtheta(ithetaRow,ithetaCol)*(L+1)/(2*L+3)
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



end subroutine populateMatrix


