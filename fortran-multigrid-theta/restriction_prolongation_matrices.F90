#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscmatdef.h>
#endif

subroutine restriction_prolongation_matrices(fine_level)

  use petscmat
  use indices
  use variables, only: levels, multigrid_restriction_matrices, multigrid_prolongation_matrices, numProcs, one, constraint_option, &
       restriction_option, pi, zetaMax, masterProc

  implicit none

  integer, intent(in) :: fine_level
  integer :: coarse_level
  integer :: Ntheta_fine, Ntheta_coarse
  integer :: Nzeta_fine, Nzeta_coarse
  integer :: Nxi_fine, Nxi_coarse
  PetscScalar, allocatable, dimension(:,:) :: theta_prolongation, zeta_prolongation, xi_prolongation
  integer :: itheta_fine, itheta_coarse
  integer :: izeta_fine, izeta_coarse
  integer :: ixi_fine, ixi_coarse
  PetscScalar :: theta_value, zeta_value, xi_value
  PetscErrorCode :: ierr
  Vec :: row_sums
  integer :: num_rows, num_cols, j, nnz_per_row


  coarse_level = fine_level + 1
  Ntheta_fine = levels(fine_level)%Ntheta
  Nzeta_fine  = levels(fine_level)%Nzeta
  Nxi_fine    = levels(fine_level)%Nxi
  Ntheta_coarse = levels(coarse_level)%Ntheta
  Nzeta_coarse  = levels(coarse_level)%Nzeta
  Nxi_coarse    = levels(coarse_level)%Nxi

  if (masterProc) print "(a,i4,a,i4)"," Assembling prolongation and restriction matrices between levels ",fine_level," and ",coarse_level
  
  ! First, generate prolongation matrices for each of the coordinates individually:
  allocate(theta_prolongation(Ntheta_fine, Ntheta_coarse))
  allocate(zeta_prolongation(Nzeta_fine, Nzeta_coarse))
  allocate(xi_prolongation(Nxi_fine, Nxi_coarse))

  call    periodic_interpolation(Ntheta_coarse, Ntheta_fine, 2*pi,                    levels(fine_level)%theta, theta_prolongation)
  call    periodic_interpolation( Nzeta_coarse,  Nzeta_fine, zetaMax,                 levels(fine_level)%zeta,   zeta_prolongation)
  call nonperiodic_interpolation(   Nxi_coarse,    Nxi_fine, levels(coarse_level)%xi, levels(fine_level)%xi,       xi_prolongation)

  if (masterProc) then
     print *,"Here comes theta interpolation matrix:"
     do j=1,Ntheta_fine
        print "(*(f10.6))",theta_prolongation(j,:)
     end do
     print *,"Here comes zeta interpolation matrix:"
     do j=1,Nzeta_fine
        print "(*(f10.6))",zeta_prolongation(j,:)
     end do
     print *,"Here comes xi interpolation matrix:"
     do j=1,Nxi_fine
        print "(*(f10.6))",xi_prolongation(j,:)
     end do
  end if

  ! Initialize the 3D prolongation matrix
  call MatCreate(PETSC_COMM_WORLD, multigrid_prolongation_matrices(fine_level), ierr)
  call MatSetType(multigrid_prolongation_matrices(fine_level), MATAIJ, ierr)
  call MatSetSizes(multigrid_prolongation_matrices(fine_level), PETSC_DECIDE, PETSC_DECIDE, &
       levels(fine_level)%matrixSize, levels(coarse_level)%matrixSize, ierr)
  print *,"HHH"

  ! Allocate 7 nonzero entries per row:
  ! This is a very conservative estimate: 2 off-diagonal elements in theta, zeta, and xi, plus the diagonal.
  ! This estimate could surely be made tighter if it mattered significantly, which is probably does not.
  nnz_per_row = 7
  if (numProcs == 1) then
     call MatSeqAIJSetPreallocation(multigrid_prolongation_matrices(fine_level), nnz_per_row, PETSC_NULL_INTEGER, ierr)
  else
     call MatMPIAIJSetPreallocation(multigrid_prolongation_matrices(fine_level), nnz_per_row, PETSC_NULL_INTEGER, nnz_per_row, PETSC_NULL_INTEGER, ierr)
  end if
  print *,"JJJ"
  ! Transfer the source/constraint element at the end by itself:
  if (constraint_option==1) then
     call MatSetValue(multigrid_prolongation_matrices(fine_level), levels(fine_level)%matrixSize, levels(coarse_level)%matrixSize, &
          one, ADD_VALUES, ierr)
  end if
  print *,"KKK"
  ! Now populate the 3D prolongation matrix.
  ! Parallelize the looping over rows (fine grid) but not the looping over columns (coarse grid)
  do itheta_coarse = 1,Ntheta_coarse
     do itheta_fine = levels(fine_level)%ithetaMin, levels(fine_level)%ithetaMax
        theta_value = theta_prolongation(itheta_fine, itheta_coarse)
        if (abs(theta_value)<1e-12) cycle
        do izeta_coarse = 1,Nzeta_coarse
           do izeta_fine = levels(fine_level)%izetaMin, levels(fine_level)%izetaMax
              zeta_value = zeta_prolongation(izeta_fine, izeta_coarse)
              if (abs(zeta_value)<1e-12) cycle
              do ixi_coarse = 1,Nxi_coarse
                 do ixi_fine = 1,Nxi_fine
                    xi_value = xi_prolongation(ixi_fine, ixi_coarse)
                    if (abs(xi_value)>1e-12) then
                       call MatSetValue(multigrid_prolongation_matrices(fine_level), &
                            getIndex(fine_level, itheta_fine, izeta_fine, ixi_fine), &
                            getIndex(coarse_level, itheta_coarse, izeta_coarse, ixi_coarse), &
                            theta_value*zeta_value*xi_value, ADD_VALUES, ierr)
                    end if
                 end do
              end do
           end do
        end do
     end do
  end do
  print *,"RRR"
  call MatAssemblyBegin(multigrid_prolongation_matrices(fine_level), MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(multigrid_prolongation_matrices(fine_level), MAT_FINAL_ASSEMBLY, ierr)
  print *,"TTT"
  call MatView(multigrid_prolongation_matrices(fine_level), PETSC_VIEWER_STDOUT_WORLD, ierr)
  ! Now handle the restriction matrix.
  call MatTranspose(multigrid_prolongation_matrices(fine_level), MAT_INITIAL_MATRIX, multigrid_restriction_matrices(fine_level), ierr)
  print *,"UUU"
  call MatCreateVecs(multigrid_restriction_matrices(fine_level), PETSC_NULL_OBJECT, row_sums, ierr)
  call VecGetSize(row_sums, num_rows, ierr)
  print *,"Vec initial num_rows:",num_rows
  call MatGetRowSum(multigrid_restriction_matrices(fine_level), row_sums, ierr)
  call MatGetSize(multigrid_restriction_matrices(fine_level), num_rows, num_cols, ierr)
  print *,"matrix num_rows:",num_rows,"num_cols:",num_cols
  print *,"VVV"
  select case (restriction_option)
  case (1)
     call VecReciprocal(row_sums, ierr)
     call VecGetSize(row_sums, num_rows, ierr)
     print *,"Vec final num_rows:",num_rows
     call MatDiagonalScale(multigrid_restriction_matrices(fine_level), row_sums, PETSC_NULL_OBJECT, ierr)
  case default
     print *,"Error! Invalid restriction_option:",restriction_option
     stop
  end select
  print *,"WWW"
  ! Clean up.
  call VecDestroy(row_sums, ierr)
  deallocate(theta_prolongation, zeta_prolongation, xi_prolongation)
  print *,"ZZZ"

end subroutine restriction_prolongation_matrices

