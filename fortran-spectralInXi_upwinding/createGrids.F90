#include <petsc/finclude/petscsysdef.h>

subroutine createGrids()

  use petscsys

  use variables

  implicit none

  PetscScalar, dimension(:,:), allocatable :: d2dtheta2
  PetscScalar, dimension(:,:), allocatable :: d2dzeta2
  PetscScalar :: zetaMax
  logical :: call_uniform_diff_matrices
  integer :: derivative_option_plus, derivative_option_minus, quadrature_option
  integer :: j, k
  PetscScalar, dimension(:,:), allocatable :: ddtheta_plus, ddtheta_minus
  PetscScalar, dimension(:,:), allocatable :: ddzeta_plus, ddzeta_minus
  PetscScalar, dimension(:,:), allocatable :: ddtheta_plus_preconditioner, ddtheta_minus_preconditioner
  PetscScalar, dimension(:,:), allocatable :: ddzeta_plus_preconditioner, ddzeta_minus_preconditioner


  matrixSize = Ntheta*Nzeta*Nxi
  if (constraint_option==1) matrixSize = matrixSize + 1

  ! *****************************************************
  ! Set up theta grid
  ! *****************************************************

  allocate(theta(Ntheta))
  allocate(thetaWeights(Ntheta))
  allocate(ddtheta_plus(Ntheta,Ntheta))
  allocate(ddtheta_minus(Ntheta,Ntheta))
  allocate(ddtheta_plus_preconditioner(Ntheta,Ntheta))
  allocate(ddtheta_minus_preconditioner(Ntheta,Ntheta))
  allocate(d2dtheta2(Ntheta,Ntheta))

  select case (theta_derivative_option)

  case (1)
     if (masterProc) then
        print *,"d/dtheta derivative discretized using Fourier pseudospectral method."
     end if
     derivative_option_plus = 20
     derivative_option_minus = derivative_option_plus
     
  case (2)
     if (masterProc) then
        print *,"d/dtheta derivative discretized using centered finite differences:"
        print *,"   1 point on either side."
     end if
     derivative_option_plus = 0
     derivative_option_minus = derivative_option_plus
     
  case (3)
     if (masterProc) then
        print *,"d/dtheta derivative discretized using centered finite differences:"
        print *,"   2 points on either side."
     end if
     derivative_option_plus = 10
     derivative_option_minus = derivative_option_plus
     
  case (4)
     if (masterProc) then
        print *,"d/dtheta derivative discretized using upwinded finite differences:"
        print *,"   0 points on one side, 1 point on the other side."
     end if
     derivative_option_plus  = 30
     derivative_option_minus = 40
     
  case (5)
     if (masterProc) then
        print *,"d/dtheta derivative discretized using upwinded finite differences:"
        print *,"   0 points on one side, 2 points on the other side."
     end if
     derivative_option_plus  = 50
     derivative_option_minus = 60
     
  case (6)
     if (masterProc) then
        print *,"d/dtheta derivative discretized using upwinded finite differences:"
        print *,"   1 point on one side, 2 points on the other side."
     end if
     derivative_option_plus  = 80
     derivative_option_minus = 90
     
  case (7)
     if (masterProc) then
        print *,"d/dtheta derivative discretized using upwinded finite differences:"
        print *,"   1 point on one side, 3 points on the other side."
     end if
     derivative_option_plus  = 100
     derivative_option_minus = 110
     
  case (8)
     if (masterProc) then
        print *,"d/dtheta derivative discretized using upwinded finite differences:"
        print *,"   2 point on one side, 3 points on the other side."
     end if
     derivative_option_plus  = 120
     derivative_option_minus = 130
     
  case default
     if (masterProc) then
        print *,"Error! Invalid setting for theta_derivative_option:",theta_derivative_option
     end if
     stop
  end select
  
  quadrature_option = 0
  call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_plus,  quadrature_option, theta, thetaWeights, ddtheta_plus,  d2dtheta2)
  call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_minus, quadrature_option, theta, thetaWeights, ddtheta_minus, d2dtheta2)
  
  ddtheta_sum        = (ddtheta_plus + ddtheta_minus)/two
  ddtheta_difference = (ddtheta_plus - ddtheta_minus)/two
  
  ! *******************************************************************************
  ! Handle d/dtheta for preconditioner:
  ! *******************************************************************************
  
  call_uniform_diff_matrices = .true.
  select case (abs(preconditioner_theta_derivative_option))
  case (0)
     if (masterProc) then
        print *,"d/dtheta term is dropped in the preconditioner."
     end if
     ddtheta_plus_preconditioner = 0
     ddtheta_minus_preconditioner = 0
     call_uniform_diff_matrices = .false.
     
  case (100)
     if (masterProc) then
        print *,"d/dtheta matrix is the same in the preconditioner."
     end if
     ddtheta_plus_preconditioner  = ddtheta_plus
     ddtheta_minus_preconditioner = ddtheta_minus
     call_uniform_diff_matrices = .false.
     
  case (1)
     if (masterProc) then
        print *,"Preconditioner d/dtheta derivative discretized using Fourier pseudospectral method."
     end if
     derivative_option_plus = 20
     derivative_option_minus = derivative_option_plus
     
  case (2)
     if (masterProc) then
        print *,"Preconditioner d/dtheta derivative discretized using centered finite differences:"
        print *,"   1 point on either side."
     end if
     derivative_option_plus = 0
     derivative_option_minus = derivative_option_plus
     
  case (3)
     if (masterProc) then
        print *,"Preconditioner d/dtheta derivative discretized using centered finite differences:"
        print *,"   2 points on either side."
     end if
     derivative_option_plus = 10
     derivative_option_minus = derivative_option_plus
     
  case (4)
     if (masterProc) then
        print *,"Preconditioner d/dtheta derivative discretized using upwinded finite differences:"
        print *,"   0 points on one side, 1 point on the other side."
     end if
     derivative_option_plus  = 30
     derivative_option_minus = 40
     
  case (5)
     if (masterProc) then
        print *,"Preconditioner d/dtheta derivative discretized using upwinded finite differences:"
        print *,"   0 points on one side, 2 points on the other side."
     end if
     derivative_option_plus  = 50
     derivative_option_minus = 60
     
  case (6)
     if (masterProc) then
        print *,"Preconditioner d/dtheta derivative discretized using upwinded finite differences:"
        print *,"   1 point on one side, 2 points on the other side."
     end if
     derivative_option_plus  = 80
     derivative_option_minus = 90
     
  case (7)
     if (masterProc) then
        print *,"Preconditioner d/dtheta derivative discretized using upwinded finite differences:"
        print *,"   1 point on one side, 3 points on the other side."
     end if
     derivative_option_plus  = 100
     derivative_option_minus = 110
     
  case (8)
     if (masterProc) then
        print *,"Preconditioner d/dtheta derivative discretized using upwinded finite differences:"
        print *,"   2 point on one side, 3 points on the other side."
     end if
     derivative_option_plus  = 120
     derivative_option_minus = 130
     
  case default
     if (masterProc) then
        print *,"Error! Invalid setting for preconditioner_theta_derivative_option:",preconditioner_theta_derivative_option
     end if
     stop
  end select
  
  if (call_uniform_diff_matrices) then
     quadrature_option = 0
     call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_plus,  quadrature_option, theta, thetaWeights, ddtheta_plus_preconditioner,  d2dtheta2)
     call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_minus, quadrature_option, theta, thetaWeights, ddtheta_minus_preconditioner, d2dtheta2)
  end if
  
  if (preconditioner_theta_derivative_option<0) then
     if (masterProc) then
        print *,"   But only the diagonal is kept."
     end if
     do j=1,Ntheta
        do k=1,Ntheta
           if (j .ne. k) then
              ddtheta_plus_preconditioner(j,k) = 0
              ddtheta_minus_preconditioner(j,k) = 0
           end if
        end do
     end do
  end if
  
  ddtheta_sum_preconditioner        = (ddtheta_plus_preconditioner + ddtheta_minus_preconditioner)/two
  ddtheta_difference_preconditioner = (ddtheta_plus_preconditioner - ddtheta_minus_preconditioner)/two
  
  !  call uniformDiffMatrices(Ntheta, 0.0d+0, 2*pi, thetaGridScheme, theta, thetaWeights, ddtheta, d2dtheta2)
  
  deallocate(d2dtheta2)

  ! *****************************************************
  ! Set up zeta grid
  ! *****************************************************

  allocate(zeta(Nzeta))
  allocate(zetaWeights(Nzeta))
  allocate(ddzeta_plus(Nzeta,Nzeta))
  allocate(ddzeta_minus(Nzeta,Nzeta))
  allocate(ddzeta_plus_preconditioner(Nzeta,Nzeta))
  allocate(ddzeta_minus_preconditioner(Nzeta,Nzeta))
  allocate(d2dzeta2(Nzeta,Nzeta))

  zetaMax = 2*pi/Nperiods


  select case (zeta_derivative_option)

  case (1)
     if (masterProc) then
        print *,"d/dzeta derivative discretized using Fourier pseudospectral method."
     end if
     derivative_option_plus = 20
     derivative_option_minus = derivative_option_plus
     
  case (2)
     if (masterProc) then
        print *,"d/dzeta derivative discretized using centered finite differences:"
        print *,"   1 point on either side."
     end if
     derivative_option_plus = 0
     derivative_option_minus = derivative_option_plus
     
  case (3)
     if (masterProc) then
        print *,"d/dzeta derivative discretized using centered finite differences:"
        print *,"   2 points on either side."
     end if
     derivative_option_plus = 10
     derivative_option_minus = derivative_option_plus
     
  case (4)
     if (masterProc) then
        print *,"d/dzeta derivative discretized using upwinded finite differences:"
        print *,"   0 points on one side, 1 point on the other side."
     end if
     derivative_option_plus  = 30
     derivative_option_minus = 40
     
  case (5)
     if (masterProc) then
        print *,"d/dzeta derivative discretized using upwinded finite differences:"
        print *,"   0 points on one side, 2 points on the other side."
     end if
     derivative_option_plus  = 50
     derivative_option_minus = 60
     
  case (6)
     if (masterProc) then
        print *,"d/dzeta derivative discretized using upwinded finite differences:"
        print *,"   1 point on one side, 2 points on the other side."
     end if
     derivative_option_plus  = 80
     derivative_option_minus = 90
     
  case (7)
     if (masterProc) then
        print *,"d/dzeta derivative discretized using upwinded finite differences:"
        print *,"   1 point on one side, 3 points on the other side."
     end if
     derivative_option_plus  = 100
     derivative_option_minus = 110
     
  case (8)
     if (masterProc) then
        print *,"d/dzeta derivative discretized using upwinded finite differences:"
        print *,"   2 point on one side, 3 points on the other side."
     end if
     derivative_option_plus  = 120
     derivative_option_minus = 130
     
  case default
     if (masterProc) then
        print *,"Error! Invalid setting for zeta_derivative_option:",zeta_derivative_option
     end if
     stop
  end select
  
  quadrature_option = 0
  call uniformDiffMatrices(Nzeta, zero, zetaMax, derivative_option_plus,  quadrature_option, zeta, zetaWeights, ddzeta_plus,  d2dzeta2)
  call uniformDiffMatrices(Nzeta, zero, zetaMax, derivative_option_minus, quadrature_option, zeta, zetaWeights, ddzeta_minus, d2dzeta2)
  
  ddzeta_sum        = (ddzeta_plus + ddzeta_minus)/two
  ddzeta_difference = (ddzeta_plus - ddzeta_minus)/two
  
  ! *******************************************************************************
  ! Handle d/dzeta for preconditioner:
  ! *******************************************************************************
  
  call_uniform_diff_matrices = .true.
  select case (abs(preconditioner_zeta_derivative_option))
  case (0)
     if (masterProc) then
        print *,"d/dzeta term is dropped in the preconditioner."
     end if
     ddzeta_plus_preconditioner = 0
     ddzeta_minus_preconditioner = 0
     call_uniform_diff_matrices = .false.
     
  case (100)
     if (masterProc) then
        print *,"d/dzeta matrix is the same in the preconditioner."
     end if
     ddzeta_plus_preconditioner  = ddzeta_plus
     ddzeta_minus_preconditioner = ddzeta_minus
     call_uniform_diff_matrices = .false.
     
  case (1)
     if (masterProc) then
        print *,"Preconditioner d/dzeta derivative discretized using Fourier pseudospectral method."
     end if
     derivative_option_plus = 20
     derivative_option_minus = derivative_option_plus
     
  case (2)
     if (masterProc) then
        print *,"Preconditioner d/dzeta derivative discretized using centered finite differences:"
        print *,"   1 point on either side."
     end if
     derivative_option_plus = 0
     derivative_option_minus = derivative_option_plus
     
  case (3)
     if (masterProc) then
        print *,"Preconditioner d/dzeta derivative discretized using centered finite differences:"
        print *,"   2 points on either side."
     end if
     derivative_option_plus = 10
     derivative_option_minus = derivative_option_plus
     
  case (4)
     if (masterProc) then
        print *,"Preconditioner d/dzeta derivative discretized using upwinded finite differences:"
        print *,"   0 points on one side, 1 point on the other side."
     end if
     derivative_option_plus  = 30
     derivative_option_minus = 40
     
  case (5)
     if (masterProc) then
        print *,"Preconditioner d/dzeta derivative discretized using upwinded finite differences:"
        print *,"   0 points on one side, 2 points on the other side."
     end if
     derivative_option_plus  = 50
     derivative_option_minus = 60
     
  case (6)
     if (masterProc) then
        print *,"Preconditioner d/dzeta derivative discretized using upwinded finite differences:"
        print *,"   1 point on one side, 2 points on the other side."
     end if
     derivative_option_plus  = 80
     derivative_option_minus = 90
     
  case (7)
     if (masterProc) then
        print *,"Preconditioner d/dzeta derivative discretized using upwinded finite differences:"
        print *,"   1 point on one side, 3 points on the other side."
     end if
     derivative_option_plus  = 100
     derivative_option_minus = 110
     
  case (8)
     if (masterProc) then
        print *,"Preconditioner d/dzeta derivative discretized using upwinded finite differences:"
        print *,"   2 point on one side, 3 points on the other side."
     end if
     derivative_option_plus  = 120
     derivative_option_minus = 130
     
  case default
     if (masterProc) then
        print *,"Error! Invalid setting for preconditioner_zeta_derivative_option:",preconditioner_zeta_derivative_option
     end if
     stop
  end select
  
  if (call_uniform_diff_matrices) then
     quadrature_option = 0
     call uniformDiffMatrices(Nzeta, zero, zetaMax, derivative_option_plus,  quadrature_option, zeta, zetaWeights, ddzeta_plus_preconditioner,  d2dzeta2)
     call uniformDiffMatrices(Nzeta, zero, zetaMax, derivative_option_minus, quadrature_option, zeta, zetaWeights, ddzeta_minus_preconditioner, d2dzeta2)
  end if
  
  if (preconditioner_zeta_derivative_option<0) then
     if (masterProc) then
        print *,"   But only the diagonal is kept."
     end if
     do j=1,Nzeta
        do k=1,Nzeta
           if (j .ne. k) then
              ddzeta_plus_preconditioner(j,k) = 0
              ddzeta_minus_preconditioner(j,k) = 0
           end if
        end do
     end do
  end if
  
  ddzeta_sum_preconditioner        = (ddzeta_plus_preconditioner + ddzeta_minus_preconditioner)/two
  ddzeta_difference_preconditioner = (ddzeta_plus_preconditioner - ddzeta_minus_preconditioner)/two
  
  !call uniformDiffMatrices(Nzeta, 0.0d+0, zetaMax, zetaGridScheme, zeta, zetaWeights, ddzeta, d2dzeta2)

  deallocate(d2dzeta2)

  zetaWeights = zetaWeights * Nperiods

  ! *****************************************************
  ! Compute B and its derivatives
  ! *****************************************************

  allocate(B(Ntheta,Nzeta))
  allocate(dBdtheta(Ntheta,Nzeta))
  allocate(dBdzeta(Ntheta,Nzeta))
  call computeB()

end subroutine createGrids


