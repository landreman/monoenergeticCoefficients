#include <finclude/petscsysdef.h>

subroutine createGrids()

  use petscsys

  use variables, only: Nperiods, pi, matrixSize, &
       B, dBdtheta, dBdzeta, &
       Ntheta, Nzeta, Nxi, theta, zeta, xi, &
       thetaWeights, zetaWeights, xiWeights, &
       ddtheta_positiveXi, ddtheta_negativeXi, ddzeta_positiveXi, ddzeta_negativeXi, &
       ddtheta_preconditioner_positiveXi, ddtheta_preconditioner_negativeXi, &
       ddzeta_preconditioner_positiveXi, ddzeta_preconditioner_negativeXi, &
       ddxi, LorentzOperator, &
       thetaGridScheme, thetaGridScheme_pc, zetaGridScheme, zetaGridScheme_pc

  implicit none

  PetscScalar, dimension(:), allocatable :: theta_preconditioner, thetaWeights_preconditioner
  PetscScalar, dimension(:), allocatable :: zeta_preconditioner, zetaWeights_preconditioner
  PetscScalar, dimension(:,:), allocatable :: d2dtheta2
  PetscScalar, dimension(:,:), allocatable :: d2dzeta2
  PetscScalar, dimension(:,:), allocatable :: d2dxi2
  PetscInt :: scheme, ixi
  PetscScalar :: zetaMax

  matrixSize = Ntheta*Nzeta*Nxi

  ! *****************************************************
  ! Set up theta grid
  ! *****************************************************

  allocate(theta(Ntheta))
  allocate(theta_preconditioner(Ntheta))
  allocate(thetaWeights(Ntheta))
  allocate(thetaWeights_preconditioner(Ntheta))
  allocate(ddtheta_positiveXi(Ntheta,Ntheta))
  allocate(ddtheta_negativeXi(Ntheta,Ntheta))
  allocate(ddtheta_preconditioner_positiveXi(Ntheta,Ntheta))
  allocate(ddtheta_preconditioner_negativeXi(Ntheta,Ntheta))
  allocate(d2dtheta2(Ntheta,Ntheta))

  print *,"000"
  select case (thetaGridScheme)
  case (1)
     scheme = 0
     call uniformDiffMatrices(Ntheta, 0, 2*pi, scheme, theta, thetaWeights, ddtheta_positiveXi, d2dtheta2)
     ddtheta_negativeXi = ddtheta_positiveXi
  case (2)
     scheme = 10
     print *,"123"
     call uniformDiffMatrices(Ntheta, 0, 2*pi, scheme, theta, thetaWeights, ddtheta_positiveXi, d2dtheta2)
     print *,"456"
     ddtheta_negativeXi = ddtheta_positiveXi
  case (3)
     scheme = 30
     call uniformDiffMatrices(Ntheta, 0, 2*pi, scheme, theta, thetaWeights, ddtheta_positiveXi, d2dtheta2)
     scheme = 40
     call uniformDiffMatrices(Ntheta, 0, 2*pi, scheme, theta, thetaWeights, ddtheta_negativeXi, d2dtheta2)
  case default
     stop "Error! Invalid thetaGridScheme!"
  end select
  print *,"AAA"
  select case (thetaGridScheme_pc)
  case (1)
     scheme = 0
     call uniformDiffMatrices(Ntheta, 0, 2*pi, scheme, theta_preconditioner, thetaWeights_preconditioner, ddtheta_preconditioner_positiveXi, d2dtheta2)
     ddtheta_preconditioner_negativeXi = ddtheta_preconditioner_positiveXi
  case (2)
     scheme = 10
     call uniformDiffMatrices(Ntheta, 0, 2*pi, scheme, theta_preconditioner, thetaWeights_preconditioner, ddtheta_preconditioner_positiveXi, d2dtheta2)
     ddtheta_preconditioner_negativeXi = ddtheta_preconditioner_positiveXi
  case (3)
     scheme = 30
     call uniformDiffMatrices(Ntheta, 0, 2*pi, scheme, theta_preconditioner, thetaWeights_preconditioner, ddtheta_preconditioner_positiveXi, d2dtheta2)
     scheme = 40
     call uniformDiffMatrices(Ntheta, 0, 2*pi, scheme, theta_preconditioner, thetaWeights_preconditioner, ddtheta_preconditioner_negativeXi, d2dtheta2)
  case default
     stop "Error! Invalid thetaGridScheme_pc!"
  end select

  deallocate(theta_preconditioner, thetaWeights_preconditioner, d2dtheta2)
  print *,"BBB"
  ! *****************************************************
  ! Set up zeta grid
  ! *****************************************************

  allocate(zeta(Nzeta))
  allocate(zeta_preconditioner(Nzeta))
  allocate(zetaWeights(Nzeta))
  allocate(zetaWeights_preconditioner(Nzeta))
  allocate(ddzeta_positiveXi(Nzeta,Nzeta))
  allocate(ddzeta_negativeXi(Nzeta,Nzeta))
  allocate(ddzeta_preconditioner_positiveXi(Nzeta,Nzeta))
  allocate(ddzeta_preconditioner_negativeXi(Nzeta,Nzeta))
  allocate(d2dzeta2(Nzeta,Nzeta))

  zetaMax = 2*pi/Nperiods

  select case (zetaGridScheme)
  case (1)
     scheme = 0
     call uniformDiffMatrices(Nzeta, 0, zetaMax, scheme, zeta, zetaWeights, ddzeta_positiveXi, d2dzeta2)
     ddzeta_negativeXi = ddzeta_positiveXi
  case (2)
     scheme = 10
     call uniformDiffMatrices(Nzeta, 0, zetaMax, scheme, zeta, zetaWeights, ddzeta_positiveXi, d2dzeta2)
     ddzeta_negativeXi = ddzeta_positiveXi
  case (3)
     scheme = 30
     call uniformDiffMatrices(Nzeta, 0, zetaMax, scheme, zeta, zetaWeights, ddzeta_positiveXi, d2dzeta2)
     scheme = 40
     call uniformDiffMatrices(Nzeta, 0, zetaMax, scheme, zeta, zetaWeights, ddzeta_negativeXi, d2dzeta2)
  case default
     stop "Error! Invalid zetaGridScheme!"
  end select
  print *,"CCC"
  select case (zetaGridScheme_pc)
  case (1)
     scheme = 0
     call uniformDiffMatrices(Nzeta, 0, zetaMax, scheme, zeta_preconditioner, zetaWeights_preconditioner, ddzeta_preconditioner_positiveXi, d2dzeta2)
     ddzeta_preconditioner_negativeXi = ddzeta_preconditioner_positiveXi
  case (2)
     scheme = 10
     call uniformDiffMatrices(Nzeta, 0, zetaMax, scheme, zeta_preconditioner, zetaWeights_preconditioner, ddzeta_preconditioner_positiveXi, d2dzeta2)
     ddzeta_preconditioner_negativeXi = ddzeta_preconditioner_positiveXi
  case (3)
     scheme = 30
     call uniformDiffMatrices(Nzeta, 0, zetaMax, scheme, zeta_preconditioner, zetaWeights_preconditioner, ddzeta_preconditioner_positiveXi, d2dzeta2)
     scheme = 40
     call uniformDiffMatrices(Nzeta, 0, zetaMax, scheme, zeta_preconditioner, zetaWeights_preconditioner, ddzeta_preconditioner_negativeXi, d2dzeta2)
  case default
     stop "Error! Invalid zetaGridScheme_pc!"
  end select
  print *,"DDD"
  deallocate(zeta_preconditioner, zetaWeights_preconditioner, d2dzeta2)

  ! *****************************************************
  ! Set up xi grid
  ! *****************************************************

  allocate(xi(Nxi))
  allocate(xiWeights(Nxi))
  allocate(ddxi(Nxi,Nxi))
  allocate(d2dxi2(Nxi,Nxi))
  allocate(LorentzOperator(Nxi,Nxi))

  scheme = 12
  call uniformDiffMatrices(Nxi, -1.0d+0, 1.0d+0, scheme, xi, xiWeights, ddxi, d2dxi2)

  ! (1/2) (d/dxi) (1-xi^2) (d/dxi)
  LorentzOperator = d2dxi2/2
  do ixi=1,Nxi
     LorentzOperator(ixi,:) = (1-xi(ixi)*xi(ixi)) * d2dxi2(ixi,:)/2 &
          - xi(ixi)*ddxi(ixi,:)
  end do

  deallocate(d2dxi2)

  ! *****************************************************
  ! Compute B and its derivatives
  ! *****************************************************

  allocate(B(Ntheta,Nzeta))
  allocate(dBdtheta(Ntheta,Nzeta))
  allocate(dBdzeta(Ntheta,Nzeta))
  call computeB()

end subroutine createGrids


