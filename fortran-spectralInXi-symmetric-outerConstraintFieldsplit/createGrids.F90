#include <petsc/finclude/petscsysdef.h>

subroutine createGrids()

  use petscsys
  use variables

  implicit none

  PetscScalar, dimension(:,:), allocatable :: d2dtheta2
  PetscScalar, dimension(:,:), allocatable :: d2dzeta2
  PetscScalar :: zetaMax
  integer :: L, j

  if (masterProc) print *,"Entering createGrids."

  matrixSize = Ntheta*Nzeta*Nxi + Ntheta*Nzeta + 2
  if (masterProc) print *,"matrixSize=",matrixSize

  allocate(L_scaling(Nxi))
  select case (L_scaling_option)
  case (1)
     L_scaling = 1
  case (2)
     L_scaling = [( sqrt((2*L+1.0d+0)/2), L=0, Nxi-1 )]
  case default
     print *,"Error! Invalid for L_scaling_option:",L_scaling_option
     stop
  end select
  if (masterProc) print *,"L_scaling:",L_scaling

  ! *****************************************************
  ! Set up theta grid
  ! *****************************************************

  allocate(theta(Ntheta))
  allocate(thetaWeights(Ntheta))
  allocate(ddtheta(Ntheta,Ntheta))
  allocate(d2dtheta2(Ntheta,Ntheta))

  call uniformDiffMatrices(Ntheta, 0.0d+0, 2*pi, thetaGridScheme, theta, thetaWeights, &
       ddtheta, d2dtheta2)

  deallocate(d2dtheta2)

  thetaWeight = thetaWeights(1)

  ! *****************************************************
  ! Set up zeta grid
  ! *****************************************************

  allocate(zeta(Nzeta))
  allocate(zetaWeights(Nzeta))
  allocate(ddzeta(Nzeta,Nzeta))
  allocate(d2dzeta2(Nzeta,Nzeta))

  zetaMax = 2*pi/Nperiods

  call uniformDiffMatrices(Nzeta, 0.0d+0, zetaMax, zetaGridScheme, zeta, zetaWeights, &
       ddzeta, d2dzeta2)

  deallocate(d2dzeta2)

  zetaWeight = zetaWeights(1)

  ! *****************************************************
  ! Compute B and its derivatives
  ! *****************************************************

  allocate(B(Ntheta,Nzeta))
  allocate(dBdtheta(Ntheta,Nzeta))
  allocate(dBdzeta(Ntheta,Nzeta))
  call computeB()

  allocate(sqrt_g(Ntheta,Nzeta))
  sqrt_g = (G+iota*I) / (B*B)
  VPrime = thetaWeight * zetaWeight * sum(sqrt_g)
  FSAB2 = thetaWeight * zetaWeight * sum(sqrt_g*B*B) / VPrime

  if (masterProc) then
     print *,"Here comes B:"
     do j=1,Ntheta
        print "(*(f7.2))",B(j,:)
     end do
     print *,"Here comes sqrt_g:"
     do j=1,Ntheta
        print "(*(f7.2))",sqrt_g(j,:)
     end do
     print *,"Here comes ddtheta:"
     do j=1,Ntheta
        print "(*(f7.2))",ddtheta(j,:)
     end do
     print *,"Here comes ddzeta:"
     do j=1,Nzeta
        print "(*(f7.2))",ddzeta(j,:)
     end do
  end if

  if (masterProc) print *,"Leaving createGrids."

end subroutine createGrids


