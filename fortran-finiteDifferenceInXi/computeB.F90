#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

subroutine computeB()

  use variables, only: Ntheta, Nzeta, epsilon_t, epsilon_h, helicity_l, Nperiods, B, dBdtheta, dBdzeta, theta, zeta, &
       VPrime, FSAB2, thetaWeights, zetaWeights

  implicit none
  
!!#include <finclude/petscsys.h>
!#include <finclude/petscsysdef.h>

  PetscInt :: itheta, izeta

  do itheta = 1,Ntheta
     do izeta = 1,Nzeta

        B(itheta,izeta) = 1 + epsilon_t * cos(theta(itheta)) + epsilon_h * cos(helicity_l*theta(itheta) - Nperiods*zeta(izeta))

        dBdtheta(itheta,izeta) = -epsilon_t * sin(theta(itheta)) - helicity_l * epsilon_h * sin(helicity_l*theta(itheta) - Nperiods*zeta(izeta))

        dBdzeta(itheta,izeta) = Nperiods * epsilon_h * sin(helicity_l*theta(itheta) - Nperiods*zeta(izeta))
      
     end do
  end do

  VPrime = 0
  do itheta = 1,Ntheta
     do izeta = 1,Nzeta
        VPrime = VPrime + thetaWeights(itheta)*zetaWeights(izeta)/(B(itheta,izeta) ** 2)
     end do
  end do
  FSAB2 = sum(thetaWeights)*sum(zetaWeights)/VPrime
     
end subroutine computeB
