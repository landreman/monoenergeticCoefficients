subroutine computeB()

  use variables, only: Ntheta, Nzeta, epsilon_t, epsilon_h, helicity_l, Nperiods, B, dBdtheta, dBdzeta, theta, zeta

  implicit none
  
!#include <finclude/petscsys.h>
#include <petsc/finclude/petscsysdef.h>

  PetscInt :: itheta, izeta

  do itheta = 1,Ntheta
     do izeta = 1,Nzeta

        B(itheta,izeta) = 1 + epsilon_t * cos(theta(itheta)) + epsilon_h * cos(helicity_l*theta(itheta) - Nperiods*zeta(izeta))

        dBdtheta(itheta,izeta) = -epsilon_t * sin(theta(itheta)) - helicity_l * epsilon_h * sin(helicity_l*theta(itheta) - Nperiods*zeta(izeta))

        dBdzeta(itheta,izeta) = Nperiods * epsilon_h * sin(helicity_l*theta(itheta) - Nperiods*zeta(izeta))
      
     end do
  end do

end subroutine computeB
