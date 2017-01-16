#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

subroutine computeB(Ntheta, theta, Nzeta, zeta, B, dBdtheta, dBdzeta)

  use variables, only: epsilon_t, epsilon_h, helicity_l, Nperiods

  implicit none
  
!!#include <finclude/petscsys.h>
!#include <finclude/petscsysdef.h>

  integer, intent(in) :: Ntheta, Nzeta
  PetscScalar, dimension(:), intent(in) :: theta, zeta
  PetscScalar, dimension(:,:), intent(out) :: B, dBdtheta, dBdzeta

  PetscInt :: itheta, izeta

  do itheta = 1,Ntheta
     do izeta = 1,Nzeta

        B(itheta,izeta) = 1 + epsilon_t * cos(theta(itheta)) + epsilon_h * cos(helicity_l*theta(itheta) - Nperiods*zeta(izeta))

        dBdtheta(itheta,izeta) = -epsilon_t * sin(theta(itheta)) - helicity_l * epsilon_h * sin(helicity_l*theta(itheta) - Nperiods*zeta(izeta))

        dBdzeta(itheta,izeta) = Nperiods * epsilon_h * sin(helicity_l*theta(itheta) - Nperiods*zeta(izeta))
      
     end do
  end do

end subroutine computeB
