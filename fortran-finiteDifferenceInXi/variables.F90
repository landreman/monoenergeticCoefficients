
module variables

  implicit none

!#include <finclude/petscsys.h>
#include <finclude/petscsysdef.h>

  PetscInt :: Ntheta, Nzeta, Nxi, Nperiods, helicity_l
  PetscReal :: nu, epsilon_t, epsilon_h, iota, G, I
  PetscOffset :: offset
  PetscInt :: myRank, numProcs
  logical :: masterProc
  logical :: thetaCellCentered, zetaCellCentered
  logical :: upwindTheta, upwindZeta, upwindInPCOnly

  PetscScalar, parameter :: pi = 3.14159265358979d+0
  PetscScalar, parameter :: sqrtpi = 1.77245385090552d+0

end module variables
