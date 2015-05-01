
module variables

  use stel_kinds

  implicit none

!#include <finclude/petscsysdef.h>

  integer :: Ntheta, Nzeta, Nxi, Nperiods, helicity_l, matrixSize
  real(rprec) :: nu, epsilon_t, epsilon_h, iota, G, I
  integer :: myRank, numProcs
  logical :: masterProc

  real(rprec), dimension(:), allocatable :: theta, zeta
  real(rprec), dimension(:), allocatable :: thetaWeights, zetaWeights
  real(rprec), dimension(:,:), allocatable :: ddtheta
  real(rprec), dimension(:,:), allocatable :: ddzeta
  real(rprec), dimension(:,:), allocatable :: B, dBdtheta, dBdzeta
  real(rprec) :: diagonalShift

  real(rprec), parameter :: pi = 3.14159265358979d+0
  real(rprec), parameter :: sqrtpi = 1.77245385090552d+0

  integer :: thetaGridScheme, zetaGridScheme

  REAL(rprec), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: lblk, dblk, ublk
  REAL(rprec), ALLOCATABLE, TARGET  :: brhs(:,:)

  integer :: mblock

end module variables
