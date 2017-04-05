#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

!subroutine computeB(Ntheta, theta, Nzeta, zeta, B, dBdtheta, dBdzeta)
subroutine computeB()

  use variables, only: epsilon_t, epsilon_h, helicity_l, Nperiods, levels, N_levels, VPrime

  implicit none
  
!!#include <finclude/petscsys.h>
!#include <finclude/petscsysdef.h>

!  integer, intent(in) :: Ntheta, Nzeta
!  PetscScalar, dimension(:), intent(in) :: theta, zeta
!  PetscScalar, dimension(:,:), intent(out) :: B, dBdtheta, dBdzeta
  PetscScalar, dimension(:), pointer :: theta, zeta
  PetscScalar, dimension(:,:), pointer :: B, dBdtheta, dBdzeta

  PetscInt :: itheta, izeta, level, Ntheta, Nzeta

  do level = 1,N_levels
     Ntheta = levels(level)%Ntheta
     Nzeta = levels(level)%Nzeta

     ! Allocate arrays for B and its derivatives     
     allocate(levels(level)%B(Ntheta,Nzeta))
     allocate(levels(level)%dBdtheta(Ntheta,Nzeta))
     allocate(levels(level)%dBdzeta(Ntheta,Nzeta))


     theta => levels(level)%theta
     zeta => levels(level)%zeta
     B => levels(level)%B
     dBdtheta => levels(level)%dBdtheta
     dBdzeta => levels(level)%dBdzeta

     do itheta = 1,Ntheta
        do izeta = 1,Nzeta

           B(itheta,izeta) = 1 + epsilon_t * cos(theta(itheta)) + epsilon_h * cos(helicity_l*theta(itheta) - Nperiods*zeta(izeta))
           
           dBdtheta(itheta,izeta) = -epsilon_t * sin(theta(itheta)) - helicity_l * epsilon_h * sin(helicity_l*theta(itheta) - Nperiods*zeta(izeta))

           dBdzeta(itheta,izeta) = Nperiods * epsilon_h * sin(helicity_l*theta(itheta) - Nperiods*zeta(izeta))
        end do
     end do

     if (level==1) then
        VPrime = 0
        do itheta = 1,Ntheta
           do izeta = 1,Nzeta
              VPrime = VPrime + levels(1)%thetaWeights(itheta) * levels(1)%zetaWeights(izeta) / (B(itheta,izeta) ** 2)
           end do
        end do
     end if

  end do

end subroutine computeB
