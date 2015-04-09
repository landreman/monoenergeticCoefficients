module geometry

  implicit none
  
!#include <finclude/petscsys.h>
#include <finclude/petscsysdef.h>
  
  contains

    function B(theta,zeta)

      use variables, only: epsilon_t, epsilon_h, helicity_l, Nperiods

      implicit none

      PetscReal, intent(in) :: theta, zeta
      PetscReal :: B

      B = 1 + epsilon_t * cos(theta) + epsilon_h * cos(helicity_l*theta - Nperiods*zeta)

    end function B

    ! ---------------------------------------------------------------------------------------

    function dBdtheta(theta,zeta)

      use variables, only: epsilon_t, epsilon_h, helicity_l, Nperiods

      implicit none

      PetscReal, intent(in) :: theta, zeta
      PetscReal :: dBdtheta
  
      dBdtheta = -epsilon_t * sin(theta) - helicity_l * epsilon_h * sin(helicity_l*theta - Nperiods*zeta)

    end function dBdtheta

    ! ---------------------------------------------------------------------------------------

    function dBdzeta(theta,zeta)

      use variables, only: epsilon_t, epsilon_h, helicity_l, Nperiods

      implicit none
      
      PetscReal, intent(in) :: theta, zeta
      PetscReal :: dBdzeta

      dBdzeta = Nperiods * epsilon_h * sin(helicity_l*theta - Nperiods*zeta)

    end function dBdzeta

    ! ---------------------------------------------------------------------------------------

    subroutine whereAmI(itheta,izeta,ixi,levelNtheta,levelNzeta,levelNxi,theta,zeta,xi)

      use variables, only: Nperiods, thetaCellCentered, zetaCellCentered, pi

      implicit none

      PetscInt, intent(in) :: itheta, izeta, ixi
      PetscInt, intent(in) :: levelNtheta, levelNzeta, levelNxi
      PetscReal, intent(out) :: theta, zeta, xi
      
      if (itheta<0) then
         print *,"Error! itheta cannot be <0. itheta=",itheta
         stop
      end if
      
      if (izeta<0) then
         print *,"Error! izeta cannot be <0. izeta=",izeta
         stop
      end if
      
      if (ixi<0) then
         print *,"Error! ixi cannot be <0. ixi=",ixi
         stop
      end if
      
      if (itheta>=levelNtheta) then
         print *,"Error! itheta cannot be >= levelNtheta. itheta=",itheta,", levelNtheta=",levelNtheta
         stop
      end if
      
      if (izeta>=levelNzeta) then
         print *,"Error! izeta cannot be >= levelNzeta. izeta=",izeta,", levelNzeta=",levelNzeta
         stop
      end if
      
      if (ixi>=levelNxi) then
         print *,"Error! ixi cannot be >= levelNxi. ixi=",ixi,", levelNxi=",levelNxi
         stop
      end if
      
      !globalIndex = (ixi-1)*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta
      
      ! The formulae below assume itheta, izeta, and ixi are 0-based, not 1-based!!
      
      ! theta and zeta are face-centered right now. Should they be cell-centered?
      if (thetaCellCentered) then
         theta = 2*pi*(itheta+0.5)/levelNtheta
      else
         theta = 2*pi*itheta/levelNtheta
      end if
      
      if (zetaCellCentered) then
         zeta = (2*pi/Nperiods)*(izeta+0.5)/levelNzeta
      else
         zeta = (2*pi/Nperiods)*izeta/levelNzeta
      end if
      
      ! xi is always cell-centered:
      xi = (2d+0)*(ixi+0.5d+0)/levelNxi - 1
      
    end subroutine whereAmI

  end module geometry
