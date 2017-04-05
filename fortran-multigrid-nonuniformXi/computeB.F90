#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

!subroutine computeB(Ntheta, theta, Nzeta, zeta, B, dBdtheta, dBdzeta)
subroutine computeB()

  use variables, only: epsilon_t, epsilon_h, helicity_l, Nperiods, levels, N_levels, VPrime, FSAB2, geometry_option, G, I, iota, &
       masterProc, pi

  implicit none
  
!!#include <finclude/petscsys.h>
!#include <finclude/petscsysdef.h>

!  integer, intent(in) :: Ntheta, Nzeta
!  PetscScalar, dimension(:), intent(in) :: theta, zeta
!  PetscScalar, dimension(:,:), intent(out) :: B, dBdtheta, dBdzeta
  PetscScalar, dimension(:), pointer :: theta, zeta
  PetscScalar, dimension(:,:), pointer :: B, dBdtheta, dBdzeta
  PetscInt :: itheta, izeta, level, Ntheta, Nzeta

  character(200) :: equilibrium_file = '/Users/mattland/sfincs/equilibria/w7x-sc1-ecb2.bc'
  PetscScalar :: normradius_wish = 0.9, normradius
  integer :: fileUnit, didFileAccessWork, no_of_modes_old, no_of_modes_new, j, modeind, NHarmonics, numB0s
  PetscScalar :: psiAHat, a, B0_old, B0_new, pPrimeHat, pPrimeHat_old, pPrimeHat_new, B0OverBBar, min_Bmn_to_load=0
  PetscScalar :: GHat, IHat, iota_old, iota_new, G_old, G_new, I_old, I_new, normradius_old, normradius_new
  logical :: end_of_file, proceed
  integer, parameter :: max_no_of_modes = 10000
  integer, dimension(max_no_of_modes) :: modesm_old, modesm_new, modesn_old, modesn_new
  PetscScalar, dimension(max_no_of_modes) :: modesb_old, modesb_new
  integer, dimension(:), allocatable :: BHarmonics_l, BHarmonics_n
  PetscScalar, dimension(:), allocatable :: BHarmonics_amplitudes
  character(len=200) :: lineOfFile
  integer, dimension(4) :: headerIntegers
  PetscScalar, dimension(3) :: headerReals
  PetscScalar, dimension(6) :: surfHeader
  PetscScalar, dimension(4) :: dataNumbers
  PetscScalar, dimension(8) :: data8Numbers
  integer, dimension(2) :: dataIntegers

  select case (geometry_option)
  case (1)
  case (2)

       ! Read Boozer coordinate file in .bc format used at IPP Greifswald

       fileUnit = 11
       open(unit=fileUnit, file=equilibrium_file, action="read", status="old", iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Unable to open magnetic equilibrium file ",equilibrium_file
          stop
       else
          do
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             if (lineOfFile(1:2) /= "CC") exit
          end do

          ! Read header line:
          read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) headerIntegers, headerReals
          if (didFileAccessWork /= 0) then
             print *,"Unable to read header from the magnetic equilibrium file ",equilibrium_file
             stop
          end if

          NPeriods = headerIntegers(4)
          psiAHat  = headerReals(1)/2/pi !Convert the flux from Tm^2 to Tm^2/rad
          a        = headerReals(2)      !minor radius in meters

          end_of_file = .false.

          normradius_old = 0
          no_of_modes_old = 0
          modesm_old = 0
          modesn_old = 0
          modesb_old = 0
          iota_old = 0
          G_old = 0
          I_old = 0
          B0_old = 0
          pPrimeHat_old = 0

          normradius_new = 0
          no_of_modes_new = 0
          modesm_new = 0
          modesn_new = 0
          modesb_new = 0
          iota_new = 0
          G_new = 0
          I_new = 0
          B0_new = 0
          pPrimeHat_new = 0

          ! Skip a line
          read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile

          do 
             if ((normradius_new .ge. normradius_wish) .or. end_of_file) exit

             normradius_old = normradius_new
             no_of_modes_old = no_of_modes_new
             modesm_old = modesm_new
             modesn_old = modesn_new
             modesb_old = modesb_new
             iota_old = iota_new
             G_old = G_new
             I_old = I_new
             B0_old = B0_new
             pPrimeHat_old = pPrimeHat_new
             numB0s = 0

             ! Skip a line:
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             ! Read the header for the magnetic surface:
             read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) surfHeader

             normradius_new = sqrt(surfHeader(1))       ! r/a = sqrt(psi/psi_a)
             iota_new = surfHeader(2)
             ! Note that G and I has a minus sign in the following two lines
             ! because Ampere's law comes with a minus sign in the left-handed
             ! (r,pol,tor) system.
             G_new = -surfHeader(3)*NPeriods/2/pi*(4*pi*1d-7) !Tesla*meter
             I_new = -surfHeader(4)/2/pi*(4*pi*1d-7)          !Tesla*meter
             pPrimeHat_new = surfheader(5)*(4*pi*1e-7)       ! p=pHat \bar{B}^2 / \mu_0

             ! Skip units line:
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             proceed = .true.
             modeind = 0
             do
                if (.not. proceed) exit

                read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
                if (didFileAccessWork /= 0) then
                   proceed = .false.
                   end_of_file = .true.
                else if (index(lineOfFile,"s") > 0) then
                   ! Next flux surface has been reached
                   proceed = .false.
                else
                   read(unit=lineOfFile, fmt=*) dataIntegers, dataNumbers
                   if (dataIntegers(1) == 0 .and. dataIntegers(2) == 0) then
                      B0_new = dataNumbers(4)
                      numB0s = numB0s + 1
                   else if (abs(dataNumbers(4)) > min_Bmn_to_load) then
                      modeind = modeind + 1
                      if (modeind > max_no_of_modes) then
                         print *,"The value of max_no_of_modes in geometry.F90 was insufficient."
                         print *,"Either increase this value and recompile, or else increase min_Bmn_to_load."
                         stop
                      end if
                      modesm_new(modeind) = dataIntegers(1)
                      modesn_new(modeind) = dataIntegers(2)
                      modesb_new(modeind) = dataNumbers(4)
                   end if
                end if
             end do
             if (numB0s == 0) then
                print *,"Error: no (0,0) mode found in magnetic equilibrium file ",equilibrium_file
             else if (numB0s > 1) then
                print *,"Error: more than 1 (0,0) mode found in magnetic equilibrium file ",equilibrium_file
             end if
             no_of_modes_new = modeind
          end do

       end if

       close(unit = fileUnit)
       if (masterProc) then
          print *,"Successfully read magnetic equilibrium from file ",trim(equilibrium_file)
       end if

       if (abs(normradius_old - normradius_wish) < abs(normradius_new - normradius_wish)) then
          iota = iota_old
          GHat = G_old
          IHat = I_old
          normradius = normradius_old
          B0OverBBar = B0_old
          NHarmonics = no_of_modes_old
          pPrimeHat=pPrimeHat_old
          allocate(BHarmonics_l(NHarmonics))
          allocate(BHarmonics_n(NHarmonics))
          allocate(BHarmonics_amplitudes(NHarmonics))
          BHarmonics_l = modesm_old(1:NHarmonics)
          BHarmonics_n = modesn_old(1:NHarmonics)
          BHarmonics_amplitudes = modesb_old(1:NHarmonics)
       else
          iota = iota_new
          GHat = G_new
          IHat = I_new
          normradius = normradius_new
          B0OverBBar = B0_new
          NHarmonics = no_of_modes_new
          pPrimeHat=pPrimeHat_new
          allocate(BHarmonics_l(NHarmonics))
          allocate(BHarmonics_n(NHarmonics))
          allocate(BHarmonics_amplitudes(NHarmonics))
          BHarmonics_l = modesm_new(1:NHarmonics)
          BHarmonics_n = modesn_new(1:NHarmonics)
          BHarmonics_amplitudes = modesb_new(1:NHarmonics)
       end if

       BHarmonics_amplitudes = BHarmonics_amplitudes / B0OverBBar

       if (GHat*psiAHat<0) then
          !print *,"This is a stellarator symmetric file from Joachim Geiger. It will now be turned 180 degrees around a horizontal axis <=> flip the sign of G and I, so that it matches the sign of its total toroidal flux."
          GHat=-GHat
          IHat=-IHat
       end if
       
       !Switch from a left-handed to right-handed (radial,poloidal,toroidal) system
       psiAHat=psiAHat*(-1)           !toroidal direction sign switch
       GHat = GHat*(-1)               !toroidal direction sign switch
       iota = iota*(-1)               !toroidal direction sign switch
       BHarmonics_n=BHarmonics_n*(-1) !toroidal direction sign switch

       if (masterProc) then
          print *,"This computation is for the flux surface with minor radius ",normradius*a, &
               " meters, equivalent to r/a = ",normradius
       end if

    case (12)
       ! Read Boozer coordinate file in a generalisation of the .bc format used at IPP Greifswald for non-stellarator symmetric equilibria 

       fileUnit = 11
       open(unit=fileUnit, file=equilibrium_file, action="read", status="old", iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Unable to open magnetic equilibrium file ",equilibrium_file
          stop
       else
          do
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             if (lineOfFile(1:2) /= "CC") exit
          end do

          ! Read header line:
          read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) headerIntegers, headerReals
          if (didFileAccessWork /= 0) then
             print *,"Unable to read header from the magnetic equilibrium file ",equilibrium_file
             stop
          end if

          NPeriods = headerIntegers(4)
          psiAHat  = headerReals(1)/2/pi !Convert the flux from Tm^2 to Tm^2/rad
          a        = headerReals(2)      !minor radius in meters

          end_of_file = .false.

          normradius_old = 0
          no_of_modes_old = 0
          modesm_old = 0
          modesn_old = 0
          modesb_old = 0
          iota_old = 0
          G_old = 0
          I_old = 0
          B0_old = 0
          pPrimeHat_old = 0

          normradius_new = 0
          no_of_modes_new = 0
          modesm_new = 0
          modesn_new = 0
          modesb_new = 0
          iota_new = 0
          G_new = 0
          I_new = 0
          B0_new = 0
          pPrimeHat_new = 0

          ! Skip a line
          read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile

          do 
             if ((normradius_new .ge. normradius_wish) .or. end_of_file) exit

             normradius_old = normradius_new
             no_of_modes_old = no_of_modes_new
             modesm_old = modesm_new
             modesn_old = modesn_new
             modesb_old = modesb_new
             iota_old = iota_new
             G_old = G_new
             I_old = I_new
             B0_old = B0_new
             pPrimeHat_old = pPrimeHat_new
             numB0s = 0

             ! Skip a line:
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             ! Read the header for the magnetic surface:
             read(unit=fileUnit, iostat=didFileAccessWork, fmt=*) surfHeader

             normradius_new = sqrt(surfHeader(1))       ! r/a = sqrt(psi/psi_a)
             iota_new = surfHeader(2)
             ! Note that G and I has a minus sign in the following two lines
             ! because Ampere's law comes with a minus sign in the left-handed
             ! (r,pol,tor) system.
             G_new = -surfHeader(3)*NPeriods/2/pi*(4*pi*1d-7) !Tesla*meter
             I_new = -surfHeader(4)/2/pi*(4*pi*1d-7)          !Tesla*meter
             pPrimeHat_new = surfheader(5)*(4*pi*1e-7)       ! p=pHat \bar{B}^2 / \mu_0

             ! Skip units line:
             read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
             proceed = .true.
             modeind = 0
             do
                if (.not. proceed) exit

                read(unit=fileUnit, fmt="(a)", iostat=didFileAccessWork) lineOfFile
                if (didFileAccessWork /= 0) then
                   proceed = .false.
                   end_of_file = .true.
                else if (index(lineOfFile,"s") > 0) then
                   ! Next flux surface has been reached
                   proceed = .false.
                else
                   read(unit=lineOfFile, fmt=*) dataIntegers, data8Numbers
                   if (dataIntegers(1) == 0 .and. dataIntegers(2) == 0) then
                      B0_new = data8Numbers(7)
                      numB0s = numB0s + 1
                   else if (abs(data8Numbers(7)) > min_Bmn_to_load) then
                      if (modeind + 2 > max_no_of_modes) then
                         print *,"The value of max_no_of_modes in geometry.F90 was insufficient."
                         print *,"Either increase this value and recompile, or else increase min_Bmn_to_load."
                         stop
                      end if
                      modeind = modeind + 1
                      modesm_new(modeind) = dataIntegers(1)
                      modesn_new(modeind) = dataIntegers(2)
                      modesb_new(modeind) = data8Numbers(7) !Cosinus component
                      modeind = modeind + 1
                      modesm_new(modeind) = dataIntegers(1)
                      modesn_new(modeind) = dataIntegers(2)
                      modesb_new(modeind) = data8Numbers(8) !Sinus component
                   end if
                end if
             end do
             if (numB0s == 0) then
                print *,"Error: no (0,0) mode found in magnetic equilibrium file ",equilibrium_file
             else if (numB0s > 1) then
                print *,"Error: more than 1 (0,0) mode found in magnetic equilibrium file ",equilibrium_file
             end if
             no_of_modes_new = modeind
          end do

       end if

       close(unit = fileUnit)
       if (masterProc) then
          print *," Successfully read magnetic equilibrium from file ",trim(equilibrium_file)
       end if

       if (abs(normradius_old - normradius_wish) < abs(normradius_new - normradius_wish)) then
          iota = iota_old
          GHat = G_old
          IHat = I_old
          normradius = normradius_old
          B0OverBBar = B0_old
          NHarmonics = no_of_modes_old
          allocate(BHarmonics_l(NHarmonics))
          allocate(BHarmonics_n(NHarmonics))
          allocate(BHarmonics_amplitudes(NHarmonics))
          BHarmonics_l = modesm_old(1:NHarmonics)
          BHarmonics_n = modesn_old(1:NHarmonics)
          BHarmonics_amplitudes = modesb_old(1:NHarmonics)
       else
          iota = iota_new
          GHat = G_new
          IHat = I_new
          normradius = normradius_new
          B0OverBBar = B0_new
          NHarmonics = no_of_modes_new
          allocate(BHarmonics_l(NHarmonics))
          allocate(BHarmonics_n(NHarmonics))
          allocate(BHarmonics_amplitudes(NHarmonics))
          BHarmonics_l = modesm_new(1:NHarmonics)
          BHarmonics_n = modesn_new(1:NHarmonics)
          BHarmonics_amplitudes = modesb_new(1:NHarmonics)
       end if

       BHarmonics_amplitudes = BHarmonics_amplitudes / B0OverBBar
       
       !Switch from a left-handed to right-handed (radial,poloidal,toroidal) system
       psiAHat=psiAHat*(-1)           !toroidal direction sign switch
       GHat = GHat*(-1)               !toroidal direction sign switch
       iota = iota*(-1)               !toroidal direction sign switch
       BHarmonics_n=BHarmonics_n*(-1) !toroidal direction sign switch

      if (masterProc) then
          print *," This computation is for the flux surface with minor radius ",normradius*a, &
               " meters, equivalent to r/a = ",normradius
      end if

      G = GHat
      I = IHat

  case default
     print *,"Invalid geometry_option:",geometry_option
     stop
  end select

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

     select case (geometry_option)
     case (1)
        do itheta = 1,Ntheta
           do izeta = 1,Nzeta
              
              B(itheta,izeta) = 1 + epsilon_t * cos(theta(itheta)) + epsilon_h * cos(helicity_l*theta(itheta) - Nperiods*zeta(izeta))
              
              dBdtheta(itheta,izeta) = -epsilon_t * sin(theta(itheta)) - helicity_l * epsilon_h * sin(helicity_l*theta(itheta) - Nperiods*zeta(izeta))
              
              dBdzeta(itheta,izeta) = Nperiods * epsilon_h * sin(helicity_l*theta(itheta) - Nperiods*zeta(izeta))
           end do
        end do
     case (2)
        ! Initialize arrays:
        B = B0OverBBar ! This includes the (0,0) component.
        dBdtheta = 0
        dBdzeta = 0
        
        do j = 1, NHarmonics
           do itheta = 1,Ntheta
              B(itheta,:) = B(itheta,:) + B0OverBBar * BHarmonics_amplitudes(j) * &
                   cos(BHarmonics_l(j) * theta(itheta) - NPeriods * BHarmonics_n(j) * zeta)
              
              dBdtheta(itheta,:) = dBdtheta(itheta,:) - B0OverBBar * BHarmonics_amplitudes(j) * BHarmonics_l(j) * &
                   sin(BHarmonics_l(j) * theta(itheta) - NPeriods * BHarmonics_n(j) * zeta)
              
              dBdzeta(itheta,:) = dBdzeta(itheta,:) + B0OverBBar * BHarmonics_amplitudes(j) * Nperiods * BHarmonics_n(j) * &
                   sin(BHarmonics_l(j) * theta(itheta) - NPeriods * BHarmonics_n(j) * zeta)
              
           end do
        end do

     end select

     if (level==1) then
        VPrime = 0
        do itheta = 1,Ntheta
           do izeta = 1,Nzeta
              VPrime = VPrime + levels(1)%thetaWeights(itheta) * levels(1)%zetaWeights(izeta) / (B(itheta,izeta) ** 2)
           end do
        end do
        FSAB2 = sum(levels(1)%thetaWeights)*sum(levels(1)%zetaWeights)/VPrime

        if (masterProc) then
           print *,"Final iota:",iota
           print *,"Final G:",G
           print *,"Final I:",I
           print *,"Final Nperiods:",Nperiods
           print *,"zeta:",zeta
           print *,"Here comes B:"
           do itheta=1,Ntheta
              print "(*(f5.2))",B(itheta,:)
           end do
        end if
     end if

  end do

end subroutine computeB
