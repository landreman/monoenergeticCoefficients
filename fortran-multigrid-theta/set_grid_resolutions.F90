#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif

subroutine set_grid_resolutions()

  use petscsys
  use variables

  implicit none

  PetscScalar :: temp_float
  PetscInt :: temp_int, j
  PetscScalar, parameter :: small = 1.0d-12

  N_levels = 1
  Ntheta_min = min(Ntheta_min,Ntheta)
  Nzeta_min  = min(Nzeta_min ,Nzeta)
  Nxi_min    = min(Nxi_min   ,Nxi)
!  if (coarsen_theta) N_levels = max(N_levels,ceiling(log(Ntheta*one/Ntheta_min) / log(two) - small) +1)
!  if (coarsen_zeta ) N_levels = max(N_levels,ceiling(log(Nzeta *one/ Nzeta_min) / log(two) - small) +1)
!  if (coarsen_xi   ) N_levels = max(N_levels,ceiling(log((Nxi-one)/ (Nxi_min-one)) / log(two) - small) +1)
  if (coarsen_theta) N_levels = max(N_levels,nint(log(Ntheta*one/Ntheta_min) / log(two) - small) +1)
  if (coarsen_zeta ) N_levels = max(N_levels,nint(log(Nzeta *one/ Nzeta_min) / log(two) - small) +1)
  if (coarsen_xi   ) N_levels = max(N_levels,nint(log((Nxi-one)/ (Nxi_min-one)) / log(two) - small) +1)

  if (N_levels<1) stop "Error! N_levels<1"
  allocate(levels(N_levels))

  if (N_levels==1) then
     levels(1)%Ntheta = Ntheta
     levels(1)%Nzeta = Nzeta
     levels(1)%Nxi = Nxi
  else
     ! Handle theta:
     if (coarsen_theta) then
        do j=1,N_levels
           select case (coarsen_option)
           case (1)
              temp_float = max(one*Ntheta_min, Ntheta * (0.5 ** (j-1)))
           case (2)
              temp_float = exp(log(Ntheta*one) - (j-one)/(N_levels-one)*log(Ntheta*one/Ntheta_min))
           case default
              print *,"Invalid coarsen_option:",coarsen_option
              stop
           end select

           temp_int = nint(temp_float)
           if (mod(temp_int,2)==1) then
              levels(j)%Ntheta = temp_int
           else
              ! Check whether temp_int+1 or temp_int-1 is closer to the ideal value:
              if (abs(temp_int+1 - temp_float) > abs(temp_int-1 - temp_float)) then
                 levels(j)%Ntheta = temp_int-1
              else
                 levels(j)%Ntheta = temp_int+1
              end if
           end if
        end do
     else
        do j=1,N_levels
           levels(j)%Ntheta = Ntheta
        end do
     end if

     ! Handle zeta:
     if (coarsen_zeta) then
        do j=1,N_levels
           select case (coarsen_option)
           case (1)
              temp_float = max(one*Nzeta_min, Nzeta * (0.5 ** (j-1)))
           case (2)
              temp_float = exp(log(Nzeta*one) - (j-one)/(N_levels-one)*log(Nzeta*one/Nzeta_min))
           case default
              print *,"Invalid coarsen_option:",coarsen_option
              stop
           end select

           temp_int = nint(temp_float)
           if (mod(temp_int,2)==1) then
              levels(j)%Nzeta = temp_int
           else
              ! Check whether temp_int+1 or temp_int-1 is closer to the ideal value:
              if (abs(temp_int+1 - temp_float) > abs(temp_int-1 - temp_float)) then
                 levels(j)%Nzeta = temp_int-1
              else
                 levels(j)%Nzeta = temp_int+1
              end if
           end if
        end do
     else
        do j=1,N_levels
           levels(j)%Nzeta = Nzeta
        end do
     end if

     ! Handle xi:
     if (coarsen_xi) then
        select case (coarsen_option)
        case (1)
           do j=1,N_levels
              levels(j)%Nxi = max(Nxi_min, nint((Nxi-1) * (0.5 ** (j-1)) + 1))
           end do
        case (2)
           do j=1,N_levels
              levels(j)%Nxi = nint(exp(log(Nxi-one) - (j-one)/(N_levels-one)*log((Nxi-one)/(Nxi_min-one)))) + 1
           end do
        case default
           print *,"Invalid coarsen_option:",coarsen_option
           stop
        end select

     else
        do j=1,N_levels
           levels(j)%Nxi = Nxi
        end do
     end if

  end if

  if (masterProc) then
     print *,"----- Computed parameters for multigrid: ----"
     print *," Level Ntheta  Nzeta    Nxi"
     do j=1,N_levels
        print "(i7,i7,i7,i7)",j,levels(j)%Ntheta,levels(j)%Nzeta,levels(j)%Nxi
     end do
  end if

  if (Ntheta .ne. levels(1)%Ntheta) stop "Error! Ntheta should equal Ntheta_levels(1)"
  if (Nzeta  .ne. levels(1)%Nzeta ) stop "Error! Nzeta should equal Nzeta_levels(1)"
  if (Nxi    .ne. levels(1)%Nxi   ) stop "Error! Nxi should equal Nxi_levels(1)"

end subroutine set_grid_resolutions


