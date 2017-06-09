! Main program

#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif

program mmc

  use petscksp

  use variables

  implicit none

!#include <finclude/petscsys.h>
!#include <finclude/petscvec.h>
!#include <finclude/petscmat.h>
!#include <finclude/petscksp.h>

  PetscErrorCode ierr
  PetscBool :: wasSet
  Vec :: rhs, solution
  PetscInt :: userContext ! Not used
  PetscInt :: VecLocalSize
  integer :: unit, reason
  
#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR > 6))
  PetscViewerAndFormat vf
#endif

  !external populateMatrix, populateRHS
  external apply_preconditioner

  call PETSCInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_COMM_SIZE(PETSC_COMM_WORLD, numProcs, ierr)
  call MPI_COMM_RANK(PETSC_COMM_WORLD, myRank, ierr)

  print *,"I am proc",myRank," of ",numProcs
  masterProc = (myRank==0)

  ! Set defaults:
  nu = 0.1d+0
  E = 0
  epsilon_t = -0.07053d+0
  epsilon_h = 0.05067d+0
  iota = 0.4542d+0
  G = 3.7481d+0
  I = 0d+0
  Nperiods = 10
  helicity_l = 2
  Ntheta = 13
  Nzeta = 15
  Nxi = 16

  theta_derivative_option = 8
  preconditioner_theta_derivative_option = 4
  zeta_derivative_option = 8
  preconditioner_zeta_derivative_option = 4

  constraint_option = 1

  coarsen_theta = PETSC_FALSE
  coarsen_zeta = PETSC_TRUE
  coarsen_xi = PETSC_FALSE
  Ntheta_min = 7
  Nzeta_min = 7
  Nxi_min = 9

  smoothing_option = 0 ! Default was 1 before 20170608.
  restriction_option = 1
  coarsen_option = 1
  theta_diffusion = 0
  zeta_diffusion = 0
  defect_option = 1
  shift_L0_in_smoother = .false.

  f_scaling_option = 1
  L_scaling_option = 3
  omega = 1
  upwinding_scale_factor = 1
  preconditioner_upwinding_scale_factor = 1
  preconditioning_option = 1

#if (PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR > 6))
#define new_argument PETSC_NULL_OBJECT,
#else
#define new_argument
#endif

  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Ntheta', Ntheta, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Nzeta', Nzeta, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Nxi', Nxi, wasSet, ierr)
  call PetscOptionsGetReal(new_argument PETSC_NULL_CHARACTER, '-nu', nu, wasSet, ierr)
  call PetscOptionsGetReal(new_argument PETSC_NULL_CHARACTER, '-E', E, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-theta_derivative_option', theta_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-preconditioner_theta_derivative_option', preconditioner_theta_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-zeta_derivative_option', zeta_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-preconditioner_zeta_derivative_option', preconditioner_zeta_derivative_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-constraint_option', constraint_option, wasSet, ierr)
  call PetscOptionsGetBool(new_argument PETSC_NULL_CHARACTER, '-coarsen_theta', coarsen_theta, wasSet, ierr)
  call PetscOptionsGetBool(new_argument PETSC_NULL_CHARACTER, '-coarsen_zeta', coarsen_zeta, wasSet, ierr)
  call PetscOptionsGetBool(new_argument PETSC_NULL_CHARACTER, '-coarsen_xi', coarsen_xi, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Ntheta_min', Ntheta_min, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Nzeta_min', Nzeta_min, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-Nxi_min', Nxi_min, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-smoothing_option', smoothing_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-restriction_option', restriction_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-coarsen_option', coarsen_option, wasSet, ierr)
  call PetscOptionsGetReal(new_argument PETSC_NULL_CHARACTER, '-zeta_diffusion', zeta_diffusion, wasSet, ierr)
  call PetscOptionsGetReal(new_argument PETSC_NULL_CHARACTER, '-theta_diffusion', theta_diffusion, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-defect_option', defect_option, wasSet, ierr)
  call PetscOptionsGetBool(new_argument PETSC_NULL_CHARACTER, '-shift_L0_in_smoother', shift_L0_in_smoother, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-f_scaling_option', f_scaling_option, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-L_scaling_option', L_scaling_option, wasSet, ierr)
  call PetscOptionsGetReal(new_argument PETSC_NULL_CHARACTER, '-omega', omega, wasSet, ierr)
  call PetscOptionsGetReal(new_argument PETSC_NULL_CHARACTER, '-upwinding_scale_factor', upwinding_scale_factor, wasSet, ierr)
  call PetscOptionsGetReal(new_argument PETSC_NULL_CHARACTER, '-preconditioner_upwinding_scale_factor', preconditioner_upwinding_scale_factor, wasSet, ierr)
  call PetscOptionsGetInt(new_argument PETSC_NULL_CHARACTER, '-preconditioning_option', preconditioning_option, wasSet, ierr)

  ! Make sure Ntheta and Nzeta are odd:
  if (mod(Ntheta, 2) == 0) then
     Ntheta = Ntheta + 1
  end if
  if (mod(Nzeta, 2) == 0) then
     Nzeta = Nzeta + 1
  end if

  if (masterProc) then
     print *,"Ntheta = ",Ntheta
     print *,"Nzeta = ",Nzeta
     print *,"Nxi = ",Nxi
     print *,"nu = ",nu
     print *,"E = ",E
     print *,"theta_derivative_option = ",theta_derivative_option
     print *,"preconditioner_theta_derivative_option = ",preconditioner_theta_derivative_option
     print *,"zeta_derivative_option = ",zeta_derivative_option
     print *,"preconditioner_zeta_derivative_option = ",preconditioner_zeta_derivative_option
     print *,"constraint_option = ",constraint_option
     print *,"coarsen_theta = ",coarsen_theta
     print *,"coarsen_zeta = ",coarsen_zeta
     print *,"coarsen_xi = ",coarsen_xi
     print *,"Ntheta_min = ",Ntheta_min
     print *,"Nzeta_min = ",Nzeta_min
     print *,"Nxi_min = ",Nxi_min
     print *,"smoothing_option = ",smoothing_option
     print *,"restriction_option = ",restriction_option
     print *,"theta_diffusion = ",theta_diffusion
     print *,"zeta_diffusion = ",zeta_diffusion
     print *,"defect_option = ",defect_option
     print *,"shift_L0_in_smoother = ",shift_L0_in_smoother
     print *,"f_scaling_option = ",f_scaling_option
     print *,"L_scaling_option = ",L_scaling_option
     print *,"omega = ",omega
     print *,"upwinding_scale_factor = ",upwinding_scale_factor
     print *,"preconditioner_upwinding_scale_factor = ",preconditioner_upwinding_scale_factor
     print *,"preconditioning_option = ",preconditioning_option
  end if

  if (omega<0) stop "omega must be >= 0."
  if (omega>1) stop "omega must be <= 1."

  if (constraint_option<0 .or. constraint_option>2) stop "Invalid constraint_option"

  call KSPCreate(PETSC_COMM_WORLD,outer_ksp,ierr)
  call KSPAppendOptionsPrefix(outer_ksp, 'outer_', ierr)
  call KSPSetType(outer_ksp, KSPFGMRES, ierr)
  call KSPGetPC(outer_ksp,outer_preconditioner,ierr)
  call PCSetType(outer_preconditioner, PCSHELL, ierr)
  call PCShellSetApply(outer_preconditioner, apply_preconditioner, ierr)

  call KSPCreate(PETSC_COMM_WORLD, inner_ksp, ierr)
  call KSPAppendOptionsPrefix(inner_ksp, 'inner_', ierr)
  call KSPSetType(inner_ksp, KSPPREONLY, ierr)
  call KSPGetPC(inner_ksp, inner_preconditioner, ierr)
  call PCSetType(inner_preconditioner, PCMG,ierr)

  call setup_multigrid()

  call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrixSize, rhs, ierr)
  call VecDuplicate(rhs, solution, ierr)

  call populateRHS(rhs)

!!$  ! Create the object for the multigrid preconditioner:
!!$  call KSPGetPC(ksp,preconditioner_context,ierr)
!!$  call PCSetType(preconditioner_context,PCSHELL,ierr)
!!$  call PCShellSetApply(preconditioner_context,apply_multigrid_cycle,ierr)


!  call KSPSetComputeRHS(ksp,populateRHS,userContext,ierr)
!  call KSPSetComputeOperators(ksp,populateMatrix,userContext,ierr)
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 7))
  call KSPMonitorSet(outer_ksp, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
#else
  call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, vf, ierr) 
  !call KSPMonitorSet(ksp, KSPMonitorDefault, vf, PetscViewerAndFormatDestroy, ierr)
  call KSPMonitorSet(outer_ksp, KSPMonitorTrueResidualNorm, vf, PetscViewerAndFormatDestroy, ierr)
#endif
  call KSPSetFromOptions(outer_ksp,ierr)
  call KSPSetFromOptions(inner_ksp,ierr)
  
  !call VecView(rhs, PETSC_VIEWER_STDOUT_WORLD, ierr)
  if (masterProc) then
     print *,"Beginning solve..."
  end if
  call KSPSolve(outer_ksp, rhs, solution, ierr)
  if (masterProc) then
     print *,"Done!"
  end if

  call KSPGetConvergedReason(outer_ksp, reason, ierr)
  if (reason<=0) then
     if (masterProc) print *,"Error! Outer KSP did not converge. KSPConvergedReason=",reason
  end if

  call diagnostics(solution)

  call KSPView(outer_KSP, PETSC_VIEWER_STDOUT_WORLD, ierr)

  call KSPDestroy(outer_ksp,ierr)

  call PETScFinalize(ierr)

  ! Write an ascii output file
  unit = 9
  open(unit, file='mmc_out')
  write (unit,*) Ntheta
  write (unit,*) Nzeta
  write (unit,*) Nxi
  write (unit,*) nu
  write (unit,*) E
  write (unit,*) theta_derivative_option
  write (unit,*) preconditioner_theta_derivative_option
  write (unit,*) zeta_derivative_option
  write (unit,*) preconditioner_zeta_derivative_option
  write (unit,*) constraint_option
  write (unit,*) logical_to_int(coarsen_theta)
  write (unit,*) logical_to_int(coarsen_zeta)
  write (unit,*) logical_to_int(coarsen_xi)
  write (unit,*) Ntheta_min
  write (unit,*) Nzeta_min
  write (unit,*) Nxi_min
  write (unit,*) smoothing_option
  write (unit,*) restriction_option
  write (unit,*) defect_option
  write (unit,*) flux
  write (unit,*) flow
  close(unit)


contains

  function logical_to_int(input_logical)
    logical, intent(in) :: input_logical
    integer :: logical_to_int
    if (input_logical) then
       logical_to_int = 1
    else
       logical_to_int = 0
    end if
  end function logical_to_int

end program


