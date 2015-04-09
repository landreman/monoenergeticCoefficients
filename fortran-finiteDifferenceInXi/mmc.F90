! Main program

#include <finclude/petsckspdef.h>
#include <finclude/petscdmdadef.h>

program mmc

  use petscksp
  use petscdmda

  use variables

  implicit none

!#include <finclude/petscsys.h>
!#include <finclude/petscvec.h>
!#include <finclude/petscmat.h>
!#include <finclude/petscksp.h>
!#include <finclude/petscdm.h>
!#include <finclude/petscdmda.h>

  PetscErrorCode ierr
  PetscBool :: wasSet
  KSP :: ksp
  Vec :: solution
  DM :: dmda
  PetscInt :: userContext ! Not used
  PetscInt :: dof, stencilWidth

  external populateMatrix, populateRHS

  call PETSCInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_COMM_SIZE(PETSC_COMM_WORLD, numProcs, ierr)
  call MPI_COMM_RANK(PETSC_COMM_WORLD, myRank, ierr)

  print *,"I am proc",myRank," of ",numProcs
  masterProc = (myRank==0)

  ! Set defaults:
  nu = 0.1d+0
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
  thetaCellCentered = .true.
  zetaCellCentered = .true.
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-Ntheta', Ntheta, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-Nzeta', Nzeta, wasSet, ierr)
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-Nxi', Nxi, wasSet, ierr)
  call PetscOptionsGetReal(PETSC_NULL_CHARACTER, '-nu', nu, wasSet, ierr)
  call PetscOptionsGetBool(PETSC_NULL_CHARACTER, '-thetaCellCentered', thetaCellCentered, wasSet, ierr)
  call PetscOptionsGetBool(PETSC_NULL_CHARACTER, '-zetaCellCentered', zetaCellCentered, wasSet, ierr)

  if (masterProc) then
     print *,"Ntheta = ",Ntheta
     print *,"Nzeta = ",Nzeta
     print *,"Nxi = ",Nxi
     if (thetaCellCentered) then 
        print *,"thetaCellCentered = true"
     else
        print *,"thetaCellCentered = false"
     end if
     if (zetaCellCentered) then 
        print *,"zetaCellCentered = true"
     else
        print *,"zetaCellCentered = false"
     end if
  end if

  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)

  dof = 1
  stencilWidth = 1
  call DMDACreate3d(PETSC_COMM_WORLD, &
       DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE, & ! Options: DM_BOUNDARY_NONE, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_PERIODIC
       DMDA_STENCIL_STAR,Ntheta,Nzeta,Nxi, &
       PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, &
       dof, stencilWidth, &
       PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
       dmda, ierr)

  ! I'm not totally sure what this next line does - I think it sets piecewise-constant interpolation instead of piecewise-linear:
!  call DMDASetInterpolationType(dmda, DMDA_Q0, ierr)

  call KSPSetDM(ksp,dmda,ierr)
  call KSPSetComputeRHS(ksp,populateRHS,userContext,ierr)
  call KSPSetComputeOperators(ksp,populateMatrix,userContext,ierr)
  call KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, ierr)
  call KSPSetFromOptions(ksp,ierr)
  call KSPSolve(ksp,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

  call diagnostics(ksp)

  call KSPDestroy(ksp,ierr)
  call DMDestroy(dmda,ierr)

  call PETScFinalize(ierr)

end program
