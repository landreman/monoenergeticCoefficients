#include <finclude/petsckspdef.h>
!#include <finclude/petscdm.h>
#include <finclude/petscdmdadef.h>
!#include <finclude/petscdmda.h90>

subroutine populateRHS(ksp, vec, userContext, ierr)

  use petscksp
  use petscdmda

  use geometry
  use variables, only: G, I, myRank, masterProc

  implicit none


  KSP :: ksp
  Vec :: vec
  PetscInt :: userContext
  PetscErrorCode :: ierr

  DM :: dmda
  PetscInt :: ithetaMin, izetaMin, ixiMin
  PetscInt :: ithetaMax, izetaMax, ixiMax
  PetscInt :: localNtheta, localNzeta, localNxi
  PetscInt :: levelNtheta, levelNzeta, levelNxi
  PetscInt :: itheta, izeta, ixi
  PetscReal :: theta, zeta, xi, valueToInsert
  PetscScalar, pointer :: array(:,:,:)
  PetscViewer :: viewer

  print *,"[",myRank,"] Entering populateRHS"
  
  call KSPGetDM(ksp,dmda,ierr)

  ! Get Ntheta, Nzeta, and Nxi for this level of multigrid:
  call DMDAGetInfo(dmda,PETSC_NULL_INTEGER,levelNtheta,levelNzeta,levelNxi, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER,ierr)

  print *,"[proc ",myRank,"] RHS: levelNtheta=",levelNtheta,", levelNzeta=",levelNzeta,", levelNxi=",levelNxi

  call DMDAGetCorners(dmda, &
       ithetaMin, izetaMin, ixiMin, &
       localNtheta, localNzeta, localNxi, &
       ierr)

  call DMDAVecGetArrayF90(dmda, vec, array, ierr)

  ithetaMax = ithetaMin + localNtheta - 1
  izetaMax  = izetaMin  + localNzeta  - 1
  ixiMax    = ixiMin    + localNxi    - 1

!!$  ! Switch to 1-based indices:
!!$  ithetaMin = ithetaMin + 1
!!$  ithetaMax = ithetaMin + localNtheta - 1
!!$  izetaMin = izetaMin + 1
!!$  izetaMax = izetaMin + localNzeta - 1
!!$  ixiMin = ixiMin + 1
!!$  ixiMax = ixiMin + localNxi - 1

  print *,"[proc ",myRank,"] RHS: ithetaMin=",ithetaMin,", izetaMin=",izetaMin,", ixiMin=",ixiMin,&
       ", localNtheta=",localNtheta,", localNzeta=",localNzeta,", localNxi=",localNxi
  
  if (masterProc) then
     print *,"Here comes theta:"
     izeta=1
     ixi=1
     do itheta=ithetaMin,ithetaMax
        call whereAmI(itheta,izeta,ixi,levelNtheta,levelNzeta,levelNxi,theta,zeta,xi)
        print *,theta
     end do

     print *,"Here comes zeta:"
     itheta=1
     ixi=1
     do izeta=izetaMin,izetaMax
        call whereAmI(itheta,izeta,ixi,levelNtheta,levelNzeta,levelNxi,theta,zeta,xi)
        print *,zeta
     end do

     print *,"Here comes xi:"
     itheta=1
     izeta=1
     do ixi=ixiMin,ixiMax
        call whereAmI(itheta,izeta,ixi,levelNtheta,levelNzeta,levelNxi,theta,zeta,xi)
        print *,xi
     end do
  end if

  do itheta = ithetaMin, ithetaMax
     do izeta = izetaMin, izetaMax
        do ixi = ixiMin, ixiMax
           call whereAmI(itheta,izeta,ixi,levelNtheta,levelNzeta,levelNxi,theta,zeta,xi)

           ! This next line is where the physics is:
           valueToInsert = (1+xi*xi)/B(theta,zeta)*(G*dBdtheta(theta,zeta)-I*dBdzeta(theta,zeta))

           ! This next line requires 0-based indices!
           ! See Fortran note in http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMDAVecGetArray.html
           array(itheta,izeta,ixi) = valueToInsert
        end do
     end do
  end do

  call DMDAVecRestoreArrayF90(dmda, vec, array, ierr)
  call VecAssemblyBegin(vec, ierr)
  call VecAssemblyEnd(vec, ierr)

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mmc_rhs.dat", FILE_MODE_WRITE, viewer, ierr)
  call VecView(vec, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)

end subroutine populateRHS
