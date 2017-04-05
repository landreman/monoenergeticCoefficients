#include <petsc/finclude/petsckspdef.h>
#include <petsc/finclude/petscdmdadef.h>


subroutine diagnostics(ksp)

  use petscksp
  use petscdmda

  use geometry
  use variables

  implicit none

  KSP :: ksp

  PetscErrorCode ierr
  Vec :: solution
  DM :: dmda
  PetscViewer :: viewer
  PetscInt :: ithetaMin, izetaMin, ixiMin
  PetscInt :: ithetaMax, izetaMax, ixiMax
  PetscInt :: localNtheta, localNzeta, localNxi
  PetscInt :: levelNtheta, levelNzeta, levelNxi
  PetscInt :: itheta, izeta, ixi
  PetscReal :: theta, zeta, xi
  PetscReal :: dtheta, dzeta, dxi, Vprime
  PetscScalar, pointer :: array(:,:,:)

  PetscReal :: flux, flow
  double precision :: sendBuffer(1), recvBuffer(1)

  call KSPGetSolution(ksp,solution,ierr)

  print *,"[",myRank,"] Entering diagnostics"
  
  call KSPGetDM(ksp,dmda,ierr)

  ! Get Ntheta, Nzeta, and Nxi for this level of multigrid:
  call DMDAGetInfo(dmda,PETSC_NULL_INTEGER,levelNtheta,levelNzeta,levelNxi, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER,ierr)

  print *,"[proc ",myRank,"] diagnostics: levelNtheta=",levelNtheta,", levelNzeta=",levelNzeta,", levelNxi=",levelNxi

  if ((levelNtheta .ne. Ntheta) .or. (levelNzeta .ne. Nzeta) .or. (levelNxi .ne. Nxi)) then
     print *,"Error! The final dmda has different dimensions than the requested dimensions."
     stop
  end if

  call DMDAGetCorners(dmda, &
       ithetaMin, izetaMin, ixiMin, &
       localNtheta, localNzeta, localNxi, &
       ierr)

  call DMDAVecGetArrayF90(dmda, solution, array, ierr)

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

  print *,"[proc ",myRank,"] diagnostics: ithetaMin=",ithetaMin,", izetaMin=",izetaMin,", ixiMin=",ixiMin,&
       ", localNtheta=",localNtheta,", localNzeta=",localNzeta,", localNxi=",localNxi
  
  flux = 0
  flow = 0

  do itheta = ithetaMin, ithetaMax
     do izeta = izetaMin, izetaMax
        do ixi = ixiMin, ixiMax
           call whereAmI(itheta,izeta,ixi,levelNtheta,levelNzeta,levelNxi,theta,zeta,xi)

           flux = flux + array(itheta,izeta,ixi) * (1+xi*xi) &
                * (G * dBdtheta(theta,zeta) - I * dBdzeta(theta,zeta))/(B(theta,zeta) ** 3)

           flow = flow + array(itheta,izeta,ixi) * xi / B(theta,zeta)
        end do
     end do
  end do

  ! Include the integration weights:
  dtheta = 2*pi/Ntheta
  dzeta = 2*pi/(Nperiods*Nzeta)
  dxi = (2.0d+0)/Nxi

  flux = flux*dtheta*dzeta*dxi
  flow = flow*dtheta*dzeta*dxi

  VPrime = 0
  ixi=0
  do itheta = 0,(Ntheta-1)
     do izeta = 0,(Nzeta-1)
        call whereAmI(itheta,izeta,ixi,levelNtheta,levelNzeta,levelNxi,theta,zeta,xi)
        VPrime = VPrime + dtheta*dzeta/(B(theta,zeta) ** 2)
     end do
  end do

  flow = flow * 2 / (sqrtpi*G*VPrime)
  flux = -2 / (sqrtpi*G*G*VPrime)*flux

  print *,"[",myRank,"] Results: VPrime = ",VPrime,", flux = ",flux,", flow = ",flow

  ! Sum results from all procs:
  sendBuffer = flux
  call MPI_Reduce(sendBuffer, recvBuffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  flux = recvBuffer(1)

  sendBuffer = flow
  call MPI_Reduce(sendBuffer, recvBuffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  flow = recvBuffer(1)

  if (masterProc) then
     print *,"Final results: flux = ",flux,", flow = ",flow
  end if

  call DMDAVecRestoreArrayF90(dmda, solution, array, ierr)

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mmc_soln.dat", FILE_MODE_WRITE, viewer, ierr)
  call VecView(solution, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)

end subroutine diagnostics
