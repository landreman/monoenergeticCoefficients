#include <finclude/petscdmdadef.h>
#include <finclude/petsckspdef.h>

subroutine populateMatrix(ksp, matrix, pcMatrix, userContext, ierr)

  use petscdmda
  use petscksp

  use geometry
  use variables, only: nu, iota, Nperiods, myRank, pi

  implicit none

!#include <finclude/petscsys.h>
!#include <finclude/petscvec.h>
!#include <finclude/petscmat.h>
!#include <finclude/petscdm.h>
!#include <finclude/petscdmda.h>

  KSP :: ksp
  Mat :: matrix, pcMatrix
  Vec :: inputVec
  PetscInt :: userContext
  PetscErrorCode :: ierr

  DM :: dmda
  PetscInt :: ithetaMin, izetaMin, ixiMin
  PetscInt :: ithetaMax, izetaMax, ixiMax
  PetscInt :: localNtheta, localNzeta, localNxi
  PetscInt :: levelNtheta, levelNzeta, levelNxi
  PetscInt :: itheta, izeta, ixi
  PetscReal :: theta, zeta, xi
  MatNullSpace :: nullspace
  PetscReal :: dtheta, dzeta, dxi, BHere, temp

  PetscScalar :: valuesToAdd(7)
  MatStencil :: row(4), col(4,7)

  print *,"[",myRank,"] Entering populateMatrix"

  call KSPGetDM(ksp,dmda,ierr)

  ! Get Ntheta, Nzeta, and Nxi for this level of multigrid:
  call DMDAGetInfo(dmda,PETSC_NULL_INTEGER,levelNtheta,levelNzeta,levelNxi, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER,ierr)

  call DMDAGetCorners(dmda, &
       ithetaMin, izetaMin, ixiMin, &
       localNtheta, localNzeta, localNxi, &
       ierr)

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

  print *,"[proc ",myRank,"] matrix: ithetaMin=",ithetaMin,", izetaMin=",izetaMin,", ixiMin=",ixiMin,&
       ", localNtheta=",localNtheta,", localNzeta=",localNzeta,", localNxi=",localNxi


  dtheta = 2*pi/levelNtheta
  dzeta = 2*pi/(Nperiods*levelNzeta)
  dxi = (2.0d+0)/levelNxi

  do itheta = ithetaMin, ithetaMax
     do izeta = izetaMin, izetaMax
        do ixi = ixiMin, ixiMax
           row(MatStencil_i) = itheta
           row(MatStencil_j) = izeta
           row(MatStencil_k) = ixi

           col(MatStencil_i,1) = itheta
           col(MatStencil_j,1) = izeta
           col(MatStencil_k,1) = ixi

           col(MatStencil_i,2) = itheta+1
           col(MatStencil_j,2) = izeta
           col(MatStencil_k,2) = ixi

           col(MatStencil_i,3) = itheta-1
           col(MatStencil_j,3) = izeta
           col(MatStencil_k,3) = ixi

           col(MatStencil_i,4) = itheta
           col(MatStencil_j,4) = izeta+1
           col(MatStencil_k,4) = ixi

           col(MatStencil_i,5) = itheta
           col(MatStencil_j,5) = izeta-1
           col(MatStencil_k,5) = ixi

           col(MatStencil_i,6) = itheta
           col(MatStencil_j,6) = izeta
           col(MatStencil_k,6) = ixi+1

           col(MatStencil_i,7) = itheta
           col(MatStencil_j,7) = izeta
           col(MatStencil_k,7) = ixi-1

           valuesToAdd = 0

           call whereAmI(itheta,izeta,ixi,levelNtheta,levelNzeta,levelNxi,theta,zeta,xi)
           BHere = B(theta,zeta)

           ! Add d/dtheta term:
           valuesToAdd(2) =  xi*iota*BHere/(2*dtheta)
           valuesToAdd(3) = -xi*iota*BHere/(2*dtheta)

           ! Add d/dzeta term:
           valuesToAdd(4) =  xi*BHere/(2*dzeta)
           valuesToAdd(5) = -xi*BHere/(2*dzeta)
           
           ! Add mirror term:
           temp = -(0.5d+0)*(1-xi*xi)*(dBdtheta(theta,zeta)*iota + dBdzeta(theta,zeta))
           valuesToAdd(6) =  temp/(2*dxi)
           valuesToAdd(7) = -temp/(2*dxi)

           ! Add collision operator:
           valuesToAdd(6) = valuesToAdd(6) - nu/2*(1-(xi+dxi)*(xi+dxi))/(dxi*dxi)
           valuesToAdd(1) = nu/2*(    (1-(xi+dxi)*(xi+dxi)) + (1-(xi-dxi)*(xi-dxi))   )/(dxi*dxi)
           valuesToAdd(7) = valuesToAdd(7) - nu/2*(1-(xi-dxi)*(xi-dxi))/(dxi*dxi)

           call MatSetValuesStencil(pcMatrix, 1, row, 7, col, valuesToAdd, INSERT_VALUES, ierr)
        end do
     end do
  end do

  call MatAssemblyBegin(pcMatrix, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(pcMatrix, MAT_FINAL_ASSEMBLY, ierr)

  ! The matrix has a 1D null space, with the null vector corresponding to a constant:
  call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
  call MatSetNullSpace(pcMatrix,nullspace,ierr)
  call MatNullSpaceDestroy(nullspace,ierr)

end subroutine populateMatrix

