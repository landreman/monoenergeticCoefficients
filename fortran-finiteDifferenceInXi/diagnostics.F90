#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petsckspdef.h>
#else
#include <petsc/finclude/petsckspdef.h>
#endif


subroutine diagnostics(solution)

  use petscksp

  use variables
  use indices

  implicit none

  Vec :: solution
  
  PetscErrorCode :: ierr
  VecScatter :: VecScatterContext
  Vec :: solnOnProc0
  PetscViewer :: viewer
  PetscInt :: itheta, izeta, ixi, index
  PetscReal :: flux, flow, VPrime
  PetscScalar, pointer :: solnArray(:)

  if (masterProc) then
     print *,"Entering diagnostics"
  end if

  flux = 0
  flow = 0

  ! Send the entire solution vector to the master process:
  call VecScatterCreateToZero(solution, VecScatterContext, solnOnProc0, ierr)
  call VecScatterBegin(VecScatterContext, solution, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(VecScatterContext, solution, solnOnProc0, INSERT_VALUES, SCATTER_FORWARD, ierr)
       
  if (masterProc) then
     ! Convert the PETSc vector into a normal Fortran array:
     call VecGetArrayF90(solnOnProc0, solnArray, ierr)

     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           do ixi = 1,Nxi

              ! We add 1 here to convert from petsc 0-based indices to fortran 1-based indices:
              index = getIndex(itheta,izeta,ixi)+1

              flux = flux + solnArray(index) * (1+xi(ixi)*xi(ixi)) &
                   * (G * dBdtheta(itheta,izeta) - I * dBdzeta(itheta,izeta))/(B(itheta,izeta) ** 3) &
                   * thetaWeights(itheta)*zetaWeights(izeta)*xiWeights(ixi)

              flow = flow + solnArray(index) * xi(ixi) / B(itheta,izeta) &
                   * thetaWeights(itheta)*zetaWeights(izeta)*xiWeights(ixi)
           end do
        end do
     end do
          
     call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)

     VPrime = 0
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           VPrime = VPrime + thetaWeights(itheta)*zetaWeights(izeta)/(B(itheta,izeta) ** 2)
        end do
     end do
     
     flow = flow * 2 / (sqrtpi*G*VPrime)
     flux = -2 / (sqrtpi*G*G*VPrime)*flux
     
     print *,"Results: VPrime = ",VPrime,", flux = ",flux,", flow = ",flow
  end if




  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mmc_soln.dat", FILE_MODE_WRITE, viewer, ierr)
  call VecView(solution, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)

end subroutine diagnostics
