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
  PetscScalar, pointer :: solnArray(:)
  PetscScalar, dimension(:), pointer :: xi, thetaWeights, zetaWeights, xiWeights
  PetscScalar, dimension(:,:), pointer :: B, dBdtheta, dBdzeta

  if (masterProc) then
     print *,"Entering diagnostics"
  end if

  ! For convenience, use some pointers to refer to quantities on the finest grid:
  xi => levels(1)%xi
  thetaWeights => levels(1)%thetaWeights
  zetaWeights  => levels(1)%zetaWeights
  xiWeights    => levels(1)%xiWeights
  B => levels(1)%B
  dBdtheta => levels(1)%dBdtheta
  dBdzeta => levels(1)%dBdzeta

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
              index = getIndex(1,itheta,izeta,ixi)+1

              flux = flux + solnArray(index) * (1+xi(ixi)*xi(ixi)) &
                   * (G * dBdtheta(itheta,izeta) - I * dBdzeta(itheta,izeta))/(B(itheta,izeta) ** 3) &
                   * thetaWeights(itheta)*zetaWeights(izeta)*xiWeights(ixi)

              flow = flow + solnArray(index) * xi(ixi) / B(itheta,izeta) &
                   * thetaWeights(itheta)*zetaWeights(izeta)*xiWeights(ixi)
           end do
        end do
     end do

     if (constraint_option==1) then
        print *,"Source = ",solnArray(levels(1)%matrixSize)
     end if
          
     call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)

     flow = flow * 2 / (sqrtpi*G*VPrime)
     flux = -2 / (sqrtpi*G*G*VPrime)*flux
     
     print *,"Results: VPrime = ",VPrime,", flux = ",flux,", flow = ",flow
  end if




  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mmc_soln.dat", FILE_MODE_WRITE, viewer, ierr)
  call VecView(solution, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)

end subroutine diagnostics
