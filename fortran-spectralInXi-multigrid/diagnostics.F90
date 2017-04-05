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
  PetscInt :: itheta, izeta, L, index
  PetscScalar :: spatial_part
  PetscScalar, pointer :: solnArray(:)
  PetscScalar, dimension(:), pointer :: thetaWeights, zetaWeights
  PetscScalar, dimension(:,:), pointer :: B, dBdtheta, dBdzeta

  if (masterProc) then
     print *,"Entering diagnostics"
  end if

  ! For convenience, use some pointers to refer to quantities on the finest grid:
  thetaWeights => levels(1)%thetaWeights
  zetaWeights  => levels(1)%zetaWeights
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
           spatial_part = (G * dBdtheta(itheta,izeta) - I * dBdzeta(itheta,izeta))/(B(itheta,izeta) ** 3) &
                * thetaWeights(itheta)*zetaWeights(izeta)

           L = 0
           index = getIndex(1,itheta,izeta,L+1)+1 ! We add 1 here to convert from petsc 0-based indices to fortran 1-based indices:
           flux = flux + solnArray(index) * spatial_part * 8/(3.0d+0)
           
           L = 2
           index = getIndex(1,itheta,izeta,L+1)+1 ! We add 1 here to convert from petsc 0-based indices to fortran 1-based indices:
           flux = flux + solnArray(index) * spatial_part * 4/(15.0d+0)
           
           L = 1
           index = getIndex(1,itheta,izeta,L+1)+1 ! We add 1 here to convert from petsc 0-based indices to fortran 1-based indices:
           flow = flow + solnArray(index) / B(itheta,izeta) * thetaWeights(itheta)*zetaWeights(izeta)
        end do
     end do

     if (constraint_option==1) then
        print *,"Source = ",solnArray(levels(1)%matrixSize)
     end if
          
     call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)

     flow = flow * 4 / (3*sqrtpi*G*VPrime)
     flux = -2 / (sqrtpi*G*G*VPrime)*flux
     
     print *,"VPrime = ",VPrime,", FSAB2 = ",FSAB2
     print *,"Results: flux = ",flux,", flow = ",flow
  end if




  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mmc_soln.dat", FILE_MODE_WRITE, viewer, ierr)
  call VecView(solution, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)

end subroutine diagnostics
