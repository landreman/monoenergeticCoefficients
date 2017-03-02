#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

subroutine populateRHS(vec)

  use petscvec

  use indices
  use variables, only: G, I, masterProc, B, dBdtheta, dBdzeta, xi, Ntheta, Nzeta, Nxi

  implicit none


  Vec :: vec
  PetscErrorCode :: ierr

  PetscInt :: itheta, izeta, ixi, index
  PetscReal :: valueToInsert
  PetscViewer :: viewer

  call VecSet(vec, 0.0d+0, ierr)

  if (masterProc) then
     ! For simplicity, do everything on proc 1.

     print *,"Entering populateRHS"
  
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           do ixi = 1,Nxi
              ! This next line is where the physics is:
              valueToInsert = (1+xi(ixi)*xi(ixi))/B(itheta,izeta)*(G*dBdtheta(itheta,izeta)-I*dBdzeta(itheta,izeta))
              index = getIndex(itheta,izeta,ixi)
              call VecSetValue(vec, index, valueToInsert, INSERT_VALUES, ierr)
           end do
        end do
     end do
  end if

  call VecAssemblyBegin(vec, ierr)
  call VecAssemblyEnd(vec, ierr)

  !call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mmc_rhs.dat", FILE_MODE_WRITE, viewer, ierr)
  !call VecView(vec, viewer, ierr)
  !call PetscViewerDestroy(viewer, ierr)

end subroutine populateRHS
