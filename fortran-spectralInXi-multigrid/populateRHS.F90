#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscvecdef.h>
#else
#include <petsc/finclude/petscvecdef.h>
#endif

subroutine populateRHS(vec)

  use petscvec

  use indices
  use variables, only: G, I, masterProc, levels, Ntheta, Nzeta, Nxi, L_scaling

  implicit none

  Vec :: vec
  PetscErrorCode :: ierr

  PetscInt :: itheta, izeta, L, index, vec_size
  PetscReal :: valueToInsert
  PetscViewer :: viewer

  PetscScalar, dimension(:,:), pointer :: B, dBdtheta, dBdzeta

  if (masterProc) print *,"Entering populateRHS"

  ! For convenience, use some pointers to refer to quantities on the finest grid:
  B => levels(1)%B
  dBdtheta => levels(1)%dBdtheta
  dBdzeta => levels(1)%dBdzeta

  call VecSet(vec, 0.0d+0, ierr)
  call VecGetSize(vec, vec_size, ierr)

  if (masterProc) then
     ! For simplicity, do everything on proc 1.

     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           ! This next line is where the physics is:
           valueToInsert = (1.0d+0)/B(itheta,izeta)*(G*dBdtheta(itheta,izeta)-I*dBdzeta(itheta,izeta))

           L=0
           index = getIndex(1,itheta,izeta,L+1)
           call VecSetValue(vec, index, valueToInsert*4/(3.0d+0) * L_scaling(L+1), INSERT_VALUES, ierr)

           L=2
           index = getIndex(1,itheta,izeta,L+1)
           call VecSetValue(vec, index, valueToInsert*2/(3.0d+0) * L_scaling(L+1), INSERT_VALUES, ierr)
        end do
     end do
  end if

  call VecAssemblyBegin(vec, ierr)
  call VecAssemblyEnd(vec, ierr)

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mmc_rhs.dat", FILE_MODE_WRITE, viewer, ierr)
  call VecView(vec, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)

  if (masterProc) print *,"Leaving populateRHS"

end subroutine populateRHS
