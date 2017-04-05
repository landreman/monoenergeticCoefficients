#include <petsc/finclude/petscvecdef.h>

subroutine populateRHS(vec)

  use petscvec

  use indices
  use variables, only: G, I, masterProc, B, dBdtheta, dBdzeta, Ntheta, Nzeta, Nxi

  implicit none


  Vec :: vec
  PetscErrorCode :: ierr

  PetscInt :: itheta, izeta, L, index
  PetscReal :: valueToInsert
  PetscViewer :: viewer

  call VecSet(vec, 0.0d+0, ierr)

  if (masterProc) then
     ! For simplicity, do everything on proc 1.

     print *,"Entering populateRHS"
  
     do itheta = 1,Ntheta
        do izeta = 1,Nzeta
           valueToInsert = (1.0d+0)/B(itheta,izeta)*(G*dBdtheta(itheta,izeta)-I*dBdzeta(itheta,izeta))

           L=0
           index = getIndex(itheta,izeta,L+1)
           call VecSetValue(vec, index, valueToInsert*4/(3.0d+0), INSERT_VALUES, ierr)

           L=2
           index = getIndex(itheta,izeta,L+1)
           call VecSetValue(vec, index, valueToInsert*2/(3.0d+0), INSERT_VALUES, ierr)

        end do
     end do
  end if

  call VecAssemblyBegin(vec, ierr)
  call VecAssemblyEnd(vec, ierr)

!!$  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mmc_rhs.dat", FILE_MODE_WRITE, viewer, ierr)
!!$  call VecView(vec, viewer, ierr)
!!$  call PetscViewerDestroy(viewer, ierr)

end subroutine populateRHS
