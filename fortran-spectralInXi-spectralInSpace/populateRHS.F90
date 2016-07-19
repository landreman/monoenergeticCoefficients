#include <petsc/finclude/petscvecdef.h>

subroutine populateRHS(vec)

  use petscvec

  use FourierTransformMod
  use indices
  use variables, only: G, I, masterProc, B, dBdtheta, dBdzeta, NFourier2

  implicit none


  Vec :: vec
  PetscErrorCode :: ierr

  PetscInt :: imn, L, index
  !PetscViewer :: viewer
  PetscReal, dimension(:), allocatable :: tempFourierVector

  call VecSet(vec, 0.0d+0, ierr)

  if (masterProc) then
     ! For simplicity, do everything on proc 1.

     print *,"Entering populateRHS"
  
     allocate(tempFourierVector(NFourier2))

     call FourierTransform((G*dBdtheta - I*dBdzeta)/B, tempFourierVector)

     do imn = 1,NFourier2
        !valueToInsert = (1.0d+0)/B(itheta,izeta)*(G*dBdtheta(itheta,izeta)-I*dBdzeta(itheta,izeta))

        L=0
        index = getIndex(imn,L+1)
        call VecSetValue(vec, index, tempFourierVector(imn)*4/(3.0d+0), INSERT_VALUES, ierr)

        L=2
        index = getIndex(imn,L+1)
        call VecSetValue(vec, index, tempFourierVector(imn)*2/(3.0d+0), INSERT_VALUES, ierr)

     end do
     deallocate(tempFourierVector)

  end if

  call VecAssemblyBegin(vec, ierr)
  call VecAssemblyEnd(vec, ierr)

!!$  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mmc_rhs.dat", FILE_MODE_WRITE, viewer, ierr)
!!$  call VecView(vec, viewer, ierr)
!!$  call PetscViewerDestroy(viewer, ierr)

end subroutine populateRHS
