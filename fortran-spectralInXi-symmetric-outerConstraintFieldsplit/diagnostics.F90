#include <petsc/finclude/petsckspdef.h>


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
  PetscReal :: flux, flow, VPrime_old, spatialPart
  PetscScalar, pointer :: solnArray(:)

  integer :: clockStop
  real :: elapsedTime

  character(len=6) :: filename="output"
  integer :: fileUnit=11, didFileAccessWork


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
           spatialPart = (G * dBdtheta(itheta,izeta) - I * dBdzeta(itheta,izeta))/(B(itheta,izeta) ** 3) &
                * thetaWeight*zetaWeight

           L = 0
           ! We add 1 here to convert from petsc 0-based indices to fortran 1-based indices:
           index = getIndex(itheta,izeta,L+1)+1
           flux = flux + solnArray(index) * spatialPart * (8/3.0d+0) * L_scaling(L+1)

           L = 2
           ! We add 1 here to convert from petsc 0-based indices to fortran 1-based indices:
           index = getIndex(itheta,izeta,L+1)+1
           flux = flux + solnArray(index) * spatialPart * (4/15.0d+0) * L_scaling(L+1)

           L = 1
           ! We add 1 here to convert from petsc 0-based indices to fortran 1-based indices:
           index = getIndex(itheta,izeta,L+1)+1
           flow = flow + solnArray(index) / B(itheta,izeta) * thetaWeight*zetaWeight * L_scaling(L+1)
        end do
     end do
          
     call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)

     VPrime_old = thetaWeight * zetaWeight * sum(1/(B*B))
     flow = flow * 4 / (3*sqrtpi*G*VPrime_old)
     flux = -2 / (sqrtpi*G*G*VPrime_old)*flux
     
     call system_clock(clockStop)
     elapsedTime = real(clockStop-clockStart)/clockRate

     print *,"Results:"
     print *,"  Flux = ",flux
     print *,"  Flow = ",flow
     print *,"  Time for solve (seconds) = ",elapsedTime

     open(unit=fileUnit, file=filename, action="write", iostat = didFileAccessWork)
     if (didFileAccessWork /= 0) then
        print *,"Error opening ", trim(filename)
        stop
     else
        write (unit=fileUnit,fmt="(a,i5)") "Ntheta = ", Ntheta
        write (unit=fileUnit,fmt="(a,i5)") "Nzeta  = ", Nzeta
        write (unit=fileUnit,fmt="(a,i5)") "Nxi    = ", Nxi
        write (unit=fileUnit,fmt="(a,es22.15)") "nu = ", nu
        write (unit=fileUnit,fmt="(a,i5)") "numProcs = ", numProcs
        write (unit=fileUnit,fmt="(a,es22.15)") "Flux = ", flux
        write (unit=fileUnit,fmt="(a,es22.15)") "Flow = ", flow
        write (unit=fileUnit,fmt="(a,es22.15)") "Time = ", elapsedTime
        close(unit=fileUnit)
     end if

  end if




!!$  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mmc_soln.dat", FILE_MODE_WRITE, viewer, ierr)
!!$  call VecView(solution, viewer, ierr)
!!$  call PetscViewerDestroy(viewer, ierr)

end subroutine diagnostics
