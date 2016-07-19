#include <petsc/finclude/petsckspdef.h>


subroutine diagnostics(solution)

  use petscksp
  use FourierConvolutionMatrixMod
  use FourierTransformMod
  use indices
  use variables

  implicit none

  Vec :: solution
  
  PetscErrorCode :: ierr
  VecScatter :: VecScatterContext
  Vec :: solnOnProc0
  !PetscViewer :: viewer
  PetscInt :: L, index
  PetscReal :: flux, flow, VPrime
  PetscScalar, pointer :: solnArray(:)
  PetscScalar, dimension(:), allocatable :: tempFourierVector
  PetscScalar, dimension(:,:), allocatable :: tempFourierMatrix
  integer :: imn

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
     allocate(tempFourierVector(NFourier2))
     allocate(tempFourierMatrix(NFourier2,NFourier2))

     ! Convert the PETSc vector into a normal Fortran array:
     call VecGetArrayF90(solnOnProc0, solnArray, ierr)

     L = 0
     ! We add 1 here to convert from petsc 0-based indices to fortran 1-based indices:
     call FourierTransform((G*dBdtheta-I*dBdzeta)/(B*B*B), tempFourierVector)
     call FourierConvolutionMatrix(tempFourierVector,tempFourierMatrix)
     tempFourierVector=0
     do imn = 1,NFourier2
        index = getIndex(imn,L+1)+1
        tempFourierVector(imn) = solnArray(index)
     end do
     flux = dot_product(tempFourierMatrix(1,:), tempFourierVector)*4*pi*pi * 8/3

     L = 2
     ! We add 1 here to convert from petsc 0-based indices to fortran 1-based indices:
     do imn = 1,NFourier2
        index = getIndex(imn,L+1)+1
        tempFourierVector(imn) = solnArray(index)
     end do
     flux = flux + dot_product(tempFourierMatrix(1,:), tempFourierVector)*4*pi*pi * 4/15

     L = 1
     ! We add 1 here to convert from petsc 0-based indices to fortran 1-based indices:
     call FourierTransform(1/B, tempFourierVector)
     call FourierConvolutionMatrix(tempFourierVector,tempFourierMatrix)
     do imn = 1,NFourier2
        index = getIndex(imn,L+1)+1
        tempFourierVector(imn) = solnArray(index)
     end do
     ! Take the L=1 moment of f, divide by B, and integrate over theta and zeta:
     !flow = flow + solnArray(index) / B(itheta,izeta) * thetaWeights(itheta)*zetaWeights(izeta)
     flow = dot_product(tempFourierMatrix(1,:), tempFourierVector)*4*pi*pi
          
     call VecRestoreArrayF90(solnOnProc0, solnArray, ierr)

     call FourierTransform(1/(B*B), tempFourierVector)
     call FourierConvolutionMatrix(tempFourierVector,tempFourierMatrix)
     VPrime = tempFourierVector(1)*4*pi*pi

     flow = flow * 4 / (3*sqrtpi*G*VPrime)
     flux = -2 / (sqrtpi*G*G*VPrime)*flux
     
     deallocate(tempFourierVector, tempFourierMatrix)

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
        write (unit=fileUnit,fmt="(a,i5)") "NFourier  = ", NFourier
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
