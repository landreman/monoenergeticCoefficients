! Main program

!#include <finclude/petsckspdef.h>

program mmc

!  use petscksp
  use cyclic_red, only: rank, numtasks, bcyclic_solver, clearStorage, ns0, nsn
  use stel_kinds
  use variables
 
  implicit none

  INCLUDE 'mpif.h' 
  integer :: ierr, istat
  INTEGER, POINTER :: ipivot(:,:)

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
  numTasks = numProcs
  rank = myRank


  print *,"I am proc",myRank," of ",numProcs
  masterProc = (myRank==0)

  ! Set defaults:
  nu = 0.01d+0
  diagonalShift = nu * 1.0d-5
  epsilon_t = -0.07053d+0
  epsilon_h = 0.05067d+0
  iota = 0.4542d+0
  G = 3.7481d+0
  I = 0d+0
  Nperiods = 10
  helicity_l = 2
  Ntheta = 23
  Nzeta = 17
  Nxi = 33
  thetaGridScheme = 10
  zetaGridScheme = 10

  if (masterProc) then
     print *,"Ntheta = ",Ntheta
     print *,"Nzeta = ",Nzeta
     print *,"Nxi = ",Nxi
     print *,"nu = ",nu
     print *,"diagonalShift = ",diagonalShift
     print *,"thetaGridScheme    = ",thetaGridScheme
     print *,"zetaGridScheme     = ",zetaGridScheme
  end if

  call createGrids()

  ALLOCATE (lblk(mblock,mblock,ns0:nsn), dblk(mblock,mblock,ns0:nsn),   &
       ublk(mblock,mblock,ns0:nsn), brhs(mblock,ns0:nsn),          &
       ipivot(mblock, ns0:nsn), stat=istat)

  IF (istat .ne. 0) STOP 'Allocation error!'

  call populateRHS()

  call populateMatrix()

  CALL BCYCLIC_SOLVER (lblk, dblk, ublk, ipivot, brhs, mblock, Nxi)
 
  call diagnostics()

  call ClearStorage()
  DEALLOCATE (ipivot, lblk, dblk, ublk, brhs, stat=istat)

  call MPI_Finalize(ierr)

end program
