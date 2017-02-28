#include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
#include <finclude/petscdmdadef.h>
#else
#include <petsc/finclude/petscdmdadef.h>
#endif

  subroutine create_grids(level)
    use petscdmda
       
    use variables, Ntheta_fine => Ntheta, Nzeta_fine => Nzeta, Nxi_fine => Nxi, matrixSize_fine => matrixSize


    implicit none

    integer, intent(in) :: level

    PetscScalar, dimension(:), allocatable :: theta_preconditioner, thetaWeights_preconditioner
    PetscScalar, dimension(:), allocatable :: zeta_preconditioner, zetaWeights_preconditioner
    PetscScalar, dimension(:,:), allocatable :: d2dtheta2
    PetscScalar, dimension(:,:), allocatable :: d2dzeta2, temp_matrix
    PetscInt :: scheme, quadrature_option, derivative_option_plus, derivative_option_minus, derivative_option
    PetscInt :: j, k, itheta, izeta
    PetscInt :: localNtheta, localNzeta
    DM :: myDM
    integer, parameter :: bufferLength = 200
    character(len=bufferLength) :: procAssignments
    integer :: tag, dummy(1)
    integer :: status(MPI_STATUS_SIZE)
    logical :: call_uniform_diff_matrices
    PetscErrorCode :: ierr

    ! For convenience, use some short variable names to refer to quantities on this level:
    integer, pointer :: Ntheta, Nzeta, Nxi, matrixSize
    integer, pointer :: ithetaMin, ithetaMax, izetaMin, izetaMax
    PetscScalar, dimension(:), pointer :: theta, zeta, thetaWeights, zetaWeights
    PetscScalar, dimension(:,:), pointer :: B, dBdtheta, dBdzeta
    PetscScalar, dimension(:,:), pointer :: ddtheta_plus, ddtheta_minus, ddtheta_plus_preconditioner, ddtheta_minus_preconditioner
    PetscScalar, dimension(:,:), pointer :: ddzeta_plus, ddzeta_minus, ddzeta_plus_preconditioner, ddzeta_minus_preconditioner

    ! For convenience, use some short variable names to refer to quantities on this level:
    Ntheta => levels(level)%Ntheta
    Nzeta  => levels(level)%Nzeta
    Nxi    => levels(level)%Nxi
    matrixSize => levels(level)%matrixSize
    ithetaMin => levels(level)%ithetaMin
    ithetaMax => levels(level)%ithetaMax
    izetaMin  => levels(level)%izetaMin
    izetaMax  => levels(level)%izetaMax
    
    if (masterProc) print "(a,i3,a)"," ---- Initializing grids for multigrid level",level,"----"
    matrixSize = Ntheta*Nzeta*Nxi
    if (constraint_option==1) matrixSize = matrixSize + 1
    if (masterProc) print *,"matrixSize:",matrixSize
    
    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Set up ranges of indices owned by each processor.
    !
    ! *******************************************************************************
    ! *******************************************************************************

    ! Each processor is responsible for building the rows of the matrix and rhs corresponding
    ! to its ithetaMin:ithetaMax and izetaMin:izetaMax, and each processor is resposible for all columns of the matrix.

    ! In principle we could distribute in both theta and zeta at the same time.
    ! However, this would lead to negligible increase in speed, since the bottleneck is
    ! not matrix construction but rather the solve, which is parallelized in a completely different
    ! way (determined internally by superlu_dist or mumps.)
    if (Ntheta > Nzeta) then
       ! Distribute in theta but not in zeta

       ! Assign a range of theta indices to each processor.
       ! This is done by creating a PETSc DM that is not actually used for anything else.
       call DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, Ntheta, 1, 0, PETSC_NULL_INTEGER, myDM, ierr)

       call DMDAGetCorners(myDM, ithetaMin, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
            localNtheta, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

       call DMDestroy(myDM, ierr)

       izetaMin = 0
       izetaMax = Nzeta-1
       localNzeta = Nzeta
    else
       ! Distribute in zeta but not in theta

       ! Assign a range of zeta indices to each processor.
       ! This is done by creating a PETSc DM that is not actually used for anything else.
       call DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, Nzeta, 1, 0, PETSC_NULL_INTEGER, myDM, ierr)

       call DMDAGetCorners(myDM, izetaMin, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
            localNzeta, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

       call DMDestroy(myDM, ierr)
       ithetaMin = 0
       ithetaMax = Ntheta-1
       localNtheta = Ntheta
    end if

    ! Below is some code that breaks up the theta and zeta ranges at the same time.
    ! I'm commented it out because PETSc kept giving an error when the number of
    ! procs was large compared to Ntheta and Nzeta.
!!$    ! Assign a range of theta and zeta indices to each processor.
!!$    ! This is done by creating a PETSc DM that is not actually used for anything else.
!!$    call DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, &
!!$         Ntheta, Nzeta, PETSC_DECIDE, PETSC_DECIDE, 1, 0, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, myDM, ierr)
!!$
!!$    call DMDAGetCorners(myDM, ithetaMin, izetaMin, PETSC_NULL_INTEGER, &
!!$         localNtheta, localNzeta, PETSC_NULL_INTEGER, ierr)
!!$
!!$    call DMDestroy(myDM, ierr)

    ! Switch to 1-based indices:
    ithetaMin = ithetaMin + 1
    ithetaMax = ithetaMin+localNtheta-1
    izetaMin = izetaMin + 1
    izetaMax = izetaMin+localNzeta-1

    write (procAssignments,fmt="(a,i4,a,i3,a,i3,a,i3,a,i3,a)") "Processor ",myRank," owns theta indices ",ithetaMin," to ",ithetaMax,&
         " and zeta indices ",izetaMin," to ",izetaMax

!    call PetscSynchronizedPrintf(PETSC_COMM_WORLD, procAssignments, ierr)
!    call PetscSynchronizedFlush(PETSC_COMM_WORLD, ierr)

    ! PETSc's synchronized printing functions seem buggy, so here I've implemented my own version:
    dummy = 0
    tag = 0
    if (masterProc) then
       print *,trim(procAssignments)
       do j = 1,numProcs-1
          ! To avoid a disordered flood of messages to the masterProc,
          ! ping each proc 1 at a time by sending a dummy value:
          call MPI_SEND(dummy,1,MPI_INT,j,tag,PETSC_COMM_WORLD,ierr)
          ! Now receive the message from proc i:
          call MPI_RECV(procAssignments,bufferLength,MPI_CHAR,j,MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
          print *,trim(procAssignments)
       end do
    else
       ! First, wait for the dummy message from proc 0:
       call MPI_RECV(dummy,1,MPI_INT,0,MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
       ! Now send the message to proc 0:
       call MPI_SEND(procAssignments,bufferLength,MPI_CHAR,0,tag,PETSC_COMM_WORLD,ierr)
    end if



    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Build theta grid, integration weights, and differentiation matrices:
    !
    ! *******************************************************************************
    ! *******************************************************************************

    allocate(levels(level)%theta(Ntheta))
    allocate(levels(level)%thetaWeights(Ntheta))
    allocate(levels(level)%ddtheta_plus(Ntheta,Ntheta))
    allocate(levels(level)%ddtheta_minus(Ntheta,Ntheta))
    allocate(levels(level)%ddtheta_plus_preconditioner(Ntheta,Ntheta))
    allocate(levels(level)%ddtheta_minus_preconditioner(Ntheta,Ntheta))
    allocate(levels(level)%ddtheta_sum(Ntheta,Ntheta))
    allocate(levels(level)%ddtheta_difference(Ntheta,Ntheta))
    allocate(levels(level)%ddtheta_sum_preconditioner(Ntheta,Ntheta))
    allocate(levels(level)%ddtheta_difference_preconditioner(Ntheta,Ntheta))
    allocate(d2dtheta2(Ntheta,Ntheta))

    theta => levels(level)%theta
    thetaWeights => levels(level)%thetaWeights
    ddtheta_plus => levels(level)%ddtheta_plus
    ddtheta_minus => levels(level)%ddtheta_minus
    ddtheta_plus_preconditioner => levels(level)%ddtheta_plus_preconditioner
    ddtheta_minus_preconditioner => levels(level)%ddtheta_minus_preconditioner

    ! *******************************************************************************
    ! Handle d/dtheta for the main matrix.
    ! *******************************************************************************

    select case (theta_derivative_option)

    case (1)
       if (masterProc .and. level==1) then
          print *,"d/dtheta derivatives discretized using Fourier pseudospectral method."
       end if
       derivative_option_plus = 20
       derivative_option_minus = derivative_option_plus

    case (2)
       if (masterProc .and. level==1) then
          print *,"d/dtheta derivatives discretized using centered finite differences:"
          print *,"   1 point on either side."
       end if
       derivative_option_plus = 0
       derivative_option_minus = derivative_option_plus

    case (3)
       if (masterProc .and. level==1) then
          print *,"d/dtheta derivatives discretized using centered finite differences:"
          print *,"   2 points on either side."
       end if
       derivative_option_plus = 10
       derivative_option_minus = derivative_option_plus

    case (4)
       if (masterProc .and. level==1) then
          print *,"d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   0 points on one side, 1 point on the other side."
       end if
       derivative_option_plus  = 30
       derivative_option_minus = 40

    case (5)
       if (masterProc .and. level==1) then
          print *,"d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   0 points on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 50
       derivative_option_minus = 60

    case (6)
       if (masterProc .and. level==1) then
          print *,"d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   1 point on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 80
       derivative_option_minus = 90

    case (7)
       if (masterProc .and. level==1) then
          print *,"d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   1 point on one side, 3 points on the other side."
       end if
       derivative_option_plus  = 100
       derivative_option_minus = 110

    case (8)
       if (masterProc .and. level==1) then
          print *,"d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   2 point on one side, 3 points on the other side."
       end if
       derivative_option_plus  = 120
       derivative_option_minus = 130

    case default
       print *,"Error! Invalid setting for theta_derivative_option:",theta_derivative_option
       stop
    end select

    quadrature_option = 0
    call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_plus,  quadrature_option, theta, thetaWeights, ddtheta_plus,  d2dtheta2)
    call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_minus, quadrature_option, theta, thetaWeights, ddtheta_minus, d2dtheta2)

    levels(level)%ddtheta_sum        = (ddtheta_plus + ddtheta_minus)/two
    levels(level)%ddtheta_difference = (ddtheta_plus - ddtheta_minus)/two

    ! *******************************************************************************
    ! Handle d/dtheta for the preconditioner matrix.
    ! *******************************************************************************

    call_uniform_diff_matrices = .true.
    select case (abs(preconditioner_theta_derivative_option))
    case (0)
       if (masterProc .and. level==1) then
          print *,"d/dtheta terms are dropped in the preconditioner."
       end if
       ddtheta_plus_preconditioner = 0
       ddtheta_minus_preconditioner = 0
       call_uniform_diff_matrices = .false.

    case (100)
       if (masterProc .and. level==1) then
          print *,"d/dtheta matrices are the same in the preconditioner."
       end if
       ddtheta_plus_preconditioner  = ddtheta_plus
       ddtheta_minus_preconditioner = ddtheta_minus
       call_uniform_diff_matrices = .false.

    case (1)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dtheta derivatives discretized using Fourier pseudospectral method."
       end if
       derivative_option_plus = 20
       derivative_option_minus = derivative_option_plus

    case (2)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dtheta derivatives discretized using centered finite differences:"
          print *,"   1 point on either side."
       end if
       derivative_option_plus = 0
       derivative_option_minus = derivative_option_plus

    case (3)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dtheta derivatives discretized using centered finite differences:"
          print *,"   2 points on either side."
       end if
       derivative_option_plus = 10
       derivative_option_minus = derivative_option_plus

    case (4)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   0 points on one side, 1 point on the other side."
       end if
       derivative_option_plus  = 30
       derivative_option_minus = 40

    case (5)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   0 points on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 50
       derivative_option_minus = 60

    case (6)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   1 point on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 80
       derivative_option_minus = 90

    case (7)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   1 point on one side, 3 points on the other side."
       end if
       derivative_option_plus  = 100
       derivative_option_minus = 110

    case (8)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dtheta derivatives discretized using upwinded finite differences:"
          print *,"   2 point on one side, 3 points on the other side."
       end if
       derivative_option_plus  = 120
       derivative_option_minus = 130

    case default
       print *,"Error! Invalid setting for theta_derivative_option:",theta_derivative_option
       stop
    end select

    if (call_uniform_diff_matrices) then
       quadrature_option = 0
       call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_plus,  quadrature_option, theta, thetaWeights, ddtheta_plus_preconditioner,  d2dtheta2)
       call uniformDiffMatrices(Ntheta, zero, two*pi, derivative_option_minus, quadrature_option, theta, thetaWeights, ddtheta_minus_preconditioner, d2dtheta2)
    end if

    if (theta_diffusion>0) then
       if (masterProc .and. level==1) then
          print *,"  Adding diffusion in theta:",theta_diffusion
       end if
       derivative_option = 0
       quadrature_option = 0
       allocate(temp_matrix(Ntheta,Ntheta))
       call uniformDiffMatrices(Ntheta, zero, two*pi, &
            derivative_option,  quadrature_option, theta, thetaWeights, temp_matrix, d2dtheta2)
       ddtheta_plus_preconditioner  = ddtheta_plus_preconditioner  - theta_diffusion * d2dtheta2
       ddtheta_minus_preconditioner = ddtheta_minus_preconditioner - theta_diffusion * d2dtheta2
       deallocate(temp_matrix)
    end if

    if (preconditioner_theta_derivative_option<0) then
       if (masterProc .and. level==1) then
          print *,"   But only the diagonal is kept."
       end if
       do j=1,Ntheta
          do k=1,Ntheta
             if (j .ne. k) then
                ddtheta_plus_preconditioner(j,k) = 0
                ddtheta_minus_preconditioner(j,k) = 0
             end if
          end do
       end do
    end if
    
    levels(level)%ddtheta_sum_preconditioner        = (ddtheta_plus_preconditioner + ddtheta_minus_preconditioner)/two
    levels(level)%ddtheta_difference_preconditioner = (ddtheta_plus_preconditioner - ddtheta_minus_preconditioner)/two

    ! The following arrays will not be needed:
    deallocate(d2dtheta2)


    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Build zeta grid, integration weights, and differentiation matrices:
    !
    ! *******************************************************************************
    ! *******************************************************************************

    zetaMax = 2*pi/NPeriods

    allocate(levels(level)%zeta(Nzeta))
    allocate(levels(level)%zetaWeights(Nzeta))
    allocate(levels(level)%ddzeta_plus(Nzeta,Nzeta))
    allocate(levels(level)%ddzeta_minus(Nzeta,Nzeta))
    allocate(levels(level)%ddzeta_plus_preconditioner(Nzeta,Nzeta))
    allocate(levels(level)%ddzeta_minus_preconditioner(Nzeta,Nzeta))
    allocate(levels(level)%ddzeta_sum(Nzeta,Nzeta))
    allocate(levels(level)%ddzeta_difference(Nzeta,Nzeta))
    allocate(levels(level)%ddzeta_sum_preconditioner(Nzeta,Nzeta))
    allocate(levels(level)%ddzeta_difference_preconditioner(Nzeta,Nzeta))
    allocate(d2dzeta2(Nzeta,Nzeta))

    zeta => levels(level)%zeta
    zetaWeights  => levels(level)%zetaWeights
    ddzeta_plus => levels(level)%ddzeta_plus
    ddzeta_minus => levels(level)%ddzeta_minus
    ddzeta_plus_preconditioner => levels(level)%ddzeta_plus_preconditioner
    ddzeta_minus_preconditioner => levels(level)%ddzeta_minus_preconditioner

    if (Nzeta==1) then

       ! *******************************************************************************
       ! Axisymmetry is a special case:
       ! *******************************************************************************
       zeta = 0
       zetaWeights = 2*pi
       ddzeta_plus = 0
       ddzeta_minus = 0
       ddzeta_plus_preconditioner = 0
       ddzeta_minus_preconditioner = 0
       levels(level)%ddzeta_sum = 0
       levels(level)%ddzeta_difference = 0
       levels(level)%ddzeta_sum_preconditioner = 0
       levels(level)%ddzeta_difference_preconditioner = 0
       d2dzeta2 = 0 ! d2dzeta2 is not actually used.

    else

       ! *******************************************************************************
       ! Not axisymmetric.
       ! First, handle d/dzeta for the streaming term in the main matrix:
       ! *******************************************************************************

       select case (zeta_derivative_option)

       case (2)
          if (masterProc .and. level==1) then
             print *,"d/dzeta derivatives discretized using centered finite differences:"
             print *,"   1 point on either side."
          end if
          derivative_option_plus = 0
          derivative_option_minus = derivative_option_plus
          
       case (3)
          if (masterProc .and. level==1) then
             print *,"d/dzeta derivatives discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option_plus = 10
          derivative_option_minus = derivative_option_plus
          
       case (4)
          if (masterProc .and. level==1) then
             print *,"d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 1 point on the other side."
          end if
          derivative_option_plus  = 30
          derivative_option_minus = 40
          
       case (5)
          if (masterProc .and. level==1) then
             print *,"d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 50
          derivative_option_minus = 60
          
       case (6)
          if (masterProc .and. level==1) then
             print *,"d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 80
          derivative_option_minus = 90
          
       case (7)
          if (masterProc .and. level==1) then
             print *,"d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 100
          derivative_option_minus = 110
          
       case (8)
          if (masterProc .and. level==1) then
             print *,"d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 120
          derivative_option_minus = 130
          
       case default
          print *,"Error! Invalid setting for zeta_derivative_option:",zeta_derivative_option
          stop
       end select
       
       quadrature_option = 0
       call uniformDiffMatrices(Nzeta, zero, zetaMax, &
            derivative_option_plus,  quadrature_option, zeta, zetaWeights, ddzeta_plus, d2dzeta2)
       call uniformDiffMatrices(Nzeta, zero, zetaMax, &
            derivative_option_minus, quadrature_option, zeta, zetaWeights, ddzeta_minus, d2dzeta2)

       levels(level)%ddzeta_sum        = (ddzeta_plus + ddzeta_minus)/two
       levels(level)%ddzeta_difference = (ddzeta_plus - ddzeta_minus)/two

       ! *******************************************************************************
       ! Handle d/dzeta for the streaming term in the preconditioner matrix.
       ! *******************************************************************************

       call_uniform_diff_matrices = .true.
       select case (abs(preconditioner_zeta_derivative_option))
       case (0)
          if (masterProc .and. level==1) then
             print *,"d/dzeta terms are dropped in the preconditioner."
          end if
          ddzeta_plus_preconditioner = 0
          ddzeta_minus_preconditioner = 0
          call_uniform_diff_matrices = .false.
          
       case (100)
          if (masterProc .and. level==1) then
             print *,"d/dzeta matrices are the same in the preconditioner."
          end if
          ddzeta_plus_preconditioner  = ddzeta_plus
          ddzeta_minus_preconditioner = ddzeta_minus
          call_uniform_diff_matrices = .false.          
          
       case (2)
          if (masterProc .and. level==1) then
             print *,"Preconditioner d/dzeta derivatives discretized using centered finite differences:"
             print *,"   1 point on either side."
          end if
          derivative_option_plus = 0
          derivative_option_minus = derivative_option_plus
          
       case (3)
          if (masterProc .and. level==1) then
             print *,"Preconditioner d/dzeta derivatives discretized using centered finite differences:"
             print *,"   2 points on either side."
          end if
          derivative_option_plus = 10
          derivative_option_minus = derivative_option_plus

       case (4)
          if (masterProc .and. level==1) then
             print *,"Preconditioner d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 1 point on the other side."
          end if
          derivative_option_plus  = 30
          derivative_option_minus = 40
          
       case (5)
          if (masterProc .and. level==1) then
             print *,"Preconditioner d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   0 points on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 50
          derivative_option_minus = 60
          
       case (6)
          if (masterProc .and. level==1) then
             print *,"Preconditioner d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 2 points on the other side."
          end if
          derivative_option_plus  = 80
          derivative_option_minus = 90
          
       case (7)
          if (masterProc .and. level==1) then
             print *,"Preconditioner d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   1 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 100
          derivative_option_minus = 110
          
       case (8)
          if (masterProc .and. level==1) then
             print *,"Preconditioner d/dzeta derivatives discretized using upwinded finite differences:"
             print *,"   2 point on one side, 3 points on the other side."
          end if
          derivative_option_plus  = 120
          derivative_option_minus = 130
          
       case default
          print *,"Error! Invalid setting for zeta_derivative_option:",zeta_derivative_option
          stop
       end select
       
       if (call_uniform_diff_matrices) then
          quadrature_option = 0
          call uniformDiffMatrices(Nzeta, zero, zetaMax, &
               derivative_option_plus,  quadrature_option, zeta, zetaWeights, ddzeta_plus_preconditioner, d2dzeta2)
          call uniformDiffMatrices(Nzeta, zero, zetaMax, &
               derivative_option_minus, quadrature_option, zeta, zetaWeights, ddzeta_minus_preconditioner, d2dzeta2)
       end if

       if (zeta_diffusion>0) then
          if (masterProc .and. level==1) then
             print *,"  Adding diffusion in zeta:",zeta_diffusion
          end if
          derivative_option = 0
          quadrature_option = 0
          allocate(temp_matrix(Nzeta,Nzeta))
          call uniformDiffMatrices(Nzeta, zero, zetaMax, &
               derivative_option,  quadrature_option, zeta, zetaWeights, temp_matrix, d2dzeta2)
          ddzeta_plus_preconditioner  = ddzeta_plus_preconditioner  - zeta_diffusion * d2dzeta2
          ddzeta_minus_preconditioner = ddzeta_minus_preconditioner - zeta_diffusion * d2dzeta2
          deallocate(temp_matrix)
       end if

       if (preconditioner_zeta_derivative_option<0) then
          if (masterProc .and. level==1) then
             print *,"   But only the diagonal is kept."
          end if
          do j=1,Nzeta
             do k=1,Nzeta
                if (j .ne. k) then
                   ddzeta_plus_preconditioner(j,k) = 0
                   ddzeta_minus_preconditioner(j,k) = 0
                end if
             end do
          end do
       end if
       
       levels(level)%ddzeta_sum_preconditioner        = (ddzeta_plus_preconditioner + ddzeta_minus_preconditioner)/two
       levels(level)%ddzeta_difference_preconditioner = (ddzeta_plus_preconditioner - ddzeta_minus_preconditioner)/two

    end if

    zetaWeights = zetaWeights * NPeriods

    ! The following arrays will not be needed:
    deallocate(d2dzeta2)



  end subroutine create_grids


