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
    PetscScalar, dimension(:,:), allocatable :: d2dy2, ddy, ddy_plus, ddy_minus
    PetscScalar, dimension(:), allocatable :: y_dummy, yWeights_dummy, yWeights, dxi_dy, d2xi_dy2
    PetscInt :: scheme, ixi, quadrature_option, derivative_option_plus, derivative_option_minus, derivative_option
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
    PetscScalar, dimension(:), pointer :: theta, zeta, xi, thetaWeights, zetaWeights, xiWeights, y
    PetscScalar, dimension(:,:), pointer :: B, dBdtheta, dBdzeta
    PetscScalar, dimension(:,:), pointer :: ddtheta_plus, ddtheta_minus, ddtheta_plus_preconditioner, ddtheta_minus_preconditioner
    PetscScalar, dimension(:,:), pointer :: ddzeta_plus, ddzeta_minus, ddzeta_plus_preconditioner, ddzeta_minus_preconditioner
    PetscScalar, dimension(:,:), pointer :: ddxi_plus, ddxi_minus, ddxi_plus_preconditioner, ddxi_minus_preconditioner
    PetscScalar, dimension(:,:), pointer :: pitch_angle_scattering_operator, pitch_angle_scattering_operator_preconditioner
    
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
    allocate(d2dtheta2(Ntheta,Ntheta))

    theta => levels(level)%theta
    thetaWeights => levels(level)%thetaWeights
    ddtheta_plus => levels(level)%ddtheta_plus
    ddtheta_minus => levels(level)%ddtheta_minus
    ddtheta_plus_preconditioner => levels(level)%ddtheta_plus_preconditioner
    ddtheta_minus_preconditioner => levels(level)%ddtheta_minus_preconditioner

    ! *******************************************************************************
    ! Handle d/dtheta for the ExB & magnetic drift terms in the main matrix.
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

    ! *******************************************************************************
    ! Handle d/dtheta for the ExB & magnetic drift terms in the preconditioner matrix.
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
       
    end if

    zetaWeights = zetaWeights * NPeriods

    ! The following arrays will not be needed:
    deallocate(d2dzeta2)

    ! *******************************************************************************
    ! *******************************************************************************
    !
    ! Build xi grids, integration weights, and differentiation matrices.
    !
    ! *******************************************************************************
    ! *******************************************************************************

    allocate(levels(level)%y(Nxi))
    allocate(levels(level)%xi(Nxi))
    allocate(levels(level)%xiWeights(Nxi))

    allocate(levels(level)%ddxi_plus(Nxi,Nxi))
    allocate(levels(level)%ddxi_minus(Nxi,Nxi))
    allocate(levels(level)%ddxi_plus_preconditioner(Nxi,Nxi))
    allocate(levels(level)%ddxi_minus_preconditioner(Nxi,Nxi))
    allocate(levels(level)%pitch_angle_scattering_operator(Nxi,Nxi))
    allocate(levels(level)%pitch_angle_scattering_operator_preconditioner(Nxi,Nxi))

    allocate(dxi_dy(Nxi))
    allocate(d2xi_dy2(Nxi))
    allocate(y_dummy(Nxi))
    allocate(yWeights(Nxi))
    allocate(yWeights_dummy(Nxi))
    allocate(ddy(Nxi,Nxi))
    allocate(d2dy2(Nxi,Nxi))
    allocate(ddy_plus(Nxi,Nxi))
    allocate(ddy_minus(Nxi,Nxi))

    y => levels(level)%y
    xi => levels(level)%xi
    xiWeights    => levels(level)%xiWeights

    ddxi_plus => levels(level)%ddxi_plus
    ddxi_minus => levels(level)%ddxi_minus
    ddxi_plus_preconditioner => levels(level)%ddxi_plus_preconditioner
    ddxi_minus_preconditioner => levels(level)%ddxi_minus_preconditioner

    pitch_angle_scattering_operator => levels(level)%pitch_angle_scattering_operator
    pitch_angle_scattering_operator_preconditioner => levels(level)%pitch_angle_scattering_operator_preconditioner

    ! *******************************************************************************
    ! Set up uniform y grid and the associated nonuniform xi grid
    ! *******************************************************************************

    derivative_option = 2
    call uniformDiffMatrices(Nxi, -one, one, derivative_option, xi_quadrature_option, y, yWeights, ddy, d2dy2)
    do j=1,Nxi
       xi(j) = compute_xi_from_y(y(j))
       dxi_dy(j) = compute_dxi_dy(y(j))
       d2xi_dy2(j) = compute_d2xi_dy2(y(j))
    end do
    xiWeights = dxi_dy * yWeights
    print *,"sum(xiWeights):",sum(xiWeights)
    print *,"sum(xiWeights*xi):",sum(xiWeights*xi)
    print *,"sum(xiWeights*xi*xi):",sum(xiWeights*xi*xi)
    if (abs(sum(xiWeights)-2) > 1.0e-2) then
       print *,"Error! xiWeights do not sum to 2!"
       print *,"xiWeights:",xiWeights
       print *,"Sum is",sum(xiWeights)
       stop
    end if

    ! *******************************************************************************
    ! Handle d/dxi for the pitch angle scattering operator in the main matrix.
    ! *******************************************************************************

    select case (pitch_angle_scattering_option)

    case (2)
       if (masterProc .and. level==1) then
          print *,"Pitch angle scattering operator discretized using centered finite differences:"
          print *,"   1 point on either side."
       end if
       derivative_option = 2

    case (3)
       if (masterProc .and. level==1) then
          print *,"Pitch angle scattering operator discretized using centered finite differences:"
          print *,"   2 points on either side."
       end if
       derivative_option = 12

    case default
       print *,"Error! Invalid setting for pitch_angle_scattering_option:",pitch_angle_scattering_option
       stop
    end select

    call uniformDiffMatrices(Nxi, -one, one, derivative_option, xi_quadrature_option, y_dummy, yWeights_dummy, ddy, d2dy2)
    do j=1,Nxi
       !pitch_angle_scattering_operator(j,:) = (1/two)*(1-xi(j)*xi(j))*d2dxi2(j,:) - xi(j)*ddxi(j,:)
       pitch_angle_scattering_operator(j,:) = (1/two)*(1-xi(j)*xi(j))/(dxi_dy(j)*dxi_dy(j))&
            *(d2dy2(j,:) - d2xi_dy2(j)/dxi_dy(j)*ddy(j,:)) &
            - xi(j)/dxi_dy(j)*ddy(j,:)
    end do

    ! *******************************************************************************
    ! Handle d/dxi for the pitch angle scattering operator in the preconditioner matrix.
    ! *******************************************************************************

    call_uniform_diff_matrices = .true.
    select case (abs(preconditioner_pitch_angle_scattering_option))
    case (0)
       if (masterProc .and. level==1) then
          print *,"Pitch angle scattering operator is dropped in the preconditioner."
       end if
       ddy = 0
       d2dy2 = 0
       call_uniform_diff_matrices = .false.

    case (100)
       if (masterProc .and. level==1) then
          print *,"Pitch angle scattering operator is the same in the preconditioner."
       end if
       ! ddy and d2dy2 will be carried over from the previous section then.
       call_uniform_diff_matrices = .false.

    case (2)
       if (masterProc .and. level==1) then
          print *,"Preconditioner pitch angle scattering operator is discretized using centered finite differences:"
          print *,"   1 point on either side."
       end if
       derivative_option = 2

    case (3)
       if (masterProc .and. level==1) then
          print *,"Preconditioner pitch angle scattering operator is discretized using centered finite differences:"
          print *,"   2 points on either side."
       end if
       derivative_option = 12

    case default
       print *,"Error! Invalid setting for preconditioner_pitch_angle_scattering_option:",preconditioner_pitch_angle_scattering_option
       stop
    end select

    if (call_uniform_diff_matrices) then
       call uniformDiffMatrices(Nxi, -one, one, derivative_option, xi_quadrature_option, y_dummy, yWeights_dummy, ddy, d2dy2)
    end if
    do j=1,Nxi
       !pitch_angle_scattering_operator_preconditioner(j,:) = (1/two)*(1-xi(j)*xi(j))*d2dxi2(j,:) - xi(j)*ddxi(j,:)
       pitch_angle_scattering_operator_preconditioner(j,:) = (1/two)*(1-xi(j)*xi(j))/(dxi_dy(j)*dxi_dy(j))&
            *(d2dy2(j,:) - d2xi_dy2(j)/dxi_dy(j)*ddy(j,:)) &
            - xi(j)/dxi_dy(j)*ddy(j,:)
    end do

    if (preconditioner_pitch_angle_scattering_option<0) then
       if (masterProc .and. level==1) then
          print *,"   But only the diagonal is kept."
       end if
       do j=1,Nxi
          do k=1,Nxi
             if (j .ne. k) then
                pitch_angle_scattering_operator_preconditioner(j,k) = 0
             end if
          end do
       end do
    end if

    ! *******************************************************************************
    ! Handle d/dxi for the mirror term in the main matrix.
    ! *******************************************************************************

    select case (xi_derivative_option)

    case (2)
       if (masterProc .and. level==1) then
          print *,"d/dxi derivatives discretized using centered finite differences:"
          print *,"   1 point on either side."
       end if
       derivative_option_plus = 2
       derivative_option_minus = derivative_option_plus

    case (3)
       if (masterProc .and. level==1) then
          print *,"d/dxi derivatives discretized using centered finite differences:"
          print *,"   2 points on either side."
       end if
       derivative_option_plus  = 12
       derivative_option_minus = 12

    case (4)
       if (masterProc .and. level==1) then
          print *,"d/dxi derivatives discretized using upwinded finite differences:"
          print *,"   0 points on one side, 1 point on the other side."
       end if
       derivative_option_plus  = 32
       derivative_option_minus = 42

    case (5)
       if (masterProc .and. level==1) then
          print *,"d/dxi derivatives discretized using upwinded finite differences:"
          print *,"   0 points on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 52
       derivative_option_minus = 62

    case (6)
       if (masterProc .and. level==1) then
          print *,"d/dxi derivatives discretized using upwinded finite differences:"
          print *,"   1 point on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 82
       derivative_option_minus = 92

    case (7)
       if (masterProc .and. level==1) then
          print *,"d/dxi derivatives discretized using upwinded finite differences:"
          print *,"   1 point on one side, 3 points on the other side."
       end if
       derivative_option_plus  = 102
       derivative_option_minus = 112

    case (8)
       if (masterProc .and. level==1) then
          print *,"d/dxi derivatives discretized using upwinded finite differences:"
          print *,"   2 point on one side, 3 points on the other side."
          print *,"   High accuracy at the domain ends, though upwinding breaks down there."
       end if
       derivative_option_plus  = 122
       derivative_option_minus = 132

    case (9)
       if (masterProc .and. level==1) then
          print *,"d/dxi derivatives discretized using upwinded finite differences:"
          print *,"   2 point on one side, 3 points on the other side."
          print *,"   Upwinding all the way to the domain ends, meaning lower accuracy there."
       end if
       derivative_option_plus  = 123
       derivative_option_minus = 133

    case default
       print *,"Error! Invalid setting for xi_derivative_option:",xi_derivative_option
       stop
    end select

    call uniformDiffMatrices(Nxi, -one, one, derivative_option_plus,  xi_quadrature_option, y_dummy, yWeights_dummy, ddy_plus,  d2dy2)
    call uniformDiffMatrices(Nxi, -one, one, derivative_option_minus, xi_quadrature_option, y_dummy, yWeights_dummy, ddy_minus, d2dy2)
    do j=1,Nxi
       ddxi_plus(j,:)  =  ddy_plus(j,:) / dxi_dy(j)
       ddxi_minus(j,:) = ddy_minus(j,:) / dxi_dy(j)
    end do


    ! *******************************************************************************
    ! Handle d/dxi for the mirror term in the preconditioner matrix.
    ! *******************************************************************************

    call_uniform_diff_matrices = .true.
    select case (abs(preconditioner_xi_derivative_option))
    case (0)
       if (masterProc .and. level==1) then
          print *,"d/dxi terms are dropped in the preconditioner."
       end if
       ddxi_plus_preconditioner = 0
       ddxi_minus_preconditioner = 0
       call_uniform_diff_matrices = .false.

    case (100)
       if (masterProc .and. level==1) then
          print *,"d/dxi matrices are the same in the preconditioner."
       end if
       ddxi_plus_preconditioner  = ddxi_plus
       ddxi_minus_preconditioner = ddxi_minus
       call_uniform_diff_matrices = .false.

    case (2)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dxi derivatives discretized using centered finite differences:"
          print *,"   1 point on either side."
       end if
       derivative_option_plus = 2
       derivative_option_minus = derivative_option_plus

    case (3)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dxi derivatives discretized using centered finite differences:"
          print *,"   2 points on either side."
       end if
       derivative_option_plus  = 12
       derivative_option_minus = 12

    case (4)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
          print *,"   0 points on one side, 1 point on the other side."
       end if
       derivative_option_plus  = 32
       derivative_option_minus = 42

    case (5)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
          print *,"   0 points on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 52
       derivative_option_minus = 62

    case (6)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
          print *,"   1 point on one side, 2 points on the other side."
       end if
       derivative_option_plus  = 82
       derivative_option_minus = 92

    case (7)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
          print *,"   1 point on one side, 3 points on the other side."
       end if
       derivative_option_plus  = 102
       derivative_option_minus = 112

    case (8)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
          print *,"   2 point on one side, 3 points on the other side."
          print *,"   High accuracy at the domain ends, though upwinding breaks down there."
       end if
       derivative_option_plus  = 122
       derivative_option_minus = 132

    case (9)
       if (masterProc .and. level==1) then
          print *,"Preconditioner d/dxi derivatives discretized using upwinded finite differences:"
          print *,"   2 point on one side, 3 points on the other side."
          print *,"   Upwinding all the way to the domain ends, meaning lower accuracy there."
       end if
       derivative_option_plus  = 123
       derivative_option_minus = 133

    case default
       print *,"Error! Invalid setting for xi_derivative_option:",xi_derivative_option
       stop
    end select

    if (call_uniform_diff_matrices) then
       call uniformDiffMatrices(Nxi, -one, one, derivative_option_plus,  xi_quadrature_option, y_dummy, yWeights_dummy, ddy_plus,  d2dy2)
       call uniformDiffMatrices(Nxi, -one, one, derivative_option_minus, xi_quadrature_option, y_dummy, yWeights_dummy, ddy_minus, d2dy2)
       do j=1,Nxi
          ddxi_plus_preconditioner( j,:) = ddy_plus( j,:) / dxi_dy(j)
          ddxi_minus_preconditioner(j,:) = ddy_minus(j,:) / dxi_dy(j)
       end do
    end if

    if (preconditioner_xi_derivative_option<0) then
       if (masterProc .and. level==1) then
          print *,"   But only the diagonal is kept."
       end if
       do j=1,Nxi
          do k=1,Nxi
             if (j .ne. k) then
                ddxi_plus_preconditioner(j,k) = 0
                ddxi_minus_preconditioner(j,k) = 0
             end if
          end do
       end do
    end if

!!$    print *,"Here comes ddxi_plus * xi:"
!!$    print *,matmul(ddxi_plus,xi)
!!$    print *,"Here comes ddxi_minus * xi:"
!!$    print *,matmul(ddxi_minus,xi)
!!$    print *,"Here comes (PAS * xi) + xi:"
!!$    print *,matmul(pitch_angle_scattering_operator,xi)+xi

    ! The following arrays will not be needed:
    deallocate(d2dy2,ddy,y_dummy,yWeights_dummy,yWeights,ddy_plus,ddy_minus)


    contains

      function compute_xi_from_y(yy)
        implicit none
        PetscScalar, intent(in) :: yy
        PetscScalar :: compute_xi_from_y

        compute_xi_from_y = nonuniform_xi_a * yy + nonuniform_xi_b * (yy ** 5)
      end function compute_xi_from_y

      function compute_dxi_dy(yy)
        implicit none
        PetscScalar, intent(in) :: yy
        PetscScalar :: compute_dxi_dy

        compute_dxi_dy = nonuniform_xi_a + 5*nonuniform_xi_b * (yy ** 4)
      end function compute_dxi_dy

      function compute_d2xi_dy2(yy)
        implicit none
        PetscScalar, intent(in) :: yy
        PetscScalar :: compute_d2xi_dy2

        compute_d2xi_dy2 = 20*nonuniform_xi_b * (yy ** 3)
      end function compute_d2xi_dy2

  end subroutine create_grids


