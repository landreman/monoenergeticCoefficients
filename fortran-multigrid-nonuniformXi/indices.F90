module indices

  implicit none

contains

  integer function getIndex(level, i_theta, i_zeta, i_xi)

    !use variables, only: Nxi, Ntheta, Nzeta, matrixSize
    use variables, only: levels

    implicit none

    integer, intent(in) :: level, i_xi, i_theta, i_zeta
    integer :: Ntheta, Nzeta, Nxi, matrixSize

    Ntheta = levels(level)%Ntheta
    Nzeta  = levels(level)%Nzeta
    Nxi    = levels(level)%Nxi
    matrixSize = levels(level)%matrixSize

    ! Validate inputs:

    if (i_xi < 1) then
       print *,"Error: i_xi < 1"
       stop
    end if

    if (i_xi > Nxi) then
       print *,"Error: i_xi > Nxi"
       stop
    end if

    if (i_theta < 1) then
       print *,"Error: i_theta < 1"
       stop
    end if

    if (i_theta > Ntheta) then
       print *,"Error: i_theta > Ntheta"
       stop
    end if

    if (i_zeta < 1) then
       print *,"Error: i_zeta < 1"
       stop
    end if

    if (i_zeta > Nzeta) then
       print *,"Error: i_zeta > Nzeta"
       stop
    end if

    getIndex = (i_xi-1)*Ntheta*Nzeta &
         +(i_theta-1)*Nzeta &
         +i_zeta -1

    if (getIndex < 0) then
       print *,"Error! Something went wrong, and the index came out less than 0."
       stop
    end if

    if (getIndex >= matrixSize) then
       print *,"Error! Something went wrong, and the index came out larger than the matrix size."
       stop
    end if

  end function getIndex

end module indices
