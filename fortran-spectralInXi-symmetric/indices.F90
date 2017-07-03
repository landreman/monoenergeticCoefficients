module indices

  implicit none

contains

  integer function getIndex(i_theta, i_zeta, i_xi)

    use variables, only: Nxi, Ntheta, Nzeta, matrixSize

    implicit none

    integer, intent(in) :: i_xi, i_theta, i_zeta

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
