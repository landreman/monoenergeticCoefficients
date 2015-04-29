! NOTE: Unlike getIndex in sfincs, the function getBlockIndex in this
! file does not return a 0-based index for PETSC. Instead it returns a
! fortran 1-based index!

module indices

  implicit none

contains

  integer function getBlockIndex(i_theta, i_zeta)

    use variables, only: Ntheta, Nzeta, mblock

    implicit none

    integer, intent(in) :: i_theta, i_zeta

    ! Validate inputs:

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

    getBlockIndex = (i_theta-1)*Nzeta + i_zeta

    if (getBlockIndex < 0) then
       print *,"Error! Something went wrong, and the index came out less than 0."
       stop
    end if

    if (getBlockIndex >= mblock) then
       print *,"Error! Something went wrong, and the index came out larger than the matrix size."
       stop
    end if

  end function getBlockIndex

end module indices
