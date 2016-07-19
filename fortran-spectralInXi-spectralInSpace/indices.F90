module indices

  implicit none

contains

  integer function getIndex(i_Fourier, i_xi)

    use variables, only: NFourier2, Nxi, matrixSize

    implicit none

    integer, intent(in) :: i_xi, i_Fourier

    ! Validate inputs:

    if (i_xi < 1) then
       print *,"Error: i_xi < 1"
       stop
    end if

    if (i_xi > Nxi) then
       print *,"Error: i_xi > Nxi"
       stop
    end if

    if (i_Fourier < 1) then
       print *,"Error: i_Fourier < 1"
       stop
    end if

    if (i_Fourier > NFourier2) then
       print *,"Error: i_Fourier > NFourier2"
       stop
    end if

    getIndex = (i_xi-1)*NFourier2 &
         +i_Fourier -1

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
