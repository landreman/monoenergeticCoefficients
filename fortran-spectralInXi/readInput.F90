subroutine readInput()

  use variables

  implicit none

  character(len=14) :: filename="input.namelist"
  integer :: fileUnit, didFileAccessWork

  namelist / inputs / Ntheta, Nzeta, Nxi, nu

  fileUnit=11
  open(unit=fileUnit, file=filename,    action="read", status="old", iostat=didFileAccessWork)
  if (didFileAccessWork /= 0) then
     print *,"Proc ",myRank,": Error opening ", trim(filename)
     stop
  else
     read(fileUnit, nml=inputs, iostat=didFileAccessWork)
     if (didFileAccessWork /= 0) then
        print *,"Proc ",myRank,": Error!  I was able to open the file ", trim(filename), &
             " but not read data from the inputs namelist in it."
        stop
     end if
     if (masterProc) then
        print *,"Successfully read parameters from inputs namelist in ", trim(filename), "."
     end if
  end if

  close(unit = fileUnit)

end subroutine readInput

