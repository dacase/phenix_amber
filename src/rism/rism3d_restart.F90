!<compile=optimized>

#include "../include/dprec.fh"

!> Reads/writes non-portable binary restart files for 3D-RISM
!! calculations. The file format is Fortran unformatted binary and has
!! the form:
!! 
!! numGridPoints numSolventAtomTypes
!! correlationFunction(numGridPoints,numSolventAtomTypes)
!! 
!! The number of grid points and number of solvent types, giving the
!! size of the correlation function array.
!!
!! This format is not portable and not MPI friendly in any way and is
!! deprecated.
module rism3d_restart
  use rism_report_c
  implicit none

contains

  !> Read restart file storing previous correlation function solutions
  !! in an unformatted array.
  !! @param[in] file File name.
  !! @param[out] cuv Correlation function.
  !! @param[in] numGridPoints Expected number of data points for each solvent type.
  !! @param[in] numSolventAtomTypes Expected number of solvent types.
  subroutine readRestartFile(file, cuv, numGridPoints, numSolventAtomTypes)
    implicit none
    character(*), intent(in) :: file
    integer, intent(in) :: numGridPoints, numSolventAtomTypes
    _REAL_, intent(out) :: cuv(numGridPoints, numSolventAtomTypes)
    integer :: ncha = 99, iostat, i, j
    integer :: numGridPointso, numSolventAtomTypeso
    open (ncha, file = file, form = 'unformatted', &
         status = 'old', iostat = iostat, position = "REWIND")
    if (iostat /= 0) then
       call rism_report_error("on "//file//" open")
    end if
    read (ncha, iostat = iostat)  numGridPointso, numSolventAtomTypeso
    if (iostat /= 0) then
       call rism_report_error("on "//file//" read")
    end if
    if (numGridPointso /= numGridPoints .OR. numSolventAtomTypeso /= numSolventAtomTypes)  then
       close (ncha, iostat = iostat)
       call rism_report_error( &
            '(a, i6, a, i3, a, i6, a, i3)', 'readRestartFile:  actual   NumGridPoints = ', &
            numGridPoints, ', NumberSolventAtomTypes = ', numSolventAtomTypes, &
            '; saved   NumGridPoints = ', numGridPointso, ', NumberSolventAtomTypes = ', numSolventAtomTypeso)
    end if
    do i = 1, numSolventAtomTypes
       read (ncha, iostat = iostat)  cuv(1:numGridPoints, i)
       if (iostat /= 0) then
          call rism_report_error('(a, i4)', 'readRestartFile: I / O error:', iostat)
       end if
    end do
    close (ncha, iostat = iostat)
    if (iostat /= 0) then
       call rism_report_error("on "//file//" close")
    end if
  end subroutine readRestartFile


  !> Write restart file storing previous correlation function
  !! solutions in an unformatted array.
  !! @param[in] file File name.
  !! @param[in] cuv Correlation function.
  !! @param[in] numGridPoints Number of data points for each solvent type.
  !! @param[in] numSolventAtomTypes Number of solvent types.
  subroutine writeRestartFile(file, cuv, numGridPoints, numSolventAtomTypes)
    implicit none
    character( * ), intent(in) :: file
    integer, intent(in) :: numGridPoints, numSolventAtomTypes
    _REAL_, intent(in) :: cuv(numGridPoints, numSolventAtomTypes)
    integer :: ncha = 99, iostat, i
    open (ncha, file = file, form = 'unformatted', iostat = iostat, status = "REPLACE")
    if (iostat /= 0) then
       call rism_report_error("on "//file//" open")
    end if
    write (ncha, iostat = iostat)  numGridPoints, numSolventAtomTypes
    if (iostat /= 0) then
       call rism_report_error("on "//file//" write")
    end if
    do i = 1, numSolventAtomTypes
       write (ncha, iostat = iostat)  cuv(1:numGridPoints, i)
       if (iostat /= 0) then
          call rism_report_error('(a, i4)', 'writeRestartFile:  I / O error:', iostat)
       end if
    end do
    close (ncha, iostat = iostat)
    if (iostat /= 0) then
       call rism_report_error("on "//file//" close")
    end if
  end subroutine writeRestartFile
  
end module rism3d_restart
