!<compile=optimized>

#include "../include/dprec.fh"

!> Generic object interface for various MDIIS implementations.  A
!! specific implementation is chosen by the user at run time via the
!! constructor.  All subroutines and functions should then be called
!! in the same manner with only details of execution being different.
!!
!! MDIIS uses solutions and residual (whatever needs to be minimized
!! for each data point) from previous iterations to predict a new
!! solution with a lower residual.  For all methods, working memory
!! for the data and the residual must be allocated and passed to the
!! constructor.  This memory must not be reallocated/deallocated
!! before the MDIIS object is destroyed.  Data for solutions and
!! residuals are vectors and the memory for a data vector may be of
!! any rank (in practise > 2).  The total memory size should be the
!! vector length * the number of past solutions you wish the MDIIS
!! routine to use.
!!
!! Depending on the implementation chosen, the active vector to
!! read/write data from may change location in working memory after
!! calling MDIIS_ADVANCE(). Use MDIIS_GETVECTORINDEX() to get the
!! active vector number.
!!
!! To add a new MDIIS implementation a new class must be written and
!! registered here and in the Makefile.  Within this class, the new
!! implementation must be registered with a 'USE' statement, added as
!! a type (with enumeration constant) and appropriate calls made from
!! each subroutine/function.
module mdiis_c
  use mdiis_orig_c
  use mdiis_blas_c
  use mdiis_blas2_c
  implicit none

  type mdiis
     integer :: KOVALENKO = 0, KOVALENKO_OPT = 1, KOVALENKO_OPT2 = 2
     type(mdiis_orig),pointer :: orig => null()
     type(mdiis_blas),pointer :: blas => null()
     type(mdiis_blas2),pointer :: blas2 => null()
  end type mdiis

contains

  !> Create new MDIIS object.    Provide convergence parameters and
  !! working memory.  Note that the working memory must me nVec times
  !! larger than the data array.  Further, the number of data points
  !! cannot change between iterations.  xi(:,1) and ri(:,1) are the active vectors.
  !! @param[in,out] this MDIIS object.
  !! @param[in] delta Coefficient for residual gradient.  Should be between 0 and 1.
  !! @param[in] tol Target residual.
  !! @param[in] xi Array of vector data. np X nVec.
  !! @param[in] ri Array of residual data. np X nVec.
  !! @param[in] np Number of data points (first dimension length).
  !! @param[in] nVec Length of DIIS vectors.
  !! @param[in] restart Restart threshold factor. Ratio of the current residual to the
  !!     minimum residual in the basis that causes a restart.
  subroutine mdiis_new (this, method, delta, tol, restart)
    implicit none
    type(mdiis), intent(inout) :: this
    integer, intent(in) :: method
    _REAL_, intent(in) :: delta, tol, restart

    _REAL_ :: rms
    logical :: con

    if (method == this%KOVALENKO) then
       allocate(this%orig)
       call mdiis_orig_new_serial(this%orig, delta, tol, restart)
    else if (method == this%KOVALENKO_OPT) then
       allocate(this%blas)
       call mdiis_blas_new_serial(this%blas, delta, tol, restart)
    else if (method == this%KOVALENKO_OPT2) then
       allocate(this%blas2)
       call mdiis_blas2_new_serial(this%blas2, delta, tol, restart)
    end if
  end subroutine mdiis_new

  
  !> Create new MDIIS object.  Provide convergence parameters, working
  !! memory and MPI parameters.  Note that the working memory must me
  !! nVec times larger than the data array.  Further, the number of
  !! data points cannot change between iterations.  xi(:,1) and
  !! ri(:,1) are the active vectors.
  !! @param[out] this MDIIS object.
  !! @param[in] delta Coefficient for residual gradient.  Should be between 0 and 1.
  !! @param[in] tol Target residual.
  !! @param[in] xi Array of vector data.  np X nVec.
  !! @param[in] ri Array of residual data.  np X nVec.
  !! @param[in] np Number of data points (first dimension length).
  !! @param[in] nVec Length of DIIS vectors.
  !! @param[in] restart Restart threshold factor. Ratio of the current residual to the
  !!     minimum residual in the basis that causes a restart.
  !! @param[in] rank Mpi process rank.
  !! @param[in] size Number of processes.
  !! @param[in] comm Mpi communicator.
  subroutine mdiis_new_mpi (this, method, delta, tol, restart, &
       rank, size, comm)
    implicit none
    type(mdiis), intent(out) :: this
    integer, intent(in) :: method
    _REAL_, intent(in) :: delta, tol, restart
    integer, intent(in) :: rank, size, comm

    if (method == this%KOVALENKO) then
       allocate(this%orig)
       call mdiis_orig_new_mpi(this%orig, delta, tol, restart, &
            rank, size, comm)
    else if (method == this%KOVALENKO_OPT) then
       allocate(this%blas)
       call mdiis_blas_new_mpi(this%blas, delta, tol, restart, &
            rank, size, comm)
    else if (method == this%KOVALENKO_OPT2) then
       allocate(this%blas2)
       call mdiis_blas2_new_mpi(this%blas2, delta, tol, restart, &
            rank, size, comm)
    end if
  end subroutine mdiis_new_mpi


  !! Destroy MDIIS object. All internal variables and pointers are
  !! reset. Working memory remains intact memory is swapped (if
  !! necessary) such that xi(:,1) and ri(:,1) become the active vectors.
  !! IN:
  !!    this :: mdiis object
  subroutine mdiis_destroy(this)
    implicit none
    type(mdiis), intent(inout) :: this
    if (associated(this%orig)) then
       call mdiis_orig_destroy(this%orig)
       deallocate(this%orig)
    end if
    if (associated(this%blas)) then
       call mdiis_blas_destroy(this%blas)
       deallocate(this%blas)
    end if
    if (associated(this%blas2)) then
       call mdiis_blas2_destroy(this%blas2)
       deallocate(this%blas2)
    end if
  end subroutine mdiis_destroy


  !! Reset MDIIS object. Progress is reset to a single basis
  !! vector.Working memory remains intact memory is swapped (if
  !! necessary) such that xi(:,1) and ri(:,1) become the active vectors.
  !! All other variables are untouched.
  !! IN:
  !!    this :: mdiis object
  subroutine mdiis_reset(this)
    implicit none
    type(mdiis), intent(inout) :: this
    if (associated(this%orig)) then
       call mdiis_orig_reset(this%orig)
    end if
    if (associated(this%blas)) then
       call mdiis_blas_reset(this%blas)
    end if
    if (associated(this%blas2)) then
       call mdiis_blas2_reset(this%blas2)
    end if
  end subroutine mdiis_reset


  !! Sets data vectors and resets MDIIS object. Progress is reset to a single basis
  !! vector. Working memory remains intact memory is swapped (if
  !! necessary) such that xi(:,1) and ri(:,1) become the active vectors.
  !! All other variables are untouched.
  !! IN:
  !!    this :: mdiis object
  !!    xi    :: array of vector data.  np X nvec
  !!    ri    :: array of residual data.  np X nvec
  !!    np    :: number of data points (first dimension length)
  !!    nvec   :: length of DIIS vectors
  subroutine mdiis_setData(this, xi, ri, np, nvec)
    implicit none
    type(mdiis), intent(inout) :: this
    _REAL_, target, intent(in) :: xi(np, nvec), ri(np, nvec)
    integer, intent(in) ::  np
    integer, intent(in) :: nvec
    call mdiis_resize(this, xi, ri, np, nvec)
  end subroutine mdiis_setData


  !! Resizes and resets MDIIS object. Progress is reset to a single basis
  !! vector. Working memory remains intact memory is swapped (if
  !! necessary) such that xi(:,1) and ri(:,1) become the active vectors.
  !! All other variables are untouched.
  !! IN:
  !!    this :: mdiis object
  !!    xi    :: array of vector data.  np X nvec
  !!    ri    :: array of residual data.  np X nvec
  !!    np    :: number of data points (first dimension length)
  !!    nvec   :: length of DIIS vectors
  subroutine mdiis_resize(this, xi, ri, np, nvec)
    implicit none
    type(mdiis), intent(inout) :: this
    _REAL_, target, intent(in) :: xi(np, nvec), ri(np, nvec)
    integer, intent(in) ::  np
    integer, intent(in) :: nvec
    if (associated(this%orig)) then
       call mdiis_orig_resize(this%orig, xi, ri, np, nvec)
    end if
    if (associated(this%blas)) then
       call mdiis_blas_resize(this%blas, xi, ri, np, nvec)
    end if
    if (associated(this%blas2)) then
       call mdiis_blas2_resize(this%blas2, xi, ri, np, nvec)
    end if
  end subroutine mdiis_resize

  !! Returns the number of vector currently used as a basis by the MDIIS object
  !! IN:
  !!    this :: mdiis object
  !! OUT:
  !!    The number of vectors currently used
  function getCurrentNVec (this) result(nVec)
    implicit none
    type(mdiis), intent(inout) :: this
    integer :: nVec
    if (associated(this%orig)) then
       nVec = mdiis_orig_getCurrentNVec(this%orig)
    else if (associated(this%blas)) then
       nVec = mdiis_blas_getCurrentNVec(this%blas)
    else if (associated(this%blas2)) then
       nVec = mdiis_blas2_getCurrentNVec(this%blas2)
    end if
  end function getCurrentNVec


  !! Returns the index number of the working vector
  !! IN:
  !!    this :: mdiis object
  !! OUT:
  !!     the number of the vector to read/write
  function mdiis_getWorkVector (this) result(index)
    implicit none
    type(mdiis), intent(inout) :: this
    integer :: index
    if (associated(this%orig)) then
       index = mdiis_orig_getWorkVector(this%orig)
    else if (associated(this%blas)) then
       index = mdiis_blas_getWorkVector(this%blas)
    else if (associated(this%blas2)) then
       index = mdiis_blas2_getWorkVector(this%blas2)
    end if
  end function mdiis_getWorkVector


  !! One MDIIS iteration. Query mdiis_getVectorNumer() to get the new active vector.
  !! IN:
  !!    this  :: mdiis object
  !!    rms1  :: residual at the end of the step
  !!    conver:: have we converged
  !!    tolerance :: tolerance for this calculation
  subroutine  mdiis_advance (this, rms1, conver, tolerance_o)
    implicit none
    type(mdiis), intent(inout) :: this
    _REAL_, intent(out) ::  rms1
    logical, intent(inout) ::  conver
    _REAL_, optional, intent(in) ::  tolerance_o
#ifdef RISM_DEBUG
    write(6, * )"MDIIS_advance generic", rms1, conver
#endif

    if (associated(this%orig)) then
       if (present(tolerance_o)) then
          call mdiis_orig_advance(this%orig, rms1, conver, tolerance_o)
       else
          call mdiis_orig_advance(this%orig, rms1, conver)
       end if
    else if (associated(this%blas)) then
       if (present(tolerance_o)) then
          call mdiis_blas_advance(this%blas, rms1, conver, tolerance_o)
       else
          call mdiis_blas_advance(this%blas, rms1, conver)
       end if
    else if (associated(this%blas2)) then
       if (present(tolerance_o)) then
          call mdiis_blas2_advance(this%blas2, rms1, conver, tolerance_o)
       else
          call mdiis_blas2_advance(this%blas2, rms1, conver)
       end if
    end if
  end subroutine mdiis_advance

end module mdiis_c
