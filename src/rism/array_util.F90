!<compile=optimized>

!Collection of utility functions for arrays

module array_util
  implicit none
  interface array_index
     module procedure index_c, index_i, index_r8!, index_l
  end interface
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the first element of the array that matches the given value.  If 
!!!'back' then returns the last element.
!!!IN:
!!!   a     : character array
!!!   value : character string
!!!   back  : (optional) search backwards
!!!OUT:
!!!    Index of the element matching the value, -1 if it does not exist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function index_c(a,value,back) result(i)
    implicit none
    character(len=*), intent(in) :: a(:), value
    logical, optional, intent(in) :: back
    integer :: i, start, finish, incr
    start = lbound(a,1)
    finish = ubound(a,1)
    incr = 1
    if(present(back))then
       if(back)then
          start = ubound(a,1)
          finish = lbound(a,1)
          incr = -1
       end if
    end if
    do i=start,finish,incr
       if(trim(a(i)) .eq. trim(value)) return
    end do
    i = -1
  end function index_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the first element of the array that matches the given value.  If 
!!!'back' then returns the last element.
!!!IN:
!!!   a     : integer array
!!!   value : integer value
!!!   back  : (optional) search backwards
!!!OUT:
!!!    Index of the element matching the value, -1 if it does not exist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function index_i(a,value,back) result(i)
    implicit none
    integer, intent(in) :: a(:), value
    logical, optional, intent(in) :: back
    integer :: i, start, finish, incr
    start = lbound(a,1)
    finish = ubound(a,1)
    incr = 1
    if(present(back))then
       if(back)then
          start = ubound(a,1)
          finish = lbound(a,1)
          incr = -1
       end if
    end if
    do i=start,finish,incr
       if(a(i) == value) return
    end do
    i = -1
  end function index_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the first element of the array that matches the given value.  If 
!!!'back' then returns the last element.
!!!IN:
!!!   a     : real*8 array
!!!   value : real*8 value
!!!   back  : (optional) search backwards
!!!OUT:
!!!    Index of the element matching the value, -1 if it does not exist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function index_r8(a,value,back) result(i)
    implicit none
    real*8, intent(in) :: a(:), value
    logical, optional, intent(in) :: back
    integer :: i, start, finish, incr
    start = lbound(a,1)
    finish = ubound(a,1)
    incr = 1
    if(present(back))then
       if(back)then
          start = ubound(a,1)
          finish = lbound(a,1)
          incr = -1
       end if
    end if
    do i=start,finish,incr
       if(a(i) == value) return
    end do
    i = -1
  end function index_r8
end module array_util
