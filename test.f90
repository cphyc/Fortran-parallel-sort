!> Test the parallel sort routine
program test
  use ISO_FORTRAN_ENV, only : ERROR_UNIT
  use mod_sort
  use m_mrgrnk

  integer, parameter :: N = int(1e6)

  integer :: A(N), order(N), orderControl(N)

  integer :: i
  logical :: ok

  do i = 1, N
     A(i) = N - 2*i
  end do

  call parallel_sort(A, order)
  call mrgrnk(A, orderControl)

  ok = all(A(order) == A(orderControl))

  if (.not. ok) then
     write(ERROR_UNIT, *) 'A(order) ≠ A(orderControl)'
     stop 1
  else
     stop 0
  end if

end program test
