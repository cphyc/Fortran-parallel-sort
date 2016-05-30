program sort_test
  use omp_lib
  use mod_sort
  use m_mrgrnk

  implicit none


  integer, parameter :: N = int(1e8)
  integer :: A(N), order(N), order2(N), i
  real :: before, after

  do i = 1, N
     A(i) = N - 2*i
  end do

  before = omp_get_wtime()
  call parallel_sort(A, order)
  write(*, *) 'parallel sort     :', omp_get_wtime() - before

  before = omp_get_wtime()
  call mrgrnk(A, order2)
  write(*, *) 'mrgrnk (reference):', omp_get_wtime() - before

  write(*, *) 'comparing outputs (should print "T"):', all(A(order) == A(order2))
end program sort_test
