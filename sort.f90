module mod_sort
  use omp_lib
contains

  subroutine parallel_sort (A, order)
    integer, intent(in),  dimension(:) :: A
    integer, intent(out), dimension(size(A)) :: order

    integer :: len, from, to, nthreads, thread, chunk, i

    len      = size(A)
    nthreads = omp_get_num_threads()
    chunk    = len / nthreads

    !----------------------------------------
    ! Initialize order
    !----------------------------------------
    do i = 1, len
       order(i) = i
    end do

    !----------------------------------------
    ! Sort each chunk
    !----------------------------------------
    !$OMP parallel do default(none) &
    !$OMP private(from, to) shared(A, order) &
    !$OMP firstprivate(len, nthreads)
    do thread = 0, nthreads - 1
       from = thread           * chunk + 1
       to   = min((thread + 1) * chunk, len)

       call sort(A(from:to), order(from:to))
    end do
    !$OMP end parallel

    !----------------------------------------
    ! Merge pieces together
    !----------------------------------------
    i = 1
    do while ((1. * nthreads / 2**i) >= 1.)
       !TODO: parallel
       do thread = 0, (nthreads - 1) / 2**i
          from   = thread              * chunk + 1   ! thread
          middle = (thread + 2**(i-1)) * chunk       ! thread + 1*2**i
          to     = min((thread + 2**i) * chunk, len) ! thread + 2*2**i

          if (middle < to) then
             call merge(order(from:middle), order(middle+1:to))
          end if
       end do
       i = i + 1
    end do
  end subroutine parallel_sort

  subroutine sort (A, order)
    integer, intent(in), dimension(:)        :: A
    integer, intent(inout), dimension(size(A)) :: order

    integer :: left, right, tmp

    left  = lbound(A)
    right = ubound(A)

    if (left < right + max_simple_sort_size) then
       call interchange_sort(A, order)
    else
       ref = A((left + right) / 2)
       i = left
       j = right

       do while (i <= j)
          ! find first ≥ than ref
          do while(A(order(i)) >= ref)
             i = i + 1
          end do

          ! find last ≤ ref
          do while(A(order(j)) <= ref)
             j = j - 1
          end do

          ! swap them if required
          if (i < j) then
             tmp = order(i)
             order(i) = order(j)
             order(j) = tmp
          else if (i == j) then
             i = i + 1
          end if
       end do
       ! now i >= j, recursive call with that
       if (left < j)  call sort(A(left:j), order(left:j))
       if (i < right) call sort(A(i:right), order(right:i))
    end if

  end subroutine sort

  subroutine interchange_sort (A, order)
    integer, intent(in), dimension(:) :: A
    integer, intent(inout), dimension(size(A)) :: order

    integer :: left, right
    integer :: i, j, tmp

    left = lbound(A)
    right = ubound(A)

    do i = left, right - 1
       do j = i+1, right
          if (A(order(i)) > A(order(j))) then
             tmp      = order(i)
             order(i) = order(j)
             order(j) = tmp
          end if
       end do
    end do

  end subroutine interchange_sort


end module mod_sort
