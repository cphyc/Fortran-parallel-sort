program mod_sort
  use omp_lib
  use m_mrgrnk

  integer, parameter :: N = int(1e8)
  integer :: A(N), order(N), order2(N)
  integer, parameter :: max_simple_sort_size = 20
  real :: before, after

  do i = 1, N
     A(i) = N - 2*i
  end do

  write(*, *) 'Before'
  if (N < 20) write(*, '(*(i5))') A

  before = omp_get_wtime()
  call parallel_sort(A, order)
  write(*, *) 'parallel sort:', omp_get_wtime() - before

  before = omp_get_wtime()
  call mrgrnk(A, order2)
  write(*, *) 'mrgrnk:       ', omp_get_wtime() - before

  write(*, *) 'After', all(A(order) == A(order2))
  if (N < 20) write(*, '(*(i5))') A(order)

contains

  subroutine parallel_sort (A, order)
    integer, intent(in),  dimension(:) :: A
    integer, intent(out), dimension(size(A)) :: order

    integer :: len, from, to, nthreads, thread, chunk, i

    len      = size(A)
    nthreads = omp_get_max_threads()
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
    !$OMP shared(A, order) private(from, to) &
    !$OMP shared(chunk, len, nthreads)
    do thread = 0, nthreads
       from = thread           * chunk + 1
       to   = min((thread + 1) * chunk, len)

       call mrgrnk(A(from:to), order(from:to))
       order(from:to) = order(from:to) + from - 1

    end do
    !$OMP end parallel do

    !----------------------------------------
    ! Merge pieces together
    !----------------------------------------
    i = 1
    chunk2 = chunk
    do while (chunk2 < size(A))

       !$OMP parallel do default(none)  &
       !$OMP shared(chunk2, A, order)   &
       !$OMP private(from, middle, to)
       do thread = 0, ceiling(.5 * size(A) / chunk2)
          from   = thread*2     * chunk2 + 1
          middle = (thread*2+1) * chunk2
          to     = (thread*2+2) * chunk2

          middle = min(middle, size(A))
          to     = min(to, size(A))
          if (from < to) then
             call merge(A, order, from, middle, to)
          end if
       end do
       !$OMP end parallel do
       chunk2 = chunk2 * 2
       i = i + 1
    end do
  end subroutine parallel_sort

  recursive subroutine sort (A, order, left, right)
    integer, intent(in), dimension(:)          :: A
    integer, intent(inout), dimension(size(A)) :: order

    integer, intent(in) :: left, right

    integer :: tmp

    if (left < right + max_simple_sort_size) then
       call interchange_sort(A, order, left, right)
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
       if (left < j)  call sort(A, order, left, j)
       if (i < right) call sort(A, order, i, right)
    end if

  end subroutine sort

  subroutine interchange_sort (A, order, left, right)
    integer, intent(in), dimension(:) :: A
    integer, intent(inout), dimension(size(A)) :: order

    integer, intent(in) :: left, right

    integer :: i, j, tmp

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

  !> Merge two parts of A, ordered by order from left to right
  !! around middle.
  subroutine merge (A, order, left, middle, right)
    integer, intent(in), dimension(:) :: A
    integer, intent(inout), dimension(size(A)) :: order

    integer, intent(in) :: left, middle, right

    integer :: leftA, rightA, leftB, rightB
    integer :: iA, iB
    integer :: lenA, lenB

    integer, dimension(left    :middle) :: orderA
    integer, dimension(middle+1:right ) :: orderB

    ! copy order
    orderA = order(left    :middle)
    orderB = order(middle+1:right)

    ! more explicit variables
    leftA  = left
    rightA = middle
    leftB  = middle+1
    rightB = right

    ! initialize iA, iB to their leftmost position
    iA = leftA
    iB = leftB

    i = leftA

    do while ((iA <= rightA) .and. (iB <= rightB))
       if (A(orderA(iA)) <= A(orderB(iB))) then
          order(i) = orderA(iA)
          iA = iA + 1
       else
          order(i) = orderB(iB)
          iB = iB + 1
       end if

       i = i + 1
    end do

    ! either A or B still have elements, append them to the new order
    do while (iA <= rightA)
       order(i) = orderA(iA)
       iA = iA + 1

       i  = i + 1

    end do
    do while (iB <= rightB)
       order(i) = orderB(iB)
       iB = iB + 1

       i  = i + 1
    end do

  end subroutine merge

  subroutine quick_sort(list, order)
    ! quick sort routine from:
    ! brainerd, w.s., goldberg, c.h. & adams, j.c. (1990) "programmer's guide to
    ! fortran 90", mcgraw-hill  isbn 0-07-000248-7, pages 149-150.
    ! modified by alan miller to include an associated integer array which gives
    ! the positions of the elements in the original order.
    integer, parameter :: i8b = 8

    integer, dimension (:), intent(inout)        :: list
    integer, dimension (size(list)), intent(out) :: order

    integer :: n
    ! local variable
    integer :: i

    n = size(list)
    do i = 1, n
       order(i) = i
    end do

    call quick_sort_1(list, order, 1, n)
  end subroutine quick_sort

  recursive subroutine quick_sort_1(list, order, left_end, right_end)
    integer, dimension (:), intent(inout)          :: list
    integer, dimension (size(list)), intent(inout) :: order

    integer, intent(in) :: left_end, right_end

    !     local variables
    integer             :: i, j, itemp
    integer             :: reference, temp

    if (right_end < left_end + max_simple_sort_size) then
       ! use interchange sort for small lists
       call interchange_sort_1(list, order, left_end, right_end)

    else
       ! use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1

       do
          ! scan list from left end until element >= reference is found
          do
             i = i + 1
             if (list(order(i)) >= reference) exit
          end do
          ! scan list from right end until element <= reference is found
          do
             j = j - 1
             if (list(order(j)) <= reference) exit
          end do


          if (i < j) then
             ! swap two out-of-order elements
             itemp = order(i); order(i) = order(j); order(j) = itemp
          else if (i == j) then
             i = i + 1
             exit
          else
             exit
          end if
       end do

       if (left_end < j) call quick_sort_1(list, order, left_end, j)
       if (i < right_end) call quick_sort_1(list, order, i, right_end)
    end if

  end subroutine quick_sort_1

  subroutine interchange_sort_1(list, order, left_end, right_end)
    integer, dimension (:), intent(inout)        :: list
    integer, dimension (size(list)), intent(out) :: order

    integer, intent(in) :: left_end, right_end

    !     local variables
    integer :: i, j, itemp
    integer :: temp

    do i = left_end, right_end - 1
       do j = i+1, right_end
          if (list(order(i)) > list(order(j))) then
             itemp = order(i); order(i) = order(j); order(j) = itemp
          end if
       end do
    end do

  end subroutine interchange_sort_1

end program mod_sort
