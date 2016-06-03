module mod_sort
  use omp_lib
  use m_mrgrnk

  implicit none

  integer, parameter :: max_simple_sort_size = 20
  private
  public :: parallel_sort
contains

  subroutine parallel_sort (A, order)
    integer, intent(in),  dimension(:) :: A
    integer, intent(out), dimension(size(A)) :: order

    integer :: len, from, middle, to, nthreads, thread, chunk, chunk2, i

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
    !$OMP parallel do default(none)          &
    !$OMP firstprivate(chunk, len, nthreads) &
    !$OMP shared(A, order) private(from, to)
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

  !> Merge two parts of A, ordered by order from left to right
  !! around middle.
  subroutine merge (A, order, left, middle, right)
    integer, intent(in), dimension(:) :: A
    integer, intent(inout), dimension(size(A)) :: order

    integer, intent(in) :: left, middle, right

    integer :: leftA, rightA, leftB, rightB
    integer :: iA, iB, i
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

end module mod_sort
