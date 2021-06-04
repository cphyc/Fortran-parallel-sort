module mod_sort
  use omp_lib
  use m_mrgrnk

  implicit none

  integer, parameter :: max_simple_sort_size = 20
  private
  public :: parallel_sort

  interface parallel_sort
     module procedure C_parallel_sort, D_parallel_sort, R_parallel_sort, I_parallel_sort
  end interface parallel_sort
contains

  subroutine C_parallel_sort (A, order)
    character(*), intent(in),  dimension(:) :: A
    integer, intent(out), dimension(size(A)) :: order

    integer :: ilen, from, middle, ito, nthreads, thread, chunk, chunk2, i, iremainder, extraThread 

    ilen      = size(A)
    nthreads = omp_get_max_threads()
    chunk    = ilen / nthreads

    iremainder = mod(ilen, nthreads)
    if (iremainder /= 0) then
        extraThread = 1
    else
        extraThread = 0
    endif   
    !----------------------------------------
    ! Initialize order
    !----------------------------------------
    do i = 1, ilen
       order(i) = i
    end do

    !----------------------------------------
    ! Sort each chunk
    !----------------------------------------
    !$OMP parallel do default(shared) private(thread, from, ito)
    do thread = 0, nthreads-1 + extraThread 
       from = thread*chunk + 1
       ito  = min((thread + 1)*chunk, ilen)

       call mrgrnk(A(from:ito), order(from:ito))
       order(from:ito) = order(from:ito) + from - 1
    end do
    !$OMP end parallel do

    !----------------------------------------
    ! Merge pieces together
    !----------------------------------------
    i = 1
    chunk2 = chunk
    do while (chunk2 < size(A))

       !$OMP parallel do default(shared) private(thread, from, middle, ito)
       do thread = 0, ceiling(.5 * size(A) / chunk2)
          from   = thread*2     * chunk2 + 1
          middle = (thread*2 + 1) * chunk2
          ito     = (thread*2 + 2) * chunk2

          middle = min(middle, size(A))
          ito     = min(ito, size(A))
          if (from < ito) then
             call C_merge(A, order, from, middle, ito)
          end if
       end do
       !$OMP end parallel do

       chunk2 = chunk2 * 2
       i = i + 1
    end do
  end subroutine C_parallel_sort

  !> Merge two parts of A, ordered by order from left to right
  !! around middle.
  subroutine C_merge (A, order, left, middle, right)
    character(*), intent(in), dimension(:) :: A
    integer, intent(inout), dimension(size(A)) :: order

    integer, intent(in) :: left, middle, right

    integer :: leftA, rightA, leftB, rightB
    integer :: iA, iB, i

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

  end subroutine C_merge

  subroutine D_parallel_sort (A, order)
    real(8), intent(in),  dimension(:) :: A
    integer, intent(out), dimension(size(A)) :: order

    integer :: ilen, from, middle, ito, nthreads, thread, chunk, chunk2, i, iremainder, extraThread 

    ilen      = size(A)
    nthreads = omp_get_max_threads()
    chunk    = ilen / nthreads

    iremainder = mod(ilen, nthreads)
    if (iremainder /= 0) then
        extraThread = 1
    else
        extraThread = 0
    endif   
    !----------------------------------------
    ! Initialize order
    !----------------------------------------
    do i = 1, ilen
       order(i) = i
    end do

    !----------------------------------------
    ! Sort each chunk
    !----------------------------------------
    !$OMP parallel do default(shared) private(thread, from, ito)
    do thread = 0, nthreads-1 + extraThread 
       from = thread*chunk + 1
       ito  = min((thread + 1)*chunk, ilen)

       call mrgrnk(A(from:ito), order(from:ito))
       order(from:ito) = order(from:ito) + from - 1
    end do
    !$OMP end parallel do

    !----------------------------------------
    ! Merge pieces together
    !----------------------------------------
    i = 1
    chunk2 = chunk
    do while (chunk2 < size(A))

       !$OMP parallel do default(shared) private(thread, from, middle, ito)
       do thread = 0, ceiling(.5 * size(A) / chunk2)
          from   = thread*2     * chunk2 + 1
          middle = (thread*2 + 1) * chunk2
          ito     = (thread*2 + 2) * chunk2

          middle = min(middle, size(A))
          ito     = min(ito, size(A))
          if (from < ito) then
             call D_merge(A, order, from, middle, ito)
          end if
       end do
       !$OMP end parallel do

       chunk2 = chunk2 * 2
       i = i + 1
    end do
  end subroutine D_parallel_sort

  !> Merge two parts of A, ordered by order from left to right
  !! around middle.
  subroutine D_merge (A, order, left, middle, right)
    real(8), intent(in), dimension(:) :: A
    integer, intent(inout), dimension(size(A)) :: order

    integer, intent(in) :: left, middle, right

    integer :: leftA, rightA, leftB, rightB
    integer :: iA, iB, i

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

  end subroutine D_merge

  subroutine R_parallel_sort (A, order)
    real(4), intent(in),  dimension(:) :: A
    integer, intent(out), dimension(size(A)) :: order

    integer :: ilen, from, middle, ito, nthreads, thread, chunk, chunk2, i, iremainder, extraThread 

    ilen      = size(A)
    nthreads = omp_get_max_threads()
    chunk    = ilen / nthreads

    iremainder = mod(ilen, nthreads)
    if (iremainder /= 0) then
        extraThread = 1
    else
        extraThread = 0
    endif   
    !----------------------------------------
    ! Initialize order
    !----------------------------------------
    do i = 1, ilen
       order(i) = i
    end do

    !----------------------------------------
    ! Sort each chunk
    !----------------------------------------
    !$OMP parallel do default(shared) private(thread, from, ito)
    do thread = 0, nthreads-1 + extraThread 
       from = thread*chunk + 1
       ito  = min((thread + 1)*chunk, ilen)

       call mrgrnk(A(from:ito), order(from:ito))
       order(from:ito) = order(from:ito) + from - 1
    end do
    !$OMP end parallel do

    !----------------------------------------
    ! Merge pieces together
    !----------------------------------------
    i = 1
    chunk2 = chunk
    do while (chunk2 < size(A))

       !$OMP parallel do default(shared) private(thread, from, middle, ito)
       do thread = 0, ceiling(.5 * size(A) / chunk2)
          from   = thread*2     * chunk2 + 1
          middle = (thread*2 + 1) * chunk2
          ito     = (thread*2 + 2) * chunk2

          middle = min(middle, size(A))
          ito     = min(ito, size(A))
          if (from < ito) then
             call R_merge(A, order, from, middle, ito)
          end if
       end do
       !$OMP end parallel do

       chunk2 = chunk2 * 2
       i = i + 1
    end do
  end subroutine R_parallel_sort

  !> Merge two parts of A, ordered by order from left to right
  !! around middle.
  subroutine R_merge (A, order, left, middle, right)
    real(4), intent(in), dimension(:) :: A
    integer, intent(inout), dimension(size(A)) :: order

    integer, intent(in) :: left, middle, right

    integer :: leftA, rightA, leftB, rightB
    integer :: iA, iB, i

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

  end subroutine R_merge

  subroutine I_parallel_sort (A, order)
    integer, intent(in),  dimension(:) :: A
    integer, intent(out), dimension(size(A)) :: order

    integer :: ilen, from, middle, ito, nthreads, thread, chunk, chunk2, i, iremainder, extraThread 

    ilen      = size(A)
    nthreads = omp_get_max_threads()
    chunk    = ilen / nthreads

    iremainder = mod(ilen, nthreads)
    if (iremainder /= 0) then
        extraThread = 1
    else
        extraThread = 0
    endif   
    !----------------------------------------
    ! Initialize order
    !----------------------------------------
    do i = 1, ilen
       order(i) = i
    end do

    !----------------------------------------
    ! Sort each chunk
    !----------------------------------------
    !$OMP parallel do default(shared) private(thread, from, ito)
    do thread = 0, nthreads-1 + extraThread 
       from = thread*chunk + 1
       ito  = min((thread + 1)*chunk, ilen)

       call mrgrnk(A(from:ito), order(from:ito))
       order(from:ito) = order(from:ito) + from - 1
    end do
    !$OMP end parallel do

    !----------------------------------------
    ! Merge pieces together
    !----------------------------------------
    i = 1
    chunk2 = chunk
    do while (chunk2 < size(A))

       !$OMP parallel do default(shared) private(thread, from, middle, ito)
       do thread = 0, ceiling(.5 * size(A) / chunk2)
          from   = thread*2     * chunk2 + 1
          middle = (thread*2 + 1) * chunk2
          ito     = (thread*2 + 2) * chunk2

          middle = min(middle, size(A))
          ito     = min(ito, size(A))
          if (from < ito) then
             call I_merge(A, order, from, middle, ito)
          end if
       end do
       !$OMP end parallel do

       chunk2 = chunk2 * 2
       i = i + 1
    end do
  end subroutine I_parallel_sort

  !> Merge two parts of A, ordered by order from left to right
  !! around middle.
  subroutine I_merge (A, order, left, middle, right)
    integer, intent(in), dimension(:) :: A
    integer, intent(inout), dimension(size(A)) :: order

    integer, intent(in) :: left, middle, right

    integer :: leftA, rightA, leftB, rightB
    integer :: iA, iB, i

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

  end subroutine I_merge

end module mod_sort
