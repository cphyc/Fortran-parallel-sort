! Parallel sorting routines by: Corentin Cadiou and Steve Rivkin 

!> Test the parallel sort routine
program Test_Sort_Parallel
  use ISO_FORTRAN_ENV, only : ERROR_UNIT
  use mod_sort
  use m_mrgrnk

  implicit none

  integer, parameter :: N = 97 !95 !96 ! int(1e6) ! Arbitrary length of arrays
  integer, parameter :: Cwidth = 50        ! Arbitrary width of strings

  integer :: A(N), order(N), orderSerial(N), irand(N, Cwidth)
  real(8) :: dpA(N), xrand(N), yrand(N, Cwidth)
  real(4) :: spA(N)
  character(Cwidth) :: cA(N)

  integer :: i, j

  logical :: ok
  integer :: Nerr = 0
!  ------------------------------------
  open(8, file = 'Sort_Parallel.log')
!
! Create test data for 4 kinds of sorting:
!
! Test data for: character strings
!
  call random_number(yrand(:,:))
  irand(:,:) = 48 + nint((122 - 48)*yrand(:,:))

  where ((irand(:,:) >= 58) .and. (irand(:,:) <= 64)) ! Remove some strange characters
     irand(:,:) = 65
  end where

  where ((irand(:,:) >= 91) .and. (irand(:,:) <= 96)) ! Remove some strange characters
     irand(:,:) = 97
  end where

  do i = 1, N
     cA(i) = ''
     do j = 1, Cwidth
        cA(i) = trim(cA(i))//achar(irand(i, j)) 
     end do
  end do
!
! Test data for: double precision, single precision, integer  
!
  call random_number(xrand(:))

  dpA(:) = 1.0d3 * xrand(:)
  spA(:) = real(dpA, 4)
  A(:) = nint(1.0d3 * xrand(:))
!
!  ------------------------------------
!  Parallel sort character strings
!  ------------------------------------
  call mrgrnk(cA, orderSerial)
  call parallel_sort(cA, order)

  write(8, *) 'Sort character strings:'
  do i = 1, N
     write(8, '(1x,i8, 3(1x,a))') i, cA(i), cA(order(i)), cA(orderSerial(i))
  end do

  ok = all(cA(order) == cA(orderSerial))
  if (.not. ok) then
     write(ERROR_UNIT, *) 'An error ocurred while sorting character strings'
     Nerr = Nerr + 1
  end if
!
!  ------------------------------------
!  Parallel sort double precision
!  ------------------------------------
  call mrgrnk(dpA, orderSerial)
  call parallel_sort(dpA, order)

  write(8, *) 'Sort double precision:'
  do i = 1, N
     write(8, '(1x,i8, 3(1x,G15.7))') i, dpA(i), dpA(order(i)), dpA(orderSerial(i))
  end do

  ok = all(dpA(order) == dpA(orderSerial))
  if (.not. ok) then
     write(ERROR_UNIT, *) 'An error ocurred while sorting double precision floats'
     Nerr = Nerr + 1
  end if
!
!  ------------------------------------
!  Parallel sort single precision
!  ------------------------------------
  call mrgrnk(spA, orderSerial)
  call parallel_sort(spA, order)

  write(8, *) 'Sort single precision:'
  do i = 1, N
     write(8, '(1x,i8, 3(1x,G15.7))') i, spA(i), spA(order(i)), spA(orderSerial(i))
  end do

  ok = all(spA(order) == spA(orderSerial))
  if (.not. ok) then
     write(ERROR_UNIT, *) 'An error ocurred while sorting single precision floats'
     Nerr = Nerr + 1
  end if
!
!  ------------------------------------
!  Parallel sort integer
!  ------------------------------------
  call mrgrnk(A, orderSerial)
  call parallel_sort(A, order)

  write(8, *) 'Sort integer:'
  do i = 1, N
     write(8, *) i, A(i), A(order(i)), A(orderSerial(i))
  end do

  ok = all(A(order) == A(orderSerial))
  if (.not. ok) then
     write(ERROR_UNIT, *) 'An error ocurred while sorting integers'
     Nerr = Nerr + 1
  end if

!
!  ------------------------------------
!  Check parallel sort for small arrays
!  ------------------------------------
  block
   integer, allocatable :: A(:)
   integer, allocatable :: order(:), orderSerial(:)
   integer :: N
   do N = 1, 100
      allocate(A(N), order(N), orderSerial(N))

      do i = 1, N
         A(i) = N - i
      end do

      call parallel_sort(A, order)
      call mrgrnk(A, orderSerial)

      ok = all(A(order) == A(orderSerial))
      if (.not. ok) then
         write(ERROR_UNIT, *) 'An error ocurred while sorting array with length', N
         Nerr = Nerr + 1
      end if
      deallocate(A, order, orderSerial)
   end do
  end block

!
!  ------------------------------------
!  Timings
!  ------------------------------------

  block
    use omp_lib

    real(8) :: tstart, tend, serial_time, parallel_time, speed_up, efficiency
    real(8), allocatable :: AA(:)
    integer, allocatable :: Aorder(:)
    integer :: Niter = 100
    integer :: nthreads
    integer(kind=8) :: count, count_rate, count_max

    allocate(AA(1:1000000))
    allocate(Aorder(1:size(AA)))
    call random_number(AA)

    !-------------- Time serial version --------------
    call system_clock(count, count_rate=count_rate, count_max=count_max)
    tstart = dble(count) / count_rate

    do i = 1, Niter
      call mrgrnk(AA, Aorder)
    end do
    call system_clock(count, count_rate=count_rate, count_max=count_max)
    tend = dble(count) / count_rate

    serial_time = 1000*(tend-tstart) / Niter
    print*, "mrgrnk took", serial_time, "ms/iter"

    !-------------- Time parallel version --------------
    call system_clock(count, count_rate=count_rate, count_max=count_max)
    tstart = dble(count) / count_rate

    do i = 1, Niter
      call parallel_sort(AA, Aorder)
    end do

    call system_clock(count, count_rate=count_rate, count_max=count_max)
    tend = dble(count) / count_rate

    parallel_time = 1000*(tend-tstart) / Niter
    print*, "parallel sort took", parallel_time, "ms/iter"

    !-------------- Statistics -------------------------
    speed_up = serial_time / parallel_time
    print*, "speed up factor ", speed_up

    nthreads = omp_get_max_threads()
    print*, "using ", nthreads, " threads"

    efficiency = speed_up / nthreads 
    print*, "efficiency is ", 100*efficiency, "percent"
    print*, "Note this assumes the serial version is unvectorized."
    print*, "However it probably is auto-vectorized which contributes to the efficiency appearing low."
  end block

  write(*, *) Nerr, " of 8 sorting algorithms have errors."
end program Test_Sort_Parallel
