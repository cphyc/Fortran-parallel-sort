!> Test the parallel sort routine
program Test_Sort_Parallel
  use mod_sort
  use m_mrgrnk

  integer, parameter :: N = 97 !95 !96 ! int(1e6) ! Arbitrary length of arrays
  integer, parameter :: Cwidth = 50        ! Arbitrary width of strings

  integer :: A(N), order(N), orderSerial(N), irand(N, Cwidth) 
  real(8) :: dpA(N), xrand(N), yrand(N, Cwidth)
  real(4) :: spA(N)
  character(Cwidth) :: cA(N)

  integer :: i, j
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
  spA(:) = 1.0d3 * xrand(:)
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

end program Test_Sort_Parallel
