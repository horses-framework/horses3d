!
!////////////////////////////////////////////////////////////////////////
!
!      SortingTests.f90
!      Checks that mergeSortDescendInt / mergeSortAscendInt return sorted output: ordered keys, the same keys as the
!      input, and (descending, with a second array) a permutation that keeps equal keys in their original order.
!      Inputs are pseudo-random keys from a small range (many ties), from a fixed linear congruential sequence.
!
!////////////////////////////////////////////////////////////////////////
!
      module SortingTestData
         implicit none
         integer, parameter :: SIZES(6) = [0, 1, 2, 7, 100, 1000]
         integer, parameter :: MAX_KEY  = 5

      contains

         subroutine randomKeys(n, seed, keys)
            integer, intent(in)    :: n
            integer, intent(inout) :: seed
            integer, intent(out)   :: keys(n)
            integer                :: i
            do i = 1, n
               seed    = mod(1103515245_8 * seed + 12345_8, 2147483648_8)
               keys(i) = 1 + mod(seed / 65536, MAX_KEY)
            end do
         end subroutine randomKeys

         logical function sameKeys(a, b)
!           Same values with the same multiplicity (keys are in 1..MAX_KEY)
            integer, intent(in) :: a(:), b(:)
            integer             :: k
            sameKeys = size(a) == size(b)
            do k = 1, MAX_KEY
               if ( sum(merge(1, 0, a == k)) /= sum(merge(1, 0, b == k)) ) sameKeys = .false.
            end do
         end function sameKeys

      end module SortingTestData

      SUBROUTINE testSortDescend
         USE FTAssertions
         USE Utilities, only: mergeSortDescendInt
         USE SortingTestData
         IMPLICIT NONE
         integer              :: s, n, i, seed
         integer, allocatable :: keys(:), A(:), B(:)
         logical              :: ordered, stable, isPermutation
         logical, allocatable :: seen(:)
         character(len=64)    :: msg

         seed = 2026
         do s = 1, size(SIZES)
            n = SIZES(s)
            allocate(keys(n), seen(n))
            call randomKeys(n, seed, keys)
            A = keys
            B = [(i, i = 1, n)]
            call mergeSortDescendInt(A, B)

            ordered = .true. ; stable = .true.
            do i = 2, n
               if ( A(i) > A(i-1) ) ordered = .false.
               if ( (A(i) == A(i-1)) .and. (B(i) < B(i-1)) ) stable = .false.
            end do
            seen = .false. ; isPermutation = .true.
            do i = 1, n
               if ( (B(i) < 1) .or. (B(i) > n) ) then
                  isPermutation = .false.
               else if ( seen(B(i)) ) then
                  isPermutation = .false.
               else
                  seen(B(i)) = .true.
               end if
            end do

            write(msg,'(A,I5)') "Descend: keys in descending order, n =", n
            call FTAssert(test = ordered, msg = trim(msg))
            write(msg,'(A,I5)') "Descend: second array is a permutation, n =", n
            call FTAssert(test = isPermutation, msg = trim(msg))
            if ( isPermutation ) then
               write(msg,'(A,I5)') "Descend: permutation maps to the sorted keys, n =", n
               call FTAssert(test = all(keys(B) == A), msg = trim(msg))
            end if
            write(msg,'(A,I5)') "Descend: equal keys keep their order, n =", n
            call FTAssert(test = stable, msg = trim(msg))
            deallocate(keys, seen)
         end do

      END SUBROUTINE testSortDescend

      SUBROUTINE testSortAscend
         USE FTAssertions
         USE Utilities, only: mergeSortAscendInt
         USE SortingTestData
         IMPLICIT NONE
         integer              :: s, n, i, seed
         integer, allocatable :: keys(:), A(:)
         logical              :: ordered
         character(len=64)    :: msg

         seed = 2027
         do s = 1, size(SIZES)
            n = SIZES(s)
            allocate(keys(n))
            call randomKeys(n, seed, keys)
            A = keys
            call mergeSortAscendInt(A)

            ordered = .true.
            do i = 2, n
               if ( A(i) < A(i-1) ) ordered = .false.
            end do

            write(msg,'(A,I5)') "Ascend: keys in ascending order, n =", n
            call FTAssert(test = ordered, msg = trim(msg))
            write(msg,'(A,I5)') "Ascend: same keys as the input, n =", n
            call FTAssert(test = sameKeys(A, keys), msg = trim(msg))
            deallocate(keys)
         end do

      END SUBROUTINE testSortAscend
