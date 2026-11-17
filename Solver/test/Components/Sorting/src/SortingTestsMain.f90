!
!////////////////////////////////////////////////////////////////////////
!
!      SortingTestsMain.f90
!      Tests of the integer merge sorts in Utilities: the output must be sorted
!
!////////////////////////////////////////////////////////////////////////
!
      PROGRAM SortingTestsMain
         USE TestSuiteManagerClass
         USE SMConstants
         IMPLICIT NONE

         TYPE(TestSuiteManager) :: testSuite
         INTEGER                :: numberOfFailures
         EXTERNAL               :: testSortDescend, testSortAscend

         CALL testSuite % init()

         CALL testSuite % addTestSubroutineWithName(testSortDescend,"Merge sort, descending with permutation")
         CALL testSuite % addTestSubroutineWithName(testSortAscend ,"Merge sort, ascending")

         CALL testSuite % performTests(numberOfFailures)

         CALL testSuite % finalize()

         IF(numberOfFailures > 0)   error stop 99

      END PROGRAM SortingTestsMain
