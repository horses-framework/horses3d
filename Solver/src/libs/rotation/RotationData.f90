#include "Includes.h"
module RotationData
!
!  ***********************************************************************
!    This module implements periodic revolution BCs: the rotation
!    geometry (axis, center) and the periodic sector angle.
!
!  Equation-set-agnostic: this module is built once (see libs/rotation),
!  not once per equation set. The "only supported in Navier-Stokes"
!  check therefore lives in the equation-set-specific caller (DGSEMClass),
!  not here.
!  ***********************************************************************
!
   use SMConstants
   implicit none

   private
   public :: RotationParams_t, rotationParams
   public :: RotationAxisIndex
   public :: SetPeriodicAngle
   public :: InitializeRotationFromControl, DescribeRotation
   public :: SectorAngle, SyncPeriodicAngle
   public :: CompareTwoNodesRevol, CompareTwoNodesRevolRelax

   type RotationParams_t

      ! Rotation geometry
      real(kind=RP) :: axis(NDIM) = [0.0_RP, 0.0_RP, 1.0_RP]    ! Unit rotation axis
      real(kind=RP) :: center(NDIM) = [0.0_RP, 0.0_RP, 0.0_RP]  ! Center

      ! For periodic revolution BCs
      logical       :: periodicEnabled = .false.
      real(kind=RP) :: periodicAngle = 0.0_RP  ! Sector (rad)
   end type RotationParams_t

   type(RotationParams_t), protected :: rotationParams

   contains

   subroutine DescribeRotation(nRotaryFaces)
!
!        Display the periodic revolution settings (root rank only)
!
      use Headers,          only: Section_Header
      use MPI_Process_Info, only: MPI_Process
      implicit none
      integer,       intent(in) :: nRotaryFaces    ! Global number of rotary periodic faces

      if (.not. MPI_Process % isRoot) return

      if (rotationParams % periodicEnabled) then
         call Section_Header("Periodic Revolution BCs")
         write(STD_OUT,'(30X,A,A30,F8.3)')   "->","Sector angle (deg): ", rotationParams % periodicAngle * 180.0_RP / PI
         write(STD_OUT,'(30X,A,A30,I0)') "->","Rotary periodic faces: ", nRotaryFaces
         write(STD_OUT,'(30X,A,A30,3ES10.3)') "->","Rotation axis: ", rotationParams % axis
         write(STD_OUT,'(30X,A,A30,3ES10.3)') "->","Rotation center: ", rotationParams % center
         write(STD_OUT,'(/)')
      end if

   end subroutine DescribeRotation

   pure integer function RotationAxisIndex()
!     Return the index of the rotation axis if it is aligned with x, y, or z.
      real(kind=RP) :: axisAbs(NDIM)
      integer       :: idx

      axisAbs = abs(rotationParams % axis)
      idx = maxloc(axisAbs, dim = 1)
      if (axisAbs(idx) < (1.0_RP - 1.0e-8_RP)) then
         RotationAxisIndex = 0
      else
         RotationAxisIndex = idx
      end if
   end function RotationAxisIndex

   subroutine SetPeriodicAngle(angle_)
!
!     Store the sector angle (rad) computed during periodic face matching
!
      implicit none
      real(kind=RP), intent(in) :: angle_
      rotationParams % periodicAngle = angle_
   end subroutine SetPeriodicAngle

   subroutine InitializeRotationFromControl(controlVariables)
!
!     Read rotation geometry (axis, center) and enable/disable flags
!     from the control file.  Called early during mesh construction so that
!     periodic face matching can use the rotation information.
!
      use FTValueDictionaryClass
      use FileReadingUtilities, only: getRealArrayFromString
      implicit none
      class(FTValueDictionary), intent(in) :: controlVariables
      real(kind=RP) :: axis_(NDIM), center_(NDIM)

      axis_   = rotationParams % axis
      center_ = rotationParams % center

      ! ---- rotation axis ----
      if (controlVariables % containsKey("rotation axis")) then
         axis_ = getRealArrayFromString(controlVariables % stringValueForKey("rotation axis", LINE_LENGTH))
      end if

      ! ---- rotation center ----
      if (controlVariables % containsKey("rotation center")) then
         center_ = getRealArrayFromString(controlVariables % stringValueForKey("rotation center", LINE_LENGTH))
      end if

      rotationParams % axis   = axis_ / norm2(axis_)
      rotationParams % center = center_

      ! ---- periodic enabled ----
      if (controlVariables % containsKey("periodic enabled")) then
         rotationParams % periodicEnabled = .false.
         if (controlVariables % logicalValueForKey("periodic enabled")) then
            if (RotationAxisIndex() /= 0) then
               rotationParams % periodicEnabled = .true.
            end if
         end if
      end if

   end subroutine InitializeRotationFromControl




   real(kind=RP) function SectorAngle(x1, x2)
!
!     Angle that rotates x1 into x2 around the rotation axis (wrapped into (-PI, PI]).
!     Used to get the periodic sector angle from a matched node pair.
!
      implicit none
      real(kind=RP), intent(in) :: x1(NDIM), x2(NDIM)
      real(kind=RP)             :: x0(NDIM)

      x0 = rotationParams % center
      select case (RotationAxisIndex())
      case (1) ! Rotation around x-axis, use y-z plane
         SectorAngle = atan2(x2(3) - x0(3), x2(2) - x0(2)) - atan2(x1(3) - x0(3), x1(2) - x0(2))
      case (2) ! Rotation around y-axis, use z-x plane
         SectorAngle = atan2(x2(1) - x0(1), x2(3) - x0(3)) - atan2(x1(1) - x0(1), x1(3) - x0(3))
      case (3) ! Rotation around z-axis, use x-y plane
         SectorAngle = atan2(x2(2) - x0(2), x2(1) - x0(1)) - atan2(x1(2) - x0(2), x1(1) - x0(1))
      case default
         print *, "Error: invalid rotation axis index when computing periodic sector angle."
         error stop
      end select

      SectorAngle = SectorAngle - 2.0_RP * PI * anint(SectorAngle / (2.0_RP * PI))

   end function SectorAngle

   subroutine SyncPeriodicAngle()
!
!     The sector angle is computed by the root rank only (first call to ConstructPeriodicFaces):
!     broadcast it to all the ranks.
!
#ifdef _HAS_MPI_
      use mpi
#endif
      use MPI_Process_Info, only: MPI_Process
      implicit none
#ifdef _HAS_MPI_
      real(kind=RP) :: angle
      integer       :: ierr

      if ( MPI_Process % doMPIAction ) then
         angle = rotationParams % periodicAngle
         call mpi_bcast(angle, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
         call SetPeriodicAngle(angle)
      end if
#endif
   end subroutine SyncPeriodicAngle

   SUBROUTINE CompareTwoNodesRevol(x1, x2, success, coord, x0)
   USE Utilities, only: AlmostEqual
   IMPLICIT NONE
!
!-------------------------------------------------------------------
! Similar to CompareTwoNodes, but compares with respect to revolution
! radius around axis `coord` and center `x0`.
!-------------------------------------------------------------------
!
   REAL(KIND=RP) :: x1(3)
   REAL(KIND=RP) :: x2(3)
   LOGICAL       :: success
   INTEGER       :: coord
   REAL(KIND=RP) :: x0(3)
!
   INTEGER       :: i
   REAL(KIND=RP) :: r1
   REAL(KIND=RP) :: r2

   IF (coord == 0) THEN
      success = .FALSE.
   ELSE
      r1 = 0.0_RP
      r2 = 0.0_RP
      DO i = 1,3
         IF (i /= coord) THEN
            r1 = r1 + POW2(x1(i) - x0(i))
            r2 = r2 + POW2(x2(i) - x0(i))
         ENDIF
      ENDDO

      IF ( AlmostEqual(r1, r2) .and. AlmostEqual(x1(coord), x2(coord)) ) THEN
         success = .TRUE.
      ELSE
         success = .FALSE.
      ENDIF
   ENDIF

   END SUBROUTINE CompareTwoNodesRevol
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE CompareTwoNodesRevolRelax(x1, x2, success, coord, x0, min_edge_length)
   USE Utilities, only: AlmostEqualRelax
   IMPLICIT NONE
!
!-------------------------------------------------------------------
! Similar to CompareTwoNodesRevol, but with relaxed tolerance.
!-------------------------------------------------------------------
!
   REAL(KIND=RP) :: x1(3)
   REAL(KIND=RP) :: x2(3)
   LOGICAL       :: success
   INTEGER       :: coord
   REAL(KIND=RP) :: x0(3)
   REAL(KIND=RP) :: min_edge_length
!
   INTEGER       :: i
   REAL(KIND=RP) :: r1
   REAL(KIND=RP) :: r2

   IF (coord == 0) THEN
      success = .FALSE.
   ELSE
      r1 = 0.0_RP
      r2 = 0.0_RP
      DO i = 1,3
         IF (i /= coord) THEN
            r1 = r1 + POW2(x1(i) - x0(i))
            r2 = r2 + POW2(x2(i) - x0(i))
         ENDIF
      ENDDO

      IF ( AlmostEqualRelax(r1, r2, min_edge_length) .and. &
           AlmostEqualRelax(x1(coord), x2(coord), min_edge_length) ) THEN
         success = .TRUE.
      ELSE
         success = .FALSE.
      ENDIF
   ENDIF

   END SUBROUTINE CompareTwoNodesRevolRelax

end module RotationData
