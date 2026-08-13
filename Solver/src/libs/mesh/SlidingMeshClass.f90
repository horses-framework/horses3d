!
!////////////////////////////////////////////////////////////////////////////////
!
!     SlidingMeshClass
!     ----------------
!
!     Holds everything the sliding mesh needs to describe itself: the rotating
!     cylinder, the current angle, and the tables that pair the two sides of the
!     interface. The algorithms that use them live in SlidingMeshProcedures.
!
!     Read from the control file, in a "#define slidingmesh" block:
!
!        center          two coordinates, in the plane normal to the axis
!        radius          elements entirely inside it rotate
!        angle           rotation applied per time step, radians
!        rotation axis   x, y or z
!
!     Conventions
!     -----------
!     iAx is the component index of the rotation axis, iR1 and iR2 the two that
!     span the plane normal to it. They are NOT in the usual cyclic order:
!
!        axis x  ->  (iR1, iR2) = (z, y)
!        axis y  ->  (iR1, iR2) = (x, z)
!        axis z  ->  (iR1, iR2) = (y, x)
!
!     so phi = atan2(x(iR2) - center(2), x(iR1) - center(1)) runs in the opposite
!     sense to the usual right-handed convention. This is consistent throughout
!     both files.
!
!     center is given in that same (iR1, iR2) order, which is (z,y) for axis x,
!     (x,z) for axis y and (y,x) for axis z. Only the middle one matches the
!     intuitive ordering. Invisible for a cylinder on the axis, off by ninety
!     degrees for one that is not.
!
!     omega is the cumulative angle since the start, theta the increment per
!     step, both in radians. phi is returned in (-pi, pi].
!
!     Because (iR1, iR2) is anti-cyclic for all three axes, RotatePoint turning
!     by +theta in that plane is a rotation by -theta about +e_iAx, whichever
!     axis was chosen. OmegaVector is the one place that sign is written down;
!     anything needing a grid velocity must go through it rather than rebuild
!     the vector, so that flux and geometry cannot disagree on the direction.
!
!////////////////////////////////////////////////////////////////////////////////
!

#include "Includes.h"
MODULE SlidingMeshClass
      use Utilities                       , only: toLower, almostEqual, AlmostEqualRelax
      use SMConstants
      USE MeshTypes
      USE NodeClass
      use FileReaders               , only: ReadControlFile 
      use FileReadingUtilities      , only: getFileName, getRealArrayFromString
      use ElementConnectivityDefinitions
      use FTValueDictionaryClass

      IMPLICIT NONE

      enum, bind(C)
        enumerator :: rotationAxis_X, rotationAxis_Z, rotationAxis_Y
      end enum

!     ---------------
!     Mesh definition
!     ---------------
!
    type SlidingMesh
      logical                                :: active = .false.
      logical                                :: isConfigured = .false.
      integer                                :: numPureSlidingElems = 0
      integer,                   allocatable :: mortarNeighborElems(:)!associated non-sliding neighbors across the interface (require mortars)
      integer,                   allocatable :: slidingMortarElems(:)!sliding elements at the interface (require mortars)
      integer,                   allocatable :: pureSlidingElems(:)!sliding elements fully inside the region (no mortar)
      integer,                   allocatable :: mortararr1(:,:)
      integer,                   allocatable :: mortararr2(:,:)
      integer,                   allocatable :: face_nodes(:,:)
      integer,                   allocatable :: face_othernodes(:,:)
      integer,                   allocatable :: slidingMortarConnectivity(:,:)! slidingMortarConnectivity(i,:) :
      !   Data structure describing the mortar interface configuration for each slidingMortarElems(i). 
      !   Each row encodes the connectivity between:
      !   - a sliding interface element,
      !   - its neighboring sliding element,
      !   - and the associated non-sliding (mortar) element.
      !
      !   Column description:
      !     (1)  : ID of the sliding interface element (current element)
      !     (2)  : ID of the neighboring sliding element across the interface (left of (1))
      !     (3)  : ID of the neighboring non-sliding element (mortar element) (top of (2))
      !     (4)  : local face index of (1) connected to the mortar interface
      !     (5)  : local face index of (2) connected to the mortar interface
      !     (6)  : local face index of (3) connected to the mortar interface
      !     (7)  : rotation of face (1)
      !     (8)  : rotation of face (2)
      !     (9)  : rotation of face (3)
      !     (10) : ID of the second non-sliding master element connected to (1)
      !     (11) : corresponding local face index for (10)
      !     (12) : rotation associated with (10)
      !
      !   This structure is used to fully reconstruct mortar connectivity and orientation 
      !   between sliding and non-sliding regions.
      integer,                   allocatable :: neighborConnectivity(:,:,:)
      integer,                   allocatable :: rotmortars(:)
      integer                                :: numBFacePoints
      integer                                :: numSlidingElements
      integer                                :: numSlidingInterfaceElements 
      integer                                :: currentSectorID = 0
      real(KIND=RP)                          :: center(2) = 0.0_RP
      real(KIND=RP)                          :: radius = huge(1.0_RP)
      real(KIND=RP)                          :: theta = 0.0_RP       ! rotation angle
      real(KIND=RP)                          :: omega = 0.0_RP
      real(KIND=RP)                          :: localAngle
      real(KIND=RP)                          :: angularVelocity = 0.0_RP ! signed rotation rate, rad per (nondimensional)
                                                                        ! time unit, in the same sense as theta. Must be
                                                                        ! kept equal to theta/dt by whoever advances the
                                                                        ! mesh; it is NOT set here and NOT set by
                                                                        ! AdvanceSlidingMesh.
      integer                                :: rotationAxis
      integer                                :: iAx = 2 
      integer                                :: iR1 = 1
      integer                                :: iR2 = 3
      integer                                :: numElemsPerLayer = 0
      logical                                :: initialized=.false.
      logical                                :: conforming=.false.
      integer                                :: numDuplicatedNodes = 0
      integer                                :: numLayers = 0
      integer, allocatable                   :: layerID(:)
    contains

        procedure :: read_info                          => SlidingMesh_read_info
        procedure :: GetInfo                            => SlidingMesh_GetInfo
        procedure :: Initialize                         => SlidingMesh_Initialize
        procedure :: Destruct                           => SlidingMesh_Destruct
        procedure :: phi                                => SlidingMesh_phi
        procedure :: zeta                               => SlidingMesh_zeta
        procedure :: rad                                => SlidingMesh_rad
        procedure :: RotatePoint                        => SlidingMesh_RotatePoint
        procedure :: OmegaVector                        => SlidingMesh_OmegaVector
        procedure :: RotationCenter3D                   => SlidingMesh_RotationCenter3D
        procedure :: PlaneCoords                        => SlidingMesh_PlaneCoords
        procedure :: GridVelocityAt                     => SlidingMesh_GridVelocityAt

    end type

!     ========
      CONTAINS
!     ========

! ---------------------------------------------------------------------------
!  SlidingMesh_read_info
!
!  Entry point for the control file. Does nothing unless a "#define
!  slidingmesh" block is present, in which case it defers to GetInfo.
! --------------------------------------------------------------------------
subroutine SlidingMesh_read_info( self, controlVariables )
    use FileReadingUtilities
    implicit none
        
    class(SlidingMesh), intent(inout) :: self
    class(FTValueDictionary)       :: controlVariables

    if (SlidingMeshIsDefined()) then
        call self % GetInfo( controlVariables )
    end if

end subroutine SlidingMesh_read_info

! ---------------------------------------------------------------------------
!  SlidingMesh_GetInfo
!
!  Reads center, radius, angle and rotation axis from the "#define
!  slidingmesh" block, and derives iAx, iR1 and iR2 from the axis. Every
!  entry is mandatory: a missing or malformed one stops the run.
! ---------------------------------------------------------------------------
subroutine SlidingMesh_GetInfo( self, controlVariables )
    use FileReadingUtilities
    use ParamfileRegions
    implicit none
        
    class(SlidingMesh), intent(inout) :: self
    class(FTValueDictionary)       :: controlVariables
    
    character(len=LINE_LENGTH) :: in_label, paramFile
    character(len=LINE_LENGTH) :: R_center
    real(kind=RP),allocatable :: center(:)
    real(kind=RP),allocatable :: radius
    real(kind=RP),allocatable :: angle 
    character(len=LINE_LENGTH) :: rotAxis
    integer :: i

!     **********
!     Read block
!     **********

    write(in_label , '(A)') "#define slidingmesh"
    call get_command_argument(1, paramFile) 

    call readValueInRegion( trim( paramFile ), "center", R_center, in_label, "#end" )
    call readValueInRegion( trim( paramFile ), "radius", radius, in_label, "#end" )
    call readValueInRegion( trim( paramFile ), "angle", angle, in_label, "#end" )
    call readValueInRegion( trim( paramFile ), "rotation axis", rotAxis, in_label, "#end" )

    if (R_center .ne. "") then
        center = getRealArrayFromString(R_center)
    else
        print *, "center should be specified."
        errorMessage(STD_OUT)
        error stop
    end if
    if (size(center) .ne. 2) then
        print *, "The center should have two coordinates."
        errorMessage(STD_OUT)
        error stop
    endif
    self % center = center

    if (.not. allocated(radius)) then
        print *, "radius should be specified."
        errorMessage(STD_OUT)
        error stop
    else
        self % radius = radius 
    end if 
    if (.not. allocated(angle)) then 
        print *, "angle should be specified."
        errorMessage(STD_OUT)
        error stop
    else
        self % theta = angle 
    end if

    call toLower(rotAxis)
    select case (trim(rotAxis))
    case ("x")
        self % rotationAxis = rotationAxis_X
        self % iAx = 1
    case ("y")
        self % rotationAxis = rotationAxis_Y
        self % iAx = 2
    case ("z")
        self % rotationAxis = rotationAxis_Z
        self % iAx = 3
    case default
        print *, "Unrecognized rotation axis value. Possible values are: x, y, z."
        errorMessage(STD_OUT)
        error stop
    end select   
    
    self % iR1 = mod(self % iAx + 1, 3) + 1
    self % iR2 = mod(self % iAx, 3) + 1
    self % isConfigured = .true.

end subroutine SlidingMesh_GetInfo

! ---------------------------------------------------------------------------
!  SlidingMesh_Initialize
!
!  Allocates the connectivity tables once the interface has been sized, and
!  zeroes them on the first pass only. Allocation is guarded, so calling it
!  again on a later step is harmless.
! ---------------------------------------------------------------------------
subroutine SlidingMesh_Initialize(self, numSlidingInterfaceElements, numSlidingElements, numBFacePoints)
    use Physics
    use PartitionedMeshClass
    use MPI_Process_Info

    implicit none

    class(SlidingMesh), intent(inout)        :: self
    integer, intent(in)                      :: numSlidingInterfaceElements
    integer, intent(in)                      :: numSlidingElements
    integer, intent(in)                      :: numBFacePoints
   !========================
   ! Allocation
   !========================
    if (.not. allocated(self % mortarNeighborElems))             allocate(self % mortarNeighborElems(numSlidingInterfaceElements))
    if (.not. allocated(self % slidingMortarElems))              allocate(self % slidingMortarElems(numSlidingInterfaceElements))
    if (.not. allocated(self % mortararr1))                      allocate(self % mortararr1(numSlidingInterfaceElements, 2))
    if (.not. allocated(self % mortararr2))                      allocate(self % mortararr2(numSlidingInterfaceElements, 2))
    if (.not. allocated(self % pureSlidingElems))                allocate(self % pureSlidingElems(numSlidingElements))
    if (.not. allocated(self % slidingMortarConnectivity))       allocate(self % slidingMortarConnectivity(numSlidingInterfaceElements, 12))
    if (.not. allocated(self % face_nodes))                      allocate(self % face_nodes(numSlidingInterfaceElements, 4))
    if (.not. allocated(self % face_othernodes))                 allocate(self % face_othernodes(numSlidingInterfaceElements, 4))
    if (.not. allocated(self % neighborConnectivity))            allocate(self % neighborConnectivity(numSlidingInterfaceElements, 9, 6))
    if (.not. allocated(self % rotmortars))                      allocate(self % rotmortars(2 * numSlidingInterfaceElements))
    if (.not. allocated(self % layerID))                         allocate(self % layerID(numSlidingInterfaceElements))
    
    !========================
    ! Initialization (first pass only)
    !========================
    if (.not. self % active) then
       self % mortarNeighborElems       = 0
       self % slidingMortarElems        = 0
       self % pureSlidingElems          = 0
       self % slidingMortarConnectivity = 0
       self % neighborConnectivity      = 0
       self % face_nodes                = 0
       self % face_othernodes           = 0
       self % localAngle                = 0
       self % currentSectorID           = 0
       self % numSlidingInterfaceElements = numSlidingInterfaceElements
       self % numSlidingElements          = numSlidingElements
       self % numBFacePoints              = numBFacePoints
       self % layerID                   = 0
    end if 

end subroutine SlidingMesh_Initialize

! ---------------------------------------------------------------------------
!  SlidingMesh_Destruct
!
!  Releases every table allocated by Initialize. Uses safedeallocate, so it
!  is safe on a mesh where the sliding interface was never active.
! ---------------------------------------------------------------------------
subroutine SlidingMesh_Destruct(self)
    implicit none

    class(SlidingMesh), intent(inout) :: self

    ! ---------------------------------------------------------------------
    ! Deallocate sliding mesh connectivity and mortar storage
    ! ---------------------------------------------------------------------
    safedeallocate(self % mortarNeighborElems)
    safedeallocate(self % slidingMortarElems)
    safedeallocate(self % pureSlidingElems)
    safedeallocate(self % mortararr1)
    safedeallocate(self % mortararr2)
    safedeallocate(self % face_nodes)
    safedeallocate(self % face_othernodes)
    safedeallocate(self % slidingMortarConnectivity)
    safedeallocate(self % neighborConnectivity)
    safedeallocate(self % rotmortars)
    safedeallocate(self % layerID)
    ! ---------------------------------------------------------------------
    ! Reset counters and state variables
    ! ---------------------------------------------------------------------

    self % numSlidingInterfaceElements = 0
    self % numSlidingElements          = 0
    self % numBFacePoints              = 0

    self % omega        = 0.0_RP
    self % theta        = 0.0_RP

    self % initialized  = .false.
    self % conforming   = .false.
    self % active       = .false.

end subroutine SlidingMesh_Destruct

! ---------------------------------------------------------------------------
!  SlidingMeshIsDefined
!
!  Scans the case file for a "#define slidingmesh" block and returns whether
!  one was found. Called before anything is read, so that a case without a
!  sliding interface pays nothing.
! ---------------------------------------------------------------------------
logical function SlidingMeshIsDefined()
    use ParamfileRegions
    implicit none

    ! ---------------
    ! Local variables
    ! ---------------

    character(len=LINE_LENGTH) :: case_name, line
    integer                    :: fID
    integer                    :: io

    ! Initialize
    ! ----------
    SlidingMeshIsDefined = .FALSE.

    ! Get case file name
    ! ------------------
    call get_command_argument(1, case_name)


    ! Open case file
    ! --------------
    open ( newunit = fID , file = case_name , status = "old" , action = "read" )

    ! Read the whole file to find sliding mesh region
    ! -------------------------------------------------
    readloop:do
        read ( fID , '(A)' , iostat = io ) line

        if ( io .lt. 0 ) then

            ! End of file
            ! -----------
            line = ""
            exit readloop

        elseif ( io .gt. 0 ) then

            ! Error
            ! -----
            errorMessage(STD_OUT)
            error stop "Stopped."

        else

            ! Succeeded
            ! ---------
            line = getSquashedLine( line )
            call toLower(line)

            if ( index ( line , '#defineslidingmesh' ) .gt. 0 ) then
               SlidingMeshIsDefined = .TRUE.
               close(fID)  
               return
            end if
            
        end if
         
    end do readloop

    ! Close case file
    ! ---------------
    close(fID)                             

end function SlidingMeshIsDefined

!  Azimuthal angle of x about the axis, in (-pi, pi].
pure function SlidingMesh_phi(self, x) result(p)
   class(SlidingMesh), intent(in) :: self
   real(kind=RP),      intent(in) :: x(3)
   real(kind=RP) :: p

   p = atan2( x(self % iR2) - self % center(2), x(self % iR1) - self % center(1) )

end function

!  Coordinate of x along the rotation axis.
pure function SlidingMesh_zeta(self, x) result(z)
   class(SlidingMesh), intent(in) :: self
   real(kind=RP),      intent(in) :: x(3)
   real(kind=RP) :: z

   z = x(self % iAx)

end function

!  Distance from x to the axis.
pure function SlidingMesh_rad(self, x) result(r)
   class(SlidingMesh), intent(in) :: self
   real(kind=RP),      intent(in) :: x(3)
   real(kind=RP) :: r

   r = sqrt( (x(self % iR1) - self % center(1))**2 + (x(self % iR2) - self % center(2))**2 )

end function

!  Rotates x about the axis by theta.
pure subroutine SlidingMesh_RotatePoint(self, theta, x)
   class(SlidingMesh), intent(in)    :: self
   real(kind=RP),      intent(in)    :: theta
   real(kind=RP),      intent(inout) :: x(3)

   real(kind=RP) :: u, v, c, s

   c = cos(theta) ; s = sin(theta)
   u = x(self % iR1) - self % center(1)
   v = x(self % iR2) - self % center(2)

   x(self % iR1) = self % center(1) + c*u - s*v
   x(self % iR2) = self % center(2) + s*u + c*v

end subroutine

!  -----------------------------------------------------------------------
!  Angular velocity vector of the sliding region, in the same sign
!  convention as RotatePoint. That convention is uniform across the three
!  axes: turning by +theta in the anti-cyclic (iR1, iR2) plane is a
!  rotation by -theta about +e_iAx, hence the minus sign below and the
!  absence of any per-axis branch.
!
!  This is the only place the sign lives. Any grid-velocity user must go
!  through it so that flux and geometry can never disagree on the
!  direction.
!  -----------------------------------------------------------------------
pure function SlidingMesh_OmegaVector(self) result(omg)
    implicit none
    class(SlidingMesh), intent(in) :: self
    real(kind=RP)                  :: omg(3)

    omg = 0.0_RP
    omg(self % iAx) = -self % angularVelocity

end function SlidingMesh_OmegaVector

!  -----------------------------------------------------------------------
!  The two control-file center coordinates, placed back into 3D. They are
!  given in the (iR1, iR2) order, the same one RotatePoint rotates about,
!  so a nonzero center is consistent between the geometry update and the
!  fluxes.
!  -----------------------------------------------------------------------
pure function SlidingMesh_RotationCenter3D(self) result(xc)
    implicit none
    class(SlidingMesh), intent(in) :: self
    real(kind=RP)                  :: xc(3)

    xc = 0.0_RP
    xc(self % iR1) = self % center(1)
    xc(self % iR2) = self % center(2)

end function SlidingMesh_RotationCenter3D

!  -----------------------------------------------------------------------
!  In-plane coordinate indices of the rotation plane. Kept as a named
!  accessor for code that reads better that way; it is exactly (iR1, iR2)
!  and carries the same anti-cyclic ordering, so ih and iv are (z,y) for
!  axis x, (x,z) for axis y and (y,x) for axis z. Note that for x and z
!  this is the reverse of the intuitive pairing.
!  -----------------------------------------------------------------------
pure subroutine SlidingMesh_PlaneCoords(self, ih, iv)
    implicit none
    class(SlidingMesh), intent(in)  :: self
    integer,            intent(out) :: ih, iv

    ih = self % iR1
    iv = self % iR2

end subroutine SlidingMesh_PlaneCoords

!  -----------------------------------------------------------------------
!  Grid velocity v_g = Omega x (x - xc) at a physical point x
!  -----------------------------------------------------------------------
pure function SlidingMesh_GridVelocityAt(self, x) result(vg)
    implicit none
    class(SlidingMesh), intent(in) :: self
    real(kind=RP),      intent(in) :: x(3)
    real(kind=RP)                  :: vg(3)

    real(kind=RP) :: omg(3), r(3)

    omg = self % OmegaVector()
    r   = x - self % RotationCenter3D()

    vg(1) = omg(2)*r(3) - omg(3)*r(2)
    vg(2) = omg(3)*r(1) - omg(1)*r(3)
    vg(3) = omg(1)*r(2) - omg(2)*r(1)

end function SlidingMesh_GridVelocityAt

END MODULE SlidingMeshClass
