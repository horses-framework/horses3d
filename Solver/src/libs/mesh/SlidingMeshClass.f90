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
      integer,                   allocatable :: mortarNeighborElems(:)!associated non-sliding neighbors across the interface (require mortars)
      integer,                   allocatable :: slidingMortarElems(:)!sliding elements at the interface (require mortars)
      integer,                   allocatable :: pureSlidingElems(:)!sliding elements fully inside the region (no mortar)
      integer,                   allocatable :: mortararr1(:,:)
      integer,                   allocatable :: mortararr2(:,:)
      integer,                   allocatable :: face_nodes(:,:)
      integer,                   allocatable :: face_othernodes(:,:)
      integer,                   allocatable :: slidingMortarConnectivity(:,:)! slidingMortarConnectivity(i,:) :
      !   Data structure describing the mortar interface configuration for each slidingMortarElems(i). Each row encodes the connectivity between:
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
      !   This structure is used to fully reconstruct mortar connectivity and orientation between sliding and non-sliding regions.
      integer,                   allocatable :: neighborConnectivity(:,:,:)
      integer,                   allocatable :: rotmortars(:)
      integer                                :: numBFacePoints
      integer                                :: numSlidingElements
      integer                                :: numSlidingInterfaceElements
      integer                                :: currentSectorID
      real(KIND=RP)                          :: center(2)
      real(KIND=RP)                          :: radius 
      real(KIND=RP)                          :: theta        ! rotation angle
      real(KIND=RP)                          :: omega
      real(KIND=RP)                          :: localAngle
      real(KIND=RP)                          :: angularVelocity = 0.0_RP ! rad per (nondimensional) time unit; magnitude of the
                                                                         ! rotation rate actually applied to the geometry. Must be
                                                                         ! kept equal to theta/dt by whoever advances the mesh.
      integer                                :: rotationAxis
      logical                                :: initialized=.false.
      logical                                :: conforming=.false.

    contains

        procedure :: read_info                              => SlidingMesh_read_info
        procedure :: GetInfo                                => SlidingMesh_GetInfo
        procedure :: Initialize                         => SlidingMesh_Initialize
        procedure :: Destruct                           => SlidingMesh_Destruct
        procedure :: OmegaVector                        => SlidingMesh_OmegaVector
        procedure :: RotationCenter3D                   => SlidingMesh_RotationCenter3D
        procedure :: GridVelocityAt                     => SlidingMesh_GridVelocityAt
        procedure :: PlaneCoords                        => SlidingMesh_PlaneCoords

    end type

!     ========
      CONTAINS
!     ========

!  -------------------------------------------------
!  SlidingMesh info
!  -------------------------------------------------   
subroutine SlidingMesh_read_info( self, controlVariables )
    use FileReadingUtilities
    implicit none
        
    class(SlidingMesh), intent(inout) :: self
    class(FTValueDictionary)       :: controlVariables

    if (SlidingMeshIsDefined()) then
        call self % GetInfo( controlVariables )
    end if
     
end subroutine SlidingMesh_read_info

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
    case ("y")
        self % rotationAxis = rotationAxis_Y
    case ("z")
        self % rotationAxis = rotationAxis_Z
    case default
        print *, "Unrecognized rotation axis value. Possible values are: x, y, z."
        errorMessage(STD_OUT)
        error stop
    end select

    
    self % isConfigured = .true.

end subroutine SlidingMesh_GetInfo

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
    
    end if 

end subroutine SlidingMesh_Initialize

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

!  -----------------------------------------------------------------------
!  Angular velocity vector of the sliding region, in the same sign
!  convention as the rotation matrices applied by RotateSlidingRegion:
!  a positive control-file angle corresponds to R_x(+theta), R_z(+theta)
!  and R_y(-theta) per step. Any grid-velocity user must go through this
!  function so that flux and geometry can never disagree on the direction.
!  -----------------------------------------------------------------------
pure function SlidingMesh_OmegaVector(self) result(omg)
    implicit none
    class(SlidingMesh), intent(in) :: self
    real(kind=RP)                  :: omg(3)

    omg = 0.0_RP
    select case (self % rotationAxis)
    case (rotationAxis_X)
        omg(1) =  self % angularVelocity
    case (rotationAxis_Y)
        omg(2) = -self % angularVelocity
    case (rotationAxis_Z)
        omg(3) =  self % angularVelocity
    end select

end function SlidingMesh_OmegaVector

!  -----------------------------------------------------------------------
!  The two control-file center coordinates live in the plane normal to the
!  rotation axis. NOTE: RotateSlidingRegion currently rotates about the
!  origin (no center offset), so only center = (0,0) is consistent with
!  the geometry update; this mapping exists so the fluxes stay correct the
!  day the rotation itself honors the center.
!  -----------------------------------------------------------------------
pure function SlidingMesh_RotationCenter3D(self) result(xc)
    implicit none
    class(SlidingMesh), intent(in) :: self
    real(kind=RP)                  :: xc(3)

    xc = 0.0_RP
    select case (self % rotationAxis)
    case (rotationAxis_X)      ! rotation plane (y,z)
        xc(2) = self % center(1)
        xc(3) = self % center(2)
    case (rotationAxis_Y)      ! rotation plane (x,z)
        xc(1) = self % center(1)
        xc(3) = self % center(2)
    case (rotationAxis_Z)      ! rotation plane (x,y)
        xc(1) = self % center(1)
        xc(2) = self % center(2)
    end select

end function SlidingMesh_RotationCenter3D

!  -----------------------------------------------------------------------
!  In-plane coordinate indices (horizontal ih, vertical iv) of the rotation
!  plane. All three rotation matrices in RotateSlidingRegion act as the SAME
!  counterclockwise rotation when expressed in these (h,v) coordinates
!  (R_x(+theta) in (y,z), R_y(-theta) in (x,z), R_z(+theta) in (x,y)), so the
!  direction-sensitive sliding connectivity logic written for the y-axis
!  case generalizes to any principal axis by pure index substitution. The
!  control-file center maps as center(1) <-> X(ih), center(2) <-> X(iv),
!  consistent with RotationCenter3D.
!  -----------------------------------------------------------------------
pure subroutine SlidingMesh_PlaneCoords(self, ih, iv)
    implicit none
    class(SlidingMesh), intent(in)  :: self
    integer,            intent(out) :: ih, iv

    select case (self % rotationAxis)
    case (rotationAxis_X)
        ih = 2 ; iv = 3
    case (rotationAxis_Y)
        ih = 1 ; iv = 3
    case (rotationAxis_Z)
        ih = 1 ; iv = 2
    case default
        ih = 1 ; iv = 3
    end select

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
