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

      type SlidingRingEntry
      ! communicated 
      integer       :: globID   = 0        ! global element ID
      integer       :: rank     = -1       ! owning rank
      integer       :: side     = 0        ! local face index 1..6 on the owner
      integer       :: rotation = 0        ! face rotation index
      integer       :: dir      = 0        ! azimuthal local direction, see SlidingFaceDir
      real(kind=RP) :: phi      = 0.0_RP   ! azimuth of the face centroid, IN ITS OWN FRAME
      real(kind=RP) :: zeta     = 0.0_RP   ! axial coordinate of the face centroid
      real(kind=RP) :: jacMin   = 0.0_RP   ! min jacobian of the owning element
      integer       :: Nel(2)   = 0        ! polynomial orders of the owning element on this face
      ! assigned locally, never communicated 
      integer       :: Azim_Id  = 0        ! azimuthal index, 1..numElemsPerLayer
      integer       :: Ax_Id    = 0        ! axial layer index, 1..numLayers
      integer       :: localEID = 0        ! local element ID, valid only on the owner
   end type SlidingRingEntry
   

   integer, parameter :: SM_RING_NINT  = 7      
   integer, parameter :: SM_RING_NREAL = 3

    type SlidingMesh
        logical                                :: active = .false.
        logical                                :: isConfigured = .false.
        integer                                :: numPureSlidingElems = 0
        integer, allocatable                   :: slidingElems(:)       ! All sliding elements
        integer, allocatable                   :: slidingFaces(:)       ! All sliding faces
        integer, allocatable                   :: mortarNeighborElems(:)! associated non-sliding neighbors across the interface (require mortars)
        integer, allocatable                   :: slidingMortarElems(:) ! sliding elements at the interface (require mortars)
        integer, allocatable                   :: pureSlidingElems(:)   ! sliding elements fully inside the region (no mortar)
        integer, allocatable                   :: slidingMortarConnectivity(:,:)
        !   slidingMortarConnectivity(i,:) :
        !   One row per static ring row this rank owns, in myMasterRows order.
        !   A zero element means the other side is on another rank.
        !
        !     (1),(4),(7) : element, local face and rotation of the rotating side of Mortarpos 1
        !     (2),(5),(8) : the same for Mortarpos 0
        !     (3),(6),(9) : element, local face and rotation of the static side, the master, shared by both families
        integer                                :: numBFacePoints
        integer                                :: numSlidingElements
        integer                                :: numSlidingInterfaceElements 
        integer                                :: currentSectorID = -1
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
        integer                                :: numDuplicatedNodes = 0
        integer                                :: numLayers = 0
        !   global interface directory, replicated on every rank 
        type(SlidingRingEntry), allocatable    :: rotorRing(:)    ! allocated by AllgatherRing
        type(SlidingRingEntry), allocatable    :: statorRing(:)   
        integer                                :: nRingGlobal = 0 ! numElemsPerLayer * numLayers
        !   local entries, produced by CollectLocalRingEntries 
        type(SlidingRingEntry), allocatable    :: localRotorEnt(:)
        type(SlidingRingEntry), allocatable    :: localStatorEnt(:)
        integer                                :: nLocalRotorEnt  = 0
        integer                                :: nLocalStatorEnt = 0
        !   global -> local resolution, valid for ring elements only 
        integer, allocatable                   :: rotorLocalElem(:)   ! nRingGlobal, 0 when not mine
        integer, allocatable                   :: statorLocalElem(:)  ! nRingGlobal, 0 when not mine
        integer, allocatable                   :: myMasterRows(:)     ! stator rows I own: I build the mortars
        integer, allocatable                   :: mySlaveRows(:)      ! stator rows where a rotor slave is mine
        integer                                :: nMasterRows = 0
        integer                                :: nSlaveRows  = 0
        integer                                :: nSlaveRowsMax = 0   ! max over all sectors, set once at setup
        !   global geometric anchors 
        real(kind=RP)                          :: Rint     = 0.0_RP   ! interface radius, Allreduce(max)
        real(kind=RP)                          :: zetaMin  = 0.0_RP   ! axial extent of the rotating region
        real(kind=RP)                          :: zetaMax  = 0.0_RP
        integer                                :: NfInterface(2) = 0   ! uniform polynomial order on the interface
        integer                                :: nodeGlobIDOffset = 0 ! global node count before the sliding split
        !  local counts, NOT to be confused with the global ones 
        integer                                :: nLocalSplitFaces = 0   ! of those, the ones needing a local split
        logical                                :: needsLocalSplit  = .false.
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
        procedure :: RingRow                            => SlidingMesh_RingRow
        procedure :: RingAzimId                         => SlidingMesh_RingAzimId
        procedure :: RingAxId                           => SlidingMesh_RingAxId
    end type

!     ========
      CONTAINS
!     ========
!
!////////////////////////////////////////////////////////////////////////////////
!
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
!
!////////////////////////////////////////////////////////////////////////////////
!
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
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SlidingMesh_Initialize
!
!  Allocates the connectivity tables once the interface has been sized, and
!  zeroes them on the first pass only. Allocation is guarded, so calling it
!  again on a later step is harmless.
! ---------------------------------------------------------------------------
subroutine SlidingMesh_Initialize(self, numLocalInterfaceElements, numLocalSlidingElements, numBFacePoints)
    use Physics
    use PartitionedMeshClass
    use MPI_Process_Info

    implicit none

    class(SlidingMesh), intent(inout)        :: self
    integer, intent(in)                      :: numLocalInterfaceElements
    integer, intent(in)                      :: numLocalSlidingElements
    integer, intent(in)                      :: numBFacePoints
   !========================
   ! Allocation
   !========================
    if (.not. allocated(self % mortarNeighborElems))             allocate(self % mortarNeighborElems(numLocalInterfaceElements))
    if (.not. allocated(self % slidingMortarElems))              allocate(self % slidingMortarElems(numLocalInterfaceElements))
    if (.not. allocated(self % pureSlidingElems))                allocate(self % pureSlidingElems(numLocalSlidingElements))
    if (.not. allocated(self % slidingMortarConnectivity))       allocate(self % slidingMortarConnectivity(numLocalInterfaceElements, 9))
    
    !========================
    ! Initialization (first pass only)
    !========================
    if (.not. self % active) then
       self % mortarNeighborElems       = 0
       self % slidingMortarElems        = 0
       self % pureSlidingElems          = 0
       self % slidingMortarConnectivity = 0
       self % localAngle                = 0
       self % currentSectorID           = -1
       self % numSlidingInterfaceElements = numLocalInterfaceElements
       self % numSlidingElements          = numLocalSlidingElements
       self % numBFacePoints              = numBFacePoints
    end if 

end subroutine SlidingMesh_Initialize
!
!////////////////////////////////////////////////////////////////////////////////
!
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
    safedeallocate(self % slidingMortarConnectivity)
    safedeallocate(self % rotorRing)
    safedeallocate(self % statorRing)
    safedeallocate(self % rotorLocalElem)
    safedeallocate(self % statorLocalElem)
    safedeallocate(self % myMasterRows)
    safedeallocate(self % mySlaveRows)
    safedeallocate(self % localRotorEnt)
    safedeallocate(self % localStatorEnt)
    safedeallocate(self % slidingElems)
    safedeallocate(self % slidingFaces)

    ! ---------------------------------------------------------------------
    ! Reset counters and state variables
    ! ---------------------------------------------------------------------

    self % numSlidingInterfaceElements = 0
    self % numSlidingElements          = 0
    self % numBFacePoints              = 0

    self % omega        = 0.0_RP
    self % theta        = 0.0_RP

    self % active       = .false.

    self % nRingGlobal      = 0
    self % nMasterRows      = 0
    self % nSlaveRows       = 0
    self % nLocalRotorEnt   = 0
    self % nLocalStatorEnt  = 0
    self % Rint             = 0.0_RP
    self % zetaMin          = 0.0_RP
    self % zetaMax          = 0.0_RP
    self % NfInterface      = 0
    self % nodeGlobIDOffset = 0
    self % numElemsPerLayer = 0
    self % numLayers        = 0
    self % numDuplicatedNodes = 0
    self % nLocalSplitFaces = 0
    self % needsLocalSplit  = .false.
    self % currentSectorID  = -1
    self % nSlaveRowsMax    = 0

end subroutine SlidingMesh_Destruct
!
!////////////////////////////////////////////////////////////////////////////////
!
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
!
!////////////////////////////////////////////////////////////////////////////////
!
!  Azimuthal angle of x about the axis, in (-pi, pi].
pure function SlidingMesh_phi(self, x) result(p)
   class(SlidingMesh), intent(in) :: self
   real(kind=RP),      intent(in) :: x(3)
   real(kind=RP) :: p

   p = atan2(x(self % iR2) - self % center(2), x(self % iR1) - self % center(1))

end function
!
!////////////////////////////////////////////////////////////////////////////////
!
!  Coordinate of x along the rotation axis.
pure function SlidingMesh_zeta(self, x) result(z)
   class(SlidingMesh), intent(in) :: self
   real(kind=RP),      intent(in) :: x(3)
   real(kind=RP) :: z

   z = x(self % iAx)

end function
!
!////////////////////////////////////////////////////////////////////////////////
!
!  Distance from x to the axis.
pure function SlidingMesh_rad(self, x) result(r)
   class(SlidingMesh), intent(in) :: self
   real(kind=RP),      intent(in) :: x(3)
   real(kind=RP) :: r

   r = sqrt((x(self % iR1) - self % center(1))**2 + (x(self % iR2) - self % center(2))**2)

end function
!
!////////////////////////////////////////////////////////////////////////////////
!
!  Rotates x about the axis by theta.
pure subroutine SlidingMesh_RotatePoint(self, theta, x)
   class(SlidingMesh), intent(in)    :: self
   real(kind=RP),      intent(in)    :: theta
   real(kind=RP),      intent(inout) :: x(3)

   real(kind=RP) :: u, v, c, s

   c = cos(theta) 
   s = sin(theta)
   u = x(self % iR1) - self % center(1)
   v = x(self % iR2) - self % center(2)

   x(self % iR1) = self % center(1) + c*u - s*v
   x(self % iR2) = self % center(2) + s*u + c*v

end subroutine
!
!////////////////////////////////////////////////////////////////////////////////
!
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
!
!////////////////////////////////////////////////////////////////////////////////
!
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
!
!////////////////////////////////////////////////////////////////////////////////
!
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
!
!////////////////////////////////////////////////////////////////////////////////
!
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
!
!////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Row of the ring array for azimuthal index Azim_Id of axial layer Ax_Id.
!  Azim_Id is taken modulo the number of elements per layer, so callers can
!  pass Azim_Id-1 or Azim_Id-sectorID without wrapping by hand. This is the
!  only place where the closure of the ring is expressed.
!  -----------------------------------------------------------------------
pure integer function SlidingMesh_RingRow(self, Ax_Id, Azim_Id) result(row)
   class(SlidingMesh), intent(in) :: self
   integer,            intent(in) :: Ax_Id, Azim_Id

   row = (Ax_Id - 1) * self % numElemsPerLayer &
       + modulo(Azim_Id - 1, self % numElemsPerLayer) + 1

end function SlidingMesh_RingRow
!
!////////////////////////////////////////////////////////////////////////////////
!
pure integer function SlidingMesh_RingAzimId(self, row) result(Azim_Id)
   class(SlidingMesh), intent(in) :: self
   integer,            intent(in) :: row

   Azim_Id = modulo(row - 1, self % numElemsPerLayer) + 1

end function SlidingMesh_RingAzimId
!
!////////////////////////////////////////////////////////////////////////////////
!
pure integer function SlidingMesh_RingAxId(self, row) result(Ax_Id)
   class(SlidingMesh), intent(in) :: self
   integer,            intent(in) :: row

   Ax_Id = (row - 1) / self % numElemsPerLayer + 1

end function SlidingMesh_RingAxId
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SlidingSlaveRows
!
!  The two rotating interface elements overlapping static element j indexed by
!  Mortarpos:
!
!      Azim_Id_rotor = Azim_Id_stator - sectorID - Mortarpos
!
!  After turning by omega the rotating element of index k occupies the sector of
!  effective index k + sectorID, so each static sector is straddled by exactly
!  two of them. Holds for negative sectorID.
!
!  rows(0) covers the static face as localAngle -> 0, rows(1) degenerates to
!  zero scale there.
! ---------------------------------------------------------------------------
subroutine SlidingSlaveRows(SM, j, sectorID, rows)
    implicit none
    ! =========================
    ! Arguments
    ! =========================
    type(SlidingMesh), intent(in)  :: SM
    integer,           intent(in)  :: j, sectorID
    integer,           intent(out) :: rows(0:1)          ! indexed BY Mortarpos
    ! =========================
    ! Local variables
    ! =========================
    integer :: Ax_Id, Azim_Id, Mortarpos
 
    Ax_Id   = SM % RingAxId(j)
    Azim_Id = SM % RingAzimId(j) - sectorID
 
    do Mortarpos = 0, 1
       rows(Mortarpos) = SM % RingRow(Ax_Id, Azim_Id - Mortarpos)
    end do
 
 end subroutine SlidingSlaveRows
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SlidingMaxSlaveRows
!
!  Largest number of static interface elements this rank is the slave of, over every
!  sector, so the mortar tables are sized once and never reallocated.
!
!  Sweeping the sectors gives the exact value: the trivial bound 2m is loose,
!  a contiguous arc of m rotor elements touching only m + 1 static rows.
!
!  Local, constant for the run. Call once, after the rings are built.
! ---------------------------------------------------------------------------
integer function SlidingMaxSlaveRows(SM) result(nMax)
    implicit none
    ! =========================
    ! Arguments
    ! =========================
    type(SlidingMesh), intent(in) :: SM
    ! =========================
    ! Local variables
    ! =========================
    integer :: s, j, n, rows(0:1)

    nMax = 0

    do s = 0, SM % numElemsPerLayer - 1
    n = 0
    do j = 1, SM % nRingGlobal
        call SlidingSlaveRows(SM, j, s, rows)
        if (SM % rotorLocalElem(rows(0)) > 0 .or. &
            SM % rotorLocalElem(rows(1)) > 0) n = n + 1
    end do
    nMax = max(nMax, n)
    end do

end function SlidingMaxSlaveRows
 
END MODULE SlidingMeshClass
