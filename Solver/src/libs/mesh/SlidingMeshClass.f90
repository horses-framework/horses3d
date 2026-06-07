#include "Includes.h"
MODULE SlidingMeshClass
      use Utilities                       , only: toLower, almostEqual, AlmostEqualRelax
      use SMConstants
      USE MeshTypes
      USE NodeClass
      use FileReaders               , only: ReadControlFile 
      use FileReadingUtilities      , only: getFileName
      use ElementConnectivityDefinitions
      use FTValueDictionaryClass

      IMPLICIT NONE

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
      logical                                :: initialized=.false.
      logical                                :: conforming=.false.

    contains

        procedure :: read_info                              => SlidingMesh_read_info
        procedure :: GetInfo                                => SlidingMesh_GetInfo
        procedure :: Initialize                         => SlidingMesh_Initialize
        procedure :: Destruct                           => SlidingMesh_Destruct
      
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
     
    call self % GetInfo( controlVariables )
     
end subroutine SlidingMesh_read_info

subroutine SlidingMesh_GetInfo( self, controlVariables )
    use FileReadingUtilities
    use ParamfileRegions
    implicit none
        
    class(SlidingMesh), intent(inout) :: self
    class(FTValueDictionary)       :: controlVariables
    
    character(len=LINE_LENGTH) :: in_label, paramFile
    real(kind=RP),allocatable :: centerx
    real(kind=RP),allocatable :: centery
    real(kind=RP),allocatable :: radius
    real(kind=RP),allocatable :: angle 
    integer :: i

!     **********
!     Read block
!     **********

    write(in_label , '(A)') "#define slidingmesh"
    call get_command_argument(1, paramFile) 

    call readValueInRegion( trim( paramFile ), "centerx", centerx, in_label, "#end" )
    call readValueInRegion( trim( paramFile ), "centery", centery, in_label, "#end" )
    call readValueInRegion( trim( paramFile ), "radius", radius, in_label, "#end" )
    call readValueInRegion( trim( paramFile ), "angle", angle, in_label, "#end" )   

    i=0 

    if (allocated(centerx)) then 
        self % center(1) = centerx
        i=i+1
    end if 
    if (allocated(centery)) then 
        self % center(2) = centery
        i=i+1
    end if 
    if (allocated(radius)) then 
        self % radius = radius 
        i=i+1
    end if 
    if (allocated(angle)) then 
        self % theta = angle 
        i=i+1
    end if 

    if (i==4) then 
        self % isConfigured = .true.
    else if (i .NE. 0) then
        write(*,*) 'Error, missing SlidingMesh Argument : Need centerX, centerY, radius and angle'
    end if 
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
   if (allocated(self % mortarNeighborElems)) then
       deallocate(self % mortarNeighborElems)
   end if

   if (allocated(self % slidingMortarElems)) then
       deallocate(self % slidingMortarElems)
   end if

   if (allocated(self % pureSlidingElems)) then
       deallocate(self % pureSlidingElems)
   end if

   if (allocated(self % mortararr1)) then
       deallocate(self % mortararr1)
   end if

   if (allocated(self % mortararr2)) then
       deallocate(self % mortararr2)
   end if

   if (allocated(self % face_nodes)) then
       deallocate(self % face_nodes)
   end if

   if (allocated(self % face_othernodes)) then
       deallocate(self % face_othernodes)
   end if

   if (allocated(self % slidingMortarConnectivity)) then
       deallocate(self % slidingMortarConnectivity)
   end if

   if (allocated(self % neighborConnectivity)) then
       deallocate(self % neighborConnectivity)
   end if

   if (allocated(self % rotmortars)) then
       deallocate(self % rotmortars)
   end if

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
END MODULE SlidingMeshClass