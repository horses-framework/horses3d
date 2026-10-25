!
!////////////////////////////////////////////////////////////////////////////////
!
!     SlidingMeshProcedures
!     ---------------------
!
!     Rotates part of the mesh while the solver runs. The rotating region is a
!     cylinder given by a radius, a centre and an axis, all read from the control
!     file into mesh % SlidingMesh. Elements inside it turn, the rest stay, and the
!     two sides are joined by mortar faces rebuilt at every stage.
!
!     The same code runs in sequential and in MPI. Two partitions are supported:
!        strategy 1   both sides of an axial layer on the same rank, every
!                     mortar is local
!        strategy 2   the two sides always on different ranks, each mortar is
!                     split between two ranks (see SlidingMPIMortars)
!
!     Who calls whom
!     --------------
!     AdvanceSlidingMesh is the only entry point. It runs at every Runge-Kutta
!     stage and drives everything below.
!
!        AdvanceSlidingMesh
!         |
!         +-- IdentifySlidingRegion              first call only
!         |     |  which elements rotate, plus the interface radius and axial band
!         |     +-- SlidingComputeGeometricAnchors
!         |
!         +-- SlidingBuildGlobalRings            first call only, collective
!         |     |  the interface directory, identical on every rank: one entry per
!         |     |  interface face on each side, indexed by (Azim_Id, Ax_Id)
!         |     +-- CollectLocalRingEntries
!         |     |     +-- IsInterfaceFace
!         |     |     +-- SlidingFaceDir
!         |     |     +-- SlidingFaceQuadrants   orientation tag of each face
!         |     +-- AllgatherRing
!         |     +-- BuildAzimuthalConnectivity
!         |     |     +-- SortRingEntries
!         |     +-- SlidingMaxSlaveRows
!         |
!         +-- SlidingConstructMPIMortars         first call only (SlidingMPIMortars)
!         |     buffers of the split mortars; finds nothing to do outside strategy 2
!         |
!         +-- SplitInterfaceNodes                first call only
!         |     gives the rotating side its own copy of the interface nodes
!         |
!         +-- GetSlidingTopologyState            every call
!         |     how far we have turned: which sector, how far inside it
!         |
!         +-- InitializeSlidingConnectivity      every call
!         |     +-- RotateSlidingRegion          rotates the region, offsets and scales
!         |     +-- SlidingRotateNodes           first call only, moves mesh % nodes
!         |
!         +-- UpdateSlidingMortarsConnectivity   when the sector changes
!         |     pairs each static interface face with the two rotating faces
!         |     overlapping it
!         |
!         +-- SlidingSplitInterfaceFaces         first call only, strategy 1 and sequential
!         +-- SlidingBuildGeometryLists          first call only
!         +-- SlidingPruneMPIFaces               first call only
!         |     removes the interface faces from the ordinary MPI exchanges
!         |
!         +-- ConstructSlidingMortars            every call, strategy 1 and sequential
!         |     both sides of every mortar are on this rank
!         +-- SlidingConstructSplitMortars       every call, strategy 2
!         |     full mortar on the static rank, light copy on the rotating rank
!         |     (both use SetSlidingMortarPlacement: SlidingPairRotation,
!         |      QuadAzimSign, QuadAzimDir)
!         |
!         +-- SlidingBuildCommLists              when the sector changes (SlidingMPIMortars)
!         +-- SlidingCheckCommAlignment          when the sector changes (SlidingMPIMortars)
!
!/////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"

module SlidingMeshProcedures
   use Utilities                       , only: toLower, almostEqual, AlmostEqualRelax
   use HexMeshClass
   use SlidingMeshClass
   use SMConstants
   use MeshTypes
   use NodeClass
   use ElementClass
   use FaceClass
   use ElementConnectivityDefinitions
   use MPI_Process_Info
   use MPI_Face_Class
   use SlidingMPIMortars
   use Headers
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none

   private

   public :: AdvanceSlidingMesh
   public :: IdentifySlidingRegion
   public :: InitializeSlidingConnectivity
   public :: RotateSlidingRegion
   public :: ConstructSlidingMortars
   public :: PrintMortarConnectivity
   public :: UpdateSlidingMortarsConnectivity

contains
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  AdvanceSlidingMesh
!
!  One Runge-Kutta stage. On the first call it builds the region and the pairing,
!  on every call it advances omega and the geometry by the same increment
!  (the angle argument when present, else SM % theta), reads off the new 
!  topology state, moves the nodes, and rebuilds the mortars.
!  Sets SM % active on the way out so the setup work happens once.
! ---------------------------------------------------------------------------
subroutine AdvanceSlidingMesh(mesh, numBFacePoints, nodes, useMPI, angle, rotateEntireMesh)
   use Physics
   use PartitionedMeshClass

   implicit none 
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(inout)        :: mesh
   integer, intent(inout)               :: numBFacePoints
   integer                              :: nodes
   logical, intent(in)                  :: useMPI 
   real(kind=RP), intent(in), optional  :: angle
   logical, intent(in), optional        :: rotateEntireMesh
   ! =========================
   ! Local variables
   ! =========================
   integer :: numSlidingElements
   integer :: l, eID, fID, r
   integer :: sectorID 
   logical :: isConforming 
   real(kind=RP) :: offsetParams(4)
   real(kind=RP) :: scaleParams(4)
   real(kind=RP) :: dTheta                          ! rotation increment for this call: the optional `angle`
                                                    ! argument when present (RK stage fraction), else theta.
   integer :: no_interior_faces, no_boundary_faces, no_mpi_faces
   integer :: no_sequential_elems, no_mpi_elems
   integer :: aux_array(1:3)
   logical :: sectorChanged 

   sectorChanged = .false.
   associate(SM => mesh % SlidingMesh)

      ! Setup phase, first call only
      if (.not. SM % active) then
         if (MPI_Process % isRoot) then
            write(STD_OUT,'(/)')
            call Section_Header("Sliding mesh configuration")
            write(STD_OUT,'(/)')
   
            write(STD_OUT,'(30X,A,A30,"(",1pE12.5,", ",1pE12.5,")")') "->", "Center:", SM % Center(1), SM % Center(2)
            write(STD_OUT,'(30X,A,A30,1pG10.3)') "->", "Radius: ", SM % Radius  
            write(STD_OUT,'(30X,A,A30,1pG10.3)') "->", "Angle of rotation: ", SM % theta  
            write(STD_OUT,'(/)')
         end if
      end if

      if (.not. SM % active) then
         call IdentifySlidingRegion(mesh, numSlidingElements)
         call SlidingBuildGlobalRings(mesh)
         call SM % Initialize(max(SM % nMasterRows, SM % nSlaveRowsMax), numSlidingElements, numBFacePoints)
         ! Second MPI face set, for the mortars split across ranks. Called on every
         ! rank: under strategy 1 it finds no partner and allocates nothing.
         call SlidingConstructMPIMortars(mesh)
      end if

      ! Give the rotating side its own interface nodes
      if (.not. SM % active) then
         call SplitInterfaceNodes(mesh)

         ! the local sweep is over
         safedeallocate(SM % localRotorEnt)  
         SM % nLocalRotorEnt = 0
         safedeallocate(SM % localStatorEnt) 
         SM % nLocalStatorEnt = 0
      end if

      ! Optional angle override
      if (present(angle)) then 
         dTheta = angle
      else 
         dTheta = SM % theta
      end if 
      SM % omega = SM % omega + dTheta
      
      ! Core mesh modification
      call GetSlidingTopologyState(SM % omega, SM % numElemsPerLayer, sectorID, isConforming, SM % localAngle)

      call InitializeSlidingConnectivity(mesh, dTheta, offsetParams, scaleParams)

      if (sectorID .NE. SM % currentSectorID) then
         call UpdateSlidingMortarsConnectivity(mesh, sectorID)
         SM % currentSectorID = sectorID
         sectorChanged = .true.
      end if

      if (.not. present(rotateEntireMesh)) then
         do r = 1, SM % nRingGlobal

            if (SM % rotorLocalElem(r) > 0) then
               eID = SM % rotorLocalElem(r)
               if (any(mesh % elements(eID) % MortarFaces == 1)) then
                  mesh % elements(eID) % MortarFaces = 0
               end if
            end if

            if (SM % statorLocalElem(r) > 0) then
               eID = SM % statorLocalElem(r)
               if (any(mesh % elements(eID) % MortarFaces == 1)) then
                  mesh % elements(eID) % MortarFaces = 0
               end if
            end if

         end do
      end if
      
      !  Face splitting: Split only what the partition has not split already
      if (.not. SM % active .and. SM % needsLocalSplit) then
         call SlidingSplitInterfaceFaces(mesh)
      end if

      if (.not. SM % active) then
         call SlidingBuildGeometryLists(mesh)
      end if

      !  Tag the interface faces 
      !  SlidingPruneMPIFaces then removes them from the ordinary MPI face lists
      if (.not. present(rotateEntireMesh)) then
         do r = 1, SM % nRingGlobal

            if (SM % rotorLocalElem(r) > 0) then
               fID = mesh % elements(SM % rotorLocalElem(r)) % faceIDs(SM % rotorRing(r) % side)
               mesh % faces(fID) % MortarType = MORTAR_SLIDING
               mesh % faces(fID) % faceType   = HMESH_INTERIOR
            end if

            if (SM % statorLocalElem(r) > 0) then
               fID = mesh % elements(SM % statorLocalElem(r)) % faceIDs(SM % statorRing(r) % side)
               mesh % faces(fID) % MortarType = MORTAR_SLIDING
               mesh % faces(fID) % faceType   = HMESH_INTERIOR
            end if

         end do

         !  Catch-all for the faces left undefined by the split
         do l = 1, size(mesh % faces)
            if (mesh % faces(l) % faceType == HMESH_UNDEFINED) then
               mesh % faces(l) % faceType = HMESH_INTERIOR
            end if
         end do
      end if

      ! Connectivity
      if (.not. SM % active .and. SM % needsLocalSplit) then
         call mesh % SetConnectivitiesAndLinkFaces(nodes)
      end if

      if (.not. SM % active) then
         call SlidingPruneMPIFaces(mesh)
      end if

      if (.not. SM % active) then
         do l = 1, mesh % no_of_elements
            call mesh % elements(l) % geom % destruct
         end do
         do l = 1, size(mesh % faces)
            call mesh % faces(l) % geom % destruct
         end do
         call mesh % ConstructGeometry()
      else
         do l = 1, size(SM % slidingElems)
            call mesh % elements(SM % slidingElems(l)) % geom % destruct
         end do
         do l = 1, size(SM % slidingFaces)
            call mesh % faces(SM % slidingFaces(l)) % geom % destruct
         end do
         call mesh % ConstructGeometry(facesList = SM % slidingFaces, elementList = SM % slidingElems)
      end if
      
      ! Mortar geometry cleanup 
      if (allocated(mesh % mortar_faces)) then 
         do l = 1, size(mesh % mortar_faces)
            call mesh % mortar_faces(l) % geom % destruct
         end do 
      end if

      ! Create three arrays that contain the fIDs of interior, mpi, and boundary faces
      if (.not. SM  % active) then 
         mesh % no_of_faces = size(mesh % faces)

         no_interior_faces = 0
         no_boundary_faces = 0
         no_mpi_faces      = 0

         aux_array = 0
         do fID = 1, mesh % no_of_faces
            aux_array(mesh % faces(fID) % faceType) = aux_array(mesh % faces(fID) % faceType) + 1
         end do

         no_interior_faces = aux_array(HMESH_INTERIOR)
         no_mpi_faces      = aux_array(HMESH_MPI)
         no_boundary_faces = aux_array(HMESH_BOUNDARY)

         safedeallocate(mesh % faces_interior)
         safedeallocate(mesh % faces_mpi)
         safedeallocate(mesh % faces_boundary)

         allocate(mesh % faces_interior(no_interior_faces))
         allocate(mesh % faces_mpi(no_mpi_faces))
         allocate(mesh % faces_boundary(no_boundary_faces))

         no_interior_faces = 0
         no_boundary_faces = 0
         no_mpi_faces      = 0

         do fID = 1, mesh % no_of_faces
            select case (mesh % faces(fID) % faceType)
            case(HMESH_INTERIOR)
               no_interior_faces = no_interior_faces + 1
               mesh % faces_interior(no_interior_faces) = fID

            case(HMESH_MPI)
               no_mpi_faces = no_mpi_faces + 1
               mesh % faces_mpi(no_mpi_faces) = fID

            case(HMESH_BOUNDARY)
               no_boundary_faces = no_boundary_faces + 1
               mesh % faces_boundary(no_boundary_faces) = fID

            end select
         end do

         no_sequential_elems = 0
         no_mpi_elems = 0

         do eID = 1, mesh % no_of_elements
            if (mesh % elements(eID) % hasSharedFaces) then
               no_mpi_elems = no_mpi_elems + 1 
            else
               no_sequential_elems = no_sequential_elems + 1 
            end if
         end do

         safedeallocate(mesh % elements_sequential)
         safedeallocate(mesh % elements_mpi)
         allocate(mesh % elements_sequential(no_sequential_elems))
         allocate(mesh % elements_mpi(no_mpi_elems))
         
         no_sequential_elems = 0
         no_mpi_elems = 0

         do eID = 1, mesh % no_of_elements
            if (mesh % elements(eID) % hasSharedFaces) then
               no_mpi_elems = no_mpi_elems + 1 
               mesh % elements_mpi(no_mpi_elems) = eID
            else
               no_sequential_elems = no_sequential_elems + 1 
               mesh % elements_sequential(no_sequential_elems) = eID
            end if
         end do

      end if 

      if (.not. present(rotateEntireMesh)) then

         if (.not. allocated(mesh % mortar_faces)) then
            allocate(mesh % mortar_faces(MortarSlotCount(SM)))
!           Empty slots are skipped by elementIDs <= 0; InitializeSlidingConnectivity
!           only resets them from the second call on
            do l = 1, size(mesh % mortar_faces)
               mesh % mortar_faces(l) % ID          = -1
               mesh % mortar_faces(l) % elementIDs  = -1
               mesh % mortar_faces(l) % elementSide = -1
            end do
         end if

!        Strategy 1 and sequential: every pair is local. Strategy 2: every pair
!        is split across ranks. The partition check in SlidingBuildGlobalRings
!        guarantees there is nothing in between.
         if (SM % needsLocalSplit) then
            call ConstructSlidingMortars(mesh, nodes, SM % numSlidingInterfaceElements, &
                           SM % slidingMortarConnectivity, offsetParams, scaleParams)
         else
            call SlidingConstructSplitMortars(mesh, nodes, offsetParams, scaleParams)
         end if

!        The links depend on the current pair only, and read the mortar sizes:
!        rebuilt after the construction, when the sector has changed
         if (sectorChanged) then
            call SlidingBuildCommLists(mesh, SM % currentSectorID)
            call SlidingCheckCommAlignment(mesh)   
         end if
      end if
      
      SM  % active = .true.
      mesh % slidingflux = .false.

   end associate

end subroutine AdvanceSlidingMesh
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  IdentifySlidingRegion
!
!  Splits the local elements into rotating and static, then computes the global
!  geometric anchors. It stops there on purpose: finding the interface faces
!  needs the anchors produced at the end of this routine, and is done once, for
!  both sides at once, in CollectLocalRingEntries.
!
!  First step only. COLLECTIVE.
! ---------------------------------------------------------------------------
subroutine IdentifySlidingRegion(mesh, numLocalSlidingElements)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(inout) :: mesh
   integer,        intent(out)   :: numLocalSlidingElements
   ! =========================
   ! Local variables
   ! =========================
   integer :: i, j, numNodesInsideRadius

   numLocalSlidingElements = 0

   do i = 1, size(mesh % elements)
      numNodesInsideRadius = 0
      do j = 1, 8
         if (mesh % SlidingMesh % rad(mesh % nodes(mesh % elements(i) % nodeIDs(j)) % X) &
              <= mesh % SlidingMesh % radius) then
            numNodesInsideRadius = numNodesInsideRadius + 1
         end if
      end do

      mesh % elements(i) % sliding = (numNodesInsideRadius == 8)
      if (mesh % elements(i) % sliding) then
         numLocalSlidingElements = numLocalSlidingElements+1
      end if 
   end do

   call SlidingComputeGeometricAnchors(mesh)

end subroutine IdentifySlidingRegion
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SlidingComputeGeometricAnchors
!
!  Interface radius and axial extent of the rotating region.
!
!  Every rotating element satisfies rad <= radius and the outermost layer
!  touches the interface, so the largest nodal radius found on a rotating
!  element is the interface radius. Reducing over all ranks makes it
!  independent of the partition, which is what lets IsInterfaceFace be local.
!
!  COLLECTIVE. Every rank must reach it, including those owning no rotating
!  element.
! ---------------------------------------------------------------------------
subroutine SlidingComputeGeometricAnchors(mesh)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(inout) :: mesh
   ! =========================
   ! Local variables
   ! =========================
   integer       :: eID, j, ierr
   real(kind=RP) :: xn(3)
   real(kind=RP) :: locMax(2), gloMax(2)     
   real(kind=RP) :: locMin, gloMin          

   associate(SM => mesh % SlidingMesh)

      locMax = -huge(1.0_RP)
      locMin = huge(1.0_RP)

      do eID = 1, mesh % no_of_elements
         if (.not. mesh % elements(eID) % sliding) cycle
         do j = 1, 8
            xn = mesh % nodes(mesh % elements(eID) % nodeIDs(j)) % X
            locMax(1) = max(locMax(1), SM % rad(xn))
            locMax(2) = max(locMax(2), SM % zeta(xn))
            locMin = min(locMin   , SM % zeta(xn))
         end do
      end do

      gloMax = locMax
      gloMin = locMin
#ifdef _HAS_MPI_
      if (MPI_Process % doMPIAction) then
         call MPI_Allreduce(locMax, gloMax, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
         call MPI_Allreduce(locMin, gloMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
      end if
#endif

      if (gloMax(1) <= 0.0_RP) then
         write(STD_OUT,*) 'sliding: no rotating element found on any rank, check radius and center'
         error stop
      end if

      SM % Rint = gloMax(1)
      SM % zetaMax = gloMax(2)
      SM % zetaMin = gloMin

   end associate
end subroutine SlidingComputeGeometricAnchors
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  IsInterfaceFace
!
!  True when local face j of element eID lies on the sliding interface: its
!  four nodes sit at Rint and within the axial band of the rotating region.
!
!  Side agnostic. Applied to a rotating element it returns its interface face,
!  applied to a static one it returns its own. Uses nothing but local node
!  coordinates and two global scalars, so it gives the same answer on every
!  rank whatever the partition.
! ---------------------------------------------------------------------------
logical function IsInterfaceFace(mesh, eID, j) result(onInterface)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(in) :: mesh
   integer,        intent(in) :: eID, j
   ! =========================
   ! Local variables
   ! =========================
   real(kind=RP), parameter :: tolR = 1.0e-6_RP    
   real(kind=RP), parameter :: tolZ = 1.0e-6_RP
   integer       :: fID, c
   real(kind=RP) :: xn(3), r, z

   associate(SM => mesh % SlidingMesh)

      onInterface = .false.
      fID = mesh % elements(eID) % faceIDs(j)
      if (fID <= 0) return

      do c = 1, 4
         xn = mesh % nodes(mesh % faces(fID) % nodeIDs(c)) % X
         r = SM % rad(xn)
         z = SM % zeta(xn)
         if (abs(r - SM % Rint) > tolR * SM % Rint) return
         if (z < SM % zetaMin - tolZ .or. z > SM % zetaMax + tolZ) return
      end do

      onInterface = .true.

   end associate
   
end function IsInterfaceFace
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SlidingBuildGlobalRings
!
!  Builds the two interface directories, one per side, identical on every
!  rank, and everything derived from them: the global counts, the global to
!  local tables, and the rows this rank owns or is slave of.
!
!  Row i of one ring faces row i of the other, because rotating azimuths are
!  stored in the rotating frame. The alignment check only confirms it.
!
!  Stops if the two sides are on different ranks and the polynomial order is
!  not the same on every interface face.
!
!  Stops if the partition is neither strategy 1 nor strategy 2 on this rank:
!  checked here on the whole ring, so it holds at every angle.
!
!  Called once. COLLECTIVE, on every rank.
! ---------------------------------------------------------------------------
subroutine SlidingBuildGlobalRings(mesh)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(inout) :: mesh
   ! =========================
   ! Local variables
   ! =========================
   integer :: nRingRotor,  nAzimRotor,  nAxRotor
   integer :: nRingStator, nAzimStator, nAxStator
   integer :: i, k, n0
   integer :: nWhole, nMixed, nHeld
   integer, allocatable :: nStaLayer(:), nRotLayer(:)
   real(kind=RP) :: tolz
   intrinsic :: count

   associate(SM => mesh % SlidingMesh)

      tolz = 1.0e-3_RP * SM % radius

      ! The local entries stay in SM: SplitInterfaceNodes reads them later, from another branch of the call tree.
      call CollectLocalRingEntries(mesh, SM % localRotorEnt , SM % nLocalRotorEnt, &
                                         SM % localStatorEnt, SM % nLocalStatorEnt)

      call AllgatherRing(SM % localRotorEnt , SM % nLocalRotorEnt , SM % rotorRing , nRingRotor)
      call AllgatherRing(SM % localStatorEnt, SM % nLocalStatorEnt, SM % statorRing, nRingStator)

      call BuildAzimuthalConnectivity(SM % rotorRing , nRingRotor , tolz, nAzimRotor , nAxRotor , 'rotor')
      call BuildAzimuthalConnectivity(SM % statorRing, nRingStator, tolz, nAzimStator, nAxStator, 'stator')

      if (nRingRotor /= nRingStator .or. nAzimRotor /= nAzimStator .or. nAxRotor /= nAxStator) then
         write(STD_OUT,*) 'sliding: the two sides of the interface do not match,', &
                          ' rotor', nRingRotor , nAzimRotor , nAxRotor , &
                          ' stator', nRingStator, nAzimStator, nAxStator
         error stop
      end if

      SM % numElemsPerLayer = nAzimRotor
      SM % numLayers = nAxRotor
      SM % nRingGlobal = nRingRotor

      ! global to local resolution 
      allocate(SM % rotorLocalElem (SM % nRingGlobal)) 
      SM % rotorLocalElem = 0

      allocate(SM % statorLocalElem(SM % nRingGlobal)) 
      SM % statorLocalElem = 0

      do i = 1, SM % nRingGlobal
         if (SM % rotorRing(i) % rank == MPI_Process % rank) then
            do k = 1, SM % nLocalRotorEnt
               if (SM % localRotorEnt(k) % globID == SM % rotorRing(i) % globID) then
                  SM % rotorLocalElem(i) = SM % localRotorEnt(k) % localEID 
                  exit
               end if
            end do
         end if

         if (SM % statorRing(i) % rank == MPI_Process % rank) then
            do k = 1, SM % nLocalStatorEnt
               if (SM % localStatorEnt(k) % globID == SM % statorRing(i) % globID) then
                  SM % statorLocalElem(i) = SM % localStatorEnt(k) % localEID 
                  exit
               end if
            end do
         end if
      end do

      ! static interface elements I own: build their mortars
      SM % nMasterRows = count(SM % statorLocalElem > 0)
      allocate(SM % myMasterRows(max(SM % nMasterRows,1)))
      k = 0
      do i = 1, SM % nRingGlobal
         if (SM % statorLocalElem(i) > 0) then
            k = k+1
            SM % myMasterRows(k) = i
         end if
      end do

      ! static interface elements I am the slave of: bounded once, over every sector 
      SM % nSlaveRowsMax = SlidingMaxSlaveRows(SM)
      allocate(SM % mySlaveRows(max(SM % nSlaveRowsMax,1)))
      SM % mySlaveRows = 0
      SM % nSlaveRows = 0

      ! the two rings must be aligned index by index 
      do i = 1, SM % nRingGlobal
         if (abs(WrapPi(SM % statorRing(i) % phi - SM % rotorRing(i) % phi)) &
              > 1.0e-8_RP) then
            write(STD_OUT,*) 'sliding: rings not aligned at row', i, &
                  SM % statorRing(i) % phi, SM % rotorRing(i) % phi
            error stop
         end if
      end do

      ! sides on different ranks: every interface face must have the same order
      if (any(SM % statorRing(1:SM % nRingGlobal) % rank /= &
              SM % rotorRing (1:SM % nRingGlobal) % rank)) then
         n0 = SM % statorRing(1) % Nel(1)
         if (any(SM % statorRing(1:SM % nRingGlobal) % Nel(1) /= n0) .or. &
             any(SM % statorRing(1:SM % nRingGlobal) % Nel(2) /= n0) .or. &
             any(SM % rotorRing (1:SM % nRingGlobal) % Nel(1) /= n0) .or. &
             any(SM % rotorRing (1:SM % nRingGlobal) % Nel(2) /= n0)) then
            if (MPI_Process % isRoot) then
               write(STD_OUT,*) 'sliding: static and rotating sides on different ranks require'
               write(STD_OUT,*) 'a uniform, isotropic polynomial order on the interface faces'
            end if
            error stop
         end if
      end if

      ! partition: per layer, this rank holds the whole layer (strategy 1)
      ! or one side only (strategy 2), and never both kinds of layer
      allocate(nStaLayer(SM % numLayers), nRotLayer(SM % numLayers))
      nStaLayer = 0
      nRotLayer = 0
      do i = 1, SM % nRingGlobal
         if (SM % statorRing(i) % rank == MPI_Process % rank) then
            nStaLayer(SM % statorRing(i) % Ax_Id) = nStaLayer(SM % statorRing(i) % Ax_Id) + 1
         end if
         if (SM % rotorRing(i) % rank == MPI_Process % rank) then
            nRotLayer(SM % rotorRing(i) % Ax_Id) = nRotLayer(SM % rotorRing(i) % Ax_Id) + 1
         end if
      end do

      nWhole = count(nStaLayer == SM % numElemsPerLayer .and. nRotLayer == SM % numElemsPerLayer)
      nMixed = count(nStaLayer > 0 .and. nRotLayer > 0)
      nHeld  = count(nStaLayer > 0 .or.  nRotLayer > 0)

      if (nMixed > 0 .and. nWhole /= nHeld) then
         write(STD_OUT,*) 'sliding: invalid partition on rank', MPI_Process % rank
         write(STD_OUT,*) 'a rank must hold whole interface layers (strategy 1)'
         write(STD_OUT,*) 'or a single side of each layer it touches (strategy 2)'
         error stop
      end if

      if (MPI_Process % isRoot) then
         write(STD_OUT,'(30X,A,A30,I0," x ",I0)') "->", "Interface ring: ", &
               SM % numElemsPerLayer, SM % numLayers
      end if

   end associate
end subroutine SlidingBuildGlobalRings
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  CollectLocalRingEntries
!
!  One sweep over the local elements, both sides at once, producing one ring
!  entry per interface face. Also marks MortarFaces and counts the interface
!  faces whose partner is local.
!
!  Rotating azimuths are stored in the rotating frame, phi - omega, so that
!  Azim_Id belongs to the element and does not change with time: the two
!  rings then align index by index.
! ---------------------------------------------------------------------------
subroutine CollectLocalRingEntries(mesh, rotorEnt, nRot, statorEnt, nSta)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh),                      intent(inout) :: mesh
   type(SlidingRingEntry), allocatable, intent(out)   :: rotorEnt(:), statorEnt(:)
   integer,                             intent(out)   :: nRot, nSta
   ! =========================
   ! Local variables
   ! =========================
   type(SlidingRingEntry) :: e
   integer                :: eID, j, c, fID, nFacesFound
   real(kind=RP)          :: xc(3)
   integer                :: cc, k
   real(kind=RP)          :: xCorner(3,4)
   intrinsic :: count

   associate(SM => mesh % SlidingMesh)

      allocate(rotorEnt (mesh % no_of_elements))
      allocate(statorEnt(mesh % no_of_elements))
      nRot = 0 
      nSta = 0
      SM % nLocalSplitFaces = 0

      do eID = 1, mesh % no_of_elements
         nFacesFound = 0

         do j = 1, 6
            if (.not. IsInterfaceFace(mesh, eID, j)) cycle
            nFacesFound = nFacesFound+1

            fID = mesh % elements(eID) % faceIDs(j)

            xc = 0.0_RP
            do c = 1, 4
               xc = xc + mesh % nodes(mesh % faces(fID) % nodeIDs(c)) % X
            end do
            xc = 0.25_RP * xc

            e % globID   = mesh % elements(eID) % globID
            e % rank     = MPI_Process % rank
            e % side     = j
            e % rotation = mesh % faces(fID) % rotation
            e % dir      = SlidingFaceDir(SM, mesh % faces(fID) % geom % x)
            e % zeta     = SM % zeta(xc)
            e % jacMin   = minval(mesh % elements(eID) % geom % jacobian)
            e % localEID = eID
            e % Azim_Id  = 0            ! assigned by BuildAzimuthalConnectivity
            e % Ax_Id    = 0            ! assigned by BuildAzimuthalConnectivity
            e % Nel      = mesh % elements(eID) % Nxyz(axisMap(:, j))

            !Orientation tag, see SlidingFaceQuadrants
            do cc = 1, 4
               xCorner(:,cc) = mesh % nodes( mesh % elements(eID) % nodeIDs( localFaceNode(cc, j) ) ) % X
            end do
            e % quad = SlidingFaceQuadrants(SM, xCorner, xc)

            !The four corners must fall in four distinct quadrants, otherwise
            !the tag does not identify them and the pairing breaks silently
            if (any([(count(e % quad == k), k = 0, 3)] /= 1)) then
               write(STD_OUT,*) 'sliding: degenerate corner quadrants on face', fID, e % quad
               error stop
            end if

            mesh % elements(eID) % MortarFaces(j) = 1

            if (mesh % elements(eID) % sliding) then
               ! rotating frame, so that Azim_Id never depends on time
               e % phi = WrapPi(SM % phi(xc) - SM % omega)

               if (mesh % faces(fID) % faceType /= HMESH_MPI) then 
                  SM % nLocalSplitFaces = SM % nLocalSplitFaces+1
               end if

               nRot = nRot+1 
               rotorEnt(nRot) = e
            else
               e % phi = SM % phi(xc)
               nSta = nSta+1 
               statorEnt(nSta) = e
            end if
         end do

         if (nFacesFound > 1) then
            write(STD_OUT,*) 'sliding: element', eID, 'has', nFacesFound, &
                             'interface faces (expected at most 1), closed cylinder not supported'
            error stop
         end if
      end do

      SM % needsLocalSplit = (SM % nLocalSplitFaces > 0)
   end associate

end subroutine CollectLocalRingEntries
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  AllgatherRing
!
!  Replicates the local ring entries of every rank onto every rank. Allocates
!  entAll to the global size.
!
!  Entries are flattened into one integer and one real buffer rather than sent
!  as a derived type, so the layout does not depend on compiler padding.
!  Azim_Id, Ax_Id and localEID are not packed: the first two are produced after
!  the gather, the third is meaningless outside its owner.
!
!  COLLECTIVE. Ranks with zero local entries must still call it.
! ---------------------------------------------------------------------------
subroutine AllgatherRing(entMine, nMine, entAll, nAll)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   type(SlidingRingEntry),              intent(in)  :: entMine(:)
   integer,                             intent(in)  :: nMine
   type(SlidingRingEntry), allocatable, intent(out) :: entAll(:)
   integer,                             intent(out) :: nAll
   ! =========================
   ! Local variables
   ! =========================
   integer                    :: i, p, ierr, offInt, offReal
   integer,       allocatable :: nPerRank(:)
   integer,       allocatable :: intCount(:),  intDispl(:)
   integer,       allocatable :: realCount(:), realDispl(:)
   integer,       allocatable :: intSend(:),   intRecv(:)
   real(kind=RP), allocatable :: realSend(:),  realRecv(:)

   if (.not. MPI_Process % doMPIAction) then
      nAll = nMine
      allocate(entAll(nAll))
      entAll(1:nAll) = entMine(1:nMine)
      return
   end if

#ifdef _HAS_MPI_

   ! How many entries each rank contributes
   allocate(nPerRank(MPI_Process % nProcs))
   call mpi_allgather(nMine, 1, MPI_INT, nPerRank, 1, MPI_INT, MPI_COMM_WORLD, ierr)

   nAll = sum(nPerRank)
   allocate(entAll(max(nAll,1)))

   ! Counts and displacements, one set per buffer
   allocate(intCount (MPI_Process % nProcs), intDispl (MPI_Process % nProcs))
   allocate(realCount(MPI_Process % nProcs), realDispl(MPI_Process % nProcs))

   intCount  = nPerRank * SM_RING_NINT
   realCount = nPerRank * SM_RING_NREAL

   intDispl(1)  = 0
   realDispl(1) = 0
   do p = 2, MPI_Process % nProcs
      intDispl(p) = intDispl(p-1) + intCount(p-1)
      realDispl(p) = realDispl(p-1) + realCount(p-1)
   end do

   allocate(intSend (max(nMine*SM_RING_NINT ,1)), realSend(max(nMine*SM_RING_NREAL,1)))
   allocate(intRecv (max(nAll *SM_RING_NINT ,1)), realRecv(max(nAll *SM_RING_NREAL,1)))

   ! Pack
   do i = 1, nMine
      offInt  = (i-1)*SM_RING_NINT
      offReal = (i-1)*SM_RING_NREAL

      intSend(offInt+1) = entMine(i) % globID
      intSend(offInt+2) = entMine(i) % rank
      intSend(offInt+3) = entMine(i) % side
      intSend(offInt+4) = entMine(i) % rotation
      intSend(offInt+5) = entMine(i) % dir
      intSend(offInt+6) = entMine(i) % Nel(1)
      intSend(offInt+7) = entMine(i) % Nel(2)
      intSend(offInt+8)  = entMine(i) % quad(1)
      intSend(offInt+9)  = entMine(i) % quad(2)
      intSend(offInt+10) = entMine(i) % quad(3)
      intSend(offInt+11) = entMine(i) % quad(4)

      realSend(offReal+1) = entMine(i) % phi
      realSend(offReal+2) = entMine(i) % zeta
      realSend(offReal+3) = entMine(i) % jacMin
   end do

   ! Exchange
   call mpi_allgatherv(intSend, nMine*SM_RING_NINT, MPI_INT, intRecv, intCount, intDispl, MPI_INT,&
                       MPI_COMM_WORLD, ierr)

   call mpi_allgatherv(realSend, nMine*SM_RING_NREAL, MPI_DOUBLE, realRecv, realCount, realDispl, MPI_DOUBLE, &
                       MPI_COMM_WORLD, ierr)
   ! Unpack
   do i = 1, nAll
      offInt  = (i-1)*SM_RING_NINT
      offReal = (i-1)*SM_RING_NREAL

      entAll(i) % globID   = intRecv(offInt+1)
      entAll(i) % rank     = intRecv(offInt+2)
      entAll(i) % side     = intRecv(offInt+3)
      entAll(i) % rotation = intRecv(offInt+4)
      entAll(i) % dir      = intRecv(offInt+5)
      entAll(i) % Nel(1)   = intRecv(offInt+6)
      entAll(i) % Nel(2)   = intRecv(offInt+7)
      entAll(i) % quad(1)  = intRecv(offInt+8)
      entAll(i) % quad(2)  = intRecv(offInt+9)
      entAll(i) % quad(3)  = intRecv(offInt+10)
      entAll(i) % quad(4)  = intRecv(offInt+11)
      entAll(i) % phi      = realRecv(offReal+1)
      entAll(i) % zeta     = realRecv(offReal+2)
      entAll(i) % jacMin   = realRecv(offReal+3)
      
   end do

   deallocate(nPerRank, intCount, intDispl, realCount, realDispl, intSend, realSend, intRecv, realRecv)
#endif
end subroutine AllgatherRing
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  BuildAzimuthalConnectivity
!
!  Assigns the structured interface indices Azim_Id and Ax_Id to a ring, and
!  checks that it is a regular cylindrical ring: equal layer sizes, layers
!  angularly aligned, sum(dphi) = 2*pi, uniform azimuthal spacing.
!
!  Operates on the replicated table: a rank owning only an arc sees neither a
!  closed ring nor complete layers.
!
!  Last place an angle is compared. From here on the interface is addressed by
!  (Azim_Id, Ax_Id, Mortarpos) alone, in integer arithmetic.
! ---------------------------------------------------------------------------
subroutine BuildAzimuthalConnectivity(ent, n, tolz, nPerLayer, nLayers, tag)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   type(SlidingRingEntry), intent(inout) :: ent(:)
   integer,                intent(in)    :: n
   real(kind=RP),          intent(in)    :: tolz
   integer,                intent(out)   :: nPerLayer, nLayers
   character(len=*),       intent(in)    :: tag
   ! =========================
   ! Local variables
   ! =========================
   integer                    :: gStart, gEnd, Ng, NgRef, k, km, kk
   real(kind=RP)              :: sumphi, dphi, dmin, dmax, dal
   real(kind=RP), allocatable :: phiRef(:)

   call SortRingEntries(ent, n, tolz)

   NgRef = -1
   gStart = 1
   nLayers = 0

   do while (gStart <= n)
      gEnd = gStart
      do while (gEnd < n)
         if (ent(gEnd+1) % zeta - ent(gStart) % zeta > tolz) exit
         gEnd = gEnd+1
      end do

      Ng = gEnd - gStart+1
      if (NgRef < 0) NgRef = Ng
      if (Ng /= NgRef) then
         write(STD_OUT,*) 'sliding: ', tag, ' layer sizes differ,', Ng, 'vs', NgRef
         error stop
      end if

      if (.not. allocated(phiRef)) then
         allocate(phiRef(Ng))
         phiRef(1:Ng) = ent(gStart:gEnd) % phi
      else
         do kk = 1, Ng
            dal = modulo(ent(gStart+kk-1) % phi - phiRef(kk), 2.0_RP*PI)
            dal = min(dal, 2.0_RP*PI - dal)
            if (dal > 1.0e-2_RP * (2.0_RP*PI/real(Ng,RP))) then
               write(STD_OUT,*) 'sliding: ', tag, ' axial layers are not angularly aligned,', dal
               error stop
            end if
         end do
      end if

      nLayers = nLayers+1
      ent(gStart:gEnd) % Ax_Id = nLayers
      do kk = 1, Ng
         ent(gStart+kk-1) % Azim_Id = kk 
      end do

      sumphi = 0.0_RP
      dmin = huge(1.0_RP)
      dmax = 0.0_RP

      do k = gStart, gEnd
         km = merge(gEnd, k-1, k == gStart)
         dphi = modulo(ent(k) % phi - ent(km) % phi, 2.0_RP*PI)
         sumphi = sumphi+dphi
         dmin = min(dmin, dphi)
         dmax = max(dmax, dphi)
      end do

      if (abs(sumphi - 2.0_RP*PI) > 1.0e-8_RP) then
         write(STD_OUT,*) 'sliding: ', tag, ' ring not closed, sum(dphi) =', sumphi
         error stop
      end if
      if (dmax/dmin - 1.0_RP > 1.0e-2_RP) then
         write(STD_OUT,*) 'sliding: ', tag, ' non-uniform azimuthal spacing, dphi in [', &
                          dmin, ',', dmax, ']'
         error stop
      end if

      gStart = gEnd+1
   end do

   nPerLayer = NgRef
   safedeallocate(phiRef)

end subroutine BuildAzimuthalConnectivity 
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SortRingEntries
!
!  Orders a ring by axial layer, then azimuth, then global element ID.
!
!  The last key makes the ordering total. Without it two entries with equal
!  azimuth could be ordered differently on two ranks sorting the same
!  replicated table, and the pairing would diverge silently.
! ---------------------------------------------------------------------------
subroutine SortRingEntries(ent, n, tolz)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   type(SlidingRingEntry), intent(inout) :: ent(:)
   integer,                intent(in)    :: n
   real(kind=RP),          intent(in)    :: tolz
   ! =========================
   ! Local variables
   ! =========================
   integer                :: i, j
   type(SlidingRingEntry) :: key
   logical                :: before

   do i = 2, n
      key = ent(i)
      j = i-1
      do while (j >= 1)
         if (key % zeta < ent(j) % zeta-tolz) then
            before = .true.
         else if (key % zeta > ent(j) % zeta+tolz) then
            before = .false.
         else if (key % phi  < ent(j) % phi) then
            before = .true.
         else if (key % phi  > ent(j) % phi) then
            before = .false.
         else
            before = (key % globID < ent(j) % globID)
         end if

         if (.not. before) exit
         ent(j+1) = ent(j)
         j = j-1
      end do
      ent(j+1) = key
   end do

end subroutine SortRingEntries
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SplitInterfaceNodes
!
!  Gives the rotating side its own copy of the interface nodes, node by node:
!  duplicated when a local static element holds it too, relabelled when the
!  partition has already given the static side its own copy elsewhere.
!
!  The new global ID is origGlobID + nodeGlobIDOffset in both cases, so two
!  ranks holding a copy of the same node agree without communicating.
!
!  First step only.
! ---------------------------------------------------------------------------
subroutine SplitInterfaceNodes(mesh)
   use NodeClass
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(inout) :: mesh
   ! =========================
   ! Local variables
   ! =========================
   integer, parameter :: FACE_NODES(4,6) = reshape(&
        [1,2,5,6, 4,3,7,8, 1,2,3,4, 6,2,3,7, 5,6,7,8, 5,1,4,8], [4,6])
   logical, allocatable    :: sharedWithStator(:)
   integer, allocatable    :: remap(:)
   type(Node), allocatable :: tmp(:)
   integer :: nOld, nTot, newNodeCounter, nDup
   integer :: i, j, c, ln, fID, eID, nID, ierr, nGlobNodes

   associate(SM => mesh % SlidingMesh)

      nOld = size(mesh % nodes)
      allocate(sharedWithStator(nOld)) 
      sharedWithStator = .false.
      allocate(remap(nOld))          
      remap = 0

      ! which interface nodes are also held by a LOCAL static element 
      do i = 1, SM % nLocalStatorEnt
         fID = mesh % elements(SM % localStatorEnt(i) % localEID) % faceIDs(SM % localStatorEnt(i) % side)
         do c = 1, 4
            sharedWithStator(mesh % faces(fID) % nodeIDs(c)) = .true.
         end do
      end do

      ! first pass: how many duplicates, in the reference order 
      nDup = 0
      do i = 1, SM % nLocalRotorEnt
         eID = SM % localRotorEnt(i) % localEID
         j = SM % localRotorEnt(i) % side
         do c = 1, 4
            nID = mesh % elements(eID) % nodeIDs(FACE_NODES(c,j))
            if (remap(nID) /= 0) cycle
            if (sharedWithStator(nID)) then
               nDup = nDup+1
               remap(nID) =-1           
            else
               remap(nID) = nID           
            end if
         end do
      end do

      SM % numDuplicatedNodes = nDup      

      ! global node ID offset, identical on every rank 
      nGlobNodes = maxval(mesh % nodes(:) % globID)
#ifdef _HAS_MPI_
      if (MPI_Process % doMPIAction) then
         call mpi_allreduce(MPI_IN_PLACE, nGlobNodes, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, ierr)
      end if
#endif
      SM % nodeGlobIDOffset = nGlobNodes

      ! grow the node array by the local duplicates only 
      if (nDup > 0) then
         nTot = nOld+nDup
         allocate(tmp(nTot))
         do i = 1, nOld
            call ConstructNode(tmp(i), mesh % nodes(i) % x, mesh % nodes(i) % globID)
         end do
         do i = 1, nOld
            call mesh % nodes(i) % destruct
         end do
         safedeallocate(mesh % nodes)
         call move_alloc(tmp, mesh % nodes)
      end if

      ! second pass: create or relabel, in the reference order 
      remap = 0
      newNodeCounter = nOld
      do i = 1, SM % nLocalRotorEnt
         eID = SM % localRotorEnt(i) % localEID
         j = SM % localRotorEnt(i) % side
         do c = 1, 4
            ln= FACE_NODES(c,j)
            nID = mesh % elements(eID) % nodeIDs(ln)
            if (remap(nID) /= 0) cycle

            if (sharedWithStator(nID)) then
               newNodeCounter = newNodeCounter+1
               call ConstructNode(mesh % nodes(newNodeCounter), mesh % nodes(nID) % x, &
                                  mesh % nodes(nID) % globID + SM % nodeGlobIDOffset)
               remap(nID) = newNodeCounter
            else
               mesh % nodes(nID) % globID = mesh % nodes(nID) % globID + SM % nodeGlobIDOffset
               remap(nID) = nID
            end if
         end do
      end do

      ! point every local rotating element at its own copy 
      do eID = 1, mesh % no_of_elements
         if (.not. mesh % elements(eID) % sliding) cycle
         do ln = 1, 8
            nID = mesh % elements(eID) % nodeIDs(ln)
            if (nID <= nOld) then
               if (remap(nID) > 0) mesh % elements(eID) % nodeIDs(ln) = remap(nID)
            end if
         end do
      end do

      deallocate(sharedWithStator, remap)

   end associate
end subroutine SplitInterfaceNodes
!
!////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------------------------------------
!  SlidingBuildGeometryLists
!
!  Lists the rotating elements and the faces they carry, so that the geometry
!  rebuild can be restricted to the region that actually moves.
!
!  Both lists are fixed for the whole run: the topology of the rotating region
!  never changes, only its position.
!  ---------------------------------------------------------------------------
subroutine SlidingBuildGeometryLists(mesh)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(inout) :: mesh
   ! =========================
   ! Local variables
   ! =========================
   logical, allocatable :: faceTaken(:)
   integer :: eID, j, fID, n
   intrinsic :: count

   associate(SM => mesh % SlidingMesh)

      allocate(faceTaken(size(mesh % faces)))
      faceTaken = .false.

      n = 0
      do eID = 1, mesh % no_of_elements
         if (mesh % elements(eID) % sliding) n = n+1
      end do

      safedeallocate(SM % slidingElems)
      allocate(SM % slidingElems(n))

      n = 0
      do eID = 1, mesh % no_of_elements
         if (.not. mesh % elements(eID) % sliding) cycle
         n = n+1
         SM % slidingElems(n) = eID
         do j = 1, 6
            fID = mesh % elements(eID) % faceIDs(j)
            if (fID <= 0) cycle
            faceTaken(fID) = .true.

         end do
      end do

      safedeallocate(SM % slidingFaces)
      allocate(SM % slidingFaces(count(faceTaken)))

      n = 0
      do fID = 1, size(mesh % faces)
         if (.not. faceTaken(fID)) cycle
         n = n+1
         SM % slidingFaces(n) = fID
      end do

      deallocate(faceTaken)

   end associate
   
end subroutine SlidingBuildGeometryLists
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SlidingPruneMPIFaces
!
!  Drops the interface faces from the ordinary MPI face lists and resizes the
!  buffers: their flux crosses through the mortars, and left in place they would
!  be exchanged twice.
!
!  Both ranks apply the same geometric criterion and compaction preserves the
!  order, so the packing order still agrees. No communication.
!
!  Local, no-op sequentially. Call once, after the faces are tagged and before
!  ConstructGeometry.
! ---------------------------------------------------------------------------
subroutine SlidingPruneMPIFaces(mesh)
   use PhysicsStorage
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(inout) :: mesh
   ! =========================
   ! Local variables
   ! =========================
   integer :: k, d, domain, nKeep, fID, mpifID
   integer :: MPI_NDOFS(MPI_Process % nProcs)

   if (.not. MPI_Process % doMPIAction) return

#ifdef _HAS_MPI_
   if (.not. allocated(mesh % MPIfaces % faces)) return

   ! drop the interface faces, preserving order 
   do k = 1, mesh % MPIfaces % nDomainShared
      domain = mesh % MPIfaces % listDomain(k)

      associate ( mf => mesh % MPIfaces % faces(domain) )
         nKeep = 0
         do d = 1, mf % no_of_faces
            if (mesh % faces( mf % faceIDs(d) ) % MortarType == MORTAR_SLIDING) cycle
            nKeep = nKeep+1
            mf % faceIDs(nKeep) = mf % faceIDs(d)
            mf % elementSide(nKeep) = mf % elementSide(d)
         end do
         mf % no_of_faces = nKeep
      end associate
   end do

   ! resize the ordinary buffers on the pruned lists 
   MPI_NDOFS = 0
   do domain = 1, MPI_Process % nProcs
      if (mesh % MPIfaces % faces(domain) % no_of_faces <= 0) cycle
      do mpifID = 1, mesh % MPIfaces % faces(domain) % no_of_faces
         fID = mesh % MPIfaces % faces(domain) % faceIDs(mpifID)
         MPI_NDOFS(domain) = MPI_NDOFS(domain) + product(mesh % faces(fID) % Nf+1)
      end do
   end do

#if defined(NAVIERSTOKES)
   call ConstructMPIFacesStorage(mesh % MPIfaces, NCONS, NGRAD, MPI_NDOFS)
#elif defined(INCNS)
   call ConstructMPIFacesStorage(mesh % MPIfaces, NCONS, NCONS, MPI_NDOFS)
#elif defined(CAHNHILLIARD) && !defined(MULTIPHASE)
   call ConstructMPIFacesStorage(mesh % MPIfaces, NCOMP, NCOMP, MPI_NDOFS)
#elif defined(MULTIPHASE)
   call ConstructMPIFacesStorage(mesh % MPIfaces, NCONS, NCONS, MPI_NDOFS)
#elif defined(ACOUSTIC)
   call ConstructMPIFacesStorage(mesh % MPIfaces, NCONS, NCONS, MPI_NDOFS, NCONSB_in=NCONSB)
#endif
#endif

end subroutine SlidingPruneMPIFaces
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  GetSlidingTopologyState
!
!  Turns the cumulative angle omega into three things: sectorID, the number of
!  whole element widths turned; localAngle, the remainder inside the sector; and
!  isConforming, true when that remainder is zero and the two sides line up.
!  Stops if a step lands within 1e-10 of the next sector boundary, where one of
!  the two mortars would be degenerate.
! ---------------------------------------------------------------------------
subroutine GetSlidingTopologyState(omega, nElements, sectorID, isConforming, localAngle)
   use SMConstants
   use Utilities, only: AlmostEqual
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   real(kind=RP), intent(in)  :: omega
   integer,       intent(in)  :: nElements
   integer,       intent(out) :: sectorID
   logical,       intent(out) :: isConforming
   real(kind=RP), intent(out) :: localAngle
   ! =========================
   ! Local variables
   ! =========================
   real(kind=RP) :: deltaTheta
   real(kind=RP) :: remainder

   ! Angular spacing between two conforming configurations
   deltaTheta = 2.0_RP * PI / real(nElements, RP)

   ! Determine topology sector
   sectorID = floor(omega / deltaTheta)

   ! Local angle inside current sector
   remainder = modulo(omega, deltaTheta)
   localAngle = remainder

   ! Check if configuration is conforming
   isConforming = AlmostEqual(remainder, 0.0_RP)

   if (.not. isConforming .and. deltaTheta - remainder < 1.0e-10_RP * deltaTheta) then
      write(STD_OUT,*) 'slidingmesh: step landed within tolerance of upper sector boundary, remainder =', remainder
      error stop
   end if
end subroutine GetSlidingTopologyState
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  InitializeSlidingConnectivity
!
!  Prepares the stage: rotates the region and computes the mortar offsets and
!  scales (RotateSlidingRegion), moves the nodes on the first call only
!  (SlidingRotateNodes), and clears the mortar entries so nothing from the
!  previous stage survives.
! ---------------------------------------------------------------------------
subroutine InitializeSlidingConnectivity(mesh, dTheta, offsetParams, scaleParams)
   implicit none 
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(inout) :: mesh
   real(kind=RP), intent(in)     :: dTheta                       
   real(kind=RP), intent(inout)  :: offsetParams(4)
   real(kind=RP), intent(inout)  :: scaleParams(4)
   ! =========================
   ! Local variables
   ! =========================
   integer :: l

   offsetParams = 0.0_RP
   scaleParams  = 0.0_RP

   call RotateSlidingRegion(mesh, dTheta, offsetParams, scaleParams)

   ! mesh % nodes only moves on the first call: afterwards the geometry is
   ! rebuilt from SurfInfo, which RotateSlidingRegion rotates at every stage
   if (.not. mesh % SlidingMesh % active) call SlidingRotateNodes(mesh, dTheta)

   ! Reset mortar face data
   if (allocated(mesh%mortar_faces)) then
      do l = 1, SIZE(mesh%mortar_faces)
         mesh % mortar_faces(l) % ID             = -1
         mesh % mortar_faces(l) % FaceType       = HMESH_NONE
         mesh % mortar_faces(l) % rotation       = 0
         mesh % mortar_faces(l) % NelLeft        = -1
         mesh % mortar_faces(l) % NelRight       = -1
         mesh % mortar_faces(l) % NfLeft         = -1
         mesh % mortar_faces(l) % NfRight        = -1
         mesh % mortar_faces(l) % Nf             = -1
         mesh % mortar_faces(l) % nodeIDs        = -1
         mesh % mortar_faces(l) % elementIDs     = -1
         mesh % mortar_faces(l) % elementSide    = -1
         mesh % mortar_faces(l) % projectionType = -1
         mesh % mortar_faces(l) % boundaryName   = ""
      end do
   end if

end subroutine InitializeSlidingConnectivity
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  RotateSlidingRegion
!
!  Rotates the geometry of the rotating elements by dTheta: the corners or the
!  face patches stored in SurfInfo, from which ConstructGeometry rebuilds the
!  mappings. Also computes the mortar offsets and scales from the position
!  inside the current sector. Called at every Runge-Kutta stage.
! ---------------------------------------------------------------------------
subroutine RotateSlidingRegion(mesh, dTheta, offsetParams, scaleParams)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(inout) :: mesh
   real(kind=RP), intent(in)     :: dTheta                       ! rotation increment applied by this call:
                                                                 ! the RK stage fraction of theta, not theta itself
   real(kind=RP), intent(inout)  :: offsetParams(4)              ! start of each mortar piece
   real(kind=RP), intent(inout)  :: scaleParams(4)               ! width of each mortar piece
   ! =========================
   ! Local variables
   ! =========================
   integer :: i, j, uIdx, vIdx, cornerIdx
   integer :: eID
   real(kind=RP) :: mortarShiftParam                ! position inside the sector, in [-1,1]
   real(kind=RP) :: deltaTheta
   real(KIND=RP), allocatable :: facePatchPoints(:,:,:)
   real(KIND=RP), allocatable :: uNodes(:)
   real(KIND=RP), allocatable :: vNodes(:)

   offsetParams = 0.0_RP
   scaleParams  = 0.0_RP

   ! Position inside the current sector
   deltaTheta       = 2.0_RP * PI / real(mesh % SlidingMesh % numElemsPerLayer, RP)
   mortarShiftParam = 1.0_RP - 2.0_RP * mesh % SlidingMesh % localAngle / deltaTheta

   ! Mortar geometric parameters
   offsetParams(1) = (mortarShiftParam - 1.0_RP) / 2.0_RP
   offsetParams(2) = (1.0_RP - mortarShiftParam) / 2.0_RP
   offsetParams(3) = (1.0_RP + mortarShiftParam) / 2.0_RP
   offsetParams(4) = (-mortarShiftParam - 1.0_RP) / 2.0_RP

   scaleParams(1) = offsetParams(1) + 1.0_RP
   scaleParams(2) = 1.0_RP - offsetParams(2)
   scaleParams(3) = 1.0_RP - offsetParams(3)
   scaleParams(4) = offsetParams(4) + 1.0_RP

   ! Rotate sliding element geometry
!$omp parallel do schedule(runtime) default(shared) &
!$omp private(i, eID, j, cornerIdx, uIdx, vIdx, uNodes, vNodes, facePatchPoints)
   do i = 1, mesh%no_of_elements
      if (.not. mesh % elements(i)%sliding) cycle
      eID = i
   
      ! Case 1: Hex8 element (corner-based geometry)
      if (mesh % elements(eID)%SurfInfo%IsHex8) then
         do cornerIdx = 1, 8
            call mesh % SlidingMesh % RotatePoint(dTheta, mesh % elements(eID) % SurfInfo % corners(:,cornerIdx))
         end do
      else

         ! Case 2: High-order faces
         do j = 1, 6
   
            if (.not. allocated(mesh % elements(eID)%SurfInfo%facePatches(j)%uKnots)) cycle

            ! Flat face (2x2)
            if (mesh % elements(eID)%SurfInfo%facePatches(j)%noOfKnots(1) == 2) then
   
               call mesh % SlidingMesh % RotatePoint(dTheta, mesh % elements(eID)%SurfInfo%facePatches(j)%points(:,1,1))
               call mesh % SlidingMesh % RotatePoint(dTheta, mesh % elements(eID)%SurfInfo%facePatches(j)%points(:,2,1))
               call mesh % SlidingMesh % RotatePoint(dTheta, mesh % elements(eID)%SurfInfo%facePatches(j)%points(:,2,2))
               call mesh % SlidingMesh % RotatePoint(dTheta, mesh % elements(eID)%SurfInfo%facePatches(j)%points(:,1,2))
   
            else
               
               ! Curved face (high-order patch)
               uNodes = mesh % elements(eID)%SurfInfo%facePatches(j)%uKnots
               vNodes = mesh % elements(eID)%SurfInfo%facePatches(j)%vKnots
             
               facePatchPoints = mesh % elements(eID)%SurfInfo%facePatches(j)%points
               
               do uIdx = 1, size(facePatchPoints, 3)
                  do vIdx = 1, size(facePatchPoints, 2)
                     call mesh % SlidingMesh % RotatePoint(dTheta, facePatchPoints(:,vIdx,uIdx))
                  end do
               end do
               
               if (.not. mesh % SlidingMesh % active) then
                  call mesh % elements(eID)%SurfInfo%facePatches(j)%destruct
               end if
               
               call mesh % elements(eID)%SurfInfo%facePatches(j)%construct(uNodes, vNodes, facePatchPoints)
   
            end if
         end do ! j
      end if
   end do ! i
!$omp end parallel do

end subroutine RotateSlidingRegion
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SlidingRotateNodes
!
!  Rotates the nodes of the rotating elements by dTheta, each node once.
!  First call only: afterwards nothing reads the node positions, the geometry
!  is rebuilt from SurfInfo.
! ---------------------------------------------------------------------------
subroutine SlidingRotateNodes(mesh, dTheta)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(inout) :: mesh
   real(kind=RP),  intent(in)    :: dTheta
   ! =========================
   ! Local variables
   ! =========================
   logical, allocatable :: done(:)
   integer :: eID, j, nID

   allocate(done(size(mesh % nodes)))
   done = .false.

   do eID = 1, mesh % no_of_elements
      if (.not. mesh % elements(eID) % sliding) cycle
      do j = 1, 8
         nID = mesh % elements(eID) % nodeIDs(j)
         if (done(nID)) cycle
         call mesh % SlidingMesh % RotatePoint(dTheta, mesh % nodes(nID) % X)
         done(nID) = .true.
      end do
   end do

   deallocate(done)

end subroutine SlidingRotateNodes
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  UpdateSlidingMortarsConnectivity
!
!  Refills slidingMortarConnectivity from the ring directory, for the given
!  sector, and lists the static rows this rank is the slave of.
!
!  Driven by the static interface elements this rank owns, not by its rotating 
!  elements: the static side is the master, and under a partition that splits the 
!  interface a rank may own one side and not the other. A column left at zero 
!  means the element is on another rank.
!
!  For static row j, SlidingSlaveRows gives the two rotating rows straddling it:
!  rows(0) for Mortarpos 0, rows(1) for Mortarpos 1.
!
!  Called once at setup and at every sector change.
! ---------------------------------------------------------------------------
subroutine UpdateSlidingMortarsConnectivity(mesh, sectorID)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(inout) :: mesh
   integer,        intent(in)    :: sectorID
   ! =========================
   ! Local variables
   ! =========================
   integer :: j, m, rows(0:1)

   associate(SM => mesh % SlidingMesh)

      ! static rows this rank is the slave of, for this sector 
      SM % nSlaveRows = 0
      do j = 1, SM % nRingGlobal
         call SlidingSlaveRows(SM, j, sectorID, rows)
         if (SM % rotorLocalElem(rows(0)) > 0 .or. &
            SM % rotorLocalElem(rows(1)) > 0) then
            SM % nSlaveRows = SM % nSlaveRows+1
            SM % mySlaveRows(SM % nSlaveRows) = j
         end if
      end do

      ! one entry per static row I own, both mortar families in it 
      SM % slidingMortarConnectivity = 0

      do m = 1, SM % nMasterRows
         j = SM % myMasterRows(m)
         call SlidingSlaveRows(SM, j, sectorID, rows)

         ! slave of Mortarpos 1 
         SM % slidingMortarConnectivity(m,1) = SM % rotorLocalElem(rows(1))
         SM % slidingMortarConnectivity(m,4) = SM % rotorRing(rows(1)) % side
         SM % slidingMortarConnectivity(m,7) = SM % rotorRing(rows(1)) % rotation

         ! slave of Mortarpos 0 
         SM % slidingMortarConnectivity(m,2) = SM % rotorLocalElem(rows(0))
         SM % slidingMortarConnectivity(m,5) = SM % rotorRing(rows(0)) % side
         SM % slidingMortarConnectivity(m,8) = SM % rotorRing(rows(0)) % rotation

         ! master, the static row itself 
         SM % slidingMortarConnectivity(m,3) = SM % statorLocalElem(j)
         SM % slidingMortarConnectivity(m,6) = SM % statorRing(j) % side
         SM % slidingMortarConnectivity(m,9) = SM % statorRing(j) % rotation

         !orientation tags of the three faces, for the current sector
         SM % slidingMortarConnectivity(m,10:13) = SM % rotorRing(rows(1)) % quad   ! Mortarpos 1 slave
         SM % slidingMortarConnectivity(m,14:17) = SM % rotorRing(rows(0)) % quad   ! Mortarpos 0 slave
         SM % slidingMortarConnectivity(m,18:21) = SM % statorRing(j)      % quad   ! master
      end do

   end associate
end subroutine UpdateSlidingMortarsConnectivity
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SlidingSplitInterfaceFaces
!
!  Gives every locally shared interface face a second copy, so that the static
!  and the rotating element each own one. The existing face keeps its number
!  and becomes the static element's; the rotating element gets a new face,
!  appended after the last one, built on the node copies made by
!  SplitInterfaceNodes.
!
!  Nothing else is rebuilt: every other face keeps the number, type and links
!  the mesh reader gave it, periodic pairs, MPI faces and 4:1 mortars included.
!  Interface faces shared with another rank are already one-sided and are left
!  alone.
!
!  Local. Only reached when needsLocalSplit is true, on the first call, after
!  SplitInterfaceNodes.
! ---------------------------------------------------------------------------
subroutine SlidingSplitInterfaceFaces(mesh)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(inout) :: mesh
   ! =========================
   ! Local variables
   ! =========================
   type(Face), allocatable :: oldFaces(:)
   integer :: n, k, r, c, fID
   integer :: rotorEID, rotorSide, statorEID, statorSide
   integer :: rotorNodes(4)

   associate(SM => mesh % SlidingMesh)

      !Grow the face array, keeping every existing face and its number
      n = mesh % numberOfFaces
      call move_alloc(mesh % faces, oldFaces)
      allocate(mesh % faces(n + SM % nLocalSplitFaces))
      mesh % faces(1:n) = oldFaces(1:n)
      deallocate(oldFaces)

      k = 0
      do r = 1, SM % nRingGlobal
         rotorEID = SM % rotorLocalElem(r)
         if (rotorEID <= 0) cycle
         rotorSide = SM % rotorRing(r) % side
         fID = mesh % elements(rotorEID) % faceIDs(rotorSide)

         !Already one-sided: shared with another rank, nothing to split here
         if (mesh % faces(fID) % faceType == HMESH_MPI) cycle
         if (mesh % faces(fID) % elementIDs(2) <= 0)    cycle

         !The static element on the other side of this face
         if (mesh % faces(fID) % elementIDs(1) == rotorEID) then
            statorEID  = mesh % faces(fID) % elementIDs(2)
            statorSide = mesh % faces(fID) % elementSide(2)
         else
            statorEID  = mesh % faces(fID) % elementIDs(1)
            statorSide = mesh % faces(fID) % elementSide(1)
         end if

         !The existing face keeps its number and becomes the static element's
         !alone, written in that element's own orientation
         do c = 1, 4
            mesh % faces(fID) % nodeIDs(c) = &
               mesh % elements(statorEID) % nodeIDs(localFaceNode(c, statorSide))
         end do
         mesh % faces(fID) % elementIDs  = [statorEID,  HMESH_NONE]
         mesh % faces(fID) % elementSide = [statorSide, HMESH_NONE]
         mesh % faces(fID) % rotation    = 0
         mesh % elements(statorEID) % faceSide(statorSide) = 1

         !A new face for the rotating element, on its own node copies
         k = k + 1
         do c = 1, 4
            rotorNodes(c) = mesh % elements(rotorEID) % nodeIDs(localFaceNode(c, rotorSide))
         end do
         call mesh % faces(n+k) % Construct(ID = n+k, nodeIDs = rotorNodes, &
                                            elementID = rotorEID, side = rotorSide)
         mesh % elements(rotorEID) % faceIDs(rotorSide)  = n+k
         mesh % elements(rotorEID) % faceSide(rotorSide) = 1
      end do

      !The two counts come from two different sweeps over the same faces
      if (k /= SM % nLocalSplitFaces) then
         write(STD_OUT,*) 'sliding: split', k, 'interface faces, expected', SM % nLocalSplitFaces
         error stop
      end if

      mesh % numberOfFaces = n + k

   end associate

end subroutine SlidingSplitInterfaceFaces
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  ConstructSlidingMortars
!
!  Builds the mortars when both sides are on this rank (sequential and
!  strategy 1). Two mortars per static interface face: Mortarpos 0 and 1
!  cover the two parts of the overlap. Each gets its two elements and faces,
!  its rotation, offset and scale, its projection matrices (LinkWithElements)
!  and its geometry.
!
!  At localAngle = 0, Mortarpos 0 has scale 1 and covers the whole face;
!  Mortarpos 1 has scale 0 and contributes nothing.
!
!  Called at every Runge-Kutta stage.
! ---------------------------------------------------------------------------
subroutine ConstructSlidingMortars(mesh, nodes, nelm, slidingMortarConnectivity, offsetParams, scaleParams)
   implicit none
   ! ================
   ! Arguments
   ! ================
   class(HexMesh), intent(inout) :: mesh
   integer, intent(in)           :: nodes
   integer, intent(in)           :: nelm
   integer, intent(in)           :: slidingMortarConnectivity(nelm,21)
   real(kind=RP), intent(in)     :: offsetParams(4)
   real(kind=RP), intent(in)     :: scaleParams(4)
   ! ================
   ! Local variables
   ! ================
   integer :: i, j, Mortarpos, mortarIndex
   integer :: masterFaceID, slaveFaceID
   integer :: masterElementID, slaveElementID
   integer :: masterFaceNumber, slaveFaceNumber
   integer :: mortarFaceNodeIDs(4)
   integer :: elementNodeIDs(8)
   integer :: masterQuad(4), slaveQuad(4)
   integer :: NelL(2), NelR(2)
   real(kind=RP) :: jmax

   do i = 1, nelm

      ! The master is the same static face for both families
      masterElementID = slidingMortarConnectivity(i,3)
      if (masterElementID <= 0) cycle ! unused row
      masterFaceNumber = slidingMortarConnectivity(i,6)
      masterFaceID = mesh % elements(masterElementID) % faceIDs(masterFaceNumber)

      elementNodeIDs = mesh % elements(masterElementID) % nodeIDs
      do j = 1, 4
         mortarFaceNodeIDs(j) = elementNodeIDs(localFaceNode(j, masterFaceNumber))
      end do

      if (.not. allocated(mesh % faces(masterFaceID) % Mortar)) then
         allocate(mesh % faces(masterFaceID) % Mortar(2))
         mesh % faces(masterFaceID) % Mortar = 0
      end if

      do Mortarpos = 0, 1

         mortarIndex = MortarSlot(mesh % SlidingMesh, Mortarpos, i)

         if (Mortarpos == 0) then
            slaveElementID  = slidingMortarConnectivity(i,2)
            slaveFaceNumber = slidingMortarConnectivity(i,5)
         else
            slaveElementID  = slidingMortarConnectivity(i,1)
            slaveFaceNumber = slidingMortarConnectivity(i,4)
         end if
         if (slaveElementID <= 0) cycle ! unused row

         !Orientation tags of this pair, for the current sector
         masterQuad = slidingMortarConnectivity(i,18:21)
         if (Mortarpos == 0) then
            slaveQuad = slidingMortarConnectivity(i,14:17)
         else
            slaveQuad = slidingMortarConnectivity(i,10:13)
         end if

         call mesh % mortar_faces(mortarIndex) % Construct(ID = mortarIndex, nodeIDs   = mortarFaceNodeIDs, &
                                                            elementID = masterElementID, side = masterFaceNumber)

         mesh % mortar_faces(mortarIndex) % Mortarpos = Mortarpos
         mesh % mortar_faces(mortarIndex) % FaceType  = HMESH_INTERIOR

         if (.not. allocated(mesh % mortar_faces(mortarIndex) % Mortar)) then
            allocate(mesh % mortar_faces(mortarIndex) % Mortar(2))
            mesh % mortar_faces(mortarIndex) % Mortar = 0
         end if

         mesh % mortar_faces(mortarIndex) % Mortar(1) = masterFaceID
         mesh % faces(masterFaceID) % Mortar(1+Mortarpos) = mortarIndex

         mesh % mortar_faces(mortarIndex) % elementIDs(2)  = slaveElementID
         mesh % mortar_faces(mortarIndex) % elementSide(2) = slaveFaceNumber

         mesh % elements(masterElementID) % faceSide(masterFaceNumber) = 1
         mesh % elements(slaveElementID)  % faceSide(slaveFaceNumber)  = 2

         slaveFaceID = mesh % elements(slaveElementID) % faceIDs(slaveFaceNumber)

         mesh % mortar_faces(mortarIndex) % Mortar(2) = slaveFaceID

         if (.not. allocated(mesh % faces(slaveFaceID) % Mortar)) then
            allocate(mesh % faces(slaveFaceID) % Mortar(2))
            mesh % faces(slaveFaceID) % Mortar = 0
         end if

         mesh % faces(slaveFaceID) % Mortar(2-Mortarpos) = mortarIndex

         call SetSlidingMortarPlacement(mesh % mortar_faces(mortarIndex), masterQuad, slaveQuad, &
                    mesh % SlidingMesh % statorRing(mesh % SlidingMesh % myMasterRows(i)) % dir, &
                    Mortarpos, offsetParams, scaleParams)

         mesh % elements(masterElementID) % MortarFaces(masterFaceNumber) = 3
         mesh % elements(slaveElementID)  % MortarFaces(slaveFaceNumber)  = 4

      end do
   end do

   ! Link mortar faces with neighbouring elements
!$omp parallel do schedule(runtime) if (mesh % SlidingMesh % active) private(i, NelL, NelR)
   do i = 1, size(mesh % mortar_faces)
      associate(f => mesh % mortar_faces(i))
         if (f % elementIDs(1) <= 0 .or. f % elementIDs(2) <= 0) cycle
         associate(eL => mesh % elements(f % elementIDs(1)), &
                   eR => mesh % elements(f % elementIDs(2)))
            NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
            NelR = eR % Nxyz(axisMap(:, f % elementSide(2)))
            call f % LinkWithElements(NelL, NelR, nodes, f % offset, f % s)
         end associate
      end associate
   end do
!$omp end parallel do

   ! Mortar geometrical mappings
!$omp parallel do schedule(runtime) if (mesh % SlidingMesh % active) private(i, NelL, NelR, jmax)
   do i = 1, size(mesh % mortar_faces)
      associate(f => mesh % mortar_faces(i))
         if (f % elementIDs(1) <= 0 .or. f % elementIDs(2) <= 0) cycle

         associate(eL => mesh % elements(f % elementIDs(1)), &
                   eR => mesh % elements(f % elementIDs(2)))
            NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
            NelR = eR % Nxyz(axisMap(:, f % elementSide(2)))

            if (f % s(1) > 0.0_RP) then
               call f % geom % construct(f % Nf, f % NelLeft, f % NfLeft, eL % Nxyz, &
                                         eL % geom, eL % hexMap, f % elementSide(1), &
                                         f % projectionType(1), 1, 0, .true.,        &
                                         f % MIntSliding(:,:,1), f % s(1), f % slidingDir)
            else
               call f % geom % construct(f % Nf, f % NelLeft, f % NfLeft, eL % Nxyz, &
                                         eL % geom, eL % hexMap, f % elementSide(1), &
                                         f % projectionType(1), 1, 0)
            end if
         end associate

         jmax = maxval(f % geom % jacobian)
         if (jmax <= 0.0_RP) then
            write(STD_OUT,*) 'FATAL: degenerate mortar face, fID=', f % ID
            error stop
         end if
         f % geom % h = minval(mesh % elements(f % elementIDs(2)) % geom % jacobian) &
                      / maxval(f % geom % jacobian)

      end associate
   end do
!$omp end parallel do

end subroutine ConstructSlidingMortars
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SlidingConstructSplitMortars
!
!  Strategy 2: builds the mortars whose two sides are on different ranks. Each
!  rank builds its own copy. The static rank gets the full mortar with its
!  geometry, the rotating rank a light copy without geometry
!  (elementIDs(1) = 0). What a rank does not know about the remote element
!  comes from the directory.
!
!  Called at every Runge-Kutta stage.
! ---------------------------------------------------------------------------
subroutine SlidingConstructSplitMortars(mesh, nodes, offsetParams, scaleParams)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   class(HexMesh), intent(inout) :: mesh
   integer,        intent(in)    :: nodes
   real(kind=RP),  intent(in)    :: offsetParams(4)
   real(kind=RP),  intent(in)    :: scaleParams(4)
   ! =========================
   ! Local variables
   ! =========================
   integer :: m, j, r, k, Mortarpos, mortarIndex, rows(0:1)
   integer :: eID, side, fID
   integer :: NelL(2), NelR(2)
   integer :: mortarFaceNodeIDs(4), masterQuad(4), slaveQuad(4)
   real(kind=RP) :: jmax

   if (.not. MPI_Process % doMPIAction) return

   associate(SM => mesh % SlidingMesh)

!     Static side: the full mortar
!$omp parallel do schedule(runtime) if (SM % active) &
!$omp private(m, j, rows, eID, side, fID, k, mortarFaceNodeIDs, NelL, NelR, &
!$omp masterQuad, slaveQuad, Mortarpos, r, mortarIndex, jmax)
      do m = 1, SM % nMasterRows
         j = SM % myMasterRows(m)
         call SlidingSlaveRows(SM, j, SM % currentSectorID, rows)

         eID  = SM % statorLocalElem(j)
         side = SM % statorRing(j) % side
         fID  = mesh % elements(eID) % faceIDs(side)
         do k = 1, 4
            mortarFaceNodeIDs(k) = mesh % elements(eID) % nodeIDs(localFaceNode(k, side))
         end do
         NelL       = mesh % elements(eID) % Nxyz(axisMap(:, side))
         masterQuad = SM % statorRing(j) % quad

         if (.not. allocated(mesh % faces(fID) % Mortar)) then
            allocate(mesh % faces(fID) % Mortar(2))
            mesh % faces(fID) % Mortar = 0
         end if

         do Mortarpos = 0, 1
            r = rows(Mortarpos)
            if (SM % rotorLocalElem(r) > 0) cycle      ! both sides local: ConstructSlidingMortars

            mortarIndex = MortarSlot(SM, Mortarpos, m)
            NelR        = SM % rotorRing(r) % Nel
            slaveQuad   = SM % rotorRing(r) % quad

            associate(f => mesh % mortar_faces(mortarIndex))

            call f % Construct(ID = mortarIndex, nodeIDs = mortarFaceNodeIDs, &
                               elementID = eID, side = side)
            f % Mortarpos = Mortarpos
            f % FaceType  = HMESH_INTERIOR

            if (.not. allocated(f % Mortar)) allocate(f % Mortar(2))
            f % Mortar = [fID, 0]                      ! slave face on another rank
            mesh % faces(fID) % Mortar(1 + Mortarpos) = mortarIndex

            f % elementIDs(2)  = 0                     ! slave element on another rank
            f % elementSide(2) = SM % rotorRing(r) % side

            mesh % elements(eID) % faceSide(side)    = 1
            mesh % elements(eID) % MortarFaces(side) = 3

            call SetSlidingMortarPlacement(f, masterQuad, slaveQuad, SM % statorRing(j) % dir, &
                                           Mortarpos, offsetParams, scaleParams)
            call f % LinkWithElements(NelL, NelR, nodes, f % offset, f % s)

            associate(eL => mesh % elements(eID))
            if (f % s(1) > 0.0_RP) then
               call f % geom % construct(f % Nf, f % NelLeft, f % NfLeft, eL % Nxyz, &
                                         eL % geom, eL % hexMap, f % elementSide(1), &
                                         f % projectionType(1), 1, 0, .true.,        &
                                         f % MIntSliding(:,:,1), f % s(1), f % slidingDir)
            else
               call f % geom % construct(f % Nf, f % NelLeft, f % NfLeft, eL % Nxyz, &
                                         eL % geom, eL % hexMap, f % elementSide(1), &
                                         f % projectionType(1), 1, 0)
            end if
            end associate

            jmax = maxval(f % geom % jacobian)
            if (jmax <= 0.0_RP) then
               write(STD_OUT,*) 'FATAL: degenerate mortar face, fID=', f % ID
               error stop
            end if
            ! The slave element is not local: its minimum jacobian is in the directory
            f % geom % h = SM % rotorRing(r) % jacMin / jmax

            end associate
         end do
      end do
!$omp end parallel do

      ! Rotating side: the reduced record, no geometry
      do m = 1, SM % nSlaveRows
         j = SM % mySlaveRows(m)
         if (SM % statorLocalElem(j) > 0) cycle        ! both sides local: ConstructSlidingMortars
         call SlidingSlaveRows(SM, j, SM % currentSectorID, rows)

         NelL       = SM % statorRing(j) % Nel
         masterQuad = SM % statorRing(j) % quad

         do Mortarpos = 0, 1
            r = rows(Mortarpos)
            if (SM % rotorLocalElem(r) <= 0) cycle     ! this half is not mine

            mortarIndex = MortarSlot(SM, Mortarpos, m)

            eID  = SM % rotorLocalElem(r)
            side = SM % rotorRing(r) % side
            fID  = mesh % elements(eID) % faceIDs(side)
            do k = 1, 4
               mortarFaceNodeIDs(k) = mesh % elements(eID) % nodeIDs(localFaceNode(k, side))
            end do
            NelR      = mesh % elements(eID) % Nxyz(axisMap(:, side))
            slaveQuad = SM % rotorRing(r) % quad

            if (.not. allocated(mesh % faces(fID) % Mortar)) then
               allocate(mesh % faces(fID) % Mortar(2))
               mesh % faces(fID) % Mortar = 0
            end if

            associate(f => mesh % mortar_faces(mortarIndex))

            call f % Construct(ID = mortarIndex, nodeIDs = mortarFaceNodeIDs, &
                               elementID = eID, side = side)
            f % Mortarpos = Mortarpos
            f % FaceType  = HMESH_INTERIOR

!           Mirror of the static copy: the local element sits on side 2,
!           elementIDs(1) = 0 marks the record as reduced
            f % elementIDs  = [0, eID]
            f % elementSide = [SM % statorRing(j) % side, side]

            if (.not. allocated(f % Mortar)) allocate(f % Mortar(2))
            f % Mortar = [0, fID]                      ! master face on another rank
            mesh % faces(fID) % Mortar(2 - Mortarpos) = mortarIndex

            mesh % elements(eID) % faceSide(side)    = 2
            mesh % elements(eID) % MortarFaces(side) = 4

            call SetSlidingMortarPlacement(f, masterQuad, slaveQuad, SM % statorRing(j) % dir, &
                                           Mortarpos, offsetParams, scaleParams)
            call f % LinkWithElements(NelL, NelR, nodes, f % offset, f % s)
!           No geometry: this side never solves a Riemann problem

            end associate
         end do
      end do

   end associate

end subroutine SlidingConstructSplitMortars
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SetSlidingMortarPlacement
!
!  Rotation, azimuthal direction, offsets and scales of a sliding mortar, from
!  the directory alone. Shared by the local and split constructions, so that
!  both copies of a split mortar get the same operators.
! ---------------------------------------------------------------------------
subroutine SetSlidingMortarPlacement(f, masterQuad, slaveQuad, masterDir, Mortarpos, &
   offsetParams, scaleParams)
   implicit none
   type(Face),    intent(inout) :: f
   integer,       intent(in)    :: masterQuad(4), slaveQuad(4), masterDir, Mortarpos
   real(kind=RP), intent(in)    :: offsetParams(4), scaleParams(4)

   integer :: sM

   f % rotation = SlidingPairRotation(masterQuad, slaveQuad)

   !  Crossed azimuthal directions require an axis-swapping rotation
   if (QuadAzimDir(masterQuad) /= QuadAzimDir(slaveQuad)) then
      if (all(f % rotation /= [1, 3, 4, 6])) then
         write(STD_OUT,*) 'sliding: mortar', f % ID, ' has crossed azimuthal', &
         ' directions but a non axis-swapping rotation.', masterQuad, slaveQuad
         error stop
      end if
   end if

   f % slidingDir = masterDir

   !  Master's sign on both sides: the slave is already in the master's frame (leftIndexes2Right).
   sM = QuadAzimSign(masterQuad)

   f % offset(1) = sM * offsetParams(2*Mortarpos+2)
   f % offset(2) = sM * offsetParams(2*Mortarpos+1)
   f % s(1)      = scaleParams(2*Mortarpos+1)
   f % s(2)      = scaleParams(2*Mortarpos+2)

end subroutine SetSlidingMortarPlacement
!
!////////////////////////////////////////////////////////////////////////////////
!

! ---------------------------------------------------------------------------
!  SlidingFaceDir
!
!  The contracted direction must be the azimuthal one, because that is where the
!  mortar overlap lives. Which local direction that is depends on the mesh
!  generator, not on the geometry. 
!
!                   SpecMesh                           HOPR
!             +-----------------+               +-----------------+
!             |  .   .   .   .  |               |  .   .   o   .  |
!             |                 |               |          |      |
!             |  o---o---o---o  |               |  .   .   o   .  |
!             |                 |               |          |      |
!             |  .   .   .   .  |               |  .   .   o   .  |
!      axis   |                 |        axis   |          |      |
!       ^     |  .   .   .   .  |         ^     |  .   .   o   .  |
!     i |     +-----------------+       j |     +-----------------+
!       |                                 |
!       +----> azimuth                    +----> azimuth
!          j                                 i  
!              dir 1 = axis                      dir 1 = azimuth
!              dir 2 = azimuth                   dir 2 = axis
!
!  Left, contracting dir 2 mixes along the azimuth: correct.
!  Right, it would mix along the axis and ignore the slip entirely
! ---------------------------------------------------------------------------
integer function SlidingFaceDir(SM, x) result(dir)
   use SlidingMeshClass
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   type(SlidingMesh), intent(in) :: SM
   real(kind=RP), intent(in) :: x(:,0:,0:)
   ! =========================
   ! Local variables
   ! =========================
   integer :: n1, n2
   real(kind=RP) :: phi_origin, d1, d2

   n1 = ubound(x,2)
   n2 = ubound(x,3)
   phi_origin = SM % phi(x(:,0,0))
   d1 = abs(WrapPi(SM % phi(x(:,n1,0))-phi_origin))
   d2 = abs(WrapPi(SM % phi(x(:,0,n2))-phi_origin))
   dir = merge(1, 2, d1 > d2)

end function SlidingFaceDir
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SlidingFaceQuadrants
!
!  Orientation tag of an interface face: quadrant of each corner, in
!  localFaceNode order, relative to the face centre,
!
!     quad = 2*(dphi > 0) + (dzeta > 0)          in {0,1,2,3}
!
!  Geometric and invariant under rigid rotation: computed once at setup.
! ---------------------------------------------------------------------------
function SlidingFaceQuadrants(SM, xCorner, xc) result(quad)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   type(SlidingMesh), intent(in) :: SM
   real(kind=RP),     intent(in) :: xCorner(3,4)    ! corners, localFaceNode order
   real(kind=RP),     intent(in) :: xc(3)           ! face centre
   integer                       :: quad(4)
   ! =========================
   ! Local variables
   ! =========================
   integer       :: cc
   real(kind=RP) :: pc, zc, dp, dz

   pc = SM % phi(xc)
   zc = SM % zeta(xc)

   do cc = 1, 4
      dp = WrapPi(SM % phi(xCorner(:,cc)) - pc)
      dz = SM % zeta(xCorner(:,cc)) - zc
      quad(cc) = 2*merge(1, 0, dp > 0.0_RP) + merge(1, 0, dz > 0.0_RP)
   end do

end function SlidingFaceQuadrants
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  SlidingPairRotation
!
!  Relative rotation (0..7) of two faces, computed as faceRotation (MeshTypes)
!  with the corner quadrants as corner identifiers.
! ---------------------------------------------------------------------------
integer pure function SlidingPairRotation(masterQuad, slaveQuad) result(rot)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   integer, intent(in) :: masterQuad(4), slaveQuad(4)
   ! =========================
   ! Local variables
   ! =========================
   integer, parameter  :: NEXTNODE(4) = [2, 3, 4, 1]
   integer             :: j

!  Rotate until both first corners match
   do j = 1, 4
      if (masterQuad(1) == slaveQuad(j)) exit
   end do

!  Same or opposite orientation
   if (masterQuad(2) == slaveQuad(NEXTNODE(j))) then
      rot = j - 1
   else
      rot = j + 3
   end if

end function SlidingPairRotation
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  QuadAzimSign
!
!  +1 if the azimuth increases along the face's azimuthal reference
!  coordinate, -1 otherwise.
! ---------------------------------------------------------------------------
integer pure function QuadAzimSign(quad) result(s)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   integer, intent(in) :: quad(4)

   if (quad(1)/2 /= quad(2)/2) then
      s = merge(1, -1, quad(2)/2 > quad(1)/2)     ! azimuth along the first coordinate
   else
      s = merge(1, -1, quad(4)/2 > quad(1)/2)     ! azimuth along the second
   end if

end function QuadAzimSign
!
!////////////////////////////////////////////////////////////////////////////////
!
! ---------------------------------------------------------------------------
!  QuadAzimDir
!
!  Which of the face's two reference coordinates carries the azimuth. Same
!  answer as SlidingFaceDir, from the tag rather than from the geometry, and
!  available for a face this rank does not own.
! ---------------------------------------------------------------------------
integer pure function QuadAzimDir(quad) result(d)
   implicit none
   ! =========================
   ! Arguments
   ! =========================
   integer, intent(in) :: quad(4)

   d = merge(1, 2, quad(1)/2 /= quad(2)/2)

end function QuadAzimDir
!
!////////////////////////////////////////////////////////////////////////////////
!
real(kind=RP) function WrapPi(a) result(b)
   implicit none
   real(kind=RP), intent(in) :: a

   b = a
   do while (b > PI) 
      b = b - 2.0_RP*PI 
    end do

   do while (b < -PI) 
      b = b + 2.0_RP*PI 
    end do

end function WrapPi
!
!////////////////////////////////////////////////////////////////////////////////
!  
subroutine PrintMortarConnectivity(mesh)
   implicit none
   class(HexMesh), intent(inout) :: mesh

   integer :: i

   do i=1 , size(mesh % faces)
     
      if (mesh % faces(i) % MortarType == MORTAR_SLIDING) then 
         write(*,*) 'face',i,'is connected to mortarface; connection:',mesh % faces(i) %  mortar, 'elements:', mesh % faces(i) % elementIDs 
      end if 
   end do 

end subroutine PrintMortarConnectivity 
!
!////////////////////////////////////////////////////////////////////////////////
!  
end module SlidingMeshProcedures