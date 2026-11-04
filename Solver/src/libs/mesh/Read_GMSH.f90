!//////////////////////////////////////////////////////
!
!  Module for reading hexahedral conforming meshes in GMSH (https://gmsh.info/) mesh format.
!  Current supported versions: GMSH Mesh Format 4.1 and 2.1
!  Current supported elements ID: 5,12,92, 93, 94. Respectively: P1, P2, P3, P4, P5 hexahedral.
!
!  GMSH original hexahedral ordering:
!
!         v
!  3----------2            3----13----2           3----13----2
!  |\     ^   |\           |\         |\          |\         |\
!  | \    |   | \          | 15       | 14        |15    24  | 14
!  |  \   |   |  \         9  \       11 \        9  \ 20    11 \
!  |   7------+---6        |   7----19+---6       |   7----19+---6
!  |   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
!  0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
!   \  |    \  \  |         \  17      \  18       \ 17    25 \  18
!    \ |     \  \ |         10 |        12|        10 |  21    12|
!     \|      w  \|           \|         \|          \|         \|
!      4----------5            4----16----5           4----16----5
!
!  HORSES3D corner ordering:
!         v
!  1----------4
!  |\     ^   |\
!  | \    |   | \
!  |  \   |   |  \
!  |   5------+---8
!  |   |  +-- |-- | -> u
!  2---+---\--3   |
!   \  |    \  \  |
!    \ |     \  \ |
!     \|      w  \|
!      6----------7
!
!////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
MODULE Read_GMSH
      use SMConstants
      use MeshTypes
      use ElementConnectivityDefinitions
      USE TransfiniteMapClass
      use MappedGeometryClass
      use NodeClass
      use ElementClass
      use HexMeshClass
      use sharedBCModule
      use PhysicsStorage
      use FileReadingUtilities      , only: getFileName, getRealArrayFromStringNoCommas
      use Utilities, only: UnusedUnit, toLower, my_findloc
      implicit none

      private
      public ConstructMesh_FromGMSHFile_v4_, ConstructMesh_FromGMSHFile_v2_, CheckGMSHversion, NumOfElems_GMSH_v4, NumOfElems_GMSH_v2

      public MSH_node_block_t, MSH_point_t, MSH_element_block_t, EL_MAX_ORDER, SUPPORTED_EL_TYPES, ReorderElement, MSH_LEN

!
!  ------------------------------------------------
!  Local temporary element storage.
!  ------------------------------------------------
!
   type :: MSH_BCinfo_t
!-----Variables-----------------------------------------------------------
      integer                             :: dim
      integer                             :: tag
      character(len=127)                  :: name
      integer                             :: no_of_elements
      integer                             :: no_of_nodes
      integer                             :: no_of_surfaces
      integer                             :: no_of_curves
      integer                             :: no_of_points
      integer, dimension(:), allocatable  :: element_tags
      integer, dimension(:), allocatable  :: wall_side
      integer, dimension(:), allocatable  :: node_tags
      integer, dimension(:), allocatable  :: volume_tags
      integer, dimension(:), allocatable  :: surface_tags
      integer, dimension(:), allocatable  :: curve_tags
      integer, dimension(:), allocatable  :: point_tags
      contains
!--   ---Subroutines-----------------------------------------------------------
      procedure                                  :: Destruct  => MSH_DestructBCStorage
   end type MSH_BCinfo_t

   type :: MSH_point_t
!-----Variables-----------------------------------------------------------
      integer             :: tag
      real(kind=RP)       :: x(3)
      integer             :: no_ptags
      integer             :: ptags(16)=0
   end type MSH_point_t

   type :: MSH_curve_t
!-----Variables-----------------------------------------------------------
      integer             :: tag
      real(kind=RP)       :: minX(3)
      real(kind=RP)       :: maxX(3)
      integer             :: no_ptags
      integer             :: ptags(16)=0
      integer             :: no_bps
      integer             :: bps(2)
   end type MSH_curve_t

   type :: MSH_surf_t
!-----Variables-----------------------------------------------------------
      integer             :: tag
      real(kind=RP)       :: minX(3)
      real(kind=RP)       :: maxX(3)
      integer             :: no_ptags
      integer             :: ptags(16)=0
      integer             :: no_bps
      integer             :: bps(64)
   end type MSH_surf_t

   type :: MSH_vol_t
!-----Variables-----------------------------------------------------------
      integer             :: tag
      real(kind=RP)       :: minX(3)
      real(kind=RP)       :: maxX(3)
      integer             :: no_ptags
      integer             :: ptags(16)=0
      integer             :: no_bps
      integer             :: bps(64)
   end type MSH_vol_t

   type :: MSH_node_block_t
!-----Variables-----------------------------------------------------------
      integer             :: entity_dim
      integer             :: entity_tag
      logical             :: parametric
      integer             :: no_nodes
      integer, dimension(:), allocatable         :: tags
      real(kind=RP), dimension(:,:), allocatable :: cords
      contains
!--   ---Subroutines-----------------------------------------------------------
      procedure                                  :: Construct => MSH_ConstructNodeBlock
      procedure                                  :: Destruct  => MSH_DestructNodeBlock
   end type MSH_node_block_t

   type :: MSH_element_block_t
!-----Variables-----------------------------------------------------------
      integer             :: el_type
      integer             :: no_els
      integer, dimension(:,:), allocatable   :: BCs
      integer, dimension(:), allocatable     :: tags
      integer, dimension(:,:), allocatable   :: nodes
      integer, dimension(:), allocatable   :: no_ptags
      integer, dimension(:,:), allocatable     :: ptags
      contains
!--   ---Subroutines-----------------------------------------------------------
      procedure                                  :: Construct => MSH_ConstructElementBlock
      procedure                                  :: Destruct  => MSH_DestructElementBlock
   end type MSH_element_block_t
!
   integer, parameter :: EL_MAX_ORDER = 5
   integer, parameter :: MSH_LEN = 4096
   integer, parameter :: SUPPORTED_EL_TYPES(5) = (/5,12,92,93,94/) ! GMSH HEX types, orders from 1 to 5
!     ========
      CONTAINS
!     ========
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CheckGMSHversion( fileName, ver )
!  -----------------------------------------------------------------------
!  Build mesh from GMSH file.
!  -----------------------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      character(len=*) , intent(in)    :: fileName
      integer          , intent(inout) :: ver
!-----Local-Variables-----------------------------------------------------
      integer                         :: fUnit, fileStat
      real(kind=RP)                   :: tmpd
      character(len=1024)             :: tmps
      integer                         :: tmpi
!-------------------------------------------------------------------------

!-----Check-mesh-file-----------------------------------------------------
      fUnit = UnusedUnit()
      OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
      IF ( fileStat /= 0 )     THEN
         error stop "CheckGMSHversion :: error opening mesh file."
      END IF
!-------------------------------------------------------------------------

!-----Read-header-info----------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$MeshFormat') error stop "READ_GMSH :: Wrong input file."
      read(fUnit,*) tmpd, tmpi, tmpi
!-------------------------------------------------------------------------

      ver = floor(tmpd) ! .msh version

      close(fUnit)

   end subroutine CheckGMSHversion
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ConstructMesh_FromGMSHFile_v4_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
!  ---------------------------------------------------------
!  Build mesh from GMSH file.
!  ---------------------------------------------------------
      use Physics
      use PartitionedMeshClass
      use MPI_Process_Info
      implicit none
!-----Arguments-----------------------------------------------------------
      type(HexMesh)                   :: self
      integer                         :: nodes
      character(len=*)                :: fileName
      integer                         :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
      integer                         :: dir2D
      logical                         :: periodRelative
      logical           , intent(out) :: success
!-----Local-Variables-----------------------------------------------------
      character(len=1024)             :: tmps
      real(kind=RP)                   :: tmpd
      integer                         :: tmpi, tmpi1, tmpi2, tmpi3, tmp_eltag
      logical                         :: tmpb
      integer, dimension(:), allocatable :: tmpi_vec1, tmpi_vec2
      character(len=127),    allocatable :: tmps_vec(:)
      type(MSH_BCinfo_t), dimension(:), allocatable :: msh_bcs
      integer                         :: msh_no_BCs
      integer                         :: msh_no_points
      integer                         :: msh_no_curves
      integer                         :: msh_no_surfaces
      integer                         :: msh_no_volumes
      type(MSH_point_t), dimension(:), allocatable  :: msh_points
      type(MSH_curve_t), dimension(:), allocatable  :: msh_curves
      type(MSH_surf_t) , dimension(:), allocatable  :: msh_surfaces
      type(MSH_vol_t)  , dimension(:), allocatable  :: msh_volumes
      integer                         :: msh_no_nodeblocks, msh_nodeblock, msh_nodes_per_block, msh_node, msh_global_node
      integer                    :: element_type, org_element_type
      integer                    :: msh_no_elblocks, msh_elblock, msh_els_per_block, msh_el, msh_global_el

      type(MSH_node_block_t)  , dimension(:), allocatable  :: msh_node_blocks
      type(MSH_element_block_t)  , dimension(:), allocatable  :: msh_element_blocks

      character(len=MSH_LEN) :: msh_entity
      real(kind=RP), allocatable :: msh_entity_vec(:)
      integer, dimension(EL_MAX_ORDER)      :: check_eltype

      integer                         :: numberOfElements
      integer                         :: numberOfNodes
      integer                         :: numberOfBoundaryFaces
      integer                         :: numberOfFaces

      integer                         :: bFaceOrder, numBFacePoints, innerEdgePoints
      integer                         :: i, j, k, l, jj, ii
      integer                         :: fUnit, fileStat
      integer                         :: nodeIDs(NODES_PER_ELEMENT)
      real(kind=RP)                   :: x(NDIM)
      CHARACTER(LEN=BC_STRING_LENGTH), pointer :: zoneNames(:)
      real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
      integer, allocatable            :: tags(:)
!-----Curved-patches------------------------------------------------------
      real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
      real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
!  -----------------------------------------------------------------------

!-----Check-if-a-mesh-partition-exists-----------------------------------
      if ( MPI_Process % doMPIAction ) then
         if ( mpi_partition % Constructed ) then
            call ConstructMeshPartition_FromGMSHFile_v4_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
         else
            call ConstructSimplestMesh_FromGMSHFile_v4_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
         end if
         return
      end if

      numberOfBoundaryFaces = 0
      success               = .TRUE.

      fUnit = UnusedUnit()
      OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
      IF ( fileStat /= 0 )     THEN
         PRINT *, "Error opening file: ", fileName
         success = .FALSE.
         RETURN
      END IF
!------------------------------------------------------------------------

!-----Read-header-info---------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$MeshFormat') error stop "READ_GMSH :: Wrong input file."
      read(fUnit,*) tmpd, tmpi, tmpi
      read(fUnit,*) tmps
!------------------------------------------------------------------------

!-----Read-BC-info-------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$PhysicalNames') error stop "READ_GMSH :: Wrong input file - no boundary conditions defined."
      read(fUnit,*) tmpi

      allocate(tmpi_vec1(tmpi))
      allocate(tmpi_vec2(tmpi))
      allocate(tmps_vec(tmpi))

      msh_no_BCs = 0 ! only surfaces count!
      do i=1, tmpi
         read(fUnit,*) tmpi_vec1(i), tmpi_vec2(i), tmps_vec(i)
         if(tmpi_vec1(i) .eq. 2) msh_no_BCs=msh_no_BCs+1 ! check if surface
      end do ! tmpi
      if (msh_no_BCs .eq. 0) print *, "READ_GMSH :: No boundary conditions detected."

      allocate(msh_bcs(msh_no_BCs))

      j = 0
      do i=1, tmpi
         if(tmpi_vec1(i) .eq. 2) then
            j = j + 1
            msh_bcs(j)%dim  = tmpi_vec1(i)
            msh_bcs(j)%tag  = tmpi_vec2(i)
            msh_bcs(j)%name = tmps_vec(i)
         end if
      end do ! msh_no_BCs

      deallocate(tmpi_vec1)
      deallocate(tmpi_vec2)
      deallocate(tmps_vec)

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndPhysicalNames') error stop "READ_GMSH :: Wrong input file - not all boundary conditions detected."
!------------------------------------------------------------------------

!-----Read-msh-entities--------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Entities') error stop "READ_GMSH :: Wrong input file - no entities found."
      read(fUnit,*) msh_no_points, msh_no_curves, msh_no_surfaces, msh_no_volumes

      ! allocate memory for gmsh internal entities
      allocate(msh_points(msh_no_points))
      allocate(msh_curves(msh_no_curves))
      allocate(msh_surfaces(msh_no_surfaces))
      allocate(msh_volumes(msh_no_volumes))
      allocate(msh_entity_vec(32)) ! arbitrary number

!-----Read-points--------------------------------------------------------
      do i=1, msh_no_points
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') msh_entity
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_points(i)%tag      = int(msh_entity_vec(1))
         msh_points(i)%x        = msh_entity_vec(2:4)
         msh_points(i)%no_ptags = int(msh_entity_vec(5))
         if (msh_points(i)%no_ptags .gt. 0) msh_points(i)%ptags(1:msh_points(i)%no_ptags) = int(msh_entity_vec(6:5+msh_points(i)%no_ptags))
      end do ! msh_no_points
!------------------------------------------------------------------------

!-----Read-curves--------------------------------------------------------
      do i=1, msh_no_curves
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') msh_entity
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_curves(i)%tag      = int(msh_entity_vec(1))
         msh_curves(i)%minX     = msh_entity_vec(2:4)
         msh_curves(i)%maxX     = msh_entity_vec(5:7)
         msh_curves(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_curves(i)%no_ptags .gt. 0) msh_curves(i)%ptags(1:msh_curves(i)%no_ptags) = int(msh_entity_vec(9:8+msh_curves(i)%no_ptags))
         msh_curves(i)%no_bps = msh_entity_vec(9+msh_curves(i)%no_ptags)
         msh_curves(i)%bps = 0
         msh_curves(i)%bps(1:msh_curves(i)%no_bps) = msh_entity_vec(10+msh_curves(i)%no_ptags:9+msh_curves(i)%no_ptags+msh_curves(i)%no_bps)
         msh_curves(i)%bps = abs(msh_curves(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_curves
!------------------------------------------------------------------------

!-----Read-surfaces------------------------------------------------------
      do i=1, msh_no_surfaces
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') msh_entity
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_surfaces(i)%tag      = int(msh_entity_vec(1))
         msh_surfaces(i)%minX     = msh_entity_vec(2:4)
         msh_surfaces(i)%maxX     = msh_entity_vec(5:7)
         msh_surfaces(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_surfaces(i)%no_ptags .gt. 0) msh_surfaces(i)%ptags(1:msh_surfaces(i)%no_ptags) = int(msh_entity_vec(9:8+msh_surfaces(i)%no_ptags))
         msh_surfaces(i)%no_bps = msh_entity_vec(9+msh_surfaces(i)%no_ptags)
         msh_surfaces(i)%bps = 0
         msh_surfaces(i)%bps(1:msh_surfaces(i)%no_bps) = msh_entity_vec(10+msh_surfaces(i)%no_ptags:9+msh_surfaces(i)%no_ptags+msh_surfaces(i)%no_bps)
         msh_surfaces(i)%bps = abs(msh_surfaces(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_surfaces
!------------------------------------------------------------------------

!-----Read-volumes-------------------------------------------------------
      do i=1, msh_no_volumes
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') msh_entity
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_volumes(i)%tag      = int(msh_entity_vec(1))
         msh_volumes(i)%minX     = msh_entity_vec(2:4)
         msh_volumes(i)%maxX     = msh_entity_vec(5:7)
         msh_volumes(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_volumes(i)%no_ptags .gt. 0) msh_volumes(i)%ptags(1:msh_volumes(i)%no_ptags) = int(msh_entity_vec(9:8+msh_volumes(i)%no_ptags))
         msh_volumes(i)%no_bps = msh_entity_vec(9+msh_volumes(i)%no_ptags)
         msh_volumes(i)%bps(1:msh_volumes(i)%no_bps) = msh_entity_vec(10+msh_volumes(i)%no_ptags:9+msh_volumes(i)%no_ptags+msh_volumes(i)%no_bps)
         msh_volumes(i)%bps = abs(msh_volumes(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_volumes
!------------------------------------------------------------------------

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndEntities') error stop "READ_GMSH :: Wrong input file - not all entities detected."
!------------------------------------------------------------------------

!-----Read-nodes---------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Nodes') error stop "READ_GMSH :: Wrong input file - no nodes found."

      read(fUnit,*) msh_no_nodeblocks, numberOfNodes, tmpi1, tmpi2
      if (numberOfNodes .ne. tmpi2) error stop "READ_gmsh :: Incoherent node numbering."

      ! allocate nodes storage
      allocate(msh_node_blocks(msh_no_nodeblocks))

      msh_global_node = 0
      do msh_nodeblock=1, msh_no_nodeblocks
         read(fUnit,*) tmpi1, tmpi2, tmpi3, msh_nodes_per_block
         tmpb=tmpi3
         if (tmpb) error stop "READ_gmsh :: Parametric nodes not supported."

         call msh_node_blocks(msh_nodeblock) % Construct(tmpi1, tmpi2, tmpb, msh_nodes_per_block)

         do msh_node=1, msh_nodes_per_block
            read(fUnit,*) msh_node_blocks(msh_nodeblock) % tags(msh_node)
         end do ! msh_nodes_per_block /for tags

         do msh_node=1, msh_nodes_per_block
            read(fUnit,*) msh_node_blocks(msh_nodeblock) % cords(msh_node,:)
         end do ! msh_nodes_per_block /for coordinates

         msh_global_node = msh_global_node + msh_nodes_per_block

      end do ! msh_no_nodeblocks

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndNodes') error stop "READ_GMSH :: Wrong input file - not all nodes detected."
!------------------------------------------------------------------------

!-----Read-elements------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Elements') error stop "READ_GMSH :: Wrong input file - no elements found."

      read(fUnit,*) msh_no_elblocks, numberOfElements, tmpi1, tmpi2
      if (numberOfElements .ne. tmpi2) error stop "READ_gmsh :: Incoherent element numbering."

      allocate(msh_element_blocks(msh_no_elblocks))

      msh_global_el = 0
      do msh_elblock=1, msh_no_elblocks

         read(fUnit,*) tmpi1, tmpi2, element_type, msh_els_per_block

         call msh_element_blocks(msh_elblock) % Construct(element_type,msh_els_per_block)

         do msh_el=1, msh_els_per_block
            read(fUnit,*) msh_element_blocks(msh_elblock) % tags(msh_el), msh_element_blocks(msh_elblock) % nodes(msh_el,:)
         end do ! msh_els_per_block /for tags

         msh_global_el = msh_global_el + msh_els_per_block

      end do ! msh_no_elblocks

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndElements') error stop "READ_GMSH :: Wrong input file - not all elements detected."
      close( fUnit )
!------------------------------------------------------------------------

!-----Mesh-info----------------------------------------------------------
      ! find order of elements curvature
      check_eltype = 0
      do i=1, EL_MAX_ORDER
         check_eltype(i) = count(msh_element_blocks(:) % el_type .eq.  SUPPORTED_EL_TYPES(i))
      end do
      if (sum(check_eltype) .eq. 0) error stop "READ_GMSH :: No 3D elements detected in the mesh."
      if (sum(check_eltype) .ne. maxval(check_eltype)) error stop "READ_GMSH :: More than 1 type of hexahedral detected in the mesh."
      bFaceOrder = maxloc(check_eltype,1) ! set order of the mesh
      org_element_type = SUPPORTED_EL_TYPES(bFaceOrder)
      ! find number of elements
      numberOfElements = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            numberOfElements = numberOfElements + msh_element_blocks(msh_elblock) % no_els
         end if
      end do
!------------------------------------------------------------------------

!-----Reorder-nodes-in-elements------------------------------------------
      allocate(tmpi_vec1((bFaceOrder + 1)**3))
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do j=1, msh_element_blocks(msh_elblock) % no_els
               ! re-order nodes within element to match HORSES ordering
               tmpi_vec1 = msh_element_blocks(msh_elblock) % nodes(j,:)
               call ReorderElement(tmpi_vec1,bFaceOrder)
               msh_element_blocks(msh_elblock) % nodes(j,:) = tmpi_vec1
            end do
         end if
      end do
      deallocate(tmpi_vec1)
!------------------------------------------------------------------------

!-----Assign-BCs---------------------------------------------------------
      do i=1,msh_no_BCs
         ! find no surfaces
         msh_bcs(i) % no_of_surfaces = 0
         do j=1, msh_no_surfaces
            tmpi = my_findloc(msh_surfaces(j)%ptags, msh_bcs(i)%tag, 1)
            if (tmpi .gt. 0) msh_bcs(i) % no_of_surfaces = msh_bcs(i) % no_of_surfaces + 1
         end do ! msh_no_surfaces

         ! assign surfaces and curves
         allocate(tmpi_vec1(256))
         tmpi_vec1 = 0
         allocate(msh_bcs(i) % surface_tags(msh_bcs(i) % no_of_surfaces))
         k = 0
         jj = 1
         do j=1, msh_no_surfaces
            tmpi = my_findloc(msh_surfaces(j)%ptags, msh_bcs(i)%tag, 1)
            if (tmpi .gt. 0) then
               k = k + 1
               msh_bcs(i) % surface_tags(k) = msh_surfaces(j) % tag
               tmpi_vec1(jj : jj + count(msh_surfaces(j) % bps /= 0) - 1) = msh_surfaces(j) % bps(1 : count(msh_surfaces(j) % bps /= 0))
               jj = jj + 1 + count(msh_surfaces(j) % bps /= 0)
            end if ! tmpi
         end do ! msh_no_surfaces
         call unique(tmpi_vec1,msh_bcs(i) % curve_tags)
         msh_bcs(i) % no_of_curves = size(msh_bcs(i) % curve_tags,1)
         deallocate(tmpi_vec1)

         ! assign points
         k = 0
         allocate(tmpi_vec1(2 * msh_bcs(i) % no_of_curves))
         tmpi_vec1 = 0
         do j=1, msh_no_curves
            tmpi = my_findloc(msh_bcs(i) % curve_tags, msh_curves(j)%tag, 1)
            if (tmpi .gt. 0) then
               k = k + 1
               tmpi_vec1(1 + (k-1)*2: k*2) = msh_curves(j) % bps
            end if ! tmpi
         end do ! msh_no_points
         call unique(tmpi_vec1,msh_bcs(i) % point_tags)
         msh_bcs(i) % no_of_points = size(msh_bcs(i) % point_tags,1)
         deallocate(tmpi_vec1)

         ! assign nodes
         k = 0
         allocate(tmpi_vec1(numberOfNodes)) ! approximate upper bound
         tmpi_vec1 = 0
         do j=1, msh_no_nodeblocks

            select case (msh_node_blocks(j) % entity_dim)
            case (0) ! points
               tmpi = my_findloc(msh_bcs(i) % point_tags, msh_node_blocks(j) % entity_tag, 1)
            case (1) ! curves
               tmpi = my_findloc(msh_bcs(i) % curve_tags, msh_node_blocks(j) % entity_tag, 1)
            case (2) ! surfaces
               tmpi = my_findloc(msh_bcs(i) % surface_tags, msh_node_blocks(j) % entity_tag, 1)
            case (3) ! volume
               tmpi = 0
            end select ! msh_node_blocks(j) % entity_dim

            if (tmpi .gt. 0) then
               tmpi_vec1(k+1 : k+msh_node_blocks(j) % no_nodes) = msh_node_blocks(j) % tags
               k = k + msh_node_blocks(j) % no_nodes
            end if ! tmpi

         end do ! msh_no_nodeblocks
         call unique(tmpi_vec1(1:k),msh_bcs(i) % node_tags)
         msh_bcs(i) % no_of_nodes = size(msh_bcs(i) % node_tags,1)
         deallocate(tmpi_vec1)

      end do ! msh_no_BCs
!------------------------------------------------------------------------

!-----Assign-BC-to-nodes-------------------------------------------------
      allocate(tmpi_vec1(8)) ! help array for the walls
      do i=1, msh_no_BCs
         do msh_elblock=1, msh_no_elblocks
            if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
               do j=1, msh_element_blocks(msh_elblock) % no_els

                  ! find matching nodes
                  tmpi_vec1=0
                  do l=1,8
                     tmpi = my_findloc(msh_bcs(i) % node_tags,msh_element_blocks(msh_elblock) % nodes(j,l),1)
                     if(tmpi .gt. 0) tmpi_vec1(l) = 1
                     ! TODO: Check why this doesn't work with older ifort
                     ! tmpi = any(msh_bcs(i) % node_tags .eq. msh_element_blocks(msh_elblock) % nodes(j,l))
                     ! tmpi_vec1(l)=tmpi
                  end do

                  ! assign BC face to the elements nodes
                  if ( (sum(tmpi_vec1(1:2)+tmpi_vec1(5:6))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,1) = msh_bcs(i) % tag ! 1:west
                  if ( (sum(tmpi_vec1(3:4)+tmpi_vec1(7:8))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,2) = msh_bcs(i) % tag ! 2:east
                  if ( (sum(tmpi_vec1(1:4))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,3) = msh_bcs(i) % tag ! 3:south
                  if ( (sum(tmpi_vec1(2:3)+tmpi_vec1(6:7))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,4) = msh_bcs(i) % tag ! 4:front
                  if ( (sum(tmpi_vec1(5:8))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,5) = msh_bcs(i) % tag ! 5:north
                  if ( (tmpi_vec1(1)+tmpi_vec1(4)+tmpi_vec1(5)+tmpi_vec1(8)) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,6) = msh_bcs(i) % tag ! 6:back

               end do ! msh_element_blocks(msh_elblock) % no_els
            end if
         end do ! msh_no_elblocks
      end do ! msh_no_BCs
      deallocate(tmpi_vec1)
!------------------------------------------------------------------------

!-----Assign-new-tags----------------------------------------------------
      tmp_eltag = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do j=1, msh_element_blocks(msh_elblock) % no_els
               tmp_eltag = tmp_eltag + 1
               msh_element_blocks(msh_elblock) % tags(j) = tmp_eltag
            end do
         end if
      end do
      if (.not. (tmp_eltag .eq. numberOfElements)) error stop "Read_GMSH :: Number of elements inconsistent."
!------------------------------------------------------------------------

!-----Build-nodes--------------------------------------------------------
      self % nodeType = nodes
      self % no_of_elements = numberOfElements
      self % no_of_allElements = numberOfElements
!------------------------------------------------------------------------

!-----Set-up-face-patches------------------------------------------------
      numBFacePoints = bFaceOrder + 1
      innerEdgePoints = bFaceOrder - 1
      allocate(uNodes(numBFacePoints))
      allocate(vNodes(numBFacePoints))
      allocate(values(3,numBFacePoints,numBFacePoints))
      do i = 1, numBFacePoints
         uNodes(i) = -1._RP + (i-1) * (2._RP/bFaceOrder)
         vNodes(i) = uNodes(i)
      end do
!------------------------------------------------------------------------

!-----Allocate-mem-for-elements-and-nodes--------------------------------
      allocate( self % elements(numberOfelements) )
      allocate( self % nodes(numberOfNodes) )
      allocate( self % Nx(numberOfelements) , self % Ny(numberOfelements) , self % Nz(numberOfelements) )
      self % Nx = Nx
      self % Ny = Ny
      self % Nz = Nz
!------------------------------------------------------------------------

!----Set-nodes-----------------------------------------------------------
      do msh_nodeblock=1, msh_no_nodeblocks
         do msh_node=1, msh_node_blocks(msh_nodeblock) % no_nodes
            x = msh_node_blocks(msh_nodeblock) % cords(msh_node,1:NDIM)/Lref
            call ConstructNode( self % nodes(msh_node_blocks(msh_nodeblock) % tags(msh_node)), x, msh_node_blocks(msh_nodeblock) % tags(msh_node) )
         end do ! msh_node_blocks(msh_nodeblock) % no_nodes
      end do ! msh_no_nodeblocks
!------------------------------------------------------------------------

!----Set-elements-----------------------------------------------------------
      j = 0
      call extract_bc_tags(msh_bcs, tags)
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do msh_el = 1, msh_element_blocks(msh_elblock) % no_els
               j = j + 1
               ! setting l'th element
               l =  msh_element_blocks(msh_elblock) % tags(msh_el)
               nodeIDs = msh_element_blocks(msh_elblock) % nodes(msh_el,1:8) ! only non-curved nodes

               if (bFaceOrder .eq. 1) then ! non-curved mesh
                  do k = 1, NODES_PER_ELEMENT
                     corners(:,k) = self % nodes(nodeIDs(k)) % x
                  end do
                  self % elements(l) % SurfInfo % IsHex8 = .TRUE.
                  self % elements(l) % SurfInfo % corners = corners
               else ! curved mesh
                  ! allocate tmp arrays for curved face node tags and coordinates
                  allocate(tmpi_vec1(numBFacePoints*numBFacePoints))
                  allocate(tmpi_vec2(numBFacePoints**3))

                  tmpi_vec2 = msh_element_blocks(msh_elblock) % nodes(msh_el,:)

                  do k = 1, FACES_PER_ELEMENT
                     call GetOrderedFaceNodeTags(tmpi_vec1,k,bFaceOrder,tmpi_vec2)
                     do jj = 1, numBFacePoints
                        do ii = 1, numBFacePoints
                           values(:,ii,jj) = self % nodes(tmpi_vec1(ii + (jj-1)*numBFacePoints)) % x
                        end do
                     end do
                     values = values / Lref
                     call self % elements(l) % SurfInfo % facePatches(k) % construct(uNodes, vNodes, values)

                  end do

                  deallocate(tmpi_vec1)
                  deallocate(tmpi_vec2)
               end if

               call self % elements(l) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , l, l)

               ! set element boundaries
               do k = 1, 6
                  tmpi1 = my_findloc( tags,msh_element_blocks(msh_elblock) % BCs(msh_el,k),1)
                  if (tmpi1 .gt. 0) then
                     self % elements(l) % boundaryName(k) = trim(msh_bcs(tmpi1) % name)
                  else
                     self % elements(l) % boundaryName(k) = emptyBCName
                  end if
               end do ! k

               ! set BC names to faces
               do k = 1, 6
                  if(trim(self % elements(l) % boundaryName(k)) /= emptyBCName) then
                     call toLower( self % elements(l) % boundaryName(k) )
                     numberOfBoundaryFaces = numberOfBoundaryFaces + 1
                     zoneNames => zoneNameDictionary % allKeys()
                     if ( all(trim(self % elements(l) % boundaryName(k)) .ne. zoneNames) ) then
                        call zoneNameDictionary % addValueForKey(trim(self % elements(l) % boundaryName(k)), trim(self % elements(l) % boundaryName(k)))
                     end if
                     deallocate (zoneNames)
                  end if
               end do ! k
            end do ! msh_element_blocks(msh_elblock) % no_els
         end if ! if el_type .eq. org_element_type
      end do ! msh_no_elblocks
      if (.not. (j .eq. numberOfElements)) error stop "Read_GMSH :: Not all elements assigned."
!------------------------------------------------------------------------

!-----Deallocate-msh-vairables-------------------------------------------
      do msh_nodeblock=1, msh_no_nodeblocks
         call msh_node_blocks(msh_nodeblock) % Destruct()
      end do
      deallocate(msh_node_blocks)

      do msh_elblock=1, msh_no_elblocks
         call msh_element_blocks(msh_elblock) % Destruct()
      end do
      deallocate(msh_element_blocks)

      do i=1,msh_no_BCs
         call msh_bcs(i) % Destruct()
      end do
      deallocate(msh_bcs)

      deallocate(msh_points)
      deallocate(msh_curves)
      deallocate(msh_surfaces)
      deallocate(msh_volumes)
      deallocate(msh_entity_vec)

      deallocate(uNodes) ! Check if we can do it! FIXME
      deallocate(vNodes) ! Check if we can do it! FIXME
      deallocate(values) ! Check if we can do it! FIXME
!------------------------------------------------------------------------
!
!     ---------------------------
!     Construct the element faces
!     ---------------------------
!
      numberOfFaces        = (6*numberOfElements + numberOfBoundaryFaces)/2
      self % numberOfFaces = numberOfFaces
      allocate( self % faces(self % numberOfFaces) )
      CALL ConstructFaces( self, success )
!
!     -------------------------
!     Build the different zones
!     -------------------------
!
      call self % ConstructZones()
!
!     ---------------------------
!     Construct periodic faces
!     ---------------------------
!
      CALL ConstructPeriodicFaces( self, periodRelative )
!
!     ---------------------------
!     Delete periodic- faces
!     ---------------------------
!
      CALL DeletePeriodicMinusFaces( self )
!
!     ---------------------------
!     Assign faces ID to elements
!     ---------------------------
!
      CALL getElementsFaceIDs(self)
!
!     ---------------------
!     Define boundary faces
!     ---------------------
!
      call self % DefineAsBoundaryFaces()
!
!     -----------------------------------
!     Check if this is a 2D extruded mesh
!     -----------------------------------
!
      call self % CheckIfMeshIs2D()
!
!     -------------------------------
!     Set the mesh as 2D if requested
!     -------------------------------
!
      self % dir2D_ctrl = dir2D
      if ( dir2D .ne. 0 ) then
         call SetMappingsToCrossProduct
         call self % CorrectOrderFor2DMesh(dir2D,0)
      end if
!
!     ------------------------------
!     Set the element connectivities
!     ------------------------------
!
      call self % SetConnectivitiesAndLinkFaces(nodes)
!
!     ---------------------------------------
!     Construct elements' and faces' geometry
!     ---------------------------------------
!
      call self % ConstructGeometry()
!
!     -------------------------------
!     Set the mesh as 2D if requested
!     -------------------------------
!
      if ( dir2D .ne. 0 ) then
         call self % CorrectOrderFor2DMesh(dir2D,0)
      end if
!
!     ---------------------------------
!     Describe mesh and prepare for I/O
!     ---------------------------------
!
      if (.not. self % child) CALL self % Describe( trim(fileName), bFaceOrder )
      call self % PrepareForIO

      call self % ExportBoundaryMesh (trim(fileName))

   end subroutine ConstructMesh_FromGMSHFile_v4_
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------
!  Constructor of mesh partitions
!  ------------------------------
   SUBROUTINE ConstructMeshPartition_FromGMSHFile_v4_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
!  ---------------------------------------------------------
!  Build mesh from GMSH file.
!  ---------------------------------------------------------
      USE Physics
      use PartitionedMeshClass
      use MPI_Process_Info
      use MPI_Face_Class
      implicit none
!-----Arguments-----------------------------------------------------------
      type(HexMesh)                   :: self
      integer                         :: nodes
      character(len=*)                :: fileName
      integer                         :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
      integer                         :: dir2D
      logical                         :: periodRelative
      logical           , intent(out) :: success
!-----Local-Variables-----------------------------------------------------
      character(len=1024)             :: tmps
      real(kind=RP)                   :: tmpd
      integer                         :: tmpi, tmpi1, tmpi2, tmpi3, tmp_eltag
      logical                         :: tmpb
      integer, dimension(:), allocatable :: tmpi_vec1, tmpi_vec2
      character(len=127),    allocatable :: tmps_vec(:)
      type(MSH_BCinfo_t), dimension(:), allocatable :: msh_bcs
      integer                         :: msh_no_BCs
      integer                         :: msh_no_points
      integer                         :: msh_no_curves
      integer                         :: msh_no_surfaces
      integer                         :: msh_no_volumes
      type(MSH_point_t), dimension(:), allocatable  :: msh_points
      type(MSH_curve_t), dimension(:), allocatable  :: msh_curves
      type(MSH_surf_t) , dimension(:), allocatable  :: msh_surfaces
      type(MSH_vol_t)  , dimension(:), allocatable  :: msh_volumes
      integer                         :: msh_no_nodeblocks, msh_nodeblock, msh_nodes_per_block, msh_node, msh_global_node
      integer                    :: element_type, org_element_type
      integer                    :: msh_no_elblocks, msh_elblock, msh_els_per_block, msh_el, msh_global_el

      type(MSH_node_block_t)  , dimension(:), allocatable  :: msh_node_blocks
      type(MSH_element_block_t)  , dimension(:), allocatable  :: msh_element_blocks

      type(Node)  , dimension(:), allocatable  :: local_nodes

      character(len=MSH_LEN) :: msh_entity
      real(kind=RP), allocatable :: msh_entity_vec(:)
      integer, dimension(EL_MAX_ORDER)      :: check_eltype

      integer                         :: numberOfAllElements
      integer                         :: numberOfAllNodes
      integer, allocatable            :: globalToLocalNodeID(:)
      integer, allocatable            :: globalToLocalElementID(:)

      integer                         :: numberOfElements
      integer                         :: numberOfNodes
      integer                         :: numberOfFaces

      integer                         :: bFaceOrder, numBFacePoints
      integer                         :: i, j, k, l, pNode, pElement
      integer                         :: jj, ii
      integer                         :: fUnit, fileStat
      integer                         :: nodeIDs(NODES_PER_ELEMENT)
      real(kind=RP)                   :: x(NDIM)
      CHARACTER(LEN=BC_STRING_LENGTH) :: names(FACES_PER_ELEMENT)
      CHARACTER(LEN=BC_STRING_LENGTH), pointer :: zoneNames(:)
      real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
      integer, allocatable            :: tags(:)
!-----Curved-patches------------------------------------------------------
      real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
      real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
!  -----------------------------------------------------------------------

      success               = .TRUE.

      fUnit = UnusedUnit()
      OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
      IF ( fileStat /= 0 )     THEN
         PRINT *, "Error opening file: ", fileName
         success = .FALSE.
         RETURN
      END IF
!------------------------------------------------------------------------

!-----Read-header-info---------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$MeshFormat') error stop "READ_GMSH :: Wrong input file."
      read(fUnit,*) tmpd, tmpi, tmpi
      read(fUnit,*) tmps
!------------------------------------------------------------------------

!-----Read-BC-info-------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$PhysicalNames') error stop "READ_GMSH :: Wrong input file - no boundary conditions defined."
      read(fUnit,*) tmpi

      allocate(tmpi_vec1(tmpi))
      allocate(tmpi_vec2(tmpi))
      allocate(tmps_vec(tmpi))

      msh_no_BCs = 0 ! only surfaces count!
      do i=1, tmpi
         read(fUnit,*) tmpi_vec1(i), tmpi_vec2(i), tmps_vec(i)
         if(tmpi_vec1(i) .eq. 2) msh_no_BCs=msh_no_BCs+1 ! check if surface
      end do ! tmpi
      if (msh_no_BCs .eq. 0) print *, "READ_GMSH :: No boundary conditions detected."

      allocate(msh_bcs(msh_no_BCs))

      j = 0
      do i=1, tmpi
         if(tmpi_vec1(i) .eq. 2) then
            j = j + 1
            msh_bcs(j)%dim  = tmpi_vec1(i)
            msh_bcs(j)%tag  = tmpi_vec2(i)
            msh_bcs(j)%name = tmps_vec(i)
         end if
      end do ! msh_no_BCs

      deallocate(tmpi_vec1)
      deallocate(tmpi_vec2)
      deallocate(tmps_vec)

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndPhysicalNames') error stop "READ_GMSH :: Wrong input file - not all boundary conditions detected."
!------------------------------------------------------------------------

!-----Read-msh-entities--------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Entities') error stop "READ_GMSH :: Wrong input file - no entities found."
      read(fUnit,*) msh_no_points, msh_no_curves, msh_no_surfaces, msh_no_volumes

      ! allocate memory for gmsh internal entities
      allocate(msh_points(msh_no_points))
      allocate(msh_curves(msh_no_curves))
      allocate(msh_surfaces(msh_no_surfaces))
      allocate(msh_volumes(msh_no_volumes))
      allocate(msh_entity_vec(32)) ! arbitrary number

!-----Read-points--------------------------------------------------------
      do i=1, msh_no_points
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') msh_entity
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_points(i)%tag      = int(msh_entity_vec(1))
         msh_points(i)%x        = msh_entity_vec(2:4)
         msh_points(i)%no_ptags = int(msh_entity_vec(5))
         if (msh_points(i)%no_ptags .gt. 0) msh_points(i)%ptags(1:msh_points(i)%no_ptags) = int(msh_entity_vec(6:5+msh_points(i)%no_ptags))
      end do ! msh_no_points
!------------------------------------------------------------------------

!-----Read-curves--------------------------------------------------------
      do i=1, msh_no_curves
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') msh_entity
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_curves(i)%tag      = int(msh_entity_vec(1))
         msh_curves(i)%minX     = msh_entity_vec(2:4)
         msh_curves(i)%maxX     = msh_entity_vec(5:7)
         msh_curves(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_curves(i)%no_ptags .gt. 0) msh_curves(i)%ptags(1:msh_curves(i)%no_ptags) = int(msh_entity_vec(9:8+msh_curves(i)%no_ptags))
         msh_curves(i)%no_bps = msh_entity_vec(9+msh_curves(i)%no_ptags)
         msh_curves(i)%bps = 0
         msh_curves(i)%bps(1:msh_curves(i)%no_bps) = msh_entity_vec(10+msh_curves(i)%no_ptags:9+msh_curves(i)%no_ptags+msh_curves(i)%no_bps)
         msh_curves(i)%bps = abs(msh_curves(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_curves
!------------------------------------------------------------------------

!-----Read-surfaces------------------------------------------------------
      do i=1, msh_no_surfaces
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') msh_entity
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_surfaces(i)%tag      = int(msh_entity_vec(1))
         msh_surfaces(i)%minX     = msh_entity_vec(2:4)
         msh_surfaces(i)%maxX     = msh_entity_vec(5:7)
         msh_surfaces(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_surfaces(i)%no_ptags .gt. 0) msh_surfaces(i)%ptags(1:msh_surfaces(i)%no_ptags) = int(msh_entity_vec(9:8+msh_surfaces(i)%no_ptags))
         msh_surfaces(i)%no_bps = msh_entity_vec(9+msh_surfaces(i)%no_ptags)
         msh_surfaces(i)%bps = 0
         msh_surfaces(i)%bps(1:msh_surfaces(i)%no_bps) = msh_entity_vec(10+msh_surfaces(i)%no_ptags:9+msh_surfaces(i)%no_ptags+msh_surfaces(i)%no_bps)
         msh_surfaces(i)%bps = abs(msh_surfaces(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_surfaces
!------------------------------------------------------------------------

!-----Read-volumes-------------------------------------------------------
      do i=1, msh_no_volumes
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') msh_entity
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_volumes(i)%tag      = int(msh_entity_vec(1))
         msh_volumes(i)%minX     = msh_entity_vec(2:4)
         msh_volumes(i)%maxX     = msh_entity_vec(5:7)
         msh_volumes(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_volumes(i)%no_ptags .gt. 0) msh_volumes(i)%ptags(1:msh_volumes(i)%no_ptags) = int(msh_entity_vec(9:8+msh_volumes(i)%no_ptags))
         msh_volumes(i)%no_bps = msh_entity_vec(9+msh_volumes(i)%no_ptags)
         msh_volumes(i)%bps(1:msh_volumes(i)%no_bps) = msh_entity_vec(10+msh_volumes(i)%no_ptags:9+msh_volumes(i)%no_ptags+msh_volumes(i)%no_bps)
         msh_volumes(i)%bps = abs(msh_volumes(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_volumes
!------------------------------------------------------------------------

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndEntities') error stop "READ_GMSH :: Wrong input file - not all entities detected."
!------------------------------------------------------------------------

!-----Read-nodes---------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Nodes') error stop "READ_GMSH :: Wrong input file - no nodes found."

      read(fUnit,*) msh_no_nodeblocks, numberOfNodes, tmpi1, tmpi2
      if (numberOfNodes .ne. tmpi2) error stop "READ_gmsh :: Incoherent node numbering."

      ! allocate nodes storage
      allocate(msh_node_blocks(msh_no_nodeblocks))

      msh_global_node = 0
      do msh_nodeblock=1, msh_no_nodeblocks
         read(fUnit,*) tmpi1, tmpi2, tmpi3, msh_nodes_per_block
         tmpb=tmpi3
         if (tmpb) error stop "READ_gmsh :: Parametric nodes not supported."

         call msh_node_blocks(msh_nodeblock) % Construct(tmpi1, tmpi2, tmpb, msh_nodes_per_block)

         do msh_node=1, msh_nodes_per_block
            read(fUnit,*) msh_node_blocks(msh_nodeblock) % tags(msh_node)
         end do ! msh_nodes_per_block /for tags

         do msh_node=1, msh_nodes_per_block
            read(fUnit,*) msh_node_blocks(msh_nodeblock) % cords(msh_node,:)
         end do ! msh_nodes_per_block /for coordinates

         msh_global_node = msh_global_node + msh_nodes_per_block

      end do ! msh_no_nodeblocks

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndNodes') error stop "READ_GMSH :: Wrong input file - not all nodes detected."
      numberOfAllNodes = numberOfNodes
!------------------------------------------------------------------------

!-----Read-elements------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Elements') error stop "READ_GMSH :: Wrong input file - no elements found."

      read(fUnit,*) msh_no_elblocks, numberOfElements, tmpi1, tmpi2
      if (numberOfElements .ne. tmpi2) error stop "READ_gmsh :: Incoherent element numbering."

      allocate(msh_element_blocks(msh_no_elblocks))

      msh_global_el = 0
      do msh_elblock=1, msh_no_elblocks

         read(fUnit,*) tmpi1, tmpi2, element_type, msh_els_per_block

         call msh_element_blocks(msh_elblock) % Construct(element_type,msh_els_per_block)

         do msh_el=1, msh_els_per_block
            read(fUnit,*) msh_element_blocks(msh_elblock) % tags(msh_el), msh_element_blocks(msh_elblock) % nodes(msh_el,:)
         end do ! msh_els_per_block /for tags

         msh_global_el = msh_global_el + msh_els_per_block

      end do ! msh_no_elblocks

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndElements') error stop "READ_GMSH :: Wrong input file - not all elements detected."
      close( fUnit )
!------------------------------------------------------------------------

!-----Mesh-info----------------------------------------------------------
      ! find order of elements curvature
      check_eltype = 0
      do i=1, EL_MAX_ORDER
         check_eltype(i) = count(msh_element_blocks(:) % el_type .eq.  SUPPORTED_EL_TYPES(i))
      end do
      if (sum(check_eltype) .eq. 0) error stop "READ_GMSH :: No 3D elements detected in the mesh."
      if (sum(check_eltype) .ne. maxval(check_eltype)) error stop "READ_GMSH :: More than 1 type of hexahedral detected in the mesh."
      bFaceOrder = maxloc(check_eltype,1) ! set order of the mesh
      org_element_type = SUPPORTED_EL_TYPES(bFaceOrder)
      ! find number of elements
      numberOfElements = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            numberOfElements = numberOfElements + msh_element_blocks(msh_elblock) % no_els
         end if
      end do
      numberOfAllElements = numberOfElements
!------------------------------------------------------------------------

!-----Reorder-nodes-in-elements------------------------------------------
      allocate(tmpi_vec1((bFaceOrder + 1)**3))
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do j=1, msh_element_blocks(msh_elblock) % no_els
               ! re-order nodes within element to match HORSES ordering
               tmpi_vec1 = msh_element_blocks(msh_elblock) % nodes(j,:)
               call ReorderElement(tmpi_vec1,bFaceOrder)
               msh_element_blocks(msh_elblock) % nodes(j,:) = tmpi_vec1
            end do
         end if
      end do
      deallocate(tmpi_vec1)
!------------------------------------------------------------------------

!-----Assign-BCs---------------------------------------------------------
      do i=1,msh_no_BCs
         ! find no surfaces
         msh_bcs(i) % no_of_surfaces = 0
         do j=1, msh_no_surfaces
            tmpi = my_findloc(msh_surfaces(j)%ptags, msh_bcs(i)%tag, 1)
            if (tmpi .gt. 0) msh_bcs(i) % no_of_surfaces = msh_bcs(i) % no_of_surfaces + 1
         end do ! msh_no_surfaces

         ! assign surfaces and curves
         allocate(tmpi_vec1(256))
         tmpi_vec1 = 0
         allocate(msh_bcs(i) % surface_tags(msh_bcs(i) % no_of_surfaces))
         k = 0
         jj = 1
         do j=1, msh_no_surfaces
            tmpi = my_findloc(msh_surfaces(j)%ptags, msh_bcs(i)%tag, 1)
            if (tmpi .gt. 0) then
               k = k + 1
               msh_bcs(i) % surface_tags(k) = msh_surfaces(j) % tag
               tmpi_vec1(jj : jj + count(msh_surfaces(j) % bps /= 0) - 1) = msh_surfaces(j) % bps(1 : count(msh_surfaces(j) % bps /= 0))
               jj = jj + 1 + count(msh_surfaces(j) % bps /= 0)
            end if ! tmpi
         end do ! msh_no_surfaces
         call unique(tmpi_vec1,msh_bcs(i) % curve_tags)
         msh_bcs(i) % no_of_curves = size(msh_bcs(i) % curve_tags,1)
         deallocate(tmpi_vec1)

         ! assign points
         k = 0
         allocate(tmpi_vec1(2 * msh_bcs(i) % no_of_curves))
         tmpi_vec1 = 0
         do j=1, msh_no_curves
            tmpi = my_findloc(msh_bcs(i) % curve_tags, msh_curves(j)%tag, 1)
            if (tmpi .gt. 0) then
               k = k + 1
               tmpi_vec1(1 + (k-1)*2: k*2) = msh_curves(j) % bps
            end if ! tmpi
         end do ! msh_no_points
         call unique(tmpi_vec1,msh_bcs(i) % point_tags)
         msh_bcs(i) % no_of_points = size(msh_bcs(i) % point_tags,1)
         deallocate(tmpi_vec1)

         ! assign nodes
         k = 0
         allocate(tmpi_vec1(numberOfNodes)) ! approximate upper bound
         tmpi_vec1 = 0
         do j=1, msh_no_nodeblocks

            select case (msh_node_blocks(j) % entity_dim)
            case (0) ! points
               tmpi = my_findloc(msh_bcs(i) % point_tags, msh_node_blocks(j) % entity_tag, 1)
            case (1) ! curves
               tmpi = my_findloc(msh_bcs(i) % curve_tags, msh_node_blocks(j) % entity_tag, 1)
            case (2) ! surfaces
               tmpi = my_findloc(msh_bcs(i) % surface_tags, msh_node_blocks(j) % entity_tag, 1)
            case (3) ! volume
               tmpi = 0
            end select ! msh_node_blocks(j) % entity_dim

            if (tmpi .gt. 0) then
               tmpi_vec1(k+1 : k+msh_node_blocks(j) % no_nodes) = msh_node_blocks(j) % tags
               k = k + msh_node_blocks(j) % no_nodes
            end if ! tmpi

         end do ! msh_no_nodeblocks
         call unique(tmpi_vec1(1:k),msh_bcs(i) % node_tags)
         msh_bcs(i) % no_of_nodes = size(msh_bcs(i) % node_tags,1)
         deallocate(tmpi_vec1)

      end do ! msh_no_BCs
!------------------------------------------------------------------------

!-----Assign-BC-to-nodes-------------------------------------------------
      allocate(tmpi_vec1(8)) ! help array for the walls
      do i=1, msh_no_BCs
         do msh_elblock=1, msh_no_elblocks
            if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
               do j=1, msh_element_blocks(msh_elblock) % no_els

                  ! find matching nodes
                  tmpi_vec1=0
                  do l=1,8
                     tmpi = my_findloc(msh_bcs(i) % node_tags,msh_element_blocks(msh_elblock) % nodes(j,l),1)
                     if(tmpi .gt. 0) tmpi_vec1(l) = 1
                     ! TODO: Check why this doesn't work with older ifort
                     ! tmpi = any(msh_bcs(i) % node_tags .eq. msh_element_blocks(msh_elblock) % nodes(j,l))
                     ! tmpi_vec1(l)=tmpi
                  end do

                  ! assign BC face to the elements nodes
                  if ( (sum(tmpi_vec1(1:2)+tmpi_vec1(5:6))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,1) = msh_bcs(i) % tag ! 1:west
                  if ( (sum(tmpi_vec1(3:4)+tmpi_vec1(7:8))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,2) = msh_bcs(i) % tag ! 2:east
                  if ( (sum(tmpi_vec1(1:4))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,3) = msh_bcs(i) % tag ! 3:south
                  if ( (sum(tmpi_vec1(2:3)+tmpi_vec1(6:7))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,4) = msh_bcs(i) % tag ! 4:front
                  if ( (sum(tmpi_vec1(5:8))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,5) = msh_bcs(i) % tag ! 5:north
                  if ( (tmpi_vec1(1)+tmpi_vec1(4)+tmpi_vec1(5)+tmpi_vec1(8)) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,6) = msh_bcs(i) % tag ! 6:back

               end do ! msh_element_blocks(msh_elblock) % no_els
            end if
         end do ! msh_no_elblocks
      end do ! msh_no_BCs
      deallocate(tmpi_vec1)
!------------------------------------------------------------------------

!-----Assign-new-tags----------------------------------------------------
      tmp_eltag = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do j=1, msh_element_blocks(msh_elblock) % no_els
               tmp_eltag = tmp_eltag + 1
               msh_element_blocks(msh_elblock) % tags(j) = tmp_eltag
            end do
         end if
      end do
      if (.not. (tmp_eltag .eq. numberOfElements)) error stop "Read_GMSH :: Number of elements inconsistent."
!------------------------------------------------------------------------

!-----Build-nodes--------------------------------------------------------
      self % nodeType = nodes
      self % no_of_elements = mpi_partition % no_of_elements
      self % no_of_allElements = numberOfAllElements
!------------------------------------------------------------------------

!-----Set-up-face-patches------------------------------------------------
      numBFacePoints = bFaceOrder + 1
      allocate(uNodes(numBFacePoints))
      allocate(vNodes(numBFacePoints))
      allocate(values(3,numBFacePoints,numBFacePoints))
      do i = 1, numBFacePoints
         uNodes(i) = -1._RP + (i-1) * (2._RP/bFaceOrder)
         vNodes(i) = uNodes(i)
      end do
!------------------------------------------------------------------------

!-----Allocate-mem-for-elements-and-nodes--------------------------------
      allocate( self % elements(mpi_partition % no_of_elements) )
      allocate( self % nodes(mpi_partition % no_of_nodes) )
      allocate ( self % Nx(self % no_of_elements) , self % Ny(self % no_of_elements) , self % Nz(self % no_of_elements) )
      allocate( globalToLocalNodeID(numberOfAllNodes) )
      allocate( globalToLocalElementID(numberOfAllElements) )

      globalToLocalNodeID = -1
      globalToLocalElementID = -1
!------------------------------------------------------------------------

!----Set-nodes-----------------------------------------------------------
      pNode = 1
      do msh_nodeblock=1, msh_no_nodeblocks
         do msh_node=1, msh_node_blocks(msh_nodeblock) % no_nodes
            x = msh_node_blocks(msh_nodeblock) % cords(msh_node,1:NDIM)/Lref

            if ( pNode .gt. mpi_partition % no_of_nodes ) cycle

            ! Construct only nodes that belong to the partition
            if ( msh_node_blocks(msh_nodeblock) % tags(msh_node) .eq. mpi_partition % nodeIDs(pNode) ) then
               call ConstructNode( self % nodes(pNode), x, msh_node_blocks(msh_nodeblock) % tags(msh_node) )
               globalToLocalNodeID(msh_node_blocks(msh_nodeblock) % tags(msh_node)) = pNode
               pNode = pNode + 1
            end if
         end do ! msh_node_blocks(msh_nodeblock) % no_nodes
      end do ! msh_no_nodeblocks
!------------------------------------------------------------------------

!----Set-local-nodes-----------------------------------------------------
      allocate(local_nodes(numberOfAllNodes))
      do msh_nodeblock=1, msh_no_nodeblocks
         do msh_node=1, msh_node_blocks(msh_nodeblock) % no_nodes
            x = msh_node_blocks(msh_nodeblock) % cords(msh_node,1:NDIM)/Lref
            call ConstructNode( local_nodes(msh_node_blocks(msh_nodeblock) % tags(msh_node)), x, msh_node_blocks(msh_nodeblock) % tags(msh_node) )
         end do ! msh_node_blocks(msh_nodeblock) % no_nodes
      end do ! msh_no_nodeblocks
!------------------------------------------------------------------------

!----Set-elements-----------------------------------------------------------
      j = 0
      pElement = 1
      call extract_bc_tags(msh_bcs, tags)
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do msh_el = 1, msh_element_blocks(msh_elblock) % no_els
               j = j + 1
               ! setting l'th element
               l =  msh_element_blocks(msh_elblock) % tags(msh_el)
               nodeIDs = msh_element_blocks(msh_elblock) % nodes(msh_el,1:8) ! only non-curved nodes
               nodeIDs = globalToLocalNodeID(nodeIDs)

               if ( pElement .gt. mpi_partition % no_of_elements ) then

                  ! set element boundaries
                  do k = 1, 6
                     tmpi1 = my_findloc(tags ,msh_element_blocks(msh_elblock) % BCs(msh_el,k),1)
                     if (tmpi1 .gt. 0) then
                        names(k) = trim(msh_bcs(tmpi1) % name)
                     else
                        names(k) = emptyBCName
                     end if
                  end do ! k

                  ! set BC names to faces
                  do k = 1, 6
                     if(trim(names(k)) /= emptyBCName) then
                        call toLower( names(k) )
                        zoneNames => zoneNameDictionary % allKeys()
                        if ( all(trim(names(k)) .ne. zoneNames) ) then
                           call zoneNameDictionary % addValueForKey(trim(names(k)), trim(names(k)))
                        end if
                        deallocate (zoneNames)
                     end if
                  end do ! k

                  cycle
               else if ( l .ne. mpi_partition % elementIDs(pElement) ) then

                  ! set element boundaries
                  do k = 1, 6
                     tmpi1 = my_findloc( tags,msh_element_blocks(msh_elblock) % BCs(msh_el,k),1)
                     if (tmpi1 .gt. 0) then
                        names(k) = trim(msh_bcs(tmpi1) % name)
                     else
                        names(k) = emptyBCName
                     end if
                  end do ! k

                  ! set BC names to faces
                  do k = 1, 6
                     if(trim(names(k)) /= emptyBCName) then
                        call toLower( names(k) )
                        zoneNames => zoneNameDictionary % allKeys()
                        if ( all(trim(names(k)) .ne. zoneNames) ) then
                           call zoneNameDictionary % addValueForKey(trim(names(k)), trim(names(k)))
                        end if
                        deallocate (zoneNames)
                     end if
                  end do ! k

                  cycle
               end if

               if (bFaceOrder .eq. 1) then ! non-curved mesh
                  do k = 1, NODES_PER_ELEMENT
                     corners(:,k) = self % nodes(nodeIDs(k)) % x
                  end do
                  self % elements(pElement) % SurfInfo % IsHex8 = .TRUE.
                  self % elements(pElement) % SurfInfo % corners = corners
               else ! curved mesh
                  ! allocate tmp arrays for curved face node tags and coordinates
                  allocate(tmpi_vec1(numBFacePoints*numBFacePoints))
                  allocate(tmpi_vec2(numBFacePoints**3))
                  tmpi_vec2 = msh_element_blocks(msh_elblock) % nodes(msh_el,:)

                  do k = 1, FACES_PER_ELEMENT
                     call GetOrderedFaceNodeTags(tmpi_vec1,k,bFaceOrder,tmpi_vec2)
                     do jj = 1, numBFacePoints
                        do ii = 1, numBFacePoints
                           values(:,ii,jj) = local_nodes(tmpi_vec1(ii + (jj-1)*numBFacePoints)) % x
                        end do
                     end do
                     values = values / Lref
                     call self % elements(pElement) % SurfInfo % facePatches(k) % construct(uNodes, vNodes, values)
                  end do

                  deallocate(tmpi_vec1)
                  deallocate(tmpi_vec2)
               end if

               call self % elements(pElement) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , pElement, l)

               self % Nx(pElement) = Nx(l)
               self % Ny(pElement) = Ny(l)
               self % Nz(pElement) = Nz(l)

               ! set element boundaries
               do k = 1, 6
                  tmpi1 = my_findloc( tags,msh_element_blocks(msh_elblock) % BCs(msh_el,k),1)
                  if (tmpi1 .gt. 0) then
                     self % elements(pElement) % boundaryName(k) = trim(msh_bcs(tmpi1) % name)
                  else
                     self % elements(pElement) % boundaryName(k) = emptyBCName
                  end if
               end do ! k

               ! set BC names to faces
               do k = 1, 6
                  if(trim(self % elements(pElement) % boundaryName(k)) /= emptyBCName) then
                     call toLower( self % elements(pElement) % boundaryName(k) )
                     zoneNames => zoneNameDictionary % allKeys()
                     if ( all(trim(self % elements(pElement) % boundaryName(k)) .ne. zoneNames) ) then
                        call zoneNameDictionary % addValueForKey(trim(self % elements(pElement) % boundaryName(k)), trim(self % elements(pElement) % boundaryName(k)))
                     end if
                     deallocate (zoneNames)
                  end if
               end do ! k

               globalToLocalElementID(l) = pElement
               pElement = pElement + 1
            end do ! msh_element_blocks(msh_elblock) % no_els
         end if ! if el_type .eq. org_element_type
      end do ! msh_no_elblocks
!------------------------------------------------------------------------

!-----Deallocate-msh-vairables-------------------------------------------
      do msh_nodeblock=1, msh_no_nodeblocks
         call msh_node_blocks(msh_nodeblock) % Destruct()
      end do
      deallocate(msh_node_blocks)

      do msh_elblock=1, msh_no_elblocks
         call msh_element_blocks(msh_elblock) % Destruct()
      end do
      deallocate(msh_element_blocks)

      do i=1,msh_no_BCs
         call msh_bcs(i) % Destruct()
      end do
      deallocate(msh_bcs)

      deallocate(msh_points)
      deallocate(msh_curves)
      deallocate(msh_surfaces)
      deallocate(msh_volumes)
      deallocate(msh_entity_vec)

      deallocate(uNodes) ! Check if we can do it! FIXME
      deallocate(vNodes) ! Check if we can do it! FIXME
      deallocate(values) ! Check if we can do it! FIXME

      ! deallocate local nodes
      deallocate(local_nodes)
!------------------------------------------------------------------------
!
!     ---------------------------
!     Construct the element faces
!     ---------------------------
!
      numberOfFaces        = GetOriginalNumberOfFaces(self)
      self % numberOfFaces = numberOfFaces
      allocate( self % faces(self % numberOfFaces) )
      CALL ConstructFaces( self, success )
!     --------------------------------
!     Get actual mesh element face IDs
!     --------------------------------
!
      CALL getElementsFaceIDs(self)
!
!     --------------
!     Cast MPI faces
!     --------------
!
      call ConstructMPIFaces( self % MPIfaces )
      call self % UpdateFacesWithPartition(mpi_partition, &
                                           numberOfAllElements, &
                                           globalToLocalElementID)
!
!     -------------------------
!     Build the different zones
!     -------------------------
!
      call self % ConstructZones()
!
!     ---------------------------
!     Construct periodic faces
!     ---------------------------
!
      CALL ConstructPeriodicFaces( self, periodRelative )
!
!     ---------------------------
!     Delete periodic- faces
!     ---------------------------
!
      CALL DeletePeriodicMinusFaces( self )
!
!     ---------------------------
!     Assign faces ID to elements
!     ---------------------------
!
      CALL getElementsFaceIDs(self)
!
!     ---------------------
!     Define boundary faces
!     ---------------------
!
      call self % DefineAsBoundaryFaces()
!
!     -----------------------------------
!     Check if this is a 2D extruded mesh
!     -----------------------------------
!
      call self % CheckIfMeshIs2D()
!
!     -------------------------------
!     Set the mesh as 2D if requested
!     -------------------------------
!
      if ( dir2D .ne. 0 ) then
         call SetMappingsToCrossProduct
         call self % CorrectOrderFor2DMesh(dir2D,0)
      end if
!
!     ------------------------------
!     Set the element connectivities
!     ------------------------------
!
      call self % SetConnectivitiesAndLinkFaces(nodes)
!
!     ---------------------------------------
!     Construct elements' and faces' geometry
!     ---------------------------------------
!
      call self % ConstructGeometry()

      if ( dir2D .ne. 0 ) then
         call self % CorrectOrderFor2DMesh(dir2D,0)
      end if
!
!     ---------
!     Finish up
!     ---------
!
      if (.not. self % child) then
         CALL self % Describe         ( trim(fileName), bFaceOrder )
         CALL self % DescribePartition( )
      end if
!
!     --------------------
!     Prepare mesh for I/O
!     --------------------
!
      call self % PrepareForIO

      deallocate(globalToLocalNodeID)
      deallocate(globalToLocalElementID)

   END SUBROUTINE ConstructMeshPartition_FromGMSHFile_v4_
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   subroutine ConstructSimplestMesh_FromGMSHFile_v4_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
!  ---------------------------------------------------------
!  Build mesh from GMSH file.
!  ---------------------------------------------------------
      USE Physics
      use PartitionedMeshClass
      use MPI_Process_Info
      implicit none
!-----Arguments-----------------------------------------------------------
      type(HexMesh)                   :: self
      integer                         :: nodes
      character(len=*)                :: fileName
      integer                         :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
      integer                         :: dir2D
      logical                         :: periodRelative
      logical           , intent(out) :: success
!-----Local-Variables-----------------------------------------------------
      character(len=1024)             :: tmps
      real(kind=RP)                   :: tmpd
      integer                         :: tmpi, tmpi1, tmpi2, tmpi3, tmp_eltag
      logical                         :: tmpb
      integer, dimension(:), allocatable :: tmpi_vec1, tmpi_vec2
      character(len=127),    allocatable :: tmps_vec(:)
      type(MSH_BCinfo_t), dimension(:), allocatable :: msh_bcs
      integer                         :: msh_no_BCs
      integer                         :: msh_no_points
      integer                         :: msh_no_curves
      integer                         :: msh_no_surfaces
      integer                         :: msh_no_volumes
      type(MSH_point_t), dimension(:), allocatable  :: msh_points
      type(MSH_curve_t), dimension(:), allocatable  :: msh_curves
      type(MSH_surf_t) , dimension(:), allocatable  :: msh_surfaces
      type(MSH_vol_t)  , dimension(:), allocatable  :: msh_volumes
      integer                         :: msh_no_nodeblocks, msh_nodeblock, msh_nodes_per_block, msh_node, msh_global_node
      integer                    :: element_type, org_element_type
      integer                    :: msh_no_elblocks, msh_elblock, msh_els_per_block, msh_el, msh_global_el

      type(MSH_node_block_t)  , dimension(:), allocatable  :: msh_node_blocks
      type(MSH_element_block_t)  , dimension(:), allocatable  :: msh_element_blocks

      character(len=MSH_LEN) :: msh_entity
      real(kind=RP), allocatable :: msh_entity_vec(:)
      integer, dimension(EL_MAX_ORDER)      :: check_eltype

      integer                         :: numberOfElements
      integer                         :: numberOfNodes
      integer                         :: numberOfBoundaryFaces
      integer                         :: numberOfFaces

      integer                         :: bFaceOrder, numBFacePoints
      integer                         :: i, j, k, l
      integer                         :: jj, ii
      integer                         :: fUnit, fileStat
      integer                         :: nodeIDs(NODES_PER_ELEMENT)
      real(kind=RP)                   :: x(NDIM)
      CHARACTER(LEN=BC_STRING_LENGTH), pointer :: zoneNames(:)
      real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
      integer, allocatable            :: tags(:)
!-----Curved-patches------------------------------------------------------
      real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
      real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
!  -----------------------------------------------------------------------

      numberOfBoundaryFaces = 0
      success               = .TRUE.

      fUnit = UnusedUnit()
      OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
      IF ( fileStat /= 0 )     THEN
         PRINT *, "Error opening file: ", fileName
         success = .FALSE.
         RETURN
      END IF

!------------------------------------------------------------------------

!-----Read-header-info---------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$MeshFormat') error stop "READ_GMSH :: Wrong input file."
      read(fUnit,*) tmpd, tmpi, tmpi
      read(fUnit,*) tmps
!------------------------------------------------------------------------

!-----Read-BC-info-------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$PhysicalNames') error stop "READ_GMSH :: Wrong input file - no boundary conditions defined."
      read(fUnit,*) tmpi

      allocate(tmpi_vec1(tmpi))
      allocate(tmpi_vec2(tmpi))
      allocate(tmps_vec(tmpi))

      msh_no_BCs = 0 ! only surfaces count!
      do i=1, tmpi
         read(fUnit,*) tmpi_vec1(i), tmpi_vec2(i), tmps_vec(i)
         if(tmpi_vec1(i) .eq. 2) msh_no_BCs=msh_no_BCs+1 ! check if surface
      end do ! tmpi
      if (msh_no_BCs .eq. 0) print *, "READ_GMSH :: No boundary conditions detected."

      allocate(msh_bcs(msh_no_BCs))

      j = 0
      do i=1, tmpi
         if(tmpi_vec1(i) .eq. 2) then
            j = j + 1
            msh_bcs(j)%dim  = tmpi_vec1(i)
            msh_bcs(j)%tag  = tmpi_vec2(i)
            msh_bcs(j)%name = tmps_vec(i)
         end if
      end do ! msh_no_BCs

      deallocate(tmpi_vec1)
      deallocate(tmpi_vec2)
      deallocate(tmps_vec)

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndPhysicalNames') error stop "READ_GMSH :: Wrong input file - not all boundary conditions detected."
!------------------------------------------------------------------------

!-----Read-msh-entities--------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Entities') error stop "READ_GMSH :: Wrong input file - no entities found."
      read(fUnit,*) msh_no_points, msh_no_curves, msh_no_surfaces, msh_no_volumes

      ! allocate memory for gmsh internal entities
      allocate(msh_points(msh_no_points))
      allocate(msh_curves(msh_no_curves))
      allocate(msh_surfaces(msh_no_surfaces))
      allocate(msh_volumes(msh_no_volumes))
      allocate(msh_entity_vec(32)) ! arbitrary number

!-----Read-points--------------------------------------------------------
      do i=1, msh_no_points
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') msh_entity
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_points(i)%tag      = int(msh_entity_vec(1))
         msh_points(i)%x        = msh_entity_vec(2:4)
         msh_points(i)%no_ptags = int(msh_entity_vec(5))
         if (msh_points(i)%no_ptags .gt. 0) msh_points(i)%ptags(1:msh_points(i)%no_ptags) = int(msh_entity_vec(6:5+msh_points(i)%no_ptags))
      end do ! msh_no_points
!------------------------------------------------------------------------

!-----Read-curves--------------------------------------------------------
      do i=1, msh_no_curves
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') msh_entity
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_curves(i)%tag      = int(msh_entity_vec(1))
         msh_curves(i)%minX     = msh_entity_vec(2:4)
         msh_curves(i)%maxX     = msh_entity_vec(5:7)
         msh_curves(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_curves(i)%no_ptags .gt. 0) msh_curves(i)%ptags(1:msh_curves(i)%no_ptags) = int(msh_entity_vec(9:8+msh_curves(i)%no_ptags))
         msh_curves(i)%no_bps = msh_entity_vec(9+msh_curves(i)%no_ptags)
         msh_curves(i)%bps = 0
         msh_curves(i)%bps(1:msh_curves(i)%no_bps) = msh_entity_vec(10+msh_curves(i)%no_ptags:9+msh_curves(i)%no_ptags+msh_curves(i)%no_bps)
         msh_curves(i)%bps = abs(msh_curves(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_curves
!------------------------------------------------------------------------

!-----Read-surfaces------------------------------------------------------
      do i=1, msh_no_surfaces
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') msh_entity
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_surfaces(i)%tag      = int(msh_entity_vec(1))
         msh_surfaces(i)%minX     = msh_entity_vec(2:4)
         msh_surfaces(i)%maxX     = msh_entity_vec(5:7)
         msh_surfaces(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_surfaces(i)%no_ptags .gt. 0) msh_surfaces(i)%ptags(1:msh_surfaces(i)%no_ptags) = int(msh_entity_vec(9:8+msh_surfaces(i)%no_ptags))
         msh_surfaces(i)%no_bps = msh_entity_vec(9+msh_surfaces(i)%no_ptags)
         msh_surfaces(i)%bps = 0
         msh_surfaces(i)%bps(1:msh_surfaces(i)%no_bps) = msh_entity_vec(10+msh_surfaces(i)%no_ptags:9+msh_surfaces(i)%no_ptags+msh_surfaces(i)%no_bps)
         msh_surfaces(i)%bps = abs(msh_surfaces(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_surfaces
!------------------------------------------------------------------------

!-----Read-volumes-------------------------------------------------------
      do i=1, msh_no_volumes
         msh_entity_vec=0.0_RP
         read(fUnit,'(4096a)') msh_entity
         call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

         msh_volumes(i)%tag      = int(msh_entity_vec(1))
         msh_volumes(i)%minX     = msh_entity_vec(2:4)
         msh_volumes(i)%maxX     = msh_entity_vec(5:7)
         msh_volumes(i)%no_ptags = int(msh_entity_vec(8))
         if (msh_volumes(i)%no_ptags .gt. 0) msh_volumes(i)%ptags(1:msh_volumes(i)%no_ptags) = int(msh_entity_vec(9:8+msh_volumes(i)%no_ptags))
         msh_volumes(i)%no_bps = msh_entity_vec(9+msh_volumes(i)%no_ptags)
         msh_volumes(i)%bps(1:msh_volumes(i)%no_bps) = msh_entity_vec(10+msh_volumes(i)%no_ptags:9+msh_volumes(i)%no_ptags+msh_volumes(i)%no_bps)
         msh_volumes(i)%bps = abs(msh_volumes(i)%bps) ! why is this negative in the .msh no clue
      end do ! msh_no_volumes
!------------------------------------------------------------------------

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndEntities') error stop "READ_GMSH :: Wrong input file - not all entities detected."
!------------------------------------------------------------------------

!-----Read-nodes---------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Nodes') error stop "READ_GMSH :: Wrong input file - no nodes found."

      read(fUnit,*) msh_no_nodeblocks, numberOfNodes, tmpi1, tmpi2
      if (numberOfNodes .ne. tmpi2) error stop "READ_gmsh :: Incoherent node numbering."

      ! allocate nodes storage
      allocate(msh_node_blocks(msh_no_nodeblocks))

      msh_global_node = 0
      do msh_nodeblock=1, msh_no_nodeblocks
         read(fUnit,*) tmpi1, tmpi2, tmpi3, msh_nodes_per_block
         tmpb=tmpi3
         if (tmpb) error stop "READ_gmsh :: Parametric nodes not supported."

         call msh_node_blocks(msh_nodeblock) % Construct(tmpi1, tmpi2, tmpb, msh_nodes_per_block)

         do msh_node=1, msh_nodes_per_block
            read(fUnit,*) msh_node_blocks(msh_nodeblock) % tags(msh_node)
         end do ! msh_nodes_per_block /for tags

         do msh_node=1, msh_nodes_per_block
            read(fUnit,*) msh_node_blocks(msh_nodeblock) % cords(msh_node,:)
         end do ! msh_nodes_per_block /for coordinates

         msh_global_node = msh_global_node + msh_nodes_per_block

      end do ! msh_no_nodeblocks

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndNodes') error stop "READ_GMSH :: Wrong input file - not all nodes detected."
!------------------------------------------------------------------------

!-----Read-elements------------------------------------------------------
      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$Elements') error stop "READ_GMSH :: Wrong input file - no elements found."

      read(fUnit,*) msh_no_elblocks, numberOfElements, tmpi1, tmpi2
      if (numberOfElements .ne. tmpi2) error stop "READ_gmsh :: Incoherent element numbering."

      allocate(msh_element_blocks(msh_no_elblocks))

      msh_global_el = 0
      do msh_elblock=1, msh_no_elblocks

         read(fUnit,*) tmpi1, tmpi2, element_type, msh_els_per_block

         call msh_element_blocks(msh_elblock) % Construct(element_type,msh_els_per_block)

         do msh_el=1, msh_els_per_block
            read(fUnit,*) msh_element_blocks(msh_elblock) % tags(msh_el), msh_element_blocks(msh_elblock) % nodes(msh_el,:)
         end do ! msh_els_per_block /for tags

         msh_global_el = msh_global_el + msh_els_per_block

      end do ! msh_no_elblocks

      read(fUnit,*) tmps
      if (trim(tmps) .ne. '$EndElements') error stop "READ_GMSH :: Wrong input file - not all elements detected."
      close( fUnit )
!------------------------------------------------------------------------

!-----Mesh-info----------------------------------------------------------
      ! find order of elements curvature
      check_eltype = 0
      do i=1, EL_MAX_ORDER
         check_eltype(i) = count(msh_element_blocks(:) % el_type .eq.  SUPPORTED_EL_TYPES(i))
      end do
      if (sum(check_eltype) .eq. 0) error stop "READ_GMSH :: No 3D elements detected in the mesh."
      if (sum(check_eltype) .ne. maxval(check_eltype)) error stop "READ_GMSH :: More than 1 type of hexahedral detected in the mesh."
      bFaceOrder = maxloc(check_eltype,1) ! set order of the mesh
      org_element_type = SUPPORTED_EL_TYPES(bFaceOrder)
      ! find number of elements
      numberOfElements = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            numberOfElements = numberOfElements + msh_element_blocks(msh_elblock) % no_els
         end if
      end do
!------------------------------------------------------------------------

!-----Reorder-nodes-in-elements------------------------------------------
      allocate(tmpi_vec1((bFaceOrder + 1)**3))
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do j=1, msh_element_blocks(msh_elblock) % no_els
               ! re-order nodes within element to match HORSES ordering
               tmpi_vec1 = msh_element_blocks(msh_elblock) % nodes(j,:)
               call ReorderElement(tmpi_vec1,bFaceOrder)
               msh_element_blocks(msh_elblock) % nodes(j,:) = tmpi_vec1
            end do
         end if
      end do
      deallocate(tmpi_vec1)
!------------------------------------------------------------------------

!-----Assign-BCs---------------------------------------------------------
      do i=1,msh_no_BCs
         ! find no surfaces
         msh_bcs(i) % no_of_surfaces = 0
         do j=1, msh_no_surfaces
            tmpi = my_findloc(msh_surfaces(j)%ptags, msh_bcs(i)%tag, 1)
            if (tmpi .gt. 0) msh_bcs(i) % no_of_surfaces = msh_bcs(i) % no_of_surfaces + 1
         end do ! msh_no_surfaces

         ! assign surfaces and curves
         allocate(tmpi_vec1(256))
         tmpi_vec1 = 0
         allocate(msh_bcs(i) % surface_tags(msh_bcs(i) % no_of_surfaces))
         k = 0
         jj = 1
         do j=1, msh_no_surfaces
            tmpi = my_findloc(msh_surfaces(j)%ptags, msh_bcs(i)%tag, 1)
            if (tmpi .gt. 0) then
               k = k + 1
               msh_bcs(i) % surface_tags(k) = msh_surfaces(j) % tag
               tmpi_vec1(jj : jj + count(msh_surfaces(j) % bps /= 0) - 1) = msh_surfaces(j) % bps(1 : count(msh_surfaces(j) % bps /= 0))
               jj = jj + 1 + count(msh_surfaces(j) % bps /= 0)
            end if ! tmpi
         end do ! msh_no_surfaces
         call unique(tmpi_vec1,msh_bcs(i) % curve_tags)
         msh_bcs(i) % no_of_curves = size(msh_bcs(i) % curve_tags,1)
         deallocate(tmpi_vec1)

         ! assign points
         k = 0
         allocate(tmpi_vec1(2 * msh_bcs(i) % no_of_curves))
         tmpi_vec1 = 0
         do j=1, msh_no_curves
            tmpi = my_findloc(msh_bcs(i) % curve_tags, msh_curves(j)%tag, 1)
            if (tmpi .gt. 0) then
               k = k + 1
               tmpi_vec1(1 + (k-1)*2: k*2) = msh_curves(j) % bps
            end if ! tmpi
         end do ! msh_no_points
         call unique(tmpi_vec1,msh_bcs(i) % point_tags)
         msh_bcs(i) % no_of_points = size(msh_bcs(i) % point_tags,1)
         deallocate(tmpi_vec1)

         ! assign nodes
         k = 0
         allocate(tmpi_vec1(numberOfNodes)) ! approximate upper bound
         tmpi_vec1 = 0
         do j=1, msh_no_nodeblocks

            select case (msh_node_blocks(j) % entity_dim)
            case (0) ! points
               tmpi = my_findloc(msh_bcs(i) % point_tags, msh_node_blocks(j) % entity_tag, 1)
            case (1) ! curves
               tmpi = my_findloc(msh_bcs(i) % curve_tags, msh_node_blocks(j) % entity_tag, 1)
            case (2) ! surfaces
               tmpi = my_findloc(msh_bcs(i) % surface_tags, msh_node_blocks(j) % entity_tag, 1)
            case (3) ! volume
               tmpi = 0
            end select ! msh_node_blocks(j) % entity_dim

            if (tmpi .gt. 0) then
               tmpi_vec1(k+1 : k+msh_node_blocks(j) % no_nodes) = msh_node_blocks(j) % tags
               k = k + msh_node_blocks(j) % no_nodes
            end if ! tmpi

         end do ! msh_no_nodeblocks
         call unique(tmpi_vec1(1:k),msh_bcs(i) % node_tags)
         msh_bcs(i) % no_of_nodes = size(msh_bcs(i) % node_tags,1)
         deallocate(tmpi_vec1)

      end do ! msh_no_BCs
!------------------------------------------------------------------------

!-----Assign-BC-to-nodes-------------------------------------------------
      allocate(tmpi_vec1(8)) ! help array for the walls
      do i=1, msh_no_BCs
         do msh_elblock=1, msh_no_elblocks
            if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
               do j=1, msh_element_blocks(msh_elblock) % no_els

                  ! find matching nodes
                  tmpi_vec1=0
                  do l=1,8
                     tmpi = my_findloc(msh_bcs(i) % node_tags,msh_element_blocks(msh_elblock) % nodes(j,l),1)
                     if(tmpi .gt. 0) tmpi_vec1(l) = 1
                     ! TODO: Check why this doesn't work with older ifort
                     ! tmpi = any(msh_bcs(i) % node_tags .eq. msh_element_blocks(msh_elblock) % nodes(j,l))
                     ! tmpi_vec1(l)=tmpi
                  end do

                  ! assign BC face to the elements nodes
                  if ( (sum(tmpi_vec1(1:2)+tmpi_vec1(5:6))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,1) = msh_bcs(i) % tag ! 1:west
                  if ( (sum(tmpi_vec1(3:4)+tmpi_vec1(7:8))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,2) = msh_bcs(i) % tag ! 2:east
                  if ( (sum(tmpi_vec1(1:4))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,3) = msh_bcs(i) % tag ! 3:south
                  if ( (sum(tmpi_vec1(2:3)+tmpi_vec1(6:7))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,4) = msh_bcs(i) % tag ! 4:front
                  if ( (sum(tmpi_vec1(5:8))) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,5) = msh_bcs(i) % tag ! 5:north
                  if ( (tmpi_vec1(1)+tmpi_vec1(4)+tmpi_vec1(5)+tmpi_vec1(8)) .eq. 4) msh_element_blocks(msh_elblock) % BCs(j,6) = msh_bcs(i) % tag ! 6:back

               end do ! msh_element_blocks(msh_elblock) % no_els
            end if
         end do ! msh_no_elblocks
      end do ! msh_no_BCs
      deallocate(tmpi_vec1)
!------------------------------------------------------------------------

!-----Assign-new-tags----------------------------------------------------
      tmp_eltag = 0
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do j=1, msh_element_blocks(msh_elblock) % no_els
               tmp_eltag = tmp_eltag + 1
               msh_element_blocks(msh_elblock) % tags(j) = tmp_eltag
            end do
         end if
      end do
      if (.not. (tmp_eltag .eq. numberOfElements)) error stop "Read_GMSH :: Number of elements inconsistent."
!------------------------------------------------------------------------

!-----Build-nodes--------------------------------------------------------
      self % nodeType = nodes
      self % no_of_elements = numberOfElements
!------------------------------------------------------------------------

!-----Set-up-face-patches------------------------------------------------
      numBFacePoints = bFaceOrder + 1
      allocate(uNodes(numBFacePoints))
      allocate(vNodes(numBFacePoints))
      allocate(values(3,numBFacePoints,numBFacePoints))
      do i = 1, numBFacePoints
         uNodes(i) = -1._RP + (i-1) * (2._RP/bFaceOrder)
         vNodes(i) = uNodes(i)
      end do
!------------------------------------------------------------------------

!-----Allocate-mem-for-elements-and-nodes--------------------------------
      allocate( self % elements(numberOfelements) )
      allocate( self % nodes(numberOfNodes) )
      allocate( self % Nx(numberOfelements) , self % Ny(numberOfelements) , self % Nz(numberOfelements) )
      self % Nx = Nx
      self % Ny = Ny
      self % Nz = Nz
!------------------------------------------------------------------------

!----Set-nodes-----------------------------------------------------------
      do msh_nodeblock=1, msh_no_nodeblocks
         do msh_node=1, msh_node_blocks(msh_nodeblock) % no_nodes
            x = msh_node_blocks(msh_nodeblock) % cords(msh_node,1:NDIM)/Lref
            call ConstructNode( self % nodes(msh_node_blocks(msh_nodeblock) % tags(msh_node)), x, msh_node_blocks(msh_nodeblock) % tags(msh_node) )
         end do ! msh_node_blocks(msh_nodeblock) % no_nodes
      end do ! msh_no_nodeblocks
!------------------------------------------------------------------------

!----Set-elements-----------------------------------------------------------
      j = 0
      call extract_bc_tags(msh_bcs, tags)
      do msh_elblock=1, msh_no_elblocks
         if (msh_element_blocks(msh_elblock) % el_type .eq. org_element_type) then
            do msh_el = 1, msh_element_blocks(msh_elblock) % no_els
               j = j + 1
               ! setting l'th element
               l =  msh_element_blocks(msh_elblock) % tags(msh_el)
               nodeIDs = msh_element_blocks(msh_elblock) % nodes(msh_el,1:8) ! only non-curved nodes

               if (bFaceOrder .eq. 1) then ! non-curved mesh
                  do k = 1, NODES_PER_ELEMENT
                     corners(:,k) = self % nodes(nodeIDs(k)) % x
                  end do
                  self % elements(l) % SurfInfo % IsHex8 = .TRUE.
                  self % elements(l) % SurfInfo % corners = corners
               else ! curved mesh
                  ! allocate tmp arrays for curved face node tags and coordinates
                  allocate(tmpi_vec1(numBFacePoints*numBFacePoints))
                  allocate(tmpi_vec2(numBFacePoints**3))
                  tmpi_vec2 = msh_element_blocks(msh_elblock) % nodes(msh_el,:)

                  do k = 1, FACES_PER_ELEMENT
                     call GetOrderedFaceNodeTags(tmpi_vec1,k,bFaceOrder,tmpi_vec2)
                     do jj = 1, numBFacePoints
                        do ii = 1, numBFacePoints
                           values(:,ii,jj) = self % nodes(tmpi_vec1(ii + (jj-1)*numBFacePoints)) % x
                        end do
                     end do
                     values = values / Lref
                     call self % elements(l) % SurfInfo % facePatches(k) % construct(uNodes, vNodes, values)
                  end do

                  deallocate(tmpi_vec1)
                  deallocate(tmpi_vec2)
               end if

               call self % elements(l) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , l, l)

               ! set element boundaries
               do k = 1, 6
                  tmpi1 = my_findloc( tags,msh_element_blocks(msh_elblock) % BCs(msh_el,k),1)
                  if (tmpi1 .gt. 0) then
                     self % elements(l) % boundaryName(k) = trim(msh_bcs(tmpi1) % name)
                  else
                     self % elements(l) % boundaryName(k) = emptyBCName
                  end if
               end do ! k

               ! set BC names to faces
               do k = 1, 6
                  if(trim(self % elements(l) % boundaryName(k)) /= emptyBCName) then
                     call toLower( self % elements(l) % boundaryName(k) )
                     numberOfBoundaryFaces = numberOfBoundaryFaces + 1
                     zoneNames => zoneNameDictionary % allKeys()
                     if ( all(trim(self % elements(l) % boundaryName(k)) .ne. zoneNames) ) then
                        call zoneNameDictionary % addValueForKey(trim(self % elements(l) % boundaryName(k)), trim(self % elements(l) % boundaryName(k)))
                     end if
                     deallocate (zoneNames)
                  end if
               end do ! k
            end do ! msh_element_blocks(msh_elblock) % no_els
         end if ! if el_type .eq. org_element_type
      end do ! msh_no_elblocks
      if (.not. (j .eq. numberOfElements)) error stop "Read_GMSH :: Not all elements assigned."
!------------------------------------------------------------------------

!-----Deallocate-msh-vairables-------------------------------------------
      do msh_nodeblock=1, msh_no_nodeblocks
         call msh_node_blocks(msh_nodeblock) % Destruct()
      end do
      deallocate(msh_node_blocks)

      do msh_elblock=1, msh_no_elblocks
         call msh_element_blocks(msh_elblock) % Destruct()
      end do
      deallocate(msh_element_blocks)

      do i=1,msh_no_BCs
         call msh_bcs(i) % Destruct()
      end do
      deallocate(msh_bcs)

      deallocate(msh_points)
      deallocate(msh_curves)
      deallocate(msh_surfaces)
      deallocate(msh_volumes)
      deallocate(msh_entity_vec)

      deallocate(uNodes) ! Check if we can do it! FIXME
      deallocate(vNodes) ! Check if we can do it! FIXME
      deallocate(values) ! Check if we can do it! FIXME
!------------------------------------------------------------------------
!
!     ---------------------------
!     Construct the element faces
!     ---------------------------
!
      numberOfFaces        = (6*numberOfElements + numberOfBoundaryFaces)/2
      self % numberOfFaces = numberOfFaces
      allocate( self % faces(self % numberOfFaces) )
      CALL ConstructFaces( self, success )
!
!     -------------------------
!     Build the different zones
!     -------------------------
!
      call self % ConstructZones()
!
!     ---------------------------
!     Construct periodic faces
!     ---------------------------
!
      CALL ConstructPeriodicFaces( self, periodRelative )
!
!     ---------------------------
!     Delete periodic- faces
!     ---------------------------
!
      CALL DeletePeriodicMinusFaces( self )
!
!     ---------------------------
!     Assign faces ID to elements
!     ---------------------------
!
      CALL getElementsFaceIDs(self)
!
!     ---------------------
!     Define boundary faces
!     ---------------------
!
      call self % DefineAsBoundaryFaces()
!
!     ------------------------------
!     Set the element connectivities
!     ------------------------------
!
      call self % SetConnectivitiesAndLinkFaces(nodes)

      call self % ExportBoundaryMesh (trim(fileName))
   end subroutine ConstructSimplestMesh_FromGMSHFile_v4_
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ConstructMesh_FromGMSHFile_v2_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
      !  ---------------------------------------------------------
      !  Build mesh from GMSH file.
      !  ---------------------------------------------------------
            USE Physics
            use PartitionedMeshClass
            use MPI_Process_Info
            implicit none
      !-----Arguments-----------------------------------------------------------
            type(HexMesh)                   :: self
            integer                         :: nodes
            character(len=*)                :: fileName
            integer                         :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
            integer                         :: dir2D
            logical                         :: periodRelative
            logical           , intent(out) :: success
      !-----Local-Variables-----------------------------------------------------
            character(len=1024)             :: tmps
            real(kind=RP)                   :: tmpd
            integer                         :: tmpi, tmpi1, tmpi2, tmpi3, tmp_eltag
            logical                         :: tmpb
            integer, dimension(:), allocatable :: tmpi_vec1, tmpi_vec2, el_types
            character(len=127),    allocatable :: tmps_vec(:)
            type(MSH_BCinfo_t), dimension(:), allocatable :: msh_bcs
            integer                         :: msh_no_BCs
            integer                    :: element_type, org_element_type, org_element_type_2D

            character(len=MSH_LEN) :: msh_entity
            real(kind=RP), allocatable       :: msh_entity_vec(:)
            integer, dimension(EL_MAX_ORDER) :: check_eltype

            type(MSH_node_block_t)     :: msh_nodes
            type(MSH_element_block_t)  :: msh_elements, msh_elements_3D, msh_elements_2D

            integer                         :: numberOfElements
            integer                         :: numberOfElements2D
            integer                         :: numberOfNodes
            integer                         :: numberOfBoundaryFaces
            integer                         :: numberOfFaces
            integer                         :: no_nodes_i, msh_node, msh_el

            integer                         :: bFaceOrder, numBFacePoints, innerEdgePoints
            integer                         :: i, j, k, l, jj, ii
            integer                         :: fUnit, fileStat
            integer                         :: nodeIDs(NODES_PER_ELEMENT)
            real(kind=RP)                   :: x(NDIM)
            CHARACTER(LEN=BC_STRING_LENGTH), pointer :: zoneNames(:)
            real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
            integer, allocatable            :: tags(:)
      !-----Curved-patches------------------------------------------------------
            real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
            real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
      !  -----------------------------------------------------------------------

      !-----Check-if-a-mesh-partition-exists-----------------------------------
            if ( MPI_Process % doMPIAction ) then
               if ( mpi_partition % Constructed ) then
                  call ConstructMeshPartition_FromGMSHFile_v2_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
               else
                  call ConstructSimplestMesh_FromGMSHFile_v2_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
               end if
               return
            end if

            numberOfBoundaryFaces = 0
            success               = .TRUE.

            fUnit = UnusedUnit()
            OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
            IF ( fileStat /= 0 )     THEN
               PRINT *, "Error opening file: ", fileName
               success = .FALSE.
               RETURN
            END IF
      !------------------------------------------------------------------------

      !-----Read-header-info---------------------------------------------------
            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$MeshFormat') error stop "READ_GMSH :: Wrong input file."
            read(fUnit,*) tmpd, tmpi, tmpi
            read(fUnit,*) tmps
      !------------------------------------------------------------------------

      !-----Read-BC-info-------------------------------------------------------
            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$PhysicalNames') error stop "READ_GMSH :: Wrong input file - no boundary conditions defined."
            read(fUnit,*) tmpi

            allocate(tmpi_vec1(tmpi))
            allocate(tmpi_vec2(tmpi))
            allocate(tmps_vec(tmpi))

            msh_no_BCs = 0 ! only surfaces count!
            do i=1, tmpi
               read(fUnit,*) tmpi_vec1(i), tmpi_vec2(i), tmps_vec(i)
               if(tmpi_vec1(i) .eq. 2) msh_no_BCs=msh_no_BCs+1 ! check if surface
            end do ! tmpi
            if (msh_no_BCs .eq. 0) print *, "READ_GMSH :: No boundary conditions detected."

            allocate(msh_bcs(msh_no_BCs))

            j = 0
            do i=1, tmpi
               if(tmpi_vec1(i) .eq. 2) then
                  j = j + 1
                  msh_bcs(j)%dim  = tmpi_vec1(i)
                  msh_bcs(j)%tag  = tmpi_vec2(i)
                  msh_bcs(j)%name = tmps_vec(i)
               end if
            end do ! msh_no_BCs

            deallocate(tmpi_vec1)
            deallocate(tmpi_vec2)
            deallocate(tmps_vec)

            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$EndPhysicalNames') error stop "READ_GMSH :: Wrong input file - not all boundary conditions detected."
      !------------------------------------------------------------------------

      !-----Read-nodes---------------------------------------------------------
            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$Nodes') error stop "READ_GMSH :: Wrong input file - no nodes found."

            read(fUnit,*) numberOfNodes

            call msh_nodes % Construct(0, 0, .false., numberOfNodes)

            do i=1, numberOfNodes
                  read(fUnit,*) msh_nodes % tags(i), msh_nodes % cords(i,:)
            end do ! numberOfNodes

            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$EndNodes') error stop "READ_GMSH :: Wrong input file - not all nodes detected."
      !------------------------------------------------------------------------

      !-----Read-elements------------------------------------------------------
            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$Elements') error stop "READ_GMSH :: Wrong input file - no elements found."

            read(fUnit,*) numberOfElements
            call msh_elements % Construct(element_type,numberOfElements)
            allocate(msh_entity_vec(255)) ! arbitrary number
            allocate(el_types(numberOfElements)) ! arbitrary number
            do i=1, numberOfElements

               read(fUnit,'(4096a)') msh_entity ! read row

               msh_entity_vec=0.0_RP

               call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

               msh_elements % tags(i)     = int(msh_entity_vec(1))
               el_types(i)                = int(msh_entity_vec(2))
               msh_elements % no_ptags(i) = int(msh_entity_vec(3))
               if (msh_elements % no_ptags(i) .gt. 0) msh_elements % ptags(i,1:msh_elements % no_ptags(i)) = &
                  int(msh_entity_vec(4:3+msh_elements % no_ptags(i)))

               select case (el_types(i))
               case (5) ! 3D - 1st order
                  no_nodes_i = 8
               case (12) ! 3D - 2nd order
                  no_nodes_i = 27
               case (92) ! 3D - 3rd order
                  no_nodes_i = 64
               case (93) ! 3D - 4th order
                  no_nodes_i = 125
               case (94) ! 3D - 5th order
                  no_nodes_i = 216
               case (3) ! 2D - 1st order
                  no_nodes_i = 4
               case (10) ! 2D - 2nd order
                  no_nodes_i = 9
               case (36) ! 2D - 3rd order
                  no_nodes_i = 16
               case (37) ! 2D - 4th order
                  no_nodes_i = 25
               case (38) ! 2D - 5th order
                  no_nodes_i = 36
               case default
                  no_nodes_i = 0
               end select

               msh_elements % nodes(i,1:size(msh_entity_vec(4+msh_elements % no_ptags(i):3+msh_elements % no_ptags(i)+no_nodes_i) )) = &
                int(msh_entity_vec(4+msh_elements % no_ptags(i):3+msh_elements % no_ptags(i)+no_nodes_i))

            end do ! numberOfElements

            deallocate(msh_entity_vec)

            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$EndElements') error stop "READ_GMSH :: Wrong input file - not all elements detected."
            close( fUnit )
      !------------------------------------------------------------------------

      !-----Mesh-info----------------------------------------------------------
            ! find order of elements curvature
            check_eltype = 0
            do i=1, EL_MAX_ORDER
               check_eltype(i) = count(el_types .eq. SUPPORTED_EL_TYPES(i))
            end do
            if (sum(check_eltype) .eq. 0) error stop "READ_GMSH :: No 3D elements detected in the mesh."
            if (sum(check_eltype) .ne. maxval(check_eltype)) error stop "READ_GMSH :: More than 1 type of hexahedral detected in the mesh."
            bFaceOrder = maxloc(check_eltype,1) ! set order of the mesh
            org_element_type = SUPPORTED_EL_TYPES(bFaceOrder)
            msh_elements % el_type = org_element_type

            select case (bFaceOrder)
            case (1)
               org_element_type_2D = 3
            case (2)
               org_element_type_2D = 10
            case (3)
               org_element_type_2D = 36
            case (4)
               org_element_type_2D = 37
            case (5)
               org_element_type_2D = 38
            case default
            end select

            if (numberOfElements .ne. (count(el_types .eq. org_element_type) + count(el_types .eq. org_element_type_2D)) ) &
               error stop "READ_GMSH :: Too many different types of elements."
            numberOfElements = count(el_types .eq. org_element_type)
            numberOfElements2D = count(el_types .eq. org_element_type_2D)

            call msh_elements_3D % Construct(org_element_type,numberOfElements)
            call msh_elements_2D % Construct(org_element_type_2D,numberOfElements2D)

            j = 0
            k = 0
            do i=1, size(el_types)

               if ( el_types(i) .eq. org_element_type ) then
                  j = j + 1
                  msh_elements_3D % tags(j) = msh_elements % tags(i)
                  msh_elements_3D % no_ptags(j) = msh_elements % no_ptags(i)
                  msh_elements_3D % ptags(j,:) = msh_elements % ptags(i,1)
                  msh_elements_3D % nodes(j,:) = msh_elements % nodes(i,1:size(msh_elements_3D % nodes(j,:)))
               elseif ( el_types(i) .eq. org_element_type_2D ) then
                  k = k + 1
                  msh_elements_2D % tags(k) = msh_elements % tags(i)
                  msh_elements_2D % no_ptags(k) = msh_elements % no_ptags(i)
                  msh_elements_2D % ptags(k,:) = msh_elements % ptags(i,1)
                  msh_elements_2D % nodes(k,:) = msh_elements % nodes(i,1:size(msh_elements_2D % nodes(k,:)))
               else
                  error stop "READ_GMSH :: Unknown element type in the mesh."
               end if

            end do

            deallocate(el_types)
            call msh_elements % Destruct()
      !------------------------------------------------------------------------

      !-----Reorder-nodes-in-elements------------------------------------------
            allocate(tmpi_vec1((bFaceOrder + 1)**3))
            do j=1, msh_elements_3D % no_els
               ! re-order nodes within element to match HORSES ordering
               tmpi_vec1 = msh_elements_3D % nodes(j,:)
               call ReorderElement(tmpi_vec1,bFaceOrder)
               msh_elements_3D % nodes(j,:) = tmpi_vec1
            end do
            deallocate(tmpi_vec1)
      !------------------------------------------------------------------------

      !-----Assign-BCs---------------------------------------------------------
            do i=1,msh_no_BCs
               ! find # of surfaces, in this version of the reader surfaces are 2D element faces
               msh_bcs(i) % no_of_surfaces = 0
               allocate(tmpi_vec1(numberOfNodes*3)) ! upper bound
               tmpi_vec1 = 0
               k = 0

               do j=1, numberOfElements2D

                  tmpi = my_findloc(msh_elements_2D%ptags(j,:), msh_bcs(i)%tag, 1)
                  if (tmpi .gt. 0) then
                     msh_bcs(i) % no_of_surfaces = msh_bcs(i) % no_of_surfaces + 1
                     tmpi_vec1(k+1 : k+size(msh_elements_2D % nodes(j,:))) = msh_elements_2D % nodes(j,:)
                     k = k + size(msh_elements_2D % nodes(j,:))
                  end if

               end do ! msh_no_surfaces

               call unique(tmpi_vec1(1:k),msh_bcs(i) % node_tags)
               msh_bcs(i) % no_of_nodes = size(msh_bcs(i) % node_tags,1)
               deallocate(tmpi_vec1)

            end do ! msh_no_BCs
      !------------------------------------------------------------------------

      !-----Assign-BC-to-nodes-------------------------------------------------
            allocate(tmpi_vec1(8)) ! help array for the walls
            do i=1, msh_no_BCs
               do j=1, msh_elements_3D % no_els

                  ! find matching nodes
                  tmpi_vec1=0
                  do l=1,8
                     tmpi = my_findloc(msh_bcs(i) % node_tags,msh_elements_3D % nodes(j,l),1)
                     if(tmpi .gt. 0) tmpi_vec1(l) = 1
                     ! TODO: Check why this doesn't work with older ifort
                     ! tmpi = any(msh_bcs(i) % node_tags .eq. msh_element_blocks(msh_elblock) % nodes(j,l))
                     ! tmpi_vec1(l)=tmpi
                  end do

                  ! assign BC face to the elements nodes
                  if ( (sum(tmpi_vec1(1:2)+tmpi_vec1(5:6))) .eq. 4) msh_elements_3D % BCs(j,1) = msh_bcs(i) % tag ! 1:west
                  if ( (sum(tmpi_vec1(3:4)+tmpi_vec1(7:8))) .eq. 4) msh_elements_3D % BCs(j,2) = msh_bcs(i) % tag ! 2:east
                  if ( (sum(tmpi_vec1(1:4))) .eq. 4) msh_elements_3D % BCs(j,3) = msh_bcs(i) % tag ! 3:south
                  if ( (sum(tmpi_vec1(2:3)+tmpi_vec1(6:7))) .eq. 4) msh_elements_3D % BCs(j,4) = msh_bcs(i) % tag ! 4:front
                  if ( (sum(tmpi_vec1(5:8))) .eq. 4) msh_elements_3D % BCs(j,5) = msh_bcs(i) % tag ! 5:north
                  if ( (tmpi_vec1(1)+tmpi_vec1(4)+tmpi_vec1(5)+tmpi_vec1(8)) .eq. 4) msh_elements_3D % BCs(j,6) = msh_bcs(i) % tag ! 6:back

               end do ! msh_elements_3D % no_els
            end do ! msh_no_BCs
            deallocate(tmpi_vec1)
      !------------------------------------------------------------------------

      !-----Assign-new-tags----------------------------------------------------
            tmp_eltag = 0
            do j=1, msh_elements_3D % no_els
               tmp_eltag = tmp_eltag + 1
               msh_elements_3D % tags(j) = tmp_eltag
            end do
            if (.not. (tmp_eltag .eq. numberOfElements)) error stop "Read_GMSH :: Number of elements inconsistent."
      !------------------------------------------------------------------------

      !-----Build-nodes--------------------------------------------------------
            self % nodeType = nodes
            self % no_of_elements = numberOfElements
            self % no_of_allElements = numberOfElements
      !------------------------------------------------------------------------

      !-----Set-up-face-patches------------------------------------------------
            numBFacePoints = bFaceOrder + 1
            innerEdgePoints = bFaceOrder - 1
            allocate(uNodes(numBFacePoints))
            allocate(vNodes(numBFacePoints))
            allocate(values(3,numBFacePoints,numBFacePoints))
            do i = 1, numBFacePoints
               uNodes(i) = -1._RP + (i-1) * (2._RP/bFaceOrder)
               vNodes(i) = uNodes(i)
            end do
      !------------------------------------------------------------------------

      !-----Allocate-mem-for-elements-and-nodes--------------------------------
            allocate( self % elements(numberOfelements) )
            allocate( self % nodes(numberOfNodes) )
            allocate( self % Nx(numberOfelements) , self % Ny(numberOfelements) , self % Nz(numberOfelements) )
            self % Nx = Nx
            self % Ny = Ny
            self % Nz = Nz
      !------------------------------------------------------------------------

      !----Set-nodes-----------------------------------------------------------
            do msh_node=1, msh_nodes % no_nodes
               x = msh_nodes % cords(msh_node,1:NDIM)/Lref
               call ConstructNode( self % nodes(msh_nodes % tags(msh_node)), x, msh_nodes % tags(msh_node) )
            end do ! msh_nodes % no_nodes
      !------------------------------------------------------------------------

      !----Set-elements-----------------------------------------------------------
            j = 0
            call extract_bc_tags(msh_bcs, tags)
            do msh_el = 1, msh_elements_3D % no_els
               j = j + 1
               ! setting l'th element
               l =  msh_elements_3D % tags(msh_el)
               nodeIDs = msh_elements_3D % nodes(msh_el,1:8) ! only non-curved nodes

               if (bFaceOrder .eq. 1) then ! non-curved mesh
                  do k = 1, NODES_PER_ELEMENT
                     corners(:,k) = self % nodes(nodeIDs(k)) % x
                  end do
                  self % elements(l) % SurfInfo % IsHex8 = .TRUE.
                  self % elements(l) % SurfInfo % corners = corners
               else ! curved mesh
                  ! allocate tmp arrays for curved face node tags and coordinates
                  allocate(tmpi_vec1(numBFacePoints*numBFacePoints))
                  allocate(tmpi_vec2(numBFacePoints**3))

                  tmpi_vec2 = msh_elements_3D % nodes(msh_el,:)

                  do k = 1, FACES_PER_ELEMENT
                     call GetOrderedFaceNodeTags(tmpi_vec1,k,bFaceOrder,tmpi_vec2)
                     do jj = 1, numBFacePoints
                        do ii = 1, numBFacePoints
                           values(:,ii,jj) = self % nodes(tmpi_vec1(ii + (jj-1)*numBFacePoints)) % x
                        end do
                     end do
                     values = values / Lref
                     call self % elements(l) % SurfInfo % facePatches(k) % construct(uNodes, vNodes, values)

                  end do

                  deallocate(tmpi_vec1)
                  deallocate(tmpi_vec2)
               end if

               call self % elements(l) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , l, l)

               ! set element boundaries
               do k = 1, 6
                  tmpi1 = my_findloc( tags,msh_elements_3D % BCs(msh_el,k),1)
                  if (tmpi1 .gt. 0) then
                     self % elements(l) % boundaryName(k) = trim(msh_bcs(tmpi1) % name)
                  else
                     self % elements(l) % boundaryName(k) = emptyBCName
                  end if
               end do ! k

               ! set BC names to faces
               do k = 1, 6
                  if(trim(self % elements(l) % boundaryName(k)) /= emptyBCName) then
                     call toLower( self % elements(l) % boundaryName(k) )
                     numberOfBoundaryFaces = numberOfBoundaryFaces + 1
                     zoneNames => zoneNameDictionary % allKeys()
                     if ( all(trim(self % elements(l) % boundaryName(k)) .ne. zoneNames) ) then
                        call zoneNameDictionary % addValueForKey(trim(self % elements(l) % boundaryName(k)), trim(self % elements(l) % boundaryName(k)))
                     end if
                     deallocate (zoneNames)
                  end if
               end do ! k
            end do ! msh_elements_3D % no_els
            if (.not. (j .eq. numberOfElements)) error stop "Read_GMSH :: Not all elements assigned."
      !------------------------------------------------------------------------

      !-----Deallocate-msh-vairables-------------------------------------------
            do i=1,msh_no_BCs
               call msh_bcs(i) % Destruct()
            end do
            deallocate(msh_bcs)

            call msh_nodes % Destruct()
            call msh_elements_3D % Destruct()
            call msh_elements_2D % Destruct()

            deallocate(uNodes) ! Check if we can do it! FIXME
            deallocate(vNodes) ! Check if we can do it! FIXME
            deallocate(values) ! Check if we can do it! FIXME
      !------------------------------------------------------------------------
      !
      !     ---------------------------
      !     Construct the element faces
      !     ---------------------------
      !
            numberOfFaces        = (6*numberOfElements + numberOfBoundaryFaces)/2
            self % numberOfFaces = numberOfFaces
            allocate( self % faces(self % numberOfFaces) )
            CALL ConstructFaces( self, success )
      !
      !     -------------------------
      !     Build the different zones
      !     -------------------------
      !
            call self % ConstructZones()
      !
      !     ---------------------------
      !     Construct periodic faces
      !     ---------------------------
      !
            CALL ConstructPeriodicFaces( self, periodRelative )
      !
      !     ---------------------------
      !     Delete periodic- faces
      !     ---------------------------
      !
            CALL DeletePeriodicMinusFaces( self )
      !
      !     ---------------------------
      !     Assign faces ID to elements
      !     ---------------------------
      !
            CALL getElementsFaceIDs(self)
      !
      !     ---------------------
      !     Define boundary faces
      !     ---------------------
      !
            call self % DefineAsBoundaryFaces()
      !
      !     -----------------------------------
      !     Check if this is a 2D extruded mesh
      !     -----------------------------------
      !
            call self % CheckIfMeshIs2D()
      !
      !     -------------------------------
      !     Set the mesh as 2D if requested
      !     -------------------------------
      !
            self % dir2D_ctrl = dir2D
            if ( dir2D .ne. 0 ) then
               call SetMappingsToCrossProduct
               call self % CorrectOrderFor2DMesh(dir2D,1)
            end if
      !
      !     ------------------------------
      !     Set the element connectivities
      !     ------------------------------
      !
            call self % SetConnectivitiesAndLinkFaces(nodes)
      !
      !     ---------------------------------------
      !     Construct elements' and faces' geometry
      !     ---------------------------------------
      !
            call self % ConstructGeometry()
      !
      !     -------------------------------
      !     Set the mesh as 2D if requested
      !     -------------------------------
      !
            if ( dir2D .ne. 0 ) then
               call self % CorrectOrderFor2DMesh(dir2D,0)
            end if
      !
      !     ---------------------------------
      !     Describe mesh and prepare for I/O
      !     ---------------------------------
      !
            if (.not. self % child) CALL self % Describe( trim(fileName), bFaceOrder )
            call self % PrepareForIO

            call self % ExportBoundaryMesh (trim(fileName))

         end subroutine ConstructMesh_FromGMSHFile_v2_
      !
      !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      !
      !     ------------------------------
      !     Constructor of mesh partitions
      !     ------------------------------
         SUBROUTINE ConstructMeshPartition_FromGMSHFile_v2_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
      !  ---------------------------------------------------------
      !  Build mesh from GMSH file.
      !  ---------------------------------------------------------
            USE Physics
            use PartitionedMeshClass
            use MPI_Process_Info
            use MPI_Face_Class
            implicit none
      !-----Arguments-----------------------------------------------------------
            type(HexMesh)                   :: self
            integer                         :: nodes
            character(len=*)                :: fileName
            integer                         :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
            integer                         :: dir2D
            logical                         :: periodRelative
            logical           , intent(out) :: success
      !-----Local-Variables-----------------------------------------------------
            character(len=1024)             :: tmps
            real(kind=RP)                   :: tmpd
            integer                         :: tmpi, tmpi1, tmpi2, tmpi3, tmp_eltag
            logical                         :: tmpb
            integer, dimension(:), allocatable :: tmpi_vec1, tmpi_vec2, el_types
            character(len=127),    allocatable :: tmps_vec(:)
            type(MSH_BCinfo_t), dimension(:), allocatable :: msh_bcs
            integer                         :: msh_no_BCs
            integer                    :: element_type, org_element_type, org_element_type_2D

            character(len=MSH_LEN) :: msh_entity
            real(kind=RP), allocatable :: msh_entity_vec(:)
            integer, dimension(EL_MAX_ORDER) :: check_eltype

            type(MSH_node_block_t)     :: msh_nodes
            type(MSH_element_block_t)  :: msh_elements, msh_elements_3D, msh_elements_2D

            type(Node)  , dimension(:), allocatable  :: local_nodes

            integer                         :: numberOfAllElements
            integer                         :: numberOfAllNodes
            integer, allocatable            :: globalToLocalNodeID(:)
            integer, allocatable            :: globalToLocalElementID(:)

            integer                         :: numberOfElements
            integer                         :: numberOfElements2D
            integer                         :: numberOfNodes
            integer                         :: numberOfBoundaryFaces
            integer                         :: numberOfFaces
            integer                         :: no_nodes_i, msh_node, msh_el

            integer                         :: bFaceOrder, numBFacePoints, innerEdgePoints
            integer                         :: i, j, k, l, jj, ii, pNode, pElement
            integer                         :: fUnit, fileStat
            integer                         :: nodeIDs(NODES_PER_ELEMENT)
            real(kind=RP)                   :: x(NDIM)
            CHARACTER(LEN=BC_STRING_LENGTH) :: names(FACES_PER_ELEMENT)
            CHARACTER(LEN=BC_STRING_LENGTH), pointer :: zoneNames(:)
            real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
            integer, allocatable            :: tags(:)
      !-----Curved-patches------------------------------------------------------
            real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
            real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
      !  -----------------------------------------------------------------------

            success               = .TRUE.

            fUnit = UnusedUnit()
            OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
            IF ( fileStat /= 0 )     THEN
               PRINT *, "Error opening file: ", fileName
               success = .FALSE.
               RETURN
            END IF
      !------------------------------------------------------------------------

      !-----Read-header-info---------------------------------------------------
            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$MeshFormat') error stop "READ_GMSH :: Wrong input file."
            read(fUnit,*) tmpd, tmpi, tmpi
            read(fUnit,*) tmps
      !------------------------------------------------------------------------

      !-----Read-BC-info-------------------------------------------------------
            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$PhysicalNames') error stop "READ_GMSH :: Wrong input file - no boundary conditions defined."
            read(fUnit,*) tmpi

            allocate(tmpi_vec1(tmpi))
            allocate(tmpi_vec2(tmpi))
            allocate(tmps_vec(tmpi))

            msh_no_BCs = 0 ! only surfaces count!
            do i=1, tmpi
               read(fUnit,*) tmpi_vec1(i), tmpi_vec2(i), tmps_vec(i)
               if(tmpi_vec1(i) .eq. 2) msh_no_BCs=msh_no_BCs+1 ! check if surface
            end do ! tmpi
            if (msh_no_BCs .eq. 0) print *, "READ_GMSH :: No boundary conditions detected."

            allocate(msh_bcs(msh_no_BCs))

            j = 0
            do i=1, tmpi
               if(tmpi_vec1(i) .eq. 2) then
                  j = j + 1
                  msh_bcs(j)%dim  = tmpi_vec1(i)
                  msh_bcs(j)%tag  = tmpi_vec2(i)
                  msh_bcs(j)%name = tmps_vec(i)
               end if
            end do ! msh_no_BCs

            deallocate(tmpi_vec1)
            deallocate(tmpi_vec2)
            deallocate(tmps_vec)

            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$EndPhysicalNames') error stop "READ_GMSH :: Wrong input file - not all boundary conditions detected."
      !------------------------------------------------------------------------

      !-----Read-nodes---------------------------------------------------------
            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$Nodes') error stop "READ_GMSH :: Wrong input file - no nodes found."

            read(fUnit,*) numberOfNodes

            call msh_nodes % Construct(0, 0, .false., numberOfNodes)

            do i=1, numberOfNodes
                  read(fUnit,*) msh_nodes % tags(i), msh_nodes % cords(i,:)
            end do ! numberOfNodes

            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$EndNodes') error stop "READ_GMSH :: Wrong input file - not all nodes detected."
            numberOfAllNodes = numberOfNodes
      !------------------------------------------------------------------------

      !-----Read-elements------------------------------------------------------
            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$Elements') error stop "READ_GMSH :: Wrong input file - no elements found."

            read(fUnit,*) numberOfElements
            call msh_elements % Construct(element_type,numberOfElements)
            allocate(msh_entity_vec(255)) ! arbitrary number
            allocate(el_types(numberOfElements)) ! arbitrary number
            do i=1, numberOfElements

               read(fUnit,'(4096a)') msh_entity ! read row

               msh_entity_vec=0.0_RP

               call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

               msh_elements % tags(i)     = int(msh_entity_vec(1))
               el_types(i)                = int(msh_entity_vec(2))
               msh_elements % no_ptags(i) = int(msh_entity_vec(3))
               if (msh_elements % no_ptags(i) .gt. 0) msh_elements % ptags(i,1:msh_elements % no_ptags(i)) = &
                  int(msh_entity_vec(4:3+msh_elements % no_ptags(i)))

                  select case (el_types(i))
                  case (5) ! 3D - 1st order
                     no_nodes_i = 8
                  case (12) ! 3D - 2nd order
                     no_nodes_i = 27
                  case (92) ! 3D - 3rd order
                     no_nodes_i = 64
                  case (93) ! 3D - 4th order
                     no_nodes_i = 125
                  case (94) ! 3D - 5th order
                     no_nodes_i = 216
                  case (3) ! 2D - 1st order
                     no_nodes_i = 4
                  case (10) ! 2D - 2nd order
                     no_nodes_i = 9
                  case (36) ! 2D - 3rd order
                     no_nodes_i = 16
                  case (37) ! 2D - 4th order
                     no_nodes_i = 25
                  case (38) ! 2D - 5th order
                     no_nodes_i = 36
                  case default
                     no_nodes_i = 0
                  end select

               msh_elements % nodes(i,1:size(msh_entity_vec(4+msh_elements % no_ptags(i):3+msh_elements % no_ptags(i)+no_nodes_i) )) = &
                int(msh_entity_vec(4+msh_elements % no_ptags(i):3+msh_elements % no_ptags(i)+no_nodes_i))

            end do ! numberOfElements

            deallocate(msh_entity_vec)

            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$EndElements') error stop "READ_GMSH :: Wrong input file - not all elements detected."
            close( fUnit )
      !------------------------------------------------------------------------

      !-----Mesh-info----------------------------------------------------------
            ! find order of elements curvature
            check_eltype = 0
            do i=1, EL_MAX_ORDER
               check_eltype(i) = count(el_types .eq. SUPPORTED_EL_TYPES(i))
            end do
            if (sum(check_eltype) .eq. 0) error stop "READ_GMSH :: No 3D elements detected in the mesh."
            if (sum(check_eltype) .ne. maxval(check_eltype)) error stop "READ_GMSH :: More than 1 type of hexahedral detected in the mesh."
            bFaceOrder = maxloc(check_eltype,1) ! set order of the mesh
            org_element_type = SUPPORTED_EL_TYPES(bFaceOrder)
            msh_elements % el_type = org_element_type

            select case (bFaceOrder)
            case (1)
               org_element_type_2D = 3
            case (2)
               org_element_type_2D = 10
            case (3)
               org_element_type_2D = 36
            case (4)
               org_element_type_2D = 37
            case (5)
               org_element_type_2D = 38
            case default
            end select

            if (numberOfElements .ne. (count(el_types .eq. org_element_type) + count(el_types .eq. org_element_type_2D)) ) &
               error stop "READ_GMSH :: Too many different types of elements."
            numberOfElements = count(el_types .eq. org_element_type)
            numberOfElements2D = count(el_types .eq. org_element_type_2D)

            call msh_elements_3D % Construct(org_element_type,numberOfElements)
            call msh_elements_2D % Construct(org_element_type_2D,numberOfElements2D)

            j = 0
            k = 0
            do i=1, size(el_types)

               if ( el_types(i) .eq. org_element_type ) then
                  j = j + 1
                  msh_elements_3D % tags(j) = msh_elements % tags(i)
                  msh_elements_3D % no_ptags(j) = msh_elements % no_ptags(i)
                  msh_elements_3D % ptags(j,:) = msh_elements % ptags(i,1)
                  msh_elements_3D % nodes(j,:) = msh_elements % nodes(i,1:size(msh_elements_3D % nodes(j,:)))
               elseif ( el_types(i) .eq. org_element_type_2D ) then
                  k = k + 1
                  msh_elements_2D % tags(k) = msh_elements % tags(i)
                  msh_elements_2D % no_ptags(k) = msh_elements % no_ptags(i)
                  msh_elements_2D % ptags(k,:) = msh_elements % ptags(i,1)
                  msh_elements_2D % nodes(k,:) = msh_elements % nodes(i,1:size(msh_elements_2D % nodes(k,:)))
               else
                  error stop "READ_GMSH :: Unknown element type in the mesh."
               end if

            end do

            deallocate(el_types)
            call msh_elements % Destruct()
            numberOfAllElements = numberOfElements
      !------------------------------------------------------------------------

      !-----Reorder-nodes-in-elements------------------------------------------
            allocate(tmpi_vec1((bFaceOrder + 1)**3))
            do j=1, msh_elements_3D % no_els
               ! re-order nodes within element to match HORSES ordering
               tmpi_vec1 = msh_elements_3D % nodes(j,:)
               call ReorderElement(tmpi_vec1,bFaceOrder)
               msh_elements_3D % nodes(j,:) = tmpi_vec1
            end do
            deallocate(tmpi_vec1)
      !------------------------------------------------------------------------

      !-----Assign-BCs---------------------------------------------------------
            do i=1,msh_no_BCs
               ! find # of surfaces, in this version of the reader surfaces are 2D element faces
               msh_bcs(i) % no_of_surfaces = 0
               allocate(tmpi_vec1(numberOfNodes*3)) ! upper bound
               tmpi_vec1 = 0
               k = 0

               do j=1, numberOfElements2D

                  tmpi = my_findloc(msh_elements_2D%ptags(j,:), msh_bcs(i)%tag, 1)
                  if (tmpi .gt. 0) then
                     msh_bcs(i) % no_of_surfaces = msh_bcs(i) % no_of_surfaces + 1
                     tmpi_vec1(k+1 : k+size(msh_elements_2D % nodes(j,:))) = msh_elements_2D % nodes(j,:)
                     k = k + size(msh_elements_2D % nodes(j,:))
                  end if

               end do ! msh_no_surfaces

               call unique(tmpi_vec1(1:k),msh_bcs(i) % node_tags)
               msh_bcs(i) % no_of_nodes = size(msh_bcs(i) % node_tags,1)
               deallocate(tmpi_vec1)

            end do ! msh_no_BCs
      !------------------------------------------------------------------------

      !-----Assign-BC-to-nodes-------------------------------------------------
            allocate(tmpi_vec1(8)) ! help array for the walls
            do i=1, msh_no_BCs
               do j=1, msh_elements_3D % no_els

                  ! find matching nodes
                  tmpi_vec1=0
                  do l=1,8
                     tmpi = my_findloc(msh_bcs(i) % node_tags,msh_elements_3D % nodes(j,l),1)
                     if(tmpi .gt. 0) tmpi_vec1(l) = 1
                     ! TODO: Check why this doesn't work with older ifort
                     ! tmpi = any(msh_bcs(i) % node_tags .eq. msh_element_blocks(msh_elblock) % nodes(j,l))
                     ! tmpi_vec1(l)=tmpi
                  end do

                  ! assign BC face to the elements nodes
                  if ( (sum(tmpi_vec1(1:2)+tmpi_vec1(5:6))) .eq. 4) msh_elements_3D % BCs(j,1) = msh_bcs(i) % tag ! 1:west
                  if ( (sum(tmpi_vec1(3:4)+tmpi_vec1(7:8))) .eq. 4) msh_elements_3D % BCs(j,2) = msh_bcs(i) % tag ! 2:east
                  if ( (sum(tmpi_vec1(1:4))) .eq. 4) msh_elements_3D % BCs(j,3) = msh_bcs(i) % tag ! 3:south
                  if ( (sum(tmpi_vec1(2:3)+tmpi_vec1(6:7))) .eq. 4) msh_elements_3D % BCs(j,4) = msh_bcs(i) % tag ! 4:front
                  if ( (sum(tmpi_vec1(5:8))) .eq. 4) msh_elements_3D % BCs(j,5) = msh_bcs(i) % tag ! 5:north
                  if ( (tmpi_vec1(1)+tmpi_vec1(4)+tmpi_vec1(5)+tmpi_vec1(8)) .eq. 4) msh_elements_3D % BCs(j,6) = msh_bcs(i) % tag ! 6:back

               end do ! msh_elements_3D % no_els
            end do ! msh_no_BCs
            deallocate(tmpi_vec1)
      !------------------------------------------------------------------------

      !-----Assign-new-tags----------------------------------------------------
            tmp_eltag = 0
            do j=1, msh_elements_3D % no_els
               tmp_eltag = tmp_eltag + 1
               msh_elements_3D % tags(j) = tmp_eltag
            end do
            if (.not. (tmp_eltag .eq. numberOfElements)) error stop "Read_GMSH :: Number of elements inconsistent."
      !------------------------------------------------------------------------

      !-----Build-nodes--------------------------------------------------------
            self % nodeType = nodes
            self % no_of_elements = mpi_partition % no_of_elements
            self % no_of_allElements = numberOfAllElements
      !------------------------------------------------------------------------

      !-----Set-up-face-patches------------------------------------------------
            numBFacePoints = bFaceOrder + 1
            allocate(uNodes(numBFacePoints))
            allocate(vNodes(numBFacePoints))
            allocate(values(3,numBFacePoints,numBFacePoints))
            do i = 1, numBFacePoints
               uNodes(i) = -1._RP + (i-1) * (2._RP/bFaceOrder)
               vNodes(i) = uNodes(i)
            end do
      !------------------------------------------------------------------------

      !-----Allocate-mem-for-elements-and-nodes--------------------------------
            allocate( self % elements(mpi_partition % no_of_elements) )
            allocate( self % nodes(mpi_partition % no_of_nodes) )
            allocate ( self % Nx(self % no_of_elements) , self % Ny(self % no_of_elements) , self % Nz(self % no_of_elements) )
            allocate( globalToLocalNodeID(numberOfAllNodes) )
            allocate( globalToLocalElementID(numberOfAllElements) )

            globalToLocalNodeID = -1
            globalToLocalElementID = -1
      !------------------------------------------------------------------------

      !----Set-nodes-----------------------------------------------------------
            pNode = 1
            do msh_node=1, msh_nodes % no_nodes
               x = msh_nodes % cords(msh_node,1:NDIM)/Lref

               if ( pNode .gt. mpi_partition % no_of_nodes ) cycle

               ! Construct only nodes that belong to the partition
               if ( msh_nodes % tags(msh_node) .eq. mpi_partition % nodeIDs(pNode) ) then
                  call ConstructNode( self % nodes(pNode), x, msh_nodes % tags(msh_node) )
                  globalToLocalNodeID(msh_nodes % tags(msh_node)) = pNode
                  pNode = pNode + 1
               end if
            end do ! msh_nodes % no_nodes
      !------------------------------------------------------------------------

      !----Set-local-nodes-----------------------------------------------------
            allocate(local_nodes(numberOfAllNodes))
            do msh_node=1, msh_nodes % no_nodes
               x = msh_nodes % cords(msh_node,1:NDIM)/Lref
               call ConstructNode( local_nodes(msh_nodes % tags(msh_node)), x, msh_nodes % tags(msh_node) )
            end do ! msh_nodes % no_nodes
      !------------------------------------------------------------------------

      !----Set-elements-----------------------------------------------------------
            j = 0
            pElement = 1
            call extract_bc_tags(msh_bcs, tags)
            do msh_el = 1, msh_elements_3D % no_els
               j = j + 1
               ! setting l'th element
               l =  msh_elements_3D % tags(msh_el)
               nodeIDs = msh_elements_3D % nodes(msh_el,1:8) ! only non-curved nodes
               nodeIDs = globalToLocalNodeID(nodeIDs)

               if ( pElement .gt. mpi_partition % no_of_elements ) then

                  ! set element boundaries
                  do k = 1, 6
                     tmpi1 = my_findloc( tags,msh_elements_3D % BCs(msh_el,k),1)
                     if (tmpi1 .gt. 0) then
                        names(k) = trim(msh_bcs(tmpi1) % name)
                     else
                        names(k) = emptyBCName
                     end if
                  end do ! k

                  ! set BC names to faces
                  do k = 1, 6
                     if(trim(names(k)) /= emptyBCName) then
                        call toLower( names(k) )
                        zoneNames => zoneNameDictionary % allKeys()
                        if ( all(trim(names(k)) .ne. zoneNames) ) then
                           call zoneNameDictionary % addValueForKey(trim(names(k)), trim(names(k)))
                        end if
                        deallocate (zoneNames)
                     end if
                  end do ! k

                  cycle
               else if ( l .ne. mpi_partition % elementIDs(pElement) ) then

                  ! set element boundaries
                  do k = 1, 6
                     tmpi1 = my_findloc( tags,msh_elements_3D % BCs(msh_el,k),1)
                     if (tmpi1 .gt. 0) then
                        names(k) = trim(msh_bcs(tmpi1) % name)
                     else
                        names(k) = emptyBCName
                     end if
                  end do ! k

                  ! set BC names to faces
                  do k = 1, 6
                     if(trim(names(k)) /= emptyBCName) then
                        call toLower( names(k) )
                        zoneNames => zoneNameDictionary % allKeys()
                        if ( all(trim(names(k)) .ne. zoneNames) ) then
                           call zoneNameDictionary % addValueForKey(trim(names(k)), trim(names(k)))
                        end if
                        deallocate (zoneNames)
                     end if
                  end do ! k

                  cycle
               end if

               if (bFaceOrder .eq. 1) then ! non-curved mesh
                  do k = 1, NODES_PER_ELEMENT
                     corners(:,k) = self % nodes(nodeIDs(k)) % x
                  end do
                  self % elements(pElement) % SurfInfo % IsHex8 = .TRUE.
                  self % elements(pElement) % SurfInfo % corners = corners
               else ! curved mesh
                  ! allocate tmp arrays for curved face node tags and coordinates
                  allocate(tmpi_vec1(numBFacePoints*numBFacePoints))
                  allocate(tmpi_vec2(numBFacePoints**3))
                  tmpi_vec2 = msh_elements_3D % nodes(msh_el,:)

                  do k = 1, FACES_PER_ELEMENT
                     call GetOrderedFaceNodeTags(tmpi_vec1,k,bFaceOrder,tmpi_vec2)
                     do jj = 1, numBFacePoints
                        do ii = 1, numBFacePoints
                           values(:,ii,jj) = local_nodes(tmpi_vec1(ii + (jj-1)*numBFacePoints)) % x
                        end do
                     end do
                     values = values / Lref
                     call self % elements(pElement) % SurfInfo % facePatches(k) % construct(uNodes, vNodes, values)
                  end do

                  deallocate(tmpi_vec1)
                  deallocate(tmpi_vec2)
               end if

               call self % elements(pElement) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , pElement, l)

               self % Nx(pElement) = Nx(l)
               self % Ny(pElement) = Ny(l)
               self % Nz(pElement) = Nz(l)

               ! set element boundaries
               do k = 1, 6
                  tmpi1 = my_findloc( tags,msh_elements_3D % BCs(msh_el,k),1)
                  if (tmpi1 .gt. 0) then
                     self % elements(pElement) % boundaryName(k) = trim(msh_bcs(tmpi1) % name)
                  else
                     self % elements(pElement) % boundaryName(k) = emptyBCName
                  end if
               end do ! k

               ! set BC names to faces
               do k = 1, 6
                  if(trim(self % elements(pElement) % boundaryName(k)) /= emptyBCName) then
                     call toLower( self % elements(pElement) % boundaryName(k) )
                     zoneNames => zoneNameDictionary % allKeys()
                     if ( all(trim(self % elements(pElement) % boundaryName(k)) .ne. zoneNames) ) then
                        call zoneNameDictionary % addValueForKey(trim(self % elements(pElement) % boundaryName(k)), trim(self % elements(pElement) % boundaryName(k)))
                     end if
                     deallocate (zoneNames)
                  end if
               end do ! k

               globalToLocalElementID(l) = pElement
               pElement = pElement + 1
            end do ! msh_elements_3D % no_els
      !------------------------------------------------------------------------

      !-----Deallocate-msh-vairables-------------------------------------------
            do i=1,msh_no_BCs
               call msh_bcs(i) % Destruct()
            end do
            deallocate(msh_bcs)

            call msh_nodes % Destruct()
            call msh_elements_3D % Destruct()
            call msh_elements_2D % Destruct()

            deallocate(uNodes) ! Check if we can do it! FIXME
            deallocate(vNodes) ! Check if we can do it! FIXME
            deallocate(values) ! Check if we can do it! FIXME

            ! deallocate local nodes
            deallocate(local_nodes)
      !------------------------------------------------------------------------
      !
      !     ---------------------------
      !     Construct the element faces
      !     ---------------------------
      !
            numberOfFaces        = GetOriginalNumberOfFaces(self)
            self % numberOfFaces = numberOfFaces
            allocate( self % faces(self % numberOfFaces) )
            CALL ConstructFaces( self, success )
      !        --------------------------------
      !        Get actual mesh element face IDs
      !        --------------------------------
      !
            CALL getElementsFaceIDs(self)
      !
      !     --------------
      !     Cast MPI faces
      !     --------------
      !
            call ConstructMPIFaces( self % MPIfaces )
            call self % UpdateFacesWithPartition(mpi_partition, &
                                                 numberOfAllElements, &
                                                 globalToLocalElementID)
      !
      !     -------------------------
      !     Build the different zones
      !     -------------------------
      !
            call self % ConstructZones()
      !
      !     ---------------------------
      !     Construct periodic faces
      !     ---------------------------
      !
            CALL ConstructPeriodicFaces( self, periodRelative )
      !
      !     ---------------------------
      !     Delete periodic- faces
      !     ---------------------------
      !
            CALL DeletePeriodicMinusFaces( self )
      !
      !     ---------------------------
      !     Assign faces ID to elements
      !     ---------------------------
      !
            CALL getElementsFaceIDs(self)
      !
      !     ---------------------
      !     Define boundary faces
      !     ---------------------
      !
            call self % DefineAsBoundaryFaces()
      !
      !     -----------------------------------
      !     Check if this is a 2D extruded mesh
      !     -----------------------------------
      !
            call self % CheckIfMeshIs2D()
      !
      !     -------------------------------
      !     Set the mesh as 2D if requested
      !     -------------------------------
      !
            if ( dir2D .ne. 0 ) then
               call SetMappingsToCrossProduct
               call self % CorrectOrderFor2DMesh(dir2D,1)
            end if
      !
      !     ------------------------------
      !     Set the element connectivities
      !     ------------------------------
      !
            call self % SetConnectivitiesAndLinkFaces(nodes)
      !
      !     ---------------------------------------
      !     Construct elements' and faces' geometry
      !     ---------------------------------------
      !
            call self % ConstructGeometry()

            if ( dir2D .ne. 0 ) then
               call self % CorrectOrderFor2DMesh(dir2D,0)
            end if
      !
      !     ---------
      !     Finish up
      !     ---------
      !
            if (.not. self % child) then
               CALL self % Describe         ( trim(fileName), bFaceOrder )
               CALL self % DescribePartition( )
            end if
      !
      !     --------------------
      !     Prepare mesh for I/O
      !     --------------------
      !
            call self % PrepareForIO

            deallocate(globalToLocalNodeID)
            deallocate(globalToLocalElementID)

         END SUBROUTINE ConstructMeshPartition_FromGMSHFile_v2_
      !
      !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         subroutine ConstructSimplestMesh_FromGMSHFile_v2_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
      !  ---------------------------------------------------------
      !  Build mesh from GMSH file.
      !  ---------------------------------------------------------
            USE Physics
            use PartitionedMeshClass
            use MPI_Process_Info
            implicit none
      !-----Arguments-----------------------------------------------------------
            type(HexMesh)                   :: self
            integer                         :: nodes
            character(len=*)                :: fileName
            integer                         :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
            integer                         :: dir2D
            logical                         :: periodRelative
            logical           , intent(out) :: success
      !-----Local-Variables-----------------------------------------------------
            character(len=1024)             :: tmps
            real(kind=RP)                   :: tmpd
            integer                         :: tmpi, tmpi1, tmpi2, tmpi3, tmp_eltag
            logical                         :: tmpb
            integer, dimension(:), allocatable :: tmpi_vec1, tmpi_vec2, el_types
            character(len=127),    allocatable :: tmps_vec(:)
            type(MSH_BCinfo_t), dimension(:), allocatable :: msh_bcs
            integer                         :: msh_no_BCs
            integer                    :: element_type, org_element_type, org_element_type_2D

            character(len=MSH_LEN) :: msh_entity
            real(kind=RP), allocatable :: msh_entity_vec(:)
            integer, dimension(EL_MAX_ORDER) :: check_eltype

            type(MSH_node_block_t)     :: msh_nodes
            type(MSH_element_block_t)  :: msh_elements, msh_elements_3D, msh_elements_2D

            integer                         :: numberOfElements
            integer                         :: numberOfElements2D
            integer                         :: numberOfNodes
            integer                         :: numberOfBoundaryFaces
            integer                         :: numberOfFaces
            integer                         :: no_nodes_i, msh_node, msh_el

            integer                         :: bFaceOrder, numBFacePoints, innerEdgePoints
            integer                         :: i, j, k, l, jj, ii
            integer                         :: fUnit, fileStat
            integer                         :: nodeIDs(NODES_PER_ELEMENT)
            real(kind=RP)                   :: x(NDIM)
            CHARACTER(LEN=BC_STRING_LENGTH), pointer :: zoneNames(:)
            real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
            integer, allocatable            :: tags(:)
      !-----Curved-patches------------------------------------------------------
            real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
            real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
      !  -----------------------------------------------------------------------

            numberOfBoundaryFaces = 0
            success               = .TRUE.

            fUnit = UnusedUnit()
            OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
            IF ( fileStat /= 0 )     THEN
               PRINT *, "Error opening file: ", fileName
               success = .FALSE.
               RETURN
            END IF

      !------------------------------------------------------------------------

      !-----Read-header-info---------------------------------------------------
            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$MeshFormat') error stop "READ_GMSH :: Wrong input file."
            read(fUnit,*) tmpd, tmpi, tmpi
            read(fUnit,*) tmps
      !------------------------------------------------------------------------

      !-----Read-BC-info-------------------------------------------------------
            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$PhysicalNames') error stop "READ_GMSH :: Wrong input file - no boundary conditions defined."
            read(fUnit,*) tmpi

            allocate(tmpi_vec1(tmpi))
            allocate(tmpi_vec2(tmpi))
            allocate(tmps_vec(tmpi))

            msh_no_BCs = 0 ! only surfaces count!
            do i=1, tmpi
               read(fUnit,*) tmpi_vec1(i), tmpi_vec2(i), tmps_vec(i)
               if(tmpi_vec1(i) .eq. 2) msh_no_BCs=msh_no_BCs+1 ! check if surface
            end do ! tmpi
            if (msh_no_BCs .eq. 0) print *, "READ_GMSH :: No boundary conditions detected."

            allocate(msh_bcs(msh_no_BCs))

            j = 0
            do i=1, tmpi
               if(tmpi_vec1(i) .eq. 2) then
                  j = j + 1
                  msh_bcs(j)%dim  = tmpi_vec1(i)
                  msh_bcs(j)%tag  = tmpi_vec2(i)
                  msh_bcs(j)%name = tmps_vec(i)
               end if
            end do ! msh_no_BCs

            deallocate(tmpi_vec1)
            deallocate(tmpi_vec2)
            deallocate(tmps_vec)

            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$EndPhysicalNames') error stop "READ_GMSH :: Wrong input file - not all boundary conditions detected."
      !------------------------------------------------------------------------

      !-----Read-nodes---------------------------------------------------------
            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$Nodes') error stop "READ_GMSH :: Wrong input file - no nodes found."

            read(fUnit,*) numberOfNodes

            call msh_nodes % Construct(0, 0, .false., numberOfNodes)

            do i=1, numberOfNodes
                  read(fUnit,*) msh_nodes % tags(i), msh_nodes % cords(i,:)
            end do ! numberOfNodes

            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$EndNodes') error stop "READ_GMSH :: Wrong input file - not all nodes detected."
      !------------------------------------------------------------------------

      !-----Read-elements------------------------------------------------------
            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$Elements') error stop "READ_GMSH :: Wrong input file - no elements found."

            read(fUnit,*) numberOfElements
            call msh_elements % Construct(element_type,numberOfElements)
            allocate(msh_entity_vec(255)) ! arbitrary number
            allocate(el_types(numberOfElements)) ! arbitrary number
            do i=1, numberOfElements

               read(fUnit,'(4096a)') msh_entity ! read row

               msh_entity_vec=0.0_RP

               call getRealArrayFromStringNoCommas(msh_entity,msh_entity_vec)

               msh_elements % tags(i)     = int(msh_entity_vec(1))
               el_types(i)                = int(msh_entity_vec(2))
               msh_elements % no_ptags(i) = int(msh_entity_vec(3))
               if (msh_elements % no_ptags(i) .gt. 0) msh_elements % ptags(i,1:msh_elements % no_ptags(i)) = &
                  int(msh_entity_vec(4:3+msh_elements % no_ptags(i)))

                  select case (el_types(i))
                  case (5) ! 3D - 1st order
                     no_nodes_i = 8
                  case (12) ! 3D - 2nd order
                     no_nodes_i = 27
                  case (92) ! 3D - 3rd order
                     no_nodes_i = 64
                  case (93) ! 3D - 4th order
                     no_nodes_i = 125
                  case (94) ! 3D - 5th order
                     no_nodes_i = 216
                  case (3) ! 2D - 1st order
                     no_nodes_i = 4
                  case (10) ! 2D - 2nd order
                     no_nodes_i = 9
                  case (36) ! 2D - 3rd order
                     no_nodes_i = 16
                  case (37) ! 2D - 4th order
                     no_nodes_i = 25
                  case (38) ! 2D - 5th order
                     no_nodes_i = 36
                  case default
                     no_nodes_i = 0
                  end select

               msh_elements % nodes(i,1:size(msh_entity_vec(4+msh_elements % no_ptags(i):3+msh_elements % no_ptags(i)+no_nodes_i) )) = &
                int(msh_entity_vec(4+msh_elements % no_ptags(i):3+msh_elements % no_ptags(i)+no_nodes_i))

            end do ! numberOfElements

            deallocate(msh_entity_vec)

            read(fUnit,*) tmps
            if (trim(tmps) .ne. '$EndElements') error stop "READ_GMSH :: Wrong input file - not all elements detected."
            close( fUnit )
      !------------------------------------------------------------------------

      !-----Mesh-info----------------------------------------------------------
            ! find order of elements curvature
            check_eltype = 0
            do i=1, EL_MAX_ORDER
               check_eltype(i) = count(el_types .eq. SUPPORTED_EL_TYPES(i))
            end do
            if (sum(check_eltype) .eq. 0) error stop "READ_GMSH :: No 3D elements detected in the mesh."
            if (sum(check_eltype) .ne. maxval(check_eltype)) error stop "READ_GMSH :: More than 1 type of hexahedral detected in the mesh."
            bFaceOrder = maxloc(check_eltype,1) ! set order of the mesh
            org_element_type = SUPPORTED_EL_TYPES(bFaceOrder)
            msh_elements % el_type = org_element_type

            select case (bFaceOrder)
            case (1)
               org_element_type_2D = 3
            case (2)
               org_element_type_2D = 10
            case (3)
               org_element_type_2D = 36
            case (4)
               org_element_type_2D = 37
            case (5)
               org_element_type_2D = 38
            case default
            end select

            if (numberOfElements .ne. (count(el_types .eq. org_element_type) + count(el_types .eq. org_element_type_2D)) ) &
               error stop "READ_GMSH :: Too many different types of elements."
            numberOfElements = count(el_types .eq. org_element_type)
            numberOfElements2D = count(el_types .eq. org_element_type_2D)

            call msh_elements_3D % Construct(org_element_type,numberOfElements)
            call msh_elements_2D % Construct(org_element_type_2D,numberOfElements2D)

            j = 0
            k = 0
            do i=1, size(el_types)

               if ( el_types(i) .eq. org_element_type ) then
                  j = j + 1
                  msh_elements_3D % tags(j) = msh_elements % tags(i)
                  msh_elements_3D % no_ptags(j) = msh_elements % no_ptags(i)
                  msh_elements_3D % ptags(j,:) = msh_elements % ptags(i,1)
                  msh_elements_3D % nodes(j,:) = msh_elements % nodes(i,1:size(msh_elements_3D % nodes(j,:)))
               elseif ( el_types(i) .eq. org_element_type_2D ) then
                  k = k + 1
                  msh_elements_2D % tags(k) = msh_elements % tags(i)
                  msh_elements_2D % no_ptags(k) = msh_elements % no_ptags(i)
                  msh_elements_2D % ptags(k,:) = msh_elements % ptags(i,1)
                  msh_elements_2D % nodes(k,:) = msh_elements % nodes(i,1:size(msh_elements_2D % nodes(k,:)))
               else
                  error stop "READ_GMSH :: Unknown element type in the mesh."
               end if

            end do

            deallocate(el_types)
            call msh_elements % Destruct()
      !------------------------------------------------------------------------

      !-----Reorder-nodes-in-elements------------------------------------------
            allocate(tmpi_vec1((bFaceOrder + 1)**3))
            do j=1, msh_elements_3D % no_els
               ! re-order nodes within element to match HORSES ordering
               tmpi_vec1 = msh_elements_3D % nodes(j,:)
               call ReorderElement(tmpi_vec1,bFaceOrder)
               msh_elements_3D % nodes(j,:) = tmpi_vec1
            end do
            deallocate(tmpi_vec1)
      !------------------------------------------------------------------------

      !-----Assign-BCs---------------------------------------------------------
            do i=1,msh_no_BCs
               ! find # of surfaces, in this version of the reader surfaces are 2D element faces
               msh_bcs(i) % no_of_surfaces = 0
               allocate(tmpi_vec1(numberOfNodes*3)) ! upper bound
               tmpi_vec1 = 0
               k = 0

               do j=1, numberOfElements2D

                  tmpi = my_findloc(msh_elements_2D%ptags(j,:), msh_bcs(i)%tag, 1)
                  if (tmpi .gt. 0) then
                     msh_bcs(i) % no_of_surfaces = msh_bcs(i) % no_of_surfaces + 1
                     tmpi_vec1(k+1 : k+size(msh_elements_2D % nodes(j,:))) = msh_elements_2D % nodes(j,:)
                     k = k + size(msh_elements_2D % nodes(j,:))
                  end if

               end do ! msh_no_surfaces

               call unique(tmpi_vec1(1:k),msh_bcs(i) % node_tags)
               msh_bcs(i) % no_of_nodes = size(msh_bcs(i) % node_tags,1)
               deallocate(tmpi_vec1)

            end do ! msh_no_BCs
      !------------------------------------------------------------------------

      !-----Assign-BC-to-nodes-------------------------------------------------
            allocate(tmpi_vec1(8)) ! help array for the walls
            do i=1, msh_no_BCs
               do j=1, msh_elements_3D % no_els

                  ! find matching nodes
                  tmpi_vec1=0
                  do l=1,8
                     tmpi = my_findloc(msh_bcs(i) % node_tags,msh_elements_3D % nodes(j,l),1)
                     if(tmpi .gt. 0) tmpi_vec1(l) = 1
                     ! TODO: Check why this doesn't work with older ifort
                     ! tmpi = any(msh_bcs(i) % node_tags .eq. msh_element_blocks(msh_elblock) % nodes(j,l))
                     ! tmpi_vec1(l)=tmpi
                  end do

                  ! assign BC face to the elements nodes
                  if ( (sum(tmpi_vec1(1:2)+tmpi_vec1(5:6))) .eq. 4) msh_elements_3D % BCs(j,1) = msh_bcs(i) % tag ! 1:west
                  if ( (sum(tmpi_vec1(3:4)+tmpi_vec1(7:8))) .eq. 4) msh_elements_3D % BCs(j,2) = msh_bcs(i) % tag ! 2:east
                  if ( (sum(tmpi_vec1(1:4))) .eq. 4) msh_elements_3D % BCs(j,3) = msh_bcs(i) % tag ! 3:south
                  if ( (sum(tmpi_vec1(2:3)+tmpi_vec1(6:7))) .eq. 4) msh_elements_3D % BCs(j,4) = msh_bcs(i) % tag ! 4:front
                  if ( (sum(tmpi_vec1(5:8))) .eq. 4) msh_elements_3D % BCs(j,5) = msh_bcs(i) % tag ! 5:north
                  if ( (tmpi_vec1(1)+tmpi_vec1(4)+tmpi_vec1(5)+tmpi_vec1(8)) .eq. 4) msh_elements_3D % BCs(j,6) = msh_bcs(i) % tag ! 6:back

               end do ! msh_elements_3D % no_els
            end do ! msh_no_BCs
            deallocate(tmpi_vec1)
      !------------------------------------------------------------------------

      !-----Assign-new-tags----------------------------------------------------
            tmp_eltag = 0
            do j=1, msh_elements_3D % no_els
               tmp_eltag = tmp_eltag + 1
               msh_elements_3D % tags(j) = tmp_eltag
            end do
            if (.not. (tmp_eltag .eq. numberOfElements)) error stop "Read_GMSH :: Number of elements inconsistent."
      !------------------------------------------------------------------------

      !-----Build-nodes--------------------------------------------------------
            self % nodeType = nodes
            self % no_of_elements = numberOfElements
      !------------------------------------------------------------------------

      !-----Set-up-face-patches------------------------------------------------
            numBFacePoints = bFaceOrder + 1
            allocate(uNodes(numBFacePoints))
            allocate(vNodes(numBFacePoints))
            allocate(values(3,numBFacePoints,numBFacePoints))
            do i = 1, numBFacePoints
               uNodes(i) = -1._RP + (i-1) * (2._RP/bFaceOrder)
               vNodes(i) = uNodes(i)
            end do
      !------------------------------------------------------------------------

      !-----Allocate-mem-for-elements-and-nodes--------------------------------
            allocate( self % elements(numberOfelements) )
            allocate( self % nodes(numberOfNodes) )
            allocate( self % Nx(numberOfelements) , self % Ny(numberOfelements) , self % Nz(numberOfelements) )
            self % Nx = Nx
            self % Ny = Ny
            self % Nz = Nz
      !------------------------------------------------------------------------

      !----Set-nodes-----------------------------------------------------------
            do msh_node=1, msh_nodes % no_nodes
               x = msh_nodes % cords(msh_node,1:NDIM)/Lref
               call ConstructNode( self % nodes(msh_nodes % tags(msh_node)), x, msh_nodes % tags(msh_node) )
            end do ! msh_nodes % no_nodes
      !------------------------------------------------------------------------

      !----Set-elements-----------------------------------------------------------
            call extract_bc_tags(msh_bcs, tags)
            j = 0
            do msh_el = 1, msh_elements_3D % no_els
               j = j + 1
               ! setting l'th element
               l =  msh_elements_3D % tags(msh_el)
               nodeIDs = msh_elements_3D % nodes(msh_el,1:8) ! only non-curved nodes

               if (bFaceOrder .eq. 1) then ! non-curved mesh
                  do k = 1, NODES_PER_ELEMENT
                     corners(:,k) = self % nodes(nodeIDs(k)) % x
                  end do
                  self % elements(l) % SurfInfo % IsHex8 = .TRUE.
                  self % elements(l) % SurfInfo % corners = corners
               else ! curved mesh
                  ! allocate tmp arrays for curved face node tags and coordinates
                  allocate(tmpi_vec1(numBFacePoints*numBFacePoints))
                  allocate(tmpi_vec2(numBFacePoints**3))

                  tmpi_vec2 = msh_elements_3D % nodes(msh_el,:)

                  do k = 1, FACES_PER_ELEMENT
                     call GetOrderedFaceNodeTags(tmpi_vec1,k,bFaceOrder,tmpi_vec2)
                     do jj = 1, numBFacePoints
                        do ii = 1, numBFacePoints
                           values(:,ii,jj) = self % nodes(tmpi_vec1(ii + (jj-1)*numBFacePoints)) % x
                        end do
                     end do
                     values = values / Lref
                     call self % elements(l) % SurfInfo % facePatches(k) % construct(uNodes, vNodes, values)

                  end do

                  deallocate(tmpi_vec1)
                  deallocate(tmpi_vec2)
               end if

               call self % elements(l) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , l, l)

               ! set element boundaries
               do k = 1, 6
                  tmpi1 = my_findloc(tags, msh_elements_3D % BCs(msh_el,k),1)
                  if (tmpi1 .gt. 0) then
                     self % elements(l) % boundaryName(k) = trim(msh_bcs(tmpi1) % name)
                  else
                     self % elements(l) % boundaryName(k) = emptyBCName
                  end if
               end do ! k

               ! set BC names to faces
               do k = 1, 6
                  if(trim(self % elements(l) % boundaryName(k)) /= emptyBCName) then
                     call toLower( self % elements(l) % boundaryName(k) )
                     numberOfBoundaryFaces = numberOfBoundaryFaces + 1
                     zoneNames => zoneNameDictionary % allKeys()
                     if ( all(trim(self % elements(l) % boundaryName(k)) .ne. zoneNames) ) then
                        call zoneNameDictionary % addValueForKey(trim(self % elements(l) % boundaryName(k)), trim(self % elements(l) % boundaryName(k)))
                     end if
                     deallocate (zoneNames)
                  end if
               end do ! k
            end do ! msh_elements_3D % no_els
            if (.not. (j .eq. numberOfElements)) error stop "Read_GMSH :: Not all elements assigned."
      !------------------------------------------------------------------------

      !-----Deallocate-msh-vairables-------------------------------------------
            do i=1,msh_no_BCs
               call msh_bcs(i) % Destruct()
            end do
            deallocate(msh_bcs)

            call msh_nodes % Destruct()
            call msh_elements_3D % Destruct()
            call msh_elements_2D % Destruct()

            deallocate(uNodes) ! Check if we can do it! FIXME
            deallocate(vNodes) ! Check if we can do it! FIXME
            deallocate(values) ! Check if we can do it! FIXME
      !------------------------------------------------------------------------
      !
      !     ---------------------------
      !     Construct the element faces
      !     ---------------------------
      !
            numberOfFaces        = (6*numberOfElements + numberOfBoundaryFaces)/2
            self % numberOfFaces = numberOfFaces
            allocate( self % faces(self % numberOfFaces) )
            CALL ConstructFaces( self, success )
      !
      !     -------------------------
      !     Build the different zones
      !     -------------------------
      !
            call self % ConstructZones()
      !
      !     ---------------------------
      !     Construct periodic faces
      !     ---------------------------
      !
            CALL ConstructPeriodicFaces( self, periodRelative )
      !
      !     ---------------------------
      !     Delete periodic- faces
      !     ---------------------------
      !
            CALL DeletePeriodicMinusFaces( self )
      !
      !     ---------------------------
      !     Assign faces ID to elements
      !     ---------------------------
      !
            CALL getElementsFaceIDs(self)
      !
      !     ---------------------
      !     Define boundary faces
      !     ---------------------
      !
            call self % DefineAsBoundaryFaces()
      !
      !     ------------------------------
      !     Set the element connectivities
      !     ------------------------------
      !
            call self % SetConnectivitiesAndLinkFaces(nodes)

            call self % ExportBoundaryMesh (trim(fileName))
         end subroutine ConstructSimplestMesh_FromGMSHFile_v2_
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
    function NumOfElems_GMSH_v4( fileName ) result(nelem)
        implicit none
        !----------------------------------
        CHARACTER(LEN=*), intent(in) :: fileName
        integer                      :: nelem
        !----------------------------------
        integer :: i, j, k, fUnit, nElBlocks, ElsInBlock, eltype, elem_tmp
        integer :: ierr
        logical :: success=.false.
        character (len=100) :: cline
        character (len=10)  :: cword
        !----------------------------------

        open(newunit = fUnit, FILE = trim(fileName) )

        nelem = 0
        do
            read(fUnit,'(A)',iostat=ierr) cline
            if (ierr /= 0) exit
            read (cline,*) cword ! read first word of line
            if (cword == "$Elements") then
               read(fUnit,*) nElBlocks, elem_tmp, k, k
               do i=1, nElBlocks
                  read(fUnit,*) k, k, eltype, ElsInBlock
                  if (my_findloc(SUPPORTED_EL_TYPES,eltype,1) .ge. 1) nelem = nelem + ElsInBlock
                  do j=1, ElsInBlock
                     read(fUnit,*) k
                  end do
               end do
               read(fUnit,*) cline
               if (trim(cline) .ne. '$EndElements') error stop "NumOfElems_GMSH :: Wrong input file - not all elements read."
               success=.true.
               exit
            end if
        end do

        if (success) then
        else
            print *, cline
            error stop "NumOfElems_GMSH :: Could not read number of elements."
        end if

        close(fUnit)

    end function NumOfElems_GMSH_v4
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
    function NumOfElems_GMSH_v2( fileName ) result(nelem)
      implicit none
      !----------------------------------
      CHARACTER(LEN=*), intent(in) :: fileName
      integer                      :: nelem
      !----------------------------------
      integer :: i, j, k, fUnit, nEltotal, ElsInBlock, eltype, elem_tmp
      integer :: ierr
      logical :: success=.false.
      character (len=100) :: cline
      character (len=10)  :: cword
      !----------------------------------

      open(newunit = fUnit, FILE = trim(fileName) )

      nelem = 0
      do
          read(fUnit,'(A)',iostat=ierr) cline
          if (ierr /= 0) exit
          read (cline,*) cword ! read first word of line
          if (cword == "$Elements") then
             read(fUnit,*) nEltotal
             do i=1, nEltotal
               read(fUnit,*) k, eltype
               if (my_findloc(SUPPORTED_EL_TYPES,eltype,1) .ge. 1) nelem = nelem + 1
             end do
             read(fUnit,*) cline
             if (trim(cline) .ne. '$EndElements') error stop "NumOfElems_GMSH :: Wrong input file - not all elements read."
             success=.true.
             exit
          end if
      end do

      if (success) then
      else
         print *, cline
         error stop "NumOfElems_GMSH :: Could not read number of elements."
      end if

      close(fUnit)

  end function NumOfElems_GMSH_v2
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MSH_DestructBCStorage(this)
!  ---------------------------------------------------------
!  Constructor.
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(MSH_BCinfo_t)            , intent(inout)   :: this
!-----Local-Variables---------------------------------------------
!  -----------------------------------------------------------------------
      if(allocated(this % element_tags))  deallocate(this % element_tags)
      if(allocated(this % wall_side))     deallocate(this % wall_side)
      if(allocated(this % node_tags))     deallocate(this % node_tags)
      if(allocated(this % volume_tags))   deallocate(this % volume_tags)
      if(allocated(this % surface_tags))  deallocate(this % surface_tags)
      if(allocated(this % curve_tags))    deallocate(this % curve_tags)
      if(allocated(this % point_tags))    deallocate(this % point_tags)

   end subroutine MSH_DestructBCStorage
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MSH_ConstructNodeBlock(this,edim,etag,par,no_nodes)
!  ---------------------------------------------------------
!  Constructor.
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(MSH_node_block_t)            , intent(inout)   :: this ! element block to be constructed
      integer                , intent(in)      :: edim ! dimension of the entity
      integer                , intent(in)      :: etag ! entity tag
      logical                , intent(in)      :: par ! parametric or not (TODO: parametric not supported!)
      integer                , intent(in)      :: no_nodes ! number of nodes within the block
!-----Local-Variables---------------------------------------------
!  -----------------------------------------------------------------------
      this % entity_dim = edim
      this % entity_tag = etag
      this % parametric = par
      this % no_nodes   = no_nodes
      allocate(this % tags(no_nodes))
      allocate(this % cords(no_nodes,3))

   end subroutine MSH_ConstructNodeBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MSH_DestructNodeBlock(this)
!  ---------------------------------------------------------
!  Constructor.
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(MSH_node_block_t)            , intent(inout)   :: this ! element block to be constructed
!-----Local-Variables---------------------------------------------
!  -----------------------------------------------------------------------
      deallocate(this % tags)
      deallocate(this % cords)

   end subroutine MSH_DestructNodeBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MSH_ConstructElementBlock(this,el_type,no_els)
!  ---------------------------------------------------------
!  Constructor.
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(MSH_element_block_t)            , intent(inout)   :: this ! element block to be constructed
      integer                , intent(in)      :: el_type ! type of element
      integer                , intent(in)      :: no_els ! number of elements
!-----Local-Variables---------------------------------------------
!  -----------------------------------------------------------------------
      this % el_type = el_type
      this % no_els = no_els
      allocate(this % tags(no_els))
      select case (el_type)
      case (3)  ! 2D 1st order quadrangle
         allocate(this % nodes(no_els,4))
      case (10) ! 2D 2nd order quadrangle
         allocate(this % nodes(no_els,9))
      case (36) ! 2D 3rd order quadrangle
         allocate(this % nodes(no_els,16))
      case (37) ! 2D 4th order quadrangle
         allocate(this % nodes(no_els,25))
      case (38) ! 2D 5th order quadrangle
         allocate(this % nodes(no_els,36))
      case (5)  ! 3D 1st order hexahedron
         allocate(this % nodes(no_els,8))
      case (12) ! 3D 2nd order hexahedron
         allocate(this % nodes(no_els,27))
      case (92) ! 3D 3rd order hexahedron
         allocate(this % nodes(no_els,64))
      case (93) ! 3D 4th order hexahedron
         allocate(this % nodes(no_els,125))
      case (94) ! 3D 5th order hexahedron
         allocate(this % nodes(no_els,216))
      case default
         ! print *, " READ_GMSH :: Warning! Wrong element type. Allocation for Q3 3D."
         allocate(this % nodes(no_els,216))
      end select

      ! allocate BCs only if 3D element is detected
      if (my_findloc(SUPPORTED_EL_TYPES,el_type,1) .ge. 1) then
         allocate(this % BCs(no_els,6))
         this % BCs = 0
      end if

      allocate(this % no_ptags(no_els))
      allocate(this % ptags(no_els,8))
   end subroutine MSH_ConstructElementBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MSH_DestructElementBlock(this)
!  ---------------------------------------------------------
!  Constructor.
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(MSH_element_block_t)            , intent(inout)   :: this ! element block to be constructed
!-----Local-Variables---------------------------------------------
!  -----------------------------------------------------------------------
      deallocate(this % tags)
      deallocate(this % nodes)
      deallocate(this % no_ptags)
      deallocate(this % ptags)
      if(allocated(this % BCs)) deallocate(this % BCs)

   end subroutine MSH_DestructElementBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine unique(vec,vec_unique)
      ! CREDIT: http://degenerateconic.com/unique/ email: jacob[at]degenerateconic[dot]com
      ! Return only the unique values from vec.

      implicit none

      integer,dimension(:),intent(in) :: vec
      integer,dimension(:),allocatable,intent(inout) :: vec_unique

      integer :: i,num
      logical,dimension(size(vec)) :: mask

      mask = .false.

      do i=1,size(vec)

          !count the number of occurrences of this element:
          num = count( vec(i)==vec )

          if (num==1) then
              !there is only one, flag it:
              mask(i) = .true.
          else
              !flag this value only if it hasn't already been flagged:
              if (.not. any(vec(i)==vec .and. mask) ) mask(i) = .true.
          end if

      end do

      !return only flagged elements:
      allocate( vec_unique(count(mask)) )
      vec_unique = pack( vec, mask )

      !if you also need it sorted, then do so.
      ! For example, with slatec routine:
      !call ISORT (vec_unique, [0], size(vec_unique), 1)

   end subroutine unique
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine GetOrderedFaceNodeTags(face_nodes,k,bFaceOrder,elnodes)
!  ---------------------------------------------------------
!  Routine to extract k'th face from the element.
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      integer, allocatable, dimension(:),    intent(inout)   :: face_nodes !
      integer,                               intent(in   )   :: k ! face number
      integer,                               intent(in   )   :: bFaceOrder ! polynomial order
      integer, allocatable, dimension(:),    intent(in   )   :: elnodes !
!-----Local-Variables---------------------------------------------
      integer :: numBFacePoints, innerEdgePoints
      integer :: i
!  -----------------------------------------------------------------------

      face_nodes = 0.d0
      numBFacePoints = bFaceOrder + 1
      innerEdgePoints = bFaceOrder - 1

      select case (k)
      case(1)
         ! Corners
         face_nodes(1) = elnodes( 1 )
         face_nodes(numBFacePoints) = elnodes( 2 )
         face_nodes(numBFacePoints*bFaceOrder+1) = elnodes( 5 )
         face_nodes(numBFacePoints*numBFacePoints) = elnodes( 6 )

         ! Edges
         face_nodes(2 : 1+innerEdgePoints) = &
         elnodes( 9+innerEdgePoints*1 : 9+innerEdgePoints*2-1 ) ! 10th entity

         face_nodes( (numBFacePoints + 1):(innerEdgePoints*numBFacePoints + 1):numBFacePoints ) = &
         elnodes( 9+innerEdgePoints*7 : 9+innerEdgePoints*8-1 ) ! 16th entity

         face_nodes( (numBFacePoints*2):(numBFacePoints*bFaceOrder):numBFacePoints ) = &
         elnodes( 9+innerEdgePoints*2 : 9+innerEdgePoints*3-1 ) ! 11th entity

         face_nodes( numBFacePoints*bFaceOrder + 2 : numBFacePoints*bFaceOrder + 1 + innerEdgePoints ) = &
         elnodes( 9+innerEdgePoints*9 : 9+innerEdgePoints*10-1 ) ! 18th entity

         ! Inner part of the face
         do i=1, innerEdgePoints
            face_nodes( numBFacePoints*i + 2 : numBFacePoints*i + 1 + innerEdgePoints ) = &
            elnodes( 1 +  (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) : &
                              i*innerEdgePoints + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) ! 23rd entity
         end do
      case(2)
         ! Corners
         face_nodes(1) = elnodes( 4 )
         face_nodes(numBFacePoints) = elnodes( 3 )
         face_nodes(numBFacePoints*bFaceOrder+1) = elnodes( 8 )
         face_nodes(numBFacePoints*numBFacePoints) = elnodes( 7 )

         ! Edges
         face_nodes(2 : 1+innerEdgePoints) = &
         elnodes( 9+innerEdgePoints*3 : 9+innerEdgePoints*4-1 ) ! 12th entity

         face_nodes( (numBFacePoints + 1):(innerEdgePoints*numBFacePoints + 1):numBFacePoints ) = &
         elnodes( 9+innerEdgePoints*6 : 9+innerEdgePoints*7-1 ) ! 15th entity

         face_nodes( (numBFacePoints*2):(numBFacePoints*bFaceOrder):numBFacePoints ) = &
         elnodes( 9+innerEdgePoints*4 : 9+innerEdgePoints*5-1 ) ! 13th entity

         face_nodes( numBFacePoints*bFaceOrder + 2 : numBFacePoints*bFaceOrder + 1 + innerEdgePoints ) = &
         elnodes( 9+innerEdgePoints*10 : 9+innerEdgePoints*11-1 ) ! 19th entity

         ! Inner part of the face
         do i=1, innerEdgePoints
            face_nodes( numBFacePoints*i + 2 : numBFacePoints*i + 1 + innerEdgePoints ) = &
            elnodes( 1 +  (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) : &
                              i*innerEdgePoints + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) ! 24th entity
         end do
      case(3)
         ! Corners
         face_nodes(1) = elnodes( 1 )
         face_nodes(numBFacePoints) = elnodes( 2 )
         face_nodes(numBFacePoints*bFaceOrder+1) = elnodes( 4 )
         face_nodes(numBFacePoints*numBFacePoints) = elnodes( 3 )

         ! Edges
         face_nodes(2 : 1+innerEdgePoints) = &
         elnodes( 9+innerEdgePoints*1 : 9+innerEdgePoints*2-1 ) ! 10th entity

         face_nodes( (numBFacePoints + 1):(innerEdgePoints*numBFacePoints + 1):numBFacePoints ) = &
         elnodes( 9+innerEdgePoints*5 : 9+innerEdgePoints*6-1 ) ! 14th entity

         face_nodes( (numBFacePoints*2):(numBFacePoints*bFaceOrder):numBFacePoints ) = &
         elnodes( 9+innerEdgePoints*0 : 9+innerEdgePoints*1-1 ) ! 9th entity

         face_nodes( numBFacePoints*bFaceOrder + 2 : numBFacePoints*bFaceOrder + 1 + innerEdgePoints ) = &
         elnodes( 9+innerEdgePoints*3 : 9+innerEdgePoints*4-1 ) ! 12th entity

         ! Inner part of the face
         do i=1, innerEdgePoints
            face_nodes( numBFacePoints*i + 2 : numBFacePoints*i + 1 + innerEdgePoints ) = &
            elnodes( 1 +  (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) : &
                              i*innerEdgePoints + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) ! 21st entity
         end do
      case(4)
         ! Corners
         face_nodes(1) = elnodes( 2 )
         face_nodes(numBFacePoints) = elnodes( 3 )
         face_nodes(numBFacePoints*bFaceOrder+1) = elnodes( 6 )
         face_nodes(numBFacePoints*numBFacePoints) = elnodes( 7 )

         ! Edges
         face_nodes(2 : 1+innerEdgePoints) = &
         elnodes( 9+innerEdgePoints*0 : 9+innerEdgePoints*1-1 ) ! 9th entity

         face_nodes( (numBFacePoints + 1):(innerEdgePoints*numBFacePoints + 1):numBFacePoints ) = &
         elnodes( 9+innerEdgePoints*2 : 9+innerEdgePoints*3-1 ) ! 11th entity

         face_nodes( (numBFacePoints*2):(numBFacePoints*bFaceOrder):numBFacePoints ) = &
         elnodes( 9+innerEdgePoints*4 : 9+innerEdgePoints*5-1 ) ! 13th entity

         face_nodes( numBFacePoints*bFaceOrder + 2 : numBFacePoints*bFaceOrder + 1 + innerEdgePoints ) = &
         elnodes( 9+innerEdgePoints*8 : 9+innerEdgePoints*9-1 ) ! 17th entity

         ! Inner part of the face
         do i=1, innerEdgePoints
            face_nodes( numBFacePoints*i + 2 : numBFacePoints*i + 1 + innerEdgePoints ) = &
            elnodes( 1 +  (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) : &
                              i*innerEdgePoints + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) ! 22nd entity
         end do
      case(5)
         ! Corners
         face_nodes(1) = elnodes( 5 )
         face_nodes(numBFacePoints) = elnodes( 6 )
         face_nodes(numBFacePoints*bFaceOrder+1) = elnodes( 8 )
         face_nodes(numBFacePoints*numBFacePoints) = elnodes( 7 )

         ! Edges
         face_nodes(2 : 1+innerEdgePoints) = &
         elnodes( 9+innerEdgePoints*9 : 9+innerEdgePoints*10-1 ) ! 18th entity

         face_nodes( (numBFacePoints + 1):(innerEdgePoints*numBFacePoints + 1):numBFacePoints ) = &
         elnodes( 9+innerEdgePoints*11 : 9+innerEdgePoints*12-1 ) ! 20th entity

         face_nodes( (numBFacePoints*2):(numBFacePoints*bFaceOrder):numBFacePoints ) = &
         elnodes( 9+innerEdgePoints*8 : 9+innerEdgePoints*9-1 ) ! 17th entity

         face_nodes( numBFacePoints*bFaceOrder + 2 : numBFacePoints*bFaceOrder + 1 + innerEdgePoints ) = &
         elnodes( 9+innerEdgePoints*10 : 9+innerEdgePoints*11-1 ) ! 19th entity

         ! Inner part of the face
         do i=1, innerEdgePoints
            face_nodes( numBFacePoints*i + 2 : numBFacePoints*i + 1 + innerEdgePoints ) = &
            elnodes( 1 +  (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) : &
                              i*innerEdgePoints + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) ! 26th entity
         end do
      case(6)
         ! Corners
         face_nodes(1) = elnodes( 1 )
         face_nodes(numBFacePoints) = elnodes( 4 )
         face_nodes(numBFacePoints*bFaceOrder+1) = elnodes( 5 )
         face_nodes(numBFacePoints*numBFacePoints) = elnodes( 8 )

         ! Edges
         face_nodes(2 : 1+innerEdgePoints) = &
         elnodes( 9+innerEdgePoints*5 : 9+innerEdgePoints*6-1 ) ! 14th entity

         face_nodes( (numBFacePoints + 1):(innerEdgePoints*numBFacePoints + 1):numBFacePoints ) = &
         elnodes( 9+innerEdgePoints*7 : 9+innerEdgePoints*8-1 ) ! 16th entity

         face_nodes( (numBFacePoints*2):(numBFacePoints*bFaceOrder):numBFacePoints ) = &
         elnodes( 9+innerEdgePoints*6 : 9+innerEdgePoints*7-1 ) ! 15th entity

         face_nodes( numBFacePoints*bFaceOrder + 2 : numBFacePoints*bFaceOrder + 1 + innerEdgePoints ) = &
         elnodes( 9+innerEdgePoints*11 : 9+innerEdgePoints*12-1 ) ! 20th entity


         ! Inner part of the face
         do i=1, innerEdgePoints
            face_nodes( numBFacePoints*i + 2 : numBFacePoints*i + 1 + innerEdgePoints ) = &
            elnodes( 1 +  (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) : &
                              i*innerEdgePoints + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) ! 25th entity
         end do
      case default
         error stop "Read_GMSH :: GetOrderedFaceNodeTags :: face number not in 1,6 range."
      end select

   end subroutine GetOrderedFaceNodeTags
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ReorderElement(elnodes,bFaceOrder)
!  ---------------------------------------------------------
!  Reordering routine from GMSH ordering to SpecMesh ordering.
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      integer, allocatable, dimension(:),    intent(inout)   :: elnodes !
      integer,                               intent(in   )   :: bFaceOrder ! polynomial order
!-----Local-Variables---------------------------------------------
      integer, dimension( (bFaceOrder+1)**3 ) :: buffer
      integer, dimension( (bFaceOrder-1),(bFaceOrder-1) ) :: face_buffer
      integer :: innerEdgePoints,i
!  -----------------------------------------------------------------------
      buffer = elnodes
      innerEdgePoints = bFaceOrder - 1

      ! reorder verices
      elnodes(1) = buffer(4)
      elnodes(2:4) = buffer(1:3)
      elnodes(5) = buffer(8)
      elnodes(6:8) = buffer(5:7)

      select case (bFaceOrder)
      case (1)
         ! do nothing
      case (2)
         ! do nothing
      case (3)
         ! reverse edges (reverse edge no 2,4,6,10,11,12)

         ! 10th entity (2nd edge)
         elnodes( 9+innerEdgePoints*1 : 9+innerEdgePoints*2-1 )   = buffer( 9+innerEdgePoints*2-1 : 9+innerEdgePoints*1 : -1 )
         ! 12th entity (4th edge)
         elnodes( 9+innerEdgePoints*3 : 9+innerEdgePoints*4-1 )   = buffer( 9+innerEdgePoints*4-1 : 9+innerEdgePoints*3 : -1 )
         ! 14th entity (6th edge)
         elnodes( 9+innerEdgePoints*5 : 9+innerEdgePoints*6-1 )   = buffer( 9+innerEdgePoints*6-1 : 9+innerEdgePoints*5 : -1 )
         ! 18th entity (10th edge)
         elnodes( 9+innerEdgePoints*9 : 9+innerEdgePoints*10-1 )  = buffer( 9+innerEdgePoints*10-1 : 9+innerEdgePoints*9 : -1 )
         ! 19th entity (11th edge)
         elnodes( 9+innerEdgePoints*10 : 9+innerEdgePoints*11-1 ) = buffer( 9+innerEdgePoints*11-1 : 9+innerEdgePoints*10 : -1 )
         ! 20th entity (12th edge)
         elnodes( 9+innerEdgePoints*11 : 9+innerEdgePoints*12-1 ) = buffer( 9+innerEdgePoints*12-1 : 9+innerEdgePoints*11 : -1 )

         ! 23rd entity (Face no 1)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes( innerEdgePoints + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(1,1)
         elnodes( 1 + 1*innerEdgePoints + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(2,1)
         elnodes( 2*innerEdgePoints + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(1,2)

         ! 24th entity (Face no 2)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(1,2)
         elnodes( innerEdgePoints + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(1,1)

         ! 21st entity (Face no 3)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(1,2)
         elnodes( innerEdgePoints + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(1,1)

         ! 22nd enitty (Face no 4)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 1*innerEdgePoints + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes( 2*innerEdgePoints + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(2,1)

         ! 26th enitty (Face no 5)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes( innerEdgePoints + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(1,1)
         elnodes( 1 + 1*innerEdgePoints + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(2,1)
         elnodes( 2*innerEdgePoints + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(1,2)

         ! 25th enitty (Face no 6)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(1,2)
         elnodes( innerEdgePoints + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(1,1)
      case (4)
         ! reverse edges (reverse edge no 2,4,6,10,11,12)

         ! 10th entity (2nd edge)
         elnodes( 9+innerEdgePoints*1 : 9+innerEdgePoints*2-1 )   = buffer( 9+innerEdgePoints*2-1 : 9+innerEdgePoints*1 : -1 )
         ! 12th entity (4th edge)
         elnodes( 9+innerEdgePoints*3 : 9+innerEdgePoints*4-1 )   = buffer( 9+innerEdgePoints*4-1 : 9+innerEdgePoints*3 : -1 )
         ! 14th entity (6th edge)
         elnodes( 9+innerEdgePoints*5 : 9+innerEdgePoints*6-1 )   = buffer( 9+innerEdgePoints*6-1 : 9+innerEdgePoints*5 : -1 )
         ! 18th entity (10th edge)
         elnodes( 9+innerEdgePoints*9 : 9+innerEdgePoints*10-1 )  = buffer( 9+innerEdgePoints*10-1 : 9+innerEdgePoints*9 : -1 )
         ! 19th entity (11th edge)
         elnodes( 9+innerEdgePoints*10 : 9+innerEdgePoints*11-1 ) = buffer( 9+innerEdgePoints*11-1 : 9+innerEdgePoints*10 : -1 )
         ! 20th entity (12th edge)
         elnodes( 9+innerEdgePoints*11 : 9+innerEdgePoints*12-1 ) = buffer( 9+innerEdgePoints*12-1 : 9+innerEdgePoints*11 : -1 )

         ! 23rd entity (Face no 1)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(2,1)
         elnodes( 2 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(3,2)
         elnodes( 3 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(1,1)
         elnodes( 4 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(3,1)
         elnodes( 5 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(3,3)
         elnodes( 6 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes( 7 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(1,3)
         elnodes( 8 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(2,3)
         elnodes( 9 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(1,2)

         ! 24th entity (Face no 2)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(1,2)
         elnodes( 2 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes( 3 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(1,1)
         elnodes( 4 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(2,3)
         elnodes( 5 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(3,3)
         elnodes( 6 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(3,2)
         elnodes( 7 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(1,3)
         elnodes( 8 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(3,1)
         elnodes( 9 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(2,1)

         ! 21st entity (Face no 3)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(1,2)
         elnodes( 2 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes( 3 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(1,1)
         elnodes( 4 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(2,3)
         elnodes( 5 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(3,3)
         elnodes( 6 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(3,2)
         elnodes( 7 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(1,3)
         elnodes( 8 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(3,1)
         elnodes( 9 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(2,1)

         ! 22nd enitty (Face no 4)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) )
         end do
         elnodes( 2 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes( 3 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(1,2)
         elnodes( 4 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(3,2)
         elnodes( 5 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(3,3)
         elnodes( 7 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(2,1)
         elnodes( 8 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(3,1)
         elnodes( 9 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(1,3)

         ! 26th enitty (Face no 5)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(2,1)
         elnodes( 2 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(3,2)
         elnodes( 3 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(1,1)
         elnodes( 4 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(3,1)
         elnodes( 5 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(3,3)
         elnodes( 6 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes( 7 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(1,3)
         elnodes( 8 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(2,3)
         elnodes( 9 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(1,2)

         ! 25th enitty (Face no 6)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(1,2)
         elnodes( 2 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes( 3 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(1,1)
         elnodes( 4 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(2,3)
         elnodes( 5 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(3,3)
         elnodes( 6 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(3,2)
         elnodes( 7 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(1,3)
         elnodes( 8 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(3,1)
         elnodes( 9 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(2,1)
      case (5)
         ! reverse edges (reverse edge no 2,4,6,10,11,12)

         ! 10th entity (2nd edge)
         elnodes( 9+innerEdgePoints*1 : 9+innerEdgePoints*2-1 )   = buffer( 9+innerEdgePoints*2-1 : 9+innerEdgePoints*1 : -1 )
         ! 12th entity (4th edge)
         elnodes( 9+innerEdgePoints*3 : 9+innerEdgePoints*4-1 )   = buffer( 9+innerEdgePoints*4-1 : 9+innerEdgePoints*3 : -1 )
         ! 14th entity (6th edge)
         elnodes( 9+innerEdgePoints*5 : 9+innerEdgePoints*6-1 )   = buffer( 9+innerEdgePoints*6-1 : 9+innerEdgePoints*5 : -1 )
         ! 18th entity (10th edge)
         elnodes( 9+innerEdgePoints*9 : 9+innerEdgePoints*10-1 )  = buffer( 9+innerEdgePoints*10-1 : 9+innerEdgePoints*9 : -1 )
         ! 19th entity (11th edge)
         elnodes( 9+innerEdgePoints*10 : 9+innerEdgePoints*11-1 ) = buffer( 9+innerEdgePoints*11-1 : 9+innerEdgePoints*10 : -1 )
         ! 20th entity (12th edge)
         elnodes( 9+innerEdgePoints*11 : 9+innerEdgePoints*12-1 ) = buffer( 9+innerEdgePoints*12-1 : 9+innerEdgePoints*11 : -1 )

         ! 23rd entity (Face no 1)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(1,4)
         elnodes( 2 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(3,3)
         elnodes( 3 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(3,4)
         elnodes( 4 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(1,1)
         elnodes( 5 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(3,2)
         elnodes( 6 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(4,4)
         elnodes( 7 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(4,1)
         elnodes( 8 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(2,1)
         elnodes( 9 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(3,1)
         elnodes(10 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(4,3)
         elnodes(11 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(4,2)
         elnodes(12 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes(13 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(1,3)
         elnodes(14 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(2,4)
         elnodes(15 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(2,3)
         elnodes(16 + 8 + 12*innerEdgePoints + 2*(innerEdgePoints**2) ) = face_buffer(1,2)

         ! 24th entity (Face no 2)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(1,2)
         elnodes( 2 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes( 3 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(2,1)
         elnodes( 4 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(1,1)
         elnodes( 5 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(2,3)
         elnodes( 6 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(4,2)
         elnodes( 7 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(4,1)
         elnodes( 8 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(3,4)
         elnodes( 9 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(2,4)
         elnodes(10 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(4,3)
         elnodes(11 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(4,4)
         elnodes(12 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(3,3)
         elnodes(13 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(1,3)
         elnodes(14 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(3,1)
         elnodes(15 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(3,2)
         elnodes(16 + 8 + 12*innerEdgePoints + 3*(innerEdgePoints**2) ) = face_buffer(1,4)

         ! 21st entity (Face no 3)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(1,2)
         elnodes( 2 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes( 3 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(2,1)
         elnodes( 4 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(1,1)
         elnodes( 5 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(2,3)
         elnodes( 6 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(4,2)
         elnodes( 7 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(4,1)
         elnodes( 8 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(3,4)
         elnodes( 9 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(2,4)
         elnodes(10 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(4,3)
         elnodes(11 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(4,4)
         elnodes(12 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(3,3)
         elnodes(13 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(1,3)
         elnodes(14 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(3,1)
         elnodes(15 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(3,2)
         elnodes(16 + 8 + 12*innerEdgePoints + 0*(innerEdgePoints**2) ) = face_buffer(1,4)

         ! 22nd enitty (Face no 4)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(1,1)
         elnodes( 2 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(2,1)
         elnodes( 3 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes( 4 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(1,2)
         elnodes( 5 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(3,4)
         elnodes( 6 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(4,1)
         elnodes( 7 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(4,2)
         elnodes( 8 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(2,3)
         elnodes( 9 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(3,3)
         elnodes(10 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(4,4)
         elnodes(11 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(4,3)
         elnodes(12 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(2,4)
         elnodes(13 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(1,4)
         elnodes(14 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(3,2)
         elnodes(15 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(3,1)
         elnodes(16 + 8 + 12*innerEdgePoints + 1*(innerEdgePoints**2) ) = face_buffer(1,3)

         ! 26th enitty (Face no 5)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(1,4)
         elnodes( 2 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(3,3)
         elnodes( 3 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(3,4)
         elnodes( 4 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(1,1)
         elnodes( 5 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(3,2)
         elnodes( 6 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(4,4)
         elnodes( 7 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(4,1)
         elnodes( 8 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(2,1)
         elnodes( 9 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(3,1)
         elnodes(10 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(4,3)
         elnodes(11 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(4,2)
         elnodes(12 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes(13 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(1,3)
         elnodes(14 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(2,4)
         elnodes(15 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(2,3)
         elnodes(16 + 8 + 12*innerEdgePoints + 5*(innerEdgePoints**2) ) = face_buffer(1,2)

         ! 25th enitty (Face no 6)
         do i=1, innerEdgePoints
            face_buffer(i,:) = elnodes( 1 + (i-1)*innerEdgePoints + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) : &
                                                i*innerEdgePoints + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) )
         end do
         elnodes( 1 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(1,2)
         elnodes( 2 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(2,2)
         elnodes( 3 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(2,1)
         elnodes( 4 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(1,1)
         elnodes( 5 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(2,3)
         elnodes( 6 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(4,2)
         elnodes( 7 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(4,1)
         elnodes( 8 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(3,4)
         elnodes( 9 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(2,4)
         elnodes(10 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(4,3)
         elnodes(11 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(4,4)
         elnodes(12 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(3,3)
         elnodes(13 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(1,3)
         elnodes(14 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(3,1)
         elnodes(15 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(3,2)
         elnodes(16 + 8 + 12*innerEdgePoints + 4*(innerEdgePoints**2) ) = face_buffer(1,4)
      case default
         error stop "Read_GMSH :: Curved elements with P>5 not implemented."
      end select
   end subroutine ReorderElement
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine extract_bc_tags(msh_bcs, tags)
      implicit none
!-----Arguments---------------------------------------------------
      type(MSH_BCinfo_t), intent(in)    :: msh_bcs(:)
      integer, allocatable, intent(out) :: tags(:)
!-----Local-Variables---------------------------------------------
      integer :: i, n
!  -----------------------------------------------------------------------
      
      ! 1. Get the size of the input array
      n = size(msh_bcs)
      
      ! 2. Explicitly allocate the output array
      safedeallocate(tags)
      allocate(tags(n))
      
      ! 3. Manually loop to copy the data directly (No temporaries created!)
      do i = 1, n
         tags(i) = msh_bcs(i) % tag
      end do

   end subroutine extract_bc_tags
end module Read_GMSH