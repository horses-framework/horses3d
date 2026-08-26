#include "Includes.h"
MODULE Read_SpecMesh
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
      use FileReadingUtilities      , only: getFileName
      use Utilities, only: UnusedUnit, toLower
      use SlidingMeshProcedures
      implicit none

      private
      public ConstructMesh_FromSpecMeshFile_, NumOfElems_SpecMesh

!
!     ========
      CONTAINS
!     ========
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------------------------------------
!     Subroutine that constructs mesh from SpecMesh file
!     -> Only valid for conforming meshes
!     -> If the simulation is MPI, it checks if the partitions were already made:
!        * If they were, it constructs the mesh of each partition
!        * If not, it constructs the simplest mesh only with the information needed for partitioning
!     -> If the simulation is sequential, it constructs the full sequential mesh
!     ----------------------------------------------------------------------------------------------
      SUBROUTINE ConstructMesh_FromSpecMeshFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
         USE Physics
         use PartitionedMeshClass
         use MPI_Process_Info
         IMPLICIT NONE
!
!        ---------------
!        Input variables
!        ---------------
!
         type(HexMesh)                    :: self
         integer                          :: nodes
         CHARACTER(LEN=*)                 :: fileName
         integer                          :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
         integer                          :: dir2D
         logical                          :: periodRelative
         LOGICAL           , intent(out)  :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                         :: numberOfElements
         integer                         :: numberOfNodes
         integer                         :: numberOfBoundaryFaces
         integer                         :: numberOfFaces

         integer                         :: bFaceOrder, numBFacePoints
         integer                         :: i, j, k, l
         integer                         :: fUnit, fileStat
         integer                         :: nodeIDs(NODES_PER_ELEMENT), nodeMap(NODES_PER_FACE)
         real(kind=RP)                   :: x(NDIM)
         integer                         :: faceFlags(FACES_PER_ELEMENT)
         CHARACTER(LEN=BC_STRING_LENGTH) :: names(FACES_PER_ELEMENT)
         CHARACTER(LEN=BC_STRING_LENGTH), pointer :: zoneNames(:)
         real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
!
!        ------------------
!        For curved patches
!        ------------------
!
         real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
         real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
!
!        ----------------
!        For flat patches
!        ----------------
!
         real(kind=RP)  , DIMENSION(2)     :: uNodesFlat = [-1.0_RP,1.0_RP]
         real(kind=RP)  , DIMENSION(2)     :: vNodesFlat = [-1.0_RP,1.0_RP]
         real(kind=RP)  , DIMENSION(3,2,2) :: valuesFlat
!
!        ********************************
!        Check if a mesh partition exists
!        ********************************
!
         if ( MPI_Process % doMPIAction ) then
            if ( mpi_partition % Constructed ) then
               call ConstructMeshPartition_FromSpecMeshFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
            else
               call ConstructSimplestMesh_FromSpecMeshFile_ ( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
            end if
            return
         end if

         numberOfBoundaryFaces = 0
         success               = .TRUE.
!
!        -----------------------
!        Read header information
!        -----------------------
!
         fUnit = UnusedUnit()
         OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
         IF ( fileStat /= 0 )     THEN
            PRINT *, "Error opening file: ", fileName
            success = .FALSE.
            RETURN
         END IF

         READ(fUnit,*) numberOfNodes, numberOfElements, bFaceOrder

         self % nodeType = nodes
         self % no_of_elements = numberOfElements
         self % no_of_allElements = numberOfElements
!
!        ----------------------------
!        Set up for face patches
!        Face patches are defined
!        at chebyshev- lobatto points
!        ----------------------------
!
         numBFacePoints = bFaceOrder + 1
         allocate(uNodes(numBFacePoints))
         allocate(vNodes(numBFacePoints))
         allocate(values(3,numBFacePoints,numBFacePoints))

         DO i = 1, numBFacePoints
            uNodes(i) = -COS((i-1.0_RP)*PI/(numBFacePoints-1.0_RP))
            vNodes(i) = -COS((i-1.0_RP)*PI/(numBFacePoints-1.0_RP))
         END DO

!
!        ---------------
!        Allocate memory
!        ---------------
!
         allocate( self % elements(numberOfelements) )

         allocate( self % nodes(numberOfNodes) )


         allocate ( self % Nx(numberOfelements) , self % Ny(numberOfelements) , self % Nz(numberOfelements) )
         self % Nx = Nx
         self % Ny = Ny
         self % Nz = Nz

!
!        ----------------------------------
!        Read nodes: Nodes have the format:
!        x y z
!        ----------------------------------
!
         DO j = 1, numberOfNodes
            READ( fUnit, * ) x
            x = x/Lref
            CALL ConstructNode( self % nodes(j), x, j )
         END DO
!
!        -----------------------------------------
!        Read elements: Elements have the format:
!        node1ID node2ID node3ID node4ID ... node8ID
!        b1 b2 b3 b4 b5 b6
!           (=0 for straight side, 1 for curved)
!        If curved boundaries, then for each:
!           for j = 0 to bFaceOrder
!              x_j  y_j  z_j
!           next j
!        bname1 bname2 bname3 bname4 bname5 bname6
!        -----------------------------------------
!
         DO l = 1, numberOfElements
            READ( fUnit, * ) nodeIDs
            READ( fUnit, * ) faceFlags
!
!           -----------------------------------------------------------------------------
!           If the faceFlags are all zero, then self is a straight-sided
!           hex. In that case, set the corners of the hex8Map and use that in determining
!           the element geometry.
!           -----------------------------------------------------------------------------
!
            IF(MAXVAL(faceFlags) == 0)     THEN
               DO k = 1, NODES_PER_ELEMENT
                  corners(:,k) = self % nodes(nodeIDs(k)) % x
               END DO
               self % elements(l) % SurfInfo % IsHex8 = .TRUE.
               self % elements(l) % SurfInfo % corners = corners

            ELSE
!
!              --------------------------------------------------------------
!              Otherwise, we have to look at each of the faces of the element
!              --------------------------------------------------------------
!
               DO k = 1, FACES_PER_ELEMENT
                 IF ( faceFlags(k) == 0 )     THEN
!
!                    ----------
!                    Flat faces
!                    ----------
!
                     nodeMap           = localFaceNode(:,k)
                     valuesFlat(:,1,1) = self % nodes(nodeIDs(nodeMap(1))) % x
                     valuesFlat(:,2,1) = self % nodes(nodeIDs(nodeMap(2))) % x
                     valuesFlat(:,2,2) = self % nodes(nodeIDs(nodeMap(3))) % x
                     valuesFlat(:,1,2) = self % nodes(nodeIDs(nodeMap(4))) % x

                     call self % elements(l) % SurfInfo % facePatches(k) % construct(uNodesFlat, vNodesFlat, valuesFlat)

                  ELSE
!
!                    -------------
!                    Curved faces
!                    -------------
!
                     DO j = 1, numBFacePoints
                        DO i = 1, numBFacePoints
                           READ( fUnit, * ) values(:,i,j)
                        END DO
                     END DO

                     values = values / Lref

                     call self % elements(l) % SurfInfo % facePatches(k) % construct(uNodes, vNodes, values)

                  END IF
               END DO

            END IF
!
!           -------------------------
!           Now construct the element
!           -------------------------
!
            call self % elements(l) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , l, l)

            READ( fUnit, * ) names
            CALL SetElementBoundaryNames( self % elements(l), names )

            DO k = 1, 6
               IF(TRIM(names(k)) /= emptyBCName) then
                  numberOfBoundaryFaces = numberOfBoundaryFaces + 1
                  zoneNames => zoneNameDictionary % allKeys()
                  if ( all(trim(names(k)) .ne. zoneNames) ) then
                     call zoneNameDictionary % addValueForKey(trim(names(k)), trim(names(k)))
                  end if
                  deallocate (zoneNames)
               end if
            END DO
         END DO      ! l = 1, numberOfElements

!
!        ---------------------------
!        Construct the element faces
!        ---------------------------
!
         numberOfFaces        = (6*numberOfElements + numberOfBoundaryFaces)/2

         self % numberOfFaces = numberOfFaces
         allocate( self % faces(self % numberOfFaces) )
         CALL ConstructFaces( self, success )

!
!        -------------------------
!        Build the different zones
!        -------------------------
!
         call self % ConstructZones()
!
!        ---------------------------
!        Construct periodic faces
!        ---------------------------
!


         CALL ConstructPeriodicFaces( self, periodRelative )
!
!        ---------------------------
!        Delete periodic- faces
!        ---------------------------
!
         CALL DeletePeriodicMinusFaces( self )

!
!        ---------------------------
!        Assign faces ID to elements
!        ---------------------------
!
         CALL getElementsFaceIDs(self)
!
!        ---------------------
!        Define boundary faces
!        ---------------------
!
         call self % DefineAsBoundaryFaces()
!
!        -----------------------------------
!        Check if this is a 2D extruded mesh
!        -----------------------------------
!
         call self % CheckIfMeshIs2D()
!
!        -------------------------------
!        Set the mesh as 2D if requested
!        -------------------------------
!
         self % dir2D_ctrl = dir2D
         if ( dir2D .ne. 0 ) then
            call SetMappingsToCrossProduct
            call self % CorrectOrderFor2DMesh(dir2D,0)
         end if
!
!        ------------------------------
!        Set the element connectivities
!        ------------------------------
!
         call self % SetConnectivitiesAndLinkFaces(nodes)
!
!        ---------------------------------------
!        Construct elements' and faces' geometry
!        ---------------------------------------
!
         call self % ConstructGeometry()
         CLOSE( fUnit )
!
!        -------------------------------
!        Set the mesh as 2D if requested
!        -------------------------------
!
         if ( dir2D .ne. 0 ) then
            call self % CorrectOrderFor2DMesh(dir2D,0)
         end if

!
!        ---------------------------------
!        Describe mesh and prepare for I/O
!        ---------------------------------
!
         if (.not. self % child) CALL self % Describe( trim(fileName), bFaceOrder )
         call self % PrepareForIO

         call self % ExportBoundaryMesh (trim(fileName))

         self % numBFacePoints=numBFacePoints 

      END SUBROUTINE ConstructMesh_FromSpecMeshFile_
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------
!     Constructs the simplest mesh with only the information
!     that is needed to create the MPI partitions
!     ------------------------------------------------------
      SUBROUTINE ConstructSimplestMesh_FromSpecMeshFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
         USE Physics
         use PartitionedMeshClass
         use MPI_Process_Info
         IMPLICIT NONE
         !-arguments------------------------------------------------------------------------------
         type(HexMesh)                    :: self
         integer                          :: nodes
         CHARACTER(LEN=*)                 :: fileName
         integer                          :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
         integer                          :: dir2D
         logical                          :: periodRelative
         LOGICAL           , intent(out)  :: success
         !-local-variables------------------------------------------------------------------------
         integer                         :: numberOfElements
         integer                         :: numberOfNodes
         integer                         :: numberOfBoundaryFaces
         integer                         :: numberOfFaces

         integer                         :: bFaceOrder, numBFacePoints
         integer                         :: i, j, k, l
         integer                         :: fUnit, fileStat
         integer                         :: nodeIDs(NODES_PER_ELEMENT)
         real(kind=RP)                   :: x(NDIM)
         integer                         :: faceFlags(FACES_PER_ELEMENT)
         CHARACTER(LEN=BC_STRING_LENGTH) :: names(FACES_PER_ELEMENT)
         CHARACTER(LEN=BC_STRING_LENGTH), pointer :: zoneNames(:)
         real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)

         !----------------------------------------------------------------------------------------
!
!        ***************
!        Initializations
!        ***************
!
         numberOfBoundaryFaces = 0
         success               = .TRUE.
!
!        -----------------------
!        Read header information
!        -----------------------
!
         fUnit = UnusedUnit()
         OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
         IF ( fileStat /= 0 )     THEN
            PRINT *, "Error opening file: ", fileName
            success = .FALSE.
            RETURN
         END IF

         READ(fUnit,*) numberOfNodes, numberOfElements, bFaceOrder

         self % nodeType = nodes
         self % no_of_elements = numberOfElements

!
!        ----------------------------
!        Set up for face patches
!        Face patches are defined
!        at chebyshev- lobatto points
!        ----------------------------
!
         numBFacePoints = bFaceOrder + 1

!
!        ---------------
!        Allocate memory
!        ---------------
!
         allocate( self % elements(numberOfelements) )
         allocate( self % nodes   (numberOfNodes) )

         allocate ( self % Nx(numberOfelements) , self % Ny(numberOfelements) , self % Nz(numberOfelements) )
         self % Nx = Nx
         self % Ny = Ny
         self % Nz = Nz

!
!        ----------------------------------
!        Read nodes: Nodes have the format:
!        x y z
!        ----------------------------------
!
         DO j = 1, numberOfNodes
            READ( fUnit, * ) x
            x = x / Lref
            CALL ConstructNode( self % nodes(j), x, j )
         END DO
!
!        -----------------------------------------
!        Read elements: Elements have the format:
!        node1ID node2ID node3ID node4ID ... node8ID
!        b1 b2 b3 b4 b5 b6
!           (=0 for straight side, 1 for curved)
!        If curved boundaries, then for each:
!           for j = 0 to bFaceOrder
!              x_j  y_j  z_j
!           next j
!        bname1 bname2 bname3 bname4 bname5 bname6
!        -----------------------------------------
!
         DO l = 1, numberOfElements
            READ( fUnit, * ) nodeIDs
            READ( fUnit, * ) faceFlags
!
!           -----------------------------------------------------------------------------
!           If the faceFlags are all zero, then self is a straight-sided
!           hex. In that case, set the corners of the hex8Map and use that in determining
!           the element geometry.
!           -----------------------------------------------------------------------------
!

            do k = 1, FACES_PER_ELEMENT
               if ( faceFlags(k) /= 0 ) then

!                 Skip points of curved faces
!                 ---------------------------
!
                  do j = 1, numBFacePoints
                     do i = 1, numBFacePoints
                        read( fUnit, * )
                     end do
                  end do

               end if
            end do
!
!           -------------------------
!           Now construct the element
!           -------------------------
!
            call self % elements(l) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , l, l)

            READ( fUnit, * ) names
            CALL SetElementBoundaryNames( self % elements(l), names )

            DO k = 1, 6
               IF(TRIM(names(k)) /= emptyBCName) then
                  numberOfBoundaryFaces = numberOfBoundaryFaces + 1
                  zoneNames => zoneNameDictionary % allKeys()
                  if ( all(trim(names(k)) .ne. zoneNames) ) then
                     call zoneNameDictionary % addValueForKey(trim(names(k)), trim(names(k)))
                  end if
                  deallocate (zoneNames)
               end if
            END DO
         END DO      ! l = 1, numberOfElements

!
!        ---------------------------
!        Construct the element faces
!        ---------------------------
!
         numberOfFaces        = (6*numberOfElements + numberOfBoundaryFaces)/2
         self % numberOfFaces = numberOfFaces
         allocate( self % faces(self % numberOfFaces) )
         CALL ConstructFaces( self, success )
!
!        -------------------------
!        Build the different zones
!        -------------------------
!
         call self % ConstructZones()
!
!        ---------------------------
!        Construct periodic faces
!        ---------------------------
!
         CALL ConstructPeriodicFaces( self, periodRelative )
!
!        ---------------------------
!        Delete periodic- faces
!        ---------------------------
!
         CALL DeletePeriodicMinusFaces( self )
!
!        ---------------------------
!        Assign faces ID to elements
!        ---------------------------
!
         CALL getElementsFaceIDs(self)
!
!        ---------------------
!        Define boundary faces
!        ---------------------
!
         call self % DefineAsBoundaryFaces()
!
!        ------------------------------
!        Set the element connectivities
!        ------------------------------
!
         call self % SetConnectivitiesAndLinkFaces(nodes)

         CLOSE( fUnit )

         call self % ExportBoundaryMesh (trim(fileName))

      END SUBROUTINE ConstructSimplestMesh_FromSpecMeshFile_
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------
!     Constructor of mesh partitions
!     ------------------------------
      SUBROUTINE ConstructMeshPartition_FromSpecMeshFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, periodRelative, success )
         USE Physics
         use PartitionedMeshClass
         use MPI_Process_Info
         use MPI_Face_Class
         IMPLICIT NONE
!
!        ---------------
!        Input variables
!        ---------------
!
         CLASS(HexMesh)      :: self
         integer             :: nodes
         CHARACTER(LEN=*)    :: fileName
         integer             :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
         integer, intent(inout) :: dir2D
         logical             :: periodRelative
         LOGICAL             :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                         :: numberOfAllElements
         integer                         :: numberOfAllNodes
         integer                         :: numberOfFaces
         integer, allocatable            :: globalToLocalNodeID(:)
         integer, allocatable            :: globalToLocalElementID(:)

         integer                         :: bFaceOrder, numBFacePoints
         integer                         :: i, j, k, l, pNode, pElement
         integer                         :: fUnit, fileStat
         integer                         :: nodeIDs(NODES_PER_ELEMENT), nodeMap(NODES_PER_FACE)
         real(kind=RP)                   :: x(NDIM)
         integer                         :: faceFlags(FACES_PER_ELEMENT)
         CHARACTER(LEN=BC_STRING_LENGTH) :: names(FACES_PER_ELEMENT)
         CHARACTER(LEN=BC_STRING_LENGTH), pointer :: zoneNames(:)
         real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
!
!        ------------------
!        For curved patches
!        ------------------
!
         real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
         real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
!
!        ----------------
!        For flat patches
!        ----------------
!
         real(kind=RP)  , DIMENSION(2)     :: uNodesFlat = [-1.0_RP,1.0_RP]
         real(kind=RP)  , DIMENSION(2)     :: vNodesFlat = [-1.0_RP,1.0_RP]
         real(kind=RP)  , DIMENSION(3,2,2) :: valuesFlat
         character(len=LINE_LENGTH)        :: partitionName

         
         success               = .TRUE.
!
!        -----------------------
!        Read header information
!        -----------------------
!
         fUnit = UnusedUnit()
         OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
         IF ( fileStat /= 0 )     THEN
            PRINT *, "Error opening file: ", fileName
            success = .FALSE.
            RETURN
         END IF

         READ(fUnit,*) numberOfAllNodes, numberOfAllElements, bFaceOrder

         self % nodeType = nodes
         self % no_of_elements = mpi_partition % no_of_elements
!
!        ----------------------------
!        Set up for face patches
!        Face patches are defined
!        at chebyshev- lobatto points
!        ----------------------------
!
         numBFacePoints = bFaceOrder + 1
         allocate(uNodes(numBFacePoints))
         allocate(vNodes(numBFacePoints))
         allocate(values(3,numBFacePoints,numBFacePoints))

         DO i = 1, numBFacePoints
            uNodes(i) = -COS((i-1.0_RP)*PI/(numBFacePoints-1.0_RP))
            vNodes(i) = -COS((i-1.0_RP)*PI/(numBFacePoints-1.0_RP))
         END DO
!
!        ---------------
!        Allocate memory
!        ---------------
!
         allocate( self % elements(mpi_partition % no_of_elements) )
         allocate( self % nodes(mpi_partition % no_of_nodes) )
         allocate( globalToLocalNodeID(numberOfAllNodes) )
         allocate( globalToLocalElementID(numberOfAllElements) )
         self % no_of_elements = mpi_partition % no_of_elements
         self % no_of_allElements = numberOfAllElements

         globalToLocalNodeID = -1
         globalToLocalElementID = -1

         allocate ( self % Nx(self % no_of_elements) , self % Ny(self % no_of_elements) , self % Nz(self % no_of_elements) )

!
!        ----------------------------------
!        Read nodes: Nodes have the format:
!        x y z
!        ----------------------------------
!
         pNode = 1
         DO j = 1, numberOfAllNodes
            READ( fUnit, * ) x
            x = x / Lref
            if ( pNode .gt. mpi_partition % no_of_nodes ) cycle
!
!           Construct only nodes that belong to the partition
!           -------------------------------------------------
            if ( j .eq. mpi_partition % nodeIDs(pNode) ) then
               CALL ConstructNode( self % nodes(pNode), x, j )
               globalToLocalNodeID(j) = pNode
               pNode = pNode + 1
            end if
         END DO
!
!        -----------------------------------------
!        Read elements: Elements have the format:
!        node1ID node2ID node3ID node4ID ... node8ID
!        b1 b2 b3 b4 b5 b6
!           (=0 for straight side, 1 for curved)
!        If curved boundaries, then for each:
!           for j = 0 to bFaceOrder
!              x_j  y_j  z_j
!           next j
!        bname1 bname2 bname3 bname4 bname5 bname6
!        -----------------------------------------
!
         pElement = 1
         DO l = 1, numberOfAllElements

            if ( pElement .gt. mpi_partition % no_of_elements ) then
!
!              Skip the element
!              ----------------
               READ( fUnit, * ) nodeIDs
               READ( fUnit, * ) faceFlags

               do k = 1, FACES_PER_ELEMENT
                  if ( faceFlags(k) .ne. 0 ) then
                     DO j = 1, numBFacePoints
                        DO i = 1, numBFacePoints
                           READ( fUnit, * )
                        END DO
                     END DO
                  end if
               end do
!
!              Get the boundary names to build zones (even they could be empty)
!              ----------------------------------------------------------------
               READ( fUnit, * ) names
               DO k = 1, 6
                  IF(TRIM(names(k)) /= emptyBCName) then
                     call toLower( names(k) )
                     zoneNames => zoneNameDictionary % allKeys()
                     if ( all(trim(names(k)) .ne. zoneNames) ) then
                        call zoneNameDictionary % addValueForKey(trim(names(k)), trim(names(k)))
                     end if
                     deallocate (zoneNames)
                  end if
               END DO

               cycle

            else if ( l .ne. mpi_partition % elementIDs(pElement) ) then
!
!              Skip the element
!              ----------------
               READ( fUnit, * ) nodeIDs
               READ( fUnit, * ) faceFlags

               do k = 1, FACES_PER_ELEMENT
                  if ( faceFlags(k) .ne. 0 ) then
                     DO j = 1, numBFacePoints
                        DO i = 1, numBFacePoints
                           READ( fUnit, * )
                        END DO
                     END DO
                  end if
               end do
!
!              Get the boundary names to build zones (even they could be empty)
!              ----------------------------------------------------------------
               READ( fUnit, * ) names
               DO k = 1, 6
                  IF(TRIM(names(k)) /= emptyBCName) then
                     call toLower( names(k) )
                     zoneNames => zoneNameDictionary % allKeys()
                     if ( all(trim(names(k)) .ne. zoneNames) ) then
                        call zoneNameDictionary % addValueForKey(trim(names(k)), trim(names(k)))
                     end if
                     deallocate (zoneNames)
                  end if
               END DO

               cycle
            end if
!
!           Otherwise, read the element
!           ---------------------------
            READ( fUnit, * ) nodeIDs
            READ( fUnit, * ) faceFlags
!
!           Translate nodeIDs to local nodeIDs
!           ----------------------------------
            nodeIDs = globalToLocalNodeID(nodeIDs)
!
!           -----------------------------------------------------------------------------
!           If the faceFlags are all zero, then self is a straight-sided
!           hex. In that case, set the corners of the hex8Map and use that in determining
!           the element geometry.
!           -----------------------------------------------------------------------------
!
            IF(MAXVAL(faceFlags) == 0)     THEN
               DO k = 1, NODES_PER_ELEMENT
                  corners(:,k) = self % nodes(nodeIDs(k)) % x
               END DO
               self % elements(pElement) % SurfInfo % IsHex8 = .TRUE.
               self % elements(pElement) % SurfInfo % corners = corners

            ELSE
!
!              --------------------------------------------------------------
!              Otherwise, we have to look at each of the faces of the element
!              --------------------------------------------------------------
!
               DO k = 1, FACES_PER_ELEMENT
                 IF ( faceFlags(k) == 0 )     THEN
!
!                    ----------
!                    Flat faces
!                    ----------
!
                     nodeMap           = localFaceNode(:,k)
                     valuesFlat(:,1,1) = self % nodes(nodeIDs(nodeMap(1))) % x
                     valuesFlat(:,2,1) = self % nodes(nodeIDs(nodeMap(2))) % x
                     valuesFlat(:,2,2) = self % nodes(nodeIDs(nodeMap(3))) % x
                     valuesFlat(:,1,2) = self % nodes(nodeIDs(nodeMap(4))) % x

                     call self % elements(pElement) % SurfInfo % facePatches(k) % construct(uNodesFlat, vNodesFlat, valuesFlat)

                  ELSE
!
!                    -------------
!                    Curved faces
!                    -------------
!
                     DO j = 1, numBFacePoints
                        DO i = 1, numBFacePoints
                           READ( fUnit, * ) values(:,i,j)
                        END DO
                     END DO

                     values = values / Lref

                     call self % elements(pElement) % SurfInfo % facePatches(k) % construct(uNodes, vNodes, values)

                  END IF
               END DO

            END IF
!
!           -------------------------
!           Now construct the element
!           -------------------------
!
            call self % elements(pElement) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , pElement, l)

            self % Nx(pElement) = Nx(l)
            self % Ny(pElement) = Ny(l)
            self % Nz(pElement) = Nz(l)

            READ( fUnit, * ) names
            CALL SetElementBoundaryNames( self % elements(pElement), names )
!
!           Add the boundary face to the zone Name dictionary
!           -------------------------------------------------
            DO k = 1, 6
               IF(TRIM(names(k)) /= emptyBCName) then
                  zoneNames => zoneNameDictionary % allKeys()
                  if ( all(trim(names(k)) .ne. zoneNames) ) then
                     call zoneNameDictionary % addValueForKey(trim(names(k)), trim(names(k)))
                  end if
                  deallocate (zoneNames)
               end if
            END DO
!
!           Count the partition element
!           ---------------------------
            globalToLocalElementID(l) = pElement
            pElement = pElement + 1

         END DO      ! l = 1, numberOfAllElements
!
!        ---------------------------
!        Construct the element faces
!        ---------------------------
!
         numberOfFaces        = GetOriginalNumberOfFaces(self)
         self % numberOfFaces = numberOfFaces
         allocate( self % faces(self % numberOfFaces) )
         CALL ConstructFaces( self, success )
!
!        --------------------------------
!        Get actual mesh element face IDs
!        --------------------------------
!
         CALL getElementsFaceIDs(self)
!
!        --------------
!        Cast MPI faces
!        --------------
!
         call ConstructMPIFaces( self % MPIfaces )
         call self % UpdateFacesWithPartition(mpi_partition, &
                                              numberOfAllElements, &
                                              globalToLocalElementID)

!
!        -------------------------
!        Build the different zones
!        -------------------------
!
         call self % ConstructZones()
!
!        ---------------------------
!        Construct periodic faces
!        ---------------------------
!
         CALL ConstructPeriodicFaces( self, periodRelative )
!
!        ---------------------------
!        Delete periodic- faces
!        ---------------------------
!
         CALL DeletePeriodicMinusFaces( self )
!
!        ---------------------------
!        Assign faces ID to elements
!        ---------------------------
!
         CALL getElementsFaceIDs(self)
!
!        ---------------------
!        Define boundary faces
!        ---------------------
!

         call self % DefineAsBoundaryFaces()
!
!        -----------------------------------
!        Check if this is a 2D extruded mesh
!        -----------------------------------
!
         call self % CheckIfMeshIs2D()
!
!        -------------------------------
!        Set the mesh as 2D if requested
!        -------------------------------
!
         if ( dir2D .ne. 0 ) then
            call SetMappingsToCrossProduct
            call self % CorrectOrderFor2DMesh(dir2D,1)
         end if
!
!        ------------------------------
!        Set the element connectivities
!        ------------------------------
!
         call self % SetConnectivitiesAndLinkFaces(nodes)
!
!        ---------------------------------------
!        Construct elements' and faces' geometry
!        ---------------------------------------
!
         call self % ConstructGeometry()

         if ( dir2D .ne. 0 ) then
            call self % CorrectOrderFor2DMesh(dir2D,0)
         end if


         CLOSE( fUnit )
!
!        ---------
!        Finish up
!        ---------
!
         if (.not. self % child) then
            CALL self % Describe         ( trim(fileName), bFaceOrder )
            CALL self % DescribePartition( )
         end if
!
!        --------------------
!        Prepare mesh for I/O
!        --------------------
!
         call self % PrepareForIO

         deallocate(globalToLocalNodeID)
         deallocate(globalToLocalElementID)


      END SUBROUTINE ConstructMeshPartition_FromSpecMeshFile_

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      function NumOfElems_SpecMesh( fileName ) result(nelem)
         implicit none
         !----------------------------------
         CHARACTER(LEN=*), intent(in) :: fileName
         integer                      :: nelem
         !----------------------------------
         integer :: k, fUnit
         !----------------------------------

         OPEN(newunit = fUnit, FILE = trim(fileName) )
            READ(fUnit,*) k, nelem, k                    ! Here k is used as default reader since this variables are not important now
         CLOSE(fUnit)

      end function NumOfElems_SpecMesh
end module Read_SpecMesh
