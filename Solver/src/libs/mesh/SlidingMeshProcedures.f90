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
    use FacePatchClass
    use TransfiniteMapClass
    use ElementConnectivityDefinitions
    use ZoneClass                       , only: Zone_t, ConstructZones, ReassignZones
    use MPI_Process_Info
    use MPI_Face_Class
    use FileReadingUtilities            , only: RemovePath, getFileName
    use PartitionedMeshClass            , only: mpi_partition
    use InterpolationMatrices           , only: Tset, TsetM
#ifdef _HAS_MPI_
    use mpi
#endif
   use Headers
    implicit none

    private

    public :: AdvanceSlidingMesh
    public :: IdentifySlidingRegion
    public :: InitializeSlidingConnectivity
    public :: BuildSlidingMortarConnectivity
    public :: RotateSlidingRegion
    public :: ConstructSlidingMortars
    public :: PrintMortarConnectivity
    public :: UpdateSlidingMortarsConnectivity

contains

subroutine AdvanceSlidingMesh(mesh, rotationRadius, rotationCenter, numBFacePoints, nodes, useMPI, angle, rotateEntireMesh)

   use Physics
   use PartitionedMeshClass
   use MPI_Process_Info
   use Headers
   implicit none 

   class(HexMesh), intent(inout)        :: mesh
   real(kind=RP), intent(inout)         :: rotationRadius 
   real(kind=RP), intent(inout)         :: rotationCenter(2)

   !integer, intent(inout)              :: dir2D

   integer, intent(inout)               :: numBFacePoints
   integer                              :: nodes

   logical, intent(in)                  :: useMPI 
   real(kind=RP), intent(in), optional  :: angle
   logical, intent(in), optional        :: rotateEntireMesh


   !========================
   ! Sliding / mesh data
   !========================
   integer :: numSlidingElements
   integer :: numSlidingInterfaceElements 


   !========================
   ! Temporary storage
   !========================
   type(node), allocatable :: tmp_nodes(:)
   integer,    allocatable :: tmp_SlidingMortarConnextion(:,:)

   !========================
   ! Scalars
   !========================
   integer :: inter 
   integer :: originalNodeCount
   integer :: numberOfNodes
   integer :: originalFaceCount

   integer :: l, i, j, eID, fID 
   integer :: sectorID 
   logical :: success
   logical :: isConforming = .FALSE.


   !========================
   ! Geometry 
   !========================
   real(kind=RP) :: offsetParams(4)
   real(kind=RP) :: scaleParams(4)

   real(kind=RP) :: PI 
   real(kind=RP) :: theta

   !========================
   ! Mesh Data
   !========================
   integer :: no_interior_faces, no_boundary_faces, no_mpi_faces
   integer :: no_sequential_elems, no_mpi_elems
   integer :: aux_array(1:3)

   associate(SM => mesh % SlidingMesh)

      if (.not. SM % active) then
         if ( MPI_Process % isRoot ) then
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
         call IdentifySlidingRegion(mesh, rotationRadius, rotationCenter, numSlidingElements, numSlidingInterfaceElements)
      end if
            
      if (.not. SM % active) then
      
         originalNodeCount = size(mesh % nodes)
      
         allocate(tmp_nodes(originalNodeCount))
      
         do i = 1, originalNodeCount
            call ConstructNode(tmp_nodes(i), mesh % nodes(i) % x, mesh % nodes(i) % globID)
         end do 
      
         do i = 1, originalNodeCount
            call mesh % nodes(i) % destruct
         end do 
      
         safedeallocate(mesh % nodes)
      
         numberOfNodes = size(tmp_nodes) + 2 * numSlidingInterfaceElements
         allocate(mesh % nodes(numberOfNodes))
      
         do i = 1, originalNodeCount
            call ConstructNode(mesh % nodes(i), tmp_nodes(i) % x, tmp_nodes(i) % globID)
         end do 
      
         safedeallocate(tmp_nodes)
      
      else 
      
         numberOfNodes = size(mesh % nodes)
         originalNodeCount       = numberOfNodes - 2 * numSlidingInterfaceElements
      
      end if

      sectorID=0

      if (.not. SM % active) then 

         call SM % Initialize(numSlidingInterfaceElements, numSlidingElements, numBFacePoints)

      end if 

      !========================
      ! Optional angle override
      !========================
      if (present(angle)) then 
         SM %omega = SM %omega + angle  
      else 
        SM %omega = SM %omega + SM % theta 
      end if 
      
      !========================
      ! Core mesh modification
      !========================

      call GetSlidingTopologyState(SM % omega, SM % numSlidingInterfaceElements, sectorID, isConforming, SM % localAngle)

      call InitializeSlidingConnectivity(mesh, nodes, SM % numSlidingInterfaceElements, rotationCenter, offsetParams, scaleParams, originalNodeCount, numberOfNodes, theta )

         if (sectorID .NE. SM % currentSectorID) then 

            !SM % sectorID=sectorID
            call UpdateSlidingMortarsConnectivity(mesh, sectorID)
            SM % currentSectorID=sectorID
         end if 


         if (.not. present(rotateEntireMesh)) then 

            do l = 1, size(SM % slidingMortarElems)
         
               do j = 1, 6
         
                  if (mesh % elements(SM % slidingMortarElems(l)) % MortarFaces(j) == 1) then 

                     SM  % mortararr2(l,1) = SM  % slidingMortarConnectivity(l,1)
                     SM  % mortararr2(l,2) = SM  % slidingMortarConnectivity(l,4)
         
                     mesh % elements(SM  % slidingMortarElems(l)) % MortarFaces = 0
                  end if 
         
                  if (mesh % elements(SM  % mortarNeighborElems(l)) % MortarFaces(j) == 1) then 

                     SM  % mortararr1(l,1) = SM  % slidingMortarConnectivity(l,3)
                     SM  % mortararr1(l,2) = SM  % slidingMortarConnectivity(l,6)
         
                     mesh % elements(SM  % mortarNeighborElems(l)) % MortarFaces = 0
                  end if 
         
               end do 
         
            end do 
         
         end if 
         
         originalFaceCount = size(mesh % faces)
         
         !========================
         ! Destroy old faces
         !========================

         if (.not. SM  % active) then 
         
            do i = 1, size(mesh % faces)
               call mesh % faces(i) % Destruct
            end do 
         
         end if 
         
         !========================
         ! Reallocate faces
         !========================

         if (.not. SM  % active) then 
            safedeallocate(mesh % faces)
            allocate(mesh % faces(originalFaceCount + numSlidingInterfaceElements))
         end if 
         
         
         !========================
         ! Rebuild mesh topology
         !========================

         if (.not. SM  % active) call ConstructFaces(mesh, success)
         
         if (.not. SM  % active) then 
         
            if (allocated(mesh % zones)) then 
               do i = 1, size(mesh % zones)
                  call mesh % zones(i) % destruct
               end do
               deallocate(mesh % zones)
            end if 
         
            call mesh % ConstructZones()
            call getElementsFaceIDs(mesh)
            call mesh % DefineAsBoundaryFaces()
         
            if (.not. MPI_Process % doMPIRootAction) then
               call mesh % CheckIfMeshIs2D()
            end if
         
         end if

      !if ( dir2D .ne. 0 ) then
      !   call SetMappingsToCrossProduct
         !  call mesh % CorrectOrderFor2DMesh(dir2D)
      !end if
         if (.not. present(rotateEntireMesh)) then 

            do l = 1, size(SM % slidingMortarElems)
         
               mesh % faces(mesh % elements(SM  % slidingMortarElems(l)) % faceIDs(SM  % mortararr2(l,2))) % MortarType = MORTAR_SLIDING
               mesh % faces(mesh % elements(SM  % slidingMortarElems(l)) % faceIDs(SM  % mortararr2(l,2))) % faceType = 1
         
               mesh % faces(mesh % elements(SM  % mortarNeighborElems(l)) % faceIDs(SM  % mortararr1(l,2))) % MortarType = MORTAR_SLIDING
               mesh % faces(mesh % elements(SM  % mortarNeighborElems(l)) % faceIDs(SM  % mortararr1(l,2))) % faceType = 1
         
            end do 
               
            do l = 1, size(mesh % faces)
               if (mesh % faces(l) % faceType == -1) then 
                  mesh % faces(l) % faceType = 1
               end if 
            end do
         
         end if 
         
         !========================
         ! Connectivity
         !========================

         if (.not. SM % active) then
            call mesh % SetConnectivitiesAndLinkFaces(nodes)
         end if 
         
         !========================
         ! Destroy old geometry
         !========================

         do l = 1, mesh % no_of_elements
            call mesh % elements(l) % geom % destruct
         end do
         
         do l = 1, size(mesh % faces)
            call mesh % faces(l) % geom % destruct
         end do
         
         !========================
         ! Reconstruct geometry
         !========================
         call mesh % ConstructGeometry()
         
         !========================
         ! Mortar geometry cleanup
         !========================
         
         if (allocated(mesh % mortar_faces)) then 
         
            do l = 1, size(mesh % mortar_faces)
               call mesh % mortar_faces(l) % geom % destruct
            end do 
            
         end if

         !
      !===============================================================================
      ! Create three arrays that contain the fIDs of interior, mpi, and boundary faces
      !===============================================================================

         if (.not. SM  % active ) then 
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
            allocate(mesh % faces_mpi     (no_mpi_faces     ))
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

            call ConstructSlidingMortars(mesh, nodes, SM % numSlidingInterfaceElements, &
            SM % mortarNeighborElems, SM % slidingMortarElems, SM % slidingMortarConnectivity, offsetParams, scaleParams, isConforming )

         end if
         
         SM  % active     = .true.
         mesh % slidingflux = .false.

   end associate

end subroutine AdvanceSlidingMesh


subroutine IdentifySlidingRegion (mesh, rotationRadius, rotationCenter, numSlidingElements, numSlidingInterfaceElements)
   IMPLICIT NONE 
   class(HexMesh), intent(inout)  :: mesh
   real(kind=RP), intent(in) :: rotationRadius 
   real(kind=RP), intent(in) :: rotationCenter(2) 
   integer, intent(inout) :: numSlidingElements
   integer, intent(inout) :: numSlidingInterfaceElements

   integer :: elementID, neighborElementID
   integer ::  i, j, numNodesInsideRadius
   integer :: new_nFaces


   numSlidingInterfaceElements=0
   numSlidingElements=0
   new_nFaces=SIZE(mesh % faces)

      ! =========================
      ! Detect elements fully contained
      ! inside the sliding radius
      ! =========================

   do i=1, size(mesh % elements)

      numNodesInsideRadius=0

      do j=1,8
         if (((mesh%Nodes(mesh % elements(i)%nodeIDs(j))%X(1)-rotationCenter(1))**2 +&
         (mesh%Nodes(mesh % elements(i)%nodeIDs(j))%X(3)-rotationCenter(2))**2) .le. rotationRadius**2)  then 

            numNodesInsideRadius = numNodesInsideRadius + 1

         end if 
      end do  

      if (numNodesInsideRadius==8) then

        mesh % elements(i)%sliding=.true.
         numSlidingElements=numSlidingElements+1

      end if 

   end do 

      ! =========================
      ! Detect sliding interface elements
      ! requiring duplicated nodes
      ! =========================

   do i=1, size(mesh % elements)

      if (mesh % elements(i)%sliding) then 

         elementID=mesh % elements(i)%eID
         do j=1,6

            if (mesh % faces(mesh % elements(i)%faceIDs(j))%elementIDs(1)==elementID) then 

                  neighborElementID=mesh % faces(mesh % elements(i)%faceIDs(j))%elementIDs(2)

            else   

                  neighborElementID=mesh % faces(mesh % elements(i)%faceIDs(j))%elementIDs(1)

            end if 

            if (neighborElementID .ne. 0) then 

               if (mesh % elements(neighborElementID)%sliding) then 
                  cycle
               else 

                  mesh % elements(elementID)%sliding_newnodes=.true.
                  numSlidingInterfaceElements=numSlidingInterfaceElements+1
                  new_nFaces=new_nFaces + 1

               end if 

            end if

         end do

      end if

   end do 

end subroutine IdentifySlidingRegion
  
subroutine BuildSlidingMortarConnectivity(mesh, rad, nelm, center, new_nFaces)
   IMPLICIT NONE 
   class(HexMesh), intent(inout)  :: mesh
   real(kind=RP), intent(in) :: rad 
   real(kind=RP), intent(in) :: center(2)
   integer, intent(in)    :: nelm
   integer, intent(inout) :: new_nFaces 

   integer :: elemID, neighborElemID, candidateElemID
   integer ::  i, j, faceIdx2, elemCount, nodeIdxA, nodeIdxB,ind1,ind2,ind3, iNode, iNode2, iNode3
   integer ::  iNodeCheck, isTopRegion, isBottomRegion, hasNodeBelowCenter, isAlreadyConnected, sharedNodeCounter, nNodesInsideRadius

   integer :: coordY 

   coordY      = 3
   nNodesInsideRadius      = 0
   elemCount      = 0
   new_nFaces = size(mesh % faces)


   !========================
   ! Detect sliding elements
   !========================
   do i = 1, size(mesh % elements)

      nNodesInsideRadius = 0

      do j = 1, 8

         if ( ( (mesh % Nodes(mesh % elements(i) % nodeIDs(j)) % X(1) - center(1))**2 + &
                  (mesh % Nodes(mesh % elements(i) % nodeIDs(j)) % X(coordY) - center(2))**2 ) <= rad**2 ) then 
                  nNodesInsideRadius = nNodesInsideRadius + 1
         end if 

      end do  

      if (nNodesInsideRadius == 8) then
        mesh % elements(i) % sliding = .true.
      end if 

   end do 


   !========================
   ! Detect interface (mortar) elements

   !========================
   elemCount = 0

   associate(SM => mesh % SlidingMesh)

      do i = 1, size(mesh % elements)

         if (mesh % elements(i) % sliding) then 

            elemID = mesh % elements(i) % eID

            do j = 1, 6

               if (mesh % faces(mesh % elements(i) % faceIDs(j)) % elementIDs(1) == elemID) then 
                  neighborElemID = mesh % faces(mesh % elements(i) % faceIDs(j)) % elementIDs(2)
               else   
                  neighborElemID = mesh % faces(mesh % elements(i) % faceIDs(j)) % elementIDs(1)
               end if 

               if (neighborElemID /= 0) then 

                  if (mesh % elements(neighborElemID) % sliding) then 
                     cycle
                  else 
                  mesh % elements(elemID) % MortarFaces(j)   = 1
                  mesh % elements(elemID) % sliding_newnodes = .true.

                     elemCount = elemCount + 1
                     SM % slidingMortarElems(elemCount) = elemID

                     new_nFaces = new_nFaces + 1
                  end if 

               end if 

            end do 

         end if 

      end do

      elemCount=0

      SM % slidingMortarConnectivity=0

      do i=1,size(SM % slidingMortarElems)
         SM % slidingMortarConnectivity(i,1)=SM % slidingMortarElems(i)
         isTopRegion=0
         isBottomRegion=0
         hasNodeBelowCenter=0

         ! --------------------------------------------------
         ! Classify element location relative to rotation center
         ! (top or bottom region)
         ! --------------------------------------------------

         do iNode = 1, 8

            if (mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(iNode)) % X(coordY) > 0.0_RP) then 
         
               do iNodeCheck = 1, 8 
                  if (mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(iNodeCheck)) % X(coordY) < 0.0_RP) then 
                     hasNodeBelowCenter = 1
                  end if 
               end do 
         
               if (hasNodeBelowCenter == 0) then 
                  isTopRegion = 1
                  exit 
               end if 
         
            end if 
         
         
            if (mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(iNode)) % X(coordY) < 0.0_RP) then 
               isBottomRegion = 1
               exit 
            end if 
         
         end do
         !it's on the top of the circle 
         if (isTopRegion==1) then 
            do j=1,6
               if (mesh % faces(mesh % elements(SM % slidingMortarElems(i))%faceIDs(j))%elementIDs(1)==SM % slidingMortarElems(i)) then 
                  candidateElemID=mesh % faces(mesh % elements(SM% slidingMortarElems(i))%faceIDs(j))%elementIDs(2)
               else   
                  candidateElemID=mesh % faces(mesh % elements(SM % slidingMortarElems(i))%faceIDs(j))%elementIDs(1)
               end if

               if (candidateElemID==0) cycle 

               if (.not.mesh % elements(candidateElemID)%sliding) then 

                  SM % slidingMortarConnectivity(i,4)=j 
                  SM% slidingMortarConnectivity(i,7)=mesh % faces(mesh % elements(SM % slidingMortarElems(i))%faceIDs(j))%rotation 
                  ! --------------------------------------------------
                  ! Assign local face node ordering for mortar construction
                  ! --------------------------------------------------
                  select case (j)

                  case (1)
                     SM % face_nodes(i,:)      = (/1,2,5,6/)  
                     SM % face_othernodes(i,:) = (/3,4,7,8/)

                  case (2)
                     SM % face_nodes(i,:)      = (/4,3,7,8/)
                     SM % face_othernodes(i,:) = (/1,2,6,5/)

                  case (3)
                    SM % face_nodes(i,:)      = (/1,2,3,4/)
                     SM % face_othernodes(i,:) = (/5,6,7,8/)

                  case (4)
                     SM % face_nodes(i,:)      = (/6,2,3,7/)
                     SM % face_othernodes(i,:) = (/5,1,4,8/)

                  case (5)
                     SM % face_nodes(i,:)      = (/5,6,7,8/)
                     SM % face_othernodes(i,:) = (/1,2,3,4/)

                  case (6)
                     SM % face_nodes(i,:)      = (/5,1,4,8/)
                     SM % face_othernodes(i,:) = (/6,2,3,7/)

                  end select

               end if 

               ! --------------------------------------------------
               ! Search for adjacent sliding element across the interface
               ! using geometric node-based comparison
               ! --------------------------------------------------

               if ((candidateElemID /= 0) .and. (mesh % elements(candidateElemID) % sliding_newnodes)) then
                  if (SM % slidingMortarConnectivity(i,2) == 0) then
                  
                     do nodeIdxA = 1, 8
                  
                        if (SM % slidingMortarConnectivity(i,2) == 0) then
                  
                           do nodeIdxB = 1, 8   
                  
                              if (SM % slidingMortarConnectivity(i,2) == 0) then

                                 if (mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(nodeIdxB)) % X(1) >= center(1)) then 

                                    if ( (mesh % nodes(mesh % elements(candidateElemID) % nodeIDs(nodeIdxA)) % X(1) <  &
                                       mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(nodeIdxB)) % X(1)) .and. &
                                             (mesh % nodes(mesh % elements(candidateElemID) % nodeIDs(nodeIdxB)) % X(coordY) >= center(2)) .and. &
                                             (mesh % nodes(mesh % elements(candidateElemID) % nodeIDs(nodeIdxA)) % X(coordY) >  &
                                             mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(nodeIdxB)) % X(3)) .and. &
                                             (SM % slidingMortarConnectivity(i,3) == 0) ) then 
               
                                          SM % slidingMortarConnectivity(i,2) = candidateElemID 

                                          exit
                  
                                    end if 
                                 end if 
                  
                                 if (mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(nodeIdxB)) % X(1) <= center(1)) then 
                  
                                    if ( (mesh % nodes(mesh % elements(candidateElemID) % nodeIDs(nodeIdxA)) % X(1) <  &
                                    mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(nodeIdxB)) % X(1)) .and. &
                                       (mesh % nodes(mesh % elements(candidateElemID) % nodeIDs(nodeIdxB)) % X(coordY) >= center(2)) .and. &
                                       (mesh % nodes(mesh % elements(candidateElemID) % nodeIDs(nodeIdxA)) % X(coordY) <  &
                                       mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(nodeIdxB)) % X(3)) .and. &
                                       (SM % slidingMortarConnectivity(i,3) == 0) ) then 
                  
                                          SM % slidingMortarConnectivity(i,2) = candidateElemID 

                                       exit
                  
                                    end if 

                                 end if  

                              end if   

                           end do 
                  
                        end if 
                  
                        if (SM % slidingMortarConnectivity(i,2) /= 0) exit
                  
                     end do  

                  end if 
                  
               end if

            end do 
            ! --------------------------------------------------
            ! Determine matching face on mortar element that connects
            ! back to the sliding neighbor (interface consistency)
            ! --------------------------------------------------
            do j = 1, 6

               if (SM% slidingMortarConnectivity(i,2) /= 0) then 
            
                  if (mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2)) % faceIDs(j)) % elementIDs(1) == SM % slidingMortarConnectivity(i,2)) then 
                     candidateElemID = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2)) % faceIDs(j)) % elementIDs(2)
                  else   
                     candidateElemID = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2)) % faceIDs(j)) % elementIDs(1)
                  end if
            
                  if (candidateElemID /= 0) then 
            
                     if (.not. mesh % elements(candidateElemID) % sliding) then 
            
                        SM % slidingMortarConnectivity(i,3) = candidateElemID
                        SM % mortarNeighborElems(i)  = candidateElemID
                        SM % slidingMortarConnectivity(i,5) = j 
            
                        SM % slidingMortarConnectivity(i,8) = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2)) % faceIDs(j)) % rotation
            
                     end if  
            
                  end if    
            
               end if 
            
            end do

            
            do j = 1, 6

               if (SM % slidingMortarConnectivity(i,3) /= 0) then 
            
                  if (mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,3)) % faceIDs(j)) % elementIDs(1) == SM % slidingMortarConnectivity(i,3)) then 
                     candidateElemID = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,3)) % faceIDs(j)) % elementIDs(2)
                  else   
                     candidateElemID = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,3)) % faceIDs(j)) % elementIDs(1)
                  end if
            
                  if (candidateElemID == SM % slidingMortarConnectivity(i,2)) then 
            
                     SM % slidingMortarConnectivity(i,6) = j 
            
                     mesh % elements(SM % slidingMortarConnectivity(i,3)) % MortarFaces(j) = 1 
            
                     SM % slidingMortarConnectivity(i,9) = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,3)) % faceIDs(j)) % rotation
            
                  end if 
            
               end if 
            
            end do

         end if 

         !it's on the bottom of the circle 
         if (isBottomRegion==1) then
         do j=1,6
            if (mesh % faces(mesh % elements(SM % slidingMortarElems(i))%faceIDs(j))%elementIDs(1)==SM % slidingMortarElems(i)) then 

               candidateElemID=mesh % faces(mesh % elements(SM % slidingMortarElems(i))%faceIDs(j))%elementIDs(2)

            else   

               candidateElemID=mesh % faces(mesh % elements(SM % slidingMortarElems(i))%faceIDs(j))%elementIDs(1)

            end if

               if (candidateElemID==0) cycle 

            if (.not.mesh % elements(candidateElemID)%sliding) then 

                  SM % slidingMortarConnectivity(i,4)=j 
                  SM % slidingMortarConnectivity(i,7)=mesh % faces(mesh % elements(SM % slidingMortarElems(i))%faceIDs(j))%rotation

                  ! --------------------------------------------------
                  ! Assign local face node ordering for mortar construction
                  ! --------------------------------------------------

                  select case (j)

                     case (1)
                        SM % face_nodes(i,:)      = (/1,2,5,6/)   
                        SM % face_othernodes(i,:) = (/3,4,7,8/)

                     case (2)
                        SM % face_nodes(i,:)      = (/4,3,7,8/)
                        SM % face_othernodes(i,:) = (/1,2,6,5/)

                     case (3)
                        SM % face_nodes(i,:)      = (/1,2,3,4/)
                        SM % face_othernodes(i,:) = (/5,6,7,8/)

                     case (4)
                        SM % face_nodes(i,:)      = (/6,2,3,7/)
                        SM % face_othernodes(i,:) = (/5,1,4,8/)

                     case (5)
                        SM % face_nodes(i,:)      = (/5,6,7,8/)
                        SM % face_othernodes(i,:) = (/1,2,3,4/)

                     case (6)
                        SM % face_nodes(i,:)      = (/5,1,4,8/)
                        SM % face_othernodes(i,:) = (/6,2,3,7/)

                  end select

               end if 

               if ( (candidateElemID /= 0) .and. (mesh % elements(candidateElemID) % sliding_newnodes) ) then

               if (SM % slidingMortarConnectivity(i,2) == 0) then

                  do nodeIdxA = 1, 8
            
                     if (SM % slidingMortarConnectivity(i,2) /= 0) exit
            
                     do nodeIdxB = 1, 8
            
                        if (SM % slidingMortarConnectivity(i,2) /= 0) exit
                        ! --- Left side of the rotation center ---
                        if ( mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(nodeIdxB)) % X(1) <= center(1) ) then
            
                           if (mesh % nodes(mesh % elements(candidateElemID) % nodeIDs(nodeIdxA)) % X(1) >  &
                              mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(nodeIdxB)) % X(1) .and. &
                              mesh % nodes(mesh % elements(candidateElemID) % nodeIDs(nodeIdxB)) % X(coordY) <= center(2) .and. &
                              mesh % nodes(mesh % elements(candidateElemID) % nodeIDs(nodeIdxA)) % X(coordY) <  &
                              mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(nodeIdxB)) % X(3) ) then
            
                                 SM % slidingMortarConnectivity(i,2) = candidateElemID
                              exit
            
                           end if
                        end if
                        ! --- Right side of the rotation center ---
                        if ( mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(nodeIdxB)) % X(1) >= center(1) ) then
            
                           if ( mesh % nodes(mesh % elements(candidateElemID) % nodeIDs(nodeIdxA)) % X(1) >  &
                              mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(nodeIdxB)) % X(1) .and. &
                              mesh % nodes(mesh % elements(candidateElemID) % nodeIDs(nodeIdxB)) % X(coordY) <= center(2) .and. &
                              mesh % nodes(mesh % elements(candidateElemID) % nodeIDs(nodeIdxA)) % X(coordY) >  &
                              mesh % nodes(mesh % elements(SM % slidingMortarElems(i)) % nodeIDs(nodeIdxB)) % X(3) ) then
            
                              SM % slidingMortarConnectivity(i,2) = candidateElemID
                              exit
            
                           end if

                        end if
            
                     end do

                  end do
            
               end if
            
            end if

         end do 

         ! --------------------------------------------------
         ! Identify non-sliding neighbor (mortar element) 
         ! connected to the sliding neighbor across the interface
         ! --------------------------------------------------
         do j = 1, 6

            if (SM % slidingMortarConnectivity(i,2) /= 0) then 
         
               if (mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2)) % faceIDs(j)) % elementIDs(1) == SM % slidingMortarConnectivity(i,2)) then 
                  candidateElemID = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2)) % faceIDs(j)) % elementIDs(2)
               else   
                  candidateElemID = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2)) % faceIDs(j)) % elementIDs(1)
               end if
         
               if (candidateElemID /= 0) then 
                  if (.not. mesh % elements(candidateElemID) % sliding) then 
         
                     SM % slidingMortarConnectivity(i,3) = candidateElemID
                     SM % mortarNeighborElems(i)  = candidateElemID
                     SM % slidingMortarConnectivity(i,5) = j 
         
                     SM % slidingMortarConnectivity(i,8) = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2)) % faceIDs(j)) % rotation
         
                  end if 
               end if     
         
            end if 
         
         end do

         ! --------------------------------------------------
         ! Locate the corresponding face on the mortar element 
         ! that connects back to the sliding neighbor
         ! (ensures interface consistency and orientation)
         ! --------------------------------------------------

         do j = 1, 6

            if (SM % slidingMortarConnectivity(i,3) /= 0) then 
         
               if (mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,3)) % faceIDs(j)) % elementIDs(1) == SM % slidingMortarConnectivity(i,3)) then 
                  candidateElemID = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,3)) % faceIDs(j)) % elementIDs(2)
               else   
                  candidateElemID = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,3)) % faceIDs(j)) % elementIDs(1)
               end if
         
               if (candidateElemID == SM % slidingMortarConnectivity(i,2)) then 
         
                  SM % slidingMortarConnectivity(i,6) = j 
         
                  mesh % elements(SM % slidingMortarConnectivity(i,3)) % MortarFaces(j) = 1 
         
                  SM % slidingMortarConnectivity(i,9) = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,3)) % faceIDs(j)) % rotation
         
               end if 
         
            end if 
         
         end do

      end if 
      end do 

         ! --------------------------------------------------
         ! Consistency checks:
         ! verify that sliding and mortar connectivity relations
         ! are correctly defined and physically valid
         ! --------------------------------------------------

      do i = 1, size(SM % slidingMortarElems)

         if ((SM % slidingMortarConnectivity(i,2) .NE. 0) .AND. (SM % slidingMortarConnectivity(i,3) .NE. 0)) then

            if (.not. mesh % elements(SM % slidingMortarConnectivity(i,1))%sliding_newnodes) write(*,*) 'for i=', i, '(i,1) is not sliding_newnodes'
            if (.not. mesh % elements(SM % slidingMortarConnectivity(i,2))%sliding_newnodes) write(*,*) 'for i=', i, '(i,2) is not sliding_newnodes'
            if (mesh % elements(SM % slidingMortarConnectivity(i,3))%sliding_newnodes)       write(*,*) 'for i', i, '(i,3) is sliding_newnodes'
            if (mesh % elements(SM % slidingMortarConnectivity(i,3))%sliding)                write(*,*) 'for i', i, '(i,3) is sliding'

            elemID  = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,1))%faceIDs(SM % slidingMortarConnectivity(i,4)))%elementIDs(1)
            neighborElemID = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,1))%faceIDs(SM % slidingMortarConnectivity(i,4)))%elementIDs(2)

            if (elemID == SM % slidingMortarConnectivity(i,1)) then
               if (mesh % elements(neighborElemID)%sliding) write(*,*) 'problem with (i,4)'
            else
               if (mesh % elements(elemID)%sliding) write(*,*) 'problem with (i,4)'
            end if

            if ((elemID .NE. SM % slidingMortarConnectivity(i,1)) .AND. (neighborElemID .NE. SM % slidingMortarConnectivity(i,1))) write(*,*) 'problem with (i,4) line 6195'


            elemID  = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2))%faceIDs(SM % slidingMortarConnectivity(i,5)))%elementIDs(1)
            neighborElemID = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2))%faceIDs(SM % slidingMortarConnectivity(i,5)))%elementIDs(2)

            if (elemID == SM % slidingMortarConnectivity(i,2)) then
               if (mesh % elements(neighborElemID)%sliding) write(*,*) 'problem with (i,5)'
            else
               if (mesh % elements(elemID)%sliding) write(*,*) 'problem with (i,5)'
            end if

            if ((elemID .NE. SM % slidingMortarConnectivity(i,2)) .AND. (neighborElemID .NE. SM % slidingMortarConnectivity(i,2))) write(*,*) 'problem with (i,5) line 6203'


            elemID  = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,3))%faceIDs(SM% slidingMortarConnectivity(i,6)))%elementIDs(1)
            neighborElemID = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,3))%faceIDs(SM % slidingMortarConnectivity(i,6)))%elementIDs(2)

            if (elemID == SM % slidingMortarConnectivity(i,3)) then
               if (.not. mesh % elements(neighborElemID)%sliding) write(*,*) 'problem with (i,6)'
            else
               if (.not. mesh % elements(elemID)%sliding) write(*,*) 'problem with (i,6)'
            end if

            if ((elemID .NE. SM % slidingMortarConnectivity(i,3)) .AND. (neighborElemID .NE. SM % slidingMortarConnectivity(i,3))) write(*,*) 'problem with (i,5) line 6211'

         end if

      end do

      ! --------------------------------------------------
      ! Build element-to-element connectivity graph 
      ! based on shared nodes between sliding interface elements
      ! --------------------------------------------------
      SM % neighborConnectivity = 0


      do i = 1, size(SM % slidingMortarElems)

         do j = 1, 6

            if ( (mesh % faces(mesh%Elements(SM % slidingMortarConnectivity(i,1))%faceIDs(j))%elementIDs(1) .NE. 0) .AND. &
                  (mesh % faces(mesh%Elements(SM % slidingMortarConnectivity(i,1))%faceIDs(j))%elementIDs(2) .NE. 0) ) then

                  elemID = mesh % faces(mesh%Elements(SM % slidingMortarConnectivity(i,1))%faceIDs(j))%elementIDs(1)

               if (.NOT. mesh % elements(elemID)%sliding) then

                  SM % slidingMortarConnectivity(i,10) = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,1))%faceIDs(j))%elementIDs(1)
                  SM % slidingMortarConnectivity(i,12) = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,1))%faceIDs(j))%rotation

                  do faceIdx2 = 1, 6
                     if (mesh % faces(mesh % elements(elemID)%faceIDs(faceIdx2))%elementIDs(1) == SM % slidingMortarConnectivity(i,1)) SM % slidingMortarConnectivity(i,11) = faceIdx2
                     if (mesh % faces(mesh % elements(elemID)%faceIDs(faceIdx2))%elementIDs(2) == SM % slidingMortarConnectivity(i,1)) SM % slidingMortarConnectivity(i,11) = faceIdx2
                  end do

               else

                  elemID = mesh % faces(mesh%Elements(SM % slidingMortarConnectivity(i,1))%faceIDs(j))%elementIDs(2)

                  if (.NOT. mesh % elements(elemID)%sliding) then

                     SM % slidingMortarConnectivity(i,10) = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,1))%faceIDs(j))%elementIDs(2)
                     SM % slidingMortarConnectivity(i,12) = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,1))%faceIDs(j))%rotation

                     do faceIdx2 = 1, 6
                        if (mesh % faces(mesh % elements(elemID)%faceIDs(faceIdx2))%elementIDs(1) == SM % slidingMortarConnectivity(i,1)) SM % slidingMortarConnectivity(i,11) = faceIdx2
                        if (mesh % faces(mesh % elements(elemID)%faceIDs(faceIdx2))%elementIDs(2) == SM % slidingMortarConnectivity(i,1)) SM % slidingMortarConnectivity(i,11) = faceIdx2
                     end do

                  end if

               end if

            end if

         end do

      end do


      if (MPI_Process % doMPIAction) then
         if (.NOT. mpi_partition % Constructed) then

            open(unit=10, file="slidingMortarConnectivity.txt", action="write")

            do i = 1, size(SM % slidingMortarElems)
               write(10,*) (SM % slidingMortarConnectivity(i,j), j=1,12)
            end do

            close(10)

         end if
      end if

      ! --------------------------------------------------
      !  Detect neighboring sliding elements sharing common nodes
      ! --------------------------------------------------

      do i = 1, size(SM % slidingMortarElems)

         nodeIdxA = 2
         SM % neighborConnectivity(i,1,1) = SM % slidingMortarElems(i)
         SM % neighborConnectivity(i,1,2) = i

         do elemCount = 1, size(SM % slidingMortarElems)
            do iNode2 = 1, 8
               do iNode3 = 1, 8

                  if ((mesh % elements(SM % slidingMortarElems(i))%nodeIDs(iNode2) == mesh % elements(SM % slidingMortarElems(elemCount))%nodeIDs(iNode3)) .AND. &
                  SM % slidingMortarElems(i) .ne. SM % slidingMortarElems(elemCount)) then

                     isAlreadyConnected = 0

                     do nodeIdxB = 2, 9
                        if (SM % neighborConnectivity(i,nodeIdxB,1) == SM % slidingMortarElems(elemCount)) then
                           isAlreadyConnected = 1
                        end if
                     end do

                     if (isAlreadyConnected .ne. 1) then
                        SM % neighborConnectivity(i,nodeIdxA,1) = SM % slidingMortarElems(elemCount)
                        SM % neighborConnectivity(i,nodeIdxA,2) = elemCount
                        nodeIdxA= nodeIdxA + 1
                     end if

                  end if

               end do
            end do
         end do

      end do

      ! --------------------------------------------------
      ! Store local indices of shared nodes for each neighbor pair
      ! --------------------------------------------------

      do i = 1, size(SM % slidingMortarElems)

         elemID = SM % neighborConnectivity(i,1,1)

         do iNodeCheck = 2, 9

            neighborElemID = SM % neighborConnectivity(i,iNodeCheck,1)
            sharedNodeCounter    = 3

            if ((elemID .ne. 0) .and. (neighborElemID .ne. 0)) then

               do iNode2 = 1, 8
                  do iNode3 = 1, 8

                     if (mesh % elements(elemID)%nodeIDs(iNode2) == mesh % elements(neighborElemID)%nodeIDs(iNode3)) then
                        SM % neighborConnectivity(i,iNodeCheck,sharedNodeCounter) = iNode2
                        sharedNodeCounter = sharedNodeCounter + 1
                     end if

                  end do
               end do

            end if

         end do

      end do


      if (MPI_Process % doMPIAction) then
         if (.NOT. mpi_partition % Constructed) then

            open(unit=10, file="neighborConnectivity.txt", action="write")

            do i = 1, size(SM % slidingMortarElems)

               write(10,*) (SM % neighborConnectivity(i,:,:))

            end do

            close(10)

         end if
      end if

      ! --------------------------------------------------
      ! Extract sliding elements that do not require mortars
      ! (fully internal sliding region)
      ! --------------------------------------------------
      elemCount = 0

      do i = 1, size(mesh % elements)

         if ((mesh % elements(i)%sliding) .and. .not.(mesh % elements(i)%sliding_newnodes)) then
            elemCount = elemCount + 1
            SM % pureSlidingElems(elemCount) = i
         end if

      end do


      do i = 1, size(SM % slidingMortarElems)

         if ((SM % slidingMortarConnectivity(i,2) .NE. 0) .AND. (SM % slidingMortarConnectivity(i,3) .NE. 0)) then

            if (mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,1))%faceIDs(SM % slidingMortarConnectivity(i,4)))%rotation .NE. SM % slidingMortarConnectivity(i,7)) &
               write(*,*) 'problem with rotation of slidingMortarConnectivity(i7)'

            if ((mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,1))%faceIDs(SM % slidingMortarConnectivity(i,4)))%elementIDs(1) .NE. SM % slidingMortarConnectivity(i,1)) .AND. &
                  (mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,1))%faceIDs(SM % slidingMortarConnectivity(i,4)))%elementIDs(2) .NE. SM % slidingMortarConnectivity(i,1))) &
               write(*,*) 'problem with face i1'


            if (mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2))%faceIDs(SM % slidingMortarConnectivity(i,5)))%rotation .NE. SM %slidingMortarConnectivity(i,8)) &
               write(*,*) 'problem with rotation of slidingMortarConnectivity(i8)'

            if ((mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2))%faceIDs(SM % slidingMortarConnectivity(i,5)))%elementIDs(1) .NE. SM % slidingMortarConnectivity(i,2)) .AND. &
                  (mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2))%faceIDs(SM % slidingMortarConnectivity(i,5)))%elementIDs(2) .NE. SM % slidingMortarConnectivity(i,2))) &
               write(*,*) 'problem with face i2'


            if (mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,3))%faceIDs(SM % slidingMortarConnectivity(i,6)))%rotation .NE. SM % slidingMortarConnectivity(i,9)) &
               write(*,*) 'problem with rotation of slidingMortarConnectivity(i9)'

            if ((mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,3))%faceIDs(SM % slidingMortarConnectivity(i,6)))%elementIDs(1) .NE. SM % slidingMortarConnectivity(i,3)) .AND. &
                  (mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,3))%faceIDs(SM % slidingMortarConnectivity(i,6)))%elementIDs(2) .NE. SM % slidingMortarConnectivity(i,3))) &
               write(*,*) 'problem with face i3'

         end if

      end do
      ! --------------------------------------------------
      ! Ensure uniqueness of connectivity mappings
      ! (no duplicated sliding, neighbor, or mortar assignments)
      ! --------------------------------------------------

      do i = 1, size(SM % slidingMortarElems)

         ind1 = 0
         ind2 = 0
         ind3 = 0

         do j = 1, size(SM % slidingMortarElems)

            if (SM % slidingMortarConnectivity(i,1) == SM % slidingMortarConnectivity(j,1)) ind1 = ind1 + 1
            if (SM % slidingMortarConnectivity(i,2) == SM % slidingMortarConnectivity(j,2)) ind2 = ind2 + 1
            if (SM % slidingMortarConnectivity(i,3) == SM % slidingMortarConnectivity(j,3)) ind3 = ind3 + 1

         end do

         if (ind1 .NE. 1) write(*,*) 'ind1 .NE 1 problem', ind1
         if (ind2 .NE. 1) write(*,*) 'ind2 .NE 1 problem', ind2
         if (ind3 .NE. 1) write(*,*) 'ind3 .NE 1 problem', ind3

      end do

      ! --------------------------------------------------
      ! Final validation of face connectivity and rotation consistency
      ! across sliding-mortar interfaces
      ! --------------------------------------------------
      do i = 1, size(SM % slidingMortarElems)

         if (SM % slidingMortarConnectivity(i,2) .NE. 0) then
            SM % slidingMortarConnectivity(i,8) = mesh % faces(mesh % elements(SM % slidingMortarConnectivity(i,2))%faceIDs(SM % slidingMortarConnectivity(i,5)))%rotation
         end if

      end do

   end associate
end subroutine BuildSlidingMortarConnectivity

subroutine RotateSlidingRegion(mesh, theta, omega, numElementsPerLayer, n, m, newNodeCount, new_nodes, &
                              slidingMortarElems, pureSlidingElems, neighborConnectivity, offsetParams, scaleParams, &
                              mortarFaceNodes, nonMortarFaceNodes, numBFacePoints, originalNodeCount, rotationAxis)

   implicit none

   ! =========================
   ! Arguments
   ! =========================

   class(HexMesh), intent(inout) :: mesh
   real(KIND=RP), intent(in)     :: theta        ! rotation angle
   real(KIND=RP), intent(inout)  :: omega        ! cumulative rotation parameter

   integer, intent(in)           :: numElementsPerLayer      
   integer, intent(in)           :: n, m         ! discretization / rotation counters

   integer, intent(inout)        :: newNodeCount
   type(Node), intent(inout)     :: new_nodes(newNodeCount)

   integer, intent(in) :: slidingMortarElems(numElementsPerLayer) ! sliding elements at interface (need mortars)
   integer, intent(in) :: pureSlidingElems(numElementsPerLayer)   ! sliding elements fully inside region

   integer, intent(in) :: neighborConnectivity(numElementsPerLayer, 9, 6) ! connectivity graph

   real(kind=RP), intent(inout) :: offsetParams(4) ! geometric offsets (mortar mapping)
   real(kind=RP), intent(inout) :: scaleParams(4)  ! geometric scaling (mortar mapping)

   integer, intent(inout) :: mortarFaceNodes(numElementsPerLayer,4)   ! face nodes (local indexing)
   integer, intent(inout) :: nonMortarFaceNodes(numElementsPerLayer,4)   ! complementary face nodes

   integer, intent(inout) :: numBFacePoints
   integer, intent(inout) :: originalNodeCount

   integer, intent(in) :: rotationAxis ! solver convention: 1=X, 2=Z, 3=Y

   ! =========================
   ! Local variables
   ! =========================

   real(KIND=RP) :: ROT(3,3)                     ! rotation matrix
   real(KIND=RP) :: rotatedNodeCoords(8,3)       ! rotated coordinates of element nodes
   real(KIND=RP) :: nodeCoord(3)                 ! temporary node coordinate

   real(KIND=RP) :: rotatedFlatFace(3,2,2)                 ! planar face (2x2)
   real(KIND=RP) :: rotatedFacePatch(3,numBFacePoints,numBFacePoints) ! curved face patch

   integer :: i, j, uIdx, vIdx, localNodeIdx, cornerIdx, neighborIdx
   integer :: newNodeCounter                        ! new node counter
   integer :: eID                                   ! element ID

   integer :: periodicShiftIndex                    ! periodic alignment indicator
   real(kind=RP) :: mortarShiftParam                ! geometric parameter

   integer :: nodeRemapLocal(8)                     ! local remapping of duplicated nodes

   real(KIND=RP) :: corners(3,8)                    ! element corners
   real(KIND=RP) :: faceCorners(3,2,2)              ! flat face
   real(KIND=RP) :: facePatchPoints(3,numBFacePoints,numBFacePoints)

   real(KIND=RP) :: uNodes(numBFacePoints)
   real(KIND=RP) :: vNodes(numBFacePoints)

   ! =========================
   ! Initialization
   ! =========================

   offsetParams = 0.0_RP
   scaleParams  = 0.0_RP

   rotatedNodeCoords = 0.0_RP
   nodeCoord         = 0.0_RP

   nodeRemapLocal = 0

   ! Detect potential periodic alignment after rotation
   periodicShiftIndex = MOD(m, n)

   faceCorners = 0.0_RP

   ! =========================
   ! Parameter computation
   ! =========================

   ! Detect alignment-driven geometric correction
   mortarShiftParam = 1.0_RP - periodicShiftIndex * (2.0_RP / n)

   ! Rotation-dependent scaling (linked to omega evolution)
   mortarShiftParam = (-40.0_RP / (4.0_RP * DATAN(1.0_RP))) * mesh % SlidingMesh % localAngle + 1.0_RP

   ! =========================
   ! Copy existing nodes
   ! =========================

   do i = 1, originalNodeCount
      new_nodes(i)%X      = mesh % nodes(i)%X
      new_nodes(i)%globID = mesh % nodes(i)%globID
   end do

   ! Mark nodes belonging to sliding elements
   do i = 1, size(mesh % elements)
      if (mesh % elements(i)%sliding) then
         do j = 1, 8
            new_nodes(mesh % elements(i)%nodeIDs(j))%tbrotated = .true.
         end do
      end if
   end do

   newNodeCounter = SIZE(mesh % nodes) + 1

   ! =========================
   ! Rotation matrix
   ! =========================

   select case(rotationAxis)

   case(rotationAxis_X)  ! X-axis
      ROT=0.0_RP
      ROT(1,1)=1.0_RP 
      ROT(2,2)=COS(mesh % SlidingMesh % theta)
      ROT(2,3)=-SIN(mesh % SlidingMesh % theta)
      ROT(3,2)=SIN(mesh % SlidingMesh % theta)
      ROT(3,3)=COS(mesh % SlidingMesh % theta)

   case(rotationAxis_Z)  ! Z-axis
      ROT=0.0_RP
      ROT(1,1)=COS(mesh % SlidingMesh % theta)
      ROT(2,2)=COS(mesh % SlidingMesh % theta)
      ROT(3,3)=1.0_RP 
      ROT(1,2)=-SIN(mesh % SlidingMesh % theta)
      ROT(2,1)=SIN(mesh % SlidingMesh % theta)

   case(rotationAxis_Y)  ! Y-axis (solver convention)
      ROT=0.0_RP
      ROT(1,1)=COS(-mesh % SlidingMesh % theta)
      ROT(2,2)=1.0_RP
      ROT(3,3)=COS(-mesh % SlidingMesh % theta) 
      ROT(1,3)=SIN(-mesh % SlidingMesh % theta)
      ROT(3,1)=-SIN(-mesh % SlidingMesh % theta)

   end select

   ! =========================
   ! Mortar geometric parameters
   ! =========================

   offsetParams(1) = (mortarShiftParam - 1.0_RP) / 2.0_RP
   offsetParams(2) = (1.0_RP - mortarShiftParam) / 2.0_RP
   offsetParams(3) = (1.0_RP + mortarShiftParam) / 2.0_RP
   offsetParams(4) = (-mortarShiftParam - 1.0_RP) / 2.0_RP

   scaleParams(1) = offsetParams(1) + 1.0_RP
   scaleParams(2) = 1.0_RP - offsetParams(2)
   scaleParams(3) = 1.0_RP - offsetParams(3)
   scaleParams(4) = offsetParams(4) + 1.0_RP

   ! =========================
   ! Rotate sliding element geometry
   ! - Applies rotation matrix to element corner nodes
   ! - Updates surface representation:
   !     * Linear faces (2x2 nodes)
   !     * Curved face patches (high-order surfaces)
   ! =========================
   do i = 1, mesh%no_of_elements

      if (.not. mesh % elements(i)%sliding) cycle
   
      eID = i
   
      ! =========================
      ! Case 1: Hex8 element (corner-based geometry)
      ! =========================
      if (mesh % elements(eID)%SurfInfo%IsHex8) then
   
         do cornerIdx = 1, 8
            corners(:,cornerIdx) = MATMUL(ROT, mesh % elements(eID)%SurfInfo%corners(:,cornerIdx))
         end do
   
         mesh % elements(eID)%SurfInfo%corners = corners
   
      else
   
         ! =========================
         ! Case 2: High-order faces
         ! =========================
         do j = 1, 6
   
            if (.not. allocated(mesh % elements(eID)%SurfInfo%facePatches(j)%uKnots)) cycle
   
            ! ---------------------------------
            ! Flat face (2x2)
            ! ---------------------------------
            if (mesh % elements(eID)%SurfInfo%facePatches(j)%noOfKnots(1) == 2) then
   
               faceCorners = mesh % elements(eID)%SurfInfo%facePatches(j)%points
   
               rotatedFlatFace(:,1,1) = MATMUL(ROT, faceCorners(:,1,1))
               rotatedFlatFace(:,2,1) = MATMUL(ROT, faceCorners(:,2,1))
               rotatedFlatFace(:,2,2) = MATMUL(ROT, faceCorners(:,2,2))
               rotatedFlatFace(:,1,2) = MATMUL(ROT, faceCorners(:,1,2))
   
               mesh % elements(eID)%SurfInfo%facePatches(j)%points = rotatedFlatFace
   
            else
   
               ! ---------------------------------
               ! Curved face (high-order patch)
               ! ---------------------------------
               uNodes = mesh % elements(eID)%SurfInfo%facePatches(j)%uKnots
               vNodes = mesh % elements(eID)%SurfInfo%facePatches(j)%vKnots
   
               facePatchPoints = mesh % elements(eID)%SurfInfo%facePatches(j)%points
   
               do uIdx = 1, numBFacePoints
                  do vIdx = 1, numBFacePoints
                     rotatedFacePatch(:,vIdx,uIdx) = MATMUL(ROT, facePatchPoints(:,vIdx,uIdx))
                  end do
               end do
   
               if (.not. mesh % SlidingMesh % active) then
                  call mesh % elements(eID)%SurfInfo%facePatches(j)%destruct
               end if
   
               call mesh % elements(eID)%SurfInfo%facePatches(j)%construct(uNodes, vNodes, rotatedFacePatch)
   
            end if
   
         end do ! j
   
      end if
   
   end do ! i

   ! =========================
   ! Rotate interface sliding elements (with mortars)
   ! - Rotate element nodes
   ! - Update interface nodes
   ! - Duplicate nodes to enforce non-conformity
   ! - Maintain connectivity across mortar interface
   ! =========================

   newNodeCounter = originalNodeCount + 1

   do i = 1, size(slidingMortarElems)

      ! Safety check: ensure element is marked as sliding
      if (.not. mesh % elements(slidingMortarElems(i)) % sliding) then
         write(*,*) "problem elem slidingMortarElems"
      end if

      ! -------------------------
      ! Compute rotated coordinates of element nodes
      ! -------------------------
      do j = 1, 8

         if (mesh % elements(slidingMortarElems(i)) % nodeIDs(j) .LE. originalNodeCount) then

            if (mesh % elements(slidingMortarElems(i)) % nodeIDs(j) == 0) then
               write(*,*) 'node ', j, 'of element', i, '=0'
            end if

            nodeCoord = mesh % nodes(mesh % elements(slidingMortarElems(i)) % nodeIDs(j)) % X
            rotatedNodeCoords(j,:) = MATMUL(ROT, nodeCoord)

         else
            rotatedNodeCoords(j,:) = 0.0_RP
         end if

      end do

      ! -------------------------
      ! Update nodes on non-mortar side (no duplication)
      ! -------------------------
      new_nodes(mesh % elements(slidingMortarElems(i)) % nodeIDs(nonMortarFaceNodes(i,1))) % X = rotatedNodeCoords(nonMortarFaceNodes(i,1),:)
      new_nodes(mesh % elements(slidingMortarElems(i)) % nodeIDs(nonMortarFaceNodes(i,2))) % X = rotatedNodeCoords(nonMortarFaceNodes(i,2),:)
      new_nodes(mesh % elements(slidingMortarElems(i)) % nodeIDs(nonMortarFaceNodes(i,3))) % X = rotatedNodeCoords(nonMortarFaceNodes(i,3),:)
      new_nodes(mesh % elements(slidingMortarElems(i)) % nodeIDs(nonMortarFaceNodes(i,4))) % X = rotatedNodeCoords(nonMortarFaceNodes(i,4),:)

      new_nodes(mesh % elements(slidingMortarElems(i)) % nodeIDs(nonMortarFaceNodes(i,1))) % rotated = .true.
      new_nodes(mesh % elements(slidingMortarElems(i)) % nodeIDs(nonMortarFaceNodes(i,2))) % rotated = .true.
      new_nodes(mesh % elements(slidingMortarElems(i)) % nodeIDs(nonMortarFaceNodes(i,3))) % rotated = .true.
      new_nodes(mesh % elements(slidingMortarElems(i)) % nodeIDs(nonMortarFaceNodes(i,4))) % rotated = .true.

      nodeRemapLocal = 0

      ! -------------------------
      ! Handle mortar interface (node duplication + remapping)
      ! -------------------------
      if (.not. mesh % SlidingMesh % active) then

         ! Create duplicated nodes for non-conforming interface
         do j = 1, 4

            if (mesh % elements(slidingMortarElems(i)) % nodeIDs(mortarFaceNodes(i,j)) .LE. originalNodeCount) then

               if (mesh % elements(slidingMortarElems(i)) % nodeIDs(mortarFaceNodes(i,j)) .GT. originalNodeCount) then
                  write(*,*) 'problem logistic 0'
               end if

               new_nodes(newNodeCounter) % globID = newNodeCounter
               new_nodes(newNodeCounter) % X = rotatedNodeCoords(mortarFaceNodes(i,j),:)
               new_nodes(newNodeCounter) % rotated = .true.

               nodeRemapLocal(mortarFaceNodes(i,j)) = newNodeCounter
               newNodeCounter = newNodeCounter + 1

               ! Propagate new node IDs to neighboring elements
               do neighborIdx = 2, 3

                     eID = neighborConnectivity(i,neighborIdx,1)

                     if (eID .ne. 0) then

                        do localNodeIdx = 1, 8

                           if (mesh % elements(slidingMortarElems(i)) % nodeIDs(mortarFaceNodes(i,j)) == &
                           mesh % elements(eID) % nodeIDs(localNodeIdx)) then

                            mesh % elements(eID) % nodeIDs(localNodeIdx) = nodeRemapLocal(mortarFaceNodes(i,j))

                              if (nodeRemapLocal(mortarFaceNodes(i,j)) == 0) then
                                 write(*,*) 'element', eID, 'is receiving 0 for node', localNodeIdx
                              end if

                           end if

                        end do !localNodeIdx

                     end if

               end do !neighborIdx

            end if

         end do !j

         ! Replace original node IDs with duplicated ones
         do j = 1, 4

            if ((mesh % elements(slidingMortarElems(i)) % nodeIDs(mortarFaceNodes(i,j))) .LE. originalNodeCount) then

               if (nodeRemapLocal(mortarFaceNodes(i,j)) == 0) then
                  write(*,*) 'element', slidingMortarElems(i), 'is receiving 0 for node 5', &
                  mesh % elements(slidingMortarElems(i)) % nodeIDs(mortarFaceNodes(i,j))
               end if

               mesh % elements(slidingMortarElems(i)) % nodeIDs(mortarFaceNodes(i,j)) = nodeRemapLocal(mortarFaceNodes(i,j))

            end if

         end do !j

      end if

   end do


   ! =========================
   ! Rotate interior sliding elements (no mortars)
   ! - Simple rotation of all nodes
   ! - No duplication or connectivity updates
   ! =========================

   do i = 1, size(pureSlidingElems)

      ! Safety check
      if (.not. (mesh % elements(pureSlidingElems(i)) % sliding)) then
         write(*,*) "problem elem slidingMortarElems"
      end if

      ! Rotate all nodes
      do j = 1, 8

         nodeCoord = mesh % nodes(mesh % elements(pureSlidingElems(i)) % nodeIDs(j)) % X
         rotatedNodeCoords(j,:) = MATMUL(ROT, nodeCoord)

      end do

      ! Update node coordinates if not already rotated
      do j = 1, 8

         if (.not. new_nodes(mesh % elements(pureSlidingElems(i)) % nodeIDs(j)) % rotated) then

            new_nodes(mesh % elements(pureSlidingElems(i)) % nodeIDs(j)) % X = rotatedNodeCoords(j,:)
            new_nodes(mesh % elements(pureSlidingElems(i)) % nodeIDs(j)) % rotated = .true.

         end if

      end do

   end do

end subroutine RotateSlidingRegion


subroutine InitializeSlidingConnectivity(mesh, nodes, numElementsPerLayer, center, offsetParams, scaleParams, originalNodeCount, totalNodeCount, theta)

   IMPLICIT NONE

   class(HexMesh), intent(inout) :: mesh
   integer, intent(in)           :: nodes
   integer, intent(in)           :: numElementsPerLayer

   real(kind=RP), intent(inout)  :: center(2)
   real(kind=RP), intent(inout)  :: offsetParams(4)
   real(kind=RP), intent(inout)  :: scaleParams(4)
   integer, intent(inout)        :: originalNodeCount

   integer, intent(in)           :: totalNodeCount
   real(kind=RP), intent(inout)  :: theta

   type(Node), allocatable       :: new_nodes(:)
   !type(Node), allocatable     :: tmpNodes(:)

   integer                       :: l, i, j
   integer                       :: new_nNodes, new_nFaces
   integer                       :: n, m

   integer                       :: Connect(numElementsPerLayer, 9, 6)

   integer                       :: sn

   ! =========================
   ! Initialization
   ! =========================

   allocate(new_nodes(totalNodeCount))

   new_nNodes = 0
   new_nFaces = 0

   ! =========================
   ! Rotation setup
   ! - Update cumulative rotation
   ! - Define rotation configuration
   ! =========================

   n  = 3
   m  = 1
   sn = SIZE(mesh % nodes)

   new_nNodes = SIZE(mesh % nodes)

   ! =========================
   ! Initialize mortar mapping parameters
   ! =========================

   offsetParams = 0.0_RP
   scaleParams = 0.0_RP

   ! =========================
   ! Build sliding mortar connectivity
   ! - Identify interface neighbors
   ! - Construct mortar connectivity maps
   ! =========================

   if (.not. mesh % SlidingMesh % active) then

      mesh % SlidingMesh% mortarNeighborElems = 0
      mesh % SlidingMesh% slidingMortarElems = 0

      call BuildSlidingMortarConnectivity(mesh, 1.01_RP, numElementsPerLayer, center, new_nFaces)
   end if

   ! =========================
   ! Rotate sliding mesh region
   ! - Update node coordinates
   ! - Duplicate interface nodes if needed
   ! =========================

   new_nNodes = SIZE(mesh % nodes)

   !if (.not. mesh % SlidingMesh  % active) then

   !call RotateSlidingRegion(theta, mesh % SlidingMesh%omega, numElementsPerLayer, n, m, new_nNodes, new_nodes, &
   !slidingMortarElems, pureSlidingElems, Connect, offsetParams, scaleParams, face_nodes, &
   !                        face_othernodes, numBFacePoints, originalNodeCount, rotationAxis)

   !else

   !call RotateSlidingRegion(theta, mesh % SlidingMesh%omega, numElementsPerLayer, n, m, new_nNodes, new_nodes, &
   !mesh % SlidingMesh%slidingMortarElems, mesh % SlidingMesh%pureSlidingElems, mesh % SlidingMesh%neighborConnectivity, offsetParams, scaleParams, &
   !mesh % SlidingMesh%face_nodes, mesh % SlidingMesh%face_othernodes, numBFacePoints, originalNodeCount, rotationAxis)

   call RotateSlidingRegion(mesh, theta, mesh % SlidingMesh%omega, numElementsPerLayer, n, m, new_nNodes, new_nodes, &
   mesh % SlidingMesh%slidingMortarElems, mesh % SlidingMesh%pureSlidingElems, mesh % SlidingMesh%neighborConnectivity, offsetParams, scaleParams, &
   mesh % SlidingMesh%face_nodes, mesh % SlidingMesh%face_othernodes, mesh % SlidingMesh%numBFacePoints, originalNodeCount, mesh % SlidingMesh % rotationAxis)
   !end if


   ! =========================
   ! Copy updated node coordinates
   ! =========================

   if (.not. mesh % SlidingMesh % active) then
   do i = 1, size(mesh % nodes)
      mesh % nodes(i) % X      = new_nodes(i) % X
      mesh % nodes(i) % GlobID = new_nodes(i) % GlobID
   end do
   end if

   ! =========================
   ! Reset face connectivity data
   ! =========================

   if (.not. mesh % SlidingMesh % active) then

      do l = 1, SIZE(mesh % faces)

         !call mesh % faces(l) % Destruct

         mesh % faces(l) % ID             = -1
         mesh % faces(l) % FaceType       = HMESH_NONE
         mesh % faces(l) % rotation       = 0
         mesh % faces(l) % NelLeft        = -1
         mesh % faces(l) % NelRight       = -1
         mesh % faces(l) % NfLeft         = -1
         mesh % faces(l) % NfRight        = -1
         mesh % faces(l) % Nf             = -1
         mesh % faces(l) % nodeIDs        = -1
         mesh % faces(l) % elementIDs     = -1
         mesh % faces(l) % elementSide    = -1
         mesh % faces(l) % projectionType = -1
         mesh % faces(l) % boundaryName   = ""

      end do

      end if

      ! =========================
      ! Reset mortar face data
      ! =========================

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

   deallocate(new_nodes)

end subroutine InitializeSlidingConnectivity


! -------------------------------------------------------------------------
! Construct sliding mortar interfaces between neighboring mesh elements.
!
! This routine:
!   1. Creates mortar faces for both sliding interface directions
!   2. Links mortar faces with adjacent master/slave elements
!   3. Assigns mortar offsets, scaling parameters and rotations
!   4. Builds mortar geometric mappings and connectivity
!
! Inputs:
!   - slidingMortarConnectivity : connectivity information between
!                                 neighboring sliding elements
!   - offsetParams             : parametric mortar offsets
!   - scaleParams              : parametric mortar scaling factors
!   - confor                   : conforming/non-conforming interface flag
!
! Mortar convention:
!   - MortarPos = 0 : first sliding direction
!   - MortarPos = 1 : opposite sliding direction
!
! -------------------------------------------------------------------------
subroutine ConstructSlidingMortars(mesh, nodes, nelm, mortarNeighborElems, slidingMortarElems, slidingMortarConnectivity, offsetParams, scaleParams, confor)
   implicit none
   class(HexMesh), intent(inout) :: mesh
   integer, intent(in)           :: nodes
   integer, intent(in)           :: nelm
   
   integer, intent(in)           :: mortarNeighborElems(nelm)
   integer, intent(in)           :: slidingMortarElems(nelm)
   integer, intent(in)           :: slidingMortarConnectivity(nelm,12)
   
   real(kind=RP), intent(in)     :: offsetParams(4)
   real(kind=RP), intent(in)     :: scaleParams(4)
   
   logical, intent(in)           :: confor
   
   integer                       :: i, j, mortarIndex
   integer                       :: masterFaceID, slaveFaceID
   integer                       :: masterElementID, slaveElementID
   integer                       :: masterFaceNumber
   
   integer                       :: mortarFaceNodeIDs(4)
   
   integer                       :: elementNodeIDs(8)
   
   integer                       :: NelL(2), NelR(2)

   ! Reset temporary sliding structures
   if (mesh % SlidingMesh % active) then
      call TsetM % destruct
   end if
   
   
   do i = 1, nelm
      do mortarIndex = 1, nelm
         if (mortarNeighborElems(i) == slidingMortarElems(mortarIndex)) then
            write(*,*) ' mortarNeighborElems(i)==slidingMortarElems(mortarIndex)', i, mortarIndex, mortarNeighborElems(i), slidingMortarElems(mortarIndex)
         end if
      end do
   end do
   
   ! Allocate mortar face container
   if (.not. allocated(mesh%mortar_faces)) then
      allocate(mesh % mortar_faces(2*nelm))
   end if
   
   
   mortarIndex = 1

   ! -------------------------------------------------------------------------
   ! Construct first family of mortar interfaces (MortarPos = 0).
   ! -------------------------------------------------------------------------
   
   do i = 1, size(slidingMortarElems)
   
      ! Select master element and interface face
      if (.not. confor) then
       masterElementID      = slidingMortarConnectivity(i,3)
       masterFaceNumber = slidingMortarConnectivity(i,6)
      else
       masterElementID      = slidingMortarConnectivity(i,10)
       masterFaceNumber = slidingMortarConnectivity(i,11)
      end if
   
       ! Extract face node IDs
      elementNodeIDs = mesh % elements(masterElementID) % nodeIDs
   
      do j = 1, 4
       mortarFaceNodeIDs(j) = elementNodeIDs(localFaceNode(j, masterFaceNumber))
      end do
   
   
      masterFaceID = mesh % elements(masterElementID) % faceIDs(masterFaceNumber)  
   
      if (.not. allocated(mesh % faces(masterFaceID) % Mortar)) then
         allocate(mesh % faces(masterFaceID) % Mortar(2))
      end if
      
      ! Create mortar face
      call mesh % mortar_faces(mortarIndex) % Construct( &
           ID        = mortarIndex, &
           nodeIDs   = mortarFaceNodeIDs, &
           elementID = masterElementID, &
           side      = masterFaceNumber)
      
           mesh % mortar_faces(mortarIndex) % Mortarpos = 0
           mesh % mortar_faces(mortarIndex) % FaceType  = HMESH_INTERIOR
      
      
      if (.not. allocated(mesh % mortar_faces(mortarIndex) % Mortar)) then
         allocate(mesh % mortar_faces(mortarIndex) % Mortar(2))
      end if
      
      mesh % mortar_faces(mortarIndex) % Mortar(1) = masterFaceID
      mesh % faces(masterFaceID) % Mortar(1)     = mortarIndex
      
      ! Define slave element connectivity
      slaveElementID = slidingMortarConnectivity(i,2)
      
      mesh % mortar_faces(mortarIndex) % elementIDs(2)  = slaveElementID
      mesh % mortar_faces(mortarIndex) % elementSide(2) = slidingMortarConnectivity(i,5)
      mesh % mortar_faces(mortarIndex) % FaceType       = HMESH_INTERIOR
      
      
      mesh % elements(masterElementID) % faceSide(masterFaceNumber) = 1
      mesh % elements(slaveElementID) % faceSide(slidingMortarConnectivity(i,5))   = 2
      
      
      slaveFaceID = mesh % elements(slaveElementID) % faceIDs(slidingMortarConnectivity(i,5))
      
      mesh % mortar_faces(mortarIndex) % Mortar(2) = slaveFaceID
      mesh % mortar_faces(mortarIndex) % rotation  = slidingMortarConnectivity(i,8)
      
      
      if (.not. allocated(mesh % faces(slaveFaceID) % Mortar)) then
         allocate(mesh % faces(slaveFaceID) % Mortar(2))
      end if
      
      ! Assign mortar parametric transformation
      mesh % mortar_faces(mortarIndex) % offset(1) = offsetParams(1)
      mesh % mortar_faces(mortarIndex) % offset(2) = offsetParams(2)
      
      mesh % mortar_faces(mortarIndex) % s(1) = scaleParams(1)
      mesh % mortar_faces(mortarIndex) % s(2) = scaleParams(2)

      if (confor) then
      
        mesh % mortar_faces(mortarIndex) % offset(1) = 0.0_RP
        mesh % mortar_faces(mortarIndex) % offset(2) = 0.0_RP
      
        mesh % mortar_faces(mortarIndex) % s(1) = 0.0_RP
        mesh % mortar_faces(mortarIndex) % s(2) = 0.0_RP
      
      end if
      
      mesh % faces(slaveFaceID) % Mortar(2) = mortarIndex
   
      mortarIndex = mortarIndex + 1
      
      mesh % elements(slidingMortarElems(i)) % MortarFaces(1) = 3
      mesh % elements(mortarNeighborElems(i)) % MortarFaces(2) = 3
      
   end do
   
   
   ! -------------------------------------------------------------------------
   ! Construct second mortar family (MortarPos = 1)
   ! -------------------------------------------------------------------------
   
   mortarIndex=nelm+1
   
   do i=1,size(slidingMortarElems)
   
      if (.not.confor) then 
   
           masterElementID=slidingMortarConnectivity(i,3)
           masterFaceNumber=slidingMortarConnectivity(i,6)
   
      else 
   
           masterElementID=slidingMortarConnectivity(i,10)
           masterFaceNumber=slidingMortarConnectivity(i,11)
   
      end if
   
      elementNodeIDs = mesh % elements(masterElementID) % nodeIDs
   
      do j = 1, 4
   
           mortarFaceNodeIDs(j) = elementNodeIDs(localFaceNode(j,masterFaceNumber))
   
      end do 
   
      masterFaceID=mesh % elements(masterElementID) % faceIDs(masterFaceNumber)
   
       CALL mesh % mortar_faces(mortarIndex) % Construct(ID  = mortarIndex, &
          nodeIDs = mortarFaceNodeIDs, &
          elementID = masterElementID,       &
          side = masterFaceNumber)

          mesh%mortar_faces(mortarIndex)%Mortarpos=1
          mesh % mortar_faces(mortarIndex) % FaceType       = HMESH_INTERIOR
   
       if (.not.allocated(mesh % mortar_faces(mortarIndex) % Mortar)) allocate(mesh % mortar_faces(mortarIndex) % Mortar(2)) 
   
       mesh % mortar_faces(mortarIndex) % Mortar(1)=masterFaceID
       mesh % faces(masterFaceID) % Mortar(2)=mortarIndex
   
       slaveElementID= slidingMortarConnectivity(i,1)   !3
   
       mesh % mortar_faces(mortarIndex) % elementIDs(2)  = slaveElementID
       mesh % mortar_faces(mortarIndex) % elementSide(2) = slidingMortarConnectivity(i,4)
       mesh % mortar_faces(mortarIndex) % FaceType       = HMESH_INTERIOR
   
       mesh % elements(masterElementID)%faceSide(masterFaceNumber)=1
       mesh % elements(slaveElementID)%faceSide(slidingMortarConnectivity(i,4))=2
   
       slaveFaceID=mesh % elements(slaveElementID) % faceIDs(slidingMortarConnectivity(i,4))
   
       mesh % mortar_faces(mortarIndex) % Mortar(2)=slaveFaceID
       mesh % mortar_faces(mortarIndex) % rotation=slidingMortarConnectivity(i,8)
       if (.not.allocated(mesh % faces(slaveFaceID) % Mortar)) then 
    
            if (.not.allocated(mesh % faces(slaveFaceID) % Mortar)) allocate(mesh % faces(slaveFaceID) % Mortar(2)) 
   
       end if 
   
       mesh % mortar_faces(mortarIndex)%offset(1)=offsetParams(3)
       mesh % mortar_faces(mortarIndex)%offset(2)=offsetParams(4)
   
       mesh % mortar_faces(mortarIndex)%s(1)=scaleParams(3)
       mesh % mortar_faces(mortarIndex)%s(2)= scaleParams(4)
   
       if (confor) then 
         mesh % mortar_faces(mortarIndex)%offset(1)=0.0_RP
         mesh % mortar_faces(mortarIndex)%offset(2)=0.0_RP
   
         mesh % mortar_faces(mortarIndex)%s(1)=1.0_RP
         mesh % mortar_faces(mortarIndex)%s(2)=1.0_RP
   
       end if 
   
       mesh % faces(slaveFaceID) % Mortar(1)=mortarIndex
   
       mortarIndex = mortarIndex + 1    
   
       mesh % elements(slidingMortarElems(i))% MortarFaces(1)=3 
       mesh % elements(mortarNeighborElems(i))% MortarFaces(2)=3                                                        

   end do
   
   
   ! -------------------------------------------------------------------------
   ! Link mortar faces with neighboring elements
   ! -------------------------------------------------------------------------
   
   do i=1, size(mesh%mortar_faces)  
   
       associate(f => mesh % mortar_faces(i))
           associate(eL => mesh % elements(f % elementIDs(1)), &
           eR => mesh % elements(f % elementIDs(2))   )
           NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
           NelR = eR % Nxyz(axisMap(:, f % elementSide(2)))
   
           call f % LinkWithElements(NelL, NelR, nodes, f%offset, f%s)
   
           end associate
       end associate
   end do
   
   
   ! -------------------------------------------------------------------------
   ! Construct mortar geometrical mappings
   ! -------------------------------------------------------------------------
   
   do i=1, size(mesh%mortar_faces)
   
       associate(f => mesh % mortar_faces(i))
           associate(eL => mesh % elements(f % elementIDs(1)), &
           eR => mesh % elements(f % elementIDs(2))   )
           NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
           NelR = eR % Nxyz(axisMap(:, f % elementSide(2)))
   
           call f % geom % construct(f % Nf, f % NelLeft, f % NfLeft, eL % Nxyz, &
                                    eL % geom, eL % hexMap, f % elementSide(1), &
                                    f % projectionType(1), 1, 0,.true.,f%Mortarpos, f%s(1))
   
   
         end associate
         f % geom % h = minval(mesh % elements(f % elementIDs(2)) % geom % jacobian) &
         / maxval(f % geom % jacobian)
   
       end associate
   
   end do 
   !call PrintMortarConnectivity(mesh)
end subroutine ConstructSlidingMortars

subroutine GetSlidingTopologyState(omega, nElements, sectorID, isConforming, localAngle)
   use SMConstants
   use Utilities, only: AlmostEqual
   implicit none

   real(kind=RP), intent(in)  :: omega
   integer,       intent(in)  :: nElements

   integer,       intent(out) :: sectorID
   logical,       intent(out) :: isConforming
   real(kind=RP), intent(out) :: localAngle



   real(kind=RP) :: deltaTheta
   real(kind=RP) :: remainder


   ! Angular spacing between two conforming configurations
   deltaTheta = 2.0_RP * PI / real(nElements, RP)

   ! Determine topology sector
   sectorID = floor( omega / deltaTheta )

   ! Local angle inside current sector
   remainder  = modulo(omega, deltaTheta)
   localAngle = remainder

   ! Check if configuration is conforming
   isConforming = AlmostEqual(remainder, 0.0_RP)
   !if (isConforming) write(*,*) 'its conforming'
   !if (.not.isConforming) write(*,*) 'its not conforming'
end subroutine GetSlidingTopologyState

subroutine UpdateSlidingMortarsConnectivity(mesh, sectorID)
   IMPLICIT NONE
   class(HexMesh), intent(inout) :: mesh 
   integer, intent(in)           :: sectorID 

   integer :: i,j,k
   integer, allocatable :: tmp_SlidingMortarConnetion(:,:)

   associate(SM => mesh % SlidingMesh)
      allocate (tmp_SlidingMortarConnetion(size(SM %slidingMortarElems), 12))

      do k=SM % currentSectorID+1, sectorID
         tmp_SlidingMortarConnetion=SM %slidingMortarConnectivity
         do i=1, size(SM %slidingMortarElems)
                  
            do j=1, size(SM %slidingMortarElems)

               if (tmp_SlidingMortarConnetion(i,2)==SM %slidingMortarConnectivity(j,1)) then 

                  tmp_SlidingMortarConnetion(i,3)=SM %slidingMortarConnectivity(j,3)
                  tmp_SlidingMortarConnetion(i,6)=SM %slidingMortarConnectivity(j,6)
                  tmp_SlidingMortarConnetion(i,9)=SM %slidingMortarConnectivity(j,9)
                  tmp_SlidingMortarConnetion(i,12)=SM %slidingMortarConnectivity(j,12)
                  tmp_SlidingMortarConnetion(i,10)=SM %slidingMortarConnectivity(j,10)

               end if 
            end do 

         end do
         SM %slidingMortarConnectivity=tmp_SlidingMortarConnetion
      end do 
      
      deallocate(tmp_SlidingMortarConnetion)
   end associate

end subroutine UpdateSlidingMortarsConnectivity

subroutine PrintMortarConnectivity(mesh)
   implicit none
   class(HexMesh), intent(inout) :: mesh

   INTEGER :: i,j 

   do i=1 , size(mesh % faces)
     
      if (mesh % faces(i) % MortarType == MORTAR_SLIDING) then 
         write(*,*) 'face',i,'is connected to mortarface; connection:',mesh % faces(i) %  mortar, 'elements:', mesh % faces(i) % elementIDs 
      end if 
   end do 

end subroutine PrintMortarConnectivity 

end module SlidingMeshProcedures