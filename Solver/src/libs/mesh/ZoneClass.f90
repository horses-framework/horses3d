#include "Includes.h"
module ZoneClass
   use SMConstants
   use FaceClass              , only: Face
   use SharedBCModule         , only: zoneNameDictionary
   use MeshTypes              , only: HMESH_INTERIOR, HMESH_MPI
   use BoundaryConditions     , only: ConstructBoundaryConditions
   use IntegerDataLinkedList  , only: IntegerDataLinkedList_t
   use Utilities              , only: toLower
   
   private
   public Zone_t , ConstructZones, ReassignZones, constructZoneModule, AllZoneNames

   integer, parameter      :: STR_LEN_ZONE = BC_STRING_LENGTH
   
   type Zone_t
      integer                     :: marker
      logical                     :: toBeDeleted = .false.
      character(len=STR_LEN_ZONE) :: Name
      integer                     :: no_of_faces
      integer, allocatable        :: faces(:)
      contains
         procedure   :: Initialize       => Zone_Initialize
         procedure   :: copy             => Zone_Assign
		   procedure   :: destruct         => Zone_Destruct
         generic     :: assignment(=)    => copy
         procedure   :: CreateFictitious => Zone_CreateFictitious
   end type Zone_t
   
   contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------
!     Creates the zones and assign faces to them 
!     according to values in the input file 
!     ------------------------------------------
      subroutine ConstructZones( faces , zones )
         implicit none
         type(Face), target                  :: faces(:)
         type(Zone_t), allocatable           :: zones(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                  :: zoneID
         integer                                  :: no_of_markers
         character(len=STR_LEN_ZONE), pointer     :: zoneNames(:)
!
!        Get the number of markers from the Boundary Conditions dictionary
!        -----------------------------------------------------------------         
         no_of_markers = zoneNameDictionary % COUNT() 
         if ( no_of_markers .le. 0 ) return
!
!        Gather the zone names
!        ---------------------
         zoneNames => zoneNameDictionary % allKeys()
!
!        Construct zones
!        ---------------
         allocate ( zones ( no_of_markers ) )
         do zoneID = 1 , no_of_markers
            call toLower(zoneNames(zoneID))
            call zones(zoneID) % Initialize ( zoneID , zoneNames(zoneID) )
         end do
!
!        Assign the faces
!        ----------------         
         call Zone_AssignFaces(faces,zones,no_of_markers,zoneNames)
!
!        Construct Boundary conditions
!        ----------------------------- 
         call ConstructBoundaryConditions(no_of_markers, zoneNames)
         
         deallocate (zoneNames)
		 nullify(zoneNames)
      end subroutine ConstructZones
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------------------------------------------
!     Subroutine for reassigning the zone faces (needed when periodic boundaries are prescribed)
!     ------------------------------------------------------------------------------------------
      subroutine ReassignZones( faces , zones )
         implicit none
         type(Face), target                  :: faces(:)
         type(Zone_t)                        :: zones(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                  :: zoneID, fID
         integer                                  :: no_of_markers
         character(len=STR_LEN_ZONE), pointer     :: zoneNames(:)
!
!        Get the number of markers from the Boundary Conditions dictionary
!        -----------------------------------------------------------------         
         no_of_markers = zoneNameDictionary % COUNT() 
         if ( no_of_markers .le. 0 ) return
!
!        Gather the zone names
!        ---------------------
         zoneNames => zoneNameDictionary % allKeys()
!
!        Reset faces
!        -----------
         do fID = 1, size(faces)
            faces (fID) % zone = 0
         end do
!
!        Assign the faces
!        ----------------         
         call Zone_AssignFaces(faces,zones,no_of_markers,zoneNames)
         
         deallocate (zoneNames)
		 nullify(zoneNames)
      end subroutine ReassignZones
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------
!
!     ------------------------------------------
      subroutine Zone_Initialize ( self , marker , zoneName) 
         implicit none
         class(Zone_t)           :: self
         integer, intent(in)     :: marker
         character(len=*), intent(in)  :: zoneName

         self % marker = marker
         self % Name = trim(zoneName)
         self % no_of_faces = 0

      end subroutine Zone_Initialize
!     ------------------------------------------
!     Destruct Zone
!     ------------------------------------------	  
      subroutine Zone_Destruct ( self) 
         implicit none
         class(Zone_t)           :: self

         safedeallocate(self % faces)

      end subroutine Zone_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------
!     Gather each zone faces
!     ----------------------
!
      subroutine Zone_AssignFaces(faces,zones,no_of_markers,zoneNames)
         implicit none
!        ------------------------------------------------
         integer                    , intent(in)    :: no_of_markers
         type(Face)                 , intent(inout) :: faces(:)
         type(Zone_t)               , intent(inout) :: zones(no_of_markers)
         character(len=STR_LEN_ZONE), intent(in)    :: zoneNames(no_of_markers)
!        ------------------------------------------------
         integer                          :: fID, zoneID
         real(kind=RP), allocatable       :: realArray(:)
         integer, allocatable             :: intArray(:)
         type(IntegerDataLinkedList_t)    :: zoneList(no_of_markers)
!
!        Initialize linked lists
!        -----------------------
         do zoneID=1, no_of_markers
            zoneList(zoneID) = IntegerDataLinkedList_t(.FALSE.)
         end do
!
!        Iterate over all faces and add the face ID to the zone linked list it belongs
!        -----------------------------------------------------------------------------
         do fID = 1, size(faces)
            if (faces(fID) % FaceType == HMESH_INTERIOR .or. faces(fID) % FaceType == HMESH_MPI) cycle
            
            do zoneID = 1, no_of_markers
               if (trim(zoneNames(zoneID)) == trim(faces(fID) % boundaryName)) exit
            end do

            if (zoneID > no_of_markers) cycle
            
            call zoneList(zoneID) % Add(fID)
            faces (fID) % zone = zoneID
         end do
!
!        Dump the linked list contents onto fixed-size arrays 
!        ----------------------------------------------------
         do zoneID = 1 , no_of_markers

            safedeallocate(zones(zoneID) % faces)
            safedeallocate (intArray)
            call zoneList(zoneID) % ExportToArray( intArray  )
            
            zones(zoneID) % no_of_faces = size(intArray)
            
            allocate ( zones(zoneID) % faces(zones(zoneID) % no_of_faces) )
            zones(zoneID) % faces = intArray
            deallocate(intArray)
         end do
         
!        Finish up
!        ---------
         call zoneList % destruct
         
      end subroutine Zone_AssignFaces

      function AllZoneNames(no_of_zones, zones)
         implicit none
         integer,       intent(in)  :: no_of_zones
         type(Zone_t), intent(in)  :: zones(no_of_zones)
         character(len=LINE_LENGTH) :: AllZoneNames(no_of_zones)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: zID

         do zID = 1, no_of_zones
            AllZoneNames(zID) = trim(zones(zID) % Name)
         end do

      end function AllZoneNames
      
      elemental subroutine Zone_Assign (to, from)
         implicit none
         class(Zone_t), intent(inout)  :: to
         type(Zone_t) , intent(in)     :: from
         
         to % marker = from % marker
         to % toBeDeleted = from % toBeDeleted
         to % Name = from % Name
         to % no_of_faces = from % no_of_faces
         
         safedeallocate ( to % faces )
         allocate ( to % faces ( size(from % faces) ) )
         to % faces = from % faces
      end subroutine Zone_Assign

      ! create a fictitious zone, useful to represent fictitious surfaces such as slices or FWH analogy
      Subroutine Zone_CreateFictitious(self, marker, zoneName, no_of_faces, facesID)

         implicit none
         class(Zone_t)                                 :: self
         integer, intent(in)                           :: marker, no_of_faces
         character(len=*), intent(in)                  :: zoneName
         integer, dimension(no_of_faces), intent(in)   :: facesID

         self % marker = marker
         self % Name = trim(zoneName)
         self % no_of_faces = no_of_faces
         allocate(self % faces(no_of_faces))
         self % faces = facesID

      End Subroutine Zone_CreateFictitious

end module ZoneClass