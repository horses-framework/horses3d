#include "Includes.h"
module TessellationTypes

   use SMConstants
   use Utilities
   use IntegerDataLinkedList
   use ParamfileRegions                , only: readValueInRegion

   implicit none

   public :: DescribeSTLPartitions

   integer, parameter :: FORCING_POINT = 1, NOT_FORCING_POINT = 0
   integer, parameter :: POINT_ON_PLANE = 0, POINT_IN_FRONT_PLANE = 1, POINT_BEHIND_PLANE = 2, NumOfVertices = 3

!
!  **************************************************
!  Main type for a list of points
!  **************************************************
   type PointLinkedList
   
      class(point_type), pointer :: head => null()
      integer                    :: NumOfPoints
      
      contains
         procedure :: add        => PointLinkedList_add
         procedure :: remove     => PointLinkedList_Remove
         procedure :: removeLast => PointLinkedList_RemoveLast
         procedure :: destruct   => PointLinkedList_Destruct
   
   end type PointLinkedList

!
!  **************************************************
!  Main type for a generic point point
!  **************************************************
   type point_type
   
      class(point_type), pointer :: next => null(), prev => null()
      
      real(kind=rp), dimension(NDIM) :: coords, ImagePoint_coords, normal, xi, VectorValue
      real(kind=rp)                  :: theta, dist, Rank, ScalarValue
      integer                        :: index, element_index, NumOfIntersections, &
                                        Translate = 0, partition, objIndex, isForcingPoint, &
                                        STLNum, element_in, faceID 
      integer,       dimension(NDIM) :: local_Position
      logical                        :: delete = .false., isInsideBody = .false., &
                                        forcingPoint = .false., isInsideBox = .false.
      real(kind=RP), allocatable     :: invPhi(:,:), b(:), V(:,:,:), bb(:,:), Q(:,:)
      integer,       allocatable     :: nearestPoints(:) 

      contains
         procedure :: copy => point_type_copy
   
   end type   
   
!
!  **************************************************
!  Main type for a list of points
!  **************************************************
   type ObjectLinkedList
   
      class(object_type), pointer :: head => null()
      integer                     :: NumOfObjs
      
      contains
         procedure :: add        => ObjectLinkedList_add
         procedure :: destruct   => ObjectLinkedList_Destruct
   
   end type ObjectLinkedList 
   
!
!  **************************************************
!  Main type for a generic object
!  **************************************************  
   type Object_type

      class(Object_type), pointer :: next => null(), prev => null()
      
      type(point_type), dimension(:),   allocatable :: vertices
      real(kind=rp),    dimension(NDIM)             :: normal, tangent, coords
      integer                                       :: index, NumOfVertices
      integer,          dimension(2)                :: partition

      contains
         procedure :: copy       => object_type_copy
         procedure :: build      => object_type_build
         procedure :: destruct   => object_type_destruct
                  
   end type Object_type
   
!
!  **************************************************
!  Main type for a STL file reader
!  **************************************************
   type STLfile

      type(Object_type), dimension(:), allocatable :: ObjectsList
      integer                                      :: NumOfObjs, partition,      &
                                                      motionAxis, body,          &
                                                      NumOfObjs_OLD
      real(kind=RP)                                :: angularVelocity, ds,       &
                                                      Velocity,                  & 
                                                      rotationMatrix(NDIM,NDIM), &
                                                      rotationCenter(NDIM)
      logical                                      :: move, show, construct = .false. 
      character(len=LINE_LENGTH)                   :: filename, motionType
   
       contains
          procedure :: ReadTessellation
          procedure :: getRotationaMatrix    => STLfile_getRotationaMatrix
          procedure :: getDisplacement       => STLfile_getDisplacement
          procedure :: Clip                  => STL_Clip
          procedure :: updateNormals         => STL_updateNormals
          procedure :: SetIntegration        => STL_SetIntegration
          procedure :: ComputeVectorIntegral => STL_ComputeVectorIntegral
          procedure :: ComputeScalarIntegral => STL_ComputeScalarIntegral
          procedure :: destroy               => STLfile_destroy
          procedure :: Describe              => Describe_STLfile
          procedure :: Copy                  => STLfile_copy
          procedure :: plot                  => STLfile_plot
          procedure :: SetIntegrationPoints  => STL_SetIntegrationPoints

   end type
   
   
   
   type ObjsDataLinkedList_t
      type(ObjData_t), pointer     :: head => null()
      integer                      :: no_of_entries = 0
      contains
         procedure   :: Add      => ObjsDataLinkedList_Add
         procedure   :: check    => CheckObj
         procedure   :: Destruct => ObjsDataLinkedList_Destruct
    end type ObjsDataLinkedList_t
   
    type ObjData_t
      integer                   :: value 
      type(ObjData_t), pointer  :: next 
    end type ObjData_t


    type ObjsRealDataLinkedList_t
      class(ObjRealData_t), pointer :: head => NULL()
      integer                       :: no_of_entries = 0
      contains
         procedure   :: Add      => ObjsRealDataLinkedList_Add
         procedure   :: check    => CheckReal
         procedure   :: Destruct => ObjsRealDataLinkedList_Destruct
    end type ObjsRealDataLinkedList_t
   
   type ObjRealData_t
      real(kind=RP)                 :: value
      class(ObjRealData_t), pointer :: next 
    end type ObjRealData_t



    interface ObjsDataLinkedList_t
      module procedure  ConstructObjsDataLinkedList
    end interface 
    
    interface ObjsRealDataLinkedList_t
      module procedure  ConstructObjsRealDataLinkedList
    end interface 
   
   
   
   
   
   
   
   interface PointLinkedList
      module procedure :: PointLinkedList_Construct
   end interface

   interface ObjectLinkedList
      module procedure :: ObjectLinkedList_Construct
   end interface
   
   character(len=8) :: ROTATION = "rotation"
   character(len=6) :: LINEAR = "linear"

   contains   

      function ConstructObjsDataLinkedList( )
         implicit none
         type(ObjsDataLinkedList_t) :: ConstructObjsDataLinkedList 

         ConstructObjsDataLinkedList% head => null()
         ConstructObjsDataLinkedList% no_of_entries = 0

      end function ConstructObjsDataLinkedList
      
      function ConstructObjsRealDataLinkedList( )
         implicit none
         type(ObjsRealDataLinkedList_t) :: ConstructObjsRealDataLinkedList 

         ConstructObjsRealDataLinkedList% head => null()
         ConstructObjsRealDataLinkedList% no_of_entries = 0

      end function ConstructObjsRealDataLinkedList

      subroutine ObjsDataLinkedList_Add( this, value ) 
         implicit none
         !-arguments----------------------------------------------
         class(ObjsDataLinkedList_t), intent(inout) :: this
         integer,                     intent(in)    :: value
         !-local-variables----------------------------------------
         type(ObjData_t), pointer :: current
         integer                  :: i

         if ( this% no_of_entries .eq. 0 ) then
            allocate( this% head ) 
            this% head% value = value
            this% no_of_entries = 1
         else
            current => this% head    
            do i = 1, this% no_of_entries-1
               current => current% next
            end do
            allocate(current% next)
            current% next% value = value
            this% no_of_entries = this% no_of_entries + 1 
         end if

      end subroutine ObjsDataLinkedList_Add
      
      subroutine ObjsRealDataLinkedList_Add( this, value ) 
         implicit none
         !-arguments---------------------------------------------
         class(ObjsRealDataLinkedList_t), intent(inout) :: this
         real(kind=RP),                   intent(in)    :: value
         !-local-variables---------------------------------------
         type(ObjRealData_t), pointer  :: current
         integer                       :: i

         if ( this% no_of_entries .eq. 0 ) then
            allocate( this% head ) 
            this% head% value = value
            this% no_of_entries = 1
         else
            current => this% head    
            do i = 1, this% no_of_entries-1
               current => current% next
            end do
            allocate(current% next)
            current% next% value = value
            this% no_of_entries = this% no_of_entries + 1 
         end if

      end subroutine ObjsRealDataLinkedList_Add

      logical function CheckObj( this, value ) result( found )
         implicit none
         !-arguments---------------------------------------------
         class(ObjsDataLinkedList_t), intent(inout) :: this
         integer,                     intent(in)    :: value
         !-local-variables---------------------------------------
         type(ObjData_t), pointer :: current
         integer                  :: i
        
         found = .false.
        
         if( this% no_of_entries .eq.0 ) return
        
         current => this% head
        
         do i = 1, this% no_of_entries
            if( current% value .eq. value ) then
               found = .true.
               exit
            end if
            current => current% next
         end do
     
      end function CheckObj

      logical function CheckReal( this, value ) result( found )
         implicit none
         !-arguments----------------------------------------------
         class(ObjsRealDataLinkedList_t), intent(inout) :: this
         real(kind=RP),                   intent(in)    :: value
         !-local-variables----------------------------------------
         type(ObjRealData_t), pointer :: current
         integer                      :: i
        
         found = .false.
        
         if( this% no_of_entries .eq.0 ) return
        
         current => this% head
        
         do i = 1, this% no_of_entries
            if( almostEqual(current% value,value) ) then
               found = .true.
               exit
            end if
            current => current% next
         end do
     
      end function CheckReal

      elemental subroutine ObjsDataLinkedList_Destruct(this)
         implicit none
         !-arguments---------------------------------------------
         class(ObjsDataLinkedList_t), intent(inout) :: this
         !-local-variables---------------------------------------
         type(ObjData_t), pointer :: data, nextdata
         integer                  :: i

         data => this% head
         do i = 1, this% no_of_entries
            nextdata => data% next

            deallocate(data)
            data => nextdata
         end do
         
         this% no_of_entries = 0

      end subroutine ObjsDataLinkedList_Destruct
      
      elemental subroutine ObjsRealDataLinkedList_Destruct(this)
         implicit none
         !-arguments---------------------------------------------
         class(ObjsRealDataLinkedList_t), intent(inout) :: this
         !-local-variables---------------------------------------
         type(ObjRealData_t), pointer :: data, nextdata
         integer                      :: i

         data => this% head
         do i = 1, this% no_of_entries
            nextdata => data% next

            deallocate(data)
            data => nextdata
         end do
         
         this% no_of_entries = 0

      end subroutine ObjsRealDataLinkedList_Destruct

!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine initializes a points list
!  ------------------------------------------------ 
   function PointLinkedList_Construct(  )
   
      implicit none
      
      type(PointLinkedList) :: PointLinkedList_Construct
      
      PointLinkedList_Construct% head        => null()
      PointLinkedList_Construct% NumOfPoints =  0
   
   end function PointLinkedList_Construct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine adds a point to a points list
!  ------------------------------------------------ 
   subroutine PointLinkedList_Add( this, point )
      implicit none
      !-arguments-----------------------------------
      class(PointLinkedList) :: this
      type(point_type)       :: point
      !-local-variables-----------------------------
      type(point_type), pointer :: current,    &
                                   currentNext 
   
      if( .not. associated(this% head) ) then 
         allocate(this% head)
         call this% head% copy(point)
         this% head% next => this% head
         this% head% prev => this% head
      else
         current => this% head% prev
         allocate(currentNext) 
         call currentNext% copy(point)
         currentNext% next => this% head
         currentNext% prev => current
         nullify(current% next)
         current% next => currentNext
         nullify(this% head% prev)
         this% head% prev => currentNext
      end if
      
      nullify(current, currentNext)

      this% NumOfPoints = this% NumOfPoints + 1
   
   end subroutine PointLinkedList_Add
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine deletes the point p from a points list
!  ------------------------------------------------  
   subroutine PointLinkedList_Remove( this, p )
      implicit none
      !-arguments--------------------------------------------------------------
      class(PointLinkedList),   intent(inout) :: this
      type(point_type), target, intent(inout) :: p
      !-local-variables--------------------------------------------------------
      type(point_type), pointer :: dataNext => null(), dataPrev => null(), &
                                   dataDel => null()
      
      dataDel => p
      
      dataPrev => p% prev; dataNext => p% next
      
      if( associated(dataDel, this% head) ) then
         error stop ":: Head of list can't be deleted here!"
         return
      end if
      
      dataPrev% next => null()
      dataPrev% next => dataNext
      dataNext% prev => null()
      dataNext% prev => dataPrev
   
      nullify(dataPrev, dataNext)

      this% NumOfPoints = this% NumOfPoints - 1
   
   end subroutine PointLinkedList_Remove
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine deletes the last point from a points list
!  ------------------------------------------------  
   subroutine PointLinkedList_RemoveLast( this )
      implicit none
      !-arguments--------------------------------------------------------
      class(PointLinkedList), intent(inout) :: this
      !-local-variables--------------------------------------------------
      type(point_type), pointer :: data => null(), dataPrev => null()
      
      data => this% head% prev
   
      if( associated(data, this% head) ) then
         error stop ":: Head of list can't be deleted here!"
         return
      end if
      
      dataPrev => data% prev
      
      deallocate(data)
      
      dataPrev% next => null()
      dataPrev% next => this% head
      this% head% prev => null() 
      this% head% prev => dataPrev 
   
      nullify(dataPrev)

      this% NumOfPoints = this% NumOfPoints - 1
   
   end subroutine PointLinkedList_RemoveLast
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine deletes a point from a points list
!  ------------------------------------------------ 
   subroutine PointLinkedList_Destruct( this )
      implicit none
      !-arguments--------------------------------------
      class(PointLinkedList), intent(inout) :: this
      !-local-variables--------------------------------
      class(point_type), pointer :: current, next 
      integer                    :: i
      
      if( this% NumOfPoints .eq. 0 ) then
         nullify(this% head)
         return
      end if
      
      current => this% head

      do i = 1, this% NumOfPoints - 1
         next => current% next
         deallocate(current)
         current => next
      end do
      deallocate(current)
      
      this% NumOfPoints = 0
      nullify(this% head)
      
   end subroutine PointLinkedList_Destruct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine copies a point
!  ------------------------------------------------  

   subroutine point_type_copy( this, point )
      implicit none
      !-arguments-------------------------------
      class(point_type), intent(inout) :: this
      type(point_type),  intent(in)    :: point
    
      this% coords            = point% coords
      this% ImagePoint_coords = point% ImagePoint_coords
      this% theta             = point% theta
      this% index             = point% index
      this% element_index     = point% element_index
      this% local_position    = point% local_position
      this% Translate         = point% Translate
      this% partition         = point% partition
   
   end subroutine point_type_copy
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine initializes an objects list
!  ------------------------------------------------  
   function ObjectLinkedList_Construct(  )
      implicit none
      !-local-variables-------------------------------------
      type(ObjectLinkedList) :: ObjectLinkedList_Construct
      
      ObjectLinkedList_Construct% head => null()
      ObjectLinkedList_Construct% NumOfObjs = 0
   
   end function ObjectLinkedList_Construct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine adds an object to an objects list
!  ------------------------------------------------  

   subroutine ObjectLinkedList_Add( this, object )
      implicit none
      !arguemnts--------------------------------------------
      class(ObjectLinkedList) :: this
      type(Object_type)       :: object
      !-local-variables-------------------------------------
      type(object_type), pointer :: current => null(), &
                                    currentNext => null()
   
      if( .not. associated(this% head) ) then 
         allocate(this% head)
         call this% head% copy(object)
         this% head% next => this% head
         this% head% prev => this% head
      else
         current => this% head% prev
         allocate(currentNext) 
         call currentNext% copy(object)
         currentNext% next => this% head
         currentNext% prev => current
         nullify(current% next)
         current% next => currentNext
         nullify(this% head% prev)
         this% head% prev => currentNext
         nullify(current, currentNext)
      end if
      
      this% NumOfObjs = this% NumOfObjs + 1
   
   end subroutine ObjectLinkedList_Add
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine destroys an objects list
!  ------------------------------------------------ 
   subroutine ObjectLinkedList_Destruct( this )
      implicit none
      !arguemnts----------------------------------------
      class(ObjectLinkedList), intent(inout) :: this
      !-local-variables---------------------------------
      type(object_type), pointer :: data => null(), &
                                    dataPrev => null()
      
      if( .not. associated(this% head) ) return
      
      data => this% head% prev
      
      if( .not. associated(data, this% head) ) then
         do
            dataPrev => data% prev   
            call data% destruct()
            deallocate(data)
            data => dataPrev
            if( associated(data, this% head) ) exit
         end do
      end if
      
      deallocate(data)
      
      this% head      => null()
      this% NumOfObjs =  0
      
      nullify(dataPrev)
      
   end subroutine ObjectLinkedList_Destruct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine builds an object
!  -----------------------------------------------
   subroutine object_type_build( this, Points, normal, NumOfVertices, index )
      implicit none
      !-arguments---------------------------
      class(Object_type), intent(inout) :: this
      real(kind=RP),      intent(in)    :: Points(:,:)
      real(kind=RP),      intent(in)    :: normal(:)
      integer,            intent(in)    :: NumOfVertices, index
      !-local-variables--------------------------------------
      integer :: i
      
      if( allocated(this% vertices) ) deallocate(this% vertices)
      allocate(this% vertices(NumOfVertices))
      
      do i = 1, NumOfVertices
         this% vertices(i)% coords = Points(:,i)
      end do
      
      this% index            = index 
      this% normal           = normal 
      this% NumOfVertices    = NumOfVertices
      
   end subroutine object_type_build
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine copies an object
!  -----------------------------------------------
   subroutine object_type_copy( this, Object )
    
      implicit none
      !-arguments---------------------------
      class(Object_type), intent(inout) :: this
      type(Object_type),  intent(in)    :: Object
      !-local-variables-----------------------
      integer :: i
      
      allocate(this% vertices(Object% NumOfVertices))
      
      do i = 1, Object% NumOfVertices
         call this% vertices(i)% copy(Object% vertices(i))
      end do
      
      this% index            = Object% index 
      this% normal           = Object% normal 
      this% NumOfVertices    = Object% NumOfVertices
      this% partition        = Object% partition
      
   end subroutine object_type_copy
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine copies an object
!  -----------------------------------------------
   subroutine object_type_destruct( this )
    
      implicit none
      !-arguments---------------------------
      class(Object_type), intent(inout) :: this
      
      deallocate(this% vertices)
      
      this% index  = 0
      this% normal = 0.0_RP 
      this% NumOfVertices = 0
      
   end subroutine object_type_destruct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine adds an object to a list of objects
!  -----------------------------------------------
   subroutine addObj( ObjList, Object )
   
      implicit none
      !-arguments----------------------------------------------------
      type(Object_type), target, intent(inout) :: ObjList
      type(Object_type),         intent(in)    :: Object
      !-local-variables----------------------------------------------
      type(Object_type), pointer :: Obj=>null(), Objprev => null() 
     
      if( .not. associated(ObjList% next)) then
         call ObjList% copy( Object )
         ObjList% next => ObjList
         ObjList% prev => ObjList
         return
      end if
      
      Obj => ObjList% prev
!            
!     previous point
!     ------------
      Objprev => ObjList% prev
      
      Obj => Obj% next
      allocate(Obj)
      
      Obj% next => ObjList
      Obj% prev => Objprev
      
      Objprev% next => Obj
      
      call Obj% copy( Object )
      
      ObjList% prev => null()
      ObjList% prev => Obj
      
      if( associated(Obj% prev, ObjList) ) then
         ObjList% next => Obj 
      end if
      
      nullify(Obj, Objprev)
      
   end subroutine addObj
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine describes the .stl file
!  -----------------------------------------------
   subroutine Describe_STLfile( this, filename )
      use Headers
      use MPI_Process_Info
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------
      class(STLfile),   intent(inout) :: this
      character(len=*), intent(in)    :: filename

      if( MPI_Process% isRoot) then
       
         write(STD_OUT,'(/)')
         call Section_Header("Reading stl file")
         write(STD_OUT,'(/)')
         call SubSection_Header('Stl file "' // trim(fileName) // '"')
      
         write(STD_OUT,'(30X,A,A35,I10)') "->" , "Number of objects: " , this% NumOfObjs
         if( this% move ) then
            write(STD_OUT,'(30X,A,A35,A10)') "->" , "Motion: " , this% motionType
            if( this% motionType .eq. ROTATION ) then 
               write(STD_OUT,'(30X,A,A35,F10.3,A)') "->" , "Angular Velocity: " , this% angularVelocity, " rad/s."
               write(STD_OUT,'(30X,A,A35,F10.3,A,F10.3,A,F10.3,A)') "->" , "Center of rotation: [" , this% rotationCenter(1), ',', & 
                                                                                                     this% rotationCenter(2), ',', & 
                                                                                                     this% rotationCenter(3), ']'
            elseif( this% motionType .eq. LINEAR ) then
               write(STD_OUT,'(30X,A,A35,F10.3,A)') "->" , "Translation Velocity: " , this% Velocity, " m/s."
            end if
            write(STD_OUT,'(30X,A,A35,I10)') "->" , "Axis of motion: " , this% motionAxis
         end if
         
      end if
   
   end subroutine Describe_STLfile
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine describes the .stl file
!  -----------------------------------------------
   subroutine DescribeSTLPartitions( partition, NumOfObjs )
      use Headers
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------------
      integer, intent(in) :: partition, NumOfObjs
      !-local-variables--------------------------
      character(len=LINE_LENGTH) :: myString
      
      if( .not. MPI_Process% doMPIAction ) return

      write(myString,'(i100)') partition+1
      
      if( partition .eq. 0 ) then
         write(STD_OUT,'(/)')
         call Section_Header("stl partitions")
         write(STD_OUT,'(/)')
      end if
      
      call SubSection_Header('partition ' // trim(adjustl(myString)))
      write(STD_OUT,'(30X,A,A32,I10)') "->" , "Number of objects: " , NumOfObjs
      
   end subroutine DescribeSTLPartitions
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine reads the .stl file
!  -----------------------------------------------
   subroutine  ReadTessellation( this, filename )
      use PhysicsStorage
      implicit none
      !-arguments---------------------------------------
      class(STLfile),             intent(inout) :: this
      character(len=*),           intent(in)    :: filename
      !-local-variables---------------------------------
      integer              :: i, j, funit, NumOfTri,   &
                              fileStat, NumOfVertices
      integer*2            :: padding
      real*4, dimension(3) :: norm, vertex
      character(len=80)    :: header

      NumOfVertices = 3 
      
      this% partition = 1
      
      funit = UnusedUnit()
      open(unit=funit,file='MESH/'//trim(filename)//'.stl',status='old',access='stream',form='unformatted', iostat=fileStat)
      
      if( fileStat .ne. 0 ) then
         print *, "Read Tessellation: file '",trim(filename),"' not found"
         error stop
      end if
      
      this% filename = filename

      read(funit) header
      read(funit) NumOfTri
      
      this% NumOfObjs = NumOfTri

      allocate(this% ObjectsList(NumOfTri))

      associate( Objs => this% ObjectsList )
      
      do i = 1, NumOfTri
         read(funit) norm(1), norm(2), norm(3)
         Objs(i)% normal = norm
         allocate(Objs(i)% vertices(NumOfVertices))
         do j = 1, NumOfVertices
            read(funit) vertex(1), vertex(2), vertex(3)
            Objs(i)% vertices(j)% coords = vertex  !/Lref -> always 1           
         end do
         read(funit) padding
         Objs(i)% index = i
         Objs(i)% NumOfVertices = NumOfVertices
         Objs(i)% partition = 1
      end do
        
      end associate   
       
      close(unit=funit)

   end subroutine  ReadTessellation

   subroutine STL_SetIntegrationPoints( this )

      implicit none

      class(STLfile), intent(inout) :: this 

      real(kind=RP) :: tmp(NDIM,NumOfVertices)
      integer       :: i, j, indeces(3)

      do i = 1, this% NumOfObjs
         associate( obj => this% ObjectsList(i) )
         do j = 1, NumOfVertices
            tmp(:,j) = obj% vertices(j)% coords
         end do 
         deallocate(obj% vertices)
         allocate(obj% vertices(NumOfVertices+4))
         obj% vertices(NumOfVertices+4)% coords = 0.0_RP
         indeces =(/ 2, 3, 1 /)
         do j = 1, NumOfVertices
            obj% vertices(j)% coords = tmp(:,j)
            obj% vertices(NumOfVertices+j)% coords = 0.5_RP*(tmp(:,j) + tmp(:,indeces(j)))
            obj% vertices(NumOfVertices+4)% coords = obj% vertices(NumOfVertices+4)% coords + &
                                                     tmp(:,j)
         end do
         obj% vertices(NumOfVertices+4)% coords = obj% vertices(NumOfVertices+4)% coords/3.0_RP
         end associate
      end do

   end subroutine STL_SetIntegrationPoints  


   subroutine STLfile_plot( this, iter )
      use MPI_Process_Info
      use PhysicsStorage
      implicit none
      !-arguments----------------------------------------------
      class(STLfile), intent(inout) :: this
      integer,        intent(in)    :: iter
      !-local-variables----------------------------------------
      character(len=LINE_LENGTH)     :: filename
      integer                        :: i, j, funit
      integer*2                      :: padding = 0
      real*4, dimension(3)           :: norm, vertex
      character(len=80)              :: header = repeat(' ',80)
            
      funit = UnusedUnit()
      
      write(filename,'(A,A,I10.10)') trim(this% filename),'_', iter
      
      open(funit,file='MESH/'//trim(filename)//'.stl', status='unknown',access='stream',form='unformatted')
 
      write(funit) header, this% NumOfObjs 
      
      do i = 1, this% NumOfObjs
         norm = this% ObjectsList(i)% normal
         write(funit) norm(1), norm(2), norm(3)
         do j = 1, NumOfVertices      
            vertex = this% ObjectsList(i)% vertices(j)% coords * Lref
            write(funit) vertex(1), vertex(2), vertex(3)
         end do
         write(funit) padding      
      end do

      close(funit)
   
   end subroutine STLfile_plot

   subroutine STLfile_copy( this, STL )

      implicit none 

      class(STLfile), intent(inout) :: this
      type(STLfile),  intent(in)    :: STL

      integer :: NumOfObjs, i, j

      NumOfObjs = STL% NumOfObjs

      this% NumOfObjs = NumOfObjs

      allocate( this% ObjectsList(NumOfObjs) )

!$omp parallel
!$omp do schedule(runtime) private(j)
      do i = 1, NumOfObjs
         allocate(this% ObjectsList(i)% vertices(NumOfVertices))
         do j = 1, NumOfVertices
            this% ObjectsList(i)% vertices(j)% coords = STL% ObjectsList(i)% vertices(j)% coords
         end do
         this% ObjectsList(i)% normal        = STL% ObjectsList(i)% normal
         this% ObjectsList(i)% index         = STL% ObjectsList(i)% index
         this% ObjectsList(i)% NumOfVertices = NumOfVertices
         this% ObjectsList(i)% partition     = 1
      end do
!$omp end do
!$omp end parallel
   end subroutine STLfile_copy
  
   subroutine STLfile_GetMotionInfo( this, STLfilename, NumOfSTL )
      use FileReadingUtilities
      use FTValueDictionaryClass
      use PhysicsStorage
      implicit none
      !-arguments-------------------------------------------------------
      type(STLfile),           intent(inout) :: this
      character(len=*),        intent(in)    :: STLfilename
      integer,                 intent(in)    :: NumOfSTL
      !-local-arguments-------------------------------------------------
      integer                    :: i
      integer,       allocatable :: motionAxis_STL
      real(kind=RP), allocatable :: angularVelocity_STL, Velocity_STL, &
                                    RC_x_STL, RC_y_STL, RC_z_STL
      character(len=LINE_LENGTH) :: in_label, paramFile,  &
                                    motion_STL, STL_name 

      this% move = .false.

      do i = 1, NumOfSTL
      
         write(in_label , '(A,I0)') "#define stl motion ", i

         call get_command_argument(1, paramFile)
         call readValueInRegion ( trim ( paramFile )  , "stl name"          , STL_name           , in_label, "#end" ) 
         call readValueInRegion ( trim ( paramFile )  , "type"              , motion_STL         , in_label, "#end" ) 
         call readValueInRegion ( trim ( paramFile )  , "angular velocity"  , angularVelocity_STL, in_label, "#end" ) 
         call readValueInRegion ( trim ( paramFile )  , "velocity"          , Velocity_STL       , in_label, "#end" ) 
         call readValueInRegion ( trim ( paramFile )  , "motion axis"       , motionAxis_STL     , in_label, "#end" ) 
         call readValueInRegion ( trim ( paramFile )  , "rotation center x" , RC_x_STL           , in_label, "#end" ) 
         call readValueInRegion ( trim ( paramFile )  , "rotation center y" , RC_y_STL           , in_label, "#end" ) 
         call readValueInRegion ( trim ( paramFile )  , "rotation center z" , RC_z_STL           , in_label, "#end" ) 

         if( trim(STLfilename) .ne. trim(STL_name) ) cycle

         this% move = .true.

         select case( trim(motion_STL) )
           case( "rotation" )
               this% motionType = ROTATION
            case( "linear" )
               this% motionType = LINEAR
            case default
               print *, "STLfile_GetMotionInfo: motion not recognized. Available motions are ", ROTATION," or ",LINEAR,"."
               error stop
         end select
         
         if( allocated(angularVelocity_STL) ) then
            this% angularVelocity = angularVelocity_STL
         elseif( this% motionType .eq. ROTATION ) then
            print *, "STLfile_GetMotionInfo: 'angular velocity' must be specified for ", ROTATION, " motion."
            error stop            
         end if
         
         if( allocated(Velocity_STL) ) then
            this% Velocity = Velocity_STL
         elseif( this% motionType .eq. LINEAR ) then
            print *, "STLfile_GetMotionInfo: 'velocity' must be specified for ", LINEAR, " motion."
            error stop            
         end if
         
         if( allocated(motionAxis_STL) ) then
            this% motionAxis = motionAxis_STL
            if( this% motionAxis .gt. 3 .or. this% motionAxis .lt. 1 ) then
               print *, "STLfile_GetMotionInfo: 'motion axis' =", this% motionAxis, " not valid:"
               print *, "select 1 for x-axis, 2 for y-axis or 3 for z-axis."
               error stop                
            end if
         elseif( this% move ) then
            print *, "STLfile_GetMotionInfo: 'motion axis' must be specified."
            error stop            
         end if

         if( this% motionType .eq. ROTATION ) then 
            if( .not. allocated(RC_x_STL) ) then 
               print *, "STLfile_GetMotionInfo: 'rotation center' must be specified."
               error stop 
            end if 
            this% rotationCenter(1) = RC_x_STL
            this% rotationCenter(2) = RC_y_STL
            this% rotationCenter(3) = RC_z_STL
         end if 
         
         return

      end do
 
   end subroutine STLfile_GetMotionInfo
   
   subroutine STLfile_getRotationaMatrix( this, t, angle )
      use PhysicsStorage
      use FluidData
      implicit none
      !-arguments-----------------------------
      class(STLfile),           intent(inout):: this
      real(kind=RP),            intent(in)   :: t
      real(kind=RP),  optional, intent(in)   :: angle 
      !-local-variables-----------------------
      real(kind=RP) :: time, theta

      if( present(angle) ) then 

         this% rotationMatrix = 0.0_RP
         theta = PI/180.0_RP*angle 

         select case( this% motionAxis )
            case( IX )
               this% rotationMatrix(1,1) = 1.0_RP
               this% rotationMatrix(2,2) = cos(theta)
               this% rotationMatrix(2,3) = -sin(theta)
               this% rotationMatrix(3,2) = sin(theta)
               this% rotationMatrix(3,3) = cos(theta)
            case( IY ) 
               this% rotationMatrix(2,2) = 1.0_RP
               this% rotationMatrix(1,1) = cos(theta)
               this% rotationMatrix(1,3) = sin(theta)
               this% rotationMatrix(3,1) = -sin(theta)
               this% rotationMatrix(3,3) = cos(theta)
            case( IZ )
               this% rotationMatrix(3,3) = 1.0_RP
               this% rotationMatrix(1,1) = cos(theta)
               this% rotationMatrix(1,2) = -sin(theta)
               this% rotationMatrix(2,1) = sin(theta)
               this% rotationMatrix(2,2) = cos(theta)
         end select
         return
      end if     
#if defined(NAVIERSTOKES)
      time = t * Lref/refValues% V
   
      theta = this% angularVelocity * time
   
      this% rotationMatrix = 0.0_RP
   
      select case( this% motionAxis )
         case( IX )
            this% rotationMatrix(1,1) = 1.0_RP
            this% rotationMatrix(2,2) = cos(theta)
            this% rotationMatrix(2,3) = -sin(theta)
            this% rotationMatrix(3,2) = sin(theta)
            this% rotationMatrix(3,3) = cos(theta)
         case( IY ) 
            this% rotationMatrix(2,2) = 1.0_RP
            this% rotationMatrix(1,1) = cos(theta)
            this% rotationMatrix(1,3) = sin(theta)
            this% rotationMatrix(3,1) = -sin(theta)
            this% rotationMatrix(3,3) = cos(theta)
         case( IZ )
            this% rotationMatrix(3,3) = 1.0_RP
            this% rotationMatrix(1,1) = cos(theta)
            this% rotationMatrix(1,2) = -sin(theta)
            this% rotationMatrix(2,1) = sin(theta)
            this% rotationMatrix(2,2) = cos(theta)
      end select
#endif       
   end subroutine STLfile_getRotationaMatrix
   
   subroutine STLfile_getDisplacement( this, t )
      use FluidData
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------
      class(STLfile), intent(inout):: this
      real(kind=RP),  intent(in)   :: t

      real(kind=RP) :: time 
#if defined(NAVIERSTOKES)
      time = t * Lref/refValues% V
   
      this% ds = this% Velocity * time

      this% ds = this% ds/Lref
#endif
   end subroutine STLfile_getDisplacement

   subroutine STL_updateNormals( this )
      use MappedGeometryClass
      implicit none
      !-arguments------------------------------------
      class(STLfile), intent(inout):: this
      !-local-variables------------------------------
      real(kind=rp) :: du(NDIM), dv(NDIM), dw(NDIM)
      integer       :: i 
!$omp parallel
!$omp do schedule(runtime) private(du,dv,dw)
      do i = 1, this% NumOfObjs
         du = this% ObjectsList(i)% vertices(2)% coords - this% ObjectsList(i)% vertices(1)% coords
         dw = this% ObjectsList(i)% vertices(3)% coords - this% ObjectsList(i)% vertices(1)% coords
         call vcross(du,dw,dv)
         this% ObjectsList(i)% normal = dv/norm2(dv)      
      end do
!$omp end do
!$omp end parallel 
   end subroutine STL_updateNormals
   
   subroutine STLfile_destroy( this )
   
      implicit none
      !-arguments------------------------------
      class(STLfile), intent(inout) :: this
      !-local-variables------------------------
      integer :: i
      
      do i = 1, this% NumOfObjs
         deallocate(this% ObjectsList(i)% vertices)
      end do
      
      deallocate(this% ObjectsList)  

      this% NumOfObjs = 0
   
   end subroutine STLfile_destroy
!//////////////////////////////////////////////
!
!   Procedures form STL triangles splitting   
!   
!//////////////////////////////////////////////
! SPLITTING TRIANGLES 
   subroutine STL_Clip( this, minplane, maxplane, axis, describe )
      implicit none
      !-arguments--------------------------------------------------------------------
      class(STLfile), intent(inout) :: this
      real(kind=RP),  intent(in)    :: minplane, maxplane
      integer,        intent(in)    :: axis
      logical,        intent(in)    :: describe
      !-local-variables--------------------------------------------------------------
      type(ObjectLinkedList)     :: ObjectsLinkedList, ObjectsLinkedListFinal
      type(object_type), pointer :: obj 
      real(kind=RP)              :: minplane_point(NDIM), maxplane_point(NDIM),   & 
                                    minplane_normal(NDIM), maxplane_normal(NDIM), &
                                    Objmax, Objmin, vertices(NDIM,3)
      integer                    :: i, j

      minplane_point  = 0.0_RP
      maxplane_point  = 0.0_RP
      minplane_normal = 0.0_RP
      maxplane_normal = 0.0_RP
      
      minplane_point(axis)  = minplane
      maxplane_point(axis)  = maxplane
      minplane_normal(axis) = -1.0_RP
      maxplane_normal(axis) =  1.0_RP
 
      ObjectsLinkedList      = ObjectLinkedList_Construct()
      ObjectsLinkedListFinal = ObjectLinkedList_Construct()

      do i = 1, this% NumOfObjs
          call ClipPloy( this% ObjectsList(i), maxplane_normal, maxplane_point, ObjectsLinkedList )
      end do

      obj => ObjectsLinkedList% head 

      do i = 1, ObjectsLinkedList% NumOfObjs
         call ClipPloy( obj, minplane_normal, minplane_point, ObjectsLinkedListFinal )
         obj => obj% next 
      end do  
      
      call ObjectsLinkedList% destruct()       

      call this% destroy()

      this% partition = 1
      this% NumOfObjs = ObjectsLinkedListFinal% NumOfObjs
  
      allocate(this% ObjectsList(this% NumOfObjs))
 
      obj => ObjectsLinkedListFinal% head 

      do i = 1, this% NumOfObjs
         allocate(this% ObjectsList(i)% vertices(obj% NumOfVertices))
         do j = 1, obj% NumOfVertices
            this% ObjectsList(i)% vertices(j)% coords = obj% vertices(j)% coords
         end do
         this% ObjectsList(i)% normal        = obj% normal
         this% ObjectsList(i)% index         = i
         this% ObjectsList(i)% NumOfVertices = obj% NumOfVertices
         this% ObjectsList(i)% partition     = 1
         obj => obj% next 
      end do 

      call ObjectsLinkedListFinal% destruct()

      if( describe ) call this% describe(this% filename)

      call this% plot(0)

   end subroutine STL_Clip

   subroutine ClipPloy( obj, plane_normal, plane_point, ObjectsLinkedList )
      use MappedGeometryClass
      implicit none
      !-arguments--------------------------------------------------------------------
      type(object_type),      intent(in)    :: obj
      real(kind=rp),          intent(in)    :: plane_normal(:), plane_point(:)
      type(ObjectLinkedList), intent(inout) :: ObjectsLinkedList
      !-local-variables--------------------------------------------------------------
      real(kind=RP)     :: PointFront(NDIM,4), PointBack(NDIM,4)
      type(object_type) :: objBack
      real(kind=RP)     :: PointA(NDIM), PointB(NDIM), Point_inters(NDIM), v(NDIM),u(NDIM),W(NDIM)
      integer           :: PointA_Is, PointB_Is, n_front, n_back, i 

      n_front = 0; n_back = 0
      
      pointA = obj% vertices(obj% NumOfVertices)% coords
      
      PointA_Is = Point_wrt_Plane( plane_normal, plane_point, pointA )
   
      do i = 1, obj% NumOfVertices
         PointB    = obj% vertices(i)% coords
         PointB_Is = Point_wrt_Plane( plane_normal, plane_point, pointB )
         if( PointB_Is .eq. POINT_IN_FRONT_PLANE ) then
            if( PointA_Is .eq. POINT_BEHIND_PLANE ) then
               Point_inters = EdgePlaneIntersection( plane_normal, plane_point, PointA, PointB )
               n_front = n_front + 1
               n_back  = n_back + 1
               PointFront(:,n_front) = Point_Inters
               PointBack(:,n_back) = Point_Inters
            end if
            n_front = n_front + 1
            PointFront(:,n_front) = PointB
         elseif( PointB_Is .eq. POINT_BEHIND_PLANE ) then
            if( PointA_Is .eq. POINT_IN_FRONT_PLANE ) then
               Point_inters = EdgePlaneIntersection( plane_normal, plane_point, PointA, PointB )
               n_front = n_front + 1
               n_back  = n_back + 1
               PointFront(:,n_front) = Point_Inters
               PointBack(:,n_back) = Point_Inters
            elseif( PointA_Is .eq. POINT_ON_PLANE ) then
               n_back  = n_back + 1
               PointBack(:,n_back) = PointA
            end if
            n_back  = n_back + 1
            PointBack(:,n_back) = PointB 
         else
            n_front = n_front + 1
            PointFront(:,n_front) = PointB
            if( PointA_Is .eq. POINT_BEHIND_PLANE ) then
               n_back  = n_back + 1
               PointBack(:,n_back) = PointB 
            end if
         end if
         PointA    = PointB
         PointA_Is = PointB_Is 
      end do
      
      ! take only back elements !! 
      if( n_back .eq. 3 ) then
         call objBack% build( PointBack(:,1:n_back), obj% normal, obj% NumOfVertices, obj% index )
         call ObjectsLinkedList% add(objBack)
         call objBack% destruct()
      elseif( n_back .eq. 4 ) then
         call objBack% build( PointBack(:,1:n_back-1), obj% normal, obj% NumOfVertices, obj% index )
         call ObjectsLinkedList% add(objBack)
         call objBack% destruct()
         call objBack% build( PointBack(:,(/n_back-1,n_back,1/)), obj% normal, obj% NumOfVertices, obj% index )
         call ObjectsLinkedList% add(objBack)
         call objBack% destruct()
      elseif( n_back .eq. 0 ) then 
      else
         print *, "ClipPloy:: wrong number of vertices: ", n_back
         error stop
      end if 
   
   end subroutine ClipPloy

   subroutine STL_SetIntegration( this, NumOfInterPoints )
      use MPI_Process_Info
      implicit none 

      class(STLfile), intent(inout) :: this 
      integer,        intent(in)    :: NumOfInterPoints

      integer :: i, j

      if( .not. MPI_Process% isRoot ) return 
 
      do i = 1, this% NumOfObjs
         do j = 1, NumOfVertices + 4
            allocate( this% ObjectsList(i)% vertices(j)% nearestPoints(NumOfInterPoints),           &
                      this% ObjectsList(i)% vertices(j)% invPhi(NumOfInterPoints,NumOfInterPoints), &
                      this% ObjectsList(i)% vertices(j)% b(NumOfInterPoints)                        )
         end do 
      end do 
 
   end subroutine STL_SetIntegration

   function STL_ComputeScalarIntegral( this )

      implicit none 

      class(STLfile), intent(inout) :: this 
      real(kind=RP)                 :: STL_ComputeScalarIntegral

      real(kind=RP) :: ScalarVar(NumOfVertices+4)
      integer       :: i, j

      STL_ComputeScalarIntegral = 0.0_RP

      do i = 1, this% NumOfObjs
         associate( obj => this% ObjectsList(i) )
         do j = 1, NumOfVertices+4
            ScalarVar(j) = obj% vertices(j)% ScalarValue
         end do
         STL_ComputeScalarIntegral = STL_ComputeScalarIntegral + TriangleScalarIntegral( obj, ScalarVar )
         end associate
      end do

   end function STL_ComputeScalarIntegral

   function STL_ComputeVectorIntegral( this )

      implicit none 

      class(STLfile), intent(inout) :: this 
      real(kind=RP)                 :: STL_ComputeVectorIntegral(NDIM)

      real(kind=RP) :: VectorVar(NDIM,NumOfVertices+4)
      integer       :: i, j

      STL_ComputeVectorIntegral = 0.0_RP

      do i = 1, this% NumOfObjs
         associate( obj => this% ObjectsList(i) )
         do j = 1, NumOfVertices+4
            VectorVar(:,j) = obj% vertices(j)% VectorValue
         end do
         STL_ComputeVectorIntegral = STL_ComputeVectorIntegral + TriangleVectorIntegral( obj, VectorVar )
         end associate
      end do

   end function STL_ComputeVectorIntegral

   function TriangleScalarIntegral( obj, ScalarVar ) result( Val )
      use MappedGeometryClass
      implicit none 

      type(object_type), intent(in) :: obj 
      real(kind=RP),     intent(in) :: ScalarVar(NumOfVertices+4)
      real(kind=RP)                 :: Val

      real(kind=RP) :: AB(NDIM), AC(NDIM), S(NDIM), A

      AB = obj% vertices(2)% coords - obj% vertices(1)% coords 
      AC = obj% vertices(3)% coords - obj% vertices(1)% coords 

      call vcross(AB,AC,S)

      A = 0.5_RP * norm2(S)

      Val = A/60.0_RP * ( 27.0_RP * ScalarVar(NumOfVertices+4) +                    &
                           3.0_RP * sum(ScalarVar(1:NumOfVertices)) +               &
                           8.0_RP * sum(ScalarVar(NumOfVertices+1:NumOfVertices+3)) )

   end function TriangleScalarIntegral

   function TriangleVectorIntegral( obj, VectorVar ) result( Val )
      use MappedGeometryClass
      implicit none 

      type(object_type), intent(in) :: obj 
      real(kind=RP),     intent(in) :: VectorVar(NDIM,NumOfVertices+4)
      real(kind=RP)                 :: Val(NDIM)

      real(kind=RP) :: AB(NDIM), AC(NDIM), S(NDIM), A

      AB = obj% vertices(2)% coords - obj% vertices(1)% coords 
      AC = obj% vertices(3)% coords - obj% vertices(1)% coords 

      call vcross(AB,AC,S)

      A = 0.5_RP * norm2(S)

      Val = A/60.0_RP * ( 27.0_RP * VectorVar(:,NumOfVertices+4) +                    &
                           3.0_RP * sum(VectorVar(:,1:NumOfVertices)) +               &
                           8.0_RP * sum(VectorVar(:,NumOfVertices+1:NumOfVertices+3)) )

   end function TriangleVectorIntegral
   
   integer function Point_wrt_Plane( plane_normal, plane_point, point ) result( PointIs )
      use MappedGeometryClass
      implicit none
      !-arguments-----------------------------------------------------------------
      real(kind=RP), dimension(:), intent(in) :: plane_normal, plane_point, point
      !-local-variables-----------------------------------------------------------
      real(kind=RP) :: d, eps
   
      d = dot_product(plane_normal,point) - dot_product(plane_normal,plane_point)
      
      if( d > EPSILON(eps) ) then
         PointIs = POINT_IN_FRONT_PLANE
      elseif( d < -EPSILON(eps) ) then
         PointIs = POINT_BEHIND_PLANE
      else 
         PointIs = POINT_ON_PLANE
      end if
   
   end function Point_wrt_Plane
  
   function EdgePlaneIntersection( plane_normal, plane_point, PointA, PointB ) result( Point_inters )
      use MappedGeometryClass
      implicit none
      !-arguments-----------------------------------------------------------------
      real(kind=RP), intent(in) :: plane_normal(:), plane_point(:), &
                                   PointA(:), PointB(:)
      real(kind=RP)             :: Point_inters(NDIM)
      !-local-variables-----------------------------------------------------------
      real(kind=RP) :: B_A(NDIM), d, t
 
      B_A = PointB - PointA
      
      d = dot_product(plane_normal,plane_point)
      
      t = ( d - dot_product(plane_normal,PointA) )/dot_product(plane_normal,B_A)
 
      Point_inters = PointA + t*B_A 
 
   end function EdgePlaneIntersection
!//////////////////////////////////////////////
   
   subroutine TecFileHeader( FileName, Title, I, J, K, funit, DATAPACKING, ZONETYPE )
   
      implicit none
      !-arguments------------------------------------------------------
      character(len=*), intent(in)  :: FileName, Title, DATAPACKING
      integer,          intent(in)  :: I, J, K
      character(len=*), optional    :: ZONETYPE
      integer,          intent(out) :: funit
   
      funit = UnusedUnit()
   
      open(funit,file=trim(FileName)//'.tec', status='unknown')
      
      write(funit,"(A9,A,A)") 'TITLE = "',trim(Title),'"'
      write(funit,"(A23)") 'VARIABLES = "x","y","z"'
      if( present(ZONETYPE) ) then
         write(funit,"(A7,I0,A3,I0,A3,I0,A14,A,A11,A)") 'ZONE I=',I,',J=',J,',K=',K,', DATAPACKING=',trim(DATAPACKING),', ZONETYPE=',trim(ZONETYPE)
      else
         write(funit,"(A7,I0,A3,I0,A3,I0,A14,A)") 'ZONE I=',I,',J=',J,',K=',K,', DATAPACKING=',trim(DATAPACKING)
      end if
      
   end subroutine TecFileHeader

   subroutine PolynomialVector( x, y, z, v )
   
      implicit none 

      real(kind=rp), intent(in)    :: x, y, z 
      real(kind=rp), intent(inout) :: v(:)

      integer :: n_of_points 

      n_of_points = size(v)

      select case( n_of_points )  
         case( 3 )                    !1st order
            v(1) = 1.0_RP 
            v(2) = x 
            v(3) = y 
            !v(4) = z 
         case( 6 )                    !2nd order
            v(1)  = 1.0_RP 
            v(2)  = x 
            v(3)  = y 
            !v(4)  = z
            v(4)  = POW2(x) 
            v(5)  = POW2(y) 
            !v(7)  = POW2(z)
            v(6)  = x*y          
            !v(9)  = x*z          
            !v(10) = z*y 
         case( 10 )                    !3rd order
            v(1)  = 1.0_RP 
            v(2)  = x 
            v(3)  = y 
            !v(4)  = z
            v(4)  = POW2(x) 
            v(5)  = POW2(y) 
            !v(7)  = POW2(z)
            v(6)  = x*y          
            !v(9)  = x*z          
            !v(10) = z*y 
            v(7) = POW3(x) 
            v(8) = POW3(y) 
            !v(13) = POW3(z)
            v(9) = POW2(x)*y          
            !v(15) = POW2(x)*z          
            v(10) = POW2(y)*x         
            !v(17) = POW2(y)*z         
            !v(18) = POW2(z)*x         
            !v(19) = POW2(z)*y         
            !v(20) = x*y*z         
         case( 17 )                    !4th order
            v(1)  = 1.0_RP 
            v(2)  = x 
            v(3)  = y 
            !v(4)  = z
            v(4)  = POW2(x) 
            v(5)  = POW2(y) 
            !v(7)  = POW2(z)
            v(6)  = x*y          
            !v(9)  = x*z          
            !v(10) = z*y 
            v(7) = POW3(x) 
            v(8) = POW3(y) 
            !v(13) = POW3(z)
            v(9) = POW2(x)*y          
            !v(15) = POW2(x)*z          
            v(10) = POW2(y)*x         
            !v(17) = POW2(y)*z         
            v(11) = POW2(z)*x         
            v(12) = POW2(z)*y         
            !v(20) = x*y*z
            v(13) = POW2(x)*POW2(x)         
            v(14) = POW2(y)*POW2(y)         
            !v(23) = POW2(z)*POW2(z)
            v(15) = POW3(x)*y         
            !v(25) = POW3(x)*z
            v(16) = POW3(y)*x         
            !v(27) = POW3(y)*z
            !v(28) = POW3(z)*x         
            !v(29) = POW3(z)*y               
            !v(30) = POW2(x)*y*z          
            !v(31) = POW2(y)*x*z          
            !v(32) = POW2(z)*x*y 
            v(17) = POW2(x)*POW2(y) 
            !v(34) = POW2(x)*POW2(z)          
            !v(35) = POW2(y)*POW2(z)            
      end select 

   end subroutine PolynomialVector

   subroutine buildMatrixPolySpline( Phi, P, P_T, V )

      implicit none 

      real(kind=rp), intent(in)    :: Phi(:,:), P(:,:), P_T(:,:)
      real(kind=rp), intent(inout) :: V(:,:)

      integer :: start, final 

      V = 0.0_RP 
      
      V(1:size(Phi,1),1:size(Phi,2)) = Phi 
      start = size(Phi,2)+1; final = start + size(P,2)-1
      V(1:size(Phi,1),start:final) = P 
      start = size(Phi,1)+1; final = start + size(P_T,1)-1
      V(start:final,1:size(Phi,2)) = P_T  

   end subroutine buildMatrixPolySpline

   real(kind=rp) function EvaluateInterp( coeff, x, y, z ) result( value )

      implicit none
      
      real(kind=rp), intent(in) :: coeff(:), x, y, z
      
      real(kind=rp) :: v(size(coeff))

      call PolynomialVector( x, y, z, v )

      value = dot_product(coeff,v)

   end function EvaluateInterp

end module TessellationTypes
