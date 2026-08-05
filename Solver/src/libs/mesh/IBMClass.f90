#include "Includes.h"
module IBMClass

   use SMConstants
   use Utilities
   use FTValueDictionaryClass
   use ElementClass
   use FaceClass
   use TessellationTypes
   use OrientedBoundingBox
   use KDClass
   use IntegerDataLinkedList
   use MPI_IBMUtilities
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none

   type :: Integral_t

      logical              :: ListComputed = .false., &
                              compute      = .false., &
                              constructed  = .false.
   end type

   type IBM_type

      type(STLfile),              allocatable :: stl(:), stlSurfaceIntegrals(:), stlMove(:)
      type(KDtree),               allocatable :: root(:), rootDistance(:), rootPoints(:)
      type(PointLinkedList)                   :: BandPoints
      type(IBMpoints),            allocatable :: BandRegion(:), BandRegion4Distance(:)
      type(point_type),           allocatable :: ImagePoints(:)
      type(Integral_t),           allocatable :: Integral(:)
      character(len=LINE_LENGTH), allocatable :: STLfilename(:)
      character(len=LINE_LENGTH)              :: filename
      logical                                 :: plotOBB              = .false., &
                                                 plotKDtree           = .false., &
                                                 active               = .false., & 
                                                 TimePenal            = .false., &
                                                 semiImplicit         = .false., &
                                                 ComputeBandRegion    = .false., &
                                                 plotBandPoints       = .false., &
                                                 plotMask             = .false., &
                                                 ComputeInterpolation = .false., &
                                                 Wallfunction         = .false., &
                                                 ComputeDistance      = .false., &
                                                 AAB                  = .false.                                                
      real(kind=rp)                           :: eta, BandRegionCoeff, IP_Distance = 0.0_RP, &
                                                 y_plus_target, minCOORDS, maxCOORDS,        &
                                                 penalCoeff
      real(kind=rp),              allocatable :: penalization(:)
      integer                                 :: KDtree_Min_n_of_Objs, NumOfInterPoints,     &
                                                 n_of_INpoints,  rank, lvl = 0, NumOfSTL,    &
                                                 NumOfForcingPoints, Clipaxis = 0,           &
                                                 Nx, Ny, Nz, LocClipAxis = 0,                &
                                                 InterpolationType
      integer,                    allocatable :: ImagePoint_NearestPoints(:,:)
      
      contains  
         procedure :: read_info                           => IBM_read_info
         procedure :: construct                           => IBM_construct
         procedure :: constructMask                       => IBM_constructMask
         procedure :: constructSTL_KDtree                 => IBM_constructSTL_KDtree
         procedure :: CheckPoint                          => IBM_checkPoint
         procedure :: constructBandRegion                 => IBM_constructBandRegion
         procedure :: constructBandRegion4Distance        => IBM_constructBandRegion4Distance
         procedure :: build                               => IBM_build
         procedure :: SetPolynomialOrder                  => IBM_SetPolynomialOrder
         procedure :: GetMask                             => IBM_GetMask
         procedure :: MPI_sendOBB                         => IBM_MPI_sendOBB
         procedure :: MPI_sendSTLpartitions               => IBM_MPI_sendSTLpartitions
         procedure :: MPI_sendMask2Root                   => IBM_MPI_sendMask2Root
         procedure :: MPI_sendMask2Partitions             => IBM_MPI_sendMask2Partitions
         procedure :: MPI_sendNormals2Root                => IBM_MPI_sendNormals2Root
         procedure :: MPI_sendDistNormals2partitions      => IBM_MPI_sendDistNormals2partitions
         procedure :: BandRegionPoints                    => IBM_bandRegionPoints
         procedure :: GetForcingPointsGeom                => IBM_GetForcingPointsGeom
         procedure :: GetInfo                             => IBM_GetInfo 
         procedure :: SourceTerm                          => IBM_SourceTerm
         procedure :: ComputeIBMWallDistance              => IBM_ComputeIBMWallDistance
         procedure :: GetDistanceInsideBox                => IBM_GetDistanceInsideBox
         procedure :: GetDistanceOutsideBox               => IBM_GetDistanceOutsideBox
         procedure :: SemiImplicitCorrection              => IBM_SemiImplicitCorrection
         procedure :: GetImagePoint_nearest               => IBM_GetImagePoint_nearest
         procedure :: GetBandRegionStates                 => IBM_GetBandRegionStates
         procedure :: GetDomainExtreme                    => IBM_GetDomainExtreme
         procedure :: SourceTermTurbulence                => IBM_SourceTermTurbulence
         procedure :: semiImplicitShiftJacobian           => IBM_semiImplicitShiftJacobian
         procedure :: semiImplicitTurbulenceShiftJacobian => IBM_semiImplicitTurbulenceShiftJacobian
         procedure :: semiImplicitJacobian                => IBM_semiImplicitJacobian
         procedure :: semiImplicitTurbulenceJacobian      => IBM_semiImplicitTurbulenceJacobian  
         procedure :: GetSemiImplicitStep                 => IBM_GetSemiImplicitStep   
         procedure :: GetSemiImplicitStepTurbulence       => IBM_GetSemiImplicitStepTurbulence   
         procedure :: SetIntegration                      => IBM_SetIntegration
         procedure :: copy                                => IBM_copy
         procedure :: MoveBody                            => IBM_MoveBody
         procedure :: CleanMask                           => IBM_CleanMask
         procedure :: BandPoint_state                     => IBM_BandPoint_state
         procedure :: Describe                            => IBM_Describe
         procedure :: plot_Mask                           => IBM_plot_Mask
         procedure :: Destruct                            => IBM_Destruct
         procedure :: DestroyKDtree                       => IBM_DestroyKDtree
         procedure :: constructDistance_KDtree            => IBM_constructDistance_KDtree
         procedure :: MPI_PointsListOperations            => IBM_MPI_PointsListOperations
         procedure :: MaskVelocity                        => IBM_MaskVelocity
   end type

   public :: expCoeff, EXPONENTIAL

   real(kind=RP)      :: expCoeff
   integer, parameter :: EXPONENTIAL = 1, IDW = 2, POLYHARMONIC_SPLINE = 3, POLYNOMIAL = 4, MLS = 5

   contains
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  Immersed boundary info
!  -------------------------------------------------   
   subroutine IBM_read_info( this, controlVariables )
      use FileReadingUtilities
      implicit none
      
      class(IBM_type), intent(inout) :: this
      class(FTValueDictionary)       :: controlVariables
 
      call this% GetInfo( controlVariables )

   end subroutine IBM_read_info
   

   subroutine IBM_GetInfo( this, controlVariables )
      use FileReadingUtilities
#if defined(NAVIERSTOKES)
      use WallFunctionDefinitions  
#endif
      implicit none
      !-arguments----------------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      class(FTValueDictionary)       :: controlVariables
      !-local-variables----------------------------------------------------------
      logical,       allocatable :: active_in, plotOBB_in,              &
                                    plotKDtree_in, semiImplicit_in,     &
                                    plotBandPoints_in, plotMask_in,     &
                                    BandRegion_in, Distance_in, AAB_in
      real(kind=rp), allocatable :: penalization_in, y_plus_target_in,  &
                                    BandRegionCoeff_in,                 & 
                                    penalization_coeff_in
      integer,       allocatable :: n_of_Objs_in, n_of_interpoints_in,  &
                                    Nx_in, Ny_in, Nz_in, Clipaxis_in
      character(len=LINE_LENGTH) :: in_label, paramFile, name_in, tmp,  &
                                    InterpolationType_in
      real(kind=rp), allocatable :: coords(:)   
      logical                    :: correct
      
      character(len=LINE_LENGTH), parameter :: NumberOfSTL = "number of stl" 
      
!     Read block
!     **********
      write(in_label , '(A)') "#define ibm"
      call get_command_argument(1, paramFile) 
      call readValueInRegion( trim( paramFile ), "name",                           name_in,               in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "active",                         active_in,             in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "penalization",                   penalization_in,       in_label, "#end" )   
      call readValueInRegion( trim( paramFile ), "penalization coeff",             penalization_coeff_in, in_label, "#end" )   
      call readValueInRegion( trim( paramFile ), "semi implicit",                  semiImplicit_in,       in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "nx",                             Nx_in,                 in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "ny",                             Ny_in,                 in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "nz",                             Nz_in,                 in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "target y plus",                  y_plus_target_in,      in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "number of objects",              n_of_Objs_in,          in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "number of interpolation points", n_of_interpoints_in,   in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "clip axis",                      Clipaxis_in,           in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "band region",                    BandRegion_in,         in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "distance",                       Distance_in,           in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "band region coeff",              BandRegionCoeff_in,    in_label, "#end" )    
      call readValueInRegion( trim( paramFile ), "aab",                            AAB_in,                in_label, "#end" )     
      call readValueInRegion( trim( paramFile ), "intepolation",                   InterpolationType_in,  in_label, "#end" )     
      call readValueInRegion( trim( paramFile ), "plot obb",                       plotOBB_in,            in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "plot kdtree" ,                   plotKDtree_in,         in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "plot mask",                      plotMask_in,           in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "plot band points",               plotBandPoints_in,     in_label, "#end" )

      this% filename = trim(name_in)

      if( allocated(active_in) ) then
         this% active = active_in
      else
         this% active = .FALSE.
      end if
      
      if( .not. this% active) return
      
      if( allocated(semiImplicit_in) ) then
         this% semiImplicit = semiImplicit_in
      else
         this% semiImplicit = .FALSE.
      end if       

      if( allocated(Nx_in) ) then
         this% Nx = Nx_in 
      else
         this% Nx = 0
      end if   
           
      if( allocated(Ny_in) ) then
         this% Ny = Ny_in
      else
         this% Ny = 0
      end if 
             
      if( allocated(Nz_in) ) then
         this% Nz = Nz_in
      else
         this% Nz = 0
      end if           
   
      if( allocated(plotOBB_in) ) then
         this% plotOBB = plotOBB_in
      else
         this% plotOBB = .FALSE.
      end if

      if( allocated(plotKDtree_in) ) then
         this% plotKDtree = plotKDtree_in
      else
         this% plotKDtree = .FALSE.
      end if

      if( allocated(penalization_in) ) then
         this% eta = penalization_in
         this% TimePenal = .false.
      else
         this% eta = 0.1_RP
         this% TimePenal = .true.
      end if

      if( allocated(penalization_coeff_in) ) then
         this% penalCoeff = penalization_coeff_in
      else
         this% penalCoeff = 1.0_RP
      end if

      if( allocated(n_of_Objs_in) ) then
         this% KDtree_Min_n_of_Objs = n_of_Objs_in
      else
         this% KDtree_Min_n_of_Objs = 5
      end if
      
      if( allocated(n_of_interpoints_in) ) then
         this% NumOfInterPoints = n_of_interpoints_in
      else
         this% NumOfInterPoints = 15
      end if 
      
      if( allocated(plotBandPoints_in) ) then
         this% plotBandPoints = plotBandPoints_in
      else
         this% plotBandPoints = .false.
      end if
      
      if( allocated(plotMask_in) ) then
         this% plotMask = plotMask_in
      else
         this% plotMask = .false.
      end if
      
      if( allocated(BandRegion_in) ) then
         this% ComputeBandRegion = BandRegion_in
      else
         this% ComputeBandRegion = .false.
      end if
      
      if( allocated(Distance_in) ) then
         this% ComputeDistance = Distance_in
      else
         this% ComputeDistance = .false.
      end if

      if( allocated(Clipaxis_in) ) then 
         this% ClipAxis = Clipaxis_in
      else
         this% ClipAxis = 0
      endif

     if( controlVariables% containsKey("wall function") ) then
        this% Wallfunction = .true.
#if defined(NAVIERSTOKES)
        call Initialize_Wall_Function(controlVariables, correct)   !TO BE REMOVED
#endif
        if( allocated(y_plus_target_in) ) then
           this% y_plus_target = y_plus_target_in
        else
           this% y_plus_target = 50.0_RP
        end if
        
        this% ComputeBandRegion = .true.
        this% ComputeDistance   = .true.
   
     else
        this% Wallfunction = .false.
     end if

      if( allocated(BandRegionCoeff_in) ) then
         this% BandRegionCoeff = BandRegionCoeff_in 
      else
         this% BandRegionCoeff = 2.0_RP      
      end if
 
      if( controlVariables% containsKey(trim(NumberOfSTL)) ) then
         tmp = controlVariables% StringValueForKey(trim(NumberOfSTL),LINE_LENGTH) 
         this% NumOfSTL = GetIntValue(tmp)
      else
         this% NumOfSTL = 1
      end if          

      if( allocated(AAB_in) ) then 
         this% AAB = AAB_in
      else
         this% AAB = .false.
      end if

      select case(trim(InterpolationType_in))
         case("exp")
            this% InterpolationType =  EXPONENTIAL
         case("idw")
            this% InterpolationType =  IDW
         case("spline")
            this% InterpolationType =  POLYHARMONIC_SPLINE
         case("polynomial")
            this% InterpolationType =  POLYNOMIAL
         case("mls")
            this% InterpolationType =  MLS
         case default 
            this% InterpolationType =  IDW
      end select 

   end subroutine IBM_GetInfo
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  KDtree procedures 
!  -------------------------------------------------   
   subroutine IBM_construct( this, controlVariables )   
      use mainKeywordsModule
      use FileReadingUtilities
      use MPI_Process_Info
      implicit none
      !-arguments----------------------------------------
      class(IBM_type), intent(inout) :: this  
      class(FTValueDictionary)       :: controlVariables
      !-local-variables----------------------------------
      character(len=LINE_LENGTH) :: filename, MyString
      real(kind=RP)              :: axis(NDIM)
      integer                    :: STLNum, j, k  

      call this% describe()         

      allocate( this% stl(this% NumOfSTL),                 & 
                OBB(this% NumOfSTL),                       &
                this% root(this% NumOfSTL),                &
                this% integral(this% NumOfSTL),            &
                this% STLfilename(this% NumOfSTL),         &
                this% stlSurfaceIntegrals(this% NumOfSTL), &
                this% stlMove(this% NumOfSTL)              )

      if( this% ComputeBandRegion ) then
         allocate( this% rootPoints(this% NumOfSTL),   & 
                   this% BandRegion(this% NumOfSTL)    )
      end if

      if( this% ComputeDistance ) allocate(this% rootDistance(this% NumOfSTL))
      
      do STLNum = 1, this% NumOfSTL 
         write(MyString, '(i100)') STLNum
         if( STLNum .eq. 1 ) then
            filename = stlFileNameKey
         else
            filename = trim(stlFileNameKey)//trim(adjustl(MyString))
         end if          
         this% STLfilename(STLNum) = controlVariables% stringValueForKey(trim(filename), requestedLength = LINE_LENGTH)
         OBB(STLNum)% filename     = this% STLfilename(STLNum)
         call STLfile_GetMotionInfo( this% stl(STLNum), this% STLfilename(STLNum), this% NumOfSTL ) 
         if( MPI_Process% isRoot ) then
            this% stl(STLNum)% show = .true. 
            call this% stl(STLNum)% ReadTessellation( this% STLfilename(STLNum) ) 
            if( this% ClipAxis .ne. 0 ) then 
               call this% stl(STLNum)% Clip( this% minCOORDS, this% maxCOORDS, this% ClipAxis, .true. ) 
            else
               call this% stl(STLNum)% describe(this% STLfilename(STLNum))
            end if
            if( this% stl(STLNum)% move ) call this% stlMove(STLNum)% copy( this% stl(STLNum) )
            call OBB(STLNum)% construct( this% stl(STLNum), this% plotOBB, this% AAB )
            call OBB(STLNum)% ChangeObjsRefFrame( this% stl(STLNum)% ObjectsList, LOCAL )  
         end if
         call this% MPI_sendOBB(STLNum)  
         call this% constructSTL_KDtree( STLNum )                
      end do

   end subroutine IBM_Construct
   
   subroutine IBM_constructSTL_KDtree( this, STLNum )
      use MPI_Process_Info 
      implicit none
      !-arguments-----------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum
      !-local-variables-----------------------------
      real(kind=RP) :: vertices(NDIM,8)

      this% root(STLNum)% STLNum       = STLNum
      this% root(STLNum)% which_KDtree = TRIANGLES_KDTREE_SAH

      vertices = OBB(STLNum)% LocVertices

      call this% MPI_sendSTLpartitions( STLNum, vertices )

      call this% root(STLNum)% construct( stl           = this% stl(STLNum),         &
                                          vertices      = vertices,                  &
                                          isPlot        = this% plotKDtree,          & 
                                          Min_n_of_Objs = this% KDtree_Min_n_of_Objs )

      if( this% ComputeDistance ) call this% constructDistance_KDtree( STLNum ) 

      call this% stl(STLNum)% destroy() 

   end subroutine IBM_constructSTL_KDtree
   
   subroutine IBM_constructDistance_KDtree( this, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments----------------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum
      !-local-variables----------------------------------
      real(kind=RP) :: BRvertices(NDIM,8)

      this% rootDistance(STLNum)% STLNum       = STLNum      
      this% rootDistance(STLNum)% which_KDtree = TRIANGLES_KDTREE_MEDIAN      

      call GetBRvertices( this% root(STLNum)% vertices, this% BandRegionCoeff, this% root(STLNum)% MaxAxis, STLNum, BRvertices )

      call this% rootDistance(STLNum)% construct( stl           = this% stl(STLNum),         &
                                                  Vertices      = BRvertices,                &
                                                  isPlot        = .false.,                   & 
                                                  Min_n_of_Objs = this% KDtree_Min_n_of_Objs )

   end subroutine IBM_constructDistance_KDtree
!  
!  Copy a KD tree   
!  ---------------
   subroutine IBM_copy( this, parent, lvl )
   
      implicit none
      !-arguments----------------------------------------------------------------
      class(IBM_type),         intent(inout) :: this
      type(IBM_type),  target, intent(in)    :: parent 
      integer,                 intent(in)    :: lvl
      !-local-variables----------------------------------------------------------
      integer :: STLNum
   
      allocate( this% root(parent% NumOfSTL),       &
                this% STLfilename(parent% NumOfSTL) )

      if( parent% ComputeBandRegion ) then
         allocate( this% rootPoints(parent% NumOfSTL), & 
                   this% BandRegion(parent% NumOfSTL)  )
      end if
      
      if( parent% ComputeDistance ) allocate(this% rootDistance(parent% NumOfSTL))

      do STLNum = 1, parent% NumOfSTL
         this% STLfilename(STLNum) = parent% STLfilename(STLNum)
         this% root(STLNum)        = parent% root(STLNum)
         if( parent% ComputeDistance ) this% rootDistance(STLNum) = parent% rootDistance(STLNum)
      end do
      
      this% ClipAxis = parent% ClipAxis

      this% lvl = lvl

   end subroutine IBM_copy
!  
!  Destroy the KD tree   
!  --------------------
   subroutine IBM_DestroyKDtree( this, isChild, DistanceKDtree )
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------------------------------
      class(IBM_type),           intent(inout) :: this
      logical,                   intent(in)    :: isChild
      logical,         optional, intent(in)    :: DistanceKDtree
      !-local-variables--------------------------------------------
      integer :: STLNum
      
      if(this% ComputeDistance ) then
         do STLNum = 1, this% NumOfSTL
            call this% rootDistance(STLNum)% destruct( isChild )
         end do 
         deallocate(this% rootDistance)
         return      
      end if
        
      do STLNum = 1, this% NumOfSTL
         call this% root(STLNum)% destruct( isChild )
      end do
   
   end subroutine IBM_DestroyKDtree
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  Mask procedures 
!  -------------------------------------------------   
  subroutine IBM_GetMask( this, elements, no_of_elements, no_of_DoFs, STLNum, iter )
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: no_of_elements, no_of_DoFs, STLNum, iter

      call GetMaskCandidates( elements, no_of_elements, no_of_DoFs, STLNum, this% NumOfSTL )

      call this% MPI_PointsListOperations( Mask )

      call this% constructmask( elements, STLNum )

      if( this% plotMask ) call this% plot_Mask( iter, STLNum )
 
      deallocate( Mask% x )

   end subroutine IBM_GetMask
!  
! 
!  mask construction
!  -----------------------------------   
   subroutine IBM_constructmask( this, elements, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments---------------------------------------------------- 
      class(IBM_type),  intent(inout) :: this
      integer,          intent(in)    :: STLNum
      type(element),    intent(inout) :: elements(:)
      !-local-variables----------------------------------------------
      real(kind=RP) :: Point(NDIM) 
      integer       :: eID, n, i, j, k
!$omp parallel 
!$omp do schedule(runtime) private(Point)
      do n = 1, Mask% NumOfObjs    
         call OBB(STLNum)% ChangeRefFrame( Mask% x(n)% coords, LOCAL, Point )
         Mask% x(n)% isInsideBody       = .false.
         Mask% x(n)% NumOfIntersections = 0
         if( isInsideBox( Point, this% root(STLNum)% vertices ) ) then
            call this% CheckPoint( Point, STLNum, Mask% x(n)% NumOfIntersections )  
         end if   
      end do
!$omp end do
!$omp end parallel
      call this% MPI_sendMask2Root()       
      call this% MPI_sendMask2Partitions()
!$omp parallel 
!$omp do schedule(runtime) private(eID,i,j,k)
       do n = 1, Mask% NumOfObjs
          if( Mask% x(n)% partition .eq. MPI_Process% rank ) then
             eID = Mask% x(n)% element_index
             i = Mask% x(n)% local_Position(1)
             j = Mask% x(n)% local_Position(2)
             k = Mask% x(n)% local_Position(3)
             elements(eID)% isInsideBody(i,j,k) = Mask% x(n)% isInsideBody
             if( elements(eID)% isInsideBody(i,j,k) ) then 
                 elements(eID)% STL(i,j,k) = STLNum
             end if 
         end if    
      end do
!$omp end do
!$omp end parallel
      if( MPI_Process% isRoot ) then 
         if( Mask% NumOfObjs .eq. 0 .and. this% lvl .gt. 0 ) then
            print *, "The mask for the multigrid level ", this% lvl, " is made of 0 points."
            print *, "Try to increase the polynomial order or to refine the mesh."
            error stop
         elseif( Mask% NumOfObjs .eq. 0 ) then
            print *, "The mask is made of 0 points."
            print *, "Try to increase the polynomial order or to refine the mesh."
            error stop         
         end if
      end if
   end subroutine IBM_constructmask
!  
!  Mask plot
!  -----------------------------------------------------------
   subroutine IBM_plot_Mask(  this, iter, STLNum)
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: iter, STLNum
      !-local-variables------------------------------------------------
      character(len=LINE_LENGTH) :: filename
      integer                    :: i, funit, NumOfObjs
      logical                    :: add_Point = .false.
      
      if( .not. MPI_Process% isRoot ) return

      NumOfObjs = 0
      do i = 1, Mask% NumOfObjs
         if( Mask% x(i)% isInsideBody ) NumOfObjs = NumOfObjs + 1
      end do

      if( NumOfObjs .eq. 0 ) then 
         print*, "Mask is made of 0 points"
         error stop 
      end if

      if( this% lvl .gt. 0 ) then
         write(filename,'(A,A,I1,A,I10.10)') trim(this% STLfilename(STLNum)),'_MGlevel',this% lvl,'_',iter
      else
         write(filename,'(A,A,I10.10)') trim(this% STLfilename(STLNum)),'_',iter
      end if

      call TecFileHeader( 'IBM/Mask_'//trim(filename), 'Mask Points', NumOfObjs/2+mod(NumOfObjs,2),2,1, funit, 'POINT')
     
      if( mod(NumOfObjs,2) .ne. 0 )  add_Point  = .true.

      do i = 1, Mask% NumOfObjs
         if( Mask% x(i)% isInsideBody ) then
            write(funit,'(3E13.5)')  Mask% x(i)% coords(1), Mask% x(i)% coords(2), Mask% x(i)% coords(3)
            if( add_Point ) then
               write(funit,'(3E13.5)')  Mask% x(i)% coords(1), Mask% x(i)% coords(2), Mask% x(i)% coords(3)
               add_Point = .false.
            end if
         end if 
      end do

      close(funit)
      
   end subroutine IBM_plot_Mask
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -----------------------------------------------------------------------
!  Getting extreme of the domain mesh. These are used to clip the STL file
!  -----------------------------------------------------------------------   
   
   subroutine IBM_GetDomainExtreme( this, elements )
      use MPI_Process_Info
      implicit none
      !-arguments---------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(in)   :: elements(:)
      !-local-variables---------------------------------------
      real(kind=rp) :: ElemMax, ElemMin
      integer       :: eID, i, axis, j, k 
#ifdef _HAS_MPI_
      real(kind=RP) :: localmax, localmin
      integer       :: ierr
#endif
      this% maxCOORDS = -huge(1.0_RP); this% minCOORDS = huge(1.0_RP)
      axis = this% ClipAxis
      if( axis .eq. 0 ) return
      do eID = 1, size(elements)
         do k = 0, elements(eID)% Nxyz(3)   ; do j = 0, elements(eID)% Nxyz(2) ; do i = 0, elements(eID)% Nxyz(1)
            this% maxCOORDS = max(this% maxCOORDS,elements(eID)% geom% x(axis,i,j,k)); this% minCOORDS = min(this% minCOORDS,elements(eID)% geom% x(axis,i,j,k))
         end do; end do; end do
      end do
      this% maxCOORDS = this% maxCOORDS + 1.0e-8
      this% minCOORDS = this% minCOORDS - 1.0e-8
#ifdef _HAS_MPI_
      localmax = this% maxCOORDS; localmin = this% minCOORDS
      call mpi_allreduce(localmax, this% maxCOORDS, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
      call mpi_allreduce(localmin, this% minCOORDS, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif
   end subroutine IBM_GetDomainExtreme
!
!  Allocation and coordinates for the integration points
!  -----------------------------------------------------  
   subroutine IBM_SetIntegration( this, STLNum )
      use PhysicsStorage
      implicit none
      !-arguments------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum
      !-local-variables------------------------
      integer :: i, j 

      if( this% Integral(STLNum)% compute ) return

      allocate(this% BandRegion(STLNum)% U_x(NCONS,this% BandRegion(STLNum)% NumOfObjs))
      allocate(this% BandRegion(STLNum)% U_y(NCONS,this% BandRegion(STLNum)% NumOfObjs))
      allocate(this% BandRegion(STLNum)% U_z(NCONS,this% BandRegion(STLNum)% NumOfObjs))

   end subroutine IBM_SetIntegration
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  building the immersed boundary 
!  -------------------------------------------------    
   subroutine IBM_build( this, elements, no_of_elements, no_of_DoFs, isChild, movingSTL, iter )
      use MPI_Process_Info
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------------------------------------
      class(IBM_type),           intent(inout) :: this
      type(element),             intent(inout) :: elements(:)
      integer,                   intent(in)    :: no_of_elements, no_of_DoFs
      logical,                   intent(in)    :: isChild       
      integer,         optional, intent(in)    :: movingSTL, iter       
      !-local-variables-----------------------------------------------------------
      integer :: MaskPoints, STLNum
#ifdef _HAS_MPI_
      integer :: localVal, ierr
#endif
      do STLNum = 1, this% NumOfSTL
         if( .not. present(movingSTL) )then 
            call this% GetMask( elements, no_of_elements, no_of_DoFs, STLNum, 0 )
         else
            if( STLNum .eq. movingSTL ) call this% GetMask( elements, no_of_elements, no_of_DoFs, STLNum, iter )
         end if
      end do  

      if( .not. present(movingSTL) )then
         if (.not. allocated(this% penalization)) then
            allocate( this% penalization(no_of_elements) )
            this% penalization = this% eta
         end if 
      end if 
#if defined(NAVIERSTOKES)
      if( this% ComputeDistance ) then
         if( .not. present(movingSTL) ) then
            call this% ComputeIBMWallDistance( elements )
         else
            call this% ComputeIBMWallDistance( elements, movingSTL )
         end if
      end if
#endif   
      if( this% ComputeBandRegion ) then
         do STLNum = 1, this% NumOfSTL
            if( .not. present(movingSTL) ) then
               call this% constructBandRegion( elements, no_of_elements, STLNum ) 
            else
               if( STLNum .eq. movingSTL ) then
                  call this% constructBandRegion( elements, no_of_elements, STLNum ) 
               end if
            end if
         end do 
      end if

      if( this% Wallfunction ) then
#if defined(NAVIERSTOKES)
         call this% GetImagePoint_nearest( elements )
#endif
      end if

   end subroutine IBM_build
!  
!  Immersed boundary description
!  -----------------------------
   subroutine IBM_Describe( this )
      use Headers
      use mainKeywordsModule
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------
      class(IBM_type),  intent(inout) :: this
      
      if ( .not. MPI_Process % isRoot ) return
      write(STD_OUT,'(/)')
      call Section_Header("IBM parameters")
      write(STD_OUT,'(/)')
      
      call SubSection_Header('IBM info')
      write(STD_OUT,'(30X,A,A35,L10)') "->" , "Semi implicit treatment: ", this% semiImplicit
      if( .not. this% TimePenal ) then
         write(STD_OUT,'(30X,A,A35,ES14.2)') "->" , "Penalization term: " , this% eta
      else
         write(STD_OUT,'(30X,A,A35,A10)') "->" , "Penalization term: ", " Dt"
      end if

      write(STD_OUT,'(30X,A,A35,I10)') "->" , "Minimum number of objects: ", this% KDtree_Min_n_of_Objs
      write(STD_OUT,'(30X,A,A35,I10)') "->" , "Number of interpolation points: ", this% NumOfInterPoints

   end subroutine IBM_Describe
!
!  Destroy immersed boundary
!  -------------------------
   subroutine IBM_Destruct( this, isChild )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------
      class(IBM_type), intent(inout) :: this
      logical,         intent(in)    :: isChild
      !----------------------------------------
      integer :: STLNum, i, j 
      
      do STLNum = 1, this% NumOfSTL
         call this% root(STLNum)% destruct( isChild )
         if( this% ComputeDistance ) call this% rootDistance(STLNum)% destruct( isChild ) 
         if( this% ComputeBandRegion ) then
            call this% rootPoints(STLNum)% Destruct( .false. )
            deallocate( this% BandRegion(STLNum)% x, &
                        this% BandRegion(STLNum)% Q  )
            if( this% Integral(STLNum)% compute ) then
               deallocate( this% BandRegion(STLNum)% U_x, &
                           this% BandRegion(STLNum)% U_y, &
                           this% BandRegion(STLNum)% U_z  )
            end if
         end if
         if( this% Integral(STLNum)% compute ) call this% stlSurfaceIntegrals(STLNum)% destroy()
         if( this% stl(STLNum)% move ) call this% stlMove(STLNum)% destroy()
      end do
      
      if( this% Wallfunction ) then 
         do i = 1, this% NumOfForcingPoints
            deallocate( this% ImagePoints(i)% invPhi, &
                        this% ImagePoints(i)% b       )
         end do
         deallocate( this% ImagePoints ) 
      end if 

      deallocate( this% penalization,        &
                  this% stl,                 &
                  this% stlSurfaceIntegrals, &
                  this% stlMove,             &
                  this% Integral,            &
                  this% STLfilename          )

      if( this% ComputeBandRegion ) deallocate( this% BandRegion )

   end subroutine IBM_Destruct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  ------------------------------------------------------------
!  Creating the .omesh file for p-adaptation close to the STL  
!  ------------------------------------------------------------
   subroutine IBM_SetPolynomialOrder( this, elements, corners )
      implicit none
      
      !-arguments----------------------------------------------
      class(IBM_type),         intent(inout) :: this
      type(element),           intent(inout) :: elements(:)
      real(kind=RP), optional, intent(in)    :: corners(:,:) 
      !-local-variables----------------------------------------
      real(kind=RP) :: corner(NDIM)
      integer       :: STLNum, eID, k

      if( present(corners) ) then 
         do eID = 1, size(elements)
            do k = 1, 8 
               if( isInsideBox( elements(eID)% SurfInfo% corners(:,k), corners ) ) then
                  elements(eID)% Nxyz(1) = this% Nx 
                  elements(eID)% Nxyz(2) = this% Ny
                  elements(eID)% Nxyz(3) = this% Nz 
                  exit
               end if 
            end do 
         end do  
      else 
         do eID = 1, size(elements)
            !loop over the stl files
            do STLNum = 1, this% NumOfSTL
               loop: do k = 1, 8 
                  call OBB(STLNum)% ChangeRefFrame( elements(eID)% SurfInfo% corners(:,k), LOCAL, corner )
                  if( isInsideBox( corner, this% BandRegionCoeff*OBB(STLNum)% LocVertices ) ) then
                     if( this% Nx .ne. 0 ) elements(eID)% Nxyz(1) = this% Nx
                     if( this% Ny .ne. 0 ) elements(eID)% Nxyz(2) = this% Ny
                     if( this% Nz .ne. 0 ) elements(eID)% Nxyz(3) = this% Nz
                     exit loop
                  end if
               end do loop
            end do
         end do     
      end if 

   end subroutine IBM_SetPolynomialOrder
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -----------------------------------------------------------------------
!  Band region for the computation of the forces and turbulent quantities
!  -----------------------------------------------------------------------

   subroutine IBM_constructBandRegion( this, elements, no_of_elements, STLNum )
      use MPI_Process_Info
      use PhysicsStorage
      implicit none
      !-arguments------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: no_of_elements, &
                                        STLNum
      !-local-variables------------------------------------------------
      type(STLfile)             :: stl 
      type(point_type), pointer :: p 
      integer                   :: i
#ifdef _HAS_MPI_
      integer                   :: localVal, ierr
#endif        
      call this% BandRegionPoints( elements, no_of_elements, STLNum )

      this% BandRegion(STLNum)% LocNumOfObjs = this% BandPoints% NumOfPoints
#ifdef _HAS_MPI_
      call mpi_allreduce(this% BandPoints% NumOfPoints, this% BandRegion(STLNum)% NumOfObjs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
      this% BandRegion(STLNum)% NumOfObjs = this% BandPoints% NumOfPoints 
#endif

      if( this% BandRegion(STLNum)% NumOfObjs .eq. 0 ) then
         print *, "IBM_bandRegionPoints: Number of points in the band region is 0"
         error stop
      end if

      allocate(this% BandRegion(STLNum)% x(this% BandRegion(STLNum)% NumOfObjs))

      p => this% BandPoints% head

      do i = 1, this% BandPoints% NumOfPoints
         this% BandRegion(STLNum)% x(i)% index          = i
         this% BandRegion(STLNum)% x(i)% coords         = p% coords
         this% BandRegion(STLNum)% x(i)% local_Position = p% local_Position
         this% BandRegion(STLNum)% x(i)% element_index  = p% element_index
         this% BandRegion(STLNum)% x(i)% partition      = p% partition
         this% BandRegion(STLNum)% x(i)% STLNum         = STLNum
         this% BandRegion(STLNum)% x(i)% normal         = 0.0_RP
         this% BandRegion(STLNum)% x(i)% Dist           = huge(1.0_RP)
         p => p% next 
      end do

      call this% MPI_PointsListOperations( this% BandRegion(STLNum) )
  
      ! STL not used for rootpoints
      this% rootPoints(STLNum)% STLNum       = STLNum
      this% rootPoints(STLNum)% which_KDtree = POINTS_KDTREE

      call this% rootPoints(STLNum)% construct( stl           = this% stl(STLNum),                              &
                                                Vertices      = this% BandRegionCoeff*OBB(STLNum)% LocVertices, &
                                                isPlot        = this% plotBandPoints,                           &
                                                Min_n_of_Objs = this% NumOfInterPoints,                         &
                                                pointList     = this% BandRegion(STLNum)% x,                    &
                                                lvl           = this% lvl                                       )

      allocate( this% BandRegion(STLNum)% Q(NCONS,this% BandRegion(STLNum)% NumOfObjs) )

      call this% BandPoints% destruct()
 
   end subroutine IBM_constructBandRegion
!
!  Building band region for the computation of the STL-DoFs distance
!  ------------------------------------------------------------------
   subroutine IBM_constructBandRegion4Distance( this, elements, no_of_elements, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: no_of_elements, &
                                        STLNum
      !-local-variables------------------------------------
      type(point_type), pointer :: p 
      integer                   :: i
#ifdef _HAS_MPI_
      integer                   :: localVal, ierr
#endif
      call this% BandRegionPoints( elements, no_of_elements, STLNum )

      this% BandRegion4Distance(STLNum)% LocNumOfObjs = this% BandPoints% NumOfPoints
#ifdef _HAS_MPI_
      call mpi_allreduce(this% BandPoints% NumOfPoints, this% BandRegion4Distance(STLNum)% NumOfObjs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
      this% BandRegion4Distance(STLNum)% NumOfObjs  = this% BandPoints% NumOfPoints
#endif
      allocate(this% BandRegion4Distance(STLNum)% x(this% BandRegion4Distance(STLNum)% NumOfObjs))

      p => this% BandPoints% head

      do i = 1, this% BandPoints% NumOfPoints
         this% BandRegion4Distance(STLNum)% x(i)% index          = i
         this% BandRegion4Distance(STLNum)% x(i)% coords         = p% coords
         this% BandRegion4Distance(STLNum)% x(i)% local_Position = p% local_Position
         this% BandRegion4Distance(STLNum)% x(i)% element_index  = p% element_index
         this% BandRegion4Distance(STLNum)% x(i)% partition      = p% partition
         this% BandRegion4Distance(STLNum)% x(i)% STLNum         = STLNum
         this% BandRegion4Distance(STLNum)% x(i)% normal         = 0.0_RP
         this% BandRegion4Distance(STLNum)% x(i)% Dist           = huge(1.0_RP)
         p => p% next 
      end do

      call this% MPI_PointsListOperations( this% BandRegion4Distance(STLNum) )
        
      call this% BandPoints% destruct()

   end subroutine IBM_constructBandRegion4Distance
!
!  Band region points are found and stored
!  ---------------------------------------   
   subroutine IBM_bandRegionPoints( this, elements, n_of_elements, STLNum )
      use ElementClass
      use FaceClass
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: n_of_elements, &
                                        STLNum
      !-local-variables-----------------------------------
      type(point_type)  :: p          
      real(kind=RP)     :: Point(NDIM)     
      integer           :: eID, n, i, j, k, NumOfPoints

      this% BandPoints = PointLinkedList()      

      n = 0; NumOfPoints = 0

      do eID = 1, n_of_elements
         associate( e => elements(eID) )
         do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3) 
            if( .not. e% isInsideBody(i,j,k) .and. .not. e% isForcingPoint(i,j,k) ) then
               call OBB(STLNum)% ChangeRefFrame(e% geom% x(:,i,j,k),LOCAL,Point)
               if( isInsideBox(Point, this% BandRegionCoeff*OBB(STLNum)% LocVertices, .true. ) ) then
                  p% coords = Point
                  p% element_index = eID; p% local_Position = (/ i,j,k /)
                  n = n + 1; NumOfPoints = NumOfPoints + 1 
                  p% index = n; p% partition = MPI_Process% rank
                  call this% BandPoints% add(p)
               end if
            end if
         end do; end do; end do
         end associate 
      end do

      this% BandPoints% NumOfPoints = NumOfPoints

   end subroutine IBM_bandRegionPoints
!
!  Band region points' state is stored
!  -----------------------------------
   subroutine IBM_BandPoint_state( this, elements, STLNum, gradients )
      use PhysicsStorage
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------
      class(IBM_type),           intent(inout) :: this
      type(element),             intent(in)    :: elements(:)
      integer,                   intent(in)    :: STLnum
      logical,                   intent(in)    :: gradients
      !-local-variables------------------------------------------
      integer                    :: n, i, j, k, eID 
#ifdef _HAS_MPI_
      real(kind=rp), allocatable :: local_sum(:), global_sum(:)
      integer                    :: ierr, NumOfObjs
#endif
      this% BandRegion(STLNum)% Q = 0.0_RP
      if( gradients ) then 
         this% BandRegion(STLNum)% U_x = 0.0_RP
         this% BandRegion(STLNum)% U_y = 0.0_RP
         this% BandRegion(STLNum)% U_z = 0.0_RP
      end if
!$omp parallel
!$omp do schedule(runtime) private(i,j,k,eID)
      do n = 1, this% BandRegion(STLNum)% NumOfObjs
         if( this% BandRegion(STLNum)% x(n)% partition .eq. MPI_Process% rank ) then
            i   = this% BandRegion(STLNum)% x(n)% local_Position(1)
            j   = this% BandRegion(STLNum)% x(n)% local_Position(2)
            k   = this% BandRegion(STLNum)% x(n)% local_Position(3)
            eID = this% BandRegion(STLNum)% x(n)% element_index
            this% BandRegion(STLNum)% Q(:,n) = elements(eID)% storage% Q(:,i,j,k)   
            if( gradients ) then
               this% BandRegion(STLNum)% U_x(:,n) = elements(eID)% storage% U_x(:,i,j,k)
               this% BandRegion(STLNum)% U_y(:,n) = elements(eID)% storage% U_y(:,i,j,k)
               this% BandRegion(STLNum)% U_z(:,n) = elements(eID)% storage% U_z(:,i,j,k)
            end if
         end if
      end do
!$omp end do
!$omp end parallel
#ifdef _HAS_MPI_
      if( MPI_Process% doMPIAction ) then     
         NumOfObjs = this% BandRegion(STLNum)% NumOfObjs
         allocate(local_sum(NumOfObjs),global_sum(NumOfObjs))
         do i = 1, NCONS
            local_sum = this% BandRegion(STLNum)% Q(i,:)
            call mpi_allreduce(local_sum, global_sum, NumOfObjs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr) 
            this% BandRegion(STLNum)% Q(i,:) = global_sum  
            if( gradients ) then
               local_sum = this% BandRegion(STLNum)% U_x(i,:)
               call mpi_allreduce(local_sum, global_sum, NumOfObjs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
               this% BandRegion(STLNum)% U_x(i,:) = global_sum
               local_sum = this% BandRegion(STLNum)% U_y(i,:)
               call mpi_allreduce(local_sum, global_sum, NumOfObjs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
               this% BandRegion(STLNum)% U_y(i,:) = global_sum
               local_sum = this% BandRegion(STLNum)% U_z(i,:)
               call mpi_allreduce(local_sum, global_sum, NumOfObjs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
               this% BandRegion(STLNum)% U_z(i,:) = global_sum  
            end if
         end do 
         deallocate(local_sum,global_sum)
      end if
#endif
   end subroutine IBM_BandPoint_state
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  ------------------------------------------------------------
!  MPI routines  
!  ------------------------------------------------------------
   subroutine IBM_MPI_sendOBB( this, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum
      
      if ( MPI_Process% doMPIAction ) then   
         call recvOBB( STLNum )
      end if

      if( MPI_Process% doMPIRootAction ) then
         call SendOBB( STLNum )
      end if
   
   end subroutine IBM_MPI_sendOBB
  
   subroutine IBM_MPI_sendSTLpartitions( this, STLNum, vertices  )
      use MPI_Process_Info
      implicit none
      !-arguments---------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum
      real(kind=RP),   intent(inout) :: vertices(:,:)
      
      if( MPI_Process% doMPIAction ) then 
         call receiveSTLpartitions( this% stl(STLNum), STLNum, vertices, this% root(STLNum)% MaxAxis )
      end if 

      if( MPI_Process% isRoot ) then
         call sendSTL2Partitions( this% stl(STLNum), STLNum, vertices, this% root(STLNum)% MaxAxis )
      end if

      this% stl(STLNum)% construct = .true.
  
   end subroutine IBM_MPI_sendSTLpartitions

   subroutine IBM_MPI_PointsListOperations( this, pointsList )
      use MPI_Process_Info
      implicit none 
      !-arguments----------------------------
      class(IBM_type),           intent(inout) :: this 
      type(IBMpoints),           intent(inout) :: pointsList

      if( MPI_Process% doMPIAction ) then
         call SendPointsList2Root( pointsList )
      end if

      if( MPI_Process% doMPIRootAction ) then 
         call RecvPointsListRoot( pointsList )
      end if

      if( MPI_Process% doMPIRootAction ) then 
         call SendPointsList2partitions( pointsList )
      end if

      if( MPI_Process% doMPIAction ) then
         call RecvPointsListpartitions( pointsList )
      end if

   end subroutine IBM_MPI_PointsListOperations
!
!  MASK: root sends and receives its mask points 
!  ----------------------------------------------
   subroutine IBM_MPI_sendMask2Root( this )
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------
      class(IBM_type), intent(inout) :: this
   
      if( MPI_Process% isRoot ) then 
         call recvPointsMaskRoot()
      end if 

      if( MPI_Process% doMPIAction ) then 
         call sendPointsMask2Root()
      endif  
      
   end subroutine IBM_MPI_sendMask2Root
!
!  MASK: all points go from root to partitions
!  -------------------------------------------
   subroutine IBM_MPI_sendMask2Partitions( this )
      use MPI_Process_Info
      implicit none
      !-arguments----------------------------
      class(IBM_type), intent(inout) :: this
 
      if( MPI_Process% doMPIAction ) then 
         call recvPointsMaskPartitions()
      end if 

      if( MPI_Process% doMPIRootAction ) then 
         call sendPointsMask2Partitions()
      endif

   end subroutine IBM_MPI_sendMask2Partitions

   subroutine IBM_MPI_sendNormals2Root( this, pointsList, ranks )
      use MPI_Process_Info
      implicit none 
      !-arguments----------------------------
      class(IBM_type), intent(inout) :: this 
      type(IBMpoints), intent(inout) :: pointsList
      real(kind=RP),   intent(in)    :: ranks(:)
      
      if( MPI_Process% doMPIAction ) then
         call sendNormals2Root( pointsList )
      end if
      if( MPI_Process% doMPIRootAction ) then
         call recvNormalsRoot( PointsList, ranks )
     end if
   
   end subroutine IBM_MPI_sendNormals2Root
        
   subroutine IBM_MPI_sendDistNormals2partitions( this, PointsList )
      use MPI_Process_Info
      implicit none
      !-arguments-------------------------------
      class(IBM_type), intent(inout) :: this
      type(IBMpoints), intent(inout) :: pointsList
     
      if( MPI_Process% doMPIRootAction ) then
         call sendDistanceANDNormals2partitions( PointsList )
      end if
      if( MPI_Process% doMPIAction ) then
         call recvDistancesANDNormalspartitions( PointsList )
      end if
   
   end subroutine IBM_MPI_sendDistNormals2partitions
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  ------------------------------------------------------------
!  Moving bodies  
!  ------------------------------------------------------------      
   subroutine IBM_MoveBody( this, elements, no_of_elements, no_of_DoFs, isChild, t, iter, autosave )
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------------------------------------------
      class(IBM_type),           intent(inout) :: this
      type(element),             intent(inout) :: elements(:)
      integer,                   intent(in)    :: no_of_elements, no_of_DoFs
      logical,                   intent(in)    :: isChild 
      real(kind=RP),             intent(in)    :: t
      integer,         optional, intent(in)    :: iter
      logical,         optional, intent(in)    :: autosave
      !-local-variables-----------------------------------------------------------
      integer :: STLNum
 
      do STLNum = 1, this% NumOfSTL
         if( this% stl(STLNum)% move ) then
            call this% root(STLNum)% Destruct( isChild )
            if( this% ComputeBandRegion ) then
               deallocate( this% BandRegion(STLNum)% x )
               call this% rootPoints(STLNum)% destruct( isChild )
            end if
            if( this% ComputeDistance ) call this% rootDistance(STLNum)% destruct( isChild )
            call this% CleanMask( elements, no_of_elements, STLNum ) 
            if( MPI_Process% isRoot .and. .not. isChild ) then
               call this% stl(STLNum)% ReadTessellation( this% STLfilename(STLNum) )
               if( this% stl(STLNum)% motionType .eq. ROTATION ) then
                  call this% stl(STLNum)% getRotationaMatrix( t )
                  call OBB(STLNum)% STL_rotate( this% stl(STLNum), .false. )
               elseif( this% stl(STLNum)% motionType .eq. LINEAR ) then
                  call this% stl(STLNum)% getDisplacement( t )
                  call OBB(STLNum)% STL_translate( this% stl(STLNum), .false. )
               end if
               call this% stl(STLNum)% updateNormals()
               this% stl(STLNum)% show = .false.
               this% plotMask          = .false.
               if( present(autosave) .and. autosave ) this% plotMask = .true.
               if( present(autosave) .and. autosave ) call this% stl(STLNum)% plot( iter )
               call OBB(STLNum)% construct( this% stl(STLNum), this% plotOBB, this% AAB )
               call OBB(STLNum)% ChangeObjsRefFrame( this% stl(STLNum)% ObjectsList, LOCAL ) 
            end if
            call this% MPI_sendOBB(STLNum) 
            this% plotKDtree = .false.
            call this% constructSTL_KDtree( STLNum ) 
         end if
      end do
      
      do STLNum = 1, this% NumOfSTL
         if( this% stl(STLNum)% move ) call this% build( elements, no_of_elements, no_of_DoFs, isChild, STLNum, iter )
      end do

   end subroutine IBM_MoveBody
!
!  Mask is cleaned
!  ---------------
   subroutine IBM_CleanMask( this, elements, no_of_elements, STLNum )
   
      implicit none
      !-arguments------------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: no_of_elements, STLNum
      !-local-variables------------------------------------------
      integer :: eID, i, j, k
!$omp parallel
!$omp do schedule(runtime) private(i,j,k)
      do eID = 1, no_of_elements
         do i = 0, elements(eID)% Nxyz(1); do j = 0, elements(eID)% Nxyz(2); do k = 0, elements(eID)% Nxyz(3) 
            if( elements(eID)% STL(i,j,k) .eq. STLNum .and. elements(eID)% isInsideBody(i,j,k) ) then
               elements(eID)% STL(i,j,k) = 0
               elements(eID)% isInsideBody(i,j,k) = .false.
            end if
         end do; end do; end do
      end do 
!$omp end do
!$omp end parallel
   end subroutine IBM_CleanMask
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  ----------------------------------------------------
!  Geometrical quantities for the turbulence wall model  
!  ----------------------------------------------------   
   subroutine IBM_GetForcingPointsGeom( this, elements )
      use MPI_Process_Info
      use MeshTypes
      implicit none
      !-arguments----------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      !-local-variables----------------------------------
      real(kind=RP) :: dx, dy, dz, d, d_min, Dist
      integer       :: eID, i, j, k, n, Total
#ifdef _HAS_MPI_
      integer       :: ierr
#endif
      this% NumOfForcingPoints = 0
      d = huge(1.0_RP) 

      do eID = 1, size(elements) 
         dx = maxval(elements(eID)% SurfInfo% corners(1,:)) - minval(elements(eID)% SurfInfo% corners(1,:))
         dy = maxval(elements(eID)% SurfInfo% corners(2,:)) - minval(elements(eID)% SurfInfo% corners(2,:))
         dz = maxval(elements(eID)% SurfInfo% corners(3,:)) - minval(elements(eID)% SurfInfo% corners(3,:))
         d  = min(d,min(abs(dx),abs(dy),abs(dz)))   
      end do    
      d_min = sqrt(3.0_RP)*d 
#ifdef _HAS_MPI_
      call mpi_allreduce(d_min, this% IP_Distance, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
      Dist = this% IP_Distance
      this% IP_Distance = 1.5_RP*this% IP_Distance
#else
      this% IP_Distance = d_min
      Dist = this% IP_Distance
      this% IP_Distance = 1.5_RP*this% IP_Distance
#endif        
      do eID = 1, size(elements) 
         do k = 0, elements(eID)% Nxyz(3); do j = 0, elements(eID)% Nxyz(2); do i = 0, elements(eID)% Nxyz(1)   
            if( elements(eID)% geom% dWall(i,j,k) .gt. Dist ) cycle
            elements(eID)% isForcingPoint(i,j,k) = .true.
            this% NumOfForcingPoints = this% NumOfForcingPoints + 1
         end do; end do; end do     
      end do  
#ifdef _HAS_MPI_
      call mpi_allreduce(this% NumOfForcingPoints, Total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
      Total = this% NumOfForcingPoints
#endif
      if( Total .eq. 0 .and. MPI_Process% isRoot ) then 
         print*, "No forcing points found, change y+ target"
         error stop 
      endif

      allocate(this% ImagePoints(this% NumOfForcingPoints))
      n = 0
 
      do eID = 1, size(elements) 
         do k = 0, elements(eID)% Nxyz(3); do j = 0, elements(eID)% Nxyz(2); do i = 0, elements(eID)% Nxyz(1) 
            if( elements(eID)% isForcingPoint(i,j,k) ) then 
               Dist = this% IP_distance - elements(eID)% geom% dWall(i,j,k)
               n = n + 1
               this% ImagePoints(n)% coords = elements(eID)% geom% x(:,i,j,k) + Dist * elements(eID)% geom% normal(:,i,j,k)
               this% ImagePoints(n)% element_index = eID 
               this% ImagePoints(n)% local_Position = (/i,j,k/)
               this% ImagePoints(n)% index = n
            end if 
         end do; end do; end do 
      end do

      call Plot_Forcing_Imagepoints( this, elements )

   end subroutine IBM_GetForcingPointsGeom
      
   subroutine Plot_Forcing_Imagepoints( IBM, elements )
      use MPI_Process_Info
      implicit none
      !-arguments----------------------------------------
      type(IBM_type), intent(in) :: IBM
      type(element),  intent(in) :: elements(:)
      !-local-variables----------------------------------
      real(kind=RP)              :: ImagePoint(NDIM), &
                                    Dist
      character(len=LINE_LENGTH) :: filenameFP,       &
                                    filenameIP,       &
                                    lvl, rank
      integer                    :: eID, i, j, k,     &
                                    funit
      logical                    :: add_Point = .false.
      
      if( IBM% NumOfForcingPoints .eq. 0 ) return
 
      if( MPI_Process% nProcs .gt. 1 ) then 
         write(rank,*) MPI_Process% rank
         if( IBM% lvl .gt. 0 ) then 
            write(lvl,*)  IBM% lvl
            filenameFP = 'ForcingPoints_'//trim(IBM% filename)//'_MGlevel'//trim(adjustl(lvl))//'_Process'//trim(adjustl(rank))
            filenameIP = 'ImagePoints_'//trim(IBM% filename)//'_MGlevel'//trim(adjustl(lvl))//'_Process'//trim(adjustl(rank))
         else
            filenameFP = 'ForcingPoints_'//trim(IBM% filename)//'_Process'//trim(adjustl(rank))
            filenameIP = 'ImagePoints_'//trim(IBM% filename)//'_Process'//trim(adjustl(rank))        
         end if    
      else
         if( IBM% lvl .gt. 0 ) then 
            write(lvl,*)  IBM% lvl
            filenameFP = 'ForcingPoints_'//trim(IBM% filename)//'_MGlevel'//trim(adjustl(lvl))
            filenameIP = 'ImagePoints_'//trim(IBM% filename)//'_MGlevel'//trim(adjustl(lvl))
         else
            filenameFP = 'ForcingPoints_'//trim(IBM% filename)
            filenameIP = 'ImagePoints_'//trim(IBM% filename)       
         end if
      endif

      call TecFileHeader( 'IBM/'//trim(filenameFP), 'Forcing Points', IBM% NumOfForcingPoints/2+mod(IBM% NumOfForcingPoints,2),2,1, funit, 'POINT' )

      if( mod(IBM% NumOfForcingPoints,2) .ne. 0 ) add_Point  = .true.

      do eID = 1, size(elements)
         associate ( e => elements(eID) )
         do k = 0, e% Nxyz(3); do j = 0, e% Nxyz(2) ; do i = 0, e% Nxyz(1)
            if( e% isForcingPoint(i,j,k) ) then
               write(funit,'(3E13.5)')  e% geom% x(1,i,j,k), e% geom% x(2,i,j,k), e% geom% x(3,i,j,k)
               if( add_Point ) then
                  write(funit,'(3E13.5)') e% geom% x(1,i,j,k), e% geom% x(2,i,j,k), e% geom% x(3,i,j,k)
                  add_Point = .false.
               end if
            end if
         end do; end do; end do
         end associate
      end do
         
      close(funit)

      call TecFileHeader( 'IBM/'//trim(filenameIP), 'Image Points', IBM% NumOfForcingPoints/2+mod(IBM% NumOfForcingPoints,2),2,1, funit, 'POINT' )      

      do i = 1, IBM% NumOfForcingPoints
         write(funit,'(3E13.5)')  IBM% ImagePoints(i)% coords(1), IBM% ImagePoints(i)% coords(2), IBM% ImagePoints(i)% coords(3)
         if( add_Point ) then
            write(funit,'(3E13.5)') IBM% ImagePoints(i)% coords(1), IBM% ImagePoints(i)% coords(2), IBM% ImagePoints(i)% coords(3)
            add_Point = .false.
         end if
      end do
         
      close(funit)   

   end subroutine Plot_Forcing_Imagepoints
   
!
!  DoFs closer to each image point
!  -------------------------------
   subroutine IBM_GetImagePoint_nearest( this, elements )
      use MappedGeometryClass
      use PhysicsStorage
      use MPI_Process_Info
      use ElementClass
      implicit none
      !-arguments-------------------------------------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      !-local-variables-------------------------------------------------------------------------------
      real(kind=RP) :: iP(NDIM), Dist 
      integer       :: i, n, STLNum  
!$omp parallel
!$omp do schedule(runtime) private(Dist,STLNum,iP,n)
         do i = 1, this% NumOfForcingPoints 
            STLNum = elements(this% ImagePoints(i)% element_index)% STL( this% ImagePoints(i)% local_Position(1), &
                                                                         this% ImagePoints(i)% local_Position(2), &
                                                                         this% ImagePoints(i)% local_Position(3)  )

            call OBB(STLNum)% ChangeRefFrame(this% ImagePoints(i)% coords,LOCAL,iP)

            allocate( this% ImagePoints(i)% nearestPoints(this% NumOfInterPoints),                 &
                      this% ImagePoints(i)% invPhi(this% NumOfInterPoints,this% NumOfInterPoints), &
                      this% ImagePoints(i)% b(this% NumOfInterPoints)                              )

            this% ImagePoints(i)% nearestPoints = 0

            do n = 1, this% NumOfInterPoints 
               call MinimumDistancePoints( iP, this% rootPoints(STLNum), this% BandRegion(STLNum), &
                                           Dist, n, this% ImagePoints(i)% nearestPoints            ) 
            end do
                  
            if(any(this% ImagePoints(i)% nearestPoints .eq. 0) ) then
               print *, "Can't fine nearest point for the image point: ", i
               error stop
            end if

            call GetMatrixInterpolationSystem( iP,                                                               &
                                               this% BandRegion(STLNum)% x(this% ImagePoints(i)% nearestPoints), &
                                               this% ImagePoints(i)% invPhi,                                     &
                                               this% ImagePoints(i)% b, this% InterpolationType                  )
         end do
!$omp end do
!$omp end parallel
   end subroutine IBM_GetImagePoint_nearest
!
!  Distance from the wall needed 
!  -----------------------------
   subroutine IBM_ComputeIBMWallDistance( this, elements, movingSTL )
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------------------
      class(IBM_type),    intent(inout) :: this
      type(element),      intent(inout) :: elements(:)
      integer, optional , intent(in)    :: movingSTL
      !-local-variables--------------------------------
      integer :: STLNum

      allocate(this% BandRegion4Distance(this% NumOfSTL))

      do STLNum = 1, this% NumOfSTL
         if( .not. present(movingSTL) )then
            call this% GetDistanceInsideBox( elements, STLNum )
         else
            if( STLNum .eq. movingSTL) call this% GetDistanceInsideBox( elements, STLNum )
         end if
      end do

      do STLNum = 1, this% NumOfSTL
         if( .not. present(movingSTL) )then
            call this% GetDistanceOutsideBox( STLNum )
         else
            if( STLNum .eq. movingSTL) call this% GetDistanceOutsideBox( STLNum )
         end if
      end do

      do STLNum = 1, this% NumOfSTL
         if( .not. present(movingSTL) )then
            call SetDistances( this% BandRegion4Distance(STLNum), elements )
         else
            if( STLNum .eq. movingSTL) call SetDistances( this% BandRegion4Distance(STLNum), elements )
         endif
         deallocate( this% BandRegion4Distance(STLNum)% x)
      end do

      deallocate(this% BandRegion4Distance)
 
      if( this% Wallfunction ) then
#if defined(NAVIERSTOKES)
         call this% GetForcingPointsGeom( elements )
#endif
      endif

   end subroutine IBM_ComputeIBMWallDistance
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------------------------------------------
!  Computing the distance between a point and the triangles inside the box it lies in 
!  -------------------------------------------------------------------------------------    
   subroutine IBM_GetDistanceInsideBox( this, elements, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: STLNum
      !-local-variables------------------------------
      integer       :: i, no_of_elements
      real(kind=RP) :: xP(NDIM), Dist, normal(NDIM)

      no_of_elements = size(elements)

      call this% constructBandRegion4Distance( elements, no_of_elements, STLNum )
!$omp parallel  
!$omp do schedule(runtime) private(xP,Dist,normal)
      do i = 1, this% BandRegion4Distance(STLNum)% NumOfObjs
         this% BandRegion4Distance(STLNum)% x(i)% Dist = huge(1.0_RP)   
         this% BandRegion4Distance(STLNum)% x(i)% STLNum = STLNum      
         this% BandRegion4Distance(STLNum)% x(i)% normal = 0.0_RP      
         xP = this% BandRegion4Distance(STLNum)% x(i)% coords                            
         if( isInsideBox(xP,this% rootDistance(STLNum)% vertices) ) then                              
            call MinimumDistance( xP, this% rootDistance(STLNum), Dist, normal )   
            if( Dist .lt. this% BandRegion4Distance(STLNum)% x(i)% Dist ) then 
               this% BandRegion4Distance(STLNum)% x(i)% Dist   = Dist
               this% BandRegion4Distance(STLNum)% x(i)% normal = normal  
            end if
         end if      
      end do   
!$omp end do
!$omp end parallel    
   end subroutine IBM_GetDistanceInsideBox
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------------------------------------------------
!  Computing the distance between a point and the triangles inside the box it does not lay in 
!  -------------------------------------------------------------------------------------------       
   subroutine IBM_GetDistanceOutsideBox( this, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum
      !-local-variables-----------------------------
      integer       :: i, no_of_elements
      real(kind=RP) :: xP(NDIM), Dist, normal(NDIM)
#ifdef _HAS_MPI_ 
      real(kind=rp), allocatable :: input(:,:), output(:,:)
      integer                    :: ierr, rank
#endif
#ifdef _HAS_MPI_
      allocate( input(2,this% BandRegion4Distance(STLNum)% NumOfObjs), &
                output(2,this% BandRegion4Distance(STLNum)% NumOfObjs) )
#endif     
!$omp parallel  
!$omp do schedule(runtime) private(xP,Dist,normal)
      do i = 1, this% BandRegion4Distance(STLNum)% NumOfObjs         
         xP = this% BandRegion4Distance(STLNum)% x(i)% coords        
         this% BandRegion4Distance(STLNum)% x(i)% STLNum = STLNum                   
         if( .not. isInsideBox(xP,this% rootDistance(STLNum)% vertices) ) then 
           call MinimumDistanceSphereKDtree( xP, this% rootDistance(STLNum), this% BandRegion4Distance(STLNum)% x(i)% Dist, Dist, normal )        
            if( Dist .lt. this% BandRegion4Distance(STLNum)% x(i)% Dist ) then 
               this% BandRegion4Distance(STLNum)% x(i)% Dist   = Dist
               this% BandRegion4Distance(STLNum)% x(i)% normal = normal      
            end if
         end if        
#ifdef _HAS_MPI_ 
         input(1,i) = this% BandRegion4Distance(STLNum)% x(i)% Dist
         input(2,i) = MPI_Process% rank
#endif    
      end do 
!$omp end do
!$omp end parallel
#ifdef _HAS_MPI_ 
      call mpi_allreduce( input, output, this% BandRegion4Distance(STLNum)% NumOfObjs, &
                          MPI_2DOUBLE_PRECISION, MPI_MINLOC,                           &
                          MPI_COMM_WORLD, ierr                                         )
       
      if( MPI_Process% isRoot ) then
         do i = 1, this% BandRegion4Distance(STLNum)% NumOfObjs
            this% BandRegion4Distance(STLNum)% x(i)% Dist = output(1,i)
         end do
      end if

      call this% MPI_sendNormals2Root( this% BandRegion4Distance(STLNum), output(2,:) )
      call this% MPI_sendDistNormals2partitions( this% BandRegion4Distance(STLNum)  )    
 
      deallocate(input,output)
#endif       
   end subroutine IBM_GetDistanceOutsideBox
   
   subroutine SetDistances( BandRegion4Distance, elements )
      use MPI_Process_Info
      implicit none
      !-arguments---------------------------------------------------
      type(IBMpoints), intent(inout) :: BandRegion4Distance
      type(element),   intent(inout) :: elements(:)
      !-local-variables---------------------------------------------
      integer :: i
!$omp parallel
!$omp do schedule(runtime)
      do i = 1, BandRegion4Distance% NumOfObjs
         BandRegion4Distance% x(i)% Dist = sqrt(BandRegion4Distance% x(i)% Dist)
         if( BandRegion4Distance% x(i)% partition .eq. MPI_Process% rank ) then           
            elements(BandRegion4Distance% x(i)% element_index)% geom%  &
               dWall(BandRegion4Distance% x(i)% local_Position(1),     &
                     BandRegion4Distance% x(i)% local_Position(2),     &
                     BandRegion4Distance% x(i)% local_Position(3)  ) = BandRegion4Distance% x(i)% Dist
            elements(BandRegion4Distance% x(i)% element_index)% geom%        &
                     normal(:, BandRegion4Distance% x(i)% local_Position(1), &
                               BandRegion4Distance% x(i)% local_Position(2), &
                               BandRegion4Distance% x(i)% local_Position(3)  ) = BandRegion4Distance% x(i)% normal
            elements(BandRegion4Distance% x(i)% element_index)%                &
                             STL(BandRegion4Distance% x(i)% local_Position(1), &
                                 BandRegion4Distance% x(i)% local_Position(2), &
                                 BandRegion4Distance% x(i)% local_Position(3)  ) = BandRegion4Distance% x(i)% STLNum
         end if
      end do
!$omp end do
!$omp end parallel
   end subroutine SetDistances
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  Source terms for the immersed boundary 
!  ------------------------------------------------   
   subroutine IBM_SourceTerm( this, eID, Q, Q_target, Source, Wallfunction )
      use PhysicsStorage
      implicit none
      !-arguments--------------------------------------------
      class(IBM_type),           intent(inout) :: this
      integer,                   intent(in)    :: eID
      real(kind=rp),             intent(in)    :: Q(:)
      real(kind=rp),   optional, intent(in)    :: Q_target(:)
      real(kind=rp),             intent(inout) :: Source(:)
      logical,                   intent(in)    :: Wallfunction
      !-local-variables--------------------------------------
      real(kind=rp) :: rho, rho_s, u, u_s, v, v_s, w, w_s
#if defined(SPALARTALMARAS)
      real(kind=rp) :: theta, theta_s
#endif
      Source = 0.0_RP
#if defined(NAVIERSTOKES)
      if( present(Q_target) ) then
         rho_s = Q_target(IRHO)
         u_s   = Q_target(IRHOU)/rho_s
         v_s   = Q_target(IRHOV)/rho_s
         w_s   = Q_target(IRHOW)/rho_s   
#if defined(SPALARTALMARAS)
         theta_s = Q_target(IRHOTHETA)/rho_s  
#endif
      else
         rho_s = Q(IRHO)
         u_s   = 0.0_RP
         v_s   = 0.0_RP
         w_s   = 0.0_RP
#if defined(SPALARTALMARAS)
         theta_s = 0.0_RP  
#endif
      end if   
        
      rho = Q(IRHO)
      u   = Q(IRHOU)/rho
      v   = Q(IRHOV)/rho
      w   = Q(IRHOW)/rho
      
      Source(IRHOU) = rho*u-rho_s*u_s  
      Source(IRHOV) = rho*v-rho_s*v_s 
      Source(IRHOW) = rho*w-rho_s*w_s
      Source(IRHOE) = 0.5_RP * rho*( POW2(u) + POW2(v) + POW2(w) )  & 
                     -0.5_RP * rho_s*( POW2(u_s) + POW2(v_s) + POW2(w_s) )
#if defined(SPALARTALMARAS)
      theta = Q(IRHOTHETA)/rho
      Source(IRHOTHETA) = rho*theta - rho_s*theta_s
#endif                     
      if( Wallfunction ) then 
         Source = -1.0_RP/(this% penalCoeff * this% penalization(eID)) * Source 
      else 
         Source = -1.0_RP/this% penalization(eID) * Source
      end if 
#endif   
#if defined(INCNS)
      rho_s = Q(INSRHO)
      u_s   = 0.0_RP
      v_s   = 0.0_RP
      w_s   = 0.0_RP

      rho = Q(INSRHO)
      u   = Q(INSRHOU)/rho
      v   = Q(INSRHOV)/rho
      w   = Q(INSRHOW)/rho
      
      ! no penalty for pressure
      Source(INSRHOU) = rho*u-rho_s*u_s  
      Source(INSRHOV) = rho*v-rho_s*v_s 
      Source(INSRHOW) = rho*w-rho_s*w_s

      Source = -1.0_RP/this% penalization(eID) * Source
#endif   
!
#if defined(MULTIPHASE)
      ! this is not u, is sqrt(rho)*, we are not computing rho for now to save time as is not really use without velocity
      u_s   = 0.0_RP
      v_s   = 0.0_RP
      w_s   = 0.0_RP

      ! no penalty for pressure
      ! the source for MULTIPHASE is sqrt(rho)*u - sqrt(rho_s)*u_s

      Source(IMSQRHOU) = Q(IMSQRHOU) - u_s
      Source(IMSQRHOV) = Q(IMSQRHOV) - v_s
      Source(IMSQRHOW) = Q(IMSQRHOW) - w_s

      ! no extra scaling for multiphase source term as is already in sqrt(rho)*u and no rho*u
      Source = -1.0_RP/this% penalization(eID) * Source
#endif   
   end subroutine IBM_SourceTerm

   function IBM_MaskVelocity( this, Q, nEqn, STLNum, x, t ) result( Q_target )
      use PhysicsStorage
      use FluidData
      use TessellationTypes
      use OrientedBoundingBox
      implicit none

      class(IBM_type), intent(in) :: this 
      integer,         intent(in)    :: nEqn, STLNum
      real(kind=RP),   intent(in)    :: Q(nEqn), x(NDIM), t
      real(kind=RP)                  :: Q_target(nEqn)

      real(kind=RP) :: R, v_plane, time, theta, &
                       rho, u_s, v_s, w_s
      
      Q_target = Q

      u_s = 0.0_RP; v_s = 0.0_RP; w_s = 0.0_RP
#if defined(NAVIERSTOKES)
      if( this% stl(STLNum)% motionType .eq. LINEAR ) then
         select case( this% stl(STLNum)% motionAxis )
         case( IX )
            u_s = this% stl(STLNum)% Velocity/refValues% V
         case( IY )
            v_s = this% stl(STLNum)% Velocity/refValues% V
         case( IZ ) 
            w_s = this% stl(STLNum)% Velocity/refValues% V
         end select 
      elseif( this% stl(STLNum)% motionType .eq. ROTATION ) then
         R = norm2( x - this% stl(STLNum)% rotationCenter ) 

         v_plane = this% stl(STLNum)% angularVelocity * R * Lref
         v_plane = v_plane/refValues% V

         time    = t * Lref/refValues% V
         theta   = this% stl(STLNum)% angularVelocity * time

         select case( this% stl(STLNum)% motionAxis )
         case( IX ) 
            v_s = -sin(theta) * v_plane 
            w_s =  cos(theta) * v_plane 
         case( IY )
            u_s =  cos(theta) * v_plane  
            w_s = -sin(theta) * v_plane 
         case( IZ )
            u_s = -sin(theta) * v_plane
            v_s =  cos(theta) * v_plane 
         end select
      end if 

      rho = Q_target(IRHO)

      Q_target(IRHOU) = rho * u_s      
      Q_target(IRHOV) = rho * v_s      
      Q_target(IRHOW) = rho * w_s      
#endif
   end function IBM_MaskVelocity
!
!  Analytical jacobian: dS/dQ 
!  ---------------------------    
   subroutine IBM_semiImplicitJacobian( this, eID, Q, dS_dQ )
      use PhysicsStorage
      implicit none
      !-arguments----------------------------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: eID
      real(kind=rp),   intent(in)    :: Q(:)
      real(kind=rp),   intent(inout) :: dS_dQ(:,:)
      !-local-variables----------------------------------------------
      real(kind=rp) :: rho, u, v, w

      dS_dQ = 0.0_RP       
#if defined(NAVIERSTOKES) 
      rho = Q(IRHO)
      u   = Q(IRHOU)/rho 
      v   = Q(IRHOV)/rho 
      w   = Q(IRHOW)/rho 

      dS_dQ(IRHOU,IRHOU) = 1.0_RP
      dS_dQ(IRHOV,IRHOV) = 1.0_RP
      dS_dQ(IRHOW,IRHOW) = 1.0_RP
      dS_dQ(IRHOE,IRHO)  = -0.5_RP*( POW2(u) + POW2(v) + POW2(w) )
      dS_dQ(IRHOE,IRHOU) = u 
      dS_dQ(IRHOE,IRHOV) = v 
      dS_dQ(IRHOE,IRHOW) = w
#if defined(SPALARTALMARAS)
       dS_dQ(IRHOTHETA,IRHOTHETA) = 1.0_RP
#endif
#endif
       dS_dQ = -1.0_RP/this% penalization(eID) * dS_dQ
 
    end subroutine IBM_semiImplicitJacobian

   subroutine IBM_semiImplicitTurbulenceJacobian( this, ImagePoint, Q, normal, dWall, STLNum, dS_dQ )
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------------------------------------
      class(IBM_type),  intent(inout) :: this
      type(point_type), intent(in)    :: ImagePoint  
      real(kind=rp),    intent(in)    :: Q(:), dWall, normal(:)
      integer,          intent(in)    :: STLNum
      real(kind=rp),    intent(inout) :: dS_dQ(:,:)
      !-local-variables-----------------------------------------------------------
      real(kind=rp) :: Q_IP(NCONS), Q_FP(NCONS), Qnext(NCONS), &
                       Q_IPnext(NCONS), Q_FPnext(NCONS),       &
                       dS_dQns(NCONS,NCONS), iP(NDIM), eps 
      integer       :: i

      do i = 1, NCONS 
         Q_IP(i) = GetInterpolatedValue( this% BandRegion(STLNum)% Q(i,ImagePoint% nearestPoints), &
                                         ImagePoint% invPhi,                                       &
                                         ImagePoint% b,                                            &
                                         this% InterpolationType                                   ) 
      end do

      Q_FP = Q
#if defined(NAVIERSTOKES) 
      call ForcingPointState( Q_IP, this% IP_Distance, dWall, normal, Q_FP )
#endif
      eps = sqrt(EPSILON(eps))

      do i = 1, NCONS
         Q_FPnext = Q; Q_IPnext = Q_IP       
         Q_IPnext(i) = Q_IP(i) + eps 
#if defined(NAVIERSTOKES)
         call ForcingPointState( Q_IPnext, this% IP_Distance, dWall, normal, Q_FPnext )
#endif
         dS_dQ(:,i) = (Q_FPnext - Q_FP)/eps 
      end do 

      dS_dQ = -1.0_RP/(this% penalCoeff*this% penalization(ImagePoint% element_index)) * dS_dQ 

      call this% semiImplicitJacobian( ImagePoint% element_index, Q, dS_dQns )

      dS_dQ = dS_dQns - dS_dQ

   end subroutine IBM_semiImplicitTurbulenceJacobian
!
!  Analytical inverse jacobian: (1 - dt*dS/dQ)^(-1) 
!  --------------------------------------------------   
   subroutine IBM_semiImplicitShiftJacobian( this, eID, Q, dt, invdS_dQ )
      use PhysicsStorage
      implicit none
         !-arguments----------------------------------------------------------------
         class(IBM_type), intent(inout) :: this
         integer,         intent(in)    :: eID
         real(kind=rp),   intent(in)    :: Q(:)
         real(kind=rp),   intent(in)    :: dt
         real(kind=rp),   intent(inout) :: invdS_dQ(:,:)
         !-local-variables----------------------------------------------------------
         real(kind=rp) :: rho, u, v, w

         invdS_dQ = 0.0_RP  
         
         associate( eta => this% penalization(eID) )     
#if defined(NAVIERSTOKES) 
         rho = Q(IRHO)
         u   = Q(IRHOU)/rho 
         v   = Q(IRHOV)/rho 
         w   = Q(IRHOW)/rho 
         
         invdS_dQ(IRHO,IRHO)   = 1.0_RP
         invdS_dQ(IRHOU,IRHOU) = eta/( dt + eta )
         invdS_dQ(IRHOV,IRHOV) = eta/( dt + eta )
         invdS_dQ(IRHOW,IRHOW) = eta/( dt + eta )
         invdS_dQ(IRHOE,IRHO)  = 0.5_RP*dt/eta * (POW2(u) + POW2(v) + POW2(w))
         invdS_dQ(IRHOE,IRHOU) = -dt*u/( dt + eta )
         invdS_dQ(IRHOE,IRHOV) = -dt*v/( dt + eta )
         invdS_dQ(IRHOE,IRHOW) = -dt*w/( dt + eta )
         invdS_dQ(IRHOE,IRHOE) =  1.0_RP 
#if defined(SPALARTALMARAS)
         invdS_dQ(IRHOTHETA,IRHOTHETA) = eta/( dt + eta )
#endif                                                                           
#endif
       end associate
 
   end subroutine IBM_SemiImplicitShiftJacobian
!
!  (1/dt - dS_dQt)
!
   subroutine IBM_SemiImplicitTurbulenceShiftJacobian( this, dt, invdS_dQ )
      use PhysicsStorage
      implicit none 
      !-arguments-----------------------------------
      class(IBM_type), intent(inout) :: this 
      real(kind=RP),   intent(in)    :: dt 
      real(kind=RP),   intent(inout) :: invdS_dQ(:,:)
      !-local-variables-----------------------------
      integer :: i

      invdS_dQ = dt*invdS_dQ

      do i = 1, NCONS 
         invdS_dQ(i,i) = 1.0_RP + invdS_dQ(i,i) 
      end do

   end subroutine IBM_SemiImplicitTurbulenceShiftJacobian
!
!  Second order Strang splitting correction Q^*. (1/dt - dS/dQ)^(-1)*Q^* = Q + dt*(S - dS/dQ*Q^*)
!  ------------------------------------------------------------------------------------------------
   subroutine IBM_GetSemiImplicitStep( this, eID, dt, Q, Q_target )
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------------------------
      class(IBM_type),           intent(inout) :: this
      integer,                   intent(in)    :: eID
      real(kind=rp),             intent(in)    :: dt
      real(kind=rp),             intent(inout) :: Q(NCONS)
      real(kind=rp),   optional, intent(in)    :: Q_target(NCONS)
      !-local-variables-----------------------------------------------
      real(kind=rp) :: dS_dQ(NCONS,NCONS), invdS_dQ(NCONS,NCONS), &
                       IBMSource(NCONS)

      call this% semiImplicitJacobian( eID, Q, dS_dQ )

      call this% semiImplicitShiftJacobian( eID, Q, dt, invdS_dQ ) 
      if( present(Q_target) ) then 
         call this% SourceTerm(eID = eID, Q = Q, Q_target = Q_target, Source = IBMSource, Wallfunction  =.false.)
      else
         call this% SourceTerm(eID = eID, Q = Q, Source = IBMSource, Wallfunction  =.false.)          
      end if 

      Q = matmul(invdS_dQ, Q + dt*( IBMSource - matmul(dS_dQ,Q) ))

   end subroutine IBM_GetSemiImplicitStep
!
!  Second order Strang splitting correction Q^*. (1/dt - dS/dQ)^(-1)*Q^* = Q + dt*(S - dS/dQ*Q^*)
!  ------------------------------------------------------------------------------------------------
   subroutine IBM_GetSemiImplicitStepTurbulence( this, ImagePoint, dt, Q, normal, dWall, STLNum )
      use PhysicsStorage
      use DenseMatUtilities
      implicit none
      !-arguments-------------------------------------------------------
      class(IBM_type),  intent(inout) :: this
      type(point_type), intent(in)    :: ImagePoint
      real(kind=rp),    intent(in)    :: dt, normal(:), dWall
      integer,          intent(in)    :: STLNum 
      real(kind=rp),    intent(inout) :: Q(:)
      !-local-variables-------------------------------------------------
      real(kind=rp) :: dS_dQ(NCONS,NCONS), invdS_dQ(NCONS,NCONS), &
                       TurbulenceSource(NCONS)

      call this% semiImplicitTurbulenceJacobian( ImagePoint, Q, normal, dWall, STLNum, dS_dQ )

      invdS_dQ = -dS_dQ

      call this% SemiImplicitTurbulenceShiftJacobian( dt, invdS_dQ )

      invdS_dQ = inverse(invdS_dQ)

      call this% SourceTermTurbulence( ImagePoint, Q, normal, dWall, STLNum, TurbulenceSource )       
 
      Q = matmul(invdS_dQ, Q + dt*( TurbulenceSource - matmul(dS_dQ,Q) ))
  
    end subroutine IBM_GetSemiImplicitStepTurbulence
!    
!   Splitting correction is applied. In the splitting, dt = dt/2
!   ------------------------------------------------------------
   subroutine IBM_SemiImplicitCorrection( this, elements, t, dt )
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      real(kind=RP),   intent(in)    :: t, dt
      !-local-variables-----------------------------
      real(kind=RP) :: Q_target(NCONS)
      integer       :: eID, i, j, k, iP

      if( .not. this% semiImplicit ) return

      if( this% Wallfunction ) call this% GetBandRegionStates( elements )
!$omp parallel
!$omp do schedule(runtime) private(i,j,k,Q_target)
      do eID = 1, SIZE( elements )
         associate(e => elements(eID))
         do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3)
            if( e% isInsideBody(i,j,k) ) then 
               if( this% stl(e% STL(i,j,k))% move ) then 
#if defined(NAVIERSTOKES)
                  Q_target = this% MaskVelocity( e% storage% Q(:,i,j,k), NCONS, e% STL(i,j,k), e% geom% x(:,i,j,k), t )
#endif
                  call this% GetSemiImplicitStep( eID, 0.5_RP*dt, e% storage% Q(:,i,j,k), Q_target )
               else
                  call this% GetSemiImplicitStep( eID, 0.5_RP*dt, e% storage% Q(:,i,j,k) ) 
               end if 
            end if 
         end do; end do; end do
         end associate
      end do
!$omp end do 
      if( this% Wallfunction ) then
!$omp do schedule(runtime) private(i,j,k)
         do iP = 1, this% NumOfForcingPoints
            associate( e => elements(this% ImagePoints(iP)% element_index))
            i = this% ImagePoints(iP)% local_position(1)
            j = this% ImagePoints(iP)% local_position(2)
            k = this% ImagePoints(iP)% local_position(3)
            call this% GetSemiImplicitStepTurbulence( this% ImagePoints(iP), 0.5_RP*dt, e% storage% Q(:,i,j,k),      &
                                                      e% geom% normal(:,i,j,k), e% geom% dWall(i,j,k), e% STL(i,j,k) )
            end associate 
         end do
!$omp end do 
      end if 
!$omp end parallel    
    end subroutine IBM_SemiImplicitCorrection
    
   subroutine IBM_GetBandRegionStates( this, elements )
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------------------
      class(IBM_type),            intent(inout) :: this
      type(element),              intent(in)    :: elements(:)
      !-local-variables-----------------------------------------
      integer :: STLNum

      do STLNum = 1, this% NumOfSTL
         call this% BandPoint_state( elements, STLNum, .false. )
      end do
    
   end subroutine IBM_GetBandRegionStates
!
!   Turbulent source term when wall function is applied
!   ---------------------------------------------------
   subroutine IBM_SourceTermTurbulence( this, ImagePoint, Q, normal, dWall, STLNum, TurbulenceSource )
      use PhysicsStorage
      use NodalStorageClass, only: NodalStorage
      implicit none
      !-arguments-----------------------------------------------------------------
      class(IBM_type),  intent(inout) :: this
      type(point_type), intent(in)    :: ImagePoint  
      real(kind=rp),    intent(in)    :: Q(:), dWall, normal(:)
      integer,          intent(in)    :: STLNum
      real(kind=rp),    intent(inout) :: TurbulenceSource(NCONS)
      !-local-variables-----------------------------------------------------------
      real(kind=rp) :: Q_IP(NCONS), Q_FP(NCONS), iP(NDIM)
      integer       :: i

      do i = 1, NCONS 
         Q_IP(i) = GetInterpolatedValue( this% BandRegion(STLNum)% Q(i,ImagePoint% nearestPoints), &
                                         ImagePoint% invPhi,                                       &
                                         ImagePoint% b,                                            &
                                         this% InterpolationType                                   )
      end do

      Q_FP = Q
#if defined(NAVIERSTOKES) 
      call ForcingPointState( Q_IP, this% IP_Distance, dWall, normal, Q_FP )
#endif
      call this% SourceTerm( ImagePoint% element_index, Q, Q_FP, TurbulenceSource, .true. )

      !TurbulenceSource = -1.0_RP/(1.0_RP*this% penalization(ImagePoint% element_index)) * (Q - Q_FP)

   end subroutine IBM_SourceTermTurbulence
!
!   State to be imposed on the forcing points due to the wall model
!   ---------------------------------------------------------------
#if defined(NAVIERSTOKES)  
   subroutine ForcingPointState( Q_IP, y_IP, y_FP, normal, Q_FP )
      use PhysicsStorage
      use WallFunctionDefinitions
      use WallFunctionBC
      use VariableConversion
      use FluidData   
#if defined(SPALARTALMARAS)
      use SpallartAlmarasTurbulence
#endif
      implicit none
      !-arguments--------------------------------------------------------------
      real(kind=rp), intent(in)    :: Q_IP(:), normal(:)
      real(kind=rp), intent(in)    :: y_IP, y_FP
      real(kind=rp), intent(inout) :: Q_FP(:)
      !-local-variables--------------------------------------------------------
      real(kind=rp)            :: nu_IP, nu_FP, mu_IP, mu_FP, T_IP, T_FP, &
                                  u_tau, u_IP(NDIM), u_IPt, u_FPt,        &
                                  u_FPn, u_FP(NDIM), u_IP_t(NDIM),        &
                                  tangent(NDIM), yplus_FP, nu_tilde,      &
                                  kappa_IP, kappa_FP
#if defined(SPALARTALMARAS)
      real(kind=rp)            :: Dump, chi, fv1, nu_t
#endif     
      call get_laminar_mu_kappa(Q_IP,mu_IP,kappa_IP)      
      nu_IP = mu_IP/Q_IP(IRHO)
            
      u_IP   = Q_IP(IRHOU:IRHOW)/Q_IP(IRHO)
       
      u_IP_t = u_IP - (dot_product(u_IP,normal) * normal)
   
      if( almostEqual(norm2(u_IP_t),0.0_RP) ) return 
      
      tangent = u_IP_t/norm2(u_IP_t)
   
      u_IPt = dot_product(u_IP,tangent)
 
      u_tau = u_tau_f(u_IPt, y_IP, nu_IP, u_tau0=.1_RP)

      call get_laminar_mu_kappa(Q_FP,mu_FP,kappa_FP)
      nu_FP = mu_FP/Q_FP(IRHO)

      u_FPt = u_plus_f(y_plus_f(y_FP, u_tau, nu_FP)) * u_tau

      u_FPn = dot_product(u_IP,normal) * y_FP/y_IP 
         
      u_FP = u_FPt*tangent + u_FPn*normal 

      T_FP = Temperature(Q_IP) + (dimensionless% Pr)**(1._RP/3._RP)/(2.0_RP*thermodynamics% cp) * (POW2(u_IPt) - POW2(u_FPt))      
#if defined(SPALARTALMARAS)      
      !yplus_FP = y_plus_f( y_FP, u_tau, nu_FP) 
      !Dump     = POW2(1.0_RP - exp(-yplus_FP/19.0_RP))      
      
      !chi  = QuarticRealPositiveRoot( 1.0_RP, -kappa*u_tau*y_FP*Dump/nu_FP, 0.0_RP, 0.0_RP, -kappa*u_tau*y_FP*Dump/nu_FP*POW3(SAmodel% cv1) )

      !nu_t = nu_FP * chi
      nu_t = kappa * u_tau * y_FP
#endif
     ! Q_FP(IRHO)  = Pressure(Q_IP)*refvalues% p/(thermodynamics% R * refValues% T * T_FP)
     ! Q_FP(IRHO)  = Q_FP(IRHO)/refvalues% rho   

      Q_FP(IRHOU:IRHOW) = Q_FP(IRHO) * u_FP

     ! Q_FP(IRHOE) = pressure(Q_IP)/( thermodynamics% gamma - 1._RP) + 0.5_RP*Q_FP(IRHO)*(POW2(u_FP(1)) + POW2(u_FP(2)) + POW2(u_FP(3)))
      Q_FP(IRHOE) = 0.5_RP*Q_FP(IRHO)*(POW2(u_FP(1)) + POW2(u_FP(2)) + POW2(u_FP(3)))
#if defined(SPALARTALMARAS)
      Q_FP(IRHOTHETA) = Q_FP(IRHO) * nu_t
#endif              
   end subroutine ForcingPointState
#endif      
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------------------
!  Ray tracing 
!  -------------------------------------------------------------
   subroutine IBM_CheckPoint( this, Point, STLNum, NumOfIntersections )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      real(kind=rp),   intent(inout) :: Point(:)
      integer,         intent(in)    :: STLNum  
      integer,         intent(inout) :: NumOfIntersections     
      !-local-variables------------------------------------------------------------
      type(KDtree),  pointer     :: tree 
      real(kind=rp)              :: RayDirection(NDIM), vecAxis(NDIM)
      integer                    :: Axis, minAxis(1), &
                                    lastIndex
      logical                    :: Upward, OnSurface, OnTriBound, delete
      type(ObjsDataLinkedList_t) :: Integer_List  

      RayDirection = 0.0_RP

      OnSurface = .false.; delete = .false.

      vecAxis = (/OBB(STLNum)% MBR% Length,OBB(STLNum)% MBR% Width,abs(OBB(STLNum)% nMax) + abs(OBB(STLNum)% nMin)/)
      
      minAxis = minloc(vecAxis)
      axis = minAxis(1)

      call RayAxis( this% root(STLNum)% maxAxis, this% ClipAxis, axis )

      if( Point(axis) .le. 0.0_RP ) then
         RayDirection(axis) = -1.0_RP
         Upward = .true.
            vecAxis(1) = OBB(STLNum)% LocVertices(1,1)
            vecAxis(2) = OBB(STLNum)% LocVertices(2,1)
            vecAxis(3) = OBB(STLNum)% LocVertices(3,1)
      else
         RayDirection(axis) = 1.0_RP
         Upward = .false.
            vecAxis(1) = OBB(STLNum)% LocVertices(1,7)
            vecAxis(2) = OBB(STLNum)% LocVertices(2,7)
            vecAxis(3) = OBB(STLNum)% LocVertices(3,7)
      end if

      Integer_List = ObjsDataLinkedList_t()

      NumOfIntersections = 0; lastIndex = -1

      do 

         call this% root(STLNum)% FindLeaf( Point, tree, .false. )
        
         if( tree% index .eq. lastIndex ) then
            call this% root(STLNum)% FindLeaf( Point, tree, .true. )
         endif

         lastIndex = tree% index

         call isPointInside( Point, RayDirection, this% root(STLNum)% ObjectsList, tree, Integer_List, NumOfIntersections, OnTriBound )

         if( OnTriBound ) delete = .true.
 
         if( Upward ) then
            Point(axis) = tree% vertices(axis,1) 
            if( Point(axis) .le. vecAxis(axis) ) exit
         elseif( .not. Upward ) then
            Point(axis) = tree% vertices(axis,7) 
            if( Point(axis) .ge. vecAxis(axis) ) exit
         end if

      end do

      if( delete ) NumOfIntersections = NumOfIntersections - 1

      call integer_List% Destruct()

   end subroutine IBM_CheckPoint
!
!  Intersection between a ray an a set of triangles 
!  ------------------------------------------------
   subroutine isPointInside( Point, RayDirection, ObjectsList, tree, Integer_List, NumOfIntersections, OnTriBound )
      use RealDataLinkedList
      use omp_lib
      implicit none
      !-arguments----------------------------------------------------------------
      real(kind=rp),              intent(in)    :: Point(:), RayDirection(:)
      type(object_type),          intent(in)    :: ObjectsList(:)
      type(KDtree),               intent(inout) :: tree
      type(ObjsDataLinkedList_t), intent(inout) :: Integer_List
      integer,                    intent(inout) :: NumOfIntersections
      logical,                    intent(inout) :: OnTriBound
      !-local-variables----------------------------------------------------------
      logical                        :: Intersect, found
      integer                        :: i, index

      OnTriBound = .false.

      if( tree% NumOfObjs .eq. 0 ) then
         return
      end if 
      
      do i = 1, tree% NumOfObjs
         index = tree% ObjsIndeces(i)
         found = integer_List% Check( index )
         if( .not. found ) then
            call Integer_List% Add( index )
            
            call PointIntersectTriangle( Point,ObjectsList(index)% vertices(1)% coords, &
                                               ObjectsList(index)% vertices(2)% coords, &
                                               ObjectsList(index)% vertices(3)% coords, &
                                               RayDirection, Intersect, OnTriBound      )  
            if( Intersect ) then
               NumOfIntersections = NumOfIntersections + 1  
            end if
         end if
      end do

   end subroutine isPointInside
!  
! This subroutine checks if a ray (RayDirection) starting from a point (Point) intersects
! a triangle in 3D space. If present, the intersection point Q is Q = P + t*RayDirection.
! If there is more than one intersections, the second one is not counted if it has the same t
! as the one previously found.
! See Fast, Minimum Storage Ray/Trinagle INtersection,  Moller TRumbore
!  --------------------------------------------------------------------------------------------   
   subroutine PointIntersectTriangle( Point, TriangleVertex1, TriangleVertex2, &
                                      TriangleVertex3, RayDirection,           &
                                      Intersect, OnTriBound                    )
      use MappedGeometryClass
      implicit none
      !-arguments---------------------------------------------------------------------------------------
      real(kind=rp), intent(in)    :: Point(NDIM), RayDirection(NDIM)
      real(kind=rp), intent(in)    :: TriangleVertex1(NDIM), TriangleVertex2(NDIM), TriangleVertex3(NDIM)
      logical,       intent(inout) :: Intersect, OnTriBound 
      !-local-variables----------------------------------------------------------------------------------
      real(kind=rp) :: E1vec(NDIM), E2vec(NDIM), Pvec(NDIM), &
                       Qvec(NDIM), Tvec(NDIM), N(NDIM),      &
                       Det, u, v, invDet, t
      logical       :: isInside
      
      Intersect  = .false.
      OnTriBound = .false.
      t          = -1.0_RP
      
      E1vec = TriangleVertex2 - TriangleVertex1 
      E2vec = TriangleVertex3 - TriangleVertex1
      Tvec   = Point - TriangleVertex1
      
      call vcross(RayDirection,E2vec,Pvec)
      Det = dot_product( E1vec, Pvec )

      If( almostEqual(Det,0.0_RP) ) return
      
      call vcross(Tvec,E1vec,Qvec)
      
      invDet = 1.0_RP/Det
      
      u = dot_product( Tvec, Pvec )*invDet

      if( u < 0.0_RP .or. u > 1.0_RP ) return
      
      v = dot_product( RayDirection, Qvec )*invDet

      if( v < 0.0_RP .or. u+v > 1.0_RP ) return

      t  = dot_product( E2vec, Qvec )*invDet

      if( almostequal(t,0.0_RP) ) return 
   
      ! Check if the point lies on the boundaries of the triangle
      !----------------------------------------------------------
      if( almostEqual(u,0.0_RP) .and. ( v .ge. 0.0_RP .and. v .le. 1.0_RP ) ) OnTriBound = .true. 
      if( almostEqual(v,0.0_RP) .and. ( u .ge. 0.0_RP .and. u .le. 1.0_RP ) ) OnTriBound = .true. 
      if( almostEqual(u+v,1.0_RP) )  OnTriBound = .true.

      if( t .gt. 0.0_RP ) Intersect = .true.


   end subroutine PointIntersectTriangle
!  
! Point-triangle intersection (point and triangle are on the same plane)
!  ------------------------------------------------   
   logical function isPointInsideTri( Point, TriangleVertex1, TriangleVertex2, &
                                      TriangleVertex3 ) result( isInside )
      use MappedGeometryClass
      implicit none
      !-arguments-----------------------------------------------------------------
      real(kind=rp), intent(in) :: Point(:), TriangleVertex1(:),          &
                                   TriangleVertex2(:), TriangleVertex3(:)
      !-local-variables-----------------------------------------------------------
      real(kind=rp) :: bb(NDIM), E0(NDIM), E1(NDIM), dd(NDIM), &
                       a, b, c, d, e, f, det, s, t
      integer       :: region
      
      isInside = .false.
      
      bb = TriangleVertex1
      E0 = TriangleVertex1 - TriangleVertex2
      E1 = TriangleVertex1 - TriangleVertex3
      dd = bb - Point
   
      a = dot_product(E0,E0)  
      b = dot_product(E0,E1)
      c = dot_product(E1,E1)
      d = dot_product(E0,dd)
      e = dot_product(E1,dd)
      f = dot_product(dd,dd)
      
      det = abs( a*c - b*b)
      s    = (b*e - c*d)/det
      t    = (b*d - a*e)/det
      
      region = FindRegion( det, s, t )
      
      if( region .eq. 0 ) isInside = .true.
      
   end function isPointInsideTri

   subroutine RayAxis( maxAxis, ClipAxis, axis )
      implicit none 
      !-arguments----------------------------------
      integer, intent(in)    :: maxAxis, ClipAxis
      integer, intent(inout) :: axis
      !-local-variables----------------------------
      integer :: vecAxis(NDIM), i

      vecAxis = (/1,2,3/)

      vecAxis(maxAxis)  = 0
      if( ClipAxis .ne. 0 ) vecAxis(ClipAxis) = 0

      if( axis .eq. maxAxis .or. axis .eq. ClipAxis ) then 
         do i = 1, NDIM
            if(vecAxis(i) .eq. 0) cycle
            axis = i; exit
         end do
      end if

   end subroutine RayAxis
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  Procedures for computing the point-STL distance
!  -------------------------------------------------
! 
! This function computes the minimum distance from a point to a triangle in 3D. 
! for more ditails see https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
!  ------------------------------------------------   
   subroutine MinimumPointTriDistance( Point, TriangleVertex1, TriangleVertex2, &
                                       TriangleVertex3, dist, IntersectionPoint ) 
      use MappedGeometryClass
      implicit none
      !-arguments--------------------------------------------------------------------
      real(kind=rp), intent(in)  :: Point(:), TriangleVertex1(:), &
                                    TriangleVertex2(:),        &
                                    TriangleVertex3(:)
      real(kind=rp), intent(out) :: IntersectionPoint(NDIM) 
      real(kind=rp), intent(out) :: dist     
      !-local-variables--------------------------------------------------------------
      real(kind=rp) :: bb(NDIM), E0(NDIM), E1(NDIM), dd(NDIM), & 
                       a, b, c, d, e, f, det, s, t,    &
                       tmp1, tmp0, numer, denom
      integer       :: region
      
      bb = TriangleVertex1
      E0 = TriangleVertex2 - bb
      E1 = TriangleVertex3 - bb
      dd = bb - Point
   
      a = dot_product(E0,E0) 
      b = dot_product(E0,E1) 
      c = dot_product(E1,E1) 
      d = dot_product(E0,dd) 
      e = dot_product(E1,dd) 
      f = dot_product(dd,dd)
 
      det = a*c - b*b  
      s   = b*e - c*d      
      t   = b*d - a*e

      if( (s + t) <= det ) then 
         if( s < 0.0_RP ) then 
            if( t < 0.0_RP ) then 
               region = 4 
            else 
               region = 3 
            end if 
         else if ( t < 0.0_RP ) then 
            region = 5 
         else 
            region = 0 
         end if 
      else 
         if( s < 0.0_RP ) then 
            region = 2
         elseif( t < 0.0_RP ) then 
            region = 6
         else 
            region = 1
         end if 
      end if 

      select case( region )
      case( 0 )
         s = s/det 
         t = t/det 
      case( 1 )
         numer = (c + e) - (b + d)
         if( numer <= 0.0_RP ) then 
            s = 0.0_RP 
         else 
            denom = a - 2.0_RP * b + c
            if( numer >= denom ) then 
               s = 1.0_RP 
            else 
               s = numer/denom
            end if 
         end if 
         t = 1.0_RP - s 
      case( 2 ) 
         tmp0 = b + d 
         tmp1 = c + e 
         if( tmp1 > tmp0 ) then 
            numer = tmp1 - tmp0 
            denom = a - 2.0_RP*b + c 
            if( numer >= denom ) then 
               s = 1.0_RP 
            else 
               s = numer/denom 
            end if 
            t = 1.0_RP - s 
         else
            s = 0.0_RP 
            if ( tmp1 <= 0.0_RP ) then 
               t = 1.0_RP 
            elseif( e >= 0.0_RP ) then 
               t = 0.0_RP 
            else 
               t = -e/c 
            end if 
         end if 
      case( 3 )
         s = 0.0_RP 
         if( e >= 0.0_RP ) then 
            t = 0.0_RP
         elseif( -e >= c ) then 
            t = 1.0_RP 
         else 
            t = -e/c 
         end if 
      case( 4 )
         if( d < 0.0_RP ) then 
            t = 0.0_RP 
            if( -d >= a ) then 
               s = 1.0_RP 
            else 
               s = -d/a 
            end if 
         else 
            s = 0.0_RP 
            if( e >= 0.0_RP ) then 
               t = 0.0_RP 
            elseif( -e >= c ) then 
               t = 1.0_RP 
            else 
               t = -e/c 
            end if 
         end if 
      case( 5 )
         t = 0.0_RP 
         if( d >= 0.0_RP ) then 
            s = 0.0_RP 
         elseif( -d >= a ) then 
            s = 1.0_RP 
         else 
            s = -d/a 
         end if 
      case( 6 )
         tmp0 = b + e 
         tmp1 = a + d 
         if( tmp1 > tmp0 ) then 
            numer = tmp1 - tmp0 
            denom = a - 2.0_RP*b + c 
            if( numer >= denom ) then 
               t = 1.0_RP 
            else 
               t = numer/denom 
            end if 
            s = 1.0_RP - t
         else  
            t = 0.0_RP 
            if( tmp1 <= 0.0_RP ) then 
               s = 1.0_RP 
            elseif( d >= 0.0_RP ) then 
               s = 0.0_RP 
            else 
               s = -d/a 
            end if 
         end if 
      end select

      IntersectionPoint = TriangleVertex1 + s*E0 + t*E1    
        
      dist = norm2(Point - IntersectionPoint)

   end subroutine MinimumPointTriDistance
!  
!  Identification of the section, see below.
!  --------------------------------------------------------------------------
!      \ region2 ^t
!       \        |
!        \       |
!         \      | 
!          \     |
!           \    |
!            \   |
!             \  |
!              \ |
!               \|
!                *
!                |\    P2
!                | \
!                |  \
!                |   \
!                |    \
!                |     \
!                |      \
!                |       \
!                |        \
!       region3  |         \ region1
!                |          \
!                | region0   \
!                |            \
!                |             \ P1
! ---------------*--------------*------->s
!                |P0             \
!  region4       |    region5     \ region6
   integer function FindRegion( det, s, t ) result( region )
      implicit none
      !-arguments-----------------------------
      real(kind=rp), intent(in) :: det, s, t
      
      if( s+t <= det ) then
         if( s < 0 ) then
            if( t < 0 ) then
               region = 4
            else
               region = 3
            end if
         elseif( t < 0 ) then
            region = 5
         else
            region = 0
         end if
      else
         if( s < 0 ) then
            region = 2
         elseif ( t < 0 ) then
            region = 6
         else
            region = 1
         end if
      end if 
      
   end function FindRegion
!
! This subroutine computes the minumum distance from x (global ref. frame) to the body. First,
! the box (tree) where x lies is found, then the min distance between the objects inside the tree and
! the point is computed. If the sphere whose radius is minDist, is enclosed in the box, the subroutine stops.
! If the latter condition is not true, all the boxes that intersects the sphere are checked in the same way as
! the initial one, getting new_minDist. If a lower distance is found, minDist is updated.
! minDist is the SQUARE of the actual distance.
!  -------------------------------------------------------------------------------------------------------------   
   subroutine MinimumDistance( Point, root, minDist, normal )
      implicit none
      !-arguments----------------------------------------------------------
      real(kind=rp), intent(in)    :: Point(:)
      type(KDtree),  intent(inout) :: root
      real(kind=rp), intent(inout) :: minDist
      real(kind=rp), intent(out)   :: normal(NDIM)    
      !-local-variables-----------------------------------------------------
      real(kind=rp)         :: IntersPoint(NDIM), new_IntersPoint(NDIM), &
                               dsvec(NDIM), x(NDIM), P(NDIM), IP(NDIM),  &
                               IntersectionPoint(NDIM), Dist,            &
                               New_minDist, Radius, ds
      logical               :: Intersect
      type(KDtree), pointer :: tree, newGuess
      integer               :: i, index, LeafIndex,             &
                               TriangleIndex, New_TriangleIndex
                                        
      minDist = huge(1.0_RP)

      call root% FindLeaf( Point, tree, .false. )

      LeafIndex = tree% index

      do i = 1, tree% NumOfObjs
         index = tree% ObjsIndeces(i)

         call MinimumPointTriDistance( Point, root% ObjectsList(index)% vertices(1)% coords, &
                                       root% ObjectsList(index)% vertices(2)% coords,        &
                                       root% ObjectsList(index)% vertices(3)% coords, Dist,  &
                                       IntersPoint                                           )
          
         if( Dist .lt. minDist ) then
            minDist           = Dist
            IntersectionPoint = IntersPoint   
            TriangleIndex     = index           
         end if
      end do  

      if( tree% NumOfObjs .gt. 0 ) then
      ! Check the sphere
      !-----------------
         Radius    = POW2(minDist)
         Intersect = CheckHypersphere( tree, Point, Radius )
      else
         print *, "IBM:: MinimumDistance: "
         print *, "Can't find triangles in leaf ", LeafIndex
         error stop
      end if

      nullify(tree) 

      if( Intersect ) then
         New_minDist = huge(1.0_RP)
         call MinimumDistOtherBoxes( Point, root, root, Radius, LeafIndex, New_minDist, New_IntersPoint, New_TriangleIndex )
         if( New_minDist .lt. minDist ) then 
            minDist           = New_minDist
            IntersectionPoint = New_IntersPoint
            TriangleIndex     = New_TriangleIndex
         end if
      end if    

      call OBB(root% STLNum)% ChangeRefFrame( Point, GLOBAL, P ) 
      call OBB(root% STLNum)% ChangeRefFrame( IntersectionPoint, GLOBAL, IP ) 
      
      normal = (P - IP)/norm2(P - IP)
              
   end subroutine MinimumDistance
!  
!  Distance between the point and triangles in other boxes, i.e. not the one containing it 
!  -----------------------------------------------------------------------------------------
   recursive subroutine MinimumDistOtherBoxes( Point, root, tree, Radius, LeafIndex,           &
                                               New_minDist, New_IntersPoint, New_TriangleIndex )
                                                                              
      implicit none
      !-arguments-----------------------------------------------
      real(kind=rp), intent(in)    :: Point(:)
      type(KDtree),  intent(in)    :: root
      type(KDtree),  intent(inout) :: tree
      real(kind=rp), intent(in)    :: Radius
      integer,       intent(in)    :: LeafIndex
      real(kind=rp), intent(inout) :: New_minDist
      real(kind=rp), intent(inout) :: New_IntersPoint(:)
      integer,       intent(inout) :: New_TriangleIndex 
      !-local-variables-----------------------------------------
      real(kind=rp)                  :: IntersPoint(NDIM), Dist
      logical                        :: Intersect
      integer                        :: i, index
      
      Intersect = BoxIntersectSphere( Radius, Point, tree% vertices )
      
      if( Intersect ) then
         if( tree% isLast ) then
            if( LeafIndex .ne. tree% index ) then 
               do i = 1, tree% NumOfObjs
                  index = tree% ObjsIndeces(i)
                  call MinimumPointTriDistance( Point, root% ObjectsList(index)% vertices(1)% coords, &
                                                root% ObjectsList(index)% vertices(2)% coords,        &
                                                root% ObjectsList(index)% vertices(3)% coords, Dist,  &
                                                IntersPoint                                           )
                  if( Dist .lt. New_minDist ) then
                     New_minDist       = Dist
                     New_IntersPoint   = IntersPoint
                     New_TriangleIndex = index      
                  end if
               end do    
            end if
         else
            call MinimumDistOtherBoxes( Point, root, tree% child_L,        & 
                                        Radius, LeafIndex, New_minDist,    &
                                        New_IntersPoint, New_TriangleIndex )    
            call MinimumDistOtherBoxes( Point, root, tree% child_R,        &
                                        Radius, LeafIndex, New_minDist,    &
                                        New_IntersPoint, New_TriangleIndex )
         end if
      end if
          
   end subroutine MinimumDistOtherBoxes
   
   subroutine MinimumDistanceSphereKDtree( Point, root, Radius, minDist, normal )
      implicit none
      !-arguments----------------------------------------
      real(kind=rp), intent(in)    :: Point(:), Radius
      type(KDtree),  intent(inout) :: root
      real(kind=rp), intent(inout) :: minDist
      real(kind=rp), intent(out)   :: normal(NDIM)    
      !-local-variables----------------------------------
      real(kind=RP) :: P(NDIM), IP(NDIM),       &
                       IntersectionPoint(NDIM)
      integer       :: TriangleIndex
      
      minDist = huge(1.0_RP)

      if( .not. BoxIntersectSphere( Radius, Point, root% vertices ) ) return

      call MinimumDistOtherBoxes( Point, root, root, Radius, 0, minDist, &
                                  IntersectionPoint, TriangleIndex       )
   
      call OBB(root% STLNum)% ChangeRefFrame( Point, GLOBAL, P ) 
      call OBB(root% STLNum)% ChangeRefFrame( IntersectionPoint, GLOBAL, IP ) 
      
      normal = (P - IP)/norm2(P - IP)      
   
   end subroutine MinimumDistanceSphereKDtree
   
   
!
!  Once the first distance is found, a sphere whose radius is equal to the distance is computed. If 
!  it lies inside the starting box, nothing is done; otherwise the boxes intersecting the sphere are
!  selected
!  ---------------------------------------------------------------------------------------------------
   logical function CheckHypersphere( tree, Point, minDist) result( Intersect )
   
      implicit none
      !-arguments-----------------------------------------------------
      real(kind=rp),         intent(in)    :: Point(:)
      type(KDtree),  target, intent(inout) :: tree
      real(kind=rp),         intent(inout) :: minDist
      
      Intersect = BoxIntersectSphere( minDist, Point, tree% vertices )
   
   end function CheckHypersphere
!  
!  Circle-triangle intersection
!  ---------------------------------------------------------
   logical function CircleRectIntersection( RectLength, RectWidth, RectCenter, Center, Radius ) result( Intersect )

      implicit none
      !-arguments------------------------------------------------
      real(kind=rp), intent(in) :: RectCenter(:), Center(:)
      real(kind=rp), intent(in) :: RectLength, RectWidth, Radius 
      !-local-variables------------------------------------------
      real(kind=rp) :: d(2)

      Intersect = .false.
      
      d = Center - RectCenter

      If( abs(d(1)) .gt. 0.5_RP*RectLength + Radius ) return
      If( abs(d(2)) .gt. 0.5_RP*RectWidth + Radius ) return

      if( abs(d(1)) .le. 0.5_RP*RectLength .and. &
          abs(d(2)) .le. 0.5_RP*RectWidth ) then
         Intersect = .true.
         return
      end if
      
      if( POW2(abs(d(1)) - 0.5_RP*RectLength) + POW2(abs(d(2)) - 0.5_RP*RectWidth) .le. &
          POW2(Radius) ) Intersect = .true.

   end function CircleRectIntersection
!  
!  Point-box distance
!  ------------------------------------------------
   real(kind=rp) function PointBoxDistance( Point, vertices ) result( sqrDist )
   
      implicit none
      !-arguments--------------------------------
      real(kind=rp), intent(in) :: Point(:)
      real(kind=rp), intent(in) :: vertices(:,:)
      !-local-variables--------------------------
      integer :: i
   
      sqrDist = 0.0_RP
   
      ! if the point's x-coordinate (xp) is less than the x-min coord. of the box, then
      ! the minimum distance on the x-dir is the distance between the xp and x-min, or
      ! vertices(i,1) - x-coordinate. The following loop avoids the computation of the point
      ! with the minimum distance from P. If the point P is inside the box, sqrDist = 0.  
      !--------------------------------------------------------------------------------------      
      do i = 1, NDIM
         if( Point(i) .lt. minval(vertices(i,:)) ) sqrDist = sqrDist + POW2(minval(vertices(i,:))-Point(i))
         if( Point(i) .gt. maxval(vertices(i,:)) ) sqrDist = sqrDist + POW2(Point(i)-maxval(vertices(i,:)))
      end do
   
   end function PointBoxDistance
!  
!  Sphere-box intersection, the radius must be SQUARED
!  ------------------------------------------------
   logical function BoxIntersectSphere( Radius, Center, vertices ) result( Intersect )
   
      implicit none
      !-arguments-------------------------------------------      
      real(kind=rp), intent(in) :: vertices(:,:)
      real(kind=rp), intent(in) :: Center(:)
      real(kind=rp), intent(in) :: Radius
      !-local-variables-------------------------------------
      real(kind=rp) :: sqrDist
      
      Intersect = .false.
      sqrDist = PointBoxDistance( Center, vertices )
      
      if( sqrDist .le. Radius ) Intersect = .true.
      
   end function BoxIntersectSphere
!
!  Nearest-neighbor algorithm for selecting the DoFs closer to Point. 
!  --------------------------------------------------------------------
   subroutine MinimumDistancePoints( Point, root, BandRegion, minDist, actualIndex, PointsIndex )
      use ElementClass
      implicit none
      !-arguments-------------------------------------------------------
      real(kind=rp),   intent(in)    :: Point(:)
      type(KDtree),    intent(inout) :: root
      type(IBMpoints), intent(in)    :: BandRegion
      real(kind=rp),   intent(inout) :: minDist
      integer,         intent(in)    :: actualIndex
      integer,         intent(inout) :: PointsIndex(:)
      !-local-variables-------------------------------------------------
      real(kind=rp)         :: BandPoint(NDIM), sqrDist, &
                               New_sqrDist, Radius
      logical               :: Intersect
      type(KDtree), pointer :: tree
      integer               :: i, LeafIndex, new_index, k
           
      minDist = huge(1.0_RP)

      call root% FindLeaf( Point, tree, .false. ) 
      
      LeafIndex = tree% index

      do i = 1, tree% NumOfObjs
               
         if( any(PointsIndex .eq. tree% ObjsIndeces(i)) ) cycle

         BandPoint = BandRegion% x(tree% ObjsIndeces(i))% coords
         
         sqrDist = 0.0_RP
         do k = 1, NDIM
            sqrDist = sqrDist + POW2(Point(k) - BandPoint(k))
         end do

         if( sqrDist .lt. minDist .or. almostEqual(sqrDist,minDist) ) then
               minDist                  = sqrDist
               PointsIndex(actualIndex) = tree% ObjsIndeces(i)
         end if 
         
      end do

      if( tree% NumOfObjs .gt. 0 ) then
      ! Check the sphere
      !-----------------
         if( PointsIndex(actualIndex) .eq. 0 ) then 
            BandPoint = BandRegion% x(PointsIndex(actualIndex-1))% coords
         else
            BandPoint = BandRegion% x(PointsIndex(actualIndex))% coords
         end if

         sqrDist = 0.0_RP
         do k = 1, NDIM
            sqrDist = sqrDist + POW2(Point(k) - BandPoint(k))
         end do
         
         Radius = sqrDist  !minDist
         Intersect = CheckHypersphere( tree, Point, Radius )
      else 
         print *, "IBM:: MinimumDistance: "
         print *, "Can't find triangles in leaf ", LeafIndex
         error stop
      end if
      
      nullify(tree)

      if( Intersect ) then
         New_sqrDist = huge(1.0_RP)
         call MinimumDistOtherBoxesPoints( Point, root, BandRegion, Radius, &
                                           New_sqrDist, PointsIndex,        &
                                           LeafIndex, new_index             )         
         if( New_sqrDist .le. minDist ) then
            minDist = New_sqrDist; PointsIndex(actualIndex) = new_index 
         end if
      end if    

      minDist = sqrt(minDist)
      
   end subroutine MinimumDistancePoints
!  
!  Distance between the point and DoFs in other boxes, i.e. not the one containing it 
!  -----------------------------------------------------------------------------------------
   recursive subroutine MinimumDistOtherBoxesPoints( Point, tree, BandRegion, Radius, &
                                                     New_sqrDist,  PointsIndex,       &
                                                     LeafIndex,  new_index            )
      use elementClass                                                                 
      implicit none
      !-arguments---------------------------------------------------
      real(kind=rp),   intent(in)    :: Point(:) 
      type(KDtree),    intent(inout) :: tree
      type(IBMpoints), intent(in)    :: BandRegion
      real(kind=rp),   intent(in)    :: Radius
      real(kind=rp),   intent(inout) :: New_sqrDist
      integer,         intent(inout) :: new_index
      integer,         intent(in)    :: PointsIndex(:)
      integer,         intent(in)    :: LeafIndex
      !-local-variables---------------------------------------------
      real(kind=rp) :: sqrDist, BandPoint(NDIM)
      logical       :: Intersect
      integer       :: i, k
      
      Intersect = BoxIntersectSphere( Radius, Point, tree% vertices )
      
      if( Intersect ) then
         if( tree% isLast ) then
            if( LeafIndex .ne. tree% index ) then
               do i = 1, tree% NumOfObjs
           
                  if( any(PointsIndex .eq. tree% ObjsIndeces(i)) ) cycle

                  BandPoint = BandRegion% x(tree% ObjsIndeces(i))% coords
                 
                  sqrDist = 0.0_RP
                  do k = 1, NDIM
                     sqrDist = sqrDist + POW2(Point(k) - BandPoint(k))
                  end do
                  
                  if( sqrDist .lt. New_sqrDist .or. almostEqual(sqrDist,New_sqrDist)  )then
                     New_sqrDist = sqrDist
                     new_index   = tree% ObjsIndeces(i)          
                  end if

               end do    
            end if
         else
            call MinimumDistOtherBoxesPoints( Point, tree% child_L, BandRegion, &
                                              Radius, New_sqrDist,PointsIndex,  &
                                              LeafIndex, new_index              )    
            call MinimumDistOtherBoxesPoints( Point, tree% child_R, BandRegion, &
                                              Radius, New_sqrDist, PointsIndex, &
                                              LeafIndex, new_index              )
         end if
      end if
          
   end subroutine MinimumDistOtherBoxesPoints


   subroutine GetMatrixInterpolationSystem( Point, x, invPhi, b, INTERPOLATION )
      use DenseMatUtilities
      implicit none
      !-arguments------------------------------------
      real(kind=RP),     intent(in)    :: Point(:)
      type(point_type),  intent(in)    :: x(:)
      real(kind=RP),     intent(inout) :: invPhi(:,:), b(:)
      integer,           intent(in)    :: INTERPOLATION 
      !-local-variables------------------------------
      real(kind=RP) :: Phi(size(x),size(x)),   &
                       dist(size(x),size(x)),  &
                       L(size(x)), d,          &
                       P(size(x),size(x)),     &
                       P_T(size(x),size(x)),   &
                       V(2*size(x),2*size(x)), &
                       lambda(size(x)),        &
                       pp(size(x)), w(size(b))         
      integer       :: i, j, k

      select case( INTERPOLATION )
         case( EXPONENTIAL )

            do i = 1, size(x)
               do j = i, size(x)
                  dist(i,j) = norm2(x(i)% coords - x(j)% coords)
                  dist(j,i) = dist(i,j)
               end do 
            end do

            expCoeff = 0.001_RP!norm2(Point - x(1)% coords) 

            do i = 1, size(x)
               do j = i, size(x)
                  Phi(i,j) = interpolationfunction(dist(i,j), EXPONENTIAL )
                  Phi(j,i) = Phi(i,j)
               end do
               d    = norm2(Point - x(i)% coords)
               b(i) = interpolationfunction(d, EXPONENTIAL)
            end do

            invPhi = inverse(Phi) 

         case( IDW )

            b = 0.0_RP; invPhi = 0.0_RP
            do i = 1, size(x)
               d           = norm2(Point - x(i)% coords)
               invPhi(i,i) = 1.0_RP/d
               b           = b + invPhi(i,i)
            end do  
            do i = 1, size(x)
               b(i) = 1.0_RP/b(i)
            end do 

         case( POLYHARMONIC_SPLINE )

            b = 0.0_RP; invPhi = 0.0_RP; Phi = 0.0_RP
            
            do i = 1, size(x)
               do j = i, size(x)
                  dist(i,j) = norm2(x(i)% coords - x(j)% coords)
                  dist(j,i) = dist(i,j)
               end do 
            end do
           
            do i = 1, size(x) 
               do j = i, size(x) 
                  Phi(i,j) = interpolationfunction(dist(i,j), POLYHARMONIC_SPLINE )
                  Phi(j,i) = Phi(i,j)
               end do 
               d    = norm2(Point - x(i)% coords)
               b(i) = interpolationfunction(d, POLYHARMONIC_SPLINE)
               call PolynomialVector( x(i)% coords(1), x(i)% coords(2), x(i)% coords(3), P(i,:) ) 
            end do  

            P_T = transpose(P)

            call buildMatrixPolySpline( Phi, P, P_T, V )

            invPhi = inverse(V)

         case( POLYNOMIAL )

            Phi = 0.0_RP; invPhi = 0.0_RP; b = 0.0_RP
            do i = 1, size(x)
               call PolynomialVector( x(i)% coords(1), x(i)% coords(2), x(i)% coords(3), Phi(i,:) )
            end do 

            call PolynomialVector( Point(1), Point(2), Point(3), b )

            invPhi = inverse(Phi)

         case( MLS )

            Phi = 0.0_RP; invPhi = 0.0_RP; b = 0.0_RP
            expCoeff = norm2(Point - x(1)% coords) 
            do i = 1, size(x)
               call PolynomialVector( x(i)% coords(1)-Point(1), x(i)% coords(2)-Point(2), x(i)% coords(3)-Point(3), P(i,:) )
               d    = norm2(Point - x(i)% coords)
               W(i) = interpolationfunction( d, MLS )
            end do 

            call PolynomialVector( 0.0_RP, 0.0_RP, 0.0_RP, pp )

            do j = 1, size(p,2)
               do k = 1, size(p,2)
                  do i = 1, size(x)
                     Phi(j,k) = Phi(j,k) + P(i,j)*P(i,k)*W(i)
                  end do
               end do 
            end do

            lambda = matmul(inverse(Phi),pp)

            do i = 1, size(p,2)
               do j = 1, size(p,2)
                  b(i) = b(i) + lambda(j) * p(i,j)
               end do 
               b(i) = W(i)*b(i)
            end do 

      end select 

   end subroutine GetMatrixInterpolationSystem

   real(kind=RP) function interpolationfunction( x, INTERPOLATION, k ) result(f)
      implicit none 
      real(kind=RP),           intent(in) :: x 
      integer,                 intent(in) :: INTERPOLATION 
      real(kind=rp), optional, intent(in) :: k 

      select case( INTERPOLATION )
         case( EXPONENTIAL )
         
            f = exp(-POW2(x/expCoeff))
         
         case( POLYHARMONIC_SPLINE )
         
            if( abs(x) .lt. 1.0d-12 ) then 
               f = 0.0_RP 
               return 
            end if 
            if( present(k) ) then 
               f = x**k * log(x)
            else 
               f = POW2(x) * log(x)
            end if 
         
         case( MLS )

            f = 1.0_RP/(POW2(x) + 1.0d-12)

      end select 
   end function interpolationfunction

   real(kind=RP) function GetInterpolatedValue( forcing, invPhi, b, INTERPOLATION, x, y, z ) result( value )
   implicit none 
   !-arguments------------------------------------------------------------
   real(kind=RP),           intent(in) :: forcing(:), invPhi(:,:), b(:)
   integer,                 intent(in) :: INTERPOLATION
   real(kind=RP), optional, intent(in) :: x, y, z
   !-local-variables------------------------------------------------------
   real(kind=RP), allocatable :: weights(:), f(:)

   select case( INTERPOLATION )
      case( POLYHARMONIC_SPLINE )
         
         allocate( weights(2*size(b)), &
                   f(2*size(b))        )
         f            = 0.0_RP 
         f(1:size(b)) = forcing 
         weights = matmul(invPhi,f)
         value = dot_product(weights(1:size(b)),b) +                   &
                 EvaluateInterp( weights(size(b)+1:2*size(b)), x, y, z )
         deallocate(weights,f)
      case( MLS )

         value = dot_product(forcing,b)

      case default 
         allocate(weights(size(b)))
         weights = matmul(invPhi,forcing)
         value = dot_product(weights,b)
         deallocate(weights)
   end select 

   end function GetInterpolatedValue
!  
!  Inverse Distance Weighting interpolation
!  ----------------------------------------
   subroutine GetIDW_value( Point, x, Q, value )
      use PhysicsStorage
      use MappedGeometryClass
      implicit none
      !-arguments----------------------------------------------------
      real(kind=RP),    intent(in)  :: Point(:)
      type(point_type), intent(in)  :: x(:)
      real(kind=RP),    intent(in)  :: Q(:,:)
      real(kind=RP),    intent(out) :: value(NCONS)
      !-local-variables----------------------------------------------
      real(kind=RP)            :: DistanceNormal(NDIM), d,  &
                                  distanceSqr, sqrd1,       &
                                  num(NCONS), den, xP(NDIM)
      real(kind=RP), parameter :: eps = 1.0d-10
      integer                  :: i, k
   
      num = 0.0_RP; den = 0.0_RP
   
      do i = 1, size(x)
           
         d = norm2(x(i)% coords - Point)

         num = num + Q(:,i)/(d+eps)
         den = den + 1.0_RP/(d+eps)

      end do 
   
      value = num/den
   
   end subroutine GetIDW_value
   
! estimate the yplus for a flat plate 
   
   real(kind=RP) function InitializeDistance( y_plus ) result( y )
      use FluidData
      implicit none
      !-arguments--------------------------
      real(kind=rp), intent(in) :: y_plus
      !-local-varirables-------------------
      real(kind=RP) :: nu, Cf
      
      y = 0.0_RP
#if defined(NAVIERSTOKES)     
      nu = refValues% mu / refValues% rho
#endif
      Cf = Estimate_Cf()
#if defined(NAVIERSTOKES)      
      y = sqrt(2.0_RP) * y_plus * nu /(refValues% V * sqrt(Cf)) 
#endif  
   end function InitializeDistance
   
   real(kind=RP) function Estimate_Cf() result( Cf )
      use FluidData
      implicit none

      Cf = 0.0_RP      
#if defined(NAVIERSTOKES)   
     ! Schlichting, Hermann (1979), Boundary Layer Theory... ok if Re < 10^9
!~       Cf = (2.0_RP*log10(dimensionless% Re) - 0.65_RP)**(-2.3_RP)     
      Cf = 0.058_RP * dimensionless% Re**(-0.2_RP)     
#endif        
   end function Estimate_Cf
   
   real(kind=RP) function Estimate_u_tau( ) result( u_tau )
      use FluidData
      implicit none
      !-local-variables--------------------------------
      real(kind=RP) :: Cf
      
      u_tau = 0.0_RP
      
      Cf = Estimate_Cf()
#if defined(NAVIERSTOKES) 
      u_tau = sqrt( 0.5_RP * POW2(refValues% V) * Cf ) !sqrt(tau_w/rho)
#endif      
   end function Estimate_u_tau
   
   real(kind=RP) function GetEstimated_y( y_plus, nu, u_tau ) result( y )
      use PhysicsStorage 
      implicit none
      !-arguments--------------------------------------
      real(kind=RP), intent(in) :: y_plus, nu, u_tau
   
      y = y_plus * nu / u_tau
      
      y = y/Lref 
   
   end function GetEstimated_y


   real(kind=RP) function QuarticRealPositiveRoot(a0, b0, c0, d0, e0) result( value )
      !https://math.stackexchange.com/questions/115631/quartic-equation-solution-and-conditions-for-real-roots
      implicit none
      !-arguments--------------------------------------
      real(kind=rp), intent(in) :: a0, b0, c0, d0, e0
      !-local-variables--------------------------------
      real(kind=rp) :: a, b, c, d, x, m2, m, n,  &
                       alpha, beta, psi,         & 
                       solution(2)

      a = b0/a0; b = c0/a0; c = d0/a0; d = e0/a0
      
      x = cubic_poly( 1.0_RP, -b, a*c-4.0_RP*d, d*(4.0_RP*b-POW2(a)) - POW2(c) )

      m2 = 0.25_RP*POW2(a)-b+x
      
      if( m2 .gt. 0.0_RP ) then
         m = sqrt(m2)
         n = (a*x-2.0_RP*c)/(4.0_RP*m)
      elseif( almostEqual(m2,0.0_RP) ) then
         m = 0.0_RP
         n = sqrt(0.25_RP*POW2(x) - d)
      else
         print *, 'no real roots found'
         error stop
      end if
      
      alpha = 0.5_RP*POW2(a) - x - b
      
      beta = 4.0_RP * n - a * m
      
      if( (alpha + beta) .ge. 0.0_RP ) then
         psi = sqrt(alpha+beta)
         solution(1) = -0.25_RP*a + 0.5_RP*m + 0.5_RP*psi
         solution(2) = -0.25_RP*a + 0.5_RP*m - 0.5_RP*psi
      elseif( (alpha - beta) .ge. 0.0_RP ) then
         psi = sqrt(alpha-beta)
         solution(1) = -0.25_RP*a - 0.5_RP*m + 0.5_RP*psi
         solution(2) = -0.25_RP*a - 0.5_RP*m - 0.5_RP*psi        
      end if
      !take positive solution
      value = maxval(solution) 

   end function QuarticRealPositiveRoot
   
   real(KIND=rp) function cubic_poly( a0, b0, c0, d0 ) result( x )
   
      implicit none
      !-arguments----------------------------------------------
      real(kind=RP),      intent(in) :: a0, b0, c0, d0
      !-local-variables----------------------------------------
      real(kind=RP) :: a, b, c, Jac, x0, q, r
      real(kind=RP) :: theta, m, s, txx
      integer :: i
      
      
      a = b0/a0; b = c0/a0; c = d0/a0
            
      x0 = 1.0_RP
      
      do  i = 1, 100
      
         Jac = 3.0_RP*a0*POW2(x0) + 2.0_RP*b0*x0 + c0
   
         x = x0 - eval_cubic(a0, b0, c0, d0, x0)/Jac
   
         if( abs(x-x0) .lt. 1.0E-10_RP ) return
   
         x0 = x
   
      end do
    
      print *, "cubic_poly:: Newtown method doesn't coverge"
      error stop
    
   end function cubic_poly

   real(kind=RP) function eval_cubic( a0, b0, c0, d0, x0 ) result( value )
   
      implicit none
      !-arguments----------------------------------
      real(kind=RP), intent(in) :: a0, b0, c0, d0, x0
      
      value = a0*POW3(x0) + b0*POW2(x0) + c0*x0 + d0
   
   end function eval_cubic
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  ---------------------------------------------------
!  Subourinte for generating a .tec file with forces 
!  ---------------------------------------------------   
   subroutine STLScalarTEC( x, y, z, scalar, STLNum, FileName, title, variables )
   
      implicit none
      !-arguments--------------------------------------------------
      real(kind=RP),     intent(in) :: x(:), y(:), z(:), scalar(:)
      integer,           intent(in) :: STLNum
      character(len=*),  intent(in) :: FileName, title, variables
      !-local-variables--------------------------------------------
      integer :: NumOfObjs, funit
      
      NumOfObjs = size(x)
  
      call TecHeader( funit, FileName, title, variables ) 
      call WriteScalarQuatity( funit, NumOfObjs, x, y, z, scalar )
      call TECtriangle_3points( funit, NumOfObjs )

      close(funit)

   end subroutine STLScalarTEC
   
   subroutine STLvectorTEC( x, y, z, vector_x, vector_y, vector_z, STLNum, FileName, title, variables )
   
      implicit none
      !-arguments---------------------------------------------------------
      real(kind=RP),     intent(in) :: x(:), y(:), z(:), vector_x(:), &
                                       vector_y(:), vector_z(:)
      integer,           intent(in) :: STLNum
      character(len=*),  intent(in) :: FileName, title, variables
      !-local-variables--------------------------------------------------
      integer                    :: NumOfObjs, funit
      
      NumOfObjs = size(x)
  
      call TecHeader( funit, FileName, title, variables )
      call WriteVectorQuatity( funit, NumOfObjs, x, y, z, vector_x, vector_y, vector_z )
      call TECtriangle_3points( funit, NumOfObjs )

      close(funit)
       
   end subroutine STLvectorTEC
   
   subroutine TecHeader( funit, FileName, title, variables )
   
      implicit none
      !-arguments-------------------
      integer,          intent(inout) :: funit
      character(len=*), intent(in)    :: FileName, title, &
                                         variables
   
      funit = UnusedUnit()
   
      open(funit,file='RESULTS/'//trim(FileName), status='unknown')    
   
      write(funit,'(A)') 'Title="' //trim(title)// '"'
      write(funit,'(A)') 'Variables='//trim(variables)
   
   end subroutine TecHeader
   
   subroutine WriteScalarQuatity( funit, NumOfObjs, x, y, z, scalar )
   
      implicit none
      !-arguments-------------------------------------------------------
      integer,       intent(in) :: funit, NumOfObjs
      real(kind=RP), intent(in) :: x(:), y(:), z(:), scalar(:)
      !-local-variables-------------------------------------------------
      integer :: i, j
     
      write(funit,'(a,i6,a,i6,a)') 'ZONE NODES =', NumOfObjs, ',ELEMENTS =', NumOfObjs/3,', DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
      
      do i = 1, NumOfObjs
         write(funit,'(4g15.6)') x(i), y(i), z(i), scalar(i)
      end do
   
   end subroutine WriteScalarQuatity
   
   subroutine WriteVectorQuatity( funit, NumOfObjs, x, y, z, vector_x, vector_y, vector_z )
   
      implicit none
      !-arguments---------------------------------------------------------
      integer,       intent(in) :: funit, NumOfObjs
      real(kind=RP), intent(in) :: x(:), y(:), z(:), vector_x(:), &
                                   vector_y(:), vector_z(:) 
      !-local-variables---------------------------------------------------
      integer :: i, j
     
      write(funit,'(a,i6,a,i6,a)') 'ZONE NODES=', NumOfObjs, ',ELEMENTS =', NumOfObjs/3,', DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
      
      do i = 1, NumOfObjs
         write (funit,'(6g15.6)') x(i), y(i), z(i), vector_x(i), vector_y(i), vector_z(i)
      end do
   
   end subroutine WriteVectorQuatity
   
   subroutine WriteTimeFile( t, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments---------------------------------------------------------
      real(kind=RP), intent(in) :: t
      integer ,      intent(in) :: STLNum  
      !-local-variables---------------------------------------------------
      integer :: funit 
     
      if( .not. MPI_Process% isRoot ) return

      funit = UnusedUnit()

      open(funit,file='RESULTS/'//trim(OBB(STLNum)% FileName)//'_time.dat', action = "write" , access = "append" , status = "unknown")

      write(funit,'(1g15.6)') t

      close(funit)
   
   end subroutine WriteTimeFile

   subroutine TECtriangle_3points( funit, NumOfObjs )
      
      implicit none
      !-arguments------------------------------------
      integer, intent(in) :: funit, NumOfObjs
      !-local-variables------------------------------
      integer :: i, index
   
      index = 0

      do i = 1, NumOfObjs/3
         write(funit,'(i0,a,i0,a,i0)') index + 1,' ', index + 2,' ', index + 3
         index = index + 3
      end do
   
   end subroutine TECtriangle_3points

end module IBMClass
