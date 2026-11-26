!
!//////////////////////////////////////////////////////
!
!      Class cotaining routines for adapting polynomial orders based on the Truncation Error.
!        -> The adaptation procedure is performed with  
!           truncation error estimations (TruncationErrorClass.f90).
!        -> The current implementation is only valid for shared memory 
!           parallelization (OpenMP). TODO: Update to MPI.
!
!////////////////////////////////////////////////////////////////////////
!
module pAdaptationClassTE
   use SMConstants
#ifdef NAVIERSTOKES
   use PhysicsStorage                  , only: CTD_IGNORE_MODE, flowIsNavierStokes
#else
   use PhysicsStorage                  , only: CTD_IGNORE_MODE
#endif
   use ElementClass
   use DGSEMClass                      , only: DGSem, ComputeTimeDerivative_f, MaxTimeStep, ComputeMaxResiduals
   use TruncationErrorClass            , only: NMINest, TruncationError_t, ISOLATED_TE, NON_ISOLATED_TE, GenerateExactTEmap, OLD_TE, NEW_TE
   use FTValueDictionaryClass          , only: FTValueDictionary
   use StorageClass
   use FileReadingUtilities            , only: getFileName, getIntArrayFromString, getCharArrayFromString, GetRealValue, GetIntValue
   use ParamfileRegions                , only: readValueInRegion
   use Utilities                       , only: toLower
   use ReadMeshFile                    , only: NumOfElemsFromMeshFile
   use ExplicitMethods                 , only: TakeRK3Step
   use InterpolationMatrices           , only: Interp3DArrays
   use MultiTauEstimationClass         , only: MultiTauEstim_t
   use StopwatchClass                  , only: Stopwatch
   use pAdaptationClass                
   implicit none
   
#include "Includes.h"
   private
   public pAdaptationTE_t
   
   integer, parameter :: MAX_STEPS_SMOOTHING = 1000
   
   integer, parameter :: SMOOTH_RK3 = 0
   integer, parameter :: SMOOTH_FAS = 1
   !---------------------------------
   ! Type for storing the source term
   !---------------------------------
   type :: SourceStorage_t
      integer                    :: Nold(NDIM)
      real(kind=RP), allocatable :: S(:,:,:,:)
   end type SourceStorage_t
   
   !--------------------------------------------------
   ! Main type for performing a p-adaptation procedure
   !--------------------------------------------------
   type, extends(pAdaptation_t) :: pAdaptationTE_t
      real(kind=RP)                     :: reqTEc                          ! Requested truncation error for coarsening
      logical                           :: TEFiles = .FALSE.               ! Write truncation error files?
      logical                           :: increasing      = .FALSE.       ! Performing an increasing adaptation procedure?
      logical                           :: ErrorEstimFromFiles = .FALSE.   ! Read error estimation from files?
      character(len=LINE_LENGTH)        :: EstimFilesName
      integer                           :: EstimFilesNumber(2)
      type(MultiTauEstim_t)             :: MultiTauEstim
      integer                           :: TruncErrorType                  ! Truncation error type (either ISOLATED_TE or NON_ISOLATED_TE)
      integer                           :: TruncErrorForm                  ! Truncation error form (either OLD_TE or NEW_TE)
      integer, allocatable              :: maxNdecrease
      ! Variables for post-smoothing
      logical                           :: postSmoothing   = .FALSE.       ! Signals if a smoothing operation must be performed post p-adaptation
      logical                           :: Compute_dt      = .FALSE.
      integer                           :: postSmoothMethod
      real(kind=RP)                     :: postSmoothRes                   ! Post smoothing residual
      real(kind=RP)                     :: cfl, dcfl
      real(kind=RP)                     :: dt
      type(SourceStorage_t),allocatable :: Source(:)
      type(TruncationError_t), allocatable :: TE(:)         ! Truncation error for every element(:)
      
      contains
         ! Base class procedures
         procedure :: pAdaptation_Construct
         procedure :: pAdaptation_Destruct
         procedure :: pAdaptation_pAdapt
         ! Extended procedures
         procedure :: StoreQdot  => pAdaptation_storeQdot
         procedure :: postSmooth => pAdaptation_postSmooth
   end type pAdaptationTE_t
!
!  ----------------
!  Module variables
!  ----------------
!
   integer    :: NMIN(NDIM) = 1
!!    integer               :: dN_Inc = 3 
   integer    :: fN_Inc = 2
   integer    :: nelem           ! number of elements in mesh
   
   ! Here we define the input variables that can be changed after p-adaptation
   character(len=18), parameter :: ReplaceableVars(4) = (/'mg sweeps         ', &
                                                          'mg sweeps pre     ', &
                                                          'mg sweeps post    ', &
                                                          'mg sweeps coarsest'/)
   type(FTValueDictionary) :: ReplacedVariables
!========
 contains
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ROUTINES FOR ADAPTATION
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------
!  Routine for constructing the p-adaptator
!  ----------------------------------------
   subroutine pAdaptation_Construct(this, controlVariables, t0, mesh, adaptiveTimeStep)
      implicit none
      !--------------------------------------
      class(pAdaptationTE_t) , intent(inout) :: this             !>  P-Adaptator
      type(FTValueDictionary), intent(in)    :: controlVariables !<  Input values
      real(kind=RP)          , intent(in)    :: t0
      class(HexMesh)         , intent(inout) :: mesh
      type(adaptiveTimeStep_t), intent(inout):: adaptiveTimeStep
      !--------------------------------------
      ! For block reading
      character(LINE_LENGTH)         :: paramFile
      character(LINE_LENGTH)         :: in_label
      character(20*BC_STRING_LENGTH) :: confBoundaries
      character(LINE_LENGTH)         :: simulationType
      character(LINE_LENGTH)         :: R_Nmax, R_Nmin, R_TruncErrorType, R_TruncErrorForm, R_OrderAcrossFaces, replacedValue, R_mode, R_interval, R_pSmoothingMethod, R_EstimFiles, R_EstimFilesNum
      logical      , allocatable     :: R_increasing, R_TEFiles, reorganize_z, R_restart
      real(kind=RP), allocatable     :: TruncError, R_pSmoothing, R_cTruncError
      ! Extra vars
      integer                        :: i      ! Element counter
      integer                        :: no_of_overen_boxes
      !--------------------------------------
      
      this % Adapt  = pAdaptationIsDefined()
      
      if (this % Adapt) then
         this % Constructed = .TRUE.
      else
         this % Constructed = .FALSE.
         return
      end if
      
      nelem = NumOfElemsFromMeshFile( controlVariables % stringValueForKey("mesh file name", requestedLength = LINE_LENGTH) )  ! THIS DOES NOT WORK IN PARALLEL!!!
!
!     **************************************************
!     * p-adaptation is defined - Proceed to construct *     
!     **************************************************
!

!     Read block
!     **********
      write(in_label , '(A)') "#define p-adaptation"
      
      call get_command_argument(1, paramFile) !
      
      call readValueInRegion ( trim ( paramFile )  , "conforming boundaries"  , confBoundaries     , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "increasing"             , R_increasing       , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "truncation error"       , TruncError         , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "coarse truncation error", R_cTruncError      , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "write error files"      , R_TEFiles          , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "adjust nz"              , reorganize_z       , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "nmax"                   , R_Nmax             , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "nmin"                   , R_Nmin             , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "truncation error type"  , R_TruncErrorType   , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "truncation error formulation"  , R_TruncErrorForm   , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "order across faces"     , R_OrderAcrossFaces , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "max n decrease"         , this % maxNdecrease, in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "mode"                   , R_mode             , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "interval"               , R_interval         , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "restart files"          , R_restart          , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "post smoothing residual", R_pSmoothing       , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "post smoothing method"  , R_pSmoothingMethod , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "estimation files"       , R_EstimFiles       , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "estimation files number", R_EstimFilesNum    , in_label , "# end" )
      
!     Conforming boundaries?
!     ----------------------
      if ( confBoundaries /= "" ) then
         call getCharArrayFromString (confBoundaries,BC_STRING_LENGTH,this % conformingBoundaries)
         do i=1, size(this % conformingBoundaries)
            call toLower(this % conformingBoundaries(i))
         end do
      end if
      
!     Increasing adaptation?
!     ----------------------
      if ( allocated(R_increasing) ) then
         this % increasing = R_increasing
      end if
      
!     Truncation error threshold
!     --------------------------
      if ( allocated(TruncError) ) then
         this % reqTE = TruncError
      else
         error stop 'A truncation error must be specified for p-adapt'
      end if

!     Truncation error threshold for coarsening
!     -----------------------------------------
      if ( allocated(R_cTruncError) ) then
         this % reqTEc = R_cTruncError
      else
         this % reqTEc = this % reqTE
      end if
      
!     Write truncation error files?
!     -----------------------------
      if ( allocated(R_TEFiles) ) then
         this % TEFiles = R_TEFiles
      end if
      
!     Adjust Nz?
!     ----------
      if ( allocated(reorganize_z) ) then
         reorganize_Nz = reorganize_z
      end if
      
!     Nmax
!     ----
      if ( R_Nmax /= "" ) then
         this % NxyzMax = getIntArrayFromString(R_Nmax)
      else
         error stop 'Keyword Nmax is mandatory for p-adaptation'
      end if
      
!     Nmin -> If this is a p-nonconforming 3D case, it should be 2
!     ----
      if ( R_Nmin /= "" ) then
         NMIN = getIntArrayFromString(R_Nmin)
      end if
      
!     Truncation error type
!     ---------------------
      if ( R_TruncErrorType /= "" ) then
         call toLower (R_TruncErrorType)
         select case ( trim(R_TruncErrorType) )
            case ('isolated')
               this % TruncErrorType = ISOLATED_TE
            case ('non-isolated')
               this % TruncErrorType = NON_ISOLATED_TE
            case default
               write(STD_OUT,*) "Not recognized 'truncation error type'. Defaulting to isolated."
               this % TruncErrorType = ISOLATED_TE
         end select
      else
         this % TruncErrorType = ISOLATED_TE
      end if

!     Truncation error formulation, the selecion is between old and new (Andre's Rueda Ramirez thesis for more details)
!     ---------------------
      if ( R_TruncErrorForm /= "" ) then
         call toLower (R_TruncErrorForm)
         select case ( trim(R_TruncErrorForm) )
            case ('old')
               this % TruncErrorForm = OLD_TE
            case ('new')
               this % TruncErrorForm = NEW_TE
            case default
               write(STD_OUT,*) "Not recognized 'truncation error formulation'. Defaulting to old."
               this % TruncErrorForm = OLD_TE
         end select
      else
         this % TruncErrorForm = OLD_TE
      end if

      
!     Polynomial order jump
!     ---------------------
      if ( R_OrderAcrossFaces /= "" ) then
         select case ( trim (R_OrderAcrossFaces) )
            case ("n*2/3")
               GetOrderAcrossFace => NumberN23
            case default
               GetOrderAcrossFace => NumberN_1
         end select
      else
         GetOrderAcrossFace => NumberN_1
      end if
      
!     Adaptation mode: Steady(default) or unsteady
!     --------------------------------------------
      
      if ( R_mode == "" ) R_mode = "static"
      select case ( trim(R_mode) )
         case ("time")
            this % adaptation_mode = ADAPT_DYNAMIC_TIME
            this % time_interval   = GetRealValue(R_interval)
            this % iter_interval   = huge(this % iter_interval)
            this % nextAdaptationTime = t0   ! + this % time_interval
         case ("iteration")
            this % adaptation_mode = ADAPT_DYNAMIC_ITER
            this % time_interval   = huge(this % time_interval)
            this % iter_interval   = GetIntValue(R_interval)
         case ("static")
            this % adaptation_mode = ADAPT_STATIC
            this % time_interval   = huge(this % time_interval)
            this % iter_interval   = huge(this % iter_interval)
         case default
            WRITE(STD_OUT,*) 'Not recognized adaptation mode. Options are:'
            WRITE(STD_OUT,*) '   * time'
            WRITE(STD_OUT,*) '   * iteration'
            WRITE(STD_OUT,*) '   * steady'
      end select
      
!     Restart files
!     -------------
      if ( allocated(R_restart) ) then
         this % restartFiles = R_restart
      end if
      
!     Post smoothing residual and related variables
!     ---------------------------------------------
      if ( allocated(R_pSmoothing) ) then
         this % postSmoothing = .TRUE.
         this % postSmoothRes = R_pSmoothing
         
         select case (R_pSmoothingMethod)
            case('RK3')  ; this % postSmoothMethod = SMOOTH_RK3
            case('FAS')  ; this % postSmoothMethod = SMOOTH_FAS
            case default ; this % postSmoothMethod = SMOOTH_RK3
         end select
         
         allocate ( this % Source(nelem) )
      
         if (controlVariables % containsKey("cfl")) then
#ifdef FLOW
            this % Compute_dt = .TRUE.
            this % cfl        = controlVariables % doublePrecisionValueForKey("cfl")
#if defined(NAVIERSTOKES)
            if (flowIsNavierStokes) then
               if (controlVariables % containsKey("dcfl")) then
                  this % dcfl       = controlVariables % doublePrecisionValueForKey("dcfl")
               else
                  error stop '"cfl" and "dcfl", or "dt" keyword must be specified for the time integrator'
               end if
            end if
#endif
#elif defined(CAHNHILLIARD)
            print*, "Error, use fixed time step to solve Cahn-Hilliard equations"
            errorMessage(STD_OUT)
            error stop
#endif
         elseif (controlVariables % containsKey("dt")) then
            this % Compute_dt = .FALSE.
            this % dt         = controlVariables % doublePrecisionValueForKey("dt")
         else
            error stop '"cfl" or "dt" keyword must be specified for the time integrator'
         end if
      end if
      
!     Truncation error estimation in files
!     ------------------------------------
      if ( R_EstimFiles /= "" ) then
         this % ErrorEstimFromFiles = .TRUE.
         this % EstimFilesName = trim (R_EstimFiles)
         
         if ( R_EstimFilesNum == "" ) then
            error stop 'pAdaptation :: The user must provide the file numbers to read the tau-estimation'
         end if
         
         this % EstimFilesNumber = getIntArrayFromString(R_EstimFilesNum)
         
         call this % MultiTauEstim % constructForPAdaptator ( this % TruncErrorType, this % NxyzMax, nelem, this % EstimFilesName )
      end if
      
!     Check replaceable control variables
!     -----------------------------------
      call ReplacedVariables % initWithSize(16)
      
      do i=1, size(ReplaceableVars)
         call readValueInRegion ( trim ( paramFile )  , "padapted " // trim(ReplaceableVars(i)), replacedValue, in_label , "# end" )
         if (replacedValue == "") cycle
         call ReplacedVariables % addValueForKey(trim(replacedValue),TRIM(ReplaceableVars(i)))
      end do
      
!
!     Stopwatch events
!     ****************
!  
      call Stopwatch % CreateNewEvent("pAdapt: Error estimation")
      call Stopwatch % CreateNewEvent("pAdapt: PolOrder selection")
      call Stopwatch % CreateNewEvent("pAdapt: Adaptation")
      
      
!
!     Some things are read from the control file
!     ******************************************
!      
      this % solutionFileName = trim(getFileName(controlVariables % stringValueForKey("solution file name", requestedLength = LINE_LENGTH)))
      this % saveGradients    = controlVariables % logicalValueForKey("save gradients with solution")
      this % saveSensor       = controlVariables % logicalValueForKey("save sensor with solution")
      simulationType = controlVariables % StringValueForKey("simulation type",LINE_LENGTH)
      call toLower(simulationType)
      if ( trim(simulationType) == 'time-accurate' ) this % UnSteady = .TRUE.
      
      
!     Construct the truncation error
!     ******************************    
      
      allocate (this % TE(nelem))
      do i = 1, nelem
         call this % TE(i) % construct(this % NxyzMax)
      end do
      
      
!     Adaptation overenriching
!     *************************
      
      call getNoOfOverEnrichingBoxes(no_of_overen_boxes)
         
      if (no_of_overen_boxes > 0) then
         allocate ( this % overenriching(no_of_overen_boxes) )
         
         do i = 1, no_of_overen_boxes
            call this % overenriching(i) % initialize (i)
         end do
      end if

      safedeallocate(R_increasing)
      safedeallocate(R_TEFiles)
      safedeallocate(reorganize_z)
      safedeallocate(R_restart)
      safedeallocate(TruncError)
      safedeallocate(R_pSmoothing)
      safedeallocate(R_cTruncError)

   end subroutine pAdaptation_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------
!  Routine for destructing the p-adaptator
!  ----------------------------------------
   subroutine pAdaptation_Destruct(this)
      implicit none
      !--------------------------------------
      class(pAdaptationTE_t) :: this
      !--------------------------------------
      integer                :: iEl
      !--------------------------------------
      
      ! Truncation error
      do iEl = 1, nelem
         call this % TE(iEl) % destruct()
      end do
      
      deallocate (this % TE)
      
      safedeallocate  (this % conformingBoundaries)
      safedeallocate  (this % overenriching)
      
      call this % MultiTauEstim % destruct
   end subroutine pAdaptation_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------------
!  Main routine for adapting the polynomial order in all elements based on 
!  the truncation error estimation
!  ------------------------------------------------------------------------
   subroutine pAdaptation_pAdapt(this,sem,itera,t, computeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables, adaptiveTimeStep)
      use AnisFASMultigridClass
      implicit none
      !-arguments----------------------------
      class(pAdaptationTE_t)     :: this              !<> Adaptation class
      type(DGSem)                :: sem               !<> sem
      integer                    :: itera             !<  iteration
      real(kind=RP)              :: t                 !< time!!
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivativeIsolated
      type(FTValueDictionary)    :: controlVariables  !<> Input variables (that can be modified depending on the user input)
      type(adaptiveTimeStep_t)   :: adaptiveTimeStep 
      !-local-variables----------------------
      integer                    :: eID               !   Element counter
      integer                    :: Dir               !   Direction
      integer                    :: NNew(3,nelem)     !   New polynomial orders of mesh (after adaptation!)
      integer                    :: Error(3,nelem)    !   Stores (with ==1) elements where the truncation error has a strange behavior of the truncation error (in every direction)
      integer                    :: Warning(nelem)    !   Stores (with ==1) elements where the specified truncation error was not achieved
      integer, save              :: Stage = 0         !   Stage of p-adaptation for the increasing method
      CHARACTER(LEN=LINE_LENGTH) :: newInput          !   Variable used to change the input in controlVariables after p-adaptation 
      character(len=LINE_LENGTH) :: RegfileName
      integer                    :: i             !   Counters
      TYPE(AnisFASMultigrid_t)   :: AnisFASpAdaptSolver
      character(len=LINE_LENGTH) :: AdaptedMeshFile
      logical                    :: last
      !--------------------------------------
#if defined(NAVIERSTOKES)
      
      write(STD_OUT,*)
      write(STD_OUT,*)
      select case (this % TruncErrorType)
         case (ISOLATED_TE)
            write(STD_OUT,'(A)') '****     Performing p-Adaptation with isolated truncation error estimates      ****'
         case (NON_ISOLATED_TE)
            write(STD_OUT,'(A)') '****     Performing p-Adaptation with non-isolated truncation error estimates      ****'
         case default
            error stop ':: Truncation error type was not defined'
      end select
      write(STD_OUT,*)
      
      Error = 0
      Warning= 0
      Stage = Stage + 1
      
!
!     -------------------------------------------
!     Get exact truncation error map if requested
!     -------------------------------------------
!
      if ( controlVariables % containsKey("get exact temap elem") ) then
         eID = controlVariables % integerValueForKey("get exact temap elem")
         call GenerateExactTEmap(sem, NMIN, this % NxyzMax, t, computeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables, eID, TruncErrorType=this % TruncErrorType, TruncErrorForm =this % TruncErrorForm)
      end if
      
!
!     --------------------------------------
!     Write pre-adaptation mesh and solution
!     --------------------------------------
!
      if (this % restartFiles) then
         write(AdaptedMeshFile,'(A,A,I2.2,A)')  trim( this % solutionFileName ), '_pre-Adapt_Stage_', Stage, '.hsol'
         call sem % mesh % Export(AdaptedMeshFile)         
         call sem % mesh % SaveSolution(itera,t,trim(AdaptedMeshFile),this % saveGradients,this % saveSensor)
      end if
      
!
!     -------------------------------------------------------------------------
!     Estimate the truncation error using the anisotropic multigrid (if needed)
!     -------------------------------------------------------------------------
!
      
      call Stopwatch % Start("pAdapt: Error estimation")
      if (.not. this % ErrorEstimFromFiles) then

         CALL AnisFASpAdaptSolver % construct(controlVariables,sem,estimator=.TRUE.,NMINestim = NMINest)
         CALL AnisFASpAdaptSolver % solve(itera,t,computeTimeDerivative,ComputeTimeDerivativeIsolated,this % TE, TEType = this % TruncErrorType, TEForm = this % TruncErrorForm)
         CALL AnisFASpAdaptSolver % destruct
         
         ! Write truncation error files if specified
         if (this % TEFiles) then
            do eID = 1, nelem
               write(RegfileName,'(A,I8.8,A)') 'TauEstimation/Elem_',eID,'.dat'
               call this % TE(eID) % ExportToFile(RegfileName)
            end do
         end if
         
      end if
      call Stopwatch % Pause("pAdapt: Error estimation")
!
!     ------------------------
!     Post-smoothing procedure
!     ------------------------
!
      if ( this % postSmoothing) then
         call this % StoreQdot(sem, t, ComputeTimeDerivative)
      end if
!
!     -------------------------------------------------------------
!     Find the polynomial order that fulfills the error requirement
!     -------------------------------------------------------------
!
      
      call Stopwatch % Start("pAdapt: PolOrder selection")
      
!$omp parallel do schedule(runtime)
      do eID = 1, nelem
         call pAdaptation_pAdaptTE_SelectElemPolorders (this, sem % mesh % elements(eID), NNew(:,eID), error(:,eID), warning(eID) )
      end do
!$omp end parallel do
      
      call Stopwatch % Pause("pAdapt: PolOrder selection")
      
!
!     ----------------------------------------------------------------------------
!     In case of increasing adaptator, modify polynomial orders according to stage
!      And decide if it is necessary to continue adapting
!     ----------------------------------------------------------------------------
!
      if (this % increasing) then
!!          Stage = Stage + 1
!!          NInc = NInc + dN_Inc
         NInc = NInc * fN_Inc
         
         if (MAXVAL(NNew) > NInc) then
            where(NNew > NInc) NNew = NInc
         else
            this % Adapt = .FALSE.
         end if
         
      else ! Only adapt once
         this % Adapt = .FALSE.
      end if
      
!
!     ----------------------------
!     Overenrich specified regions
!     ----------------------------
!
      call OverEnrichRegions(this % overenriching,sem % mesh,NNew, this % NxyzMax, NMIN)
!
!     ---------------------------------------------------
!     Restrict polynomial order decrease in every element
!     ---------------------------------------------------
!
      if ( allocated(this % maxNdecrease) ) then
         do eID=1, nelem
            associate (e => sem % mesh % elements(eID) )
            do dir=1, NDIM
               
               if ( NNew(dir,eID) < (e % Nxyz(dir) - this % maxNdecrease) ) NNew(dir,eID) = e % Nxyz(dir) - this % maxNdecrease
               
            end do
            end associate
         end do
      end if
!
!     ------------------------------
!     Restrict polynomial order jump
!     ------------------------------
!
      last = .FALSE.
      do while (.not. last)
         last = .TRUE.
         call this % makeBoundariesPConforming(sem % mesh,NNew,last)
         call ReorganizePolOrders(sem % mesh % faces,NNew,last)
      end do
!
!     -----------------------
!     Plot files if requested
!     -----------------------
!
!~      if (this % PlotInfo) call this % plot(sem,TE,Stage,NNew,Error,Warning)
!
!     ----------------------------------
!     Adapt sem to new polynomial orders
!     ----------------------------------
!
!~      call Stopwatch % Pause("Solver")
!~      call Stopwatch % Start("Preprocessing")
      
      call Stopwatch % Start("pAdapt: Adaptation")
      
      call sem % mesh % pAdapt (NNew, controlVariables)
      
      call Stopwatch % Pause("pAdapt: Adaptation")
      
      ! Reconstruct probes
      do i=1, sem % monitors % no_of_probes
         call sem % monitors % probes(i) % Initialization (sem % mesh, i, trim(sem % monitors % solution_file), .FALSE.)
      end do
      
!~      call Stopwatch % Pause("Preprocessing")
!~      call Stopwatch % Start("Solver")
      
!
!     ------------------------------------------------------
!     Perform the post p-adaptation smoothing (if requested)
!     ------------------------------------------------------
!
      if ( this % postSmoothing ) then
         call this % postSmooth(sem, t, ComputeTimeDerivative, controlVariables)
      end if
!
!     ---------------------------------------------------
!     Write post-adaptation mesh, solution and order file
!     ---------------------------------------------------
!
      if ( this % UnSteady) then
         write(AdaptedMeshFile,'(A,A,I10.10,A)')  trim( this % solutionFileName ), '_', itera+1, '.hsol'
      else
         write(AdaptedMeshFile,'(A,A,I2.2,A)')  trim( this % solutionFileName ), '_p-Adapted_Stage_', Stage, '.hsol'
      end if
      
      call sem % mesh % Export(AdaptedMeshFile)
      call sem % mesh % ExportOrders(AdaptedMeshFile)
      
      if (this % restartFiles) call sem % mesh % SaveSolution(itera,t,trim(AdaptedMeshFile),this % saveGradients,this % saveSensor)
      
!
!     ------------------------------------------------
!     Rewrite controlVariables according to user input
!     ------------------------------------------------
!
      do i = 1, size(ReplaceableVars)
         if ( ReplacedVariables % containsKey( trim(ReplaceableVars(i)) ) ) then
            newInput = ReplacedVariables % stringValueForKey(trim(ReplaceableVars(i)), requestedLength = LINE_LENGTH)
            call controlVariables % removeObjectForKey( trim(ReplaceableVars(i)) )
            call controlVariables % addValueForKey( TRIM(newInput) , trim(ReplaceableVars(i)) )
         end if
      end do
!
!     ----------------
!     Update residuals
!     ----------------
!
      call ComputeTimeDerivative(sem % mesh, sem % particles, t, CTD_IGNORE_MODE)

!
!     ---------------------
!     Update HO arrays
!     ---------------------
!
      call sem % mesh % UpdateHOArrays()
      
      write(STD_OUT,*) '****    p-Adaptation done, DOFs=', SUM((NNew(1,:)+1)*(NNew(2,:)+1)*(NNew(3,:)+1)), '****'

#endif
   end subroutine pAdaptation_pAdapt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------------------
!  pAdaptation_pAdaptTE_SelectElemPolorders:
!  Select the polynomial orders for one element based on the truncation error estimation
!  -------------------------------------------------------------------------------------
   subroutine pAdaptation_pAdaptTE_SelectElemPolorders (this, e, NNew, error, warning)
      implicit none
      !-arguments----------------------------------
      class(pAdaptationTE_t), intent(inout) :: this        !<> Adaptation class: Only needed as inout to extrapolate TE (change?)
      type(Element)       , intent(in)    :: e
      integer             , intent(out)   :: NNew(NDIM)
      integer             , intent(out)   :: error(NDIM) !>  0 if everything is right, 1 if no extrapolation was possible
      integer             , intent(out)   :: warning     !>  0 if everything is right, 1 if reqTE was not achieved in TEmap
      !-local-variables----------------------------
      logical        :: notenough(NDIM)  ! Not enough points to extrapolate
      integer        :: Pxyz(3)     ! Estimation polynomial order
      integer        :: P_1 (3)     ! Estimation polynomial order - 
      integer        :: i, j, k     ! Coordinate counters
      integer        :: dir         ! Coordinate direction
      integer        :: L(NDIM)     ! Temporary variable to store polynomial orders
      integer        :: DOFs        ! DOFs of a configuration
      integer        :: NewDOFs     ! DOFs of the currently selected configuration
      real(kind=RP)  :: TE_L        ! Truncation error of the currently selected configuration
      real(kind=RP)  :: TEmap(NMINest:this % NxyzMax(1),NMINest:this % NxyzMax(2),NMINest:this % NxyzMax(3))
      !--------------------------------------------
      
!     Initialization of output
!     ------------------------
      error   = 0
      warning = 0
      
!     Initialization of P and P-1
!     ---------------------------
      Pxyz = e % Nxyz ! Assumed as Nxyz (not necessarily true for static p-adaptation, but makes no difference - P_1 is corrected)
      P_1  = Pxyz - 1
      
      ! Correction of P-1
      do dir=1, NDIM
         if ( (P_1(dir) < NMIN(dir)) .and. (P_1(dir) < NMINest+1) ) P_1(dir) = NMIN(dir)
      end do
      
!     Initialization of TEmap and NNew
!     -------------------------------
      TEmap = 0._RP
      NNew = -1 ! Initialized to negative value
!
!     --------------------------------------------------------------------------------
!     If ErrorEstimFromFiles, get complete TEmap, otherwise compute the inner map only
!     --------------------------------------------------------------------------------
!
      if (this % ErrorEstimFromFiles) then
         call this % MultiTauEstim % GetTEmap(this % EstimFilesNumber, e % globID,this % NxyzMax,NMIN,P_1,TEmap)
      else
         do k = NMINest, P_1(3) ; do j = NMINest, P_1(2) ; do i = NMINest, P_1(1) 
                  ! 1. Generate a TEmap entry
                  TEmap(i,j,k) = this % TE(e % eID) % Dir(1) % maxTE(i) + &  !xi   contribution
                                 this % TE(e % eID) % Dir(2) % maxTE(j) + &  !eta  contribution
                                 this % TE(e % eID) % Dir(3) % maxTE(k)      !zeta contribution
         end do                 ; end do                 ; end do
      end if
      
!
!     -----------------------------------------------------------------------
!     Check if the desired TE can be obtained using the "inner" map (N_i<P_i)
!        Accomplished in 1 merged steps:
!           1. Check every entry of the map 
!           2. select the one that
!              fulfills the requirement with the lowest DOFs
!     ----------------------------------------------------------------------
!
      NewDOFs = (P_1(1) + 1) * (P_1(2) + 1) * (P_1(3) + 1) ! Initialized with maximum number of DOFs (for N<P)
      TE_L = huge(1._RP)
      
      do k = NMIN(3), P_1(3) ; do j = NMIN(2), P_1(2) ; do i = NMIN(1), P_1(1)

         ! 2. Check if it fulfills the requirement
         if (TEmap(i,j,k) < this % reqTEc) then
            DOFs = (i+1) * (j+1) * (k+1)
            
            !  3. Select the entry if it minimizes the DOFs
            if (DOFs < NewDOFs) then
               NNew = [i,j,k]
               NewDOFs = DOFs
               TE_L    = TEmap(i,j,k)
            elseif ( (DOFs == NewDOFs) .and. (TEmap(i,j,k) < TE_L) ) then
               NNew = [i,j,k]
               NewDOFs = DOFs
               TE_L    = TEmap(i,j,k)
            end if
         end if
      end do                 ; end do                 ; end do
!
!     -----------------------------------------------------------------------
!     Extra check:
!        When NMIN > P_1, it is not always necessary to extrapolate.
!        If Tau(P_1) < threshold, then NMIN should fulfill it too
!     -----------------------------------------------------------------------
!
      if ( any(NMIN > P_1) .and. all( P_1 >= NMINest ) ) then
         do dir=1, NDIM
            if (P_1(dir) <= NMIN(dir)) then
               L(dir) = P_1(dir)
            else
               L(dir) = NMIN(dir)
            end if
         end do
         
         if (TEmap(L(1),L(2),L(3)) < this % reqTEc) NNew = NMIN
      end if
            
!
!     -----------------------------------------------------------------------
!     If the desired TE could NOT be achieved using the inner map (N_i<P_i),
!     find it with a higher order
!        Accomplished in 3 steps:
!           1. Perform regression analysis and extrapolate the decoupled TE
!           2. Obtain the extended TE map for higher orders
!           3. Check every entry of the extended map
!           4. Select the N>P that fulfills the requirement with lowest DOFs
!     -----------------------------------------------------------------------
!
      if (any(NNew<NMIN)) then
         
         if (.not. this % ErrorEstimFromFiles) then
!
!              1. Regression analysis
!              ----------------------
            do Dir = 1, 3
               call this % TE(e % eID) % ExtrapolateInOneDir(P_1       = P_1(Dir)             , & 
                                                             NMax      = this % NxyzMax(Dir), &
                                                             Dir       = Dir                  , &
                                                             notenough = notenough(Dir)   , &
                                                             error     = error(Dir)   )
            end do
         end if
         
         ! If the truncation error behaves as expected, continue, otherwise skip steps 2-3-4 and select maximum N. TODO: Change? extrapolation can still be done in some directions...
         if (all(error < 1)) then
!
!           2. Generate outer map
!           -> ">=" Pxyz(dir) can be changed by ">" to avoid enriching 2
!           ---------------------------------------------------------------
            if (.not. this % ErrorEstimFromFiles) then
               do k = NMIN(3), this % NxyzMax(3) ; do j = NMIN(2), this % NxyzMax(2) ; do i = NMIN(1), this % NxyzMax(1)
                        ! cycle if it is not necessary/possible to compute the TEmap
                        if (k <= P_1(3) .AND. j <= P_1(2) .AND. i <= P_1(1)) cycle ! This is the inner map (already computed)
                        if ( (notenough(1) .AND. i >= Pxyz(1)) .or. & ! The regression was not possible in xi   (too few points), hence only checking with known values
                             (notenough(2) .AND. j >= Pxyz(2)) .or. & ! The regression was not possible in eta  (too few points), hence only checking with known values
                             (notenough(3) .AND. k >= Pxyz(3)) ) then ! The regression was not possible in zeta (too few points), hence only checking with known values
                           TEmap(i,j,k) = huge(1._RP)
                        else
                           TEmap(i,j,k) = this % TE(e % eID) % Dir(1) % maxTE(i) + &  !x contribution
                                          this % TE(e % eID) % Dir(2) % maxTE(j) + &  !y contribution
                                          this % TE(e % eID) % Dir(3) % maxTE(k)      !z contribution
                        end if
               end do                            ; end do                              ; end do
            end if
!
!           3. Check the outer map
!           ----------------------
            NewDOFs = (this % NxyzMax(1) + 1) * (this % NxyzMax(2) + 1) * (this % NxyzMax(3) + 1) ! Initialized as maximum DOFs possible
            TE_L = huge(1._RP)
            
            do k = NMIN(3), this % NxyzMax(3) ; do j = NMIN(2), this % NxyzMax(2) ; do i = NMIN(1), this % NxyzMax(1)
               
               ! cycle if it is not necessary/possible to compute the TEmap
               if (k <= P_1(3) .AND. j <= P_1(2) .AND. i <= P_1(1)) cycle ! This is the inner map (already checked)
               
               ! Check if TE was achieved
               if (TEmap(i,j,k) < this % reqTE) then
                  DOFs = (i+1) * (j+1) * (k+1)
                  ! 4. Select the entry if it minimizes the DOFs
                  if (DOFs < NewDOFs) then
                     NNew = [i,j,k]
                     NewDOFs = DOFs
                     TE_L    = TEmap(i,j,k)
                  elseif ( (DOFs == NewDOFs) .and. (TEmap(i,j,k) < TE_L) ) then
                     NNew = [i,j,k]
                     NewDOFs = DOFs
                     TE_L    = TEmap(i,j,k)
                  end if
               end if
            end do                            ; end do                            ; end do
         else
            write(STD_OUT,'(A)')  'p-Adaptation ERROR: Unexpected behavior of truncation error.'
            write(STD_OUT,'(A,I0)')  ' -> Element: ', e % globID
            write(STD_OUT,'(A,3I3)') ' -> Direction (1,2,3): ', error
            write(STD_OUT,'(A,3I3)') ' -> Using maximum polynomial order, N = ', this % NxyzMax
            NNew = this % NxyzMax
         end if
      end if
      
!~      if (e % eID==1) call PrintTEmap(NMIN,TEmap,1,this % solutionFileName)
!
!     ---------------------------------------------------------------------------
!     If the maximum polynomial order was not found, select the maximum available
!     ---------------------------------------------------------------------------
!
      if (any(NNew < NMIN)) then
         write(STD_OUT,'(A,I0)')  'p-Adaptation WARNING: Truncation error not achieved within limits.'
         write(STD_OUT,'(A,I0)')  ' -> Element: ', e % globID
         write(STD_OUT,'(A,3I3)') ' -> Using max polynomial order instead, N = ', this % NxyzMax
         warning = 1
         NNew = this % NxyzMax
      end if
      
   end subroutine pAdaptation_pAdaptTE_SelectElemPolorders
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------------------
!  Subroutine to post-smooth the solution after p-adaptation
!  ---------------------------------------------------------
   subroutine pAdaptation_storeQdot(this, sem, t, ComputeTimeDerivative)
      implicit none
      !-arguments--------------------------------------------------
      class(pAdaptationTE_t), intent(inout) :: this            !<> Adaptation class
      type(DGSem)         , intent(inout) :: sem
      real(kind=RP)       , intent(in)    :: t
      procedure(ComputeTimeDerivative_f)  :: ComputeTimeDerivative
      !-local-variables--------------------------------------------
      integer :: eID
      integer :: N(NDIM)
      integer :: nEqn
      !------------------------------------------------------------
      
      nEqn = size(sem % mesh % storage % elements(1) % Qdot,1)
      
      call ComputeTimeDerivative(sem % mesh, sem % particles, t, CTD_IGNORE_MODE)
      
      do eID = 1, sem % mesh % no_of_elements
         N = sem % mesh % elements(eID) % Nxyz
         
         this % Source(eID) % Nold = N
         allocate (this % Source(eID) % S ( nEqn, N(1), N(2), N(3) ) )
         
         this % Source(eID) % S = - sem % mesh % storage % elements(eID) % Qdot
      end do
      
   end subroutine pAdaptation_storeQdot
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------------------
!  Subroutine to post-smooth the solution after p-adaptation
!  ---------------------------------------------------------
   subroutine pAdaptation_postSmooth(this, sem, t, ComputeTimeDerivative, controlVariables)
      use FASMultigridClass, only: FASMultigrid_t
      implicit none
      !-arguments--------------------------------------------------
      class(pAdaptationTE_t), intent(inout) :: this            !<> Adaptation class
      type(DGSem)         , intent(inout) :: sem
      real(kind=RP)       , intent(in)    :: t
      procedure(ComputeTimeDerivative_f)  :: ComputeTimeDerivative
      type(FTValueDictionary), intent(in) :: controlVariables
      !-local-variables--------------------------------------------
      integer       :: k
      integer       :: nEqn
      real(kind=RP),allocatable :: maxResidual(:)
      type(FASMultigrid_t) :: FASSolver
      !------------------------------------------------------------
      
      nEqn = size(sem % mesh % storage % elements(1) % Qdot,1)
      
      allocate ( maxResidual(nEqn) )
      
!     Assign source term
!     ------------------
#if defined(NAVIERSTOKES) || defined(INCNS)
      do k=1, sem % mesh % no_of_elements
         
         call Interp3DArrays  (  nEqn, &
                                 this % Source(k) % Nold        , this % Source(k) % S, &
                                 sem % mesh % elements(k) % Nxyz, sem % mesh % storage % elements(k) % S_NS )
         
      end do
#endif
      
      if (this % postSmoothMethod == SMOOTH_FAS) call FASSolver % construct(controlVariables,sem)
      
!     Smooth solution
!     ---------------
      do k=1, MAX_STEPS_SMOOTHING
         if ( this % Compute_dt ) call MaxTimeStep( self=sem, cfl=this % cfl, dcfl=this % dcfl , MaxDt=this % dt)
            
         select case (this % postSmoothMethod)
            case(SMOOTH_RK3)
               call TakeRK3Step( sem % mesh, sem % particles, t, this % dt, ComputeTimeDerivative )
            case(SMOOTH_FAS)
               call FASSolver % solve(0, t, ComputeTimeDerivative, ComputeTimeDerivative)
         end select
         
         maxResidual = ComputeMaxResiduals(sem % mesh)
         if (maxval(maxResidual) <= this % postSmoothRes )  then
            print*, 'post-smoothed in', k, 'iterations'
            exit
         end if
      end do
      
!     Remove source term
!     ------------------
#if defined(NAVIERSTOKES) || defined(INCNS)
      do k=1, sem % mesh % no_of_elements
         sem % mesh % storage % elements(k) % S_NS = 0._RP
         deallocate ( this % Source(k) % S )
      end do
#endif
      
      if (this % postSmoothMethod == SMOOTH_FAS) call FASSolver % destruct
   end subroutine pAdaptation_postSmooth   
   
end module pAdaptationClassTE