!
!//////////////////////////////////////////////////////
!
!      Class cotaining routines for adapting polynomial orders based on Reinforcement Learning.
!        -> The adaptation procedure is performed with a RL agent trined with a Value Iteration algorithm.
!        -> The current implementation is compatible with OpenMP and MPI.
!
!////////////////////////////////////////////////////////////////////////
!
module pAdaptationClassRL
   use SMConstants
#ifdef NAVIERSTOKES
   use PhysicsStorage                  , only: CTD_IGNORE_MODE, flowIsNavierStokes
   use FluidData                       , only: thermodynamics
#elif defined(MULTIPHASE)
   use PhysicsStorage                  , only: CTD_IGNORE_MODE, IMP
#else
   use PhysicsStorage                  , only: CTD_IGNORE_MODE
#endif
   use ElementClass
   use DGSEMClass                      , only: DGSem, ComputeTimeDerivative_f, MaxTimeStep, ComputeMaxResiduals
   use FTValueDictionaryClass          , only: FTValueDictionary
   use StorageClass
   use MPI_Process_Info
   use FileReadingUtilities            , only: getFileName, getIntArrayFromString, getCharArrayFromString, GetRealValue, GetIntValue
   use ParamfileRegions                , only: readValueInRegion
   use Utilities                       , only: toLower
   use ReadMeshFile                    , only: NumOfElemsFromMeshFile
   use ExplicitMethods                 , only: TakeRK3Step
   use InterpolationMatrices           , only: Interp3DArrays
   use MultiTauEstimationClass         , only: MultiTauEstim_t
   use StopwatchClass                  , only: Stopwatch
   use pAdaptationClass      
   use ReinforcementLearning           , only: pAdaptationAgent_t   
   
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none
   
#include "Includes.h"
   private
   public pAdaptationRL_t
   
   !--------------------------------------------------
   ! Main type for performing a p-adaptation procedure
   !--------------------------------------------------
   type, extends(pAdaptation_t) :: pAdaptationRL_t
      type(pAdaptationAgent_t) :: agent
      real(kind=RP)            :: threshold
      real(kind=RP)            :: tol = 1e-2_RP
      logical                  :: error_estimation = .false.
      logical                  :: avg_error_type = .true. !True for average error, false for max error
      integer                  :: error_variable !1:u, 2:v, 3:w, 4:rho*u, 5:rho*v, 6:rho*w, 7:p (only for Navier-Stokes), 8:rho
	  integer                  :: pJump
      logical                  :: acoustics = .false.
      real(kind=RP)            :: acoustic_tol = 1e-4_RP
      integer                  :: acoustic_variable !7:p, 8:rho
      real(kind=RP)            :: acoustic_distance = 1_RP
      real(kind=RP)            :: observer(NDIM)
      character(len=BC_STRING_LENGTH), allocatable :: acoustic_sources(:)
      
      contains
         ! Base class procedures
         procedure :: pAdaptation_Construct
         procedure :: pAdaptation_Destruct
         procedure :: pAdaptation_pAdapt
   end type pAdaptationRL_t
!
!  ----------------
!  Module variables
!  ----------------
!
   integer    :: NMIN(NDIM) = 1

!========
 contains
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------
!  Routine for constructing the p-adaptator
!  ----------------------------------------
   subroutine pAdaptation_Construct(this, controlVariables, t0, mesh, adaptiveTimeStep)
      use SurfaceMesh, only: surfacesMesh
      implicit none
      !--------------------------------------
      class(pAdaptationRL_t) , intent(inout) :: this             !>  P-Adaptator
      type(FTValueDictionary), intent(in)    :: controlVariables !<  Input values
      real(kind=RP)          , intent(in)    :: t0
      class(HexMesh)         , intent(inout) :: mesh
      type(adaptiveTimeStep_t), intent(inout):: adaptiveTimeStep
      !--------------------------------------
      ! For block reading
      character(LINE_LENGTH)         :: paramFile
      character(LINE_LENGTH)         :: in_label
      character(LINE_LENGTH)         :: agentFile
      character(LINE_LENGTH)         :: simulationType
      character(20*BC_STRING_LENGTH) :: confBoundaries, R_acoustic_sources
      character(LINE_LENGTH)         :: R_Nmax, R_Nmin, R_OrderAcrossFaces, replacedValue, R_mode, R_interval, R_Jump, cwd, R_ErrorType, R_ErrorVariable, R_observer, R_acoustic_variable
      logical      , allocatable     :: R_increasing, reorganize_z, R_restart, R_ErrorEstimation, R_acoustics
      real(kind=RP), allocatable     :: R_tolerance, R_threshold, R_acoustic_tol, R_acoustic_distance
      ! Extra vars
      integer                        :: i      ! Element counter
      integer                        :: no_of_overen_boxes
      logical                        :: adaptive_dt = .FALSE.
      !--------------------------------------
      
      this % Adapt  = pAdaptationIsDefined()
      
      if (this % Adapt) then
         this % Constructed = .TRUE.
      else
         this % Constructed = .FALSE.
         return
      end if
 
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
      call readValueInRegion ( trim ( paramFile )  , "tolerance"              , R_tolerance        , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "adjust nz"              , reorganize_z       , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "nmax"                   , R_Nmax             , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "nmin"                   , R_Nmin             , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "order across faces"     , R_OrderAcrossFaces , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "mode"                   , R_mode             , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "interval"               , R_interval         , in_label , "# end" )
	   call readValueInRegion ( trim ( paramFile )  , "max polynomial diff"    , R_Jump             , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "restart files"          , R_restart          , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "agent file"             , agentFile          , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "threshold"              , R_threshold        , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "error estimation"       , R_ErrorEstimation  , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "error type"             , R_ErrorType        , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "error variable"         , R_ErrorVariable    , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "acoustics"              , R_acoustics        , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "acoustic tolerance"     , R_acoustic_tol     , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "acoustic distance"      , R_acoustic_distance, in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "acoustic observer"      , R_observer         , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "acoustic sources"       , R_acoustic_sources , in_label , "# end" )
      call readValueInRegion ( trim ( paramFile )  , "acoustic variable"      , R_acoustic_variable, in_label , "# end" )
      
!     Conforming boundaries
!     ----------------------
      if ( confBoundaries /= "" ) then
         call getCharArrayFromString (confBoundaries,BC_STRING_LENGTH,this % conformingBoundaries)
         do i=1, size(this % conformingBoundaries)
            call toLower(this % conformingBoundaries(i))
         end do
      end if
      
!     Truncation error tolerance
!     --------------------------
      if ( allocated(R_tolerance) ) then
         this % reqTE = R_tolerance !Required for static adaptation mode
         this % tol = R_tolerance
      end if

!     Adaptation threshold
!     ---------------------
      if ( allocated(R_threshold) ) then
         this % threshold = R_threshold
      else
         this % threshold = 0.0_RP
      end if
      
!     Adjust Nz
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

!     Error estimator
!     ---------------
      if ( allocated(R_ErrorEstimation) ) then
         this % error_estimation = R_ErrorEstimation
      end if

      if (this % error_estimation) then
         if ( R_ErrorType /= "" ) then
            select case ( trim (R_ErrorType) )
            case ("avg")
               this % avg_error_type = .true.
            case ("max")
               this % avg_error_type = .false.
            case default
               WRITE(STD_OUT,*) 'Not recognized adaptation mode. Options are:'
               WRITE(STD_OUT,*) '   * avg'
               WRITE(STD_OUT,*) '   * max'
               error stop ' '
            end select
         end if

         if ( R_ErrorVariable /= "" ) then
            select case ( trim (R_ErrorVariable) )
            case ("u")
               this % error_variable = 1
            case ("v")
               this % error_variable = 2
            case ("w")
               this % error_variable = 3
            case ("rhou")
               this % error_variable = 4
            case ("rhov")
               this % error_variable = 5
            case ("rhow")
               this % error_variable = 6
            case ("p")
#ifdef NAVIERSTOKES
               this % error_variable = 7
#elif defined(MULTIPHASE)
               this % error_variable = 7
#endif
            case ("rho")
               this % error_variable = 8
            case default
               WRITE(STD_OUT,*) 'Not recognized error variable. Using u velocity by default.'
               this % error_variable = 1
            end select
         else
            error stop 'Keyword error variable is mandatory for error estimation'
         end if
      end if
      
!     Adaptation mode: Steady(default) or unsteady
!     --------------------------------------------
      
      if ( R_mode == "" ) R_mode = "static"
      select case ( trim(R_mode) )
         case ("time")
            this % adaptation_mode = ADAPT_DYNAMIC_TIME
            if ( R_interval /= "" ) then
               this % time_interval = GetRealValue(R_interval)
            else
               error stop 'Keyword interval is mandatory for p-adaptation if adaptation mode is time'
            end if
            this % iter_interval   = huge(this % iter_interval)
            this % nextAdaptationTime = t0   ! + this % time_interval
         case ("iteration")
            this % adaptation_mode = ADAPT_DYNAMIC_ITER
            this % time_interval   = huge(this % time_interval)
            if ( R_interval /= "" ) then
               this % iter_interval   = GetIntValue(R_interval)
            else
               error stop 'Keyword interval is mandatory for p-adaptation if adaptation mode is iteration'
            end if
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
	  
!     Maximum polynomial difference between neighbour elements
!     --------------------------------------------------------
      if ( R_Jump /= "" ) then
	    this % pJump = GetIntValue(R_Jump)
	  else
		this % pJump = 10 
	  end if 

!     Acoustics
!     ---------------
      if ( allocated(R_acoustics) ) then
         this % acoustics = R_acoustics
      end if

      if (this % acoustics) then
         if ( allocated(R_acoustic_tol) ) then
            this % acoustic_tol = R_acoustic_tol
         end if
         if ( allocated(R_acoustic_distance) ) then
            this % acoustic_distance = R_acoustic_distance
         end if
         if ( R_observer /= "" ) then
            this % observer = getRealArrayFromString(R_observer)
         else
            error stop 'Keyword observer is mandatory for p-adaptation with acoustics'
         end if
         if ( R_acoustic_sources /= "" ) then
            call getCharArrayFromString (R_acoustic_sources,BC_STRING_LENGTH,this % acoustic_sources)
            do i=1, size(this % acoustic_sources)
               call toLower(this % acoustic_sources(i))
            end do
         else
            error stop 'Keyword acoustic sources is mandatory for p-adaptation with acoustics'
         end if

         if ( R_acoustic_variable /= "" ) then
            select case ( trim (R_acoustic_variable) )
            case ("p")
#ifdef NAVIERSTOKES
               this % acoustic_variable = 7
#elif defined(MULTIPHASE)
               this % acoustic_variable = 7
#endif
            case ("rho")
               this % acoustic_variable = 8
            case default
               this % acoustic_variable = 7
			   if ( MPI_Process % isRoot ) WRITE(STD_OUT,*) 'Undefined acoustic variable. Using pressure by default.'
            end select
         else
            this % acoustic_variable = 7
			if ( MPI_Process % isRoot ) WRITE(STD_OUT,*) 'Undefined acoustic variable. Using pressure by default.'
         end if

         call mesh % DefineAcousticElements(this % observer, this % acoustic_sources, this % acoustic_distance, surfacesMesh % zones)

      end if

!     Adaptive dt
!     -----------
      if (controlVariables % containsKey("adaptive dt")) then
         adaptive_dt = controlVariables % logicalValueForKey("adaptive dt")
         if (adaptive_dt) then
            if (this % error_estimation) then
               call adaptiveTimeStep % construct(controlVariables, t0) ! Construct the adaptive time step
               adaptiveTimeStep % error_variable = this % error_variable ! Set the error variable for the adaptive time step
            else
               error stop 'Adaptive dt is only available with error estimation'
            end if
         end if
      end if 
      
!
!     Stopwatch events
!     ****************
!  
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
      
      
!     Adaptation overenriching
!     *************************
      
      call getNoOfOverEnrichingBoxes(no_of_overen_boxes)
         
      if (no_of_overen_boxes > 0) then
         allocate ( this % overenriching(no_of_overen_boxes) )
         
         do i = 1, no_of_overen_boxes
            call this % overenriching(i) % initialize (i)
         end do
      end if
	  
!     Adaptation based on variable value
!     **********************************
      
      call getNoOfpAdaptVariables(no_of_overen_boxes)
         
      if (no_of_overen_boxes > 0) then
         allocate ( this % adaptVariable(no_of_overen_boxes) )
         
         do i = 1, no_of_overen_boxes
            call this % adaptVariable(i) % initialize (i)
         end do
      end if

!     Policy definition for the RL agent
!     **********************************
      if ( agentFile /= "" ) then
         call this % agent % construct(trim(agentFile))
         NMIN(1) = max(NMIN(1), this % agent % pmin)
         NMIN(2) = max(NMIN(2), this % agent % pmin)
         NMIN(3) = max(NMIN(3), this % agent % pmin)
         this % NxyzMax(1) = min(this % NxyzMax(1), this % agent % pmax)
         this % NxyzMax(2) = min(this % NxyzMax(2), this % agent % pmax)
         this % NxyzMax(3) = min(this % NxyzMax(3), this % agent % pmax)
      else
         error stop 'Keyword agentFile is mandatory for p-adaptation with RL'
      end if

      safedeallocate(R_increasing)
      safedeallocate(reorganize_z)
      safedeallocate(R_restart)
      safedeallocate(R_ErrorEstimation)
      safedeallocate(R_acoustics)
      safedeallocate(R_tolerance)
      safedeallocate(R_threshold)
      safedeallocate(R_acoustic_tol)
      safedeallocate(R_acoustic_distance)

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
      class(pAdaptationRL_t) :: this
      !--------------------------------------
      
      call this % agent % destruct()
      
      safedeallocate  (this % conformingBoundaries)
      safedeallocate  (this % overenriching)
      safedeallocate  (this % acoustic_sources)
      
   end subroutine pAdaptation_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------------
!  Main routine for adapting the polynomial order in all elements based on 
!  the Value Iteration RL agent
!  ------------------------------------------------------------------------
   subroutine pAdaptation_pAdapt(this, sem, itera, t, computeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables, adaptiveTimeStep)
      use AnisFASMultigridClass
#if defined(FLOW) 
      use SpongeClass, only: sponge, creatRampSponge
#endif
      implicit none
      !-arguments----------------------------
      class(pAdaptationRL_t)     :: this              !<> Adaptation class
      type(DGSem)                :: sem               !<> sem
      integer                    :: itera             !<  iteration
      real(kind=RP)              :: t                 !< time!!
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivativeIsolated
      type(FTValueDictionary)    :: controlVariables  !<> Input variables (that can be modified depending on the user input)
      type(adaptiveTimeStep_t)   :: adaptiveTimeStep  !<> Adaptive time step class
      !-local-variables----------------------
      integer                    :: eID                                     !   Element counter
      integer                    :: Dir                                     !   Direction
      integer                    :: NNew(3,sem % mesh % no_of_elements)     !   New polynomial orders of mesh (after adaptation!)
      integer, save              :: Stage = 0                               !   Stage of p-adaptation for the increasing method
      integer                    :: i                                       !   Counters
      character(len=LINE_LENGTH) :: AdaptedMeshFile
      logical                    :: last
      integer                    :: Ndir = 3, Ndir_acoustics = 4
	   integer                    :: maxPGlob, minPGlob, maxP, minP
      ! integer                    :: pressure_var = 7, rho_var = 8
      !-mpi-variables-------------------------
      integer                    :: ierr
      integer                    :: local_DOFs, global_DOFs
      integer                    :: adaptedElements, allAdaptedElements
      real(kind=RP)              :: adaptationPercentage
      !--------------------------------------
#if (defined(NAVIERSTOKES) || defined(MULTIPHASE))
      
      Stage = Stage + 1    
!
!     -------------------------------------------------------------
!     Find the polynomial order that fulfills the error requirement
!     -------------------------------------------------------------
!
      ! if (this % error_estimation) then
      !    Ndir = 4
      ! end if
      
      call Stopwatch % Start("pAdapt: PolOrder selection")
      adaptedElements = 0
      if (this % acoustics) then
!$omp parallel do schedule(runtime) private(eID)
         do i = 1, size(sem % mesh % elements_acoustics)
            eID = sem % mesh % elements_acoustics(i)
            call pAdaptation_pAdaptRL_SelectElemPolorders (this, sem % mesh % elements(eID) , NNew(:,eID), this % acoustic_tol, Ndir_acoustics, this % acoustic_variable) 
            if ( .not. all( sem % mesh % elements(eID)  % Nxyz == NNew(:,eID)) ) then
!$omp critical
               adaptedElements = adaptedElements + 1
!$omp end critical
            end if
         end do
!$omp end parallel do

!$omp parallel do schedule(runtime) private(eID)
         do i = 1, size(sem % mesh % elements_aerodynamics)
            eID = sem % mesh % elements_aerodynamics(i)
            call pAdaptation_pAdaptRL_SelectElemPolorders (this, sem % mesh % elements(eID) , NNew(:,eID), this % tol, Ndir)
            if ( .not. all( sem % mesh % elements(eID)  % Nxyz == NNew(:,eID)) ) then
!$omp critical
               adaptedElements = adaptedElements + 1
!$omp end critical
            end if
         end do
!$omp end parallel do

      else
!$omp parallel do schedule(runtime)
         do eID = 1, sem % mesh % no_of_elements
            call pAdaptation_pAdaptRL_SelectElemPolorders (this, sem % mesh % elements(eID) , NNew(:,eID), this % tol, Ndir)
            if ( .not. all( sem % mesh % elements(eID)  % Nxyz == NNew(:,eID)) ) then
!$omp critical
               adaptedElements = adaptedElements + 1
!$omp end critical
            end if
         end do
!$omp end parallel do
      end if
      call Stopwatch % Pause("pAdapt: PolOrder selection")

      if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
         call mpi_allreduce ( adaptedElements, allAdaptedElements , 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
         adaptationPercentage = 100.0_RP * allAdaptedElements / sem % mesh % no_of_allElements
#endif
      else
         adaptationPercentage = 100.0_RP * adaptedElements / sem % mesh % no_of_elements
      end if  
      
!
!     -----------------------------------------------------------------------------
!     Adapt only if the percentage of elements to be adapted is above the threshold
!     -----------------------------------------------------------------------------
!
      if (adaptationPercentage > this % threshold .and. this % Adapt) then

         adaptiveTimeStep % pAdapted = .TRUE. ! Set the pAdapted flag to true for the adaptive time step

         if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,*)
               write(STD_OUT,*)
               write(STD_OUT,'(A)') '****     Performing p-Adaptation with a VI Reinforcement Learning agent    ****'
               write(STD_OUT,*)
            end if
#endif
         else
            write(STD_OUT,*)
            write(STD_OUT,*)
            write(STD_OUT,'(A)') '****     Performing p-Adaptation with a VI Reinforcement Learning agent    ****'
            write(STD_OUT,*)
         end if

         if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,*) '****    Percentage of elements to be adapted: ', adaptationPercentage, '% ****'
               write(STD_OUT,*)
            end if
#endif
         else
            write(STD_OUT,*) '****    Percentage of elements to be adapted: ', adaptationPercentage, '% ****'
            write(STD_OUT,*)
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
!     ----------------------------
!     Overenrich specified regions
!     ----------------------------
!
      call OverEnrichRegions(this % overenriching, sem % mesh, NNew, this % NxyzMax, NMIN)
!
!     ----------------------------
!     Adaptation based on variable
!     ----------------------------
!
      call pAdaptVariableRange(this % adaptVariable, sem % mesh, NNew)
!
!     --------------------------------------------------------------------------
!     Restrict polynomial order jump between elements to be pJump (Default is 1)
!     --------------------------------------------------------------------------
!
	  maxP=maxval(NNew)
	  minP=minval(NNew) 
	  maxPGlob = maxP
	  minPGlob = minP
	  
#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
          call mpi_allreduce(maxP, maxPGlob, 1, MPI_INT, MPI_MAX, &
                            MPI_COMM_WORLD, ierr)
          call mpi_allreduce(minP, minPGlob, 1, MPI_INT, MPI_MIN, &
                            MPI_COMM_WORLD, ierr)
      end if
#endif

	  do i=maxPGlob,minPGlob+2,-1
		call pAdapt_CheckNeighbour(sem % mesh, i, this % pJump, NNew)
	  end do 

!
!     ---------------------------------------------------------------
!     Restrict polynomial order jump and make boundaries p-conforming
!     ---------------------------------------------------------------
!
      last = .FALSE.
      do while (.not. last)
         last = .TRUE.
         call this % makeBoundariesPConforming(sem % mesh, NNew, last)
         ! call ReorganizePolOrders(sem % mesh % faces, NNew, last)  #Not implemented for MPI, but not required if pmax<=6 and pmin>=2
      end do

!
!     ----------------------------------
!     Adapt sem to new polynomial orders
!     ----------------------------------
!
      call Stopwatch % Start("pAdapt: Adaptation")
      call sem % mesh % pAdapt_MPI (NNew, controlVariables)
      call Stopwatch % Pause("pAdapt: Adaptation")

!
!     ----------------------------------
!     Reconstruct sponge
!     ----------------------------------
!
#if defined(FLOW) 
      call creatRampSponge(sponge, sem % mesh)
#endif
      
      ! Reconstruct probes
      do i=1, sem % monitors % no_of_probes
         call sem % monitors % probes(i) % Initialization (sem % mesh, i, trim(sem % monitors % solution_file), .FALSE.)
      end do
          
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

!
!     --------------------------------------------------------------------------
!     Perform a reduction to know how many DOFs are in each process
!     --------------------------------------------------------------------------
      local_DOFs = SUM((NNew(1,:)+1)*(NNew(2,:)+1)*(NNew(3,:)+1))

      if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
      call mpi_reduce ( local_DOFs, global_DOFs , 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
      if ( MPI_Process % isRoot ) then
         write(STD_OUT,*) '****    p-Adaptation done, DOFs=', global_DOFs, '****'
         write(STD_OUT,*)
         write(STD_OUT,*)
      end if
#endif
      else
         write(STD_OUT,*) '****    p-Adaptation done, DOFs=', local_DOFs, '****'
         write(STD_OUT,*)
         write(STD_OUT,*)
      end if

   ! End of adapting
   end if

   ! Only adapt once
   this % Adapt = .FALSE.

!End NAVIERSTOKES
#endif 
   end subroutine pAdaptation_pAdapt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------
!  pAdaptation_pAdaptRL_SelectElemPolorders:
!  Select the polynomial orders for one element based on the VI Reinforcement Learning agent
!  -----------------------------------------------------------------------------------------
   subroutine pAdaptation_pAdaptRL_SelectElemPolorders (this, e, NNew, tolerance, Ndir, adaptation_var)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      !-arguments----------------------------------
      class(pAdaptationRL_t), intent(in)    :: this        !<> Adaptation class
      type(Element)         , intent(in)    :: e
      integer               , intent(out)   :: NNew(NDIM)
      real(kind=RP)         , intent(in)    :: tolerance
      integer, optional     , intent(in)    :: adaptation_var
      !-local-variables----------------------------
      integer                    :: Pxyz(NDIM)  ! Initial polynomial order
      integer                    :: i, j, k     ! Coordinate counters
      integer                    :: dir         ! Coordinate direction
      integer                    :: Ndir        ! Number of variables used to adapt
      integer(kind=1)            :: action      ! Action
      integer                    :: indices_dir1(e % Nxyz(1)+1), indices_dir2(e % Nxyz(2)+1), indices_dir3(e % Nxyz(3)+1)
      integer                    :: indices_error1(e % Nxyz(1)+1), indices_error2(e % Nxyz(2)+1), indices_error3(e % Nxyz(3)+1)
      real(kind=RP)              :: Q_dir1(Ndir, 0:e % Nxyz(1)), Q_dir2(Ndir, 0:e % Nxyz(2)), Q_dir3(Ndir, 0:e % Nxyz(3))
      real(kind=RP)              :: Q_error1(0:e % Nxyz(1)), Q_error2(0:e % Nxyz(2)), Q_error3(0:e % Nxyz(3))
      real(kind=RP)              :: minQ, maxQ, min_errorQ, max_errorQ
      real(kind=RP)              :: error_sensor, min_error_sensor
      !--------------------------------------------
      min_error_sensor = 1e-10_RP
      
!     Initialization of P
!     ---------------------------
      Pxyz = e % Nxyz 
      
!     Initialization of NNew
!     -------------------------------
      NNew = -1 ! Initialized to negative value
	  
	  if (maxval(Pxyz).gt.6) then ! The RL only allow pMax .le.6 - This is to bypass the process should higher p is assign by other refinement method
		NNew = Pxyz
	  else
!     --------------------------------
!     Select the polynomial order in Direction 1
!     --------------------------------
      action = -1 ! Initialized to minimum value
      error_sensor = min_error_sensor ! Initialized to minimum value

      do i = 0, Pxyz(3)
         do j = 0, Pxyz(2)
            !Compute the state variables for each Gauss-node in axis 1
            do k = 0, Pxyz(1)
               do dir = 1, Ndir
                  Q_dir1(dir, k) = e % storage % Q(dir+1, k, j, i) !RHOU, RHOV, RHOW, RHOE (or error variable)
                  if (this % error_estimation) then 
                     if(mod(this % error_variable, 3) == mod(dir, 3) .and. dir <= 3) then
                        if (this % error_variable <= 6) then
                           Q_error1(k) = e % storage % Q(dir+1, k, j, i)
                           if (this % error_variable <= 3) then
                              Q_error1(k) = Q_error1(k) / e % storage % Q(1, k, j, i)
                           end if
                        else if (this % error_variable == 7) then
#ifdef NAVIERSTOKES
                           !Pressure
                           Q_error1(k) = thermodynamics % gammaMinus1*(e % storage % Q(5,k,j,i) - 0.5_RP*(e % storage % Q(2,k,j,i)**2 + e % storage % Q(3,k,j,i)**2 + e % storage % Q(4,k,j,i)**2)/e % storage % Q(1,k,j,i))
#elif defined(MULTIPHASE)
                           !Pressure
                           Q_error1(k) = e % storage % QNS(IMP,k,j,i)
#endif
                        else if (this % error_variable == 8) then
                           !Density
                           Q_error1(k) = e % storage % Q(1,k,j,i)
                        end if
                     end if
                     if (dir == Ndir) then 
                        Q_dir1(dir, k) = Q_error1(k) !RHOE is replaced by the error variable (if any)
                     end if
                  end if
                  if (present(adaptation_var) .and. dir==Ndir) then
                     if (adaptation_var == 7) then
#ifdef NAVIERSTOKES
                        !Pressure
                        Q_dir1(dir, k) = thermodynamics % gammaMinus1*(e % storage % Q(5,k,j,i) - 0.5_RP*(e % storage % Q(2,k,j,i)**2 + e % storage % Q(3,k,j,i)**2 + e % storage % Q(4,k,j,i)**2)/e % storage % Q(1,k,j,i))
#elif defined(MULTIPHASE)
                        !Pressure
                        Q_dir1(dir, k) = e % storage % QNS(IMP,k,j,i)
#endif
                     else if (adaptation_var == 8) then
                        !Density
                        Q_dir1(dir, k) = e % storage % Q(1,k,j,i)
                     end if
                  end if

               enddo
            enddo
            !Select the most restrictive action among RHOU, RHOV, RHOW, RHOE (or error variable)
            do dir = 1, Ndir
               ! Find minimum and maximum values in axis 1
               minQ = minval(Q_dir1(dir, :))
               maxQ = maxval(Q_dir1(dir, :))

               !p-adaptation
               if (any(.not. ieee_is_finite(Q_dir1(dir, :))) .or. (maxQ - minQ < tolerance)) then
                  indices_dir1(:) = this % agent % smax + 1
                  ! Choose the best action: +1 increase polynomial order, -1 decrease polynomial order, 0 do nothing
                  if (Pxyz(1) > 2) then
                     action = max(action, this % agent % policy(Pxyz(1) - this % agent % pmin + 1) % matrix % getData(indices_dir1))
                  else if (Pxyz(1) > 1) then
                     action = max(action, -1)
                  else
                     action = max(action, 0)
                  end if
               else
                  !Compute non-dimensional state variables for each Gauss-node in axis 1
                  indices_dir1(:) = nint(max(-real(this % agent % smax, RP), min(real(this % agent % smax, RP), &
                                         2.0_RP * this % agent % smax * (Q_dir1(dir, :) - minQ) / (maxQ - minQ) - this % agent % smax))) + &
                                         this % agent % smax + 1
                  ! Choose the best action: +1 increase polynomial order, -1 decrease polynomial order, 0 do nothing
                  if (Pxyz(1) > 2) then
                     action = max(action, this % agent % policy(Pxyz(1) - this % agent % pmin + 1) % matrix % getData(indices_dir1))
                  else if (Pxyz(1) > 1) then
                     action = max(action, this % agent % policy(Pxyz(1) - this % agent % pmin + 1) % matrix % getData(indices_dir1), 0)
                  else
                     action = 1
                  end if
               end if           
            enddo

            !Error estimator
            if (this % error_estimation) then
               min_errorQ = minval(Q_error1(:))
               max_errorQ = maxval(Q_error1(:))

               if (any(.not. ieee_is_finite(Q_error1(:))) .or. (max_errorQ - min_errorQ < tolerance)) then
                  indices_error1(:) = this % agent % smax + 1
                  ! Compute the error
                  if (Pxyz(1) > 2) then
                     error_sensor = error_sensor + max(this % agent % policy(Pxyz(1) - this % agent % pmin + 1) % matrix % getError(indices_error1), min_error_sensor) * (max_errorQ - min_errorQ)**2.0_RP / 4.0_RP
                  ! else if (Pxyz(1) > 1) then
                  !    error_sensor = error_sensor + 0.0_RP
                  ! else
                  !    error_sensor = error_sensor + 0.0_RP
                  end if
               else
                  !Compute non-dimensional state variables for each Gauss-node in axis 1
                  indices_error1(:) = nint(max(-real(this % agent % smax, RP), min(real(this % agent % smax, RP), &
                                          2.0_RP * this % agent % smax * (Q_error1(:) - min_errorQ) / (max_errorQ - min_errorQ) - this % agent % smax))) + &
                                          this % agent % smax + 1
                  ! Compute the error
                  if (Pxyz(1) > 2) then
                     error_sensor = error_sensor + max(this % agent % policy(Pxyz(1) - this % agent % pmin + 1) % matrix % getError(indices_error1), min_error_sensor) * (max_errorQ - min_errorQ)**2.0_RP / 4.0_RP
                  else if (Pxyz(1) > 1) then
                     error_sensor = error_sensor + max(this % agent % policy(Pxyz(1) - this % agent % pmin + 1) % matrix % getError(indices_error1), min_error_sensor) * (max_errorQ - min_errorQ)**2.0_RP / 4.0_RP
                  else
                     error_sensor = error_sensor + 1.0_RP * (max_errorQ - min_errorQ)**2.0_RP / 4.0_RP
                  end if
               end if
            end if

         enddo
      enddo
      ! Update NNew and make sure it is within the boundaries
      NNew(1) = min(max(Pxyz(1) + action, NMIN(1)), this % NxyzMax(1))
      if (this % error_estimation) then
         e % storage % sensor = sqrt(error_sensor / ((Pxyz(2)+1) * (Pxyz(3)+1))) 
      end if

!     --------------------------------
!     Select the polynomial order in Direction 2
!     --------------------------------
      action = -1 ! Initialized to minimum value
      error_sensor = min_error_sensor

      do i = 0, Pxyz(3)
         do j = 0, Pxyz(1)
            !Compute the state variables for each Gauss-node in axis 2
            do k = 0, Pxyz(2)
               do dir = 1, Ndir
                  Q_dir2(dir, k) = e % storage % Q(dir+1, j, k, i) !RHOU, RHOV, RHOW, RHOE (or error variable)
                  if (this % error_estimation) then 
                     if(mod(this % error_variable, 3) == mod(dir, 3) .and. dir <= 3) then
                        if (this % error_variable <= 6) then
                           Q_error2(k) = e % storage % Q(dir+1, j, k, i)
                           if (this % error_variable <= 3) then
                              Q_error2(k) = Q_error2(k) / e % storage % Q(1, j, k, i)
                           end if
                        else if (this % error_variable == 7) then
#ifdef NAVIERSTOKES
                           !Pressure
                           Q_error2(k) = thermodynamics % gammaMinus1*(e % storage % Q(5,j,k,i) - 0.5_RP*(e % storage % Q(2,j,k,i)**2 + e % storage % Q(3,j,k,i)**2 + e % storage % Q(4,j,k,i)**2)/e % storage % Q(1,j,k,i))
#elif defined(MULTIPHASE)
                           !Pressure
                           Q_error2(k) = e % storage % QNS(IMP,j,k,i)
#endif
                        else if (this % error_variable == 8) then
                           !Density
                           Q_error2(k) = e % storage % Q(1,j,k,i)
                        end if
                     end if
                     if (dir == Ndir) then 
                        Q_dir2(dir, k) = Q_error2(k) !RHOE is replaced by the error variable (if any)
                     end if
                  end if

                  if (present(adaptation_var) .and. dir==Ndir) then
                     if (adaptation_var == 7) then
#ifdef NAVIERSTOKES
                        !Pressure
                        Q_dir2(dir, k) = thermodynamics % gammaMinus1*(e % storage % Q(5,j,k,i) - 0.5_RP*(e % storage % Q(2,j,k,i)**2 + e % storage % Q(3,j,k,i)**2 + e % storage % Q(4,j,k,i)**2)/e % storage % Q(1,j,k,i))
#elif defined(MULTIPHASE)
                        !Pressure
                        Q_dir2(dir, k) = e % storage % QNS(IMP,j,k,i)
#endif
                     else if (adaptation_var == 8) then
                        !Density
                        Q_dir2(dir, k) = e % storage % Q(1,j,k,i)
                     end if
                  end if
               enddo
            enddo
            !Select the most restrictive action among RHOU, RHOV, RHOW, RHOE (or error variable)
            do dir = 1, Ndir
               ! Find minimum and maximum values in axis 2
               minQ = minval(Q_dir2(dir, :))
               maxQ = maxval(Q_dir2(dir, :))

               !p-adaptation
               if (any(.not. ieee_is_finite(Q_dir2(dir, :))) .or. (maxQ - minQ < tolerance)) then
                  indices_dir2(:) = this % agent % smax + 1
                  ! Choose the best action: +1 increase polynomial order, -1 decrease polynomial order, 0 do nothing
                  if (Pxyz(2) > 2) then
                     action = max(action, this % agent % policy(Pxyz(2) - this % agent % pmin + 1) % matrix % getData(indices_dir2))
                  else if (Pxyz(2) > 1) then
                     action = max(action, -1)
                  else
                     action = max(action, 0)
                  end if
               else
                  !Compute non-dimensional state variables for each Gauss-node in axis 2
                  indices_dir2(:) = nint(max(-real(this % agent % smax, RP), min(real(this % agent % smax, RP), &
                                         2.0_RP * this % agent % smax * (Q_dir2(dir, :) - minQ) / (maxQ - minQ) - this % agent % smax))) + &
                                         this % agent % smax + 1
                  ! Choose the best action: +1 increase polynomial order, -1 decrease polynomial order, 0 do nothing
                  if (Pxyz(2) > 2) then
                     action = max(action, this % agent % policy(Pxyz(2) - this % agent % pmin + 1) % matrix % getData(indices_dir2))
                  else if (Pxyz(2) > 1) then
                     action = max(action, this % agent % policy(Pxyz(2) - this % agent % pmin + 1) % matrix % getData(indices_dir2), 0)
                  else
                     action = 1
                  end if
               end if
            enddo

            !Error estimator
            if (this % error_estimation) then
               min_errorQ = minval(Q_error2(:))
               max_errorQ = maxval(Q_error2(:))

               if (any(.not. ieee_is_finite(Q_error2(:))) .or. (max_errorQ - min_errorQ < tolerance)) then
                  indices_error2(:) = this % agent % smax + 1
                  ! Compute the error
                  if (Pxyz(2) > 2) then
                     error_sensor = error_sensor + max(this % agent % policy(Pxyz(2) - this % agent % pmin + 1) % matrix % getError(indices_error2), min_error_sensor) * (max_errorQ - min_errorQ)**2.0_RP / 4.0_RP
                  ! else if (Pxyz(2) > 1) then
                  !    error_sensor = error_sensor + 0.0_RP
                  ! else
                  !    error_sensor = error_sensor + 0.0_RP
                  end if
               else
                  !Compute non-dimensional state variables for each Gauss-node in axis 2
                  indices_error2(:) = nint(max(-real(this % agent % smax, RP), min(real(this % agent % smax, RP), &
                                         2.0_RP * this % agent % smax * (Q_error2(:) - min_errorQ) / (max_errorQ - min_errorQ) - this % agent % smax))) + &
                                         this % agent % smax + 1
                  ! Compute the error
                  if (Pxyz(2) > 2) then
                     error_sensor = error_sensor + max(this % agent % policy(Pxyz(2) - this % agent % pmin + 1) % matrix % getError(indices_error2), min_error_sensor) * (max_errorQ - min_errorQ)**2.0_RP / 4.0_RP
                  else if (Pxyz(2) > 1) then
                     error_sensor = error_sensor + max(this % agent % policy(Pxyz(2) - this % agent % pmin + 1) % matrix % getError(indices_error2), min_error_sensor) * (max_errorQ - min_errorQ)**2.0_RP / 4.0_RP
                  else
                     error_sensor = error_sensor + 1.0_RP * (max_errorQ - min_errorQ)**2.0_RP / 4.0_RP
                  end if
               end if
            end if
         enddo
      enddo
      ! Update NNew and make sure it is within the boundaries
      NNew(2) = min(max(Pxyz(2) + action, NMIN(2)), this % NxyzMax(2))
      if (this % error_estimation) then
         if (this % avg_error_type) then
            e % storage % sensor = e % storage % sensor + sqrt(error_sensor / ((Pxyz(1)+1) * (Pxyz(3)+1)))
         else
            e % storage % sensor = max(e % storage % sensor, sqrt(error_sensor / ((Pxyz(1)+1) * (Pxyz(3)+1))))
         end if
      end if

!     --------------------------------
!     Select the polynomial order in Direction 3
!     --------------------------------
      action = -1 ! Initialized to minimum value
      error_sensor = min_error_sensor

      do i = 0, Pxyz(2)
         do j = 0, Pxyz(1)
            !Compute the state variables for each Gauss-node in axis 3
            do k = 0, Pxyz(3)
               do dir = 1, Ndir
                  Q_dir3(dir, k) = e % storage % Q(dir+1, j, i, k) !RHOU, RHOV, RHOW, RHOE (or error variable)
                  if (this % error_estimation) then 
                     if(mod(this % error_variable, 3) == mod(dir, 3) .and. dir <= 3) then
                        if (this % error_variable <= 6) then
                           Q_error3(k) = e % storage % Q(dir+1, j, i, k)
                           if (this % error_variable <= 3) then
                              Q_error3(k) = Q_error3(k) / e % storage % Q(1, j, i, k)
                           end if
                        else if (this % error_variable == 7) then
#ifdef NAVIERSTOKES
                           !Pressure
                           Q_error3(k) = thermodynamics % gammaMinus1*(e % storage % Q(5,j,i,k) - 0.5_RP*(e % storage % Q(2,j,i,k)**2 + e % storage % Q(3,j,i,k)**2 + e % storage % Q(4,j,i,k)**2)/e % storage % Q(1,j,i,k))
#elif defined(MULTIPHASE)
                           !Pressure
                           Q_error3(k) = e % storage % QNS(IMP,j,i,k)
#endif
                        else if (this % error_variable == 8) then
                           !Density
                           Q_error3(k) = e % storage % Q(1,j,i,k)
                        end if
                     end if
                     if (dir == Ndir) then 
                        Q_dir3(dir, k) = Q_error3(k) !RHOE is replaced by the error variable (if any)
                     end if
                  end if

                  if (present(adaptation_var) .and. dir==Ndir) then
                     if (adaptation_var == 7) then
#ifdef NAVIERSTOKES
                        !Pressure
                        Q_dir3(dir, k) = thermodynamics % gammaMinus1*(e % storage % Q(5,j,i,k) - 0.5_RP*(e % storage % Q(2,j,i,k)**2 + e % storage % Q(3,j,i,k)**2 + e % storage % Q(4,j,i,k)**2)/e % storage % Q(1,j,i,k))
#elif defined(MULTIPHASE)
                        !Pressure
                        Q_dir3(dir, k) = e % storage % QNS(IMP,j,i,k)
#endif
                     else if (adaptation_var == 8) then
                        !Density
                        Q_dir3(dir, k) = e % storage % Q(1,j,i,k)
                     end if
                  end if
               enddo
            enddo
            !Select the most restrictive action among RHOU, RHOV, RHOW, RHOE (or error variable)
            do dir = 1, Ndir
               ! Find minimum and maximum values in axis 3
               minQ = minval(Q_dir3(dir, :))
               maxQ = maxval(Q_dir3(dir, :))

               !p-adaptation
               if (any(.not. ieee_is_finite(Q_dir3(dir, :))) .or. (maxQ - minQ < tolerance)) then
                  indices_dir3(:) = this % agent % smax + 1
                  ! Choose the best action: +1 increase polynomial order, -1 decrease polynomial order, 0 do nothing
                  if (Pxyz(3) > 2) then
                     action = max(action, this % agent % policy(Pxyz(3) - this % agent % pmin + 1) % matrix % getData(indices_dir3))
                  else if (Pxyz(3) > 1) then
                     action = max(action, -1)
                  else
                     action = max(action, 0)
                  end if
               else
                  !Compute non-dimensional state variables for each Gauss-node in axis 3
                  indices_dir3(:) = nint(max(-real(this % agent % smax, RP), min(real(this % agent % smax, RP), &
                                         2.0_RP * this % agent % smax * (Q_dir3(dir, :) - minQ) / (maxQ - minQ) - this % agent % smax))) + &
                                         this % agent % smax + 1
                  ! Choose the best action: +1 increase polynomial order, -1 decrease polynomial order, 0 do nothing
                  if (Pxyz(3) > 2) then
                     action = max(action, this % agent % policy(Pxyz(3) - this % agent % pmin + 1) % matrix % getData(indices_dir3))
                  else if (Pxyz(3) > 1) then
                     action = max(action, this % agent % policy(Pxyz(3) - this % agent % pmin + 1) % matrix % getData(indices_dir3), 0)
                  else
                     action = 1
                  end if
               end if
            enddo

            !Error estimator
            if (this % error_estimation) then
               min_errorQ = minval(Q_error3(:))
               max_errorQ = maxval(Q_error3(:))

               if (any(.not. ieee_is_finite(Q_error3(:))) .or. (max_errorQ - min_errorQ < tolerance)) then
                  indices_error3(:) = this % agent % smax + 1
                  ! Compute the error
                  if (Pxyz(3) > 2) then
                     error_sensor = error_sensor + max(this % agent % policy(Pxyz(3) - this % agent % pmin + 1) % matrix % getError(indices_error3), min_error_sensor) * (max_errorQ - min_errorQ)**2.0_RP / 4.0_RP
                  ! else if (Pxyz(3) > 1) then
                  !    error_sensor = error_sensor + 0.0_RP
                  ! else
                  !    error_sensor = error_sensor + 0.0_RP
                  end if
               else
                  !Compute non-dimensional state variables for each Gauss-node in axis 3
                  indices_error3(:) = nint(max(-real(this % agent % smax, RP), min(real(this % agent % smax, RP), &
                                          2.0_RP * this % agent % smax * (Q_error3(:) - min_errorQ) / (max_errorQ - min_errorQ) - this % agent % smax))) + &
                                          this % agent % smax + 1
                  ! Compute the error
                  if (Pxyz(3) > 2) then
                     error_sensor = error_sensor + max(this % agent % policy(Pxyz(3) - this % agent % pmin + 1) % matrix % getError(indices_error3), min_error_sensor) * (max_errorQ - min_errorQ)**2.0_RP / 4.0_RP
                  else if (Pxyz(3) > 1) then
                     error_sensor = error_sensor + max(this % agent % policy(Pxyz(3) - this % agent % pmin + 1) % matrix % getError(indices_error3), min_error_sensor) * (max_errorQ - min_errorQ)**2.0_RP / 4.0_RP
                  else
                     error_sensor = error_sensor + 1.0_RP * (max_errorQ - min_errorQ)**2.0_RP / 4.0_RP
                  end if
               end if
            end if

         enddo
      enddo
      ! Update NNew and make sure it is within the boundaries
      NNew(3) = min(max(Pxyz(3) + action, NMIN(3)), this % NxyzMax(3)) 
      if (this % error_estimation) then
         if (this % avg_error_type) then  
            e % storage % sensor = e % storage % sensor + sqrt(error_sensor / ((Pxyz(1)+1) * (Pxyz(2)+1)))
            e % storage % sensor = e % storage % sensor / 3.0_RP
         else
            e % storage % sensor = max(e % storage % sensor, sqrt(error_sensor / ((Pxyz(1)+1) * (Pxyz(2)+1))))
         end if
      end if
	  
	  end if 
      
   end subroutine pAdaptation_pAdaptRL_SelectElemPolorders 
   
end module pAdaptationClassRL
