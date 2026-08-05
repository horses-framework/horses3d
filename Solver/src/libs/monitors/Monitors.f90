#include "Includes.h"
module MonitorsClass
   use SMConstants
   use NodalStorageClass
   use HexMeshClass
   use MonitorDefinitions
   use ResidualsMonitorClass
   use VolumeMonitorClass
   use LoadBalancingMonitorClass
   use FileReadingUtilities      , only: getFileName
#ifdef FLOW
   use ProbeClass
#endif
#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE) || defined(ACOUSTIC)
   use StatisticsMonitor
   use SurfaceMonitorClass
#endif
   implicit none
!

   private
   public      Monitor_t
!
!  *****************************
!  Main monitor class definition
!  *****************************
!  
   type Monitor_t
      character(len=LINE_LENGTH)                 :: solution_file
      integer                                    :: no_of_probes
      integer                                    :: no_of_surfaceMonitors
      integer                                    :: no_of_volumeMonitors
      integer                                    :: no_of_loadBalancingMonitors
      integer                                    :: bufferLine
      integer                      , allocatable :: iter(:)
      integer                                    :: dt_restriction
      logical                                    :: write_dt_restriction
      real(kind=RP)                , allocatable :: t(:)
      real(kind=RP)                , allocatable :: SolverSimuTime(:)
      real(kind=RP)                , allocatable :: TotalSimuTime(:)
      type(Residuals_t)                          :: residuals
      class(VolumeMonitor_t)       , allocatable :: volumeMonitors(:)
      class(LoadBalancingMonitor_t), allocatable :: loadBalancingMonitors(:)
#ifdef FLOW
      class(Probe_t)               , allocatable :: probes(:)
#endif
#if defined(NAVIERSTOKES) || defined(INCNS)
      class(SurfaceMonitor_t)      , allocatable :: surfaceMonitors(:)
#endif
#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE) || defined(ACOUSTIC)
      type(StatisticsMonitor_t)                  :: stats
#endif
      contains
         procedure   :: Construct       => Monitors_Construct
         procedure   :: WriteLabel      => Monitor_WriteLabel
         procedure   :: WriteUnderlines => Monitor_WriteUnderlines
         procedure   :: WriteValues     => Monitor_WriteValues
         procedure   :: UpdateValues    => Monitor_UpdateValues
         procedure   :: WriteToFile     => Monitor_WriteToFile
         procedure   :: destruct        => Monitor_Destruct
         procedure   :: copy            => Monitor_Assign
         generic     :: assignment(=)   => copy
   end type Monitor_t
!
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Monitors_Construct( Monitors, mesh, controlVariables ) 
         use FTValueDictionaryClass
         use mainKeywordsModule
         implicit none
         class(Monitor_t)                     :: Monitors
         class(HexMesh), intent(in)           :: mesh
         class(FTValueDictionary), intent(in) :: controlVariables
         
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                         :: fID , io
         integer                         :: i
         character(len=STR_LEN_MONITORS) :: line
         character(len=STR_LEN_MONITORS) :: solution_file                                            
         logical, save                   :: FirstCall = .TRUE.
         logical                         :: saveGradients, saveLambVector
         logical                         :: saveSoundVelocitySquared
         logical                         :: saveGradSoundVelocitySquared
         logical                         :: savePressureSquared
		   logical                         :: monitorMem =.FALSE.

         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: saveLambVectorToStatsKey = "save lamb vector stats"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: saveSoundVelocitySquaredToStatsKey = "save sound velocity squared stats"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: saveGradSoundVelocitySquaredToStatsKey = "save gradient sound velocity squared stats"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: savePressureSquaredToStatsKey = "save pressure squared stats"
!
!        Setup the buffer
!        ----------------
         if (controlVariables % containsKey("monitors flush interval") ) then
            BUFFER_SIZE = controlVariables % integerValueForKey("monitors flush interval")
         end if
         
         allocate ( Monitors % TotalSimuTime(BUFFER_SIZE), &
                    Monitors % SolverSimuTime(BUFFER_SIZE), &
                    Monitors % t(BUFFER_SIZE), &
                    Monitors % iter(BUFFER_SIZE) )
!
!        Memory Monitor
!        ---------------------
         if ( controlVariables % ContainsKey("monitor memory") ) then
			monitorMem = controlVariables % LogicalValueForKey ("monitor memory")
		 end if 
!
!        Get the solution file name
!        --------------------------
         solution_file = controlVariables % stringValueForKey( solutionFileNameKey, requestedLength = STR_LEN_MONITORS )
!
!        Remove the *.hsol termination
!        -----------------------------
         solution_file = trim(getFileName(solution_file))
         Monitors % solution_file = trim(solution_file)
!
!        Search in case file for probes, surface monitors, and volume monitors
!        ---------------------------------------------------------------------
         if (mesh % child) then ! Return doing nothing if this is a child mesh
            Monitors % no_of_probes = 0
            Monitors % no_of_surfaceMonitors = 0
            Monitors % no_of_volumeMonitors = 0
            Monitors % no_of_loadBalancingMonitors = 0
         else
            call getNoOfMonitors( Monitors % no_of_probes, Monitors % no_of_surfaceMonitors, Monitors % no_of_volumeMonitors, Monitors % no_of_loadBalancingMonitors )
         end if
!
!        Initialize
!        ----------
         call Monitors % residuals % Initialization( solution_file , FirstCall, monitorMem)

         allocate ( Monitors % volumeMonitors ( Monitors % no_of_volumeMonitors )  )
         do i = 1 , Monitors % no_of_volumeMonitors
            call Monitors % volumeMonitors(i) % Initialization ( mesh , i, solution_file , FirstCall  )
         end do

         allocate ( Monitors % loadBalancingMonitors ( Monitors % no_of_loadBalancingMonitors )  )
         do i = 1 , Monitors % no_of_loadBalancingMonitors
            call Monitors % loadBalancingMonitors(i) % Initialization ( mesh , i, solution_file , FirstCall )
         end do

#ifdef FLOW
         allocate ( Monitors % probes ( Monitors % no_of_probes )  )
         do i = 1 , Monitors % no_of_probes
            call Monitors % probes(i) % Initialization ( mesh , i, solution_file , FirstCall )
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE) || defined(ACOUSTIC)
         saveGradients    = controlVariables % logicalValueForKey(saveGradientsToSolutionKey)
         saveLambVector   = controlVariables % logicalValueForKey(saveLambVectorToStatsKey)
         saveSoundVelocitySquared = controlVariables % logicalValueForKey(saveSoundVelocitySquaredToStatsKey)
         saveGradSoundVelocitySquared = controlVariables % logicalValueForKey(saveGradSoundVelocitySquaredToStatsKey)
         savePressureSquared = controlVariables % logicalValueForKey(savePressureSquaredToStatsKey)
         call Monitors % stats     % Construct(mesh, saveGradients, saveLambVector, saveSoundVelocitySquared, saveGradSoundVelocitySquared, savePressureSquared)
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
         allocate ( Monitors % surfaceMonitors ( Monitors % no_of_surfaceMonitors )  )
         do i = 1 , Monitors % no_of_surfaceMonitors
            call Monitors % surfaceMonitors(i) % Initialization ( mesh , i, solution_file , FirstCall )
         end do
#endif

         Monitors % write_dt_restriction = controlVariables % logicalValueForKey( "write dt restriction" )
         
         Monitors % bufferLine = 0
         
         FirstCall = .FALSE.
      end subroutine Monitors_Construct

      subroutine Monitor_WriteLabel ( self )
!
!        ***************************************************
!           This subroutine prints the labels for the time
!         integrator Display procedure.
!        ***************************************************
!
         use MPI_Process_Info
         implicit none
         class(Monitor_t)              :: self
         integer                       :: i 
      
         if ( .not. MPI_Process % isRoot ) return
!
!        Write "Iteration" and "Time"
!        ----------------------------
         write ( STD_OUT , ' ( A10    ) ' , advance = "no" ) "Iteration"
         write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) "Time"
!
!        Write residuals labels
!        ----------------------
         call self % residuals % WriteLabel
!
!        Write volume monitors labels
!        -----------------------------
         do i = 1 , self % no_of_volumeMonitors
            call self % volumeMonitors(i) % WriteLabel
         end do
!
!        Write load balancing monitor labels
!        ------------------------------------
         do i = 1 , self % no_of_loadBalancingMonitors
            call self % loadBalancingMonitors(i) % WriteLabel
         end do

#ifdef FLOW
!
!        Write probes labels
!        -------------------
         do i = 1 , self % no_of_probes
            call self % probes(i) % WriteLabel
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
!
!        Write surface monitors labels
!        -----------------------------
         do i = 1 , self % no_of_surfaceMonitors
            call self % surfaceMonitors(i) % WriteLabel
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE)
         call self % stats % WriteLabel
#endif
!
!        Write label for dt restriction
!        ------------------------------
         if (self % write_dt_restriction) write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) "dt restr."

         write(STD_OUT , *) 

      end subroutine Monitor_WriteLabel

      subroutine Monitor_WriteUnderlines( self ) 
!
!        ********************************************************
!              This subroutine displays the underlines for the
!           time integrator Display procedure.
!        ********************************************************
!
         use PhysicsStorage
         use MPI_Process_Info
         implicit none
         class(Monitor_t)                         :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                  :: i, j
         character(len=MONITOR_LENGTH), parameter :: dashes = "----------"

         if ( .not. MPI_Process % isRoot ) return
!
!        Print dashes for "Iteration" and "Time"
!        ---------------------------------------
         write ( STD_OUT , ' ( A10    ) ' , advance = "no" ) trim ( dashes ) 
         write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) trim ( dashes ) 
!
!        Print dashes for residuals
!        --------------------------
         do i = 1 , NCONS
            write(STD_OUT , '(3X,A10)' , advance = "no" ) trim(dashes)
         end do
!
!        Print dashes for memory
!        --------------------------
         if (self % residuals % memory) then
            write(STD_OUT , '(3X,A10)' , advance = "no" ) trim(dashes)
         end if
!
!        Print dashes for volume monitors
!        --------------------------------
         do i = 1 , self % no_of_volumeMonitors  ; do j=1, size ( self % volumeMonitors(i) % values, 1 )
            write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % volumeMonitors(i) % monitorName ) + 2 ) )
         end do                                  ; end do
!
!        Print dashes for load balancing monitor
!        --------------------------------------
         do i = 1 , self % no_of_loadBalancingMonitors ; do j=1, size ( self % loadBalancingMonitors(i) % values, 1 )
            write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % loadBalancingMonitors(i) % monitorName ) + 2 ) )
         end do                                  ; end do

#ifdef FLOW
!
!        Print dashes for probes
!        -----------------------
         do i = 1 , self % no_of_probes
            if ( self % probes(i) % active ) then
               write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % probes(i) % monitorName ) + 2 ) )
            end if
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
!
!        Print dashes for surface monitors
!        ---------------------------------
         do i = 1 , self % no_of_surfaceMonitors
            write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % surfaceMonitors(i) % monitorName ) + 2 ) )
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE)
         if ( self % stats % state .ne. 0 ) write(STD_OUT,'(3X,A10)',advance="no") trim(dashes)
#endif

         
!
!        Print dashes for dt restriction
!        -------------------------------
         if (self % write_dt_restriction) write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) trim ( dashes ) 
         
         write(STD_OUT , *) 

      end subroutine Monitor_WriteUnderlines

      subroutine Monitor_WriteValues ( self )
!
!        *******************************************************
!              This subroutine prints the values for the time
!           integrator Display procedure.
!        *******************************************************
!
         use MPI_Process_Info
         implicit none
         class(Monitor_t)           :: self
         integer                    :: i

         if ( .not. MPI_Process % isRoot ) return
!
!        Print iteration and time
!        ------------------------
         write ( STD_OUT , ' ( I10            ) ' , advance = "no" ) self % iter    ( self % bufferLine ) 
         write ( STD_OUT , ' ( 1X,A,1X,ES10.3 ) ' , advance = "no" ) "|" , self % t ( self % bufferLine ) 
!
!        Print residuals
!        ---------------
         call self % residuals % WriteValues( self % bufferLine )
!
!        Print volume monitors
!        ---------------------
         do i = 1 , self % no_of_volumeMonitors
            call self % volumeMonitors(i) % WriteValues ( self % bufferLine )
         end do
!
!        Print load balancing monitors
!        -----------------------------
         do i = 1 , self % no_of_loadBalancingMonitors
            call self % loadBalancingMonitors(i) % WriteValues ( self % bufferLine )
         end do

#ifdef FLOW
!
!        Print probes
!        ------------
         do i = 1 , self % no_of_probes
            call self % probes(i) % WriteValues ( self % bufferLine )
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
!
!        Print surface monitors
!        ----------------------
         do i = 1 , self % no_of_surfaceMonitors
            call self % surfaceMonitors(i) % WriteValues ( self % bufferLine )
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE)
         call self % stats % WriteValue
#endif
!
!        Print dt restriction
!        --------------------
         if (self % write_dt_restriction) then
            select case (self % dt_restriction)
               case (DT_FIXED) ; write ( STD_OUT , ' ( 1X,A,1X,A10) ' , advance = "no" ) "|" , 'Fixed'
               case (DT_DIFF)  ; write ( STD_OUT , ' ( 1X,A,1X,A10) ' , advance = "no" ) "|" , 'Diffusive'
               case (DT_CONV)  ; write ( STD_OUT , ' ( 1X,A,1X,A10) ' , advance = "no" ) "|" , 'Convective'
            end select
         end if

         write(STD_OUT , *) 

      end subroutine Monitor_WriteValues

      subroutine Monitor_UpdateValues ( self, mesh, t , iter, maxResiduals, Autosave, dt )
!
!        ***************************************************************
!              This subroutine updates the values for the residuals,
!           for the probes, surface and volume monitors.
!        ***************************************************************
!        
         use PhysicsStorage
         use StopwatchClass
         implicit none
         class(Monitor_t)    :: self
         class(HexMesh)      :: mesh
         real(kind=RP)       :: t
         integer             :: iter
         real(kind=RP)       :: maxResiduals(NCONS), dt
         logical             :: Autosave
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i 
!
!        Move to next buffer line
!        ------------------------
         self % bufferLine = self % bufferLine + 1
!
!        Save time, iteration and CPU-time
!        -----------------------
         self % t       ( self % bufferLine )  = t
         self % iter    ( self % bufferLine )  = iter
         self % SolverSimuTime ( self % bufferLine )  = Stopwatch % ElapsedTime("Solver")
         self % TotalSimuTime ( self % bufferLine )  = Stopwatch % ElapsedTime("TotalTime")
!
!        Compute current residuals
!        -------------------------
         call self % residuals % Update( mesh, maxResiduals, self % bufferLine )
!
!        Update volume monitors
!        ----------------------
         do i = 1 , self % no_of_volumeMonitors
            call self % volumeMonitors(i) % Update( mesh , self % bufferLine, t )
         end do
!
!        Update load balancing monitors
!        ------------------------------
         do i = 1 , self % no_of_loadBalancingMonitors
            call self % loadBalancingMonitors(i) % Update( mesh , self % bufferLine )
         end do

#ifdef FLOW
!
!        Update probes
!        -------------
         do i = 1 , self % no_of_probes
            call self % probes(i) % Update( mesh , self % bufferLine )
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
!
!        Update surface monitors
!        -----------------------
         do i = 1 , self % no_of_surfaceMonitors
            call self % surfaceMonitors(i) % Update( mesh , self % bufferLine, iter, autosave, dt )
         end do
#endif
#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE)
!
!        Update statistics
!        -----------------
         call self % stats % Update(mesh, iter, t, trim(self % solution_file) )
#endif

!
!        Update dt restriction
!        ---------------------
         if (self % write_dt_restriction) self % dt_restriction = mesh % dt_restriction 
         
      end subroutine Monitor_UpdateValues

      subroutine Monitor_WriteToFile ( self , mesh, force) 
!
!        ******************************************************************
!              This routine has a double behaviour:
!           force = .true.  -> Writes to file and resets buffers
!           force = .false. -> Just writes to file if the buffer is full
!        ******************************************************************
!
         use MPI_Process_Info
         implicit none
         class(Monitor_t)        :: self
         class(HexMesh)          :: mesh
         logical, optional       :: force
!        ------------------------------------------------
         integer                 :: i 
         logical                 :: forceVal

         if ( present ( force ) ) then
            forceVal = force

         else
            forceVal = .false.

         end if

         if ( forceVal ) then 
!
!           In this case the monitors are exported to their files and the buffer is reset
!           -----------------------------------------------------------------------------
            call self % residuals % WriteToFile ( self % iter , self % t, self % TotalSimuTime, self % SolverSimuTime , self % bufferLine )
   
            do i = 1 , self % no_of_volumeMonitors
               call self % volumeMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do

            do i = 1 , self % no_of_loadBalancingMonitors
               call self % loadBalancingMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do

#ifdef FLOW
            do i = 1 , self % no_of_probes
               call self % probes(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do
#endif
   
#if defined(NAVIERSTOKES) || defined(INCNS)
            do i = 1 , self % no_of_surfaceMonitors
               call self % surfaceMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do
#endif
!
#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE) || defined(ACOUSTIC)
!
!              Write statistics
!              ----------------
            if ( self % bufferLine .eq. 0 ) then
               i = 1
            else
               i = self % bufferLine
            end if
            call self % stats % WriteFile(mesh, self % iter(i), self % t(i), self % solution_file)
#endif
!
!           Reset buffer
!           ------------
            self % bufferLine = 0

         else
!
!           The monitors are exported just if the buffer is full
!           ----------------------------------------------------
            if ( self % bufferLine .eq. BUFFER_SIZE ) then

               call self % residuals % WriteToFile ( self % iter , self % t, self % TotalSimuTime, self % SolverSimuTime, BUFFER_SIZE )

               do i = 1 , self % no_of_volumeMonitors
                  call self % volumeMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
               end do

               do i = 1 , self % no_of_loadBalancingMonitors
                  call self % loadBalancingMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
               end do

#ifdef FLOW
               do i = 1 , self % no_of_probes
                  call self % probes(i) % WriteToFile ( self % iter , self % t , self % bufferLine ) 
               end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
               do i = 1 , self % no_of_surfaceMonitors
                  call self % surfaceMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
               end do
#endif
!
!              Reset buffer
!              ------------
               self % bufferLine = 0

            end if
         end if

      end subroutine Monitor_WriteToFile
      
      subroutine Monitor_Destruct (self)
         implicit none
         class(Monitor_t)        :: self
         
         deallocate (self % iter)
         deallocate (self % t)
         deallocate (self % TotalSimuTime)
         deallocate (self % SolverSimuTime)
         
         call self % residuals % destruct
         
         call self % volumeMonitors % destruct
         safedeallocate(self % volumeMonitors)

         call self % loadBalancingMonitors % destruct
         safedeallocate(self % loadBalancingMonitors)
         
#ifdef FLOW
         call self % probes % destruct
         safedeallocate (self % probes)
#endif
         
#if defined(NAVIERSTOKES) || defined(INCNS)
         call self % surfaceMonitors % destruct
         safedeallocate (self % surfaceMonitors)
#endif         
         
#if defined(NAVIERSTOKES) || defined(INCNS)
         !call self % stats % destruct
#endif         
      end subroutine
      
      impure elemental subroutine Monitor_Assign ( to, from )
         implicit none
         !-arguments--------------------------------------
         class(Monitor_t), intent(inout)  :: to
         type(Monitor_t) , intent(in)     :: from
         !-local-variables--------------------------------
         !------------------------------------------------
         
         to % solution_file               = from % solution_file
         to % no_of_probes                = from % no_of_probes
         to % no_of_surfaceMonitors       = from % no_of_surfaceMonitors
         to % no_of_volumeMonitors        = from % no_of_volumeMonitors
         to % no_of_loadBalancingMonitors = from % no_of_loadBalancingMonitors
         to % bufferLine                  = from % bufferLine
         
         safedeallocate ( to % iter )
         allocate ( to % iter ( size(from % iter) ) )
         to % iter = from % iter
         
         to % dt_restriction        = from % dt_restriction
         to % write_dt_restriction  = from % write_dt_restriction
         
         safedeallocate (to % t)
         allocate (to % t (size (from % t) ) ) 
         to % t = from % t
         
         safedeallocate ( to % TotalSimuTime )
         allocate ( to % TotalSimuTime ( size(from % TotalSimuTime) ) )
         to % TotalSimuTime = from % TotalSimuTime
         
         safedeallocate ( to % SolverSimuTime )
         allocate ( to % SolverSimuTime ( size(from % SolverSimuTime) ) )
         to % SolverSimuTime = from % SolverSimuTime
         
         to % residuals = from % residuals
         
         safedeallocate ( to % volumeMonitors )
         allocate ( to % volumeMonitors ( size(from % volumeMonitors) ) )
         to % volumeMonitors = from % volumeMonitors

         safedeallocate ( to % loadBalancingMonitors )
         allocate ( to % loadBalancingMonitors ( size(from % loadBalancingMonitors) ) )
         to % loadBalancingMonitors = from % loadBalancingMonitors
      
#ifdef FLOW
         safedeallocate ( to % probes )
         allocate ( to % probes ( size(from % probes) ) )
         to % probes = from % probes
#endif
         
#if defined(NAVIERSTOKES) || defined(INCNS)
         safedeallocate ( to % surfaceMonitors )
         allocate ( to % surfaceMonitors ( size(from % surfaceMonitors) ) )
         to % surfaceMonitors = from % surfaceMonitors
#endif
         
#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE)
         to % stats = from % stats
#endif
         
      end subroutine Monitor_Assign
      
!
!//////////////////////////////////////////////////////////////////////////////
!
!        Auxiliars
!
!//////////////////////////////////////////////////////////////////////////////
!
   subroutine getNoOfMonitors(no_of_probes, no_of_surfaceMonitors, no_of_volumeMonitors, no_of_loadBalancingMonitors)
      use ParamfileRegions
      implicit none
      integer, intent(out)    :: no_of_probes
      integer, intent(out)    :: no_of_surfaceMonitors
      integer, intent(out)    :: no_of_volumeMonitors
      integer, intent(out)    :: no_of_loadBalancingMonitors
!
!     ---------------
!     Local variables
!     ---------------
!
      character(len=LINE_LENGTH) :: case_name, line
      integer                    :: fID
      integer                    :: io
!
!     Initialize
!     ----------
      no_of_probes = 0
      no_of_surfaceMonitors = 0
      no_of_volumeMonitors = 0
      no_of_loadBalancingMonitors = 0
!
!     Get case file name
!     ------------------
      call get_command_argument(1, case_name)

!
!     Open case file
!     --------------
      open ( newunit = fID , file = case_name , status = "old" , action = "read" )

!
!     Read the whole file to find monitors
!     ------------------------------------
readloop:do 
         read ( fID , '(A)' , iostat = io ) line

         if ( io .lt. 0 ) then
!
!           End of file
!           -----------
            line = ""
            exit readloop

         elseif ( io .gt. 0 ) then
!
!           Error
!           -----
            errorMessage(STD_OUT)
            error stop "Stopped."

         else
!
!           Succeeded
!           ---------
            line = getSquashedLine( line )

            if ( index ( line , '#defineprobe' ) .gt. 0 ) then
               no_of_probes = no_of_probes + 1

            elseif ( index ( line , '#definesurfacemonitor' ) .gt. 0 ) then
               no_of_surfaceMonitors = no_of_surfaceMonitors + 1 

            elseif ( index ( line , '#definevolumemonitor' ) .gt. 0 ) then
               no_of_volumeMonitors = no_of_volumeMonitors + 1 

            elseif ( index ( line , '#defineloadbalancingmonitor' ) .gt. 0 ) then
               no_of_loadBalancingMonitors = no_of_loadBalancingMonitors + 1

            end if
            
         end if

      end do readloop
!
!     Close case file
!     ---------------
      close(fID)                             

end subroutine getNoOfMonitors

end module MonitorsClass
!
!///////////////////////////////////////////////////////////////////////////////////
!
