      PROGRAM HORSES3DMainNSTPT

      USE SMConstants
      use FTValueDictionaryClass
      USE PhysicsStorage
      USE SharedBCModule
      USE DGSEMClass
      USE TimeIntegratorClass
      USE mainKeywordsModule
      USE Headers
      use StopwatchClass
      use MPI_Process_Info
      use SpatialDiscretization
      use pAdaptationClass          , only: GetMeshPolynomialOrders
      use NodalStorageClass
      use ManufacturedSolutionsNS
      use FluidData
      use FileReaders               , only: ReadControlFile
      use FileReadingUtilities      , only: getFileName
      use InterpolationMatrices     , only: Initialize_InterpolationMatrices, Finalize_InterpolationMatrices
      use ProblemFileFunctions
      use BoundaryConditions        , only: DestructBoundaryConditions
	  use PartitionedMeshClass
#ifdef _HAS_MPI_
      use mpi
#endif

      IMPLICIT NONE

      procedure(UserDefinedStartup_f)     :: UserDefinedStartup
      procedure(UserDefinedFinalSetup_f)  :: UserDefinedFinalSetup
      procedure(UserDefinedFinalize_f)    :: UserDefinedFinalize
      procedure(UserDefinedTermination_f) :: UserDefinedTermination
      TYPE( FTValueDictionary)            :: controlVariables
      TYPE( DGSem )                       :: sem
      TYPE( TimeIntegrator_t )            :: timeIntegrator
      LOGICAL                             :: success, saveGradients, saveSensor, saveLES, saveSource
	  logical                             :: generateMonitor = .TRUE.
	  logical                             :: optimizePartitionLevel=.FALSE.
	  logical                             :: isMLRK=.FALSE.
      integer                             :: initial_iteration
      INTEGER                             :: ierr
      real(kind=RP)                       :: initial_time
      character(len=LINE_LENGTH)          :: solutionFileName
      integer, allocatable                :: Nx(:), Ny(:), Nz(:)
      integer                             :: Nmax


      call SetSolver(NAVIERSTOKES_SOLVER)
!
!     -----------------------------------------
!     Start measuring the total simulation time
!     -----------------------------------------
!
      call Stopwatch % CreateNewEvent("TotalTime")
      call Stopwatch % Start("TotalTime")
!
!     ---------------
!     Initializations
!     ---------------
!
      call MPI_Process % Init
      call CheckIfTheVersionIsRequested
!
!     ----------------------------------------------------------------------------------
!     The main is always compiled, so that __DATE__ and __TIME__ are updated accordingly
!     ----------------------------------------------------------------------------------
!
      if ( MPI_Process % doMPIAction ) then
         CALL Main_Header("HORSES3D High-Order (DG) Spectral Element Parallel Navier-Stokes Solver",__DATE__,__TIME__)

      else
         CALL Main_Header("HORSES3D High-Order (DG) Spectral Element Sequential Navier-Stokes Solver",__DATE__,__TIME__)

      end if

      CALL controlVariables % initWithSize(16)
      CALL UserDefinedStartup
      CALL ConstructSharedBCModule

      CALL ReadControlFile( controlVariables )
      CALL CheckInputIntegrity(controlVariables, success)
      IF(.NOT. success)   error stop "Control file reading error"

!
!     ----------------
!     Set up the DGSEM
!     ----------------
!
      CALL ConstructPhysicsStorage( controlVariables, success )
      IF(.NOT. success)   error stop "Physics parameters input error"

      ! Initialize manufactured solutions if necessary
      sem % ManufacturedSol = controlVariables % containsKey("manufactured solution")

      IF (sem % ManufacturedSol) THEN
         CALL InitializeManufacturedSol(controlVariables % StringValueForKey("manufactured solution",LINE_LENGTH))
      END IF

      call GetMeshPolynomialOrders(controlVariables,Nx,Ny,Nz,Nmax)
      call InitializeNodalStorage (controlVariables ,Nmax)
      call Initialize_InterpolationMatrices(Nmax)

	  ! MPI Partition option for MLRK time integeration
	  if ((MPI_Process % doMPIAction).and.(trim(controlVariables % stringValueForKey('explicit method', requestedLength = LINE_LENGTH)) == 'multi level rk3')) then
		  
		  isMLRK = .TRUE. ! The time integration is Multi-Level RK
		  
		  if ( controlVariables % ContainsKey("optimized partition") ) then
			 optimizePartitionLevel = controlVariables % LogicalValueForKey ("optimized partition")
		  end if
		  
		  if (optimizePartitionLevel) generateMonitor =.FALSE.  ! Do not generate monitor for the first construction of sem as it will be reconstruct
	  end if
	  
      ! Construct DGSEM library
      call sem % construct (  controlVariables  = controlVariables,                                         &
                                 Nx_ = Nx,     Ny_ = Ny,     Nz_ = Nz,                                                 &
                                 success           = success, generateMonitor =generateMonitor)
!
!     -------------------------
!     Set the initial condition
!     -------------------------
!
      call UserDefinedFinalSetup(sem % mesh, thermodynamics, dimensionless, refValues)		
      call sem % SetInitialCondition(controlVariables, initial_iteration, initial_time)  
!
!     ----------------------------------------------
!     Reconstruct for MLRK explicit time step method
!     ----------------------------------------------
!
      if (isMLRK .and. optimizePartitionLevel .and. MPI_Process % doMPIAction) then
         call sem % reconstruct (  controlVariables  = controlVariables,                                         &
                                 Nx_ = Nx,     Ny_ = Ny,     Nz_ = Nz,                                                 &
                                 success           = success)
		 call UserDefinedFinalSetup(sem % mesh, thermodynamics, dimensionless, refValues)		
		 call sem % SetInitialCondition(controlVariables, initial_iteration, initial_time)
	  end if 
!
!     -----------------------------
!     Initialize the discretization
!     -----------------------------
!
      call Initialize_SpaceAndTimeMethods(controlVariables, sem)

      IF(.NOT. success)   error stop "Mesh reading error"
      IF(.NOT. success)   error stop "Boundary condition specification error"
      																   
      !
      !     -------------------
      !     Build the particles
      !     -------------------
      !
      sem % particles % active = controlVariables % logicalValueForKey("lagrangian particles")
      if ( sem % particles % active ) then
            call sem % particles % construct(sem % mesh, controlVariables, sem % monitors % solution_file)
      endif
!
!     -----------------------------
!     Construct the time integrator
!     -----------------------------
!
      CALL timeIntegrator % construct (controlVariables, sem, initial_iteration, initial_time)
!
!     -----------------
!     Integrate in time
!     -----------------
!
      CALL timeIntegrator % integrate(sem, controlVariables, sem % monitors, sem % samplings, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
!
!     ------------------------------------------
!     Finish measuring the total simulation time
!     ------------------------------------------
!
      call Stopwatch % Pause("TotalTime")
!
!     --------------------------
!     Show simulation statistics
!     --------------------------
!
      call DisplaySimulationStatistics(sem % numberOftimeSteps, sem % mesh)
!
!     -----------------------------------------------------
!     Let the user perform actions on the computed solution
!     -----------------------------------------------------
!
      CALL UserDefinedFinalize(sem % mesh, timeIntegrator % time, sem % numberOfTimeSteps, &
                              sem % maxResidual, thermodynamics, dimensionless, refValues, &
                              sem % monitors, Stopwatch % ElapsedTime("Solver"), &
                              Stopwatch % CPUTime("Solver"))
#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
         call mpi_barrier(MPI_COMM_WORLD, ierr)
      end if
#endif
!
!     -------------------------------------
!     Save the results to the solution file
!     -------------------------------------
!
      IF(controlVariables % stringValueForKey(solutionFileNameKey,LINE_LENGTH) /= "none")     THEN
         solutionFileName = trim(getFileName(controlVariables % stringValueForKey(solutionFileNameKey,LINE_LENGTH))) // ".hsol"
         saveGradients    = controlVariables % logicalValueForKey(saveGradientsToSolutionKey)
         saveSensor       = controlVariables % logicalValueForKey(saveSensorToSolutionKey)
         saveLES          = controlVariables % logicalValueForKey(saveLESToSolutionKey)
         saveSource       = controlVariables % logicalValueForKey(saveSourceToSolutionKey)
         CALL sem % mesh % SaveSolution(sem % numberOfTimeSteps, timeIntegrator % time, solutionFileName, saveGradients, saveSensor, saveLES, saveSource)
         if (controlVariables % logicalValueForKey(saveLambVectorKey)) then
            call sem % mesh % saveLambVector(sem % numberOfTimeSteps, timeIntegrator % time, solutionFileName, controlVariables)
         endif
         if ( sem % particles % active ) then
            call sem % particles % ExportToVTK ( sem % numberOfTimeSteps, sem % monitors % solution_file )
         end if
      END IF
      call Stopwatch % WriteSummaryFile(getFileName(controlVariables % stringValueForKey(solutionFileNameKey,LINE_LENGTH)))

!
!     ---------
!     Finish up
!     ---------
!
      call Stopwatch % destruct
      CALL timeIntegrator % destruct()
      CALL sem % destruct()
      call DestructBoundaryConditions
      call Finalize_SpaceAndTimeMethods
      call Finalize_InterpolationMatrices
      call DestructGlobalNodalStorage()
      CALL destructSharedBCModule

      CALL UserDefinedTermination

      call MPI_Process % Close

      END PROGRAM HORSES3DMainNSTPT
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE CheckInputIntegrity( controlVariables, success )
         use SMConstants
         use Utilities, only: toLower
         USE FTValueDictionaryClass
         USE mainKeywordsModule
         use FTValueClass
         use MPI_Process_Info
         use SpatialDiscretization, only: viscousDiscretizationKey
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(FTValueDictionary) :: controlVariables
         LOGICAL                 :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         CLASS(FTObject), POINTER :: obj
         INTEGER                  :: i
         character(len=LINE_LENGTH)    :: inviscidDiscretization, discretizationNodes

         success = .TRUE.
!
!        Control variables with default value
!        ------------------------------------
         obj => controlVariables % objectForKey(saveGradientsToSolutionKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey(".false.",saveGradientsToSolutionKey)
         end if

         obj => controlVariables % objectForKey(saveSensorToSolutionKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey(".false.",saveSensorToSolutionKey)
         end if

         obj => controlVariables % objectForKey(discretizationNodesKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey("Gauss",discretizationNodesKey)
         end if

         obj => controlVariables % objectForKey(inviscidDiscretizationKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey("Standard",inviscidDiscretizationKey)
         end if

         obj => controlVariables % objectForKey(viscousDiscretizationKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey("BR1",viscousDiscretizationKey)
         end if

         obj => controlVariables % objectForKey(splitFormKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey("Ducros",splitFormKey)
         end if
!
!        Check for inconsistencies in the input variables
!        ------------------------------------------------
         inviscidDiscretization = trim(controlVariables % stringValueForKey(inviscidDiscretizationKey, LINE_LENGTH))
         discretizationNodes = trim(controlVariables % stringValueForKey(discretizationNodesKey, LINE_LENGTH))

         call toLower(inviscidDiscretization)
         call toLower(discretizationNodes)

         if ( (trim(inviscidDiscretization) .eq. "split-form") .and. (trim(discretizationNodes) .eq. "gauss") ) then
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,'(A)') "*** WARNING:    Only Gauss-Lobatto nodes are available for Split-Form discretizations"
               write(STD_OUT,'(A)') "*** WARNING:    Automatically switched to Gauss-Lobatto points"
            end if
            call controlVariables % removeObjectForKey(discretizationNodesKey)
            call controlVariables % addValueForKey("Gauss-Lobatto",discretizationNodesKey)
         end if
!
!        Check the controlVariables created
!        ----------------------------------
         DO i = 1, SIZE(mainKeywords)
            obj => controlVariables % objectForKey(mainKeywords(i))
            IF ( .NOT. ASSOCIATED(obj) )     THEN
               PRINT *, "Input file is missing entry for keyword: ",mainKeywords(i)
               success = .FALSE.
            END IF
         END DO

      END SUBROUTINE checkInputIntegrity

      subroutine DisplaySimulationStatistics(iter,mesh)
         use SMConstants
         use HexMeshClass
         use StopwatchClass
         use Headers
         use MPI_Process_Info
#ifdef _HAS_MPI_
         use mpi
#endif
         implicit none
         integer,    intent(in)      :: iter
         type(HexMesh),   intent(in) :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: eID
         integer                    :: NDOF, localNDOF, ierr
         real(kind=RP)              :: Naverage, localNaverage
         real(kind=RP)              :: t_elaps, t_cpu

         if ( MPI_Process % isRoot ) write(STD_OUT,'(/)')
         call Section_Header("Simulation statistics")
         if ( MPI_Process % isRoot ) write(STD_OUT,'(/)')
!
!        Get mesh-related quantities
!        ---------------------------
         NDOF = 0
         Naverage = 0

         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            NDOF = NDOF + (e % Nxyz(1) + 1)*(e % Nxyz(2) + 1)*(e % Nxyz(3) + 1)
            Naverage = Naverage + e % Nxyz(1) + e % Nxyz(2) + e % Nxyz(3)
            end associate
         end do

         Naverage = Naverage / (3.0_RP * mesh % no_of_elements)
!
!        Perform a broadcast for the MPI solver
!        --------------------------------------
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
            localNDOF = NDOF
            localNaverage = Naverage * 3.0_RP * mesh % no_of_elements
            call mpi_allreduce(localNDOF, NDOF, 1, MPI_INT, MPI_SUM, &
                               MPI_COMM_WORLD, ierr)

            call mpi_allreduce(localNaverage, Naverage, 1, MPI_DOUBLE, MPI_SUM, &
                               MPI_COMM_WORLD, ierr)

            Naverage = Naverage / (3.0_RP * mesh % no_of_allElements)

         end if
#endif

         if ( .not. MPI_Process % isRoot ) return
!
!        Show preprocessing time
!        -----------------------
         t_elaps = Stopwatch % Elapsedtime("Preprocessing")
         t_cpu   = Stopwatch % CPUTime("Preprocessing")

         call Subsection_Header("Preprocessing")

         write(STD_OUT,'(30X,A,I0,A,F5.2,A,I0,A)')      "->   ", mesh % no_of_allElements, &
                                                      " elements with polynomial order ",Naverage," (NDOF = ",NDOF,")."
         write(STD_OUT,'(30X,A,A30,ES10.3,A,ES10.3,A)') "->", "Preprocessing time: ",t_elaps," seconds (total CPU time: ",t_cpu,")."

!
!        Show simulation time
!        --------------------
         write(STD_OUT,'(/)')
         call Subsection_Header("Solver")
         if ( iter .le. 0 ) return

         t_elaps = Stopwatch % ElapsedTime("Solver")
         t_cpu   = Stopwatch % CPUTime("Solver")

         write(STD_OUT,'(30X,A,A30,ES10.3,A)') "->", "Simulation elapsed time: ",t_elaps," seconds."
         write(STD_OUT,'(30X,A,A30,ES10.3,A,ES10.3,A)') "->", "Simulation CPU time: ",t_cpu," seconds (ratio is ",t_cpu/t_elaps ,")."
         write(STD_OUT,'(30X,A,A30,ES10.3,A)') "->", "Solver efficiency: " , t_elaps/(NDOF * iter)*1.0e6_RP, " seconds/(1 Million DOF·iter)."

      end subroutine DisplaySimulationStatistics

      subroutine CheckIfTheVersionIsRequested
         use SMConstants
         use MPI_Process_Info
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: nArgs, i
         character(len=128)    :: arg

         if ( .not. MPI_Process % isRoot ) return

         nArgs = command_argument_count()

         do i = 1, nArgs
            call get_command_argument(i, arg)

            if ( trim(arg) .eq. "--version" ) then
               print*, "Current HORSES version: ", trim(VERSION)
               error stop
            end if
         end do

      end subroutine CheckIfTheVersionIsRequested
