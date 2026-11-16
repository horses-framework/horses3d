!
!////////////////////////////////////////////////////////////////////////
!
!      The Problem File contains user defined procedures
!      that are used to "personalize" i.e. define a specific
!      problem to be solved. These procedures include initial conditions,
!      exact solutions (e.g. for tests), etc. and allow modifications 
!      without having to modify the main code.
!
!      The procedures, *even if empty* that must be defined are
!
!      UserDefinedSetUp
!      UserDefinedInitialCondition(mesh)
!      UserDefinedPeriodicOperation(mesh)
!      UserDefinedFinalize(mesh)
!      UserDefinedTermination
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedStartup
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            IMPLICIT NONE  
         END SUBROUTINE UserDefinedStartup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalSetup(mesh &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
                                        )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            class(HexMesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedInitialCondition(mesh &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
                                        )
!
!           ------------------------------------------------
!           called to set the initial condition for the flow
!              - by default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            use smconstants
            use physicsstorage
            use hexmeshclass
            use fluiddata
            implicit none
            class(hexmesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            REAL(KIND=RP) :: x(3)        
            INTEGER       :: i, j, k, eID
            REAL(KIND=RP) :: rho , u , v , w , p
            REAL(KIND=RP) :: L, u_0, rho_0, p_0
            integer       :: Nx, Ny, Nz

#if defined(NAVIERSTOKES)
            associate( gamma => thermodynamics_ % gamma )
            do eID = 1, size(mesh % elements)
               associate( e => mesh % elements(eID) % storage )
               e % Q(IRHO,:,:,:)          = 1.0_RP
               e % Q(IRHOU:IRHOW,:,:,:)   = 0.0_RP
               e % Q(IRHOE,:,:,:)         = 100.0_RP / (gamma - 1.0_RP)
               end associate
            end do
            end associate
#endif

!
!           ---------------------------------------
!           Cahn-Hilliard default initial condition
!           ---------------------------------------
!
#if defined(CAHNHILLIARD)
            call random_seed()
         
            do eid = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eid) % Nxyz(1), &
                          Ny => mesh % elements(eid) % Nxyz(2), &
                          Nz => mesh % elements(eid) % Nxyz(3) )
               associate(e => mesh % elements(eID) % storage)
               call random_number(e % c) 
               e % c = 2.0_RP * (e % c - 0.5_RP)
               end associate
               end associate
            end do
#endif

         end subroutine UserDefinedInitialCondition
#if defined(NAVIERSTOKES)
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
!
!           -------------------------------------------------
!           Used to define an user defined boundary condition
!           -------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
         end subroutine UserDefinedState1

         subroutine UserDefinedGradVars1(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            use VariableConversion, only: GetGradientValues_f
            implicit none
            real(kind=RP), intent(in)          :: x(NDIM)
            real(kind=RP), intent(in)          :: t
            real(kind=RP), intent(in)          :: nHat(NDIM)
            real(kind=RP), intent(in)          :: Q(NCONS)
            real(kind=RP), intent(inout)       :: U(NGRAD)
            procedure(GetGradientValues_f)     :: GetGradients
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedGradVars1

         subroutine UserDefinedNeumann1(x, t, nHat, U_x, U_y, U_z)
!
!           --------------------------------------------------------
!           Used to define a Neumann user defined boundary condition
!           --------------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: U_x(NGRAD)
            real(kind=RP), intent(inout)  :: U_y(NGRAD)
            real(kind=RP), intent(inout)  :: U_z(NGRAD)
         end subroutine UserDefinedNeumann1
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
#if defined(NAVIERSTOKES)
            use MonitorsClass
#endif
            IMPLICIT NONE
            class(HexMesh)               :: mesh
            real(kind=RP)                :: time
            real(kind=RP)                :: dt
#if defined(NAVIERSTOKES)
            type(Monitor_t), intent(in) :: monitors
#else
            logical, intent(in) :: monitors
#endif
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
#if defined(NAVIERSTOKES)
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: Q(NCONS)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(inout) :: S(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
            S = 0.0_RP
   
         end subroutine UserDefinedSourceTermNS
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual &
#if defined(NAVIERSTOKES)
                                                    , thermodynamics_ &
                                                    , dimensionless_  &
                                                    , refValues_ & 
#endif   
#if defined(CAHNHILLIARD)
                                                    , multiphase_ &
#endif
                                                    , monitors, &
                                                      elapsedTime, &
                                                      CPUTime   )
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            use SMConstants
            use FTAssertions
            USE HexMeshClass
#ifdef _HAS_MPI_
            use mpi
#endif
            use PhysicsStorage
            use FluidData
            use MonitorsClass
            IMPLICIT NONE
            class(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)    :: multiphase_
#endif
            type(Monitor_t),        intent(in)    :: monitors
            real(kind=RP),             intent(in) :: elapsedTime
            real(kind=RP),             intent(in) :: CPUTime
!
!           ---------------
!           Local variables
!           ---------------
!
#if defined(NAVIERSTOKES)
!           With "optimized partition", every rank should own about 1/nProcs of the elements of
!           each level (counted with MLevelwN, the level the partition is built from).
            integer,       parameter            :: N_LEVELS = 2
            integer,       parameter            :: expectedTotal(N_LEVELS) = [320, 192]
            real(kind=RP), parameter            :: BALANCE_TOL = 0.1_RP   ! allowed deviation from the mean, per level
            character(len=*), parameter         :: testName = "MLRK level-aware partition"
            type(FTAssertionsManager), pointer  :: sharedManager
            integer                             :: eID, lev, r, ierr
            integer                             :: myRank, nRanks   ! MPI_Process is not shared with this library: ask MPI
            integer                             :: localCount(N_LEVELS)
            integer, allocatable                :: allCount(:,:)
            real(kind=RP)                       :: mean
            logical                             :: balanced(N_LEVELS)
            character(len=64)                   :: msg

            myRank = 0
            nRanks = 1
#ifdef _HAS_MPI_
            call mpi_comm_rank(MPI_COMM_WORLD, myRank, ierr)
            call mpi_comm_size(MPI_COMM_WORLD, nRanks, ierr)
#endif

            localCount = 0
            do eID = 1, size(mesh % elements)
               lev = mesh % elements(eID) % MLevelwN
               localCount(lev) = localCount(lev) + 1
            end do

            allocate(allCount(N_LEVELS, 0:nRanks-1))
            allCount = 0
            allCount(:, myRank) = localCount
#ifdef _HAS_MPI_
            if ( nRanks > 1 ) then
               call mpi_allreduce(MPI_IN_PLACE, allCount, size(allCount), MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
            end if
#endif

            do lev = 1, N_LEVELS
               mean = real(sum(allCount(lev,:)), RP) / nRanks
               balanced(lev) = all( abs(allCount(lev,:) - mean) <= BALANCE_TOL * mean + 1.0_RP )
            end do

            if ( myRank == 0 ) then
               write(6,'(A)') "Elements per rank and level (MLevelwN):"
               write(6,'(A8,*(A10,I1))') "rank", ("    level ", lev, lev = 1, N_LEVELS)
               do r = 0, nRanks-1
                  write(6,'(I8,*(I11))') r, allCount(:, r)
               end do
               write(6,'(A8,*(I11))') "total", (sum(allCount(lev,:)), lev = 1, N_LEVELS)
            end if

            call initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            do lev = 1, N_LEVELS
               write(msg,'(A,I2)') "Total elements with MLevelwN =", lev
               call FTAssertEqual(expectedValue = expectedTotal(lev), actualValue = sum(allCount(lev,:)), msg = trim(msg))
               write(msg,'(A,I2,A)') "Level", lev, " balanced over the ranks"
               call FTAssert(test = balanced(lev), msg = trim(msg))
            end do

            if ( myRank == 0 ) call sharedManager % summarizeAssertions(title = testName, iUnit = 6)

            if ( sharedManager % numberOfAssertionFailures() == 0 ) then
               if ( myRank == 0 ) write(6,*) testName, " ... Passed"
            else
               if ( myRank == 0 ) write(6,*) testName, " ... Failed"
            end if
!           The root's report must be out before any rank stops: an error stop on one rank kills the others
            flush(6)
#ifdef _HAS_MPI_
            if ( nRanks > 1 ) call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
            if ( sharedManager % numberOfAssertionFailures() /= 0 ) error stop 99

            call finalizeSharedAssertionsManager
            call detachSharedAssertionsManager
#endif



         END SUBROUTINE UserDefinedFinalize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedTermination
!
!        -----------------------------------------------
!        Called at the the end of the main driver after 
!        everything else is done.
!        -----------------------------------------------
!
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
      