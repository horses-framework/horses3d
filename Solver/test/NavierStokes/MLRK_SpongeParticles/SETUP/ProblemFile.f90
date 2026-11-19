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
            REAL(KIND=RP) :: rho , u , v , w , p, p_0

#if defined(NAVIERSTOKES)
!           Taylor-Green vortex
            associate( gamma => thermodynamics_ % gamma )
            p_0 = 1.0_RP / (gamma * dimensionless_ % Mach**2)
            do eID = 1, size(mesh % elements)
               associate( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  x   = e % geom % x(:,i,j,k)
                  rho = 1.0_RP
                  u   =  sin(x(1)) * cos(x(2)) * cos(x(3))
                  v   = -cos(x(1)) * sin(x(2)) * cos(x(3))
                  w   = 0.0_RP
                  p   = p_0 + (cos(2.0_RP*x(1)) + cos(2.0_RP*x(2))) * (cos(2.0_RP*x(3)) + 2.0_RP) / 16.0_RP
                  e % storage % Q(:,i,j,k) = [rho, rho*u, rho*v, rho*w, p/(gamma - 1.0_RP) + 0.5_RP*rho*(u*u + v*v + w*w)]
               end do ; end do ; end do
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
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
!           The references are per-level sums of |Q| over all nodes, at the final time. A sponge or particle source
!           missing in, or leaking into, any level changes them.
            integer,       parameter            :: N_LEVELS = 4
            integer,       parameter            :: expectedCount(N_LEVELS) = [192, 128, 128, 64]
            real(kind=RP), parameter            :: TOL = 1.0e-10_RP
            real(kind=RP), parameter            :: refSums(NCONS,N_LEVELS) = reshape( [ &
               1.5360023362154382E+03_RP, 2.9517207585522686E+02_RP, 3.8984217634476374E+02_RP, &
               2.8518480959374397E+01_RP, 4.2884502960687043E+05_RP, &
               1.0239946460282946E+03_RP, 3.2271032302764138E+02_RP, 1.3399680550875107E+02_RP, &
               1.8830332240856908E+01_RP, 2.8567895450999722E+05_RP, &
               1.0239830088911115E+03_RP, 2.4683023719963697E+02_RP, 2.3276660301521116E+02_RP, &
               2.1799397543263407E+01_RP, 2.8583852788984956E+05_RP, &
               2.1953264432484775E+04_RP, 3.2523190700025948E+03_RP, 7.2339358571605790E+03_RP, &
               2.8864693374861503E+02_RP, 6.1321429594813148E+06_RP ], [NCONS,N_LEVELS] )
            character(len=*), parameter         :: testName = "MLRK with sponge and particles"
            type(FTAssertionsManager), pointer  :: sharedManager
            real(kind=RP)                       :: sums(NCONS,N_LEVELS)
            integer                             :: eID, lev, ivar, nElemLevel(N_LEVELS)
            logical                             :: spongeInLevel(N_LEVELS), particlesInLevel(N_LEVELS)
            character(len=64)                   :: msg

            sums             = 0.0_RP
            nElemLevel       = 0
            spongeInLevel    = .false.
            particlesInLevel = .false.
            do eID = 1, size(mesh % elements)
               associate( e => mesh % elements(eID) )
               lev = e % MLevel
               nElemLevel(lev) = nElemLevel(lev) + 1
               do ivar = 1, NCONS
                  sums(ivar,lev) = sums(ivar,lev) + sum(abs(e % storage % Q(ivar,:,:,:)))
               end do
               if ( maxval(e % storage % intensitySponge) > 0.0_RP ) spongeInLevel(lev) = .true.
               if ( maxval(abs(e % storage % S_NSP)) > 0.0_RP ) particlesInLevel(lev) = .true.
               end associate
            end do

            do lev = 1, N_LEVELS
               write(STD_OUT,'(A,I2,A,I5,A,L2,A,L2)') "Level", lev, ": elements =", nElemLevel(lev), &
                     ", sponge:", spongeInLevel(lev), ", particle source:", particlesInLevel(lev)
               write(STD_OUT,'(A,5(ES24.16,A))') "   sum |Q|: ", (sums(ivar,lev), "_RP, ", ivar = 1, NCONS)
            end do

            call initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            call FTAssertEqual(expectedValue = N_LEVELS, actualValue = mesh % MLRK % maxLevel, &
                               msg = "Number of MLRK levels")
            do lev = 1, N_LEVELS
               write(msg,'(A,I2)') "Elements in level", lev
               call FTAssertEqual(expectedValue = expectedCount(lev), actualValue = nElemLevel(lev), msg = trim(msg))
               write(msg,'(A,I2)') "Test setup: sponge elements in level", lev
               call FTAssert(test = spongeInLevel(lev), msg = trim(msg))
               write(msg,'(A,I2)') "Test setup: particle source in level", lev
               call FTAssert(test = particlesInLevel(lev), msg = trim(msg))
               do ivar = 1, NCONS
                  write(msg,'(A,I2,A,I2)') "Sum of |Q| in level", lev, ", variable", ivar
                  call FTAssertEqual(expectedValue = refSums(ivar,lev), actualValue = sums(ivar,lev), &
                                     tol = TOL, msg = trim(msg))
               end do
            end do

            call sharedManager % summarizeAssertions(title = testName, iUnit = 6)

            if ( sharedManager % numberOfAssertionFailures() == 0 ) then
               write(6,*) testName, " ... Passed"
            else
               write(6,*) testName, " ... Failed"
               error stop 99
            end if
            write(6,*)

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
      