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
#include "Includes.h"
module ProblemFileFunctions
   implicit none

   abstract interface
      subroutine UserDefinedStartup_f
      end subroutine UserDefinedStartup_f

      subroutine UserDefinedFinalSetup_f(mesh &
#ifdef FLOW
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                     , multiphase_ &
#endif
                                     )
         use HexMeshClass
         use FluidData
         implicit none
         class(HexMesh)                      :: mesh
#ifdef FLOW
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      end subroutine UserDefinedFinalSetup_f

      subroutine UserDefinedInitialCondition_f(mesh &
#ifdef FLOW
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                     , multiphase_ &
#endif
                                     )
         use smconstants
         use physicsstorage
         use hexmeshclass
         use fluiddata
         implicit none
         class(hexmesh)                      :: mesh
#ifdef FLOW
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      end subroutine UserDefinedInitialCondition_f
#ifdef FLOW

      subroutine UserDefinedState_f(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         implicit none
         real(kind=RP), intent(in)          :: x(NDIM)
         real(kind=RP), intent(in)          :: t
         real(kind=RP), intent(in)          :: nHat(NDIM)
         real(kind=RP), intent(inout)       :: Q(NCONS)
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
      end subroutine UserDefinedState_f

      subroutine UserDefinedGradVars_f(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
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
      end subroutine UserDefinedGradVars_f

      subroutine UserDefinedNeumann_f(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         implicit none
         real(kind=RP), intent(in)    :: x(NDIM)
         real(kind=RP), intent(in)    :: t
         real(kind=RP), intent(in)    :: nHat(NDIM)
         real(kind=RP), intent(in)    :: Q(NCONS)
         real(kind=RP), intent(in)    :: U_x(NGRAD)
         real(kind=RP), intent(in)    :: U_y(NGRAD)
         real(kind=RP), intent(in)    :: U_z(NGRAD)
         real(kind=RP), intent(inout) :: flux(NCONS)
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
      end subroutine UserDefinedNeumann_f

#endif
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine UserDefinedPeriodicOperation_f(mesh, time, dt, Monitors)
         use SMConstants
         use HexMeshClass
         use MonitorsClass
         implicit none
         class(HexMesh)               :: mesh
         real(kind=RP)                :: time
         real(kind=RP)                :: dt
         type(Monitor_t), intent(in) :: monitors
      end subroutine UserDefinedPeriodicOperation_f
!
!////////////////////////////////////////////////////////////////////////
!
#ifdef FLOW
      subroutine UserDefinedSourceTermNS_f(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
,multiphase_ &
#endif
)
         use SMConstants
         use HexMeshClass
         use FluidData
         use PhysicsStorage
         implicit none
         real(kind=RP),             intent(in)  :: x(NDIM)
         real(kind=RP),             intent(in)  :: Q(NCONS)
         real(kind=RP),             intent(in)  :: time
         real(kind=RP),             intent(inout) :: S(NCONS)
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      end subroutine UserDefinedSourceTermNS_f
#endif
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine UserDefinedFinalize_f(mesh, time, iter, maxResidual &
#ifdef FLOW
                                                 , thermodynamics_ &
                                                 , dimensionless_  &
                                                 , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                                 , multiphase_ &
#endif
                                                 , monitors, &
                                                   elapsedTime, &
                                                   CPUTime   )
         use SMConstants
         use HexMeshClass
         use FluidData
         use MonitorsClass
         implicit none
         class(HexMesh)                        :: mesh
         real(kind=RP)                         :: time
         integer                               :: iter
         real(kind=RP)                         :: maxResidual
#ifdef FLOW
         type(Thermodynamics_t), intent(in)    :: thermodynamics_
         type(Dimensionless_t),  intent(in)    :: dimensionless_
         type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)    :: multiphase_
#endif
         type(Monitor_t),        intent(in)    :: monitors
         real(kind=RP),             intent(in) :: elapsedTime
         real(kind=RP),             intent(in) :: CPUTime
      end subroutine UserDefinedFinalize_f

      subroutine UserDefinedTermination_f
         implicit none
      end subroutine UserDefinedTermination_f
   end interface

end module ProblemFileFunctions

         subroutine UserDefinedStartup
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            implicit none
         end subroutine UserDefinedStartup
!
!////////////////////////////////////////////////////////////////////////
!
         subroutine UserDefinedFinalSetup(mesh &
#ifdef FLOW
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ &
#endif
#ifdef CAHNHILLIARD
                                        , multiphase_ &
#endif
                                        )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            use HexMeshClass
            use PhysicsStorage
            use FluidData
            implicit none
            class(HexMesh)                      :: mesh
#ifdef FLOW
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
         end subroutine UserDefinedFinalSetup
!
!////////////////////////////////////////////////////////////////////////
!
         subroutine UserDefinedInitialCondition(mesh &
#ifdef FLOW
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ &
#endif
#ifdef CAHNHILLIARD
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
#ifdef FLOW
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            real(KIND=RP) :: x(3)
            integer       :: i, j, k, eID
            real(KIND=RP) :: rho , u , v , w , p
            real(KIND=RP) :: L, u_0, rho_0, p_0
            integer       :: Nx, Ny, Nz

!
!           ---------------------------------------
!           Navier-Stokes default initial condition
!           ---------------------------------------
!
#if defined(NAVIERSTOKES)
            L     = 1.0_RP
            u_0   = 1.0_RP
            rho_0 = 1.0_RP
            p_0   = 100.0_RP

            associate( gamma => thermodynamics_ % gamma )
            do eID = 1, size(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               do k = 0, Nz
                  do j = 0, Ny
                     do i = 0, Nx

                         x = mesh % elements(eID) % geom % x(:,i,j,k)

                         rho =  rho_0
                         u   =  u_0 * sin(x(1)/L) * cos(x(2)/L) * cos(x(3)/L)
                         v   = -u_0 * cos(x(1)/L) * sin(x(2)/L) * cos(x(3)/L)
                         w   =  0.0_RP
                         p   =  p_0 + rho_0 / 16.0_RP * (                               &
                               cos(2.0_RP*x(1)/L)*cos(2.0_RP*x(3)/L) +                  &
                               2.0_RP*cos(2.0_RP*x(2)/L) + 2.0_RP*cos(2.0_RP*x(1)/L) +  &
                               cos(2.0_RP*x(2)/L)*cos(2.0_RP*x(3)/L)                    &
                               )

                         mesh % elements(eID) % storage % Q(1,i,j,k) = rho
                         mesh % elements(eID) % storage % Q(2,i,j,k) = rho*u
                         mesh % elements(eID) % storage % Q(3,i,j,k) = rho*v
                         mesh % elements(eID) % storage % Q(4,i,j,k) = rho*w
                         mesh % elements(eID) % storage % Q(5,i,j,k) = p / (gamma - 1.0_RP) + 0.5_RP * rho * (u*u + v*v + w*w)

                     end do
                  end do
               end do

            end do
            end associate
#endif
!
!           ------------------------------------------------------
!           Incompressible Navier-Stokes default initial condition
!           ------------------------------------------------------
!
#if defined(INCNS)
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx
                  mesh % elements(eID) % storage % q(:,i,j,k) = [1.0_RP, 1.0_RP,0.0_RP,0.0_RP,0.0_RP]
               end do;        end do;        end do
               end associate
            end do
#endif

!
!           ---------------------------------------
!           Cahn-Hilliard default initial condition
!           ---------------------------------------
!
#ifdef CAHNHILLIARD
            call random_seed()

            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          Ny => mesh % elements(eID) % Nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               associate(e => mesh % elements(eID) % storage)
               call random_number(e % c)
               e % c = 2.0_RP * (e % c - 0.5_RP)
               end associate
               end associate
            end do
#endif
         end subroutine UserDefinedInitialCondition

#ifdef FLOW
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
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

         subroutine UserDefinedNeumann1(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)    :: x(NDIM)
            real(kind=RP), intent(in)    :: t
            real(kind=RP), intent(in)    :: nHat(NDIM)
            real(kind=RP), intent(in)    :: Q(NCONS)
            real(kind=RP), intent(in)    :: U_x(NGRAD)
            real(kind=RP), intent(in)    :: U_y(NGRAD)
            real(kind=RP), intent(in)    :: U_z(NGRAD)
            real(kind=RP), intent(inout) :: flux(NCONS)
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedNeumann1
#endif
!
!////////////////////////////////////////////////////////////////////////
!
         subroutine UserDefinedPeriodicOperation(mesh, time, dt, Monitors)
!
!           ----------------------------------------------------------
!           Called before every time-step to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use MonitorsClass
            implicit none
            class(HexMesh)               :: mesh
            real(kind=RP)                :: time
            real(kind=RP)                :: dt
            type(Monitor_t), intent(in) :: monitors

         end subroutine UserDefinedPeriodicOperation
!
!////////////////////////////////////////////////////////////////////////
!
#ifdef FLOW
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
                                          , multiphase_ &
#endif
)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            use HexMeshClass
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: Q(NCONS)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(inout) :: S(NCONS)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            integer  :: i, j, k, eID
!
!           Usage example
!           -------------
!           S(:) = x(1) + x(2) + x(3) + time
            S  = 0.0_RP

         end subroutine UserDefinedSourceTermNS
#endif
!
!////////////////////////////////////////////////////////////////////////
!
         subroutine UserDefinedFinalize(mesh, time, iter, maxResidual &
#ifdef FLOW
                                                    , thermodynamics_ &
                                                    , dimensionless_  &
                                                    , refValues_ &
#endif
#ifdef CAHNHILLIARD
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
            use HexMeshClass
            use PhysicsStorage
            use FluidData
            use MonitorsClass

            implicit none
            class(HexMesh)                     :: mesh
            real(kind=RP)                      :: time
            integer                            :: iter
            real(kind=RP)                      :: maxResidual
#ifdef FLOW
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in) :: multiphase_
#endif
            type(Monitor_t),        intent(in) :: monitors
            real(kind=RP),          intent(in) :: elapsedTime
            real(kind=RP),          intent(in) :: CPUTime
!
!           ---------------
!           Local variables
!           ---------------
!
#if defined(NAVIERSTOKES)
            CHARACTER(LEN=29)                  :: testName = "Taylor-Green vortex (SVV-LES)"
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            real(kind=RP), parameter           :: kinEn = 0.12499999968280391_RP
            real(kind=RP), parameter           :: kinEnRate = -7.6138351130710094e-8_RP
            real(kind=RP), parameter           :: SVVdiss = 7.7873340678919037e-8_RP
            real(kind=RP), parameter           :: res(5) = [ 9.6089299201416065e-5_RP,   &
                                                             0.12704234777396151_RP,   &
                                                             0.12704234731965755_RP,   &
                                                             0.25000646761649736_RP,   &
                                                             0.62824891369699054_RP    ]


            call initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            call FTAssertEqual(expectedValue = res(1) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(1,1) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "Continuity residual")

            call FTAssertEqual(expectedValue = res(2), &
                               actualValue   = monitors % residuals % values(2,1), &
                               tol           = 1.0e-7_RP, &
                               msg           = "X-momentum residual")

            call FTAssertEqual(expectedValue = res(3), &
                               actualValue   = monitors % residuals % values(3,1), &
                               tol           = 1.0e-7_RP, &
                               msg           = "Y-momentum residual")

            call FTAssertEqual(expectedValue = res(4), &
                               actualValue   = monitors % residuals % values(4,1), &
                               tol           = 1.0e-7_RP, &
                               msg           = "Z-momentum residual")

            call FTAssertEqual(expectedValue = res(5), &
                               actualValue   = monitors % residuals % values(5,1), &
                               tol           = 1.0e-7_RP, &
                               msg           = "Energy residual")

            call FTAssertEqual(expectedValue = kinEn, &
                               actualValue   = monitors % volumeMonitors(1) % values(1,1), &
                               tol           = 1.0e-7_RP, &
                               msg           = "Kinetic energy")

            call FTAssertEqual(expectedValue = kinEnRate + 1.0_RP, &
                               actualValue   = monitors % volumeMonitors(2) % values(1,1) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "Kinetic energy rate")

            call FTAssertEqual(expectedValue = svvDiss + 1.0_RP, &
                               actualValue   = monitors % volumeMonitors(3) % values(1,1) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "SVV dissipation")

            call sharedManager % summarizeAssertions(title = testName,iUnit = 6)

            if ( sharedManager % numberOfAssertionFailures() == 0 ) then
               write(6,*) testName, " ... Passed"
               write(6,*) "This test checks if:"
               write(6,*) "   - The residuals are exact to 1e-7"
               write(6,*) "   - The kinetic energy and the SVV dissipation are also exact (1e-7)"
            else
               write(6,*) testName, " ... Failed"
               error stop 99
            end if
            write(6,*)

            call finalizeSharedAssertionsManager
            call detachSharedAssertionsManager
#endif

         end subroutine UserDefinedFinalize
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine UserDefinedTermination
!
!        -----------------------------------------------
!        Called at the the end of the main driver after
!        everything else is done.
!        -----------------------------------------------
         implicit none
      end subroutine UserDefinedTermination