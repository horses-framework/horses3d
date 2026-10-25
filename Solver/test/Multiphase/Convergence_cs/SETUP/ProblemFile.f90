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
   
      SUBROUTINE UserDefinedFinalSetup_f(mesh &
#ifdef FLOW
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ & 
#endif
#ifdef CAHNHILLIARD
                                     , multiphase_ &
#endif
                                     )
         USE HexMeshClass
         use FluidData
         IMPLICIT NONE
         CLASS(HexMesh)                      :: mesh
#ifdef FLOW
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      END SUBROUTINE UserDefinedFinalSetup_f

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
      end subroutine UserDefinedState_f

      subroutine UserDefinedNeumann_f(x, t, nHat, U_x, U_y, U_z)
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
      end subroutine UserDefinedNeumann_f
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedPeriodicOperation_f(mesh, time, dt, Monitors)
         use SMConstants
         USE HexMeshClass
         use MonitorsClass
         IMPLICIT NONE
         CLASS(HexMesh)               :: mesh
         REAL(KIND=RP)                :: time
         REAL(KIND=RP)                :: dt
         type(Monitor_t), intent(in) :: monitors
      END SUBROUTINE UserDefinedPeriodicOperation_f
!
!//////////////////////////////////////////////////////////////////////// 
! 
#ifdef FLOW
      subroutine UserDefinedSourceTermNS_f(x, Q, time, S, thermodynamics_, dimensionless_, refValues_&
#ifdef CAHNHILLIARD
,multiphase_ &
#endif
)
         use SMConstants
         USE HexMeshClass
         use FluidData
         use PhysicsStorage
         IMPLICIT NONE
         real(kind=RP),             intent(in)  :: x(NDIM)
         real(kind=RP),             intent(in)  :: Q(NCONS)
         real(kind=RP),             intent(in)  :: time
         real(kind=RP),             intent(inout) :: S(NCONS)
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#ifdef CAHNHILLIARD
         type(Multiphase_t),  intent(in)  :: multiphase_
#endif
      end subroutine UserDefinedSourceTermNS_f
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedFinalize_f(mesh, time, iter, maxResidual &
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
         USE HexMeshClass
         use FluidData
         use MonitorsClass
         IMPLICIT NONE
         CLASS(HexMesh)                        :: mesh
         REAL(KIND=RP)                         :: time
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
      END SUBROUTINE UserDefinedFinalize_f

      SUBROUTINE UserDefinedTermination_f
         implicit none
      END SUBROUTINE UserDefinedTermination_f
   end interface
   
end module ProblemFileFunctions

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
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            CLASS(HexMesh)                      :: mesh
#ifdef FLOW
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
         END SUBROUTINE UserDefinedFinalSetup
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
!           local variables
!           ---------------
!
            integer        :: eid, i, j, k
            real(kind=RP)  :: c, u, v, w, p, phi
            real(kind=RP)  :: x, y, z
!
!           ------------------------------------------------------
!           Incompressible Navier-Stokes default initial condition
!           ------------------------------------------------------
!
#if defined(MULTIPHASE)
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  x = mesh % elements(eID) % geom % x(IX,i,j,k)               
                  y = mesh % elements(eID) % geom % x(IY,i,j,k)               
                  z = mesh % elements(eID) % geom % x(IZ,i,j,k)               

                  phi = 0.0_RP
                  c = 0.5_RP
      
                  u = 0.0_RP
                  v = 0.0_RP
                  w = 0.0_RP

                  p = 2*sin(PI*x)*sin(PI*y)
   
                  mesh % elements(eID) % storage % q(:,i,j,k) = [c,u,v,w,p] 
               end do;        end do;        end do
               end associate
            end do
#endif

         end subroutine UserDefinedInitialCondition
#ifdef FLOW
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

         subroutine UserDefinedGradVars1(x, t, nHat, Q, U, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)          :: x(NDIM)
            real(kind=RP), intent(in)          :: t
            real(kind=RP), intent(in)          :: nHat(NDIM)
            real(kind=RP), intent(in)          :: Q(NCONS)
            real(kind=RP), intent(inout)       :: U(NGRAD)
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
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors)
!
!           ----------------------------------------------------------
!           Called before every time-step to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)               :: mesh
            REAL(KIND=RP)                :: time
            REAL(KIND=RP)                :: dt
            type(Monitor_t), intent(in) :: monitors
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
#ifdef FLOW
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_&
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
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: Q(NCONS)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(inout) :: S(NCONS)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
!
!           ---------------
!           Local variables
!           ---------------
!
            real(kind=RP)  :: rho0, eta = 0.001_RP
            real(kind=RP)  :: rho1, rho2, cs1, cs2

            rho1 = dimensionless_ % rho(2)
            rho2 = dimensionless_ % rho(1)
            cs1 = sqrt(thermodynamics_ % c02(2))
            cs2 = sqrt(thermodynamics_ % c02(1))

            S(IMC) = 0.5_RP*cos(time)*cos(PI*x(IX))*cos(PI*x(IY)) &  !diff(c,t)                          
                     +PI* &
                     (4.0_RP*sin(time)*(cos(PI*x(IX)))**2*(cos(PI*x(IY)))**2 &
                     -1.0_RP*sin(time)*(cos(PI*x(IX)))**2  &
                     -1.0_RP*sin(time)*(cos(PI*x(IY)))**2  &
                     +2.0_RP*cos(PI*x(IX))*cos(PI*x(IY))   &
                     )*(sin(time)) &  !+diff(c*u,x)+diff(c*v,y)
                     +(PI**2)*(multiphase_ % M0*multiphase_ % sigma/multiphase_ % eps) &
                     *(                                 &
                      3.0_RP*(PI*multiphase_ % eps)**2 &
                      +108.0_RP*(sin(time)*sin (PI*x(IX))*sin(PI*x(IY)))**2 &
                      -72.0_RP*(sin(time)*sin(PI*x(IX)))**2 &
                      -72.0_RP*(sin(time)*sin(PI*x(IY)))**2 &
                      +36.0_RP*(sin(time)**2) &
                      -12.0_RP &
                      ) &
                     *sin(time)*cos(PI*x(IX))*cos(PI*x(IY))  !-M0*(diff(μ,x,2)+diff(μ,y,2)) 

            S(IMSQRHOU) = 1.0_RP*( &
                           -1.5_RP*rho1*sin(time)*cos(PI*x(IX))*cos(PI*x(IY)) + rho1 & 
                           +1.5_RP*rho2*sin(time)*cos(PI*x(IX))*cos(PI*x(IY)) + rho2 &
                          ) &
                          *sin(PI*x(IX))*cos(time)*cos(PI*x(IY)) & !sqrt(ρ)*diff(sqrt(ρ)*u,t)                                                                                                                                                        
                         + PI * ( &
                           - 6.0_RP * rho1 * sin(time) * (cos(PI*x(IX)))**2 * (cos(PI*x(IY)))**3 &
                           + 2.0_RP * rho1 * sin(time) * (cos(PI*x(IX)))**2 * cos(PI*x(IY)) &
                           + 1.0_RP * rho1 * sin(time) * cos(PI*x(IY))**3 &
                           + 4.0_RP * rho1 * cos(PI*x(IX)) * (cos(PI*x(IY)))**2 &
                           - 1.0_RP * rho1 * cos(PI*x(IX)) &
                           + 6.0_RP * rho2 * sin(time) * (cos(PI*x(IX)))**2 * (cos(PI*x(IY)))**3 &
                           - 2.0_RP * rho2 * sin(time) * (cos(PI*x(IX)))**2 * cos(PI*x(IY)) &
                           - 1.0_RP * rho2 * sin(time) * cos(PI*x(IY))**3 &
                           + 4.0_RP * rho2 * cos(PI*x(IX)) * (cos(PI*x(IY)))**2 &
                           - 1.0_RP * rho2 * cos(PI*x(IX)) &
                             ) * (sin(time))**2 * sin(PI*x(IX)) & !0.5*(diff(ρ*u*u,x)+diff(ρ*u*v,y))
                          -PI*( &
                           rho1*(sin(time)*cos(PI*x(IX))*cos(PI*x(IY)) - 1.0_RP) &
                          -rho2*(sin(time)*cos(PI*x(IX))*cos(PI*x(IY)) + 1.0_RP) &
                           ) &
                          *(sin(time))**2*sin(PI*x(IX))*cos(PI*x(IX))*cos(2*PI*x(IY)) &    !0.5*ρ*(u*diff(u,x)+v*diff(u,y))
                          -0.5_RP*PI*multiphase_ % sigma / multiphase_ % eps  & 
                          *(sin(time)*cos(PI*x(IX))*cos(PI*x(IY)) + 1.0_RP) &
                          *(1.5_RP*(PI*multiphase_ % eps)**2  + 18.0_RP*(sin(time)*cos(PI*x(IX))*cos(PI*x(IY)))**2 - 6.0_RP) &
                          *sin(time)*sin(PI*x(IX))*cos(PI*x(IY)) &  !+ c *diff(μ,x)
                          + 2.0_RP*PI*sin(PI*x(IY))*cos(time)*cos(PI*x(IX)) & !+diff(p,x)
                          + 8.0_RP*(PI**2) *eta*sin(time)*sin(PI*x(IX))*cos(PI*x(IY)) !- diff(η*(diff(u,x) + diff(u,x)),x) - diff(η*(diff(u,y) + diff(v,x)),y)

            S(IMSQRHOV) = 1.0_RP*( &
                           -1.5_RP*rho1*sin(time)*cos(PI*x(IX))*cos(PI*x(IY)) + rho1 & 
                           +1.5_RP*rho2*sin(time)*cos(PI*x(IX))*cos(PI*x(IY)) + rho2 &
                          ) &
                          *sin(PI*x(IY))*cos(time)*cos(PI*x(IX)) & !sqrt(ρ)*diff(sqrt(ρ)*v,t)        
                          +PI * ( &
                             -6.0_RP * rho1 * sin(time) * (cos(PI*x(IX)))**3 * (cos(PI*x(IY)))**2 &
                             + 1.0_RP * rho1 * sin(time) * (cos(PI*x(IX)))**3 &
                             + 2.0_RP * rho1 * sin(time) * cos(PI*x(IX)) * (cos(PI*x(IY)))**2 &
                             + 4.0_RP * rho1 * (cos(PI*x(IX)))**2 * cos(PI*x(IY)) &
                             - 1.0_RP * rho1 * (cos(PI*x(IY))) &
                             + 6.0_RP * rho2 * sin(time) * (cos(PI*x(IX)))**3 * (cos(PI*x(IY)))**2 &
                             - 1.0_RP * rho2 * sin(time) * (cos(PI*x(IX)))**3 &
                             - 2.0_RP * rho2 * sin(time) * cos(PI*x(IX)) * (cos(PI*x(IY)))**2 &
                             + 4.0_RP * rho2 * (cos(PI*x(IX)))**2 * cos(PI*x(IY)) &
                             - 1.0_RP * rho2 * (cos(PI*x(IY))) ) &
                           * (sin(time))**2 * sin(PI*x(IY)) & !0.5*(diff(ρ*v*u,x)+diff(ρ*v*v,y))
                          -PI*( &
                           rho1*(sin(time)*cos(PI*x(IX))*cos(PI*x(IY)) - 1.0_RP) &
                          -rho2*(sin(time)*cos(PI*x(IX))*cos(PI*x(IY)) + 1.0_RP) &
                           ) &
                          *(sin(time))**2*sin(PI*x(IY))*cos(2.0_RP*PI*x(IX))*cos(PI*x(IY)) &    !0.5*ρ*(u*diff(v,x)+v*diff(v,y))
                          -0.5_RP*PI*multiphase_ % sigma / multiphase_ % eps  & 
                          *(sin(time)*cos(PI*x(IX))*cos(PI*x(IY)) + 1.0_RP) &
                          *(1.5_RP*(PI*multiphase_ % eps)**2  + 18.0_RP*(sin(time)*cos(PI*x(IX))*cos(PI*x(IY)))**2 - 6.0_RP) &
                          *sin(time)*sin(PI*x(IY))*cos(PI*x(IX)) &  !+ c *diff(μ,y)
                          + 2.0_RP*PI*sin(PI*x(IX))*cos(time)*cos(PI*x(IY)) & !+diff(p,y)
                          + 8.0_RP *(PI**2)*eta*sin(time)*sin(PI*x(IY))*cos(PI*x(IX)) !- diff(η*(diff(v,x) + diff(u,y)),x) - diff(η*(diff(v,y) + diff(v,y)),y)

            
            S(IMSQRHOW) = 0.0_RP
            
            S(IMP) = sin(PI*x(IX))*sin(PI*x(IY))*sin(time)*(-2.0_RP) & !diff(p,t) 
                     -0.5_RP * PI * ( &
                         cs1 * (sin(time) * cos(PI*x(IX)) * cos(PI*x(IY)) - 1.0_RP) &
                       - cs2 * (sin(time) * cos(PI*x(IX)) * cos(PI*x(IY)) + 1.0_RP) )**2 &
                       * ( &
                         rho1 * (sin(time) * cos(PI*x(IX)) * cos(PI*x(IY)) - 1.0_RP) &
                       - rho2 * (sin(time) * cos(PI*x(IX)) * cos(PI*x(IY)) + 1.0_RP) ) &
                       * sin(time) * cos(PI*x(IX)) * cos(PI*x(IY)) !+ρ*cs^2*(diff(u,x)+diff(v,y))
#endif                  
   
         end subroutine UserDefinedSourceTermNS
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual &
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
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use MonitorsClass
            use FTAssertions
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
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
            real(kind=RP)  :: x, y, z, c, locErr(5), phi, u, v, w, p, rho, rho0
            real(kind=RP), parameter  :: saved_errors(5) = [4.8755376700349018E-06_RP, &
                                          6.8771002939893577E-05_RP, &
                                          6.8771002939893577E-05_RP, &
                                          7.4217194511428453E-17_RP, & 
                                          6.1903872140001148E-04_RP]
            integer        :: i, j,k, eID 
            CHARACTER(LEN=49)                  :: testName           = "Multiphase convergence non constant sound speed"
            real(kind=RP)  :: error(5)
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
            real(kind=RP), parameter   :: w_LGL(0:5) = [0.066666666666667_RP, &
                                                      0.378474956297847_RP, &
                                                      0.554858377035486_RP, &
                                                      0.554858377035486_RP, &
                                                      0.378474956297847_RP, &
                                                      0.066666666666667_RP ]
!
!           *********************************
!           Check the L-inf norm of the error
!           *********************************
!
#ifdef MULTIPHASE
            rho0 = dimensionless_ % rho(2)
#endif
            error = 0.0_RP
            do eID = 1, mesh % no_of_elements
               associate(e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  x = e % geom % x(IX,i,j,k)
                  y = e % geom % x(IY,i,j,k)
                  z = e % geom % x(IZ,i,j,k)

                  phi = cos(PI*x)*cos(PI*y)*sin(time)
                  c = 0.5_RP*(1.0_RP + phi)

                  rho = c + rho0*(1.0_RP-c)

                  u = 2*sin(PI*x)*cos(PI*y)*sin(time)
                  v = 2*cos(PI*x)*sin(PI*y)*sin(time)
                  w = 0.0_RP
                  p = 2*sin(PI*x)*sin(PI*y)*cos(time)

                  locErr = e % storage % Q(:,i,j,k) - [c,sqrt(rho)*u,sqrt(rho)*v,sqrt(rho)*w,p]

                  error = error + e % geom % jacobian(i,j,k)*locErr**2*w_LGL(i)*w_LGL(j)*2.0_RP

               end do            ; end do ; end do
               end associate
            end do

            error = sqrt(error)

            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
            CALL FTAssertEqual(expectedValue = saved_errors(1)+1.0_RP, &
                               actualValue   = error(1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Concentration error")


            CALL FTAssertEqual(expectedValue = saved_errors(2)+1.0_RP, &
                               actualValue   = error(2)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "X-Momentum error")

            CALL FTAssertEqual(expectedValue = saved_errors(3)+1.0_RP, &
                               actualValue   = error(3)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Y-Momentum error")

            CALL FTAssertEqual(expectedValue = saved_errors(4)+1.0_RP, &
                               actualValue   = error(4)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Z-Momentum error")

            CALL FTAssertEqual(expectedValue = saved_errors(5)+1.0_RP, &
                               actualValue   = error(5)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Pressure error")

            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
               WRITE(6,*) "This test case has no expected solution yet, only checks the residual after 100 iterations."
            ELSE
               WRITE(6,*) testName, " ... Failed"
               WRITE(6,*) "NOTE: Failure is expected when the max eigenvalue procedure is changed."
               WRITE(6,*) "      If that is done, re-compute the expected values and modify this procedure"
                error stop 99
            END IF 
            WRITE(6,*)
            
            CALL finalizeSharedAssertionsManager
            CALL detachSharedAssertionsManager




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
      