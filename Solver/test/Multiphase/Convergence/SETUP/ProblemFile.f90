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

                  p = 2*sin(PI*x)*sin(PI*z)
   
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

            rho0 = dimensionless_ % rho(2)

            S(IMC) =  (multiphase_ % eps*PI*cos(PI*X(IX))**2-multiphase_ % eps*PI*cos(PI*X(IZ))**2-multiphase_ % eps*PI*cos(PI*&
                      X(IX))**2*cos(time)**2+multiphase_ % eps*PI*cos(PI*X(IZ))**2*cos(time)**2+multiphase_ % eps*cos(PI*X(IX))*cos&
                      (PI*X(IZ))*cos(time)*(1.0D0/2.0D0)-multiphase_ % M0*PI**2*multiphase_ % sigma*cos(PI*X(IX))*cos(PI*X(IZ))&
                      *sin(time)*1.2D1-multiphase_ % M0*PI**2*multiphase_ % sigma*cos(PI*X(IX))*cos(PI*X(IZ))**3*sin(time)*3.6D1&
                      -multiphase_ % M0*PI**2*multiphase_ % sigma*cos(PI*X(IX))**3*cos(PI*X(IZ))*sin(time)*3.6D1+multiphase_ % M0*PI**2*&
                      multiphase_ % sigma*cos(PI*X(IX))**3*cos(PI*X(IZ))**3*sin(time)*1.08D2-multiphase_ % M0*PI**2*multiphase_ % sigma*cos&
                      (PI*X(IX))**3*cos(PI*X(IZ))**3*cos(time)**2*sin(time)*1.08D2+multiphase_ % M0*PI**2*multiphase_ % sigma*&
                      cos(PI*X(IX))*cos(PI*X(IZ))**3*cos(time)**2*sin(time)*3.6D1+multiphase_ % M0*PI**2*multiphase_ % sigma*cos&
                      (PI*X(IX))**3*cos(PI*X(IZ))*cos(time)**2*sin(time)*3.6D1+multiphase_ % M0*multiphase_ % eps**2*PI**4*multiphase_ % sigma&
                      *cos(PI*X(IX))*cos(PI*X(IZ))*sin(time)*3.0D0)/multiphase_ % eps 

            S(IMSQRHOU) = cos(PI*X(IZ))*sin(PI*X(IX))*cos(time)+PI*cos(PI*X(IX))*sin(PI*X(IZ))*cos(time)&
                          *2.0D0+rho0*cos(PI*X(IZ))*sin(PI*X(IX))*cos(time)-PI*cos(PI*X(IZ))**3*sin(PI&
                          *X(IX))**3*sin(time)**3+PI*cos(PI*X(IX))*cos(PI*X(IZ))**2*sin(PI*X(IX))*sin(time)**2&
                          *2.0D0+cos(PI*X(IX))*cos(PI*X(IZ))**2*sin(PI*X(IX))*cos(time)*sin(time)*(3.0D0/2.0D0&
                          )+PI*cos(PI*X(IX))*sin(PI*X(IX))*sin(PI*X(IZ))**2*sin(time)**2*2.0D0+PI*&
                          rho0*cos(PI*X(IZ))**3*sin(PI*X(IX))**3*sin(time)**3+eta*PI**2*cos(PI*X(IZ))*sin&
                          (PI*X(IX))*sin(time)*4.0D0+PI*cos(PI*X(IX))**2*cos(PI*X(IZ))**3*sin(PI*X(IX))*&
                          sin(time)**3*2.0D0-rho0*cos(PI*X(IX))*cos(PI*X(IZ))**2*sin(PI*X(IX))*cos(time)*sin&
                          (time)*(3.0D0/2.0D0)+PI*rho0*cos(PI*X(IX))*sin(PI*X(IX))*sin(PI*X(IZ))**2*sin&
                          (time)**2*2.0D0+(PI*multiphase_ % sigma*cos(PI*X(IZ))*sin(PI*X(IX))*sin(time)*3.0D0)/multiphase_ % eps-multiphase_ % eps&
                          *PI**3*multiphase_ % sigma*cos(PI*X(IZ))*sin(PI*X(IX))*sin(time)*(3.0D0/4.0D0)-PI*rho0&
                          *cos(PI*X(IX))**2*cos(PI*X(IZ))**3*sin(PI*X(IX))*sin(time)**3*2.0D0+PI*cos(pi&
                          *X(IX))**2*cos(PI*X(IZ))*sin(PI*X(IX))*sin(PI*X(IZ))**2*sin(time)**3*3.0D0+PI*&
                          rho0*cos(PI*X(IX))*cos(PI*X(IZ))**2*sin(PI*X(IX))*sin(time)**2*2.0D0-PI*rho0*&
                          cos(PI*X(IX))**2*cos(PI*X(IZ))*sin(PI*X(IX))*sin(PI*X(IZ))**2*sin(time)**3*3.0D0-(&
                          PI*multiphase_ % sigma*cos(PI*X(IX))**2*cos(PI*X(IZ))**3*sin(PI*X(IX))*sin(time)**3*9.0D0)/&
                          multiphase_ % eps-(PI*multiphase_ % sigma*cos(PI*X(IX))**3*cos(PI*X(IZ))**4*sin(PI*X(IX))*sin(time)**4*&
                          9.0D0)/multiphase_ % eps+(PI*multiphase_ % sigma*cos(PI*X(IX))*cos(PI*X(IZ))**2*sin(PI*X(IX))*sin(time)**2*&
                          3.0D0)/multiphase_ % eps-multiphase_ % eps*PI**3*multiphase_ % sigma*cos(PI*X(IX))*cos(PI*X(IZ))**2*sin(PI*X(IX))*sin&
                          (time)**2*(3.0D0/4.0D0) 

            S(IMSQRHOV) = 0.0_RP 

            S(IMSQRHOW) = -cos(PI*X(IX))*sin(PI*X(IZ))*cos(time)+PI*cos(PI*X(IZ))*sin(PI*X(IX))*cos(time&
                           )*2.0D0-rho0*cos(PI*X(IX))*sin(PI*X(IZ))*cos(time)-PI*cos(PI*X(IX))**3*sin(PI&
                           *X(IZ))**3*sin(time)**3+PI*cos(PI*X(IX))**2*cos(PI*X(IZ))*sin(PI*X(IZ))*sin(time)**&
                           2*2.0D0-cos(PI*X(IX))**2*cos(PI*X(IZ))*sin(PI*X(IZ))*cos(time)*sin(time)*(3.0D0/2.0D0&
                           )+PI*cos(PI*X(IZ))*sin(PI*X(IX))**2*sin(PI*X(IZ))*sin(time)**2*2.0D0+PI*&
                           rho0*cos(PI*X(IX))**3*sin(PI*X(IZ))**3*sin(time)**3-eta*PI**2*cos(PI*X(IX))*&
                           sin(PI*X(IZ))*sin(time)*4.0D0+PI*cos(PI*X(IX))**3*cos(PI*X(IZ))**2*sin(PI*X(IZ))*&
                           sin(time)**3*2.0D0+rho0*cos(PI*X(IX))**2*cos(PI*X(IZ))*sin(PI*X(IZ))*cos(time)*&
                           sin(time)*(3.0D0/2.0D0)+PI*rho0*cos(PI*X(IZ))*sin(PI*X(IX))**2*sin(PI*X(IZ))*&
                           sin(time)**2*2.0D0+(PI*multiphase_ % sigma*cos(PI*X(IX))*sin(PI*X(IZ))*sin(time)*3.0D0)/multiphase_ % eps-&
                           multiphase_ % eps*PI**3*multiphase_ % sigma*cos(PI*X(IX))*sin(PI*X(IZ))*sin(time)*(3.0D0/4.0D0)-PI*&
                           rho0*cos(PI*X(IX))**3*cos(PI*X(IZ))**2*sin(PI*X(IZ))*sin(time)**3*2.0D0+PI*cos(&
                           PI*X(IX))*cos(PI*X(IZ))**2*sin(PI*X(IX))**2*sin(PI*X(IZ))*sin(time)**3*3.0D0+PI*&
                           rho0*cos(PI*X(IX))**2*cos(PI*X(IZ))*sin(PI*X(IZ))*sin(time)**2*2.0D0-PI*rho0*&
                           cos(PI*X(IX))*cos(PI*X(IZ))**2*sin(PI*X(IX))**2*sin(PI*X(IZ))*sin(time)**3*3.0D0-(&
                           PI*multiphase_ % sigma*cos(PI*X(IX))**3*cos(PI*X(IZ))**2*sin(PI*X(IZ))*sin(time)**3*9.0D0)/&
                           multiphase_ % eps-(PI*multiphase_ % sigma*cos(PI*X(IX))**4*cos(PI*X(IZ))**3*sin(PI*X(IZ))*sin(time)**4*&
                           9.0D0)/multiphase_ % eps+(PI*multiphase_ % sigma*cos(PI*X(IX))**2*cos(PI*X(IZ))*sin(PI*X(IZ))*sin(time)**2*&
                           3.0D0)/multiphase_ % eps-multiphase_ % eps*PI**3*multiphase_ % sigma*cos(PI*X(IX))**2*cos(PI*X(IZ))*sin(PI*X(IZ))*&
                           sin(time)**2*(3.0D0/4.0D0) 

            S(IMP) = sin(PI*X(IX))*sin(PI*x(IZ))*sin(time)*(-2.0_RP)
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
            real(kind=RP), parameter  :: saved_errors(5) = [1.6902480371161803E-06_RP, &
                                          1.9409444210673717E-05_RP, &
                                          7.4217194511428453E-17_RP, &
                                          1.9836864099170690E-05_RP, & 
                                          1.7108477257107393E-04_RP]
            integer        :: i, j,k, eID 
            CHARACTER(LEN=29)                  :: testName           = "Multiphase convergence"
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

                  phi = cos(PI*x)*cos(PI*z)*sin(time)
                  c = 0.5_RP*(1.0_RP + phi)

                  rho = c + rho0*(1.0_RP-c)

                  u = 2*sin(PI*x)*cos(PI*z)*sin(time)
                  v = 0.0_RP
                  w = -2*cos(PI*x)*sin(PI*z)*sin(time)
                  p = 2*sin(PI*x)*sin(PI*z)*cos(time)

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
      