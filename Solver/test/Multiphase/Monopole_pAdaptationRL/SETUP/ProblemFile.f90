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
!           Local variables
!           ---------------
!
            REAL(KIND=RP) :: x(3)        
            INTEGER       :: i, j, k, eID
            REAL(KIND=RP) :: rho , u , v , w , p
            REAL(KIND=RP) :: L, u_0, rho_0, p_0, theta, phi, eps,c
            integer       :: Nx, Ny, Nz
#if defined(NAVIERSTOKES)
            
            L     = 1.0_RP
            u_0   = 0.0_RP
            rho_0 = 1.0_RP 
            !p_0   = 0.714285714_RP !original test case definition
            p_0    = 1.0_RP/(dimensionless_%gammaM2)
            ! p_0   = rho_0 * u_0**2 / (dimensionless_ % gammaM2) 

            theta = refvalues_ % AOAtheta*(pi/180.0_RP)
            phi   = refvalues_ % AOAphi*(pi/180.0_RP)
            u  = u_0*cos(theta)*cos(phi)
            v  = u_0*sin(theta)*cos(phi)
            w  = u_0*sin(phi)
            p  = p_0

            associate( gamma => thermodynamics_ % gamma ) 
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx 

                         mesh % elements(eID) % storage % Q(1,i,j,k) = rho_0
                         mesh % elements(eID) % storage % Q(2,i,j,k) = rho_0*u
                         mesh % elements(eID) % storage % Q(3,i,j,k) = rho_0*v
                         mesh % elements(eID) % storage % Q(4,i,j,k) = rho_0*w
                         mesh % elements(eID) % storage % Q(5,i,j,k) = p / (gamma - 1.0_RP) + 0.5_RP * rho_0 * (u*u + v*v + w*w)

                     END DO
                  END DO
               END DO 
               
            END DO 
            end associate
#endif   

#if defined(INCNS)

            eps =0.25
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx 
                         
                         c  = 0.5*(1.0+tanh(2.0*(( mesh % elements(eID) % geom % x(IY,i,j,k) + 2.0))/eps))
                        

                         mesh % elements(eID) % storage % Q(1,i,j,k) = dimensionless_ % rho(1)*c + dimensionless_ % rho(2)*(1-c)  
                         mesh % elements(eID) % storage % Q(2,i,j,k) = 0.0
                         mesh % elements(eID) % storage % Q(3,i,j,k) = 0.0
                         mesh % elements(eID) % storage % Q(4,i,j,k) = 0.0
                         mesh % elements(eID) % storage % Q(5,i,j,k) = 0.0

                     END DO
                  END DO
               END DO 
               
            END DO 

#endif

#if defined(MULTIPHASE)

            
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx 
                         
                         c  = 0.5*(1.0+tanh(2.0*(( mesh % elements(eID) % geom % x(IY,i,j,k) + 2.0))/multiphase_ % eps))
                        
                         rho = dimensionless_ % rho(1)*c + dimensionless_ % rho(2)*(1.0_RP - c)
                         if (c<1e-5) then
                            c = 1e-5
                         endif

                         mesh % elements(eID) % storage % Q(1,i,j,k) = c 
                         mesh % elements(eID) % storage % Q(2,i,j,k) = 0.0
                         mesh % elements(eID) % storage % Q(3,i,j,k) = 0.0
                         mesh % elements(eID) % storage % Q(4,i,j,k) = 0.0
                         mesh % elements(eID) % storage % Q(5,i,j,k) = 0.0

                        !  if( mesh % elements(eID) % geom % x(IY,i,j,k)>0) then
                        !    mesh % elements(eID) % storage % Q(5,i,j,k) = 0.0
                        !  else
                        !    mesh % elements(eID) % storage % Q(5,i,j,k) = 1.0
                        !  endif

                     END DO
                  END DO
               END DO 
               
            END DO 

#endif

         end subroutine UserDefinedInitialCondition


         subroutine UserDefinedState1(x, t, nHat, Q &
#if defined(FLOW)
         ,thermodynamics_, dimensionless_, refValues_ &
#endif
#if defined(CAHNHILLIARD)
	,multiphase_ &
#endif
	)
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
#ifdef FLOW
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif	
	    real(kind=RP)     :: lambda, period, amplitude, depth, omega, wave_n
	    real(kind=RP)     :: pos, u_air,v_air,w_air, u_water, v_water, w_water  ,c, u, v, w, p,rho 
	         
         end subroutine UserDefinedState1

#if defined(FLOW)
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
            real(kind=RP), intent(inout)    :: U_x(NGRAD)
            real(kind=RP), intent(inout)    :: U_y(NGRAD)
            real(kind=RP), intent(inout)    :: U_z(NGRAD)
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
#endif
            real(kind=RP)                          :: f,c0,b,w
            real(kind=RP)                          :: x0(NDIM-1),r(NDIM-1)
            real(kind=RP)                          :: freqTerm
            integer, parameter                     :: Nfreq=0
            integer                                :: i
            real(kind=RP)                          :: fMax,fMin,df,freqVector(0:Nfreq)
            real(kind=RP)                          :: phi(0:Nfreq),xwrap(0:Nfreq),dummy(0:Nfreq)
            
#if defined(NAVIERSTOKES)
!
!           Usage example
!           -------------
            S = 0.0_RP
            b = 0.2_RP
            ! w = 2.0_RP*PI
            x0 = 0.0_RP
            ! fMax = 5.0_RP
            ! fMin = 0.5_RP
            fMin = 500.0_RP
            df = 0.5_RP

            c0 = 1.0 / dimensionless_ % Mach

            freqVector(0:Nfreq) = [(fMin+i*df,i=0, Nfreq)]
            ! phase using parabolic distribution to avoid over increases
            dummy = 1.0
            phi = [(i,i=1,Nfreq+1)]
            phi = 1 - phi * phi
            phi = -PI/(real(Nfreq,RP)+1)*phi
            xwrap = mod(phi, 2.0_RP*PI)
            phi = xwrap + merge(-2.0_RP*PI,0.0_RP,abs(xwrap)>PI)*sign(dummy, xwrap)

            ! s of p
            freqTerm = 0.0_RP
            do i = 0,Nfreq
                w = 2.0_RP*PI*freqVector(i)
                freqTerm = freqTerm + cos(w*time + phi(i))
            end do
            ! r = x-x0
            r = x(IX:IY)-x0
            ! f = 1.0_RP * exp(-log(2.0_RP)/(b*b)*sum(r*r) ) * cos(w*time)
            f = 1.0_RP * exp(-log(2.0_RP)/(b*b)*sum(r*r) ) * freqTerm

            !S(1) = f /(c0*c0)
            S(5) = f / (thermodynamics_ % gamma - 1.0_RP)
#endif
#if defined(INCNS)
!
!           Usage example
!           -------------
            S = 0.0_RP
            b = 0.2_RP
            ! w = 2.0_RP*PI
            x0 = 0.0_RP
            ! fMax = 5.0_RP
            ! fMin = 0.5_RP
            fMin = 500.0_RP
            df = 0.5_RP

            c0 = sqrt(thermodynamics_ % rho0c02 /  dimensionless_ % rho(1))

            freqVector(0:Nfreq) = [(fMin+i*df,i=0, Nfreq)]
            ! phase using parabolic distribution to avoid over increases
            dummy = 1.0
            phi = [(i,i=1,Nfreq+1)]
            phi = 1 - phi * phi
            phi = -PI/(real(Nfreq,RP)+1)*phi
            xwrap = mod(phi, 2.0_RP*PI)
            phi = xwrap + merge(-2.0_RP*PI,0.0_RP,abs(xwrap)>PI)*sign(dummy, xwrap)

            ! s of p
            freqTerm = 0.0_RP
            do i = 0,Nfreq
                w = 2.0_RP*PI*freqVector(i)
                freqTerm = freqTerm + cos(w*time + phi(i))
            end do
            ! r = x-x0
            r = x(IX:IY)-x0
            ! f = 1.0_RP * exp(-log(2.0_RP)/(b*b)*sum(r*r) ) * cos(w*time)
            f = 1.0_RP * exp(-log(2.0_RP)/(b*b)*sum(r*r) ) * freqTerm

            !S(1) = f /(c0*c0)
            S(5) = f 
#endif
#if defined(MULTIPHASE)
!
!           Usage example
!           -------------
            S = 0.0_RP
            b = 0.2_RP
            ! w = 2.0_RP*PI
            x0 = 0.0_RP
            ! fMax = 5.0_RP
            ! fMin = 0.5_RP
            fMin = 500.0_RP
            df = 0.5_RP

            c0 = sqrt(thermodynamics_ % c02(1)) 

            freqVector(0:Nfreq) = [(fMin+i*df,i=0, Nfreq)]
            ! phase using parabolic distribution to avoid over increases
            dummy = 1.0
            phi = [(i,i=1,Nfreq+1)]
            phi = 1 - phi * phi
            phi = -PI/(real(Nfreq,RP)+1)*phi
            xwrap = mod(phi, 2.0_RP*PI)
            phi = xwrap + merge(-2.0_RP*PI,0.0_RP,abs(xwrap)>PI)*sign(dummy, xwrap)

            ! s of p
            freqTerm = 0.0_RP
            do i = 0,Nfreq
                w = 2.0_RP*PI*freqVector(i)
                freqTerm = freqTerm + cos(w*time + phi(i))
            end do
            ! r = x-x0
            r = x(IX:IY)-x0
            ! f = 1.0_RP * exp(-log(2.0_RP)/(b*b)*sum(r*r) ) * cos(w*time)
            f = 1.0_RP * exp(-log(2.0_RP)/(b*b)*sum(r*r) ) * freqTerm

            !S(1) = f /(c0*c0)
            S(5) = f 
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
            use FTAssertions
            USE HexMeshClass
            use PhysicsStorage
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
!
!           ---------------
!           Local variables
!           ---------------
!
            CHARACTER(LEN=29)                  :: testName           = "Monopole P-adaptation RL"
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, N
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
!
!           -----------------------------------------------------------------------------------------
!           Expected solutions. 
!           InnerCylinder 0.0 NoSlipAdiabaticWall
!           Front 0.0 Inflow
!           bottom 0.0 FreeSlipWall
!           top 0.0 FreeSlipWall
!           Back 0.0 Inflow
!           Left 0.0 Inflow
!           Right 0.0 OutflowSpecifyP 
!           -----------------------------------------------------------------------------------------
!
!
!           -----------------------------------------------------------------
!           Number of iterations are for dt = 0.0000005, for the exact riemann solver,
!           air-water multiphase and a monopole with pressure-based RL p-adaptation each 100 iterations in the acoustic path
!           -----------------------------------------------------------------
!
#if defined(MULTIPHASE)
            real(kind=RP), parameter :: residuals(5) = [ 6.7612400583858490E+00_RP, &
                                                         1.7330535121898896E+00_RP, &
                                                         1.4330967844664411E+00_RP, &
                                                         4.5116531519710653E-15_RP, &
                                                         5.1783390542919660E+03_RP]

            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            CALL FTAssertEqual(expectedValue = residuals(1)+1.0_RP, &
                               actualValue   = monitors % residuals % values(1,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Continuity residual")

            CALL FTAssertEqual(expectedValue = residuals(2)+1.0_RP, &
                               actualValue   = monitors % residuals % values(2,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "X-Momentum residual")

            CALL FTAssertEqual(expectedValue = residuals(3)+1.0_RP, &
                               actualValue   = monitors % residuals % values(3,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Y-Momentum residual")

            CALL FTAssertEqual(expectedValue = residuals(4)+1.0_RP, &
                               actualValue   = monitors % residuals % values(4,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Z-Momentum residual")

            CALL FTAssertEqual(expectedValue = residuals(5)+1.0_RP, &
                               actualValue   = monitors % residuals % values(5,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Energy residual")


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
      


