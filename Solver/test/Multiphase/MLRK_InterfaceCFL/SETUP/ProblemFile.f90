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
            REAL(KIND=RP) :: L, u_0, rho_0, p_0, theta, phi, eps,c,angle
            integer       :: Nx, Ny, Nz


#if defined(MULTIPHASE)

            
angle = -10.0*PI/180.0 

DO eID = 1, SIZE(mesh % elements)
   Nx = mesh % elements(eID) % Nxyz(1)
   Ny = mesh % elements(eID) % Nxyz(2)
   Nz = mesh % elements(eID) % Nxyz(3)

   DO k = 0, Nz
      DO j = 0, Ny
         DO i = 0, Nx 
             
            
            c  = 1.0 - 0.5*(1.0+tanh(2.0*(( mesh % elements(eID) % geom % x(IX,i,j,k)*cos(angle) + mesh % elements(eID) % geom % x(IY,i,j,k)*sin(angle) + 0.0))/multiphase_ % eps))
             
                          
            if (c<1e-14) then
              c=1e-14
            endif
            

            mesh % elements(eID) % storage % Q(1,i,j,k) = c
!           Fast block well inside phase 1, far from the interface: the fastest elements are not interface elements
            if (mesh % elements(eID) % geom % x(IX,i,j,k) < -0.6_RP) then
               mesh % elements(eID) % storage % Q(2,i,j,k) = sqrt(dimensionless_ % rho(1)) * 1000.0_RP
            else
               mesh % elements(eID) % storage % Q(2,i,j,k) = 1e-14_RP
            end if
            mesh % elements(eID) % storage % Q(3,i,j,k) = 0.0
            mesh % elements(eID) % storage % Q(4,i,j,k) = 0.0
            mesh % elements(eID) % storage % Q(5,i,j,k) = 1e-14_RP  !+ exp(-1000.0*(mesh % elements(eID) % geom % x(IX,i,j,k)+0.70)**2.0) 

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
            integer, parameter                     :: Nfreq=1
            integer                                :: i
            real(kind=RP)                          :: fMax,fMin,df,freqVector(0:Nfreq)
            real(kind=RP)                          :: phi(0:Nfreq),xwrap(0:Nfreq),dummy(0:Nfreq)
            

#if defined(MULTIPHASE)
!           Usage example
!           -------------
            S = 0.0_RP
            b = 0.01_RP
            ! w = 2.0_RP*PI
            x0 = -0.55_RP
            ! fMax = 5.0_RP
            ! fMin = 0.5_RP
            !fMin = 500.0_RP
            !df = 0.5_RP

            !c0 = sqrt(thermodynamics_ % c02(1)) 

            !freqVector(0:Nfreq) = [(fMin+i*df,i=0, Nfreq)]
            !freqVector(0) = 1000.0_RP

            ! phase using parabolic distribution to avoid over increases
            ! dummy = 1.0
            ! phi = [(i,i=1,Nfreq+1)]
            ! phi = 1 - phi * phi
            ! phi = -PI/(real(Nfreq,RP)+1)*phi
            ! xwrap = mod(phi, 2.0_RP*PI)
            ! phi = xwrap + merge(-2.0_RP*PI,0.0_RP,abs(xwrap)>PI)*sign(dummy, xwrap)

            ! s of p
            ! freqTerm = 0.0_RP
            ! do i = 0,Nfreq
            !     w = 2.0_RP*PI*freqVector(i)
            !     freqTerm = freqTerm + cos(w*time + phi(i))
            ! end do
            ! r = x-x0
            r = x(IX:IY)-x0
            ! f = 1.0_RP * exp(-log(2.0_RP)/(b*b)*sum(r*r) ) * cos(w*time)
            ! f = 1.0_RP * exp(-log(2.0_RP)/(b*b)*sum(r*r) ) * freqTerm
            freqTerm = 1000.0_RP
            f = 1.0_RP * exp(-log(2.0_RP)/(b*b)*sum((x(IX)-x0)*(x(IX)-x0)) ) * cos(2.0_RP*PI*time*freqTerm)
            !S(1) = f /(c0*c0)
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
#if defined(MULTIPHASE)
!           Interface elements are assigned a CFL of 100 x (largest CFL among interface elements).
!           The case is tuned so that this lands in level 2, while 100 x (largest CFL of any element),
!           which is what the old interface test produced (true in all of phase 1), lands in level 3.
            real(kind=RP), parameter            :: TOL_EXPECTED = 1.0e-3_RP   ! "mlrk interface tolerance" in the control file
            real(kind=RP), parameter            :: DT_RATIO = 2.5_RP          ! level ratio used in MultiLevel_RK_Update
            character(len=*), parameter         :: testName = "MLRK multiphase interface CFL"
            type(FTAssertionsManager), pointer  :: sharedManager
            real(kind=RP)                       :: tol, maxc, minc, interfMax, globalMax
            integer                             :: eID, interfLevel, nInterf, nWrong
            real(kind=RP), allocatable          :: c0(:,:,:)

            tol       = mesh % MLRK % interfaceTol
            interfMax = 0.0_RP
            globalMax = 0.0_RP
            do eID = 1, size(mesh % elements)
               associate( e => mesh % elements(eID) )
               globalMax = max(globalMax, e % ML_CFL)
               c0 = initialC(e % geom % x)
               if ( any( (c0 >= tol) .and. (c0 <= 1.0_RP - tol) ) ) then
                  interfMax = max(interfMax, e % ML_CFL)
               end if
               end associate
            end do

            interfLevel = levelOf(100.0_RP * interfMax)
            nInterf = 0
            nWrong  = 0
            do eID = 1, size(mesh % elements)
               associate( e => mesh % elements(eID) )
               c0   = initialC(e % geom % x)
               maxc = min(maxval(c0), 1.0_RP)
               minc = max(minval(c0), 0.0_RP)
               if ( (maxc > tol) .and. (minc < 1.0_RP - tol) ) then
                  nInterf = nInterf + 1
                  if ( e % MLevel /= interfLevel ) nWrong = nWrong + 1
               end if
               end associate
            end do

            write(6,'(A,ES24.16)') "Interface tolerance:                          ", tol
            write(6,'(A,ES24.16)') "CFL cut-off:                                  ", mesh % MLRK % CFL_CutOff
            write(6,'(A,ES24.16,A,I2)') "100 x largest CFL among interface elements:   ", 100.0_RP * interfMax, &
                                        "  -> level", interfLevel
            write(6,'(A,ES24.16,A,I2)') "100 x largest CFL among all elements:         ", 100.0_RP * globalMax, &
                                        "  -> level", levelOf(100.0_RP * globalMax)
            write(6,'(A,I6,A,I6)') "Interface elements:", nInterf, ", in the wrong level:", nWrong

            call initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            call FTAssertEqual(expectedValue = TOL_EXPECTED, actualValue = tol, tol = 1.0e-15_RP, &
                               msg = "Interface tolerance read from the control file")
            call FTAssertEqual(expectedValue = 2, actualValue = interfLevel, &
                               msg = "Test setup: interface CFL x 100 falls in level 2")
            call FTAssertEqual(expectedValue = 3, actualValue = levelOf(100.0_RP * globalMax), &
                               msg = "Test setup: largest CFL x 100 falls in level 3")
            call FTAssert(test = nInterf > 0, msg = "Test setup: there are interface elements")
            call FTAssertEqual(expectedValue = 0, actualValue = nWrong, &
                               msg = "Interface elements in the level of the interface CFL")

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

         contains

            integer function levelOf(cfl)
               real(kind=RP), intent(in) :: cfl
               real(kind=RP)             :: cutoff
               integer                   :: i
               levelOf = mesh % MLRK % nLevel
               cutoff  = mesh % MLRK % CFL_CutOff
               do i = 1, mesh % MLRK % nLevel - 1
                  if (cfl < cutoff) then
                     levelOf = i
                     return
                  end if
                  cutoff = cutoff * DT_RATIO
               end do
            end function levelOf

!           Levels were assigned from the initial state, and the fast block changes c within one step,
!           so the interface is identified from the initial condition (same formula as UserDefinedInitialCondition)
            function initialC(x) result(c)
               real(kind=RP), intent(in)  :: x(:,0:,0:,0:)
               real(kind=RP)              :: c(0:size(x,2)-1, 0:size(x,3)-1, 0:size(x,4)-1)
               real(kind=RP), parameter   :: angle = -10.0_RP*PI/180.0_RP
               c = 1.0_RP - 0.5_RP*(1.0_RP+tanh(2.0_RP*(x(IX,:,:,:)*cos(angle) + x(IY,:,:,:)*sin(angle))/multiphase_ % eps))
               c = max(c, 1.0e-14_RP)
            end function initialC
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
      


