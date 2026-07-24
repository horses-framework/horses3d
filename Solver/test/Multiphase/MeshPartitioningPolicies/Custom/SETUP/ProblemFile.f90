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

module UserDefinedMeshPartitioning
   implicit none
   public
   contains
   
   subroutine mypartitioning(mesh, no_of_domains, elementsDomain, nodesDomain, controlVariables)
      use SMConstants
      use HexMeshClass
      use FTValueDictionaryClass
      implicit none
      type(HexMesh), intent(in)              :: mesh
      integer,       intent(in)              :: no_of_domains
      integer,       intent(out)             :: elementsDomain(mesh % no_of_elements)
      integer,       intent(out)             :: nodesDomain(size(mesh % nodes))
      type(FTValueDictionary), intent(in)    :: controlVariables

      integer :: ea, pe, dom
      real(kind=rp) :: centroid, centroid_d(3)

      ! Bounds in x direction are: [-1,2]
      do ea = 1, mesh % no_of_elements
         centroid_d = sum(mesh % elements(ea) % SurfInfo % corners, dim=2) / 8
         centroid = centroid_d(1)
         do dom = 1, no_of_domains
            if (centroid .lt. -1.0_rp + 3.0_rp / no_of_domains * real(dom, kind=rp)) then
               elementsDomain(ea) = dom
               exit
            end if
         end do
      end do
      ! Assign the nodes to each MPI domain
      do ea = 1, mesh % no_of_elements
         do pe = 1, 8
            nodesDomain(mesh % elements(ea) % nodeIDs(pe)) = elementsDomain(ea)
         end do
      end do
   end subroutine mypartitioning

end module UserDefinedMeshPartitioning



         SUBROUTINE UserDefinedStartup
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            use MeshPartitioningPolicies, only: setPointerGetRegion, setPointerUserDefinedMeshPartitioning
            use UserDefinedMeshPartitioning
            IMPLICIT NONE  

            call setPointerUserDefinedMeshPartitioning(mypartitioning)

         END SUBROUTINE UserDefinedStartup


!                   SUBROUTINE UserDefinedStartup
! !
! !        --------------------------------
! !        Called before any other routines
! !        --------------------------------
! !
! #ifdef MULTIPHASE
!             ! use, intrinsic :: iso_c_binding
!             use MeshPartitioningPolicies
!             use ajrmodule
!             IMPLICIT NONE

!             procedure(getRegion_f), pointer :: aaa
!             aaa => getRegion_fun
!             call setTimeStepCallback(aaa)
! #endif
!          END SUBROUTINE UserDefinedStartup
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
            mesh % elements(eID) % storage % Q(2,i,j,k) =  1e-14_RP !+ (1.0/sqrt(dimensionless_ % rho(1)))*(1.0/sqrt(thermodynamics_ % c02(1))) * exp(-1000.0*(mesh % elements(eID) % geom % x(IX,i,j,k)+0.70)**2.0) 
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
      


