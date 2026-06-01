#include "Includes.h"
module RiemannSolvers_NSKeywordsModule

     integer,                       parameter :: KEYWORD_LENGTH           = 132
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_SOLVER_NAME_KEY  = "riemann solver"
     character(len=KEYWORD_LENGTH), parameter :: LAMBDA_STABILIZATION_KEY = "lambda stabilization"
     character(len=KEYWORD_LENGTH), parameter :: AVG_NAME_KEY             = "averaging"
!
!    --------------------------
!    Riemann solver definitions
!    --------------------------
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_CENTRAL_NAME    = "central"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_ROE_NAME        = "roe"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_RUSANOV_NAME    = "rusanov"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_LXF_NAME        = "lax-friedrichs"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_STDROE_NAME     = "standard roe"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_ROEPIKE_NAME    = "roe-pike"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_LOWDISSROE_NAME = "low dissipation roe"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_MATRIXDISS_NAME = "matrix dissipation"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_VISCOUSNS_NAME  = "viscous ns"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_UDISS_NAME      = "u-diss"

     enum, bind(C)
        enumerator :: RIEMANN_ROE = 1, RIEMANN_LXF, RIEMANN_RUSANOV
        enumerator :: RIEMANN_STDROE, RIEMANN_CENTRAL, RIEMANN_ROEPIKE
        enumerator :: RIEMANN_LOWDISSROE, RIEMANN_VISCOUSNS, RIEMANN_MATRIXDISS
        enumerator :: RIEMANN_UDISS
     end enum
!
!    -----------------------------
!    Available averaging functions
!    -----------------------------
!
     character(len=KEYWORD_LENGTH), parameter :: STANDARD_AVG_NAME      = "standard"
     character(len=KEYWORD_LENGTH), parameter :: MORINISHI_AVG_NAME     = "morinishi"
     character(len=KEYWORD_LENGTH), parameter :: DUCROS_AVG_NAME        = "ducros"
     character(len=KEYWORD_LENGTH), parameter :: KENNEDYGRUBER_AVG_NAME = "kennedy-gruber"
     character(len=KEYWORD_LENGTH), parameter :: PIROZZOLI_AVG_NAME     = "pirozzoli"
     character(len=KEYWORD_LENGTH), parameter :: ENTROPYCONS_AVG_NAME   = "entropy conserving"
     character(len=KEYWORD_LENGTH), parameter :: CHANDRASEKAR_AVG_NAME  = "chandrasekar"

     enum, bind(C)
        enumerator :: STANDARD_AVG = 1, MORINISHI_AVG
        enumerator :: DUCROS_AVG, KENNEDYGRUBER_AVG
        enumerator :: PIROZZOLI_AVG, ENTROPYCONS_AVG
        enumerator :: CHANDRASEKAR_AVG
     end enum

end module RiemannSolvers_NSKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!
module RiemannSolvers_NS
   use SMConstants
   use Physics_NS
   use PhysicsStorage_NS
   use VariableConversion_NS
   use FluidData_NS

   implicit none

   private
   public whichAverage, whichRiemannSolver
   public SetRiemannSolver, DescribeRiemannSolver
   public RiemannSolver, AveragedStates, TwoPointFlux, RiemannSolver_dFdQ

   abstract interface
      subroutine RiemannSolverFCN(QLeft, QRight, nHat, t1, t2, flux)
         use SMConstants
         use PhysicsStorage_NS
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
      end subroutine RiemannSolverFCN
      subroutine AveragedStatesFCN(QLeft, QRight, pL, pR, invRhoL, invRhoR, flux)
         use SMConstants
         use PhysicsStorage_NS
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
      end subroutine AveragedStatesFCN
      subroutine TwoPointFluxFCN(QLeft, QRight, JaL, JaR, fSharp)
         use SMConstants
         use PhysicsStorage_NS
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
      end subroutine TwoPointFluxFCN
      subroutine RiemannSolver_dFdQFCN(ql, qr, nHat, dfdq_num, side)
         use SMConstants
         use PhysicsStorage_NS
         real(kind=RP), intent (in)  :: ql(NCONS)                 !<  Current solution on the left
         real(kind=RP), intent (in)  :: qr(NCONS)                 !<  Current solution on the right
         real(kind=RP), intent (in)  :: nHat(NDIM)                !<  Normal vector
         real(kind=RP), intent(out)  :: dfdq_num(NCONS,NCONS)     !>  Numerical flux Jacobian
         integer      , intent (in)  :: side
      end subroutine RiemannSolver_dFdQFCN
   end interface

   procedure(RiemannSolverFCN),      protected, pointer  :: RiemannSolver      => NULL()
   procedure(AveragedStatesFCN),     protected, pointer  :: AveragedStates     => NULL()
   procedure(TwoPointFluxFCN),       protected, pointer  :: TwoPointFlux       => NULL()
   procedure(RiemannSolver_dFdQFCN), protected, pointer  :: RiemannSolver_dFdQ => NULL()

   integer, protected :: whichRiemannSolver = -1
   integer, protected :: whichAverage = -1
   real(RP)           :: lambdaStab = 1.0_RP
!
!  ========
   contains
!  ========
!
      subroutine SetRiemannSolver(controlVariables)
!
!        -------
!        Modules
!        -------
         use Utilities, only: toLower
         use FTValueDictionaryClass
         use RiemannSolvers_NSKeywordsModule
!
!        ---------
!        Interface
!        ---------
         type(FTValueDictionary), intent(in) :: controlVariables
!
!        ---------------
!        Local variables
!        ---------------
         character(len=KEYWORD_LENGTH) :: keyword


         RiemannSolver_dFdQ => BaseClass_RiemannSolver_dFdQ
!
!        ---------------------------------------------
!        Choose the Riemann solver (by default is Roe)
!        ---------------------------------------------
         if (controlVariables % containsKey(RIEMANN_SOLVER_NAME_KEY)) then

            keyword = controlVariables % stringValueForKey(RIEMANN_SOLVER_NAME_KEY, KEYWORD_LENGTH)
            call toLower(keyword)

            select case (keyword)
            case (RIEMANN_ROE_NAME)
               RiemannSolver => RoeRiemannSolver
               whichRiemannSolver = RIEMANN_ROE

            case (RIEMANN_LXF_NAME)
               RiemannSolver => LxFRiemannSolver
               RiemannSolver_dFdQ => LxFRiemannSolver_dFdQ
               whichRiemannSolver = RIEMANN_LXF

            case (RIEMANN_RUSANOV_NAME)
              RiemannSolver => RusanovRiemannSolver
               whichRiemannSolver = RIEMANN_RUSANOV

            case (RIEMANN_STDROE_NAME)
              RiemannSolver => StdRoeRiemannSolver
               whichRiemannSolver = RIEMANN_STDROE

            case (RIEMANN_CENTRAL_NAME)
               RiemannSolver => CentralRiemannSolver
               whichRiemannSolver = RIEMANN_CENTRAL

            case (RIEMANN_ROEPIKE_NAME)
               RiemannSolver => RoePikeRiemannSolver
               whichRiemannSolver = RIEMANN_ROEPIKE

            case (RIEMANN_LOWDISSROE_NAME)
               RiemannSolver => LowDissipationRoeRiemannSolver
               whichRiemannSolver = RIEMANN_LOWDISSROE

            case (RIEMANN_MATRIXDISS_NAME)
               RiemannSolver => MatrixDissipationRiemannSolver
               whichRiemannSolver = RIEMANN_MATRIXDISS

            case (RIEMANN_VISCOUSNS_NAME)
               RiemannSolver => ViscousNSRiemannSolver
               whichRiemannSolver = RIEMANN_VISCOUSNS

            case (RIEMANN_UDISS_NAME)
              RiemannSolver => u_dissRiemannSolver
               whichRiemannSolver = RIEMANN_UDISS

            case default
               print*, "Riemann Solver not recognized."
               errorMessage(STD_OUT)
               error stop

            end select

         else
!
!           Select Roe by default
!           ---------------------
            RiemannSolver => RoeRiemannSolver
            whichRiemannSolver = RIEMANN_ROE

         end if
!
!        --------------------
!        Lambda stabilization
!        --------------------
         if (controlVariables % containsKey(LAMBDA_STABILIZATION_KEY)) then
            lambdaStab = controlVariables % doublePrecisionValueForKey(LAMBDA_STABILIZATION_KEY)

         else
!
!           By default, lambda is 1 (full upwind stabilization)
!           ---------------------------------------------------
            lambdaStab = 1.0_RP

         end if
!
!        If central fluxes are used, set lambdaStab to zero
!        --------------------------------------------------
         if (whichRiemannSolver .eq. RIEMANN_CENTRAL) lambdaStab = 0.0_RP
!
!        ----------------------------
!        Set up an averaging function
!        ----------------------------
         if (controlVariables % containsKey(AVG_NAME_KEY)) then

            keyword = controlVariables % stringValueForKey(AVG_NAME_KEY, KEYWORD_LENGTH)
            call toLower(keyword)

            select case (keyword)
            case (STANDARD_AVG_NAME)
               AveragedStates => StandardAverage
               TwoPointFlux => StandardDG_TwoPointFlux
               whichAverage = STANDARD_AVG

            case (DUCROS_AVG_NAME)
               AveragedStates => DucrosAverage
               TwoPointFlux => Ducros_TwoPointFlux
               whichAverage = DUCROS_AVG

            case (MORINISHI_AVG_NAME)
               AveragedStates => MorinishiAverage
               TwoPointFlux => Morinishi_TwoPointFlux
               whichAverage = MORINISHI_AVG

            case (PIROZZOLI_AVG_NAME)
               AveragedStates => PirozzoliAverage
               TwoPointFlux => Pirozzoli_TwoPointFlux
               whichAverage = PIROZZOLI_AVG

            case (KENNEDYGRUBER_AVG_NAME)
               AveragedStates => KennedyGruberAverage
               TwoPointFlux => KennedyGruber_TwoPointFlux
               whichAverage = KENNEDYGRUBER_AVG

            case (ENTROPYCONS_AVG_NAME)
               AveragedStates => EntropyConservingAverage
               TwoPointFlux => EntropyConserving_TwoPointFlux
               whichAverage = ENTROPYCONS_AVG

            case (CHANDRASEKAR_AVG_NAME)
               AveragedStates => ChandrasekarAverage
               TwoPointFlux => Chandrasekar_TwoPointFlux
               whichAverage = CHANDRASEKAR_AVG

            case default
               print*, "Averaging not recognized."
               errorMessage(STD_OUT)
               error stop

            end select

         else
!
!           Select standard by default
!           --------------------------
            AveragedStates => StandardAverage
            TwoPointFlux => StandardDG_TwoPointFlux
            whichAverage = STANDARD_AVG

         end if

      end subroutine SetRiemannSolver

      subroutine DescribeRiemannSolver
!
!        -------
!        Modules
!        -------
         use RiemannSolvers_NSKeywordsModule


         select case (whichAverage)
         case (STANDARD_AVG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Averaging function: ","Standard"

         case (MORINISHI_AVG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Averaging function: ","Morinishi"

         case (DUCROS_AVG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Averaging function: ","Ducros"

         case (KENNEDYGRUBER_AVG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Averaging function: ","Kennedy-Gruber"

         case (PIROZZOLI_AVG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Averaging function: ","Pirozzoli"

         case (ENTROPYCONS_AVG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Averaging function: ","Entropy conserving"

         case (CHANDRASEKAR_AVG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Averaging function: ","Chandrasekar"

         end select

         select case (whichRiemannSolver)
         case (RIEMANN_ROE)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Roe"

         case (RIEMANN_LXF)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Lax-Friedrichs"

         case (RIEMANN_RUSANOV)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Rusanov"

         case (RIEMANN_STDROE)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Standard Roe"

         case (RIEMANN_CENTRAL)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Central"

         case (RIEMANN_ROEPIKE)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Roe-Pike"

         case (RIEMANN_MATRIXDISS)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Matrix dissipation"

         case (RIEMANN_LOWDISSROE)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Low dissipation Roe"

         case (RIEMANN_UDISS)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","u-diss"

         end select

         write(STD_OUT,'(30X,A,A30,F10.3)') "->","Lambda stabilization: ", lambdaStab

      end subroutine DescribeRiemannSolver
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!        Riemann solvers
!        ---------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BaseClass_RiemannSolver_dFdQ(ql,qr,nHat,dfdq_num,side)
         implicit none
         !--------------------------------------------
         real(kind=RP), intent (in)  :: ql(NCONS)                 !<  Current solution on the left
         real(kind=RP), intent (in)  :: qr(NCONS)                 !<  Current solution on the right
         real(kind=RP), intent (in)  :: nHat(NDIM)                !<  Normal vector
         real(kind=RP), intent(out)  :: dfdq_num(NCONS,NCONS)     !>  Numerical flux Jacobian
         integer      , intent (in)  :: side                      !<  Either LEFT or RIGHT
         !--------------------------------------------

         error stop 'Requested Riemann solver not implemented for implicit time-integration'
      end subroutine BaseClass_RiemannSolver_dFdQ

      subroutine CentralRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: rhoL, rhouL, rhovL, rhowL, rhoeL, pL, rhoV2L
         real(kind=RP) :: rhoR, rhouR, rhovR, rhowR, rhoeR, pR, rhoV2R
         real(kind=RP) :: invRhoL, invRhoR
         real(kind=RP) :: QLRot(5), QRRot(5), oflux(NCONS)

         associate(gm1 => thermodynamics % GammaMinus1)
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(1)                  ; rhoR = QRight(1)
         invRhoL = 1.0_RP / rhoL          ; invRhoR = 1.0_RP / rhoR

         rhouL = QLeft(2) * nHat(1) + QLeft(3) * nHat(2) + QLeft(4) * nHat(3)
         rhovL = QLeft(2) * t1(1)   + QLeft(3) * t1(2)   + QLeft(4) * t1(3)
         rhowL = QLeft(2) * t2(1)   + QLeft(3) * t2(2)   + QLeft(4) * t2(3)

         rhouR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)
         rhovR = QRight(2) * t1(1)   + QRight(3) * t1(2)   + QRight(4) * t1(3)
         rhowR = QRight(2) * t2(1)   + QRight(3) * t2(2)   + QRight(4) * t2(3)

         rhoeL = QLeft(5) ; rhoeR = QRight(5)

         rhoV2L = (POW2(rhouL) + POW2(rhovL) + POW2(rhowL)) * invRhoL
         rhoV2R = (POW2(rhouR) + POW2(rhovR) + POW2(rhowR)) * invRhoR

         pL = gm1 * (rhoeL - 0.5_RP * rhoV2L)
         pR = gm1 * (rhoeR - 0.5_RP * rhoV2R)
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         QLRot = (/ rhoL, rhouL, rhovL, rhowL, rhoeL /)
         QRRot = (/ rhoR, rhouR, rhovR, rhowR, rhoeR /)
         call AveragedStates(QLRot, QRRot, pL, pR, invRhoL, invRhoR, flux)
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

         end associate

      end subroutine CentralRiemannSolver

      subroutine StdRoeRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
         use RiemannSolvers_NSKeywordsModule
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i
         real(kind=RP)  :: QLRot(5), QRRot(5), VL(NPRIM), VR(NPRIM), aL, aR
         real(kind=RP)  :: dQ(5), lambda(5), K(5,5), V2abs, alpha(5), dLambda
         real(kind=RP)  :: rho, u, v, w, V2, H, a
         real(kind=RP)  :: stab(5)     ! Careful with this variable

         associate(gm1 => thermodynamics % gammaMinus1)
!
!        ********************
!        Perform the rotation
!        ********************
!
         QLRot(1) = QLeft(1)  ; QRRot(1) = QRight(1)

         QLRot(2) = QLeft (2) * nHat(1) + QLeft (3) * nHat(2) + QLeft (4) * nHat(3)
         QRRot(2) = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

         QLRot(3) = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
         QRRot(3) = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

         QLRot(4) = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
         QRRot(4) = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

         QLRot(5) = QLeft(5) ; QRRot(5) = QRight(5)
!
!        ***************************
!        Compute primitive variables
!        ***************************
!
         call getPrimitiveVariables(QLRot, VL)
         call getPrimitiveVariables(QRRot, VR)

         aL = sqrt(VL(IPA2))  ; aR = sqrt(VR(IPA2))
!
!        *********************
!        Compute Roe variables: [rho, u, v, w, H, a]
!        *********************
!
         call getRoeVariables(QLRot, QRRot, VL, VR, rho, u, v, w, V2, H, a)
!
!        Eigenvalues
!        -----------
         lambda(1)   = u-a
         lambda(2:4) = u
         lambda(5)   = u+a
!
!        Eigenvectors
!        ------------
         K(:,1) = (/ 1.0_RP, u-a, v, w, H-u*a /)
         K(:,2) = (/ 1.0_RP, u, v, w, 0.5_RP*V2 /)
         K(:,3) = (/ 0.0_RP, 0.0_RP, 1.0_RP, 0.0_RP, v /)
         K(:,4) = (/ 0.0_RP, 0.0_RP, 0.0_RP, 1.0_RP, w /)
         K(:,5) = (/ 1.0_RP, u+a, v, w, H+u*a /)
!
!        Projections
!        -----------
         dQ = QRRot - QLRot

         alpha(3) = dQ(3) - v * dQ(1)  ; alpha(4) = dQ(4) - w * dQ(1)

         dQ(5) = dQ(5) - alpha(3) * v - alpha(4) * w

         alpha(2) = gm1 * ( dQ(1)*(H - u*u) + u * dQ(2) - dQ(5) ) / (POW2(a))
         alpha(1) = 0.5_RP * (dQ(1)*lambda(5) - dQ(2) - a*alpha(2)) / a
         alpha(5) = dQ(1) - alpha(1) - alpha(2)
!
!        **********************
!        Perform an entropy fix. Here we use Van Leer's modification of Harten's entropy fix, derived
!        in: A. Harten, "High resolution schemes for hyperbolic conservation laws". To recover the
!        Harten entropy fix, set dLambda to 0.5
!        **********************
!
!        Wave #1
!        -------
         dLambda = max((VR(IPU)-aR) - (VL(IPU)-aL), 0.0_RP)
         if ( abs(lambda(1)) .ge. 2.0_RP * dLambda ) then
            lambda(1) = abs(lambda(1))

         else
            lambda(1) = POW2(lambda(1)) / (4.0_RP * dLambda) + dLambda

         end if
!
!        Wave #5
!        -------
         dLambda = max((VR(IPU)+aR) - (VL(IPU)+aL), 0.0_RP)
         if ( abs(lambda(5)) .ge. 2.0_RP * dLambda ) then
            lambda(5) = abs(lambda(5))

         else
            lambda(5) = POW2(lambda(5)) / (4.0_RP * dLambda) + dLambda

         end if

!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         call AveragedStates(QLRot, QRRot, VL(IPP), VR(IPP), VL(IPIRHO), VR(IPIRHO), flux)
!
!        Compute the Roe stabilization
!        -----------------------------
         select case (whichAverage)
         case(PIROZZOLI_AVG, KENNEDYGRUBER_AVG)
!
!           ***************************************************************************
!           Eigenvalue matrix is corrected for PI and KG variants, see Winters et. al.
!           "A comparative study on polynomial dealiasing and split form discontinuous
!           Galerkin schemes for under-resolved turbulence computations"
!           ***************************************************************************
!
            lambda(1) = lambda(5)
         end select

         stab = 0.0_RP
         do i = 1, 5
            stab = stab + 0.5_RP * alpha(i) * abs(lambda(i)) * K(:,i)
         end do
!
!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux = flux - lambdaStab * stab
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

         end associate

      end subroutine StdRoeRiemannSolver

      subroutine MatrixDissipationRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
         use RiemannSolvers_NSKeywordsModule
         use Utilities, only: logarithmicMean
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, k
         real(kind=RP)  :: QLRot(5), QRRot(5), betaL, betaR
         real(kind=RP)  :: EVL(5), EVR(5)
         real(kind=RP)  :: dQ(5)
         real(kind=RP)  :: a_bar, h_bar
         real(kind=RP)  :: R1(5,5), T(5), Lambda(5)
         real(kind=RP)  :: stab(5)
         real(kind=RP)  :: rhoLogMean, betaLogMean, pMean, uMean, vMean, wMean, V2abs
         real(kind=RP)  :: uL, vL, wL, uR, vR, wR, vtotL, vtotR, pL, pR
         real(kind=RP)  :: invRhoL, invRhoR

         associate(gm1 => thermodynamics % gammaMinus1, gamma => thermodynamics % gamma, &
                   invGamma => thermodynamics % invGamma, cp => thermodynamics % GammaDivGammaMinus1, &
                   gammaMinus1Div2g => thermodynamics % gammaMinus1Div2g)
!
!        ********************
!        Perform the rotation
!        ********************
!
         QLRot(1) = QLeft(1)  ; QRRot(1) = QRight(1)

         QLRot(2) = QLeft (2) * nHat(1) + QLeft (3) * nHat(2) + QLeft (4) * nHat(3)
         QRRot(2) = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

         QLRot(3) = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
         QRRot(3) = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

         QLRot(4) = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
         QRRot(4) = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

         QLRot(5) = QLeft(5) ; QRRot(5) = QRight(5)
!
!        *************************
!        Compute Entropy variables
!        *************************
!
         call NSGradientVariables_ENTROPY(NCONS,NGRAD,QLRot,EVL)
         call NSGradientVariables_ENTROPY(NCONS,NGRAD,QRRot,EVR)

         invRhoL = 1.0_RP / QLRot(IRHO) ; invrhoR = 1.0_RP / QRRot(IRHO)
         uL      = QLRot(IRHOU)*invRhoL ; uR      = QRRot(IRHOU)*invRhoR
         vL      = QLRot(IRHOV)*invRhoL ; vR      = QRRot(IRHOV)*invRhoR
         wL      = QLRot(IRHOW)*invRhoL ; wR      = QRRot(IRHOW)*invRhoR
         vtotL   = uL*uL + vL*vL + wL*wL; vtotR   = uR*uR + vR*vR + wR*wR

         pL      = gm1*(QLRot(IRHOE)-0.5_RP*QLRot(IRHO)*vtotL)
         pR      = gm1*(QRRot(IRHOE)-0.5_RP*QRRot(IRHO)*vtotR)

         betaL   = -0.5_RP*EVL(IRHOE) ; betaR = -0.5_RP*EVR(IRHOE)

         call logarithmicMean(betaL      , betaR      , betaLogMean)
         call logarithmicMean(QLRot(IRHO), QRRot(IRHO), rhoLogMean)

         pMean = 0.5_RP*(QLRot(IRHO)+QRRot(IRHO))/(betaL + betaR)
         a_bar = sqrt(gamma * pMean / rhoLogMean)
         uMean = AVERAGE(uL, uR)
         vMean = AVERAGE(vL, vR)
         wMean = AVERAGE(wL, wR)
         v2Abs =   2.0_RP * ( POW2(uMean) + POW2(vMean) + POW2(wMean) )      &
                 - 0.5_RP * ( vtotL + vtotR )

         h_bar = 0.5_RP * ( cp / betaLogMean + v2Abs )
!
!        ***********************
!        Compute the eigenvalues
!        ***********************
!
         lambda(1)   = abs(uMean - a_bar)
         lambda(2:4) = abs(uMean        )
         lambda(5)   = abs(uMean + a_bar)
!
!        Eigenvectors
!        ------------
         R1(:, 1) = (/ 1.0_RP, uMean-a_bar, vMean , wMean , h_bar-uMean*a_bar/)
         R1(:, 2) = (/ 1.0_RP, uMean      , vMean , wMean , 0.5_RP*V2abs     /)
         R1(:, 3) = (/ 0.0_RP, 0.0_RP     , 1.0_RP, 0.0_RP, vMean            /)
         R1(:, 4) = (/ 0.0_RP, 0.0_RP     , 0.0_RP, 1.0_RP, wMean            /)
         R1(:, 5) = (/ 1.0_RP, uMean+a_bar, vMean , wMean , h_bar+uMean*a_bar/)
!
!        Intensities
!        -----------
         T(1) = 0.5_RP * rhoLogMean * invGamma
         T(2) = 2.0_RP * gammaMinus1Div2g * rhoLogMean
         T(3) = pMean
         T(4) = pMean
         T(5) = T(1)
!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         call AveragedStates(QLRot, QRRot, pL, pR, invRhoL, invRhoR, flux)
!
!        Compute the Roe stabilization
!        -----------------------------
         select case (whichAverage)
         case(PIROZZOLI_AVG, KENNEDYGRUBER_AVG)
!
!           ***************************************************************************
!           Eigenvalue matrix is corrected for PI and KG variants, see Winters et. al.
!           "A comparative study on polynomial dealiasing and split form discontinuous
!           Galerkin schemes for under-resolved turbulence computations"
!           ***************************************************************************
!
            lambda(1) = lambda(5)
         end select

         stab = 0.0_RP
         do i = 1, 5
            do j = 1, 5 ; do k = 1, 5
               stab(i) = stab(i) + 0.5_RP * R1(i,j) * lambda(j) * T(j) * R1(k,j) * (EVR(k) - EVL(k))
            end do      ; end do
         end do
!
!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux = flux - lambdaStab * stab
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

         end associate

      end subroutine MatrixDissipationRiemannSolver

      subroutine RoePikeRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
         use RiemannSolvers_NSKeywordsModule
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i
         real(kind=RP)  :: QLRot(5), QRRot(5), VL(NPRIM), VR(NPRIM), aL, aR
         real(kind=RP)  :: dQ(5), lambda(5), K(5,5), V2abs, alpha(5), dLambda
         real(kind=RP)  :: rho, u, v, w, V2, H, a
         real(kind=RP)  :: stab(5)

         associate(gm1 => thermodynamics % gammaMinus1)
!
!        ********************
!        Perform the rotation
!        ********************
!
         QLRot(1) = QLeft(1)  ; QRRot(1) = QRight(1)

         QLRot(2) = QLeft (2) * nHat(1) + QLeft (3) * nHat(2) + QLeft (4) * nHat(3)
         QRRot(2) = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

         QLRot(3) = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
         QRRot(3) = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

         QLRot(4) = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
         QRRot(4) = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

         QLRot(5) = QLeft(5) ; QRRot(5) = QRight(5)
!
!        ***************************
!        Compute primitive variables
!        ***************************
!
         call getPrimitiveVariables(QLRot, VL)
         call getPrimitiveVariables(QRRot, VR)

         aL = sqrt(VL(IPA2))  ; aR = sqrt(VR(IPA2))
!
!        *********************
!        Compute Roe variables: [rho, u, v, w, H, a]
!        *********************
!
         call getRoeVariables(QLRot, QRRot, VL, VR, rho, u, v, w, V2, H, a)
!
!        Eigenvalues
!        -----------
         lambda(1)   = u-a
         lambda(2:4) = u
         lambda(5)   = u+a
!
!        Eigenvectors
!        ------------
         K(:,1) = (/ 1.0_RP, u-a, v, w, H-u*a /)
         K(:,2) = (/ 1.0_RP, u, v, w, 0.5_RP*V2 /)
         K(:,3) = (/ 0.0_RP, 0.0_RP, 1.0_RP, 0.0_RP, v /)
         K(:,4) = (/ 0.0_RP, 0.0_RP, 0.0_RP, 1.0_RP, w /)
         K(:,5) = (/ 1.0_RP, u+a, v, w, H+u*a /)
!
!        Projections
!        -----------
         alpha(1) = ((VR(IPP)-VL(IPP)) - rho * a * (VR(IPU)-VL(IPU)))/(2.0_RP * a * a)
         alpha(2) = (QRight(IRHO)-QLeft(IRHO)) - (VR(IPP)-VL(IPP))/(a*a)
         alpha(3) = rho * (VR(IPV)-VL(IPV))
         alpha(4) = rho * (VR(IPW)-VL(IPW))
         alpha(5) = ((VR(IPP)-VL(IPP)) + rho * a * (VR(IPU)-VL(IPU)))/(2.0_RP * a * a)
!
!        **********************
!        Perform an entropy fix. Here we use Van Leer's modification of Harten's entropy fix, derived
!        in: A. Harten, "High resolution schemes for hyperbolic conservation laws". To recover the
!        Harten entropy fix, set dLambda to 0.5
!        **********************
!
!        Wave #1
!        -------
         dLambda = max((VR(IPU)-aR) - (VL(IPU)-aL), 0.0_RP)
         if ( abs(lambda(1)) .ge. 2.0_RP * dLambda ) then
            lambda(1) = abs(lambda(1))

         else
            lambda(1) = POW2(lambda(1)) / (4.0_RP * dLambda) + dLambda

         end if
!
!        Wave #5
!        -------
         dLambda = max((VR(IPU)+aR) - (VL(IPU)+aL), 0.0_RP)
         if ( abs(lambda(5)) .ge. 2.0_RP * dLambda ) then
            lambda(5) = abs(lambda(5))

         else
            lambda(5) = POW2(lambda(5)) / (4.0_RP * dLambda) + dLambda

         end if
!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         call AveragedStates(QLRot, QRRot, VL(IPP), VR(IPP), VL(IPIRHO), VR(IPIRHO), flux)
!
!        Compute the Roe stabilization
!        -----------------------------
         select case (whichAverage)
         case(PIROZZOLI_AVG, KENNEDYGRUBER_AVG)
!
!           ***************************************************************************
!           Eigenvalue matrix is corrected for PI and KG variants, see Winters et. al.
!           "A comparative study on polynomial dealiasing and split form discontinuous
!           Galerkin schemes for under-resolved turbulence computations"
!           ***************************************************************************
!
            lambda(1) = lambda(5)
         end select

         stab = 0.0_RP
         do i = 1, 5
            stab = stab + 0.5_RP * alpha(i) * abs(lambda(i)) * K(:,i)
         end do
!
!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux = flux - lambdaStab * stab
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

         end associate
      end subroutine RoePikeRiemannSolver

      subroutine LowDissipationRoeRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
!
!        ***********************************************************************
!           This version, presented by Oßwald et. al. [*], is a modification of
!           Roe-Pike solver, decreasing the velocity jumps intensity for low
!           Mach numbers. Normal velocities are scaled such that
!                         du <- z*du
!           where z tends to zero as the Mach number tends to zero. Precisely:
!                         z = min(1, max(ML, MR))
!
!           These normal velocities are just scaled to compute Roe dissipation's
!           intensities (alpha coefficients), not the fluxes (as in the original
!           Low dissipation method by Thornber et. al.)
!
!           * Oßwald et. al. L2Roe: a low dissipation version of Roe’s a
!              pproximate Riemann solver for low Mach numbers
!        ***********************************************************************
!
         use RiemannSolvers_NSKeywordsModule
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i
         real(kind=RP)  :: rhoL, rhouL, rhovL, rhowL, rhoeL, pL, rhoHL, rhoV2L, ML
         real(kind=RP)  :: rhoR, rhouR, rhovR, rhowR, rhoeR, pR, rhoHR, rhoV2R, MR
         real(kind=RP)  :: uL, vL, wL, uR, vR, wR, aL, aR, dLambda, z, du, dv, dw
         real(kind=RP)  :: dp
         real(kind=RP)  :: QLRot(5), QRRot(5)
         real(kind=RP)  :: sqrtRhoL, sqrtRhoR, invSumSqrtRhoLR
         real(kind=RP)  :: invSqrtRhoL, invSqrtRhoR, invRhoL, invRhoR
         real(kind=RP)  :: rho, u, v, w, H, a, dQ(5), lambda(5), K(5,5), V2abs, alpha(5)
         real(kind=RP)  :: stab(5)

         associate(gamma => thermodynamics % gamma, gm1 => thermodynamics % gammaMinus1)
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(1)                  ; rhoR = QRight(1)
         invRhoL = 1.0_RP/ rhoL           ; invRhoR = 1.0_RP / rhoR
         sqrtRhoL = sqrt(rhoL)            ; sqrtRhoR = sqrt(rhoR)
         invSqrtRhoL = 1.0_RP / sqrtRhoL  ; invSqrtRhoR = 1.0_RP / sqrtRhoR
         invSumSqrtRhoLR = 1.0_RP / (sqrtRhoL + sqrtRhoR)

         rhouL = QLeft (2) * nHat(1) + QLeft (3) * nHat(2) + QLeft (4) * nHat(3)
         rhouR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

         rhovL = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
         rhovR = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

         rhowL = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
         rhowR = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

         rhoeL = QLeft(5) ; rhoeR = QRight(5)

         uL = rhouL * invRhoL    ; uR = rhouR * invRhoR
         vL = rhovL * invRhoL    ; vR = rhovR * invRhoR
         wL = rhowL * invRhoL    ; wR = rhowR * invRhoR

         rhoV2L = (POW2(uL) + POW2(vL) + POW2(wL)) * rhoL
         rhoV2R = (POW2(uR) + POW2(vR) + POW2(wR)) * rhoR
!
!        Compute the enthalpy: here defined as rhoH = gogm1 p + 0.5 rho V^2
!        --------------------
         rhoHL = gamma*rhoeL - 0.5_RP*gm1*rhoV2L
         rhoHR = gamma*rhoeR - 0.5_RP*gm1*rhoV2R

         pL = gm1 * (rhoeL - 0.5_RP * rhoV2L)
         pR = gm1 * (rhoeR - 0.5_RP * rhoV2R)

         aL = sqrt(gamma * pL * invRhoL)
         aR = sqrt(gamma * pR * invRhoR)
!
!        Compute Roe - Pike variables
!        ----------------------------
         rho = sqrtRhoL * sqrtRhoR
         u = (invSqrtRhoL * rhouL + invSqrtRhoR * rhouR) * invSumSqrtRhoLR
         v = (invSqrtRhoL * rhovL + invSqrtRhoR * rhovR) * invSumSqrtRhoLR
         w = (invSqrtRhoL * rhowL + invSqrtRhoR * rhowR) * invSumSqrtRhoLR
         H = (invSqrtRhoL * rhoHL + invSqrtRhoR * rhoHR) * invSumSqrtRhoLR
         V2abs = POW2(u) + POW2(v) + POW2(w)
         a = sqrt(gm1*(H - 0.5_RP*V2abs))
!
!        Eigenvalues
!        -----------
         lambda(1)   = u-a
         lambda(2:4) = u
         lambda(5)   = u+a
!
!        Eigenvectors
!        ------------
         K(:,1) = (/ 1.0_RP, u-a, v, w, H-u*a /)
         K(:,2) = (/ 1.0_RP, u, v, w, 0.5_RP*V2abs /)
         K(:,3) = (/ 0.0_RP, 0.0_RP, 1.0_RP, 0.0_RP, v /)
         K(:,4) = (/ 0.0_RP, 0.0_RP, 0.0_RP, 1.0_RP, w /)
         K(:,5) = (/ 1.0_RP, u+a, v, w, H+u*a /)
!
!        Projections
!        -----------
!
!        ----------------------------------------------------------------------------
!        Low dissipation Roe-Pike Riemann solver: Reduce the dissipation associated
!        to the jump in normal velocity. See Obwald et. al. L2Roe: a low dissipation
!        version of Roe’s approximate Riemann solver for low Mach numbers
!        ----------------------------------------------------------------------------

         ML = abs(uL) / aL    ; MR = abs(uR) / aR
         z  = min(1.0_RP, max(ML,MR))

         du = z * (uR - uL)
         dv = z * (vR - vL)
         dw = z * (wR - wL)
         dp = pR - pL

         alpha(1) = (dp - rho * a * du)/(2.0_RP * a * a)
         alpha(2) = (rhoR-rhoL) - dp/(a*a)
         alpha(3) = rho * dv
         alpha(4) = rho * dw
         alpha(5) = (dp + rho * a * du)/(2.0_RP * a * a)
!
!        **********************
!        Perform an entropy fix. Here we use Van Leer's modification of Harten's entropy fix, derived
!        in: A. Harten, "High resolution schemes for hyperbolic conservation laws". To recover the
!        Harten entropy fix, set dLambda to 0.5
!        **********************
!
!        Wave #1
!        -------
         dLambda = max((uR-aR) - (uL-aL), 0.0_RP)
         if ( abs(lambda(1)) .ge. 2.0_RP * dLambda ) then
            lambda(1) = abs(lambda(1))

         else
            lambda(1) = POW2(lambda(1)) / (4.0_RP * dLambda) + dLambda

         end if
!
!        Wave #5
!        -------
         dLambda = max((uR+aR) - (uL+aL), 0.0_RP)
         if ( abs(lambda(5)) .ge. 2.0_RP * dLambda ) then
            lambda(5) = abs(lambda(5))

         else
            lambda(5) = POW2(lambda(5)) / (4.0_RP * dLambda) + dLambda

         end if
!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         QLRot = (/ rhoL, rhouL, rhovL, rhowL, rhoeL /)
         QRRot = (/ rhoR, rhouR, rhovR, rhowR, rhoeR /)
         call AveragedStates(QLRot, QRRot, pL, pR, invRhoL, invRhoR, flux)
!
!        Compute the Roe stabilization
!        -----------------------------
         select case (whichAverage)
         case(PIROZZOLI_AVG, KENNEDYGRUBER_AVG)
!
!           ***************************************************************************
!           Eigenvalue matrix is corrected for PI and KG variants, see Winters et. al.
!           "A comparative study on polynomial dealiasing and split form discontinuous
!           Galerkin schemes for under-resolved turbulence computations"
!           ***************************************************************************
!
            lambda(1) = lambda(5)
         end select

         stab = 0.0_RP
         do i = 1, 5
            stab = stab + 0.5_RP * alpha(i) * abs(lambda(i)) * K(:,i)
         end do
!
!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux = flux - lambdaStab * stab
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

         end associate

      end subroutine LowDissipationRoeRiemannSolver

      subroutine ViscousNSRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
!
!        ***********************************************************************
!           This solver is designed such that it mimics the dissipation
!           introduced by compressible Navier-Stokes terms
!        ***********************************************************************
!
         use RiemannSolvers_NSKeywordsModule
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i
         real(kind=RP)  :: rhoL, rhouL, rhovL, rhowL, rhoeL, pL, rhoHL, rhoV2L, ML, TL
         real(kind=RP)  :: rhoR, rhouR, rhovR, rhowR, rhoeR, pR, rhoHR, rhoV2R, MR, TR
         real(kind=RP)  :: uL, vL, wL, uR, vR, wR, aL, aR, dLambda, z, du, dv, dw
         real(kind=RP)  :: dp
         real(kind=RP)  :: QLRot(5), QRRot(5)
         real(kind=RP)  :: sqrtRhoL, sqrtRhoR, invSumSqrtRhoLR
         real(kind=RP)  :: invSqrtRhoL, invSqrtRhoR, invRhoL, invRhoR
         real(kind=RP)  :: rho, u, v, w, H, a, lambda(5), V2abs, tau(5,5)
         real(kind=RP)  :: stab(5), kappa
         real(kind=RP)  :: divV

         associate(gamma => thermodynamics % gamma, gm1 => thermodynamics % gammaMinus1)
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(1)                  ; rhoR = QRight(1)
         invRhoL = 1.0_RP/ rhoL           ; invRhoR = 1.0_RP / rhoR
         sqrtRhoL = sqrt(rhoL)            ; sqrtRhoR = sqrt(rhoR)
         invSqrtRhoL = 1.0_RP / sqrtRhoL  ; invSqrtRhoR = 1.0_RP / sqrtRhoR
         invSumSqrtRhoLR = 1.0_RP / (sqrtRhoL + sqrtRhoR)

         rhouL = QLeft (2) * nHat(1) + QLeft (3) * nHat(2) + QLeft (4) * nHat(3)
         rhouR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

         rhovL = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
         rhovR = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

         rhowL = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
         rhowR = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

         rhoeL = QLeft(5) ; rhoeR = QRight(5)

         uL = rhouL * invRhoL    ; uR = rhouR * invRhoR
         vL = rhovL * invRhoL    ; vR = rhovR * invRhoR
         wL = rhowL * invRhoL    ; wR = rhowR * invRhoR

         rhoV2L = (POW2(uL) + POW2(vL) + POW2(wL)) * rhoL
         rhoV2R = (POW2(uR) + POW2(vR) + POW2(wR)) * rhoR
!
!        Compute the enthalpy: here defined as rhoH = gogm1 p + 0.5 rho V^2
!        --------------------
         rhoHL = gamma*rhoeL - 0.5_RP*gm1*rhoV2L
         rhoHR = gamma*rhoeR - 0.5_RP*gm1*rhoV2R

         pL = gm1 * (rhoeL - 0.5_RP * rhoV2L)
         pR = gm1 * (rhoeR - 0.5_RP * rhoV2R)

         aL = sqrt(gamma * pL * invRhoL)
         aR = sqrt(gamma * pR * invRhoR)
!
!        Compute Roe - Pike variables
!        ----------------------------
         rho = sqrtRhoL * sqrtRhoR
         u = (invSqrtRhoL * rhouL + invSqrtRhoR * rhouR) * invSumSqrtRhoLR
         v = (invSqrtRhoL * rhovL + invSqrtRhoR * rhovR) * invSumSqrtRhoLR
         w = (invSqrtRhoL * rhowL + invSqrtRhoR * rhowR) * invSumSqrtRhoLR
         H = (invSqrtRhoL * rhoHL + invSqrtRhoR * rhoHR) * invSumSqrtRhoLR
         V2abs = POW2(u) + POW2(v) + POW2(w)
         a = sqrt(gm1*(H - 0.5_RP*V2abs))
!
!        Eigenvalues
!        -----------
         lambda(1)   = u-a
         lambda(2:4) = u
         lambda(5)   = u+a
!
!        **********************
!        Perform an entropy fix. Here we use Van Leer's modification of Harten's entropy fix, derived
!        in: A. Harten, "High resolution schemes for hyperbolic conservation laws". To recover the
!        Harten entropy fix, set dLambda to 0.5
!        **********************
!
!        Wave #1
!        -------
         dLambda = max((uR-aR) - (uL-aL), 0.0_RP)
         if ( abs(lambda(1)) .ge. 2.0_RP * dLambda ) then
            lambda(1) = abs(lambda(1))

         else
            lambda(1) = POW2(lambda(1)) / (4.0_RP * dLambda) + dLambda

         end if
!
!        Wave #5
!        -------
         dLambda = max((uR+aR) - (uL+aL), 0.0_RP)
         if ( abs(lambda(5)) .ge. 2.0_RP * dLambda ) then
            lambda(5) = abs(lambda(5))

         else
            lambda(5) = POW2(lambda(5)) / (4.0_RP * dLambda) + dLambda

         end if
!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         QLRot = (/ rhoL, rhouL, rhovL, rhowL, rhoeL /)
         QRRot = (/ rhoR, rhouR, rhovR, rhowR, rhoeR /)
         call AveragedStates(QLRot, QRRot, pL, pR, invRhoL, invRhoR, flux)
!
!        Compute the stabilization
!        -------------------------
         select case (whichAverage)
         case(PIROZZOLI_AVG, KENNEDYGRUBER_AVG)
!
!           ***************************************************************************
!           Eigenvalue matrix is corrected for PI and KG variants, see Winters et. al.
!           "A comparative study on polynomial dealiasing and split form discontinuous
!           Galerkin schemes for under-resolved turbulence computations"
!           ***************************************************************************
!
            lambda(1) = lambda(5)
         end select
!
!        Use the non-rotated velocities
!        ------------------------------
         uL = QLeft(2) / QLeft(1)   ;     uR = QRight(2) / QRight(1)
         vL = QLeft(3) / QLeft(1)   ;     vR = QRight(3) / QRight(1)
         wL = QLeft(4) / QLeft(1)   ;     wR = QRight(4) / QRight(1)

         u = 0.5_RP * (uL + uR)
         v = 0.5_RP * (vL + vR)
         w = 0.5_RP * (wL + wR)

         divV = (uR-uL)*nHat(1) + (vR-vL)*nHat(2) + (wR-wL)*nHat(3)

         tau(1,1) = 2.0_RP * (uR-uL)*nHat(1) - (1.0_RP/3.0_RP) * divV
         tau(2,2) = 2.0_RP * (vR-vL)*nHat(2) - (1.0_RP/3.0_RP) * divV
         tau(3,3) = 2.0_RP * (wR-wL)*nHat(3) - (1.0_RP/3.0_RP) * divV

         tau(1,2) = (uR-uL)*nHat(2) + (vR-vL)*nHat(1)
         tau(1,3) = (uR-uL)*nHat(3) + (wR-wL)*nHat(1)
         tau(2,3) = (vR-vL)*nHat(3) + (wR-wL)*nHat(2)

         tau(2,1) = tau(1,2)
         tau(3,1) = tau(1,3)
         tau(3,2) = tau(2,3)
!
!        Compute thermal diffusivity
!        ---------------------------
         call getThermalConductivity(1.0_RP,3.0_RP / 4.0_RP, kappa)
         call equationOfState(pL, rhoL, TL)
         call equationOfState(pR, rhoR, TR)

         stab(1) = 0.0_RP
         stab(2) = maxval(abs(lambda)) * (tau(1,1)*nHat(1) + tau(1,2)*nHat(2) + tau(1,3)*nHat(3))
         stab(3) = maxval(abs(lambda)) * (tau(2,1)*nHat(1) + tau(2,2)*nHat(2) + tau(2,3)*nHat(3))
         stab(4) = maxval(abs(lambda)) * (tau(3,1)*nHat(1) + tau(3,2)*nHat(2) + tau(3,3)*nHat(3))
         stab(5) = maxval(abs(lambda)) *(( u * tau(1,1) + v*tau(1,2) + w*tau(1,3) ) * nHat(1) &
                                       + ( u * tau(2,1) + v*tau(2,2) + w*tau(2,3) ) * nHat(2) &
                                       + ( u * tau(3,1) + v*tau(3,2) + w*tau(3,3) ) * nHat(3) + kappa * (TR-TL))
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)
!
!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux = flux - lambdaStab * stab

         end associate

      end subroutine ViscousNSRiemannSolver
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LxFRiemannSolver( QLeft, QRight, nHat, t1, t2, flux )
         implicit none
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         real(kind=RP)  :: rhoL, rhouL, rhovL, rhowL, rhoeL, pL, aL, rhoV2L
         real(kind=RP)  :: rhoR, rhouR, rhovR, rhowR, rhoeR, pR, aR, rhoV2R
         real(kind=RP)  :: QLRot(5), QRRot(5)
         real(kind=RP)  :: invRhoL, invRhoR
         real(kind=RP)  :: lambda, stab(5)

         associate(gamma => thermodynamics % gamma, gm1 => thermodynamics % gammaMinus1)
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(1)                  ; rhoR = QRight(1)
         invRhoL = 1.0_RP/ rhoL           ; invRhoR = 1.0_RP / rhoR

         rhouL = QLeft(2)  * nHat(1) + QLeft(3)  * nHat(2) + QLeft(4)  * nHat(3)
         rhouR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

         rhovL = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
         rhovR = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

         rhowL = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
         rhowR = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

         rhoV2L = (POW2(rhouL) + POW2(rhovL) + POW2(rhowL)) * invRhoL
         rhoV2R = (POW2(rhouR) + POW2(rhovR) + POW2(rhowR)) * invRhoR

         rhoeL = QLeft(5) ; rhoeR = QRight(5)

         pL = gm1 * (rhoeL - 0.5_RP * rhoV2L)
         pR = gm1 * (rhoeR - 0.5_RP * rhoV2R)

         aL = sqrt(gamma * pL * invRhoL)
         aR = sqrt(gamma * pR * invRhoR)
!
!        Eigenvalues: lambda = max(|uL|,|uR|) + max(aL,aR)
!        -----------
     !!!    lambda = max(abs(rhouL*invRhoL),abs(rhouR*invRhoR)) + max(aL, aR)   ! This is a more dissipative version (not consistent with the Jacobian below)
         lambda = max(abs(rhouL*invRhoL) + aL,abs(rhouR*invRhoR) + aR)
!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         QLRot = (/ rhoL, rhouL, rhovL, rhowL, rhoeL /)
         QRRot = (/ rhoR, rhouR, rhovR, rhowR, rhoeR /)
         call AveragedStates(QLRot, QRRot, pL, pR, invRhoL, invRhoR, flux)
!
!        Compute the Lax-Friedrichs stabilization
!        ----------------------------------------
         stab = 0.5_RP * lambda * (QRRot - QLRot)
!
!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux = flux - lambdaStab * stab
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

         end associate

      END SUBROUTINE LxFRiemannSolver

      SUBROUTINE u_dissRiemannSolver( QLeft, QRight, nHat, t1, t2, flux )
         implicit none
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         real(kind=RP)  :: rhoL, rhouL, rhovL, rhowL, rhoeL, pL, aL, rhoV2L
         real(kind=RP)  :: rhoR, rhouR, rhovR, rhowR, rhoeR, pR, aR, rhoV2R
         real(kind=RP)  :: QLRot(5), QRRot(5)
         real(kind=RP)  :: invRhoL, invRhoR
         real(kind=RP)  :: lambda, stab(5)

         associate(gamma => thermodynamics % gamma, gm1 => thermodynamics % gammaMinus1)
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(1)                  ; rhoR = QRight(1)
         invRhoL = 1.0_RP/ rhoL           ; invRhoR = 1.0_RP / rhoR

         rhouL = QLeft(2)  * nHat(1) + QLeft(3)  * nHat(2) + QLeft(4)  * nHat(3)
         rhouR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

         rhovL = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
         rhovR = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

         rhowL = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
         rhowR = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

         rhoV2L = (POW2(rhouL) + POW2(rhovL) + POW2(rhowL)) * invRhoL
         rhoV2R = (POW2(rhouR) + POW2(rhovR) + POW2(rhowR)) * invRhoR

         rhoeL = QLeft(5) ; rhoeR = QRight(5)

         pL = gm1 * (rhoeL - 0.5_RP * rhoV2L)
         pR = gm1 * (rhoeR - 0.5_RP * rhoV2R)

         aL = sqrt(gamma * pL * invRhoL)
         aR = sqrt(gamma * pR * invRhoR)
!
!        Eigenvalues: lambda = max(|uL|,|uR|)
!        -----------
         lambda = max(abs(rhouL*invRhoL),abs(rhouR*invRhoR))
!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         QLRot = (/ rhoL, rhouL, rhovL, rhowL, rhoeL /)
         QRRot = (/ rhoR, rhouR, rhovR, rhowR, rhoeR /)
         call AveragedStates(QLRot, QRRot, pL, pR, invRhoL, invRhoR, flux)
!
!        Compute the Lax-Friedrichs stabilization
!        ----------------------------------------
         stab = 0.5_RP * lambda * (QRRot - QLRot)
!
!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux = flux - lambdaStab * stab
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

         end associate

      END SUBROUTINE u_dissRiemannSolver

!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------
!     Jacobian of the Lax-Friedrichs numerical flux
!     ---------------------------------------------
      subroutine LxFRiemannSolver_dFdQ(ql,qr,nHat,dfdq_num,side)
         implicit none
         !--------------------------------------------
         real(kind=RP), intent (in)  :: ql(NCONS)                 !<  Current solution on the left
         real(kind=RP), intent (in)  :: qr(NCONS)                 !<  Current solution on the right
         real(kind=RP), intent (in)  :: nHat(NDIM)                !<  Normal vector
         real(kind=RP), intent(out)  :: dfdq_num(NCONS,NCONS)     !>  Numerical flux Jacobian
         integer      , intent (in)  :: side                      !<  Either LEFT or RIGHT
         !--------------------------------------------
         real(kind=RP), dimension(NCONS,NCONS) :: dfdq,dgdq,dhdq  ! Flux Jacobians in every direction
         real(kind=RP) :: lambda    ! Lax-Friedrichs constant (penalty term for jump)
         real(kind=RP) :: lambdaL   ! Lax-Friedrichs constant on the left
         real(kind=RP) :: lambdaR   ! Lax-Friedrichs constant on the right
         REAL(KIND=RP) :: ul, vl, wl, pL, velL, aL, srhoL, rhoV2L   ! Quantities on the left
         REAL(KIND=RP) :: ur, vr, wr, pR, velR, aR, srhoR, rhoV2R   ! Quantities on the right
         real(kind=RP) :: dvn_dq(NCONS)                  ! Derivative of the normal velocity
         real(kind=RP) :: da_dq (NCONS)                  ! Derivative of the speed of sound
         real(kind=RP) :: dlambda_dq     (1,NCONS)       ! Derivative of Lax-Friedrichs constant (row matrix)
         real(kind=RP) :: q_jump     (NCONS,1)           ! Solution jump (column matrix)
         real(kind=RP) :: lambdaTerm (NCONS,NCONS)       ! Matrix contribution of lambda terms
         integer       :: i                              ! Counter
         !--------------------------------------------

         associate( gammaMinus1 => thermodynamics % gammaMinus1, &
                    gamma => thermodynamics % gamma )

         srhoL = 1._RP / ql(IRHO)
         ul = ql(IRHOU) * srhoL
         vl = ql(IRHOV) * srhoL
         wl = ql(IRHOW) * srhoL
         rhoV2L = srhoL * (ql(IRHOU)**2 + ql(IRHOV)**2 + ql(IRHOW)**2)
         pL = gammaMinus1 * (ql(IRHOE) - 0.5d0 * rhoV2L)

         srhoR = 1._RP / qr(IRHO)
         ur = qr(IRHOU) * srhoR
         vr = qr(IRHOV) * srhoR
         wr = qr(IRHOW) * srhoR
         rhoV2R = srhoR * (qr(IRHOU)**2 + qr(IRHOV)**2 + qr(IRHOW)**2)
         pR = gammaMinus1 * (qr(IRHOE) - 0.5d0 * rhoV2R)

         velL = nHat(1)*ul + nHat(2)*vl + nHat(3)*wl
         velR = nHat(1)*ur + nHat(2)*vr + nHat(3)*wr
         aL = SQRT( gamma*pL * srhoL )
         aR = SQRT( gamma*pR * srhoR )

         lambdaL = abs(velL) + aL
         lambdaR = abs(velR) + aR
         lambda  = max(lambdaL,lambdaR)

!
!        Compute matricial terms depending on side
!        -----------------------------------------

         select case(side)
            case(LEFT)
               call InviscidJacobian(ql,dfdq,dgdq,dhdq)

               if (lambdaL .ge. lambdaR) then

                  dvn_dq = (/ -velL * srhoL, nHat(1) * srhoL, nHat(2) * srhoL, nHat(3) * srhoL, 0._RP /)
                  da_dq  = (/ ( rhoV2L -ql(IRHOE) )*srhoL**2 , -ul * srhoL , -vl * srhoL , -wl * srhoL , srhoL /) * &
                                                                                                gamma * gammaMinus1 / (2._RP * aL)
                  dlambda_dq(1,:) = sign(1._RP,velL) * dvn_dq + da_dq
                  q_jump    (:,1) = ql - qr

                  lambdaTerm = matmul(q_jump, dlambda_dq)
               else
                  lambdaTerm = 0._RP
               end if

               ! Shift with lambda
               do i = 1, NCONS
                  lambdaTerm(i,i) = lambdaTerm(i,i) + lambda
               end do

            case(RIGHT)
               call InviscidJacobian(qr,dfdq,dgdq,dhdq)

               if (lambdaR .gt. lambdaL) then

                  dvn_dq = (/ -velR * srhoR, nHat(1) * srhoR, nHat(2) * srhoR, nHat(3) * srhoR, 0._RP /)
                  da_dq  = (/ ( rhoV2R - qr(IRHOE) )*srhoR**2 , -ur * srhoR , -vr * srhoR , -wr * srhoR , srhoR /) * &
                                                                                                gamma * gammaMinus1 / (2._RP * aR)
                  dlambda_dq(1,:) = sign(1._RP,velR) * dvn_dq + da_dq
                  q_jump    (:,1) = ql - qr

                  lambdaTerm = matmul(q_jump, dlambda_dq)
               else
                  lambdaTerm = 0._RP
               end if

               ! Shift with lambda
               do i = 1, NCONS
                  lambdaTerm(i,i) = lambdaTerm(i,i) - lambda
               end do

            case default
               error stop 'LxFRiemannSolver_dFdQ: side must be LEFT or RIGHT'
         end select

!
!        Apply the lambda stabilization
!        ------------------------------
         lambdaTerm = lambdaTerm * lambdaStab

!
!        Finish up
!        ---------
         dfdq_num = 0.5_RP * (dfdq * nHat(1) + dgdq * nHat(2) + dhdq * nHat(3) + lambdaTerm)

         end associate
      end subroutine LxFRiemannSolver_dFdQ
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE RoeRiemannSolver( QLeft, QRight, nHat, t1, t2, flux )
!
!        **************************************************************
!           This Roe Riemann solver implementation does not uses the
!           averaging function, nor the lambda stabilization
!        **************************************************************
!
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         REAL(KIND=RP) :: rho , rhou , rhov , rhow  , rhoe
         REAL(KIND=RP) :: rhon, rhoun, rhovn, rhown , rhoen
         REAL(KIND=RP) :: ul  , vl   , wl   , pleft , ql  , hl  , betal
         REAL(KIND=RP) :: ur  , vr   , wr   , pright, qr  , hr  , betar
         REAL(KIND=RP) :: rtd , utd  , vtd  , wtd   , htd , atd2, atd, qtd
         REAL(KIND=RP) :: dw1 , sp1  , sp1m , hd1m  , eta1, udw1, rql
         REAL(KIND=RP) :: dw4 , sp4  , sp4p , hd4   , eta4, udw4, rqr
         REAL(KIND=RP)                   :: ds = 1.0_RP

         associate ( gamma => thermodynamics % gamma )

         rho  = QLeft(1)
         rhou = QLeft(2)
         rhov = QLeft(3)
         rhow = QLeft(4)
         rhoe = QLeft(5)

         rhon  = QRight(1)
         rhoun = QRight(2)
         rhovn = QRight(3)
         rhown = QRight(4)
         rhoen = QRight(5)

         ul = rhou/rho
         vl = rhov/rho
         wl = rhow/rho
         pleft = (gamma-1._RP)*(rhoe - 0.5_RP/rho*                        &
        &                           (rhou**2 + rhov**2 + rhow**2 ))
!
         ur = rhoun/rhon
         vr = rhovn/rhon
         wr = rhown/rhon
         pright = (gamma-1._RP)*(rhoen - 0.5_RP/rhon*                    &
        &                           (rhoun**2 + rhovn**2+ rhown**2))
!
         ql = nHat(1)*ul + nHat(2)*vl + nHat(3)*wl
         qr = nHat(1)*ur + nHat(2)*vr + nHat(3)*wr
         hl = 0.5_RP*(ul*ul + vl*vl + wl*wl) +                               &
        &                 gamma/(gamma-1._RP)*pleft/rho
         hr = 0.5_RP*(ur*ur + vr*vr + wr*wr) +                               &
        &                  gamma/(gamma-1._RP)*pright/rhon
!
!        ---------------------
!        Square root averaging
!        ---------------------
!
         rtd = sqrt(rho*rhon)
         betal = rho/(rho + rtd)
         betar = 1._RP - betal
         utd = betal*ul + betar*ur
         vtd = betal*vl + betar*vr
         wtd = betal*wl + betar*wr
         htd = betal*hl + betar*hr
         atd2 = (gamma-1._RP)*(htd - 0.5_RP*(utd*utd + vtd*vtd + wtd*wtd))
         atd = sqrt(atd2)
         qtd = utd*nHat(1) + vtd*nHat(2)  + wtd*nHat(3)
!
         IF(qtd >= 0.0_RP)     THEN

            dw1 = 0.5_RP*((pright - pleft)/atd2 - (qr - ql)*rtd/atd)
            sp1 = qtd - atd
            sp1m = min(sp1,0.0_RP)
            hd1m = ((gamma+1._RP)/4._RP*atd/rtd)*dw1
            eta1 = max(-abs(sp1) - hd1m,0.0_RP)
            udw1 = dw1*(sp1m - 0.5_RP*eta1)
            rql = rho*ql
            flux(1) = ds*(rql + udw1)
            flux(2) = ds*(rql*ul + pleft*nHat(1) + udw1*(utd - atd*nHat(1)))
            flux(3) = ds*(rql*vl + pleft*nHat(2) + udw1*(vtd - atd*nHat(2)))
            flux(4) = ds*(rql*wl + pleft*nHat(3) + udw1*(wtd - atd*nHat(3)))
            flux(5) = ds*(rql*hl + udw1*(htd - qtd*atd))

         ELSE

            dw4 = 0.5_RP*((pright - pleft)/atd2 + (qr - ql)*rtd/atd)
            sp4 = qtd + atd
            sp4p = max(sp4,0.0_RP)
            hd4 = ((gamma+1._RP)/4._RP*atd/rtd)*dw4
            eta4 = max(-abs(sp4) + hd4,0.0_RP)
            udw4 = dw4*(sp4p + 0.5_RP*eta4)
            rqr = rhon*qr
            flux(1) = ds*(rqr - udw4)
            flux(2) = ds*(rqr*ur + pright*nHat(1) - udw4*(utd + atd*nHat(1)))
            flux(3) = ds*(rqr*vr + pright*nHat(2) - udw4*(vtd + atd*nHat(2)))
            flux(4) = ds*(rqr*wr + pright*nHat(3) - udw4*(wtd + atd*nHat(3)))
            flux(5) = ds*(rqr*hr - udw4*(htd + qtd*atd))
         ENDIF

         end associate

      END SUBROUTINE RoeRiemannSolver

!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE RusanovRiemannSolver( QLeft, QRight, nHat, t1, t2, flux )

         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         REAL(KIND=RP) :: rho , rhou , rhov , rhow  , rhoe
         REAL(KIND=RP) :: rhon, rhoun, rhovn, rhown , rhoen
         REAL(KIND=RP) :: ul  , vl   , wl   , pleft , ql  , hl  , betal, al, al2
         REAL(KIND=RP) :: ur  , vr   , wr   , pright, qr  , hr  , betar, ar, ar2
         REAL(KIND=RP) :: rtd , utd  , vtd  , wtd   , htd , atd2, atd, qtd
         REAL(KIND=RP) :: dw1 , sp1  , sp1m , hd1m  , eta1, udw1, rql
         REAL(KIND=RP) :: dw4 , sp4  , sp4p , hd4   , eta4, udw4, rqr
         REAL(KIND=RP)                   :: ds = 1.0_RP

         REAL(KIND=RP) :: smax, smaxL, smaxR
         REAL(KIND=RP) :: Leigen(2), Reigen(2)

         associate ( gamma => thermodynamics % gamma )

         rho  = QLeft(1)
         rhou = QLeft(2)
         rhov = QLeft(3)
         rhow = QLeft(4)
         rhoe = QLeft(5)

         rhon  = QRight(1)
         rhoun = QRight(2)
         rhovn = QRight(3)
         rhown = QRight(4)
         rhoen = QRight(5)

         ul = rhou/rho
         vl = rhov/rho
         wl = rhow/rho
         pleft = (gamma-1._RP)*(rhoe - 0.5_RP/rho*                        &
        &                           (rhou**2 + rhov**2 + rhow**2 ))
!
         ur = rhoun/rhon
         vr = rhovn/rhon
         wr = rhown/rhon
         pright = (gamma-1._RP)*(rhoen - 0.5_RP/rhon*                    &
        &                           (rhoun**2 + rhovn**2+ rhown**2))
!
         ql = nHat(1)*ul + nHat(2)*vl + nHat(3)*wl
         qr = nHat(1)*ur + nHat(2)*vr + nHat(3)*wr
         hl = 0.5_RP*(ul*ul + vl*vl + wl*wl) +                               &
        &                 gamma/(gamma-1._RP)*pleft/rho
         hr = 0.5_RP*(ur*ur + vr*vr + wr*wr) +                               &
        &                  gamma/(gamma-1._RP)*pright/rhon
!
!        ---------------------
!        Square root averaging
!        ---------------------
!
         rtd = sqrt(rho*rhon)
         betal = rho/(rho + rtd)
         betar = 1._RP - betal
         utd = betal*ul + betar*ur
         vtd = betal*vl + betar*vr
         wtd = betal*wl + betar*wr
         htd = betal*hl + betar*hr
         atd2 = (gamma-1._RP)*(htd - 0.5_RP*(utd*utd + vtd*vtd + wtd*wtd))
         atd = sqrt(atd2)
         qtd = utd*nHat(1) + vtd*nHat(2)  + wtd*nHat(3)
         !Rusanov
         ar2 = (gamma-1.d0)*(hr - 0.5d0*(ur*ur + vr*vr + wr*wr))
         al2 = (gamma-1.d0)*(hl - 0.5d0*(ul*ul + vl*vl + wl*wl))
         ar = sqrt(ar2)
         al = sqrt(al2)
!
         rql = rho*ql
         rqr = rhon*qr
         flux(1) = ds*(rql + rqr)
         flux(2) = ds*(rql*ul + pleft*nHat(1) + rqr*ur + pright*nHat(1))
         flux(3) = ds*(rql*vl + pleft*nHat(2) + rqr*vr + pright*nHat(2))
         flux(4) = ds*(rql*wl + pleft*nHat(3) + rqr*wr + pright*nHat(3))
         flux(5) = ds*(rql*hl + rqr*hr)

         smax = MAX(ar+ABS(qr),al+ABS(ql))

         flux = (flux - ds*smax*(QRight-QLeft))/2.d0

         RETURN

         end associate

      END SUBROUTINE RusanovRiemannSolver
!
!////////////////////////////////////////////////////////////////////////////////////////////
!
!        Averaged states functions
!        -------------------------
!
!     To this averages, the states QLeft and QRight velocities have already been rotated,
!  so that only the first flux component has to be computed (Rotational invariance, see
!  E.F. Toro - Riemann solvers and numerical methods for fluid dynamics. Page 105).
!
!  Implemented two-point averages are:
!     -> Standard: Central fluxes
!     -> Ducros
!     -> Morinishi
!     -> Kennedy and Gruber
!     -> Pirozzoli
!     -> Entropy conserving
!     -> Chandrasekar
!
!////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine StandardAverage(QLeft, QRight, pL, pR, invRhoL, invRhoR, flux)
!
!        *********************************************************************
!           Computes the standard average of the two states:
!              F* = {{F}} = 0.5 * (FL + FR)
!
!           State vectors are rotated.
!        *********************************************************************
!
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: uL, vL, wL
         real(kind=RP)     :: uR, vR, wR

         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = 0.5_RP * ( QLeft(IRHOU) + QRight(IRHOU) )
         flux(IRHOU) = 0.5_RP * ( QLeft(IRHOU) * uL + QRight(IRHOU) * uR + pL + pR )
         flux(IRHOV) = 0.5_RP * ( QLeft(IRHOU) * vL + QRight(IRHOU) * vR )
         flux(IRHOW) = 0.5_RP * ( QLeft(IRHOU) * wL + QRight(IRHOU) * wR )
         flux(IRHOE) = 0.5_RP * ( uL*(QLeft(IRHOE) + pL) + uR*(QRight(IRHOE) + pR) )

      end subroutine StandardAverage

      subroutine DucrosAverage(QLeft, QRight, pL, pR, invRhoL, invRhoR, flux)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: uL, vL, wL
         real(kind=RP)     :: uR, vR, wR

         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = 0.25_RP * ( QLeft(IRHO) + QRight(IRHO) ) * (uL + uR)
         flux(IRHOU) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * (uL + uR) + 0.5_RP * (pL + pR)
         flux(IRHOV) = 0.25_RP * ( QLeft(IRHOV) + QRight(IRHOV) ) * (uL + uR)
         flux(IRHOW) = 0.25_RP * ( QLeft(IRHOW) + QRight(IRHOW) ) * (uL + uR)
         flux(IRHOE) = 0.25_RP * ( QLeft(IRHOE) + pL + QRight(IRHOE) + pR ) * (uL + uR)

      end subroutine DucrosAverage

      subroutine MorinishiAverage(QLeft,QRight, pL, pR, invRhoL, invRhoR, flux)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: uL, vL, wL, hL
         real(kind=RP)     :: uR, vR, wR, hR

         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
!
!        Here the enthalpy does not contain the kinetic energy
!        -----------------------------------------------------
         hL = dimensionless % cp * pL  ; hR = dimensionless % cp * pR
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = 0.5_RP * ( QLeft(IRHOU) + QRight(IRHOU) )
         flux(IRHOU) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * ( uL + uR ) + 0.5_RP * ( pL + pR )
         flux(IRHOV) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * ( vL + vR )
         flux(IRHOW) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * ( wL + wR )
         flux(IRHOE) = 0.5_RP * ( uL*hL + uR*hR) + 0.25_RP * ( QLeft(IRHOU)*uL + QRight(IRHOU)*uR ) * ( uL + uR ) &
                                                 + 0.25_RP * ( QLeft(IRHOU)*vL + QRight(IRHOU)*vR ) * ( vL + vR ) &
                                                 + 0.25_RP * ( QLeft(IRHOU)*wL + QRight(IRHOU)*wR ) * ( wL + wR ) &
                                                 - 0.25_RP * ( QLeft(IRHOU)*POW2(uL) + QRight(IRHOU)*POW2(uR)   ) &
                                                 - 0.25_RP * ( QLeft(IRHOU)*POW2(vL) + QRight(IRHOU)*POW2(vR)   ) &
                                                 - 0.25_RP * ( QLeft(IRHOU)*POW2(wL) + QRight(IRHOU)*POW2(wR)   )

      end subroutine MorinishiAverage

      subroutine KennedyGruberAverage(QLeft,QRight, pL, pR, invRhoL, invRhoR, flux)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: uL, vL, wL
         real(kind=RP)     :: uR, vR, wR
         real(kind=RP)     :: rho, u, v, w, e, p

         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
!
!        Compute KG averages
!        -------------------
         rho = 0.5_RP * (QLeft(IRHO) + QRight(IRHO))
         u   = 0.5_RP * (uL + uR)
         v   = 0.5_RP * (vL + vR)
         w   = 0.5_RP * (wL + wR)
         p   = 0.5_RP * (pL + pR)
         e   = 0.5_RP * (QLeft(IRHOE)*invRhoL + QRight(IRHOE)*invRhoR)
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = rho * u
         flux(IRHOU) = rho * u * u + p
         flux(IRHOV) = rho * u * v
         flux(IRHOW) = rho * u * w
         flux(IRHOE) = rho * u * e + p * u

      end subroutine KennedyGruberAverage

      subroutine PirozzoliAverage(QLeft,QRight, pL, pR, invRhoL, invRhoR, flux)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: uL, vL, wL
         real(kind=RP)     :: uR, vR, wR
         real(kind=RP)     :: rho, u, v, w, h, p

         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)

         rho = 0.5_RP * (QLeft(IRHO) + QRight(IRHO))
         u   = 0.5_RP * (uL + uR)
         v   = 0.5_RP * (vL + vR)
         w   = 0.5_RP * (wL + wR)
         p   = 0.5_RP * (pL + pR)
         h   = 0.5_RP * ((QLeft(IRHOE)+pL)*invRhoL + (QRight(IRHOE)+pR)*invRhoR)
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = rho * u
         flux(IRHOU) = rho * u * u + p
         flux(IRHOV) = rho * u * v
         flux(IRHOW) = rho * u * w
         flux(IRHOE) = rho * u * h

      end subroutine PirozzoliAverage

      subroutine EntropyConservingAverage(QLeft,QRight, pL, pR, invRhoL, invRhoR, flux)
         use SMConstants
         use Utilities, only: logarithmicMean
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: rhoL, uL, vL, wL
         real(kind=RP)     :: rhoR, uR, vR, wR
         real(kind=RP)     :: rho, u, v, w, h, p, p2
         real(kind=RP)     :: zL(NCONS), zR(NCONS), zSum(NCONS), invZ1Sum
         real(kind=RP)     :: z5Log, z1Log

         associate ( gammaPlus1Div2      => thermodynamics % gammaPlus1Div2, &
                     gammaMinus1Div2     => thermodynamics % gammaMinus1Div2, &
                     gammaDivGammaMinus1 => thermodynamics % gammaDivGammaMinus1, &
                     invGamma            => thermodynamics % invGamma )

         rhoL = QLeft(IRHO)               ; rhoR = QRight(IRHO)
         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
!
!        Compute Ismail and Roe parameter vector
!        ---------------------------------------
         zL(5) = sqrt(rhoL*pL)      ; zR(5) = sqrt(rhoR*pR)
         zL(1) = rhoL / zL(5)       ; zR(1) = rhoR / zR(5)
         zL(2) = zL(1) * uL         ; zR(2) = zR(1) * uR
         zL(3) = zL(1) * vL         ; zR(3) = zR(1) * vR
         zL(4) = zL(1) * wL         ; zR(4) = zR(1) * wR

         zSum = zL + zR
         invZ1Sum = 1.0_RP / zSum(1)

         call logarithmicMean(zL(1),zR(1), z1Log)
         call logarithmicMean(zL(5),zR(5), z5Log)

         rho = 0.5_RP * zSum(1) * z5Log
         u   = zSum(2) * invZ1Sum
         v   = zSum(3) * invZ1Sum
         w   = zSum(4) * invZ1Sum
         p   = zSum(5) * invZ1Sum
         p2  = (gammaPlus1Div2 * z5Log / z1Log + gammaMinus1Div2 * p) * invGamma
         h   = gammaDivGammaMinus1 * p2 / rho + 0.5_RP*(POW2(u) + POW2(v) + POW2(w))
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = rho * u
         flux(IRHOU) = rho * u * u + p
         flux(IRHOV) = rho * u * v
         flux(IRHOW) = rho * u * w
         flux(IRHOE) = rho * u * h

         end associate

      end subroutine EntropyConservingAverage

      subroutine ChandrasekarAverage(QLeft,QRight, pL, pR, invRhoL, invRhoR, flux)
         use SMConstants
         use Utilities, only: logarithmicMean
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: rhoL, uL, vL, wL, betaL
         real(kind=RP)     :: rhoR, uR, vR, wR, betaR
         real(kind=RP)     :: rho, u, v, w, h, p, betaLog

         associate ( gammaMinus1 => thermodynamics % gammaMinus1 )

         rhoL = QLeft(IRHO)               ; rhoR = QRight(IRHO)
         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
!
!        Compute Chandrasekar's variables
!        --------------------------------
         betaL = 0.5_RP * rhoL / pL    ; betaR = 0.5_RP * rhoR / pR
         call logarithmicMean(betaL, betaR, betaLog)

         call logarithmicMean(rhoL,rhoR,rho)
         u   = AVERAGE(uL, uR)
         v   = AVERAGE(vL, vR)
         w   = AVERAGE(wL, wR)
         p   = 0.5_RP * (rhoL + rhoR) / (betaL + betaR)
         h   =   0.5_RP/(betaLog*(gammaMinus1)) &
               - 0.5_RP*AVERAGE(POW2(uL)+POW2(vL)+POW2(wL), POW2(uR)+POW2(vR)+POW2(wR)) &
               + p/rho + POW2(u) + POW2(v) + POW2(w)
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = rho * u
         flux(IRHOU) = rho * u * u + p
         flux(IRHOV) = rho * u * v
         flux(IRHOW) = rho * u * w
         flux(IRHOE) = rho * u * h

         end associate

      end subroutine ChandrasekarAverage
!
!///////////////////////////////////////////////////////////////////////
!
!        Volumetric two-point fluxes
!        ---------------------------
!///////////////////////////////////////////////////////////////////////
!
      subroutine StandardDG_TwoPointFlux(QL,QR,JaL,JaR, fSharp)
         use SMConstants
         use PhysicsStorage_NS
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: f(NCONS), g(NCONS), h(NCONS)

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)

         pL = thermodynamics % GammaMinus1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                                                      + QL(IRHOV) * vL &
                                                                      + QL(IRHOW) * wL ))

         pR = thermodynamics % GammaMinus1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                                                      + QR(IRHOV) * vR &
                                                                      + QR(IRHOW) * wR ))
!
!        Average metrics: (Note: Here all average (1/2)s are accounted later)
!        ---------------
         Ja = (JaL + JaR)
!
!        Compute the flux
!        ----------------
         f(IRHO)  = ( QL(IRHOU) + QR(IRHOU) )
         f(IRHOU) = ( QL(IRHOU) * uL + QR(IRHOU) * uR + pL + pR )
         f(IRHOV) = ( QL(IRHOU) * vL + QR(IRHOU) * vR )
         f(IRHOW) = ( QL(IRHOU) * wL + QR(IRHOU) * wR )
         f(IRHOE) = ( uL*(QL(IRHOE) + pL) + uR*(QR(IRHOE) + pR) )

         g(IRHO)  = (QL(IRHOV) + QR(IRHOV))
         g(IRHOU) = ( QL(IRHOV) * uL + QR(IRHOV) * uR )
         g(IRHOV) = ( QL(IRHOV) * vL + QR(IRHOV) * vR + pL + pR )
         g(IRHOW) = ( QL(IRHOV) * wL + QR(IRHOV) * wR )
         g(IRHOE) = ( vL*(QL(IRHOE) + pL) + vR*(QR(IRHOE) + pR) )

         h(IRHO)  = (QL(IRHOW) + QR(IRHOW))
         h(IRHOU) = (QL(IRHOW)*uL + QR(IRHOW)*uR)
         h(IRHOV) = (QL(IRHOW)*vL + QR(IRHOW)*vR)
         h(IRHOW) = (QL(IRHOW)*wL + QR(IRHOW)*wR + pL + pR)
         h(IRHOE) = ( wL*(QL(IRHOE) + pL) + wR*(QR(IRHOE) + pR) )
!
!        Compute the sharp flux: (And account for the (1/2)^2)
!        ----------------------
         fSharp = 0.25_RP * ( f*Ja(IX) + g*Ja(IY) + h*Ja(IZ) )

      end subroutine StandardDG_TwoPointFlux

      subroutine Morinishi_TwoPointFlux(QL,QR,JaL,JaR,fSharp)
         use SMConstants
         use PhysicsStorage_NS
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL, hL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR, hR
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: f(NCONS), g(NCONS), h(NCONS)

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)

         associate(gm1 => thermodynamics % GammaMinus1)

         pL = gm1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                             + QL(IRHOV) * vL &
                                             + QL(IRHOW) * wL ))

         pR = gm1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                             + QR(IRHOV) * vR &
                                             + QR(IRHOW) * wR ))
         end associate
!
!        Here the enthalpy does not contain the kinetic energy
!        -----------------------------------------------------
         hL = dimensionless % cp * pL  ; hR = dimensionless % cp * pR
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         f(IRHO)  = 0.5_RP * ( QL(IRHOU) + QR(IRHOU) )
         f(IRHOU) = 0.25_RP * ( QL(IRHOU) + QR(IRHOU) ) * ( uL + uR ) + 0.5_RP * ( pL + pR )
         f(IRHOV) = 0.25_RP * ( QL(IRHOU) + QR(IRHOU) ) * ( vL + vR )
         f(IRHOW) = 0.25_RP * ( QL(IRHOU) + QR(IRHOU) ) * ( wL + wR )
         f(IRHOE) = 0.5_RP * ( uL*hL + uR*hR) + 0.25_RP * ( QL(IRHOU)*uL + QR(IRHOU)*uR ) * ( uL + uR ) &
                                              + 0.25_RP * ( QL(IRHOU)*vL + QR(IRHOU)*vR ) * ( vL + vR ) &
                                              + 0.25_RP * ( QL(IRHOU)*wL + QR(IRHOU)*wR ) * ( wL + wR ) &
                                              - 0.25_RP * ( QL(IRHOU)*POW2(uL) + QR(IRHOU)*POW2(uR)   ) &
                                              - 0.25_RP * ( QL(IRHOU)*POW2(vL) + QR(IRHOU)*POW2(vR)   ) &
                                              - 0.25_RP * ( QL(IRHOU)*POW2(wL) + QR(IRHOU)*POW2(wR)   )

         g(IRHO)  = 0.5_RP * ( QL(IRHOV) + QR(IRHOV) )
         g(IRHOU) = 0.25_RP * ( QL(IRHOV) + QR(IRHOV) ) * ( uL + uR )
         g(IRHOV) = 0.25_RP * ( QL(IRHOV) + QR(IRHOV) ) * ( vL + vR ) + 0.5_RP * ( pL + pR )
         g(IRHOW) = 0.25_RP * ( QL(IRHOV) + QR(IRHOV) ) * ( wL + wR )
         g(IRHOE) = 0.5_RP * ( vL*hL + vR*hR) + 0.25_RP * ( QL(IRHOV)*uL + QR(IRHOV)*uR ) * ( uL + uR ) &
                                              + 0.25_RP * ( QL(IRHOV)*vL + QR(IRHOV)*vR ) * ( vL + vR ) &
                                              + 0.25_RP * ( QL(IRHOV)*wL + QR(IRHOV)*wR ) * ( wL + wR ) &
                                              - 0.25_RP * ( QL(IRHOV)*POW2(uL) + QR(IRHOV)*POW2(uR)   ) &
                                              - 0.25_RP * ( QL(IRHOV)*POW2(vL) + QR(IRHOV)*POW2(vR)   ) &
                                              - 0.25_RP * ( QL(IRHOV)*POW2(wL) + QR(IRHOV)*POW2(wR)   )

         h(IRHO)  = 0.5_RP * ( QL(IRHOW) + QR(IRHOW) )
         h(IRHOU) = 0.25_RP * ( QL(IRHOW) + QR(IRHOW) ) * ( uL + uR )
         h(IRHOV) = 0.25_RP * ( QL(IRHOW) + QR(IRHOW) ) * ( vL + vR )
         h(IRHOW) = 0.25_RP * ( QL(IRHOW) + QR(IRHOW) ) * ( wL + wR ) + 0.5_RP * ( pL + pR )
         h(IRHOE) = 0.5_RP * ( wL*hL + wR*hR) + 0.25_RP * ( QL(IRHOW)*uL + QR(IRHOW)*uR ) * ( uL + uR ) &
                                              + 0.25_RP * ( QL(IRHOW)*vL + QR(IRHOW)*vR ) * ( vL + vR ) &
                                              + 0.25_RP * ( QL(IRHOW)*wL + QR(IRHOW)*wR ) * ( wL + wR ) &
                                              - 0.25_RP * ( QL(IRHOW)*POW2(uL) + QR(IRHOW)*POW2(uR)   ) &
                                              - 0.25_RP * ( QL(IRHOW)*POW2(vL) + QR(IRHOW)*POW2(vR)   ) &
                                              - 0.25_RP * ( QL(IRHOW)*POW2(wL) + QR(IRHOW)*POW2(wR)   )
!
!        Compute the sharp flux
!        ----------------------
         fSharp = f*Ja(IX) + g*Ja(IY) + h*Ja(IZ)

      end subroutine Morinishi_TwoPointFlux

      subroutine Ducros_TwoPointFlux(QL,QR,JaL,JaR,fSharp)
         use SMConstants
         use PhysicsStorage_NS
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: f(NCONS), g(NCONS), h(NCONS)

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)

         associate(gm1 => thermodynamics % GammaMinus1)

         pL = gm1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                             + QL(IRHOV) * vL &
                                             + QL(IRHOW) * wL ))

         pR = gm1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                             + QR(IRHOV) * vR &
                                             + QR(IRHOW) * wR ))
         end associate
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         f(IRHO)  = 0.25_RP * ( QL(IRHO) + QR(IRHO) ) * (uL + uR)
         f(IRHOU) = 0.25_RP * ( QL(IRHOU) + QR(IRHOU) ) * (uL + uR) + 0.5_RP * (pL + pR)
         f(IRHOV) = 0.25_RP * ( QL(IRHOV) + QR(IRHOV) ) * (uL + uR)
         f(IRHOW) = 0.25_RP * ( QL(IRHOW) + QR(IRHOW) ) * (uL + uR)
         f(IRHOE) = 0.25_RP * ( QL(IRHOE) + pL + QR(IRHOE) + pR ) * (uL + uR)

         g(IRHO)  = 0.25_RP * ( QL(IRHO) + QR(IRHO) ) * (vL + vR)
         g(IRHOU) = 0.25_RP * ( QL(IRHOU) + QR(IRHOU) ) * (vL + vR)
         g(IRHOV) = 0.25_RP * ( QL(IRHOV) + QR(IRHOV) ) * (vL + vR) + 0.5_RP * (pL + pR)
         g(IRHOW) = 0.25_RP * ( QL(IRHOW) + QR(IRHOW) ) * (vL + vR)
         g(IRHOE) = 0.25_RP * ( QL(IRHOE) + pL + QR(IRHOE) + pR ) * (vL + vR)

         h(IRHO)  = 0.25_RP * ( QL(IRHO) + QR(IRHO) ) * (wL + wR)
         h(IRHOU) = 0.25_RP * ( QL(IRHOU) + QR(IRHOU) ) * (wL + wR)
         h(IRHOV) = 0.25_RP * ( QL(IRHOV) + QR(IRHOV) ) * (wL + wR)
         h(IRHOW) = 0.25_RP * ( QL(IRHOW) + QR(IRHOW) ) * (wL + wR) + 0.5_RP * (pL + pR)
         h(IRHOE) = 0.25_RP * ( QL(IRHOE) + pL + QR(IRHOE) + pR ) * (wL + wR)
!
!        Compute the sharp flux
!        ----------------------
         fSharp = f*Ja(IX) + g*Ja(IY) + h*Ja(IZ)

      end subroutine Ducros_TwoPointFlux

      subroutine KennedyGruber_TwoPointFlux(QL,QR,JaL,JaR,fSharp)
         use SMConstants
         use PhysicsStorage_NS
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR
         real(kind=RP)     :: rho, u, v, w, e, p
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: f(NCONS), g(NCONS), h(NCONS)

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)

         associate(gm1 => thermodynamics % GammaMinus1)

         pL = gm1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                             + QL(IRHOV) * vL &
                                             + QL(IRHOW) * wL ))

         pR = gm1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                             + QR(IRHOV) * vR &
                                             + QR(IRHOW) * wR ))
         end associate

         rho = 0.5_RP * (QL(IRHO) + QR(IRHO))
         u   = 0.5_RP * (uL + uR)
         v   = 0.5_RP * (vL + vR)
         w   = 0.5_RP * (wL + wR)
         p   = 0.5_RP * (pL + pR)
         e   = 0.5_RP * (QL(IRHOE)*invRhoL + QR(IRHOE)*invRhoR)
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         f(IRHO)  = rho * u
         f(IRHOU) = rho * u * u + p
         f(IRHOV) = rho * u * v
         f(IRHOW) = rho * u * w
         f(IRHOE) = rho * u * e + p * u

         g(IRHO)  = rho * v
         g(IRHOU) = rho * v * u
         g(IRHOV) = rho * v * v + p
         g(IRHOW) = rho * v * w
         g(IRHOE) = rho * v * e + p * v

         h(IRHO)  = rho * w
         h(IRHOU) = rho * w * u
         h(IRHOV) = rho * w * v
         h(IRHOW) = rho * w * w + p
         h(IRHOE) = rho * w * e + p * w
!
!        Compute the sharp flux
!        ----------------------
         fSharp = f*Ja(IX) + g*Ja(IY) + h*Ja(IZ)

      end subroutine KennedyGruber_TwoPointFlux

      subroutine Pirozzoli_TwoPointFlux(QL,QR,JaL,JaR,fSharp)
         use SMConstants
         use PhysicsStorage_NS
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR
         real(kind=RP)     :: rho, u, v, w, h, p
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: ff(NCONS), gg(NCONS), hh(NCONS)

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)

         pL = thermodynamics % GammaMinus1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                                                      + QL(IRHOV) * vL &
                                                                      + QL(IRHOW) * wL ))

         pR = thermodynamics % GammaMinus1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                                                      + QR(IRHOV) * vR &
                                                                      + QR(IRHOW) * wR ))

         rho = 0.5_RP * (QL(IRHO) + QR(IRHO))
         u   = 0.5_RP * (uL + uR)
         v   = 0.5_RP * (vL + vR)
         w   = 0.5_RP * (wL + wR)
         p   = 0.5_RP * (pL + pR)
         h   = 0.5_RP * ((QL(IRHOE)+pL)*invRhoL + (QR(IRHOE)+pR)*invRhoR)
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         ff(IRHO)  = rho * u
         ff(IRHOU) = rho * u * u + p
         ff(IRHOV) = rho * u * v
         ff(IRHOW) = rho * u * w
         ff(IRHOE) = rho * u * h

         gg(IRHO)  = rho * v
         gg(IRHOU) = rho * v * u
         gg(IRHOV) = rho * v * v + p
         gg(IRHOW) = rho * v * w
         gg(IRHOE) = rho * v * h

         hh(IRHO)  = rho * w
         hh(IRHOU) = rho * w * u
         hh(IRHOV) = rho * w * v
         hh(IRHOW) = rho * w * w + p
         hh(IRHOE) = rho * w * h
!
!        Compute the sharp flux
!        ----------------------
         fSharp = ff*Ja(IX) + gg*Ja(IY) + hh*Ja(IZ)

      end subroutine Pirozzoli_TwoPointFlux

      subroutine EntropyConserving_TwoPointFlux(QL,QR,JaL,JaR,fSharp)
!
!        *******************************************************************
!           Entropy conserving split form by Ismail and Roe.
!           Parameter vector definition:
!
!           z1 = sqrt(rho / p)
!           z2 = z1 * u
!           z3 = z1 * v
!           z4 = z1 * w
!           z5 = sqrt(p rho)
!
!           Averaged states:
!
!           rho = {{z1}} {{z5:ln}}
!           u   = {{z2}} / {{z1}}
!           v   = {{z3}} / {{z1}}
!           w   = {{z4}} / {{z1}}
!           p   = {{z5}} / {{z1}}
!           h   = g*p2/(rho(g-1)) + 0.5*(u^2 + v^2 + w^2)
!
!           where
!           p2   = (g+1)/(2g) {{z5:ln}}/{{z1:ln}} + (g-1)/(2g){{z5}}/{{z1}}
!
!        *******************************************************************
!
         use SMConstants
         use PhysicsStorage_NS
         use Utilities, only: logarithmicMean
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, rhoL, uL, vL, wL, pL
         real(kind=RP)     :: invRhoR, rhoR, uR, vR, wR, pR
         real(kind=RP)     :: rho, u, v, w, h, p, p2
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: zL(NCONS), zR(NCONS), zAv(NCONS), invZ1Sum
         real(kind=RP)     :: z5Log, z1Log, invZ1Av
         real(kind=RP)     :: ff(NCONS), gg(NCONS), hh(NCONS)

         associate ( gammaPlus1Div2      => thermodynamics % gammaPlus1Div2, &
                     gammaMinus1Div2     => thermodynamics % gammaMinus1Div2, &
                     gammaDivGammaMinus1 => thermodynamics % gammaDivGammaMinus1, &
                     invGamma            => thermodynamics % invGamma )

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         rhoL = QL(IRHO)               ; rhoR = QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)

         pL = thermodynamics % GammaMinus1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                                                      + QL(IRHOV) * vL &
                                                                      + QL(IRHOW) * wL ))

         pR = thermodynamics % GammaMinus1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                                                      + QR(IRHOV) * vR &
                                                                      + QR(IRHOW) * wR ))
!
!        Compute Ismail and Roe parameter vector
!        ---------------------------------------
         zL(5) = sqrt(rhoL*pL)      ; zR(5) = sqrt(rhoR*pR)
         zL(1) = rhoL / zL(5)       ; zR(1) = rhoR / zR(5)
         zL(2) = zL(1) * uL         ; zR(2) = zR(1) * uR
         zL(3) = zL(1) * vL         ; zR(3) = zR(1) * vR
         zL(4) = zL(1) * wL         ; zR(4) = zR(1) * wR

         zAv = 0.5_RP * (zL + zR)
         invZ1Av = 1.0_RP / zAv(1)

         call logarithmicMean(zL(1),zR(1), z1Log)
         call logarithmicMean(zL(5),zR(5), z5Log)

         rho = zAv(1) * z5Log
         u   = zAv(2) * invZ1Av
         v   = zAv(3) * invZ1Av
         w   = zAv(4) * invZ1Av
         p   = zAv(5) * invZ1Av
         p2  = (gammaPlus1Div2 * z5Log / z1Log + gammaMinus1Div2 * p) * invGamma
         h   = gammaDivGammaMinus1 * p2 / rho + 0.5_RP*(POW2(u) + POW2(v) + POW2(w))
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         ff(IRHO)  = rho * u
         ff(IRHOU) = rho * u * u + p
         ff(IRHOV) = rho * u * v
         ff(IRHOW) = rho * u * w
         ff(IRHOE) = rho * u * h

         gg(IRHO)  = rho * v
         gg(IRHOU) = rho * v * u
         gg(IRHOV) = rho * v * v + p
         gg(IRHOW) = rho * v * w
         gg(IRHOE) = rho * v * h

         hh(IRHO)  = rho * w
         hh(IRHOU) = rho * w * u
         hh(IRHOV) = rho * w * v
         hh(IRHOW) = rho * w * w + p
         hh(IRHOE) = rho * w * h
!
!        Compute the sharp flux
!        ----------------------
         fSharp = ff*Ja(IX) + gg*Ja(IY) + hh*Ja(IZ)

         end associate

      end subroutine EntropyConserving_TwoPointFlux

      subroutine Chandrasekar_TwoPointFlux(QL,QR,JaL,JaR,fSharp)
         use SMConstants
         use PhysicsStorage_NS
         use Utilities, only: logarithmicMean
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, rhoL, uL, vL, wL, pL, betaL
         real(kind=RP)     :: invRhoR, rhoR, uR, vR, wR, pR, betaR
         real(kind=RP)     :: rho, u, v, w, h, p, betaLog
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: ff(NCONS), gg(NCONS), hh(NCONS)

         associate ( gammaMinus1 => thermodynamics % gammaMinus1 )

         invRhoL = 1.0_RP / QL(IRHO)   ; invRhoR = 1.0_RP / QR(IRHO)
         rhoL = QL(IRHO)               ; rhoR = QR(IRHO)
         uL = invRhoL * QL(IRHOU)      ; uR = invRhoR * QR(IRHOU)
         vL = invRhoL * QL(IRHOV)      ; vR = invRhoR * QR(IRHOV)
         wL = invRhoL * QL(IRHOW)      ; wR = invRhoR * QR(IRHOW)

         pL = thermodynamics % GammaMinus1 * ( QL(IRHOE) - 0.5_RP * (   QL(IRHOU) * uL &
                                                                      + QL(IRHOV) * vL &
                                                                      + QL(IRHOW) * wL ))

         pR = thermodynamics % GammaMinus1 * ( QR(IRHOE) - 0.5_RP * (   QR(IRHOU) * uR &
                                                                      + QR(IRHOV) * vR &
                                                                      + QR(IRHOW) * wR ))
!
!        Compute Chandrasekar's variables
!        --------------------------------
         betaL = 0.5_RP * rhoL / pL    ; betaR = 0.5_RP * rhoR / pR
         call logarithmicMean(betaL, betaR, betaLog)

         call logarithmicMean(rhoL,rhoR,rho)
         u   = AVERAGE(uL, uR)
         v   = AVERAGE(vL, vR)
         w   = AVERAGE(wL, wR)
         p   = 0.5_RP * (rhoL + rhoR) / (betaL + betaR)
         h   =   0.5_RP/(betaLog*(gammaMinus1)) &
               - 0.5_RP*AVERAGE(POW2(uL)+POW2(vL)+POW2(wL), POW2(uR)+POW2(vR)+POW2(wR)) &
               + p/rho + POW2(u) + POW2(v) + POW2(w)
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         ff(IRHO)  = rho * u
         ff(IRHOU) = rho * u * u + p
         ff(IRHOV) = rho * u * v
         ff(IRHOW) = rho * u * w
         ff(IRHOE) = rho * u * h

         gg(IRHO)  = rho * v
         gg(IRHOU) = rho * v * u
         gg(IRHOV) = rho * v * v + p
         gg(IRHOW) = rho * v * w
         gg(IRHOE) = rho * v * h

         hh(IRHO)  = rho * w
         hh(IRHOU) = rho * w * u
         hh(IRHOV) = rho * w * v
         hh(IRHOW) = rho * w * w + p
         hh(IRHOE) = rho * w * h
!
!        Compute the sharp flux
!        ----------------------
         fSharp = ff*Ja(IX) + gg*Ja(IY) + hh*Ja(IZ)

         end associate

      end subroutine Chandrasekar_TwoPointFlux

end module RiemannSolvers_NS