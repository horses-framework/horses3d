#include "Includes.h"
module RiemannSolvers_MUKeywordsModule

     integer,                       parameter :: KEYWORD_LENGTH          = 132
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_SOLVER_NAME_KEY = "riemann solver"
!
!    --------------------------
!    Riemann solver definitions
!    --------------------------
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_CENTRAL_NAME = "central"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_EXACT_NAME   = "exact"

     enum, bind(C)
        enumerator :: RIEMANN_CENTRAL = 1, RIEMANN_EXACT
     end enum

end module RiemannSolvers_MUKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!
module RiemannSolvers_MU
   use SMConstants
   use Physics_MU
   use PhysicsStorage_MU
   use VariableConversion_MU
   use FluidData_MU

   implicit none

   private
   public whichRiemannSolver
   public SetRiemannSolver, DescribeRiemannSolver
   public RiemannSolver, ExactRiemannSolver

   abstract interface
      subroutine RiemannSolverFCN(QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL,fR, invMa2L, invMa2R)
         use SMConstants
         use PhysicsStorage_MU
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: rhoL
         real(kind=RP), intent(in)       :: rhoR
         real(kind=RP), intent(in)       :: muL
         real(kind=RP), intent(in)       :: muR
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: fL(1:NCONS)
         real(kind=RP), intent(out)      :: fR(1:NCONS)
         real(kind=RP), intent(in)       :: invMa2L
         real(kind=RP), intent(in)       :: invMa2R
      end subroutine RiemannSolverFCN
   end interface

   procedure(RiemannSolverFCN), protected, pointer :: RiemannSolver  => NULL()
   integer,                     protected          :: whichRiemannSolver = -1
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
         use RiemannSolvers_MUKeywordsModule
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

!
!        --------------------------------------------
!        Choose the Riemann solver (default is exact)
!        --------------------------------------------
         if (controlVariables % containsKey(RIEMANN_SOLVER_NAME_KEY)) then

            keyword = controlVariables % stringValueForKey(RIEMANN_SOLVER_NAME_KEY, KEYWORD_LENGTH)
            call toLower(keyword)

            select case (keyword)
            case(RIEMANN_CENTRAL_NAME)
               RiemannSolver => CentralRiemannSolver
               whichRiemannSolver = RIEMANN_CENTRAL

            case(RIEMANN_EXACT_NAME)
               RiemannSolver => ExactRiemannSolver
               whichRiemannSolver = RIEMANN_EXACT

            case default
               print*, "Riemann Solver not recognized."
               errorMessage(STD_OUT)
               error stop

            end select

         else
!
!           Select exact by default
!           -----------------------
            RiemannSolver => CentralRiemannSolver
            whichRiemannSolver = RIEMANN_EXACT

         end if

      END SUBROUTINE SetRiemannSolver

      subroutine DescribeRiemannSolver
!
!        -------
!        Modules
!        -------
         use RiemannSolvers_MUKeywordsModule


         select case (whichRiemannSolver)
         case (RIEMANN_CENTRAL)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Central"

         case (RIEMANN_EXACT)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Exact"

         end select

      end subroutine DescribeRiemannSolver
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!        Riemann solvers
!        ---------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine CentralRiemannSolver(QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL, fR, invMa2L, invMa2R)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: rhoL
         real(kind=RP), intent(in)       :: rhoR
         real(kind=RP), intent(in)       :: muL
         real(kind=RP), intent(in)       :: muR
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: fL(1:NCONS)
         real(kind=RP), intent(out)      :: fR(1:NCONS)
         real(kind=RP), intent(in)       :: invMa2L
         real(kind=RP), intent(in)       :: invMa2R
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: cL, uL, vL, wL, pL, invSqrtRhoL
         real(kind=RP) :: cR, uR, vR, wR, pR, invSqrtRhoR
         

!
!        Left state variables and fluxes
!        -------------------------------
         invSqrtRhoL = 1.0_RP / sqrt(rhoL)
         cL = QLeft(IMC)
         uL = invSqrtRhoL * (QLeft(IMSQRHOU) * nHat(1) + QLeft(IMSQRHOV) * nHat(2) + QLeft(IMSQRHOW) * nHat(3))
         vL = invSqrtRhoL * (QLeft(IMSQRHOU) * t1(1)   + QLeft(IMSQRHOV) * t1(2)   + QLeft(IMSQRHOW) * t1(3))
         wL = invSqrtRhoL * (QLeft(IMSQRHOU) * t2(1)   + QLeft(IMSQRHOV) * t2(2)   + QLeft(IMSQRHOW) * t2(3))
         pL = QLeft(IMP)

         fL(IMC)      = uL*cL
         fL(IMSQRHOU) = 0.5_RP*rhoL*uL*uL + pL
         fL(IMSQRHOV) = 0.5_RP*rhoL*uL*vL
         fL(IMSQRHOW) = 0.5_RP*rhoL*uL*wL
         fL(IMP)      = 0.0_RP

!
!        Right state variables and fluxes
!        --------------------------------
         invSqrtRhoR = 1.0_RP / sqrt(rhoR)
         cR = QRight(IMC)
         uR = invSqrtRhoR * (QRight(IMSQRHOU) * nHat(1) + QRight(IMSQRHOV) * nHat(2) + QRight(IMSQRHOW) * nHat(3))
         vR = invSqrtRhoR * (QRight(IMSQRHOU) * t1(1)   + QRight(IMSQRHOV) * t1(2)   + QRight(IMSQRHOW) * t1(3))
         wR = invSqrtRhoR * (QRight(IMSQRHOU) * t2(1)   + QRight(IMSQRHOV) * t2(2)   + QRight(IMSQRHOW) * t2(3))
         pR = QRight(IMP)

         fR(IMC)      = uR*cR
         fR(IMSQRHOU) = 0.5_RP*rhoR*uR*uR + pR
         fR(IMSQRHOV) = 0.5_RP*rhoR*uR*vR
         fR(IMSQRHOW) = 0.5_RP*rhoR*uR*wR
         fR(IMP)      = 0.0_RP

!
!        Perform the average and rotation
!        --------------------------------
         fL = 0.5_RP*(fL + fR)
         fR = fL
!
!        Add the non-conservative term
!        -----------------------------
         fL(IMSQRHOU) = fL(IMSQRHOU) + 0.5_RP*cL*(muR-muL) + 0.25_RP*rhoL*uL*(uR-uL)
         fL(IMSQRHOV) = fL(IMSQRHOV) + 0.25_RP*rhoL*uL*(vR-vL)
         fL(IMSQRHOW) = fL(IMSQRHOW) + 0.25_RP*rhoL*uL*(wR-wL)
         fL(IMP)      = fL(IMP)      + 0.5_RP*invMa2L*(uR-uL)

         fR(IMSQRHOU) = fR(IMSQRHOU) + 0.5_RP*cR*(muL-muR) + 0.25_RP*rhoR*uR*(uL-uR)
         fR(IMSQRHOV) = fR(IMSQRHOV) + 0.25_RP*rhoR*uR*(vL-vR)
         fR(IMSQRHOW) = fR(IMSQRHOW) + 0.25_RP*rhoR*uR*(wL-wR)
         fR(IMP)      = fR(IMP)      + 0.5_RP*invma2R*(uL-uR)

         fL(IMSQRHOU:IMSQRHOW) = nHat*fL(IMSQRHOU) + t1*fL(IMSQRHOV) + t2*fL(IMSQRHOW)
         fR(IMSQRHOU:IMSQRHOW) = nHat*fR(IMSQRHOU) + t1*fR(IMSQRHOV) + t2*fR(IMSQRHOW)

      end subroutine CentralRiemannSolver

      subroutine ExactRiemannSolver(QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL, fR, invMa2L, invMa2R)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: rhoL
         real(kind=RP), intent(in)       :: rhoR
         real(kind=RP), intent(in)       :: muL
         real(kind=RP), intent(in)       :: muR
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: fL(1:NCONS)
         real(kind=RP), intent(out)      :: fR(1:NCONS)
         real(kind=RP), intent(in)       :: invMa2L
         real(kind=RP), intent(in)       :: invMa2R
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: cL,uL, vL, wL, pL, invRhoL, invSqrtRhoL, lambdaMinusL, lambdaPlusL
         real(kind=RP)  :: cR,uR, vR, wR, pR, invRhoR, invSqrtRhoR, lambdaMinusR, lambdaPlusR
         real(kind=RP)  :: rhoStarL, rhoStarR, uStar, pStar, rhoStar, vStar, wStar, cuStar, halfRhouStar
         real(kind=RP)  :: lambda_mu = 0.0_RP

!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         invRhoL     = 1.0_RP / rhoL
         invSqrtRhoL = sqrt(invRhoL)
         cL = QLeft(IMC)
         uL = invSqrtRhoL * (QLeft(IMSQRHOU) * nHat(1) + QLeft(IMSQRHOV) * nHat(2) + QLeft(IMSQRHOW) * nHat(3))
         vL = invSqrtRhoL * (QLeft(IMSQRHOU) * t1(1)   + QLeft(IMSQRHOV) * t1(2)   + QLeft(IMSQRHOW) * t1(3))
         wL = invSqrtRhoL * (QLeft(IMSQRHOU) * t2(1)   + QLeft(IMSQRHOV) * t2(2)   + QLeft(IMSQRHOW) * t2(3))
         pL = QLeft(IMP)

         invRhoR     = 1.0_RP / rhoR
         invSqrtRhoR = sqrt(invRhoR)
         cR = QRight(IMC)
         uR = invSqrtRhoR * (QRight(IMSQRHOU) * nHat(1) + QRight(IMSQRHOV) * nHat(2) + QRight(IMSQRHOW) * nHat(3))
         vR = invSqrtRhoR * (QRight(IMSQRHOU) * t1(1)   + QRight(IMSQRHOV) * t1(2)   + QRight(IMSQRHOW) * t1(3))
         wR = invSqrtRhoR * (QRight(IMSQRHOU) * t2(1)   + QRight(IMSQRHOV) * t2(2)   + QRight(IMSQRHOW) * t2(3))
         pR = QRight(IMP)

!
!        Compute the Star Region
!        -----------------------
         lambdaMinusR = 0.5_RP * (uR - sqrt(uR*uR + 4.0_RP*invMa2R/rhoR))
         lambdaPlusR  = 0.5_RP * (uR + sqrt(uR*uR + 4.0_RP*invMa2R/rhoR))

         lambdaMinusL = 0.5_RP * (uL - sqrt(uL*uL + 4.0_RP*invMa2L/rhoL))
         lambdaPlusL  = 0.5_RP * (uL + sqrt(uL*uL + 4.0_RP*invMa2L/rhoL))

         uStar = (pR-pL+rhoR*uR*lambdaMinusR-rhoL*uL*lambdaPlusL)/(rhoR*lambdaMinusR - rhoL*lambdaPlusL)
         pStar = pR + rhoR*lambdaMinusR*(uR-uStar)
         rhoStarL = (rhoL*lambdaPlusL)/(uStar-lambdaMinusL)
         rhoStarR = (rhoR*lambdaMinusR)/(uStar - lambdaPlusR)

         if ( uStar .ge. 0.0_RP ) then
            rhoStar = rhoStarL
            vStar   = vL
            wStar   = wL

         else
            rhoStar = rhoStarR
            vStar   = vR
            wStar   = wR

         end if

         cuStar = 0.5_RP*(cL*uL + cR*uR)
         halfRhouStar = 0.5_RP*rhoStar*uStar
!
!      - Add first the common (conservative) part
         !fL = [cuStar+lambda_mu*(muL-muR), rhoStar*uStar*uStar + pStar, rhoStar*uStar*vStar, rhoStar*uStar*wStar, dimensionless % invMa2 * uStar]
         fL = [cuStar+lambda_mu*(muL-muR), rhoStar*uStar*uStar + pStar, rhoStar*uStar*vStar, rhoStar*uStar*wStar, 0.0_RP]
         fR = fL
!
!      - Add the non--conservative part
         fL = fL + [0.0_RP, cL*0.5_RP*(muR-muL)-halfRhouStar*uL,-halfRhouStar*vL, -halfRhouStar*wL, (invMa2L)*(uStar-uL)]
         fR = fR + [0.0_RP, cR*0.5_RP*(muL-muR)-halfRhouStar*uR,-halfRhouStar*vR, -halfRhouStar*wR, (invMa2R)*(uStar-uR)]
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         fL(2:4) = nHat*fL(2) + t1*fL(3) + t2*fL(4)
         fR(2:4) = nHat*fR(2) + t1*fR(3) + t2*fR(4)

      end subroutine ExactRiemannSolver

end module RiemannSolvers_MU