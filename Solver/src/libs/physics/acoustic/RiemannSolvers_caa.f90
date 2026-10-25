#include "Includes.h"
module RiemannSolvers_CAAKeywordsModule

     integer,                       parameter :: KEYWORD_LENGTH           = 132
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_SOLVER_NAME_KEY  = "riemann solver"
     character(len=KEYWORD_LENGTH), parameter :: LAMBDA_STABILIZATION_KEY = "lambda stabilization"
!
!    --------------------------
!    Riemann solver definitions
!    --------------------------
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_CENTRAL_NAME    = "central"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_VS_NAME         = "vector-split"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_EXACT_AVG_NAME  = "exact-average"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_LXF_NAME        = "lax-friedrichs"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_EXACT_JUMP_NAME   = "exact-jump"

     enum, bind(C)
        ! enumerator :: RIEMANN_ROE = 1, RIEMANN_LXF, RIEMANN_RUSANOV
        enumerator :: RIEMANN_CENTRAL= 1, RIEMANN_LXF
        enumerator :: RIEMANN_VS, RIEMANN_EXACT_AVG, RIEMANN_EXACT_JUMP
     end enum

end module RiemannSolvers_CAAKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!
module RiemannSolvers_CAA
   use SMConstants
   use Physics_CAA
   use PhysicsStorage_CAA
   use VariableConversion_CAA
   use FluidData_CAA

   implicit none

   private
   public whichRiemannSolver
   public SetRiemannSolver, DescribeRiemannSolver
   public RiemannSolver

   abstract interface
      subroutine RiemannSolverFCN(QLeft, QRight, QbaseL, QbaseR, nHat, t1, t2, flux)
         use SMConstants
         use PhysicsStorage_CAA
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: QbaseL(1:NCONSB)
         real(kind=RP), intent(in)       :: QbaseR(1:NCONSB)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
      end subroutine RiemannSolverFCN
   end interface

   procedure(RiemannSolverFCN),      protected, pointer  :: RiemannSolver      => NULL()

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
         use RiemannSolvers_CAAKeywordsModule
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
!        ---------------------------------------------
!        Choose the Riemann solver (by default is VS)
!        ---------------------------------------------
         if (controlVariables % containsKey(RIEMANN_SOLVER_NAME_KEY)) then

            keyword = controlVariables % stringValueForKey(RIEMANN_SOLVER_NAME_KEY, KEYWORD_LENGTH)
            call toLower(keyword)

            select case (keyword)
            case (RIEMANN_VS_NAME)
               RiemannSolver => VSRiemannSolver
               whichRiemannSolver = RIEMANN_VS

            case (RIEMANN_LXF_NAME)
               RiemannSolver => LxFRiemannSolver
               whichRiemannSolver = RIEMANN_LXF

            case (RIEMANN_EXACT_AVG_NAME)
              RiemannSolver => ExactAverageRiemannSolver
               whichRiemannSolver = RIEMANN_EXACT_AVG

            case (RIEMANN_EXACT_JUMP_NAME)
              RiemannSolver => ExactJumpRiemannSolver
               whichRiemannSolver = RIEMANN_EXACT_JUMP

            case (RIEMANN_CENTRAL_NAME)
               RiemannSolver => CentralRiemannSolver
               whichRiemannSolver = RIEMANN_CENTRAL

            case default
               print*, "Riemann Solver not recognized."
               errorMessage(STD_OUT)
               error stop

            end select

         else
!
!           Select VS by default
!           ---------------------
            RiemannSolver => VSRiemannSolver
            whichRiemannSolver = RIEMANN_VS

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
      end subroutine SetRiemannSolver

      subroutine DescribeRiemannSolver
!
!        -------
!        Modules
!        -------
         use RiemannSolvers_CAAKeywordsModule

         select case (whichRiemannSolver)
         case (RIEMANN_VS)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Vector Splitting"

         case (RIEMANN_LXF)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Lax-Friedrichs"
            write(STD_OUT,'(30X,A,A30,F10.3)') "->","Lambda stabilization: ", lambdaStab

         case (RIEMANN_EXACT_AVG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Exact Average"

         case (RIEMANN_EXACT_JUMP)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Exact Jump"

         case (RIEMANN_CENTRAL)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Central"
            write(STD_OUT,'(30X,A,A30,F10.3)') "->","Lambda stabilization: ", lambdaStab

         end select

      end subroutine DescribeRiemannSolver
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!        Riemann solvers
!        ---------------
!        APE 1/4 are not Rotational Invariant. Current LxF or central should not be used for these equations.
!        ---------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine CentralRiemannSolver(QLeft, QRight, QbaseL, QbaseR, nHat, t1, t2, flux)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: QbaseL(1:NCONSB)
         real(kind=RP), intent(in)       :: QbaseR(1:NCONSB)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: rhoL, uL, vL, wL, pL, a2L
         real(kind=RP) :: rhoR, uR, vR, wR, pR, a2R
         real(kind=RP) :: QLRot(NCONS), QRRot(NCONS) 
         real(kind=RP) :: QBaseLRot(NCONSB), QbaseRRot(NCONSB) 
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(1)
         rhoR = QRight(1)

         uL = QLeft(2) * nHat(1) + QLeft(3) * nHat(2) + QLeft(4) * nHat(3)
         vL = QLeft(2) * t1(1)   + QLeft(3) * t1(2)   + QLeft(4) * t1(3)
         wL = QLeft(2) * t2(1)   + QLeft(3) * t2(2)   + QLeft(4) * t2(3)

         uR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)
         vR = QRight(2) * t1(1)   + QRight(3) * t1(2)   + QRight(4) * t1(3)
         wR = QRight(2) * t2(1)   + QRight(3) * t2(2)   + QRight(4) * t2(3)

         pL = QLeft(5)
         pR = QRight(5)
!
         QLRot = (/ rhoL, uL, vL, wL, pL /)
         QRRot = (/ rhoR, uR, vR, wR, pR /)
!
         rhoL = QbaseL(1)
         rhoR = QbaseR(1)

         uL = QbaseL(2) * nHat(1) + QbaseL(3) * nHat(2) + QbaseL(4) * nHat(3)
         vL = QbaseL(2) * t1(1)   + QbaseL(3) * t1(2)   + QbaseL(4) * t1(3)
         wL = QbaseL(2) * t2(1)   + QbaseL(3) * t2(2)   + QbaseL(4) * t2(3)

         uR = QbaseR(2) * nHat(1) + QbaseR(3) * nHat(2) + QbaseR(4) * nHat(3)
         vR = QbaseR(2) * t1(1)   + QbaseR(3) * t1(2)   + QbaseR(4) * t1(3)
         wR = QbaseR(2) * t2(1)   + QbaseR(3) * t2(2)   + QbaseR(4) * t2(3)

         pL = QbaseL(5)
         pR = QbaseR(5)

         a2L = QbaseL(6)
         a2R = QbaseR(6)
!
         QbaseLRot = (/ rhoL, uL, vL, wL, pL, a2L /)
         QbaseRRot = (/ rhoR, uR, vR, wR, pR, a2R /)
!        Perform the average using the averaging function
!        ------------------------------------------------
         call StandardAverage(QLRot, QRRot, QBaseLRot, QbaseRRot, nHat, flux)
!
!        ************************************************
!        Momentum equations are already in cartesian frame
!        ************************************************
!
         ! flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

      end subroutine CentralRiemannSolver
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LxFRiemannSolver( QLeft, QRight, QbaseL, QbaseR, nHat, t1, t2, flux )
         implicit none
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: QbaseL(1:NCONSB)
         real(kind=RP), intent(in)       :: QbaseR(1:NCONSB)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
         real(kind=RP) :: rhoL, uL, vL, wL, pL, aL
         real(kind=RP) :: rhoR, uR, vR, wR, pR, aR
         real(kind=RP) :: QLRot(NCONS), QRRot(NCONS) 
         real(kind=RP) :: QBaseLRot(NCONSB), QbaseRRot(NCONSB) 
         real(kind=RP)  :: lambda, stab(NCONS)
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QbaseL(1)
         rhoR = QbaseR(1)

         uL = QbaseL(2) * nHat(1) + QbaseL(3) * nHat(2) + QbaseL(4) * nHat(3)
         vL = QbaseL(2) * t1(1)   + QbaseL(3) * t1(2)   + QbaseL(4) * t1(3)
         wL = QbaseL(2) * t2(1)   + QbaseL(3) * t2(2)   + QbaseL(4) * t2(3)

         uR = QbaseR(2) * nHat(1) + QbaseR(3) * nHat(2) + QbaseR(4) * nHat(3)
         vR = QbaseR(2) * t1(1)   + QbaseR(3) * t1(2)   + QbaseR(4) * t1(3)
         wR = QbaseR(2) * t2(1)   + QbaseR(3) * t2(2)   + QbaseR(4) * t2(3)

         pL = QbaseL(5)
         pR = QbaseR(5)
!
         QbaseLRot = (/ rhoL, uL, vL, wL, pL, QBaseL(6) /)
         QbaseRRot = (/ rhoR, uR, vR, wR, pR, QBaseR(6) /)
!
         ! speed of sound of base flow
         aL = sqrt(QbaseL(6))
         aR = sqrt(QbaseR(6))
!
!        Eigenvalues: lambdas = u_base /cdot n +- c
!        -----------
         ! lambda = max(abs(uL), abs(uR)) + max(aL,aR)
         lambda = max(max(abs(uL-aL), abs(uL+aL)), max(abs(uR-aR), abs(uR+aR)))
!
         rhoL = QLeft(1)
         rhoR = QRight(1)

         uL = QLeft(2) * nHat(1) + QLeft(3) * nHat(2) + QLeft(4) * nHat(3)
         vL = QLeft(2) * t1(1)   + QLeft(3) * t1(2)   + QLeft(4) * t1(3)
         wL = QLeft(2) * t2(1)   + QLeft(3) * t2(2)   + QLeft(4) * t2(3)

         uR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)
         vR = QRight(2) * t1(1)   + QRight(3) * t1(2)   + QRight(4) * t1(3)
         wR = QRight(2) * t2(1)   + QRight(3) * t2(2)   + QRight(4) * t2(3)

         pL = QLeft(5)
         pR = QRight(5)
!
         QLRot = (/ rhoL, uL, vL, wL, pL /)
         QRRot = (/ rhoR, uR, vR, wR, pR /)
!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         call StandardAverage(QLRot, QRRot, QBaseLRot, QbaseRRot, nHat, flux)
!
!        Compute the Lax-Friedrichs stabilization
!        ----------------------------------------
         stab = 0.5_RP * lambda * (QRight - QLeft)
!
!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux = flux - lambdaStab * stab
!
!        ************************************************
!        Momentum equations are already in cartesian frame
!        ************************************************
!
         ! flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

      END SUBROUTINE LxFRiemannSolver
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE VSRiemannSolver( QLeft, QRight, QbaseL, QbaseR, nHat, t1, t2, flux )
!
!        **********************************************************************************************
!           This Vector Splitting Riemann solver implementation does not uses the
!           averaging function, nor the lambda stabilization
!           F* /cdot n = A^+_L*QL_L + A^-_R*Q_R
!           Based on the A matrix decomposition on positives and negative eigenvalues
!           Is an upwind Gudunov. If base flow is equal at the two states is an exact RiemannSolver
!        **********************************************************************************************
!
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: QbaseL(1:NCONSB)
         real(kind=RP), intent(in)       :: QbaseR(1:NCONSB)
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
         real(kind=RP) :: rhoL, uL, pL, aL
         real(kind=RP) :: rhoR, uR, pR, aR
         real(kind=RP) :: omega_plus, omega_minus
         real(kind=RP) :: velocity_term

         rhoL = QbaseL(1)
         rhoR = QbaseR(1)

         ! speed of sound of base flow
         aL = sqrt(QbaseL(6))
         aR = sqrt(QbaseR(6))

         uL = QLeft(2) * nHat(1) + QLeft(3) * nHat(2) + QLeft(4) * nHat(3)
         uR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)
         pL = QLeft(5)
         pR = QRight(5)
         
!        Eigenvalues: lambdas = u_base /cdot n +-c, 0,0
!        Riemann Invariants in lambdas not 0
         omega_plus = pL + rhoL*aL*uL
         omega_minus = pR - rhoR*aR*uR

         velocity_term = 0.5_RP * ( omega_plus/rhoL + omega_minus/rhoR )

         flux(1) = 0.0_RP
         flux(2) = nHat(1)*velocity_term
         flux(3) = nHat(2)*velocity_term
         flux(4) = nHat(3)*velocity_term
         flux(5) = 0.5_RP * ( aL*omega_plus - aR*omega_minus)

      END SUBROUTINE VSRiemannSolver

!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExactAverageRiemannSolver( QLeft, QRight, QbaseL, QbaseR, nHat, t1, t2, flux )
!
!        **********************************************************************************************
!           This exact Riemann solver implementation does not uses the
!           flux averaging function, nor the lambda stabilization
!           F* /cdot n = A_avg * Q*
!           Uses a linealization of A, result of the eigenvector expansion of delta F at each discontinuity
!           Is an upwind Gudunov, . If base flow is equal at the two states is an exact RiemannSolver
!        **********************************************************************************************
!
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: QbaseL(1:NCONSB)
         real(kind=RP), intent(in)       :: QbaseR(1:NCONSB)
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
         real(kind=RP) :: rhoL, uL, pL, aL, zL
         real(kind=RP) :: rhoR, uR, pR, aR, zR
         real(kind=RP) :: omega_plus, omega_minus
         real(kind=RP) :: u_star, p_star
         real(kind=RP) :: rhoa_mean, rho_inv_mean
         real(kind=RP) :: velocity_term

         rhoL = QbaseL(1)
         rhoR = QbaseR(1)

         ! speed of sound of base flow
         aL = sqrt(QbaseL(6))
         aR = sqrt(QbaseR(6))

         uL = QLeft(2) * nHat(1) + QLeft(3) * nHat(2) + QLeft(4) * nHat(3)
         uR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)
         pL = QLeft(5)
         pR = QRight(5)
         
!        Eigenvalues: lambdas = u_base /cdot n +-c, 0,0
!        Riemann Invariants in lambdas not 0
         omega_plus = pL + rhoL*aL*uL
         omega_minus = pR - rhoR*aR*uR
         zL = rhoL*aL
         zR = rhoR*aR

         u_star = ( omega_plus - omega_minus) / (zL+zR)
         p_star = ( zR*omega_plus + zL*omega_minus) / (zL+zR)

         ! mean of rho * a
         rhoa_mean = 0.5_RP * (aL*zL+aR*zR)
         ! mean of 1/rho
         rho_inv_mean = 0.5_RP * (1.0_RP/rhoL + 1.0_RP/rhoR)

         velocity_term = rho_inv_mean * p_star

         flux(1) = 0.0_RP
         flux(2) = nHat(1)*velocity_term
         flux(3) = nHat(2)*velocity_term
         flux(4) = nHat(3)*velocity_term
         flux(5) = rhoa_mean * u_star

      END SUBROUTINE ExactAverageRiemannSolver
!
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExactJumpRiemannSolver( QLeft, QRight, QbaseL, QbaseR, nHat, t1, t2, flux )
!
!        **********************************************************************************************
!           This exact Riemann solver implementation does not uses the
!           averaging function, nor the lambda stabilization
!           F* /cdot n is solved as a system of equations of the Rankine-Huginiot jump condition at each
!           eigenvalue discontinuity.
!           Is an upwind Gudunov, . If base flow is equal at the two states is an exact RiemannSolver
!        **********************************************************************************************
!
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: QbaseL(1:NCONSB)
         real(kind=RP), intent(in)       :: QbaseR(1:NCONSB)
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
         real(kind=RP) :: rhoL, uL, pL, aL, zL
         real(kind=RP) :: rhoR, uR, pR, aR, zR
         real(kind=RP) :: omega_plus, omega_minus
         real(kind=RP) :: velocity_term

         rhoL = QbaseL(1)
         rhoR = QbaseR(1)

         ! speed of sound of base flow
         aL = sqrt(QbaseL(6))
         aR = sqrt(QbaseR(6))

         uL = QLeft(2) * nHat(1) + QLeft(3) * nHat(2) + QLeft(4) * nHat(3)
         uR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)
         pL = QLeft(5)
         pR = QRight(5)
         
!        Eigenvalues: lambdas = u_base /cdot n +-c, 0,0
!        Riemann Invariants in lambdas not 0
         omega_plus = pL + rhoL*aL*uL
         omega_minus = pR - rhoR*aR*uR
         zL = rhoL*aL
         zR = rhoR*aR

         velocity_term = ( aL*omega_plus + aR*omega_minus) / (zL+zR)

         flux(1) = 0.0_RP
         flux(2) = nHat(1)*velocity_term
         flux(3) = nHat(2)*velocity_term
         flux(4) = nHat(3)*velocity_term
         flux(5) = ( zR*aL*omega_plus - zL*aR*omega_minus) / (zL+zR)

      END SUBROUTINE ExactJumpRiemannSolver
!
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
!
!////////////////////////////////////////////////////////////////////////////////////////////
      subroutine StandardAverage(QLeft, QRight, QbaseL, QbaseR, nHat, flux)
!
!        *********************************************************************
!           Computes the standard average of the two states:
!              F*_n = {{F_n}} = 0.5 * (FL_n + FR_n)
!
!           State vectors are rotated.
!        *********************************************************************
!
         use Physics_CAA, only: APEFluxNormal
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: QbaseL(1:NCONSB)
         real(kind=RP), intent(in)       :: QbaseR(1:NCONSB)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP), dimension(1:NCONS)     :: fL, fR

!
!        Compute the flux
!        ----------------
         call APEFluxNormal(QRight, fR, QbaseR, nHat)
         call APEFluxNormal(QLeft,  fL, QbaseL, nHat)

         flux = 0.5_RP*(fR+fL)

      end subroutine StandardAverage

end module RiemannSolvers_CAA
