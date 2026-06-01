#include "Includes.h"
module VariableConversion_NS
   use SMConstants
   use PhysicsStorage_NS
   use FluidData_NS
   implicit none

   private
   public   Pressure, Temperature, TemperatureDeriv, PressureDot
   public   get_laminar_mu_kappa, SutherlandsLaw
   public   NSGradientVariables_STATE
   public   NSGradientVariables_ENTROPY
   public   NSGradientVariables_ENERGY
   public   getPrimitiveVariables, getEntropyVariables
   public   getRoeVariables, GetNSViscosity, getVelocityGradients, getTemperatureGradient, getConservativeGradients
   public   getVelocityGradients_STATE
   public   set_getVelocityGradients
  

   interface getTemperatureGradient
      module procedure getTemperatureGradient_0D, getTemperatureGradient_2D, getTemperatureGradient_3D
   end interface
   
   interface getConservativeGradients
      module procedure getConservativeGradients_0D, getConservativeGradients_2D, getConservativeGradients_3D
   end interface

    abstract interface
      pure subroutine getVelocityGradients_f(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         use SMConstants
         use PhysicsStorage_NS
         implicit none
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Q_x(NGRAD), Q_y(NGRAD), Q_z(NGRAD)
         real(kind=RP), intent(out) :: U_x(NDIM), U_y(NDIM), U_z(NDIM)
      end subroutine getVelocityGradients_f
    end interface

    procedure(getVelocityGradients_f), pointer :: getVelocityGradients
   
   contains
!
! /////////////////////////////////////////////////////////////////////
!
!@mark -
!---------------------------------------------------------------------
!! Compute the pressure from the state variables
!---------------------------------------------------------------------
!
      PURE function Pressure(Q) RESULT(P)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(NCONS), INTENT(IN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: P
      
      P = thermodynamics % gammaMinus1*(Q(5) - 0.5_RP*(Q(2)**2 + Q(3)**2 + Q(4)**2)/Q(1))

      end function Pressure
!
! /////////////////////////////////////////////////////////////////////
!---------------------------------------------------------------------
!! Compute the pressure time derivate from the state variables and its time derivatives
!---------------------------------------------------------------------
!
      PURE function PressureDot(Q,QDot) RESULT(PDot)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(NCONS), INTENT(IN) :: Q
      REAL(KIND=RP), DIMENSION(NCONS), INTENT(IN) :: QDot
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: PDot
      
      PDot = QDot(5) + 0.5_RP*(Q(2)**2 + Q(3)**2 + Q(4)**2)/(Q(1)**2)*QDot(1) - dot_product(Q(2:4), QDot(2:4))/Q(1)
      PDot = thermodynamics % gammaMinus1*PDot

      end function PressureDot
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! Compute the temperature from the state variables
!---------------------------------------------------------------------
!
      PURE function Temperature(Q) RESULT(T)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(NCONS), INTENT(IN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: T
!
      T = dimensionless % gammaM2*Pressure(Q)/Q(1)

      end function Temperature
      
      pure function TemperatureDeriv (Q) result (dTdQ)
         implicit none
         !-arguments--------------------------------
         real(kind=RP), intent(in) :: Q(NCONS)
         real(kind=RP)             :: dTdQ(NCONS)
         !-local-variables--------------------------
         real(kind=RP) :: sRho
         !------------------------------------------
         
         sRho = 1._RP / Q(IRHO)
         
         dTdQ = [sRho * (-Q(IRHOE) + sRho * ( Q(IRHOU)**2 + Q(IRHOV)**2 + Q(IRHOW)**2) ), &
                 -Q(IRHOU) * sRho, &
                 -Q(IRHOV) * sRho, &
                 -Q(IRHOW) * sRho, &
                 1._RP]
         
         dTdQ = dTdQ * thermodynamics % gammaMinus1 * dimensionless % gammaM2 * sRho
         
      end function TemperatureDeriv
      
      pure subroutine GetNSViscosity(phi, mu)
         implicit none
         real(kind=RP), intent(in)   :: phi
         real(kind=RP), intent(out)  :: mu

         mu = dimensionless % mu

      end subroutine GetNSViscosity

      pure subroutine get_laminar_mu_kappa(Q,mu,kappa)
         implicit none
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(out) :: mu, kappa
!        
!        ---------------
!        Local variables        
!        ---------------
!
         real(kind=RP)  :: T, suther

         T = Temperature(Q)
         suther = SutherlandsLaw(T)

         mu    = dimensionless % mu * suther
         kappa = mu * dimensionless % mu_to_kappa

      end subroutine get_laminar_mu_kappa

      PURE FUNCTION SutherlandsLaw(T) RESULT(mu)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), INTENT(IN) :: T !! The temperature
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: mu !! The diffusivity
      real(kind=RP) :: tildeT

      tildeT = T*TemperatureReNormalization_Sutherland
!      
      mu = (1._RP + S_div_TRef_Sutherland)/(tildeT + S_div_TRef_Sutherland)*tildeT*SQRT(tildeT)


      END FUNCTION SutherlandsLaw

!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients will be taken.
!---------------------------------------------------------------------
!
      pure subroutine NSGradientVariables_STATE( nEqn, nGrad, Q, U, rho_ )
         implicit none
         integer, intent(in)        :: nEqn, nGrad
         real(kind=RP), intent(in)  :: Q(nEqn)
         real(kind=RP), intent(out) :: U(nGrad)
         real(kind=RP), intent(in), optional :: rho_
!
!        ---------------
!        Local Variables
!        ---------------
!     
         U = Q

      end subroutine NSGradientVariables_STATE

      pure subroutine NSGradientVariables_ENTROPY( nEqn, nGrad, Q, U, rho_ )
         implicit none
         integer, intent(in)        :: nEqn, nGrad
         real(kind=RP), intent(in)  :: Q(nEqn)
         real(kind=RP), intent(out) :: U(nGrad)
         real(kind=RP), intent(in), optional :: rho_
!
!        ---------------
!        Local Variables
!        ---------------
!
         real(kind=RP)  :: invRho, p, invP, rhoV2
            
         invRho = 1.0_RP / Q(IRHO)
         rhoV2 = (POW2(Q(IRHOU))+POW2(Q(IRHOV))+POW2(Q(IRHOW)))*invRho
         p = thermodynamics % gammaMinus1*(Q(IRHOE) - 0.5_RP*rhoV2)
         invP = 1.0_RP / p

         U(IRHO) =   (thermodynamics % gamma-(log(p) - thermodynamics % gamma*log(Q(IRHO))))*thermodynamics % invGammaMinus1 &
                   - 0.5_RP*rhoV2*invP
         U(IRHOU) = Q(IRHOU)*invP
         U(IRHOV) = Q(IRHOV)*invP
         U(IRHOW) = Q(IRHOW)*invP
         U(IRHOE) = -Q(IRHO)*invP


      end subroutine NSGradientVariables_ENTROPY

      pure subroutine NSGradientVariables_ENERGY( nEqn, nGrad, Q, U, rho_ )
         implicit none
         integer, intent(in)        :: nEqn, nGrad
         real(kind=RP), intent(in)  :: Q(nEqn)
         real(kind=RP), intent(out) :: U(nGrad)
         real(kind=RP), intent(in), optional :: rho_
!
!        ---------------
!        Local Variables
!        ---------------
!
         real(kind=RP)  :: invRho, p, rhoV2
            
         invRho = 1.0_RP / Q(IRHO)
         rhoV2 = (POW2(Q(IRHOU))+POW2(Q(IRHOV))+POW2(Q(IRHOW)))*invRho
         p = thermodynamics % gammaMinus1*(Q(IRHOE) - 0.5_RP*rhoV2)

         U(IRHO)  = Q(IRHO)                             ! density (only to have it)
         U(IRHOU) = Q(IRHOU)*invRho                     ! x-velocity
         U(IRHOV) = Q(IRHOV)*invRho                     ! y-velocity
         U(IRHOW) = Q(IRHOW)*invRho                     ! z-velocity
         U(IRHOE) = dimensionless % gammaM2 * p*invRho  ! Temperature

      end subroutine NSGradientVariables_ENERGY
!
! /////////////////////////////////////////////////////////////////////
!
      pure subroutine getPrimitiveVariables(U,V)
!
!        **************************************
!           Primitive variables are:
!              V = [invRho, u,v,w,p,T,a^2]
!        **************************************
!
         implicit none
         real(kind=RP), intent(in)  :: U(NCONS)
         real(kind=RP), intent(out) :: V(NPRIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: invRho

         invRho = 1.0_RP / U(IRHO)

         V(IPIRHO) = invRho
         V(IPU) = U(IRHOU) * invRho
         V(IPV) = U(IRHOV) * invRho
         V(IPW) = U(IRHOW) * invRho
         V(IPP) = thermodynamics % gammaMinus1 * ( U(IRHOE) &
                  - 0.5_RP * (V(IPU)*U(IRHOU) + V(IPV)*U(IRHOV) + V(IPW)*U(IRHOW)))
         V(IPT) = V(IPP) * dimensionless % gammaM2 * invRho 
         V(IPA2) = thermodynamics % gamma * V(IPP) * invRho

      end subroutine getPrimitiveVariables

      pure subroutine getEntropyVariables(U,p,invRho,S)
         implicit none
         real(kind=RP), intent(in)  :: U(NCONS)
         real(kind=RP), intent(in)  :: p, invRho
         real(kind=RP), intent(out) :: S(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: invP
         real(kind=RP)  :: entropy

         invP = 1.0_RP / p

         entropy = log(p * (invRho ** thermodynamics % gamma))
        
         S(1) =   (thermodynamics % gamma - entropy) / (thermodynamics % gammaMinus1) &
                - 0.5_RP * invRho * (POW2(U(IRHOU))+POW2(U(IRHOV))+POW2(U(IRHOW))) * invP

         S(2) = U(IRHOU) * invP
         S(3) = U(IRHOV) * invP
         S(4) = U(IRHOW) * invP
         S(5) = -U(IRHO) * invP

      end subroutine getEntropyVariables

      pure subroutine getRoeVariables(QL, QR, VL, VR, rho, u, v, w, V2, H, a)
!
!        ***************************************************
!           Roe variables are: [rho, u, v, w, H, a]
!        ***************************************************
!           
         implicit none
         real(kind=RP), intent(in)  :: QL(NCONS)
         real(kind=RP), intent(in)  :: QR(NCONS)
         real(kind=RP), intent(in)  :: VL(NPRIM)
         real(kind=RP), intent(in)  :: VR(NPRIM)
         real(kind=RP), intent(out) :: rho, u, v, w, V2, H, a
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: HL, HR
         real(kind=RP)  :: sqrtRhoL, sqrtRhoR
         real(kind=RP)  :: invSumSqrtRhoLR

         associate(gamma => thermodynamics % gamma, &
                   gm1   => thermodynamics % gammaMinus1)

         sqrtRhoL = sqrt(QL(IRHO))  ; sqrtRhoR = sqrt(QR(IRHO))
         invSumSqrtRhoLR = 1.0_RP / (sqrtRhoL + sqrtRhoR)
!
!        Here the enthalpy is defined as rhoH = gogm1 p + 0.5 rho V^2 = p + rhoe
!        -----------------------------------------------------------------------
         HL = (VL(IPP) + QL(IRHOE)) * VL(IPIRHO)
         HR = (VR(IPP) + QR(IRHOE)) * VR(IPIRHO)

         rho = sqrtRhoL * sqrtRhoR
         u   = (sqrtRhoL * VL(IPU) + sqrtRhoR * VR(IPU))*invSumSqrtRhoLR
         v   = (sqrtRhoL * VL(IPV) + sqrtRhoR * VR(IPV))*invSumSqrtRhoLR
         w   = (sqrtRhoL * VL(IPW) + sqrtRhoR * VR(IPW))*invSumSqrtRhoLR
         H   = (sqrtRhoL * HL      + sqrtRhoR * HR     )*invSumSqrtRhoLR
         V2  = POW2(u) + POW2(v) + POW2(w)
         a   = sqrt(gm1*(H - 0.5_RP * V2))

         end associate

      end subroutine getRoeVariables
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------------------
!     Routines to get the velocity gradients
!     --------------------------------------
!
      pure subroutine getVelocityGradients_State(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         implicit none
         !-arguments---------------------------------------------------
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Q_x(NGRAD), Q_y(NGRAD), Q_z(NGRAD)
         real(kind=RP), intent(out) :: U_x(NDIM), U_y(NDIM), U_z(NDIM)
         !-local-variables---------------------------------------------
         real(kind=RP) :: invRho, invRho2, uDivRho(NDIM)
         !-------------------------------------------------------------

         invRho  = 1._RP / Q(IRHO)
         invRho2 = invRho * invRho
         
         uDivRho = [Q(IRHOU) , Q(IRHOV) , Q(IRHOW) ] * invRho2
         
         u_x = invRho * Q_x(IRHOU:IRHOW) - uDivRho * Q_x(IRHO)
         u_y = invRho * Q_y(IRHOU:IRHOW) - uDivRho * Q_y(IRHO)
         u_z = invRho * Q_z(IRHOU:IRHOW) - uDivRho * Q_z(IRHO)

      end subroutine getVelocityGradients_State

      pure subroutine getVelocityGradients_Energy(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         implicit none
         !-arguments---------------------------------------------------
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Q_x(NGRAD), Q_y(NGRAD), Q_z(NGRAD)
         real(kind=RP), intent(out) :: U_x(NDIM), U_y(NDIM), U_z(NDIM)
         !-local-variables---------------------------------------------
         real(kind=RP) :: invRho, invRho2, uDivRho(NDIM)
         !-------------------------------------------------------------

         u_x = Q_x(IRHOU:IRHOW)
         u_y = Q_y(IRHOU:IRHOW)
         u_z = Q_z(IRHOU:IRHOW)

      end subroutine getVelocityGradients_Energy

      pure subroutine getVelocityGradients_Entropy(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         implicit none
         !-arguments---------------------------------------------------
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Q_x(NGRAD), Q_y(NGRAD), Q_z(NGRAD)
         real(kind=RP), intent(out) :: U_x(NDIM), U_y(NDIM), U_z(NDIM)
         !-local-variables---------------------------------------------
         real(kind=RP) :: pDivRho
         real(kind=RP) :: U(NDIM)
         !-------------------------------------------------------------

         pDivRho = Pressure(Q) / Q(IRHO)
         U = Q(IRHOU:IRHOW) / Q(IRHO)

         U_x = pDivRho * Q_x(IRHOU:IRHOW) + U / pDivRho * Q_x(IRHOE)
         U_y = pDivRho * Q_y(IRHOU:IRHOW) + U / pDivRho * Q_y(IRHOE)
         U_z = pDivRho * Q_z(IRHOU:IRHOW) + U / pDivRho * Q_z(IRHOE)

      end subroutine getVelocityGradients_Entropy

!
!/////////////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------
!     Routines to get the temperature gradient
!     -> Currently using the conservative and velocity gradients
!     ----------------------------------------------------------
!
!     ---
!     0D:
!     ---
      pure subroutine getTemperatureGradient_0D(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z,nablaT)
         implicit none
         !-arguments---------------------------------------------------
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Q_x(NGRAD), Q_y(NGRAD), Q_z(NGRAD)
         real(kind=RP), intent(in)  :: U_x(NDIM), U_y(NDIM), U_z(NDIM)
         real(kind=RP), intent(out) :: nablaT(NDIM)
         !-local-variables---------------------------------------------
         real(kind=RP) :: u, v, w, invRho
         !-------------------------------------------------------------
         
         invRho  = 1._RP / Q(IRHO)
         u = Q(IRHOU) / Q(IRHO)
         v = Q(IRHOV) / Q(IRHO)
         w = Q(IRHOW) / Q(IRHO)
         
         nablaT(IX) = thermodynamics % gammaMinus1*dimensionless % gammaM2*(invRho*Q_x(IRHOE) - Q(IRHOE)*invRho*invRho*Q_x(IRHO) - u*u_x(IX)-v*u_x(IY)-w*u_x(IZ))
         nablaT(IY) = thermodynamics % gammaMinus1*dimensionless % gammaM2*(invRho*Q_y(IRHOE) - Q(IRHOE)*invRho*invRho*Q_y(IRHO) - u*u_y(IX)-v*u_y(IY)-w*u_y(IZ))
         nablaT(IZ) = thermodynamics % gammaMinus1*dimensionless % gammaM2*(invRho*Q_z(IRHOE) - Q(IRHOE)*invRho*invRho*Q_z(IRHO) - u*u_z(IX)-v*u_z(IY)-w*u_z(IZ))
      
      end subroutine getTemperatureGradient_0D
!
!/////////////////////////////////////////////////////////////////////////////
!
!     ---
!     2D:
!     ---
      pure subroutine getTemperatureGradient_2D(N,Q,Q_x,Q_y,Q_z,U_x,U_y,U_z,nablaT)
         implicit none
         !-arguments---------------------------------------------------
         integer      , intent(in)  :: N(2)
         real(kind=RP), intent(in)  :: Q  ( NCONS,0:N(1), 0:N(2) )
         real(kind=RP), intent(in)  :: Q_x( NGRAD ,0:N(1), 0:N(2) )
         real(kind=RP), intent(in)  :: Q_y( NGRAD ,0:N(1), 0:N(2) )
         real(kind=RP), intent(in)  :: Q_z( NGRAD ,0:N(1), 0:N(2) )
         real(kind=RP), intent(in)  :: U_x( NDIM ,0:N(1), 0:N(2) )
         real(kind=RP), intent(in)  :: U_y( NDIM ,0:N(1), 0:N(2) )
         real(kind=RP), intent(in)  :: U_z( NDIM ,0:N(1), 0:N(2) )
         real(kind=RP), intent(out) :: nablaT( NDIM ,0:N(1), 0:N(2) )
         !-local-variables---------------------------------------------
         integer :: i,j
         !-------------------------------------------------------------
         
         do j=0, N(2) ; do i=0, N(1)
            call getTemperatureGradient_0D(Q(:,i,j),Q_x(:,i,j),Q_y(:,i,j),Q_z(:,i,j),U_x(:,i,j),U_y(:,i,j),U_z(:,i,j),nablaT(:,i,j))
         end do       ; end do
         
      end subroutine getTemperatureGradient_2D
!
!/////////////////////////////////////////////////////////////////////////////
!
!     ---
!     3D:
!     ---
      pure subroutine getTemperatureGradient_3D(N,Q,Q_x,Q_y,Q_z,U_x,U_y,U_z,nablaT)
         implicit none
         !-arguments---------------------------------------------------
         integer      , intent(in)  :: N(NDIM)
         real(kind=RP), intent(in)  :: Q  ( NCONS ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(in)  :: Q_x( NGRAD ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(in)  :: Q_y( NGRAD ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(in)  :: Q_z( NGRAD ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(out) :: U_x( NDIM  ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(out) :: U_y( NDIM  ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(out) :: U_z( NDIM  ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(out) :: nablaT(NDIM,0:N(1), 0:N(2), 0:N(3) )
         !-local-variables---------------------------------------------
         integer :: i,j,k
         !-------------------------------------------------------------
         
         do k=0, N(3) ; do j=0, N(2) ; do i=0, N(1)
            call getTemperatureGradient_0D(Q(:,i,j,k),Q_x(:,i,j,k),Q_y(:,i,j,k),Q_z(:,i,j,k),U_x(:,i,j,k),U_y(:,i,j,k),U_z(:,i,j,k),nablaT(:,i,j,k))
         end do       ; end do       ; end do
         
      end subroutine getTemperatureGradient_3D
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------------
!     Routines to get the conservative gradients from the primitive gradients
!        ( \nabla \rho must already be in Q_x(1), Q_y(1), Q_z(1) )
!     -----------------------------------------------------------------------
!
!     ---
!     0D:
!     ---
      pure subroutine getConservativeGradients_0D(Q,U_x,U_y,U_z,nablaT,Q_x,Q_y,Q_z)
         implicit none
         !-arguments---------------------------------------------------
         real(kind=RP), intent(in)     :: Q(NCONS)
         real(kind=RP), intent(in)     :: U_x(NDIM), U_y(NDIM), U_z(NDIM)
         real(kind=RP), intent(in)     :: nablaT(NDIM)
         real(kind=RP), intent(inout)  :: Q_x(NGRAD), Q_y(NGRAD), Q_z(NGRAD)
         !-local-variables---------------------------------------------
         real(kind=RP) :: u(NDIM), invRho, cons
         !-------------------------------------------------------------
         
         u = Q(IRHOU:IRHOW) / Q(IRHO)
         invRho  = 1._RP / Q(IRHO)
         cons = Q(IRHO) / (thermodynamics % gammaMinus1*dimensionless % gammaM2)
         
         Q_x(IRHOU:IRHOW) = Q(IRHO) * U_x(1:NDIM) + u(1:NDIM) * Q_x(IRHO)
         Q_y(IRHOU:IRHOW) = Q(IRHO) * U_y(1:NDIM) + u(1:NDIM) * Q_y(IRHO)
         Q_z(IRHOU:IRHOW) = Q(IRHO) * U_z(1:NDIM) + u(1:NDIM) * Q_z(IRHO)
         
         Q_x(IRHOE) = cons * nablaT(IX) + Q(IRHOE) * invRho * Q_x(IRHO) + Q(IRHOU) * U_x(IX) + Q(IRHOV) * U_x(IY) + Q(IRHOV) * U_x(IZ)
         Q_y(IRHOE) = cons * nablaT(IY) + Q(IRHOE) * invRho * Q_y(IRHO) + Q(IRHOU) * U_y(IX) + Q(IRHOV) * U_y(IY) + Q(IRHOV) * U_y(IZ)
         Q_z(IRHOE) = cons * nablaT(IZ) + Q(IRHOE) * invRho * Q_z(IRHO) + Q(IRHOU) * U_z(IX) + Q(IRHOV) * U_z(IY) + Q(IRHOV) * U_z(IZ)
         
      end subroutine getConservativeGradients_0D
!     ---
!     2D:
!     ---
      pure subroutine getConservativeGradients_2D(N,Q,U_x,U_y,U_z,nablaT,Q_x,Q_y,Q_z)
         implicit none
         !-arguments---------------------------------------------------
         integer      , intent(in)     :: N     (2)
         real(kind=RP), intent(in)     :: Q     (NCONS, 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: U_x   (NDIM , 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: U_y   (NDIM , 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: U_z   (NDIM , 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: nablaT(NDIM , 0:N(1), 0:N(2))
         real(kind=RP), intent(inout)  :: Q_x   (NGRAD, 0:N(1), 0:N(2))
         real(kind=RP), intent(inout)  :: Q_y   (NGRAD, 0:N(1), 0:N(2))
         real(kind=RP), intent(inout)  :: Q_z   (NGRAD, 0:N(1), 0:N(2))
         !-local-variables---------------------------------------------
         integer :: i,j
         !-------------------------------------------------------------
         
         do j=0, N(2) ; do i=0, N(1)
            call getConservativeGradients_0D(Q(:,i,j),U_x(:,i,j),U_y(:,i,j),U_z(:,i,j),nablaT(:,i,j),Q_x(:,i,j),Q_y(:,i,j),Q_z(:,i,j))
         end do       ; end do
         
      end subroutine getConservativeGradients_2D
!     ---
!     3D:
!     ---
      pure subroutine getConservativeGradients_3D(N,Q,U_x,U_y,U_z,nablaT,Q_x,Q_y,Q_z)
         implicit none
         !-arguments---------------------------------------------------
         integer      , intent(in)     :: N     (3)
         real(kind=RP), intent(in)     :: Q     (NCONS, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: U_x   (NDIM , 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: U_y   (NDIM , 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: U_z   (NDIM , 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: nablaT(NDIM , 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(inout)  :: Q_x   (NGRAD, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(inout)  :: Q_y   (NGRAD, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(inout)  :: Q_z   (NGRAD, 0:N(1), 0:N(2), 0:N(3))
         !-local-variables---------------------------------------------
         integer :: i,j, k
         !-------------------------------------------------------------
         
         do k=0, N(3) ; do j=0, N(2) ; do i=0, N(1)
            call getConservativeGradients_0D(Q(:,i,j,k),U_x(:,i,j,k),U_y(:,i,j,k),U_z(:,i,j,k),nablaT(:,i,j,k),Q_x(:,i,j,k),Q_y(:,i,j,k),Q_z(:,i,j,k))
         end do       ; end do       ; end do
         
      end subroutine getConservativeGradients_3D

      subroutine set_GetVelocityGradients(grad_vars_)
         implicit none
         integer, intent(in)  :: grad_vars_

         select case (grad_vars_)
         case(GRADVARS_STATE)
            getVelocityGradients => getVelocityGradients_STATE

         case(GRADVARS_ENERGY)
            getVelocityGradients => getVelocityGradients_ENERGY

         case(GRADVARS_ENTROPY)
            getVelocityGradients => getVelocityGradients_ENTROPY

         end select

      end subroutine set_getVelocityGradients
end module VariableConversion_NS