#include "Includes.h"
module ArtificialViscosity
   use SMConstants
   use PhysicsStorage_NS
   use VariableConversion_NS, only: Pressure, getVelocityGradients
   use FluidData_NS, only: thermodynamics


   private
   public ComputeShockSensor

   real(kind=RP), parameter :: sBeta0 = 0.01_RP
   real(kind=RP), parameter :: eps_omega = 1.0e-10_RP

   contains
      subroutine ComputeShockSensor(Q, Q_x, Q_y, Q_z, h, sBeta)
         implicit none
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Q_x(NGRAD)
         real(kind=RP), intent(in)  :: Q_y(NGRAD)
         real(kind=RP), intent(in)  :: Q_z(NGRAD)
         real(kind=RP), intent(in)  :: h 
         real(kind=RP), intent(out) :: sBeta
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: sBetaMax, sTheta, sOmega, c, p, divV, rotV(3)
         real(kind=RP) :: U_x(NDIM), U_y(NDIM), U_z(NDIM)
         
         call getVelocityGradients(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         
         sBetaMax = 2.0_RP / sqrt(POW2(thermodynamics % gamma) - 1.0_RP)

         p = Pressure(Q)
         c = sqrt(thermodynamics % gamma * p / Q(IRHO))

         divV = U_x(IX) + U_y(IY) + U_z(IZ)

         rotV(1) = U_y(IZ) - U_z(IY)
         rotV(2) = U_z(IX) - U_x(IZ)
         rotV(3) = U_x(IY) - U_y(IX)

         sTheta = -h * divV / c
         sOmega = (divV*divV)/(divV*divV + sum(rotV*rotV) + eps_omega)

         sBeta = sTheta * sOmega

         call SmoothedLimiter(sBeta, sBeta0, sBetaMax)

         sBeta = sBeta * sqrt(sum(POW2(Q(IRHOU:IRHOW))/(POW2(Q(IRHO)))) + c*c)
         
      end subroutine ComputeShockSensor

      subroutine SmoothedLimiter(s, s0, sMax)
         implicit none
         real(kind=RP), intent(inout) :: s
         real(kind=RP), intent(in)    :: s0
         real(kind=RP), intent(in)    :: sMax

         s = lmin(lmax(s-s0) - sMax) + sMax

      end subroutine SmoothedLimiter

      pure function Lmax(s)
         implicit none
         real(kind=RP), intent(in) :: s
         real(kind=RP)             :: Lmax
         real(kind=RP), parameter  :: b = 100.0_RP

         Lmax = (s*atan(b*s) - atan(b))/PI + 0.5_RP*(s+1)

      end function

      pure function Lmin(s)
         implicit none
         real(kind=RP), intent(in) :: s
         real(kind=RP)             :: Lmin

         Lmin = s - Lmax(s)
      
      end function


end module ArtificialViscosity