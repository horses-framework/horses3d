#include "Includes.h"
module VariableConversion
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS)) && (!(TRANSPORT))
   use VariableConversion_NS
#elif defined(NAVIERSTOKES) && (TRANSPORT)
   use VariableConversion_NSTPT
#elif defined(SPALARTALMARAS)
   USE VariableConversion_NSSA
#elif defined(INCNS)
   use VariableConversion_iNS
#elif defined(MULTIPHASE)
   use VariableConversion_MU
#endif
#if defined(CAHNHILLIARD)
   use VariableConversion_CH
#endif
#if defined(ACOUSTIC)
   use VariableConversion_CAA
#endif
   implicit none

   abstract interface
      subroutine GetGradientValues_f(nEqn, nGradEqn, Q, U, rho_)
         use SMConstants, only: RP
         implicit none
         integer, intent(in)                 :: nEqn, nGradEqn
         real(kind=RP), intent(in)           :: Q(nEqn)
         real(kind=RP), intent(out)          :: U(nGradEqn)
         real(kind=RP), intent(in), optional :: rho_
      end subroutine GetGradientValues_f
   end interface

   contains

#if (defined(CAHNHILLIARD) && defined(NAVIERSTOKES))
      pure subroutine GetNSCHViscosity(phi, mu)
         use SMConstants, only: RP
         use FluidData
         implicit none
         real(kind=RP), intent(in)     :: phi
         real(kind=RP), intent(out)    :: mu
!
!        ---------------
!        Local variables         
!        ---------------
!
         real(kind=RP)  :: cIn01, p

         cIn01 = 0.5_RP * (phi + 1.0_RP)
         p = POW3(cIn01) * (6.0_RP * POW2(cIn01) - 15.0_RP * cIn01 + 10.0_RP)

         mu = dimensionless % mu * ( (1.0_RP - p) + (p)*multiphase % viscRatio)

      end subroutine GetNSCHViscosity
#endif

#if (defined(INCNS) && defined(CAHNHILLIARD))
      pure subroutine GetiNSCHViscosity(phi, mu)
         use SMConstants, only: RP
         use FluidData
         implicit none
         real(kind=RP), intent(in)     :: phi
         real(kind=RP), intent(out)    :: mu
!
!        ---------------
!        Local variables         
!        ---------------
!
         real(kind=RP)  :: cIn01, p

         cIn01 = 0.5_RP * (phi + 1.0_RP)
         p = POW3(cIn01) * (6.0_RP * POW2(cIn01) - 15.0_RP * cIn01 + 10.0_RP)

         mu = dimensionless % mu(1) * (1.0_RP - p) + (p)*dimensionless % mu(2)
      
      end subroutine GetiNSCHViscosity
#endif

end module VariableConversion
