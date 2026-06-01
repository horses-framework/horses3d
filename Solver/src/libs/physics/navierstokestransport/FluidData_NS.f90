#include "Includes.h"
module FluidData_NS
   use SMConstants
   implicit none

   private
   public Thermodynamics_t , thermodynamics , SetThermodynamics
   public Dimensionless_t  , dimensionless  , SetDimensionless
   public RefValues_t      , refValues      , SetRefValues

   public  getThermalConductivity
   public  equationOfState

   integer, parameter   :: STR_LEN_FLUIDDATA = 128
!
!  ----------------
!  Type definitions
!  ----------------
!
   type Thermodynamics_t
      character(len=STR_LEN_FLUIDDATA) :: fluidName   ! Fluid name
      real(kind=RP)                    :: R           ! Gas constant
      real(kind=RP)                    :: gamma       ! Heat ratio
      real(kind=RP)                    :: sqrtGamma
      real(kind=RP)                    :: gammaMinus1
      real(kind=RP)                    :: gammaMinus1Div2
      real(kind=RP)                    :: gammaPlus1Div2
      real(kind=RP)                    :: gammaMinus1Div2sg
      real(kind=RP)                    :: gammaMinus1Div2g
      real(kind=RP)                    :: InvGammaPlus1Div2
      real(kind=RP)                    :: InvGammaMinus1
      real(kind=RP)                    :: InvGamma
      real(kind=RP)                    :: gammaDivGammaMinus1
      real(kind=RP)                    :: cp          ! R * gogm1
      real(kind=RP)                    :: cv          ! R * invgm1
      real(kind=RP)                    :: lambda      ! Bulk viscosity ratio
   end type Thermodynamics_t


   type RefValues_t
      real(kind=RP)        :: V
      real(kind=RP)        :: T
      real(kind=RP)        :: p
      real(kind=RP)        :: rho
      real(kind=RP)        :: mu
      real(kind=RP)        :: kappa
      real(kind=RP)        :: AoATheta
      real(kind=RP)        :: AoAPhi
   end type RefValues_t

   type Dimensionless_t
      real(kind=RP)        :: cp
      real(kind=RP)        :: cv
      real(kind=RP)        :: Re
      real(kind=RP)        :: Pr       ! Prandtl number
      real(kind=RP)        :: Prt      ! Turbulent Prandtl number
      real(kind=RP)        :: Fr
      real(kind=RP)        :: mu
      real(kind=RP)        :: kappa
      real(kind=RP)        :: mu_to_kappa
      real(kind=RP)        :: Mach
      real(kind=RP)        :: gammaM2
      real(kind=RP)        :: invFr2
      real(kind=RP)        :: gravity_dir(NDIM)
   end type Dimensionless_t
!
!  ---------
!  Instances: they are public and can be accessed from the outside.
!  ---------
!
   type(Thermodynamics_t), protected   :: thermodynamics
   type(RefValues_t),      protected   :: refValues
   type(Dimensionless_t),  protected   :: dimensionless


     interface getThermalConductivity
         module procedure getThermalConductivity0D, getThermalConductivity3D
     end interface getThermalConductivity

   contains
!
!///////////////////////////////////////////////////////////////////////
!
!        This Set* subroutines are required since the instances defined
!     here are protected. The protected argument implies that they
!     can not be modified outside this module (these quantities are
!     really important to risk accidental modification).
!
!///////////////////////////////////////////////////////////////////////
!
      subroutine SetThermodynamics( thermodynamics_ )
         implicit none
         type(Thermodynamics_t), intent(in)  :: thermodynamics_

         thermodynamics % fluidName           = trim(thermodynamics_ % fluidName)
         thermodynamics % R                   = thermodynamics_ % R
         thermodynamics % gamma               = thermodynamics_ % gamma
         thermodynamics % sqrtGamma           = thermodynamics_ % sqrtGamma
         thermodynamics % gammaMinus1         = thermodynamics_ % gammaMinus1
         thermodynamics % gammaMinus1Div2     = thermodynamics_ % gammaMinus1Div2
         thermodynamics % gammaPlus1Div2      = thermodynamics_ % gammaPlus1Div2
         thermodynamics % gammaMinus1Div2sg   = thermodynamics_ % gammaMinus1Div2sg
         thermodynamics % gammaMinus1Div2g    = thermodynamics_ % gammaMinus1Div2g
         thermodynamics % InvGammaPlus1Div2   = thermodynamics_ % InvGammaPlus1Div2
         thermodynamics % InvGammaMinus1      = thermodynamics_ % InvGammaMinus1
         thermodynamics % InvGamma            = thermodynamics_ % InvGamma
         thermodynamics % gammaDivGammaMinus1 = thermodynamics_ % gammaDivGammaMinus1
         thermodynamics % cp                  = thermodynamics_ % cp
         thermodynamics % cv                  = thermodynamics_ % cv
         thermodynamics % lambda              = thermodynamics_ % lambda

      end subroutine SetThermodynamics

      subroutine SetRefValues( refValues_ )
         implicit none
         type(RefValues_t),   intent(in)     :: refValues_

         refValues % T        = refValues_ % T
         refValues % p        = refValues_ % p
         refValues % rho      = refValues_ % rho
         refValues % V        = refValues_ % V
         refValues % mu       = refValues_ % mu
         refValues % kappa    = refValues_ % kappa
         refValues % AoATheta = refValues_ % AoATheta
         refValues % AoAPhi   = refValues_ % AoAPhi

      end subroutine SetRefValues

      subroutine SetDimensionless( dimensionless_ )
         implicit none
         type(Dimensionless_t),  intent(in)     :: dimensionless_

         dimensionless % cp              = dimensionless_ % cp
         dimensionless % cv              = dimensionless_ % cv
         dimensionless % Re              = dimensionless_ % Re
         dimensionless % Pr              = dimensionless_ % Pr
         dimensionless % Prt             = dimensionless_ % Prt
         dimensionless % Fr              = dimensionless_ % Fr
         dimensionless % mu              = dimensionless_ % mu
         dimensionless % kappa           = dimensionless_ % kappa
         dimensionless % mu_to_kappa     = dimensionless_ % mu_to_kappa
         dimensionless % Mach            = dimensionless_ % Mach
         dimensionless % gammaM2         = dimensionless_ % gammaM2
         dimensionless % invFr2 = dimensionless_ % invFr2
         dimensionless % gravity_dir     = dimensionless_ % gravity_dir

      end subroutine SetDimensionless
!
!/////////////////////////////////////////////////////////////////////////////
!
!        Equation of state
!        -----------------
!
!/////////////////////////////////////////////////////////////////////////////
!
      pure subroutine equationOfState(p, rho, T)
         implicit none
         real(kind=RP),    intent(in)  :: p
         real(kind=RP),    intent(in)  :: rho
         real(kind=RP),    intent(out) :: T
   
         T = p * dimensionless % gammaM2 / rho

      end subroutine equationOfState

!
!//////////////////////////////////////////////////////////////////////////////
!
!        Get the thermal conductivity from the viscosity and Prandtl number
!
!//////////////////////////////////////////////////////////////////////////////
!
      pure subroutine getThermalConductivity0D(mu, Pr, kappa)
         implicit none
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: Pr
         real(kind=RP), intent(out) :: kappa

         kappa = mu / (Pr * thermodynamics % gammaMinus1 * POW2(dimensionless % Mach))

      end subroutine getThermalConductivity0D

      pure subroutine getThermalConductivity3D(N, mu, Pr, kappa)
         implicit none
         integer,       intent(in)  :: N(3)
         real(kind=RP), intent(in)  :: mu(0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)  :: Pr
         real(kind=RP), intent(out) :: kappa(0:N(1), 0:N(2), 0:N(3))

         kappa = mu / (Pr * thermodynamics % gammaMinus1 * POW2(dimensionless % Mach))

      end subroutine getThermalConductivity3D

end module FluidData_NS