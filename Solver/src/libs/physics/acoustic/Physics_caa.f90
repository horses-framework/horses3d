!
!//////////////////////////////////////////////////////
!
!      CAA APE physics.
!      Modified from DSEM Code
!
!!     The variable mappings for the APE-4 Equations are
!!
!!              Q(1) = rho'
!!              Q(2) = u'
!!              Q(3) = v'
!!              Q(4) = w'
!!              Q(5) = p'
!
!////////////////////////////////////////////////////////////////////////
!    
!@mark -
!
#include "Includes.h"
!  **************
   module Physics_CAA
!  **************
!
      use SMConstants
      use PhysicsStorage_CAA
      use VariableConversion_CAA
      use FluidData_CAA
      use Utilities, only: outer_product
      implicit none

      private
      public  APEFlux, APEFlux1D, APEFluxNormal
!
!     ========
      CONTAINS 
!     ========
!
!     
!
!//////////////////////////////////////////////////////////////////////////////
!
!           INVISCID FLUXES
!           ---------------   
!
!//////////////////////////////////////////////////////////////////////////////
!
!------------------------------------------------------------------------
!! Fluxes for APE I and IV. The density is not really solved, is computed
!! from pressure, this is why we have 0 flux.
!! for pressure: F = c_0^2*rho_0*\vec{u'} + \vec{u_0}p'
!! for velocities: F = \vec{u_0} /dot /vec{u'} + p'/rho_0
!------------------------------------------------------------------------
!     
      pure subroutine APEFlux(Q, F, rho_, Qbase)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS , 1:NDIM)
         real(kind=RP), intent(in), optional :: rho_
         real(kind=RP), intent(in), optional :: Qbase(1:NCONSB)
!
!        ---------------
!        Local variables
!        ---------------
         real(kind=RP)              :: velocity_term

         velocity_term = Qbase(IBU)*Q(ICAAU) + Qbase(IBV)*Q(ICAAV) + Qbase(IBW)*Q(ICAAW) + Q(ICAAP)/Qbase(IBRHO)
!
!        X-Flux
!        ------         
         F(ICAARHO, IX ) = 0.0_RP
         F(ICAAU, IX ) = velocity_term
         F(ICAAV, IX ) = 0.0_RP
         F(ICAAW, IX ) = 0.0_RP
         F(ICAAP, IX ) = Qbase(IBA2) * Qbase(IBRHO) * Q(ICAAU) + Qbase(IBU) * Q(ICAAP)
!
!        Y-Flux
!        ------
         F(ICAARHO, IY ) = 0.0_RP
         F(ICAAU, IY ) = 0.0_RP
         F(ICAAV, IY ) = velocity_term
         F(ICAAW, IY ) = 0.0_RP
         F(ICAAP, IY ) = Qbase(IBA2) * Qbase(IBRHO) * Q(ICAAV) + Qbase(IBV) * Q(ICAAP)
!
!        Z-Flux
!        ------
         F(ICAARHO, IZ ) = 0.0_RP
         F(ICAAU, IZ ) = 0.0_RP
         F(ICAAV, IZ ) = 0.0_RP
         F(ICAAW, IZ ) = velocity_term
         F(ICAAP, IZ ) = Qbase(IBA2) * Qbase(IBRHO) * Q(ICAAW) + Qbase(IBW) * Q(ICAAP)

      end subroutine APEFlux
!     
      pure subroutine APEFlux1D(Q, F, Qbase)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS)
         real(kind=RP), intent(in), optional :: Qbase(1:NCONSB)
!
!        ---------------
!        Local variables
!        ---------------

!
!        only X-Flux
!        ------         
         F(ICAARHO) = 0.0_RP
         F(ICAAU) = Qbase(IBU)*Q(ICAAU) + Qbase(IBV)*Q(ICAAV) + Qbase(IBW)*Q(ICAAW) + Q(ICAAP)/Qbase(IBRHO)
         F(ICAAV) = 0.0_RP
         F(ICAAW) = 0.0_RP
         F(ICAAP) = Qbase(IBA2) * Qbase(IBRHO) * Q(ICAAU) + Qbase(IBU) * Q(ICAAP)

      end subroutine APEFlux1D
!     
      pure subroutine APEFluxNormal(Q, F, Qbase, nHat)
         ! Computes F_n = /vec{F} /cdot /vec{n}
         ! Q, Qbase are rotated in normal, tan1, tan2
         ! F_n is cartesian
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS)
         real(kind=RP), intent(in)   :: Qbase(1:NCONSB)
         real(kind=RP), intent(in)   :: nHat(1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
         real(kind=RP)              :: velocity_term
         real(kind=RP)              :: u_n, uBase_n

!        The dot product is mantained in the two basis, as is only a rotation
!        vec{u_0}_cartesian /dot /vec{u'}_cartesian =  vec{u_0}_rotated /dot /vec{u'}_rotated  
         velocity_term = Qbase(IBU)*Q(ICAAU) + Qbase(IBV)*Q(ICAAV) + Qbase(IBW)*Q(ICAAW) + Q(ICAAP)/Qbase(IBRHO)
!
!        Flux dot n
!        ------         
         F(ICAARHO) = 0.0_RP
         F(ICAAU) = nHat(1)*velocity_term
         F(ICAAV) = nHat(2)*velocity_term
         F(ICAAW) = nHat(3)*velocity_term
         F(ICAAP) = Qbase(IBA2) * Qbase(IBRHO) * Q(ICAAU) + Qbase(IBU) * Q(ICAAP)

      end subroutine APEFluxNormal
   END Module Physics_CAA
!
! /////////////////////////////////////////////////////////////////////
!
!----------------------------------------------------------------------
!! This routine returns the maximum eigenvalues for the APE equations 
!! for the given solution value in each spatial direction. 
!! These are to be used to compute the local time step.
!----------------------------------------------------------------------
!
      SUBROUTINE ComputeEigenvaluesForStateCAA( Q, Qbase, eigen )
      
      USE SMConstants
      USE PhysicsStorage_CAA
      use FluidData_CAA,          only: Thermodynamics
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=Rp), DIMENSION(NCONS) :: Q
      REAL(KIND=Rp), DIMENSION(NCONSB) :: Qbase
      REAL(KIND=Rp), DIMENSION(3)     :: eigen
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=Rp) :: u, v, w, p, a
!      
      u = ABS( Q(IBU) )
      v = ABS( Qbase(IBV) )
      w = ABS( Q(IBW) )
      !p = Qbase(IBP)
      a = SQRT(Qbase(IBA2))
      
      eigen(1) = u + a
      eigen(2) = v + a
      eigen(3) = w + a

      END SUBROUTINE ComputeEigenvaluesForStateCAA
