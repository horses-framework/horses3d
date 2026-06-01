!
!////////////////////////////////////////////////////////////////////////
!
!
!      Manufactured solutions definitions for 2D and 3D cases (Euler and Navier-Stokes)
!
!
!////////////////////////////////////////////////////////////////////////////////////////
MODULE ManufacturedSolutionsNS
   USE SMConstants
   USE PhysicsStorage_NS
   USE Physics_NS
   use FluidData_NS
   IMPLICIT NONE

   private
   public   InitializeManufacturedSol, ManufacturedSolP
   public   ManufacturedSolutionState, ManufacturedSolutionDeriv
   public   ManufacturedSolutionSourceEuler
   public   ManufacturedSolutionSourceNS

   REAL(KIND=RP), DIMENSION(7) :: rC, uC, vC, wC, pC

!========
 CONTAINS
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE InitializeManufacturedSol(ProblemDIM)
      IMPLICIT NONE
!
!     --------------------------------------------------------------------------------------
!     Constant definitions for the manufactured solutions. 
!
!        The primitive variable j(rho, u, v, w, p) is defined as:
!           j = jC(1) + jC(2)*Sin(pi*jC(5)*x) + jC(3)*Sin(pi*jC(6)*y) + jC(4)*Sin(pi*jC(7)*z) 
!     See:
!        |> Roy, Chris, Curt Ober, and Tom Smith. "Verification of a compressible CFD code 
!           using the method of manufactured solutions." 32nd AIAA Fluid Dynamics Conference 
!           and Exhibit. 2002.
!     --------------------------------------------------------------------------------------
!
      CHARACTER(LEN=*) :: ProblemDIM
      
      SELECT CASE (ProblemDIM)
         CASE ('2D')
!
!           -------
!           2D case
!           -------
!
            rC = (/1._RP, 0.15_RP    ,-0.1_RP    , 0._RP, 1._RP      , 0.50_RP    , 0._RP/)
            uC = (/1._RP, 0.0625_RP  ,-0.0375_RP , 0._RP, 1.5_RP     , 0.60_RP    , 0._RP/)
            vC = (/1._RP,-0.09375_RP , 0.05_RP   , 0._RP, 0.5_RP     , 2._RP/3._RP, 0._RP/)
            wC = (/0._RP, 0._RP      , 0._RP     , 0._RP, 0._RP      , 0.00_RP    , 0._RP/)
            pC = (/1._RP, 0.2_RP     , 0.5_RP    , 0._RP, 2.0_RP     , 1.00_RP    , 0._RP/)
         CASE ('3D')
!
!           -------
!           3D case
!           -------
!
            rC = (/1._RP, 0.15_RP    ,-0.10000_RP ,-0.12000_RP, 1.0_RP      , 0.50_RP    , 1.50_RP/)
            uC = (/1._RP, 0.0625_RP  ,-0.03750_RP ,-0.02250_RP, 1.5_RP      , 0.60_RP    , 0.50_RP/)
            vC = (/1._RP,-0.09375_RP , 0.05000_RP , 0.03750_RP, 0.5_RP      , 2._RP/3._RP, 1.25_RP/)
            wC = (/1._RP, 0.01875_RP ,-0.03125_RP , 0.04375_RP, 1.0_RP/3._RP, 1._RP/5._RP, 1.00_RP/)
            pC = (/1._RP, 0.2_RP     , 0.50000_RP ,-0.35000_RP, 2.0_RP      , 1.00_RP    , 1._RP/3._RP/)
!~            rC = (/1._RP,-0.15_RP   ,-0.15000_RP ,-0.15000_RP, 5.0_RP      , 5.0_RP     , 5.00_RP/)
!~            uC = (/0._RP, 0.5_RP    , 0.5_RP     , 0.5_RP    , 1.0_RP      , 1.00_RP    , 1.00_RP/)
!~            vC = (/1._RP,-0.09375_RP , 0.05000_RP , 0.03750_RP, 0.5_RP      , 2._RP/3._RP, 1.25_RP/)
!~            wC = (/1._RP, 0.01875_RP ,-0.03125_RP , 0.04375_RP, 1.0_RP/3._RP, 1._RP/5._RP, 1.00_RP/)
!~            pC = (/1._RP, 0.2_RP     , 0.50000_RP ,-0.35000_RP, 2.0_RP      , 1.00_RP    , 1._RP/3._RP/)
      END SELECT
      
   END SUBROUTINE InitializeManufacturedSol
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
   FUNCTION ManufacturedSolP(x) RESULT(p)
      IMPLICIT NONE
      REAL(KIND=RP) :: x(3), p
      
      p    = pC(1) + pC(2)*Sin(pi*pC(5)*x(1)) + pC(3)*Sin(pi*pC(6)*x(2)) + pC(4)*Sin(pi*pC(7)*x(3)) 
      
   END FUNCTION ManufacturedSolP
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ManufacturedSolutionState( x, t, Q )
      IMPLICIT NONE
!
!     ----------------------------------------------------
!     Manufactured solution state (valid for Euler and NS)
!     ----------------------------------------------------
!      
      
      REAL(KIND=RP) :: x(3), t
      REAL(KIND=RP) :: Q(NCONS)
      
      REAL(KIND=RP) :: rho, u, v, w, p
      
      associate ( gamma => thermodynamics % gamma ) 

      rho  = rC(1) + rC(2)*Sin(pi*rC(5)*x(1)) + rC(3)*Sin(pi*rC(6)*x(2)) + rC(4)*Sin(pi*rC(7)*x(3)) 
      u    = uC(1) + uC(2)*Sin(pi*uC(5)*x(1)) + uC(3)*Sin(pi*uC(6)*x(2)) + uC(4)*Sin(pi*uC(7)*x(3)) 
      v    = vC(1) + vC(2)*Sin(pi*vC(5)*x(1)) + vC(3)*Sin(pi*vC(6)*x(2)) + vC(4)*Sin(pi*vC(7)*x(3)) 
      w    = wC(1) + wC(2)*Sin(pi*wC(5)*x(1)) + wC(3)*Sin(pi*wC(6)*x(2)) + wC(4)*Sin(pi*wC(7)*x(3)) 
      p    = pC(1) + pC(2)*Sin(pi*pC(5)*x(1)) + pC(3)*Sin(pi*pC(6)*x(2)) + pC(4)*Sin(pi*pC(7)*x(3)) 
      
      Q(1) = rho
      Q(2) = rho*u
      Q(3) = rho*v
      Q(4) = rho*w
      Q(5) = p/(gamma - 1.0_RP) + 0.5_RP*rho*(u**2 + v**2 + w**2)

      end associate
      
   END SUBROUTINE ManufacturedSolutionState 
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ManufacturedSolutionDeriv( xx, t, nHat, U_x, U_y, U_z )
!
!     ------------------------------------------------------
!     Manufactured solution derivatives (only needed for NS)
!     ------------------------------------------------------
! 
      IMPLICIT NONE
      
      REAL(KIND=RP) :: xx(3), t
      REAL(KIND=RP) :: nHat(3)
      REAL(KIND=RP), INTENT(INOUT) :: U_x(NGRAD), U_y(NGRAD), U_z(NGRAD)
      
      REAL(KIND=RP) :: x, y, z
      
      associate ( gammaM2 => dimensionless % gammaM2 ) 

      x = xx(1)
      y = xx(2)
      z = xx(3)
      
      ! u (velocity)
      U_x(1) = pi*Cos(pi*x*uC(5))*uC(2)*uC(5)
      U_y(1) = pi*Cos(pi*y*uC(6))*uC(3)*uC(6)
      U_z(1) = pi*Cos(pi*z*uC(7))*uC(4)*uC(7)
      
      ! v (velocity)
      U_x(2) = pi*Cos(pi*x*vC(5))*vC(2)*vC(5)
      U_y(2) = pi*Cos(pi*y*vC(6))*vC(3)*vC(6)
      U_z(2) = pi*Cos(pi*z*vC(7))*vC(4)*vC(7)
      
      ! w (velocity)
      U_x(3) = pi*Cos(pi*x*wC(5))*wC(2)*wC(5)
      U_y(3) = pi*Cos(pi*y*wC(6))*wC(3)*wC(6)
      U_z(3) = pi*Cos(pi*z*wC(7))*wC(4)*wC(7)
      
      ! T (Temperature)
      U_x(4) = -((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
               (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
               pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))**2) +     &
               (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) + rC(4)*Sin(pi*z*rC(7)))
      U_y(4) = -((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
               (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
               pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))**2) +     &
               (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) + rC(4)*Sin(pi*z*rC(7)))
      U_z(4) = -((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
               (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
               pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))**2) +     &
               (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) + rC(4)*Sin(pi*z*rC(7)))

      end associate
      
   END SUBROUTINE ManufacturedSolutionDeriv
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ManufacturedSolutionSourceNS( xx, t, Q  )
!
!     --------------------------------
!     Source term for MS Navier Stokes
!     --------------------------------
!
      IMPLICIT NONE
      
      REAL(KIND=RP) :: xx(3), t
      REAL(KIND=RP) :: Q(NCONS)
      
      REAL(KIND=RP) :: x, y, z
      
      associate ( gamma => thermodynamics % gamma, &
                  Mach => dimensionless % Mach, &
                  Re => dimensionless % Re, &
                  Pr => dimensionless % Pr, &
                  gammaM2 => dimensionless % gammaM2 )

      x = xx(1)
      y = xx(2)
      z = xx(3)
      
!
!     -------------------------
!     Mass equation Source term
!     -------------------------
!
      Q(1) = Q(1) + pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))    &
           + pi*Cos(pi*x*uC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*uC(2)*uC(5) +     &
          pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))    &
           + pi*Cos(pi*y*vC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*vC(3)*vC(6) +     &
          pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))    &
           + pi*Cos(pi*z*wC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*wC(4)*wC(7) 

!
!     -------------------------------
!     x-momentum equation Source term
!     -------------------------------
!
      Q(2) = Q(2) + pi*Cos(pi*x*pC(5))*pC(2)*pC(5) +     &
          pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
              Sin(pi*z*uC(7))*uC(4))**2 +     &
          2*pi*Cos(pi*x*uC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*uC(2)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*    &
           uC(5) + pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))     &
           + pi*Cos(pi*y*uC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*uC(3)*uC(6)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))     &
           + pi*Cos(pi*y*vC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*(uC(1) + Sin(pi*x*uC(5))*uC(2) +     &
             Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*vC(3)*vC(6) +     &
          pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))     &
           + pi*Cos(pi*z*uC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*uC(4)*uC(7)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))     &
           + pi*Cos(pi*z*wC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*(uC(1) + Sin(pi*x*uC(5))*uC(2) +     &
             Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*wC(4)*wC(7) -     &
          ((-2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*uC(5))*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*uC(2)*uC(5))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) -     &
             (2*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*Cos(pi*x*rC(5))*Cos(pi*x*uC(5))*rC(2)*rC(5)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*uC(2)*uC(5))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (2*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*Cos(pi*x*pC(5))*Cos(pi*x*uC(5))*pC(2)*pC(5)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*uC(2)*uC(5))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*uC(5))*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*uC(2)*uC(5))/    &
              (Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (4*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*Sin(pi*x*uC(5))*uC(2)*uC(5)**2)/    &
              (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*Sin(pi*y*uC(6))*uC(3)*uC(6)**2)/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*Sin(pi*z*uC(7))*uC(4)*uC(7)**2)/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) -     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*pC(6))*pC(3)*pC(6)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) -     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*pC(7))*pC(4)*pC(7)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (2*gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                  pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
              (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) +     &
             (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                  pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
              (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*pC(5))*pC(2)*pC(5)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                  pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
              (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                  pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
              (3.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))))/Re

!
!     -------------------------------
!     y-momentum equation Source term
!     -------------------------------
!
      Q(3) = Q(3) + pi*Cos(pi*y*pC(6))*pC(3)*pC(6) +     &
       pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
        (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*    &
        (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))     &
        + pi*Cos(pi*x*uC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
          rC(4)*Sin(pi*z*rC(7)))*uC(2)*uC(5)*    &
        (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))     &
        + pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
        (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
           Sin(pi*z*vC(7))*vC(4))**2 +     &
       pi*Cos(pi*x*vC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
          rC(4)*Sin(pi*z*rC(7)))*(uC(1) + Sin(pi*x*uC(5))*uC(2) +     &
          Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*vC(2)*vC(5) +     &
       2*pi*Cos(pi*y*vC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
          rC(4)*Sin(pi*z*rC(7)))*vC(3)*    &
        (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))*    &
        vC(6) + pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
        (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))*    &
        (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))     &
        + pi*Cos(pi*z*vC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
          rC(4)*Sin(pi*z*rC(7)))*vC(4)*vC(7)*    &
        (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))     &
        + pi*Cos(pi*z*wC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
          rC(4)*Sin(pi*z*rC(7)))*(vC(1) + Sin(pi*x*vC(5))*vC(2) +     &
          Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))*wC(4)*wC(7) -     &
       (-((gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
               (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                 pC(4)*Sin(pi*z*pC(7)))*    &
               Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*Sin(pi*x*vC(5))*vC(2)*vC(5)**2)/    &
             ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7)))*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                      pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7)))))) -     &
          (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
               pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2) +     &
               (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))*    &
             (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
           ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                     pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))**2) -     &
          (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
             (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
               pC(4)*Sin(pi*z*pC(7)))*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
           ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                rC(4)*Sin(pi*z*rC(7)))**2*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))) +     &
          (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*pC(5))*pC(2)*pC(5)*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
           ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))) +     &
          (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
               pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
             (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2) +     &
               (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))*    &
             (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
           (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))) -     &
          (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*vC(6))*    &
             (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
               pC(4)*Sin(pi*z*pC(7)))*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2) +     &
               (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))*vC(3)*vC(6))/    &
           ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                     pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))**2) -     &
          (2*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*Cos(pi*y*rC(6))*Cos(pi*y*vC(6))*rC(3)*rC(6)*    &
             (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
               pC(4)*Sin(pi*z*pC(7)))*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*vC(3)*vC(6))/    &
           ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                rC(4)*Sin(pi*z*rC(7)))**2*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))) +     &
          (2*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*Cos(pi*y*pC(6))*Cos(pi*y*vC(6))*pC(3)*pC(6)*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*vC(3)*vC(6))/    &
           ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))) +     &
          (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*vC(6))*    &
             (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
               pC(4)*Sin(pi*z*pC(7)))*    &
             (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2) +     &
               (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))*vC(3)*vC(6))/    &
           (Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))) -     &
          (4*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
             (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
               pC(4)*Sin(pi*z*pC(7)))*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*Sin(pi*y*vC(6))*vC(3)*vC(6)**2)/    &
           (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))) -     &
          (gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
             (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
               pC(4)*Sin(pi*z*pC(7)))*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*Sin(pi*z*vC(7))*vC(4)*vC(7)**2)/    &
           ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))) -     &
          (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
               pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2) +     &
               (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))*    &
             (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
           ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                     pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))**2) -     &
          (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
             (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
               pC(4)*Sin(pi*z*pC(7)))*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
           ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                rC(4)*Sin(pi*z*rC(7)))**2*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))) +     &
          (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*pC(7))*pC(4)*pC(7)*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
           ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))) +     &
          (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
               pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
             (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2) +     &
               (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))*    &
             (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
           (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))) +     &
          (2*gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
               pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2) +     &
               (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))*    &
             (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
               pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
           (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                     pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))**2) +     &
          (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
             (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
               pC(4)*Sin(pi*z*pC(7)))*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
               pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
           (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                rC(4)*Sin(pi*z*rC(7)))**2*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))) -     &
          (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*pC(6))*pC(3)*pC(6)*    &
             Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
               pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
           (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))) -     &
          (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
               pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
             (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2) +     &
               (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))*    &
             (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
               pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
           (3.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                   pC(4)*Sin(pi*z*pC(7))))/    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7))))*    &
             (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
               rC(4)*Sin(pi*z*rC(7)))*    &
             (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                    pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7))))))/Re    

                
!
!     -------------------------------
!     z-momentum equation Source term
!        (2D case!)
!     -------------------------------
!
      Q(4) = Q(4) + pi*Cos(pi*z*pC(7))*pC(4)*pC(7) +     &
          pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))     &
           + pi*Cos(pi*x*uC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*uC(2)*uC(5)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))     &
           + pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))     &
           + pi*Cos(pi*y*vC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*vC(3)*vC(6)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))     &
           + pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
              Sin(pi*z*wC(7))*wC(4))**2 +     &
          pi*Cos(pi*x*wC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*(uC(1) + Sin(pi*x*uC(5))*uC(2) +     &
             Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*wC(2)*wC(5) +     &
          pi*Cos(pi*y*wC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*(vC(1) + Sin(pi*x*vC(5))*vC(2) +     &
             Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))*wC(3)*wC(6) +     &
          2*pi*Cos(pi*z*wC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*wC(4)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))*    &
           wC(7) - (-((gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                  (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                    pC(4)*Sin(pi*z*pC(7)))*    &
                  Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                        pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))*Sin(pi*x*wC(5))*wC(2)*wC(5)**2)/    &
                ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7)))*    &
                  (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                         pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7)))))) -     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) -     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*pC(5))*pC(2)*pC(5)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*Sin(pi*y*wC(6))*wC(3)*wC(6)**2)/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) -     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*pC(6))*pC(3)*pC(6)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*wC(7))*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*wC(4)*wC(7))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) -     &
             (2*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*Cos(pi*z*rC(7))*Cos(pi*z*wC(7))*rC(4)*rC(7)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*wC(4)*wC(7))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (2*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*Cos(pi*z*pC(7))*Cos(pi*z*wC(7))*pC(4)*pC(7)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*wC(4)*wC(7))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*wC(7))*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*wC(4)*wC(7))/    &
              (Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (4*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*Sin(pi*z*wC(7))*wC(4)*wC(7)**2)/    &
              (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (2*gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                  pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
              (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) +     &
             (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                  pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
              (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*pC(7))*pC(4)*pC(7)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                  pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
              (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                  pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
              (3.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))))/Re 
       
!
!     ---------------------------
!     Energy equation Source term
!     ---------------------------
!
      Q(5) = Q(5) -((-((gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                   pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                 Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                       pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                 (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                         (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                           pC(4)*Sin(pi*z*pC(7))))/    &
                       (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                          rC(4)*Sin(pi*z*rC(7)))**2) +     &
                    (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))**2)/    &
               ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))*    &
                 (S_div_TRef_Sutherland + (gammaM2*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))**2)) -     &
            (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                 pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
               Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))**2)/    &
             ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7)))*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))**2) -     &
            (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                 pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
               Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))**2)/    &
             ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7)))*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))**2) -     &
            (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
               (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                 pC(4)*Sin(pi*z*pC(7)))*    &
               Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                      (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                        pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7)))**2) +     &
                 (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7)))))/    &
             ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))**2*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                      pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))) +     &
            (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*pC(5))*pC(2)*pC(5)*    &
               Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                      (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                        pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7)))**2) +     &
                 (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7)))))/    &
             ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7)))*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                      pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))) +     &
            (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                 pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
               (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))**2)/    &
             (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7)))*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                      pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))) -     &
            (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
               (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                 pC(4)*Sin(pi*z*pC(7)))*    &
               Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                      (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                        pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7)))**2) +     &
                 (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7)))))/    &
             ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))**2*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                      pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))) +     &
            (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*pC(6))*pC(3)*pC(6)*    &
               Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                      (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                        pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7)))**2) +     &
                 (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7)))))/    &
             ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7)))*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                      pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))) +     &
            (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                 pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
               (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))**2)/    &
             (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7)))*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                      pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))) -     &
            (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
               (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                 pC(4)*Sin(pi*z*pC(7)))*    &
               Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                      (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                        pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7)))**2) +     &
                 (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7)))))/    &
             ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))**2*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                      pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))) +     &
            (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*pC(7))*pC(4)*pC(7)*    &
               Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                      (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                        pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7)))**2) +     &
                 (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7)))))/    &
             ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7)))*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                      pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))) +     &
            (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                 pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
               (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))**2)/    &
             (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7)))*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                      pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))) +     &
            (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                 pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
               Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               ((2*gammaM2*pi**2*Cos(pi*x*rC(5))**2*rC(2)**2*rC(5)**2*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**3 -     &
                 (2*gammaM2*pi**2*Cos(pi*x*pC(5))*Cos(pi*x*rC(5))*pC(2)*pC(5)*rC(2)*    &
                    rC(5))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2 +     &
                 (gammaM2*pi**2*rC(2)*rC(5)**2*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7)))*Sin(pi*x*rC(5)))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2 -     &
                 (gammaM2*pi**2*pC(2)*pC(5)**2*Sin(pi*x*pC(5)))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7)))))/    &
             ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7)))*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                      pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))) +     &
            (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                 pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
               Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               ((2*gammaM2*pi**2*Cos(pi*y*rC(6))**2*rC(3)**2*rC(6)**2*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**3 -     &
                 (2*gammaM2*pi**2*Cos(pi*y*pC(6))*Cos(pi*y*rC(6))*pC(3)*pC(6)*rC(3)*    &
                    rC(6))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2 +     &
                 (gammaM2*pi**2*rC(3)*rC(6)**2*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7)))*Sin(pi*y*rC(6)))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2 -     &
                 (gammaM2*pi**2*pC(3)*pC(6)**2*Sin(pi*y*pC(6)))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7)))))/    &
             ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7)))*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                      pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))) +     &
            (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                 pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
               Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7))))/    &
                 (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7))))*    &
               ((2*gammaM2*pi**2*Cos(pi*z*rC(7))**2*rC(4)**2*rC(7)**2*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**3 -     &
                 (2*gammaM2*pi**2*Cos(pi*z*pC(7))*Cos(pi*z*rC(7))*pC(4)*pC(7)*rC(4)*    &
                    rC(7))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2 +     &
                 (gammaM2*pi**2*rC(4)*rC(7)**2*    &
                    (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7)))*Sin(pi*z*rC(7)))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2 -     &
                 (gammaM2*pi**2*pC(4)*pC(7)**2*Sin(pi*z*pC(7)))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7)))))/    &
             ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                 rC(4)*Sin(pi*z*rC(7)))*    &
               (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                      pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))))/((-1 + gamma)*Mach**2*Pr*Re))
                    
      Q(5) = Q(5) +     &                                            ! The statement is separated here for ifort not to crash
          pi*Cos(pi*x*pC(5))*pC(2)*pC(5)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))     &
           + pi*Cos(pi*x*uC(5))*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
             pC(4)*Sin(pi*z*pC(7)))*uC(2)*uC(5) +     &
          pi*Cos(pi*y*pC(6))*pC(3)*pC(6)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))     &
           + pi*Cos(pi*y*vC(6))*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
             pC(4)*Sin(pi*z*pC(7)))*vC(3)*vC(6) +     &
          pi*Cos(pi*z*pC(7))*pC(4)*pC(7)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))     &
           + pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*    &
           ((pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                pC(4)*Sin(pi*z*pC(7)))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             ((uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))**2 +     &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))**2 +     &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))**2)/2.) +     &
          pi*Cos(pi*x*uC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*uC(2)*uC(5)*    &
           ((pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                pC(4)*Sin(pi*z*pC(7)))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             ((uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))**2 +     &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))**2 +     &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))**2)/2.) +     &
          pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))*    &
           ((pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                pC(4)*Sin(pi*z*pC(7)))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             ((uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))**2 +     &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))**2 +     &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))**2)/2.) +     &
          pi*Cos(pi*y*vC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*vC(3)*vC(6)*    &
           ((pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                pC(4)*Sin(pi*z*pC(7)))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             ((uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))**2 +     &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))**2 +     &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))**2)/2.) +     &
          pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))*    &
           ((pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                pC(4)*Sin(pi*z*pC(7)))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             ((uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))**2 +     &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))**2 +     &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))**2)/2.) +     &
          (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) + rC(4)*Sin(pi*z*rC(7)))*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*    &
           (-((pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                  (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                    pC(4)*Sin(pi*z*pC(7))))/    &
                ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2)) +     &
             (pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             (2*pi*Cos(pi*x*uC(5))*uC(2)*    &
                 (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))*uC(5) +     &
                2*pi*Cos(pi*x*vC(5))*vC(2)*    &
                 (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))*vC(5) +     &
                2*pi*Cos(pi*x*wC(5))*wC(2)*    &
                 (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))*wC(5))/2.) +     &
          (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) + rC(4)*Sin(pi*z*rC(7)))*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))*    &
           (-((pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                  (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                    pC(4)*Sin(pi*z*pC(7))))/    &
                ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2)) +     &
             (pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             (2*pi*Cos(pi*y*uC(6))*uC(3)*    &
                 (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))*uC(6) +     &
                2*pi*Cos(pi*y*vC(6))*vC(3)*    &
                 (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))*vC(6) +     &
                2*pi*Cos(pi*y*wC(6))*wC(3)*    &
                 (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))*wC(6))/2.) +     &
          pi*Cos(pi*z*wC(7))*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
             pC(4)*Sin(pi*z*pC(7)))*wC(4)*wC(7) +     &
          pi*Cos(pi*z*wC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*wC(4)*    &
           ((pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                pC(4)*Sin(pi*z*pC(7)))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             ((uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))**2 +     &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))**2 +     &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))**2)/2.)*wC(7) +     &
          (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) + rC(4)*Sin(pi*z*rC(7)))*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))*    &
           (-((pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                  (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                    pC(4)*Sin(pi*z*pC(7))))/    &
                ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2)) +     &
             (pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             (2*pi*Cos(pi*z*uC(7))*uC(4)*    &
                 (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))*uC(7) +     &
                2*pi*Cos(pi*z*vC(7))*vC(4)*    &
                 (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))*vC(7) +     &
                2*pi*Cos(pi*z*wC(7))*wC(4)*    &
                 (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))*wC(7))/2.)
          
      Q(5) = Q(5) -     &                                            ! The statement is separated here for ifort not to crash 
          (-((gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                  (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                    pC(4)*Sin(pi*z*pC(7)))*    &
                  Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                        pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))*Sin(pi*y*uC(6))*uC(3)*    &
                  (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                    Sin(pi*z*uC(7))*uC(4))*uC(6)**2)/    &
                ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7)))*    &
                  (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                         pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7)))))) -     &
             (gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*Sin(pi*z*uC(7))*uC(4)*    &
                (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                  Sin(pi*z*uC(7))*uC(4))*uC(7)**2)/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*Sin(pi*x*vC(5))*vC(2)*    &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                  Sin(pi*z*vC(7))*vC(4))*vC(5)**2)/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                  Sin(pi*z*uC(7))*uC(4))*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) -     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                  Sin(pi*z*uC(7))*uC(4))*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*pC(6))*pC(3)*pC(6)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                  Sin(pi*z*uC(7))*uC(4))*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                  Sin(pi*z*uC(7))*uC(4))*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*uC(6))*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*uC(3)*uC(6)*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                  Sin(pi*z*vC(7))*vC(4))*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) -     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                  Sin(pi*z*vC(7))*vC(4))*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*pC(5))*pC(2)*pC(5)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                  Sin(pi*z*vC(7))*vC(4))*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                  Sin(pi*z*vC(7))*vC(4))*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*vC(5))*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*vC(2)*vC(5)*    &
                (pi*Cos(pi*y*uC(6))*uC(3)*uC(6) + pi*Cos(pi*x*vC(5))*vC(2)*vC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*Sin(pi*z*vC(7))*vC(4)*    &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                  Sin(pi*z*vC(7))*vC(4))*vC(7)**2)/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*Sin(pi*x*wC(5))*wC(2)*    &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                  Sin(pi*z*wC(7))*wC(4))*wC(5)**2)/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                  Sin(pi*z*uC(7))*uC(4))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) -     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                  Sin(pi*z*uC(7))*uC(4))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*pC(7))*pC(4)*pC(7)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                  Sin(pi*z*uC(7))*uC(4))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                  Sin(pi*z*uC(7))*uC(4))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*uC(7))*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*uC(4)*uC(7)*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                  Sin(pi*z*wC(7))*wC(4))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) -     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                  Sin(pi*z*wC(7))*wC(4))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*pC(5))*pC(2)*pC(5)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                  Sin(pi*z*wC(7))*wC(4))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                  Sin(pi*z*wC(7))*wC(4))*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*wC(5))*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*wC(2)*wC(5)*    &
                (pi*Cos(pi*z*uC(7))*uC(4)*uC(7) + pi*Cos(pi*x*wC(5))*wC(2)*wC(5)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*Sin(pi*y*wC(6))*wC(3)*    &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                  Sin(pi*z*wC(7))*wC(4))*wC(6)**2)/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                  Sin(pi*z*vC(7))*vC(4))*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) -     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                  Sin(pi*z*vC(7))*vC(4))*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*pC(7))*pC(4)*pC(7)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                  Sin(pi*z*vC(7))*vC(4))*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                  Sin(pi*z*vC(7))*vC(4))*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*vC(7))*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*vC(4)*vC(7)*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) -     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                  Sin(pi*z*wC(7))*wC(4))*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                        pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                    (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7))))**2) -     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                  Sin(pi*z*wC(7))*wC(4))*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                   rC(4)*Sin(pi*z*rC(7)))**2*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*pC(6))*pC(3)*pC(6)*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                  Sin(pi*z*wC(7))*wC(4))*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*(1 + S_div_TRef_Sutherland)*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                  pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))*    &
                (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                       (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))**2) +     &
                  (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))*    &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                  Sin(pi*z*wC(7))*wC(4))*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              (2.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*    &
                (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*wC(6))*    &
                (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                  pC(4)*Sin(pi*z*pC(7)))*    &
                Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                      pC(4)*Sin(pi*z*pC(7))))/    &
                  (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                    rC(4)*Sin(pi*z*rC(7))))*wC(3)*wC(6)*    &
                (pi*Cos(pi*z*vC(7))*vC(4)*vC(7) + pi*Cos(pi*y*wC(6))*wC(3)*wC(6)))/    &
              ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))*    &
                (S_div_TRef_Sutherland + (gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                       pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7))))) +     &
             pi*Cos(pi*x*uC(5))*uC(2)*uC(5)*    &
              ((2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*uC(5))*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*uC(2)*uC(5))/    &
                 ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) -     &
                (2*gammaM2*(1 + S_div_TRef_Sutherland)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))))) +     &
             pi*Cos(pi*y*vC(6))*vC(3)*vC(6)*    &
              ((2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*vC(6))*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*vC(3)*vC(6))/    &
                 ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) -     &
                (2*gammaM2*(1 + S_div_TRef_Sutherland)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))))) +     &
             pi*Cos(pi*z*wC(7))*wC(4)*wC(7)*    &
              ((2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*wC(7))*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*wC(4)*wC(7))/    &
                 ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) -     &
                (2*gammaM2*(1 + S_div_TRef_Sutherland)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))))) +     &
             (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                Sin(pi*z*uC(7))*uC(4))*    &
              ((-2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*uC(5))*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                          (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                            pC(4)*Sin(pi*z*pC(7))))/    &
                        (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                           rC(4)*Sin(pi*z*rC(7)))**2) +     &
                     (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))*uC(2)*uC(5))/    &
                 ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                         (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                           pC(4)*Sin(pi*z*pC(7))))/    &
                       (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                         rC(4)*Sin(pi*z*rC(7))))**2) -     &
                (2*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*Cos(pi*x*rC(5))*Cos(pi*x*uC(5))*rC(2)*rC(5)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*uC(2)*uC(5))/    &
                 ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7)))**2*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) +     &
                (2*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*Cos(pi*x*pC(5))*Cos(pi*x*uC(5))*pC(2)*pC(5)*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*uC(2)*uC(5))/    &
                 ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) +     &
                (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*uC(5))*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                          (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                            pC(4)*Sin(pi*z*pC(7))))/    &
                        (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                           rC(4)*Sin(pi*z*rC(7)))**2) +     &
                     (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))*uC(2)*uC(5))/    &
                 (Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) -     &
                (4*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*Sin(pi*x*uC(5))*uC(2)*uC(5)**2)/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) +     &
                (2*gammaM2*(1 + S_div_TRef_Sutherland)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                          (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                            pC(4)*Sin(pi*z*pC(7))))/    &
                        (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                           rC(4)*Sin(pi*z*rC(7)))**2) +     &
                     (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                         (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                           pC(4)*Sin(pi*z*pC(7))))/    &
                       (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                         rC(4)*Sin(pi*z*rC(7))))**2) +     &
                (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7)))**2*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) -     &
                (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*x*pC(5))*pC(2)*pC(5)*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) -     &
                (gammaM2*(1 + S_div_TRef_Sutherland)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   (-((gammaM2*pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                          (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                            pC(4)*Sin(pi*z*pC(7))))/    &
                        (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                           rC(4)*Sin(pi*z*rC(7)))**2) +     &
                     (gammaM2*pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                         pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))))) +     &
             (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                Sin(pi*z*vC(7))*vC(4))*    &
              ((-2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*vC(6))*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                          (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                            pC(4)*Sin(pi*z*pC(7))))/    &
                        (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                           rC(4)*Sin(pi*z*rC(7)))**2) +     &
                     (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))*vC(3)*vC(6))/    &
                 ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                         (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                           pC(4)*Sin(pi*z*pC(7))))/    &
                       (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                         rC(4)*Sin(pi*z*rC(7))))**2) -     &
                (2*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*Cos(pi*y*rC(6))*Cos(pi*y*vC(6))*rC(3)*rC(6)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*vC(3)*vC(6))/    &
                 ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7)))**2*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) +     &
                (2*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*Cos(pi*y*pC(6))*Cos(pi*y*vC(6))*pC(3)*pC(6)*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*vC(3)*vC(6))/    &
                 ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) +     &
                (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*vC(6))*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                          (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                            pC(4)*Sin(pi*z*pC(7))))/    &
                        (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                           rC(4)*Sin(pi*z*rC(7)))**2) +     &
                     (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))*vC(3)*vC(6))/    &
                 (Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) -     &
                (4*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*Sin(pi*y*vC(6))*vC(3)*vC(6)**2)/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) +     &
                (2*gammaM2*(1 + S_div_TRef_Sutherland)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                          (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                            pC(4)*Sin(pi*z*pC(7))))/    &
                        (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                           rC(4)*Sin(pi*z*rC(7)))**2) +     &
                     (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                         (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                           pC(4)*Sin(pi*z*pC(7))))/    &
                       (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                         rC(4)*Sin(pi*z*rC(7))))**2) +     &
                (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7)))**2*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) -     &
                (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*y*pC(6))*pC(3)*pC(6)*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) -     &
                (gammaM2*(1 + S_div_TRef_Sutherland)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   (-((gammaM2*pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                          (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                            pC(4)*Sin(pi*z*pC(7))))/    &
                        (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                           rC(4)*Sin(pi*z*rC(7)))**2) +     &
                     (gammaM2*pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                         pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))))) +     &
             (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                Sin(pi*z*wC(7))*wC(4))*    &
              ((-2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*wC(7))*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                          (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                            pC(4)*Sin(pi*z*pC(7))))/    &
                        (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                           rC(4)*Sin(pi*z*rC(7)))**2) +     &
                     (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))*wC(4)*wC(7))/    &
                 ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                         (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                           pC(4)*Sin(pi*z*pC(7))))/    &
                       (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                         rC(4)*Sin(pi*z*rC(7))))**2) -     &
                (2*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*Cos(pi*z*rC(7))*Cos(pi*z*wC(7))*rC(4)*rC(7)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*wC(4)*wC(7))/    &
                 ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7)))**2*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) +     &
                (2*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*Cos(pi*z*pC(7))*Cos(pi*z*wC(7))*pC(4)*pC(7)*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*wC(4)*wC(7))/    &
                 ((rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) +     &
                (gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*wC(7))*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                          (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                            pC(4)*Sin(pi*z*pC(7))))/    &
                        (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                           rC(4)*Sin(pi*z*rC(7)))**2) +     &
                     (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))*wC(4)*wC(7))/    &
                 (Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) -     &
                (4*gammaM2*pi**2*(1 + S_div_TRef_Sutherland)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*Sin(pi*z*wC(7))*wC(4)*wC(7)**2)/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) +     &
                (2*gammaM2*(1 + S_div_TRef_Sutherland)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                          (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                            pC(4)*Sin(pi*z*pC(7))))/    &
                        (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                           rC(4)*Sin(pi*z*rC(7)))**2) +     &
                     (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                         (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                           pC(4)*Sin(pi*z*pC(7))))/    &
                       (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                         rC(4)*Sin(pi*z*rC(7))))**2) +     &
                (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                      rC(4)*Sin(pi*z*rC(7)))**2*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) -     &
                (2*gammaM2*pi*(1 + S_div_TRef_Sutherland)*Cos(pi*z*pC(7))*pC(4)*pC(7)*    &
                   Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                         pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))) -     &
                (gammaM2*(1 + S_div_TRef_Sutherland)*    &
                   (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                     pC(4)*Sin(pi*z*pC(7)))*    &
                   (-((gammaM2*pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                          (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                            pC(4)*Sin(pi*z*pC(7))))/    &
                        (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                           rC(4)*Sin(pi*z*rC(7)))**2) +     &
                     (gammaM2*pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7))))*    &
                   (pi*Cos(pi*x*uC(5))*uC(2)*uC(5) + pi*Cos(pi*y*vC(6))*vC(3)*vC(6) +     &
                     pi*Cos(pi*z*wC(7))*wC(4)*wC(7)))/    &
                 (3.*Sqrt((gammaM2*(pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                         pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7))))/    &
                     (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                       rC(4)*Sin(pi*z*rC(7))))*    &
                   (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))*    &
                   (S_div_TRef_Sutherland + (gammaM2*    &
                        (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                          pC(4)*Sin(pi*z*pC(7))))/    &
                      (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                        rC(4)*Sin(pi*z*rC(7)))))))/Re
   
      end associate
      
   END SUBROUTINE ManufacturedSolutionSourceNS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ManufacturedSolutionSourceEuler( xx, t, Q  )
!
!     ------------------------
!     Source term for MS Euler
!     ------------------------
!
      IMPLICIT NONE
      
      REAL(KIND=RP) :: xx(3), t
      REAL(KIND=RP) :: Q(NCONS)
      
      REAL(KIND=RP) :: x, y, z
      
      associate ( gamma => thermodynamics % gamma )
      x = xx(1)
      y = xx(2)
      z = xx(3)
      
!
!     -------------------------
!     Mass equation Source term
!     -------------------------
!
      Q(1) = Q(1) + pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))    &
            + pi*Cos(pi*x*uC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*uC(2)*uC(5) +     &
          pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))    &
            + pi*Cos(pi*y*vC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*vC(3)*vC(6) +     &
          pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))    &
            + pi*Cos(pi*z*wC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*wC(4)*wC(7) 
      
!
!     -------------------------
!     u-velocity Source term
!     -------------------------
!
      Q(2) = Q(2) + pi*Cos(pi*x*pC(5))*pC(2)*pC(5) +     &
          pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
              Sin(pi*z*uC(7))*uC(4))**2 +     &
          2*pi*Cos(pi*x*uC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*uC(2)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
             Sin(pi*z*uC(7))*uC(4))*uC(5) +     &
          pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
             Sin(pi*z*uC(7))*uC(4))*(vC(1) + Sin(pi*x*vC(5))*vC(2) +     &
             Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4)) +     &
          pi*Cos(pi*y*uC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*uC(3)*uC(6)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))    &
            + pi*Cos(pi*y*vC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*(uC(1) + Sin(pi*x*uC(5))*uC(2) +     &
             Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*vC(3)*vC(6) +     &
          pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
             Sin(pi*z*uC(7))*uC(4))*(wC(1) + Sin(pi*x*wC(5))*wC(2) +     &
             Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4)) +     &
          pi*Cos(pi*z*uC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*uC(4)*uC(7)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))    &
            + pi*Cos(pi*z*wC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*(uC(1) + Sin(pi*x*uC(5))*uC(2) +     &
             Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*wC(4)*wC(7)
      
!
!     -------------------------
!     v-velocity Source term
!     -------------------------
!
      Q(3) = Q(3) + pi*Cos(pi*y*pC(6))*pC(3)*pC(6) +     &
          pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
             Sin(pi*z*uC(7))*uC(4))*(vC(1) + Sin(pi*x*vC(5))*vC(2) +     &
             Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4)) +     &
          pi*Cos(pi*x*uC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*uC(2)*uC(5)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))    &
            + pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
              Sin(pi*z*vC(7))*vC(4))**2 +     &
          pi*Cos(pi*x*vC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*(uC(1) + Sin(pi*x*uC(5))*uC(2) +     &
             Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*vC(2)*vC(5) +     &
          2*pi*Cos(pi*y*vC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*vC(3)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
             Sin(pi*z*vC(7))*vC(4))*vC(6) +     &
          pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
             Sin(pi*z*vC(7))*vC(4))*(wC(1) + Sin(pi*x*wC(5))*wC(2) +     &
             Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4)) +     &
          pi*Cos(pi*z*vC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*vC(4)*vC(7)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))    &
            + pi*Cos(pi*z*wC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*(vC(1) + Sin(pi*x*vC(5))*vC(2) +     &
             Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))*wC(4)*wC(7)
      
!
!     -------------------------
!     w-velocity Source term
!     -------------------------
!
      Q(4) = Q(4) + pi*Cos(pi*z*pC(7))*pC(4)*pC(7) +     &
          pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
             Sin(pi*z*uC(7))*uC(4))*(wC(1) + Sin(pi*x*wC(5))*wC(2) +     &
             Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4)) +     &
          pi*Cos(pi*x*uC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*uC(2)*uC(5)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))    &
            + pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
             Sin(pi*z*vC(7))*vC(4))*(wC(1) + Sin(pi*x*wC(5))*wC(2) +     &
             Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4)) +     &
          pi*Cos(pi*y*vC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*vC(3)*vC(6)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))    &
            + pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
              Sin(pi*z*wC(7))*wC(4))**2 +     &
          pi*Cos(pi*x*wC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*(uC(1) + Sin(pi*x*uC(5))*uC(2) +     &
             Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))*wC(2)*wC(5) +     &
          pi*Cos(pi*y*wC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*(vC(1) + Sin(pi*x*vC(5))*vC(2) +     &
             Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))*wC(3)*wC(6) +     &
          2*pi*Cos(pi*z*wC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*wC(4)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
             Sin(pi*z*wC(7))*wC(4))*wC(7) 
      
!
!     -------------------------
!     Energy Source term
!     -------------------------
!
      Q(5) = Q(5) + pi*Cos(pi*x*pC(5))*pC(2)*pC(5)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) + Sin(pi*z*uC(7))*uC(4))    &
            + pi*Cos(pi*x*uC(5))*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
             pC(4)*Sin(pi*z*pC(7)))*uC(2)*uC(5) +     &
          pi*Cos(pi*y*pC(6))*pC(3)*pC(6)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) + Sin(pi*z*vC(7))*vC(4))    &
            + pi*Cos(pi*y*vC(6))*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
             pC(4)*Sin(pi*z*pC(7)))*vC(3)*vC(6) +     &
          pi*Cos(pi*z*pC(7))*pC(4)*pC(7)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) + Sin(pi*z*wC(7))*wC(4))    &
            + pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
             Sin(pi*z*uC(7))*uC(4))*((pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             ((uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))**2 +     &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))**2 +     &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))**2)/2.) +     &
          pi*Cos(pi*x*uC(5))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*uC(2)*uC(5)*    &
           ((pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                pC(4)*Sin(pi*z*pC(7)))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             ((uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))**2 +     &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))**2 +     &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))**2)/2.) +     &
          pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
             Sin(pi*z*vC(7))*vC(4))*((pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             ((uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))**2 +     &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))**2 +     &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))**2)/2.) +     &
          pi*Cos(pi*y*vC(6))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*vC(3)*vC(6)*    &
           ((pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                pC(4)*Sin(pi*z*pC(7)))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             ((uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))**2 +     &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))**2 +     &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))**2)/2.) +     &
          pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
             Sin(pi*z*wC(7))*wC(4))*((pC(1) + pC(2)*Sin(pi*x*pC(5)) +     &
                pC(3)*Sin(pi*y*pC(6)) + pC(4)*Sin(pi*z*pC(7)))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             ((uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))**2 +     &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))**2 +     &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))**2)/2.) +     &
          (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) + rC(4)*Sin(pi*z*rC(7)))*    &
           (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
             Sin(pi*z*uC(7))*uC(4))*(-((pi*Cos(pi*x*rC(5))*rC(2)*rC(5)*    &
                  (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                    pC(4)*Sin(pi*z*pC(7))))/    &
                ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2)) +     &
             (pi*Cos(pi*x*pC(5))*pC(2)*pC(5))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             (2*pi*Cos(pi*x*uC(5))*uC(2)*    &
                 (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))*uC(5) +     &
                2*pi*Cos(pi*x*vC(5))*vC(2)*    &
                 (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))*vC(5) +     &
                2*pi*Cos(pi*x*wC(5))*wC(2)*    &
                 (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))*wC(5))/2.) +     &
          (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) + rC(4)*Sin(pi*z*rC(7)))*    &
           (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
             Sin(pi*z*vC(7))*vC(4))*(-((pi*Cos(pi*y*rC(6))*rC(3)*rC(6)*    &
                  (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                    pC(4)*Sin(pi*z*pC(7))))/    &
                ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2)) +     &
             (pi*Cos(pi*y*pC(6))*pC(3)*pC(6))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             (2*pi*Cos(pi*y*uC(6))*uC(3)*    &
                 (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))*uC(6) +     &
                2*pi*Cos(pi*y*vC(6))*vC(3)*    &
                 (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))*vC(6) +     &
                2*pi*Cos(pi*y*wC(6))*wC(3)*    &
                 (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))*wC(6))/2.) +     &
          pi*Cos(pi*z*wC(7))*(pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
             pC(4)*Sin(pi*z*pC(7)))*wC(4)*wC(7) +     &
          pi*Cos(pi*z*wC(7))*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
             rC(4)*Sin(pi*z*rC(7)))*wC(4)*    &
           ((pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                pC(4)*Sin(pi*z*pC(7)))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             ((uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))**2 +     &
                (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))**2 +     &
                (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))**2)/2.)*wC(7) +     &
          (rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) + rC(4)*Sin(pi*z*rC(7)))*    &
           (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
             Sin(pi*z*wC(7))*wC(4))*(-((pi*Cos(pi*z*rC(7))*rC(4)*rC(7)*    &
                  (pC(1) + pC(2)*Sin(pi*x*pC(5)) + pC(3)*Sin(pi*y*pC(6)) +     &
                    pC(4)*Sin(pi*z*pC(7))))/    &
                ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                     rC(4)*Sin(pi*z*rC(7)))**2)) +     &
             (pi*Cos(pi*z*pC(7))*pC(4)*pC(7))/    &
              ((-1 + gamma)*(rC(1) + rC(2)*Sin(pi*x*rC(5)) + rC(3)*Sin(pi*y*rC(6)) +     &
                  rC(4)*Sin(pi*z*rC(7)))) +     &
             (2*pi*Cos(pi*z*uC(7))*uC(4)*    &
                 (uC(1) + Sin(pi*x*uC(5))*uC(2) + Sin(pi*y*uC(6))*uC(3) +     &
                   Sin(pi*z*uC(7))*uC(4))*uC(7) +     &
                2*pi*Cos(pi*z*vC(7))*vC(4)*    &
                 (vC(1) + Sin(pi*x*vC(5))*vC(2) + Sin(pi*y*vC(6))*vC(3) +     &
                   Sin(pi*z*vC(7))*vC(4))*vC(7) +     &
                2*pi*Cos(pi*z*wC(7))*wC(4)*    &
                 (wC(1) + Sin(pi*x*wC(5))*wC(2) + Sin(pi*y*wC(6))*wC(3) +     &
                   Sin(pi*z*wC(7))*wC(4))*wC(7))/2.) 

         end associate
   
   END SUBROUTINE ManufacturedSolutionSourceEuler

END MODULE ManufacturedSolutionsNS