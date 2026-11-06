!
!////////////////////////////////////////////////////////////////////////
!
!      Manufactured solution for the compressible Navier-Stokes equations on a
!      non-conforming mesh with hanging nodes (2:1 h-mortars).
!
!      Density and velocity are those of the manufactured solution of
!      test/NavierStokes/LambVector, with x+y+z replaced by kx*x+ky*y+kz*z so that
!      they are periodic in the current mesh. The energy is set so that the
!      temperature is uniform (T = 1): Sutherland's law then gives the constant
!      viscosity 1/Re that the source term assumes.
!      The source term is generated with NavierStokes_MMS_hMortars.ipynb.
!
!////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
module MMS
   use SMConstants
   implicit none

   ! Size of the periodic domain (horses3d-meshes/cartesian-band-hanging-nodes/BAND_mesh.h5)
   ! No variation in z: the domain is only 1 unit thick (the p331 case uses order 1 in z)
   real(kind=RP), parameter :: Lx = 6.0_RP, Ly = 6.0_RP
   real(kind=RP), parameter :: kx = 2.0_RP*PI/Lx, ky = 2.0_RP*PI/Ly, kz = 0.0_RP

   logical :: useGaussLobatto = .true.
   !
   ! Gauss-Lobatto quadrature weights, orders 1..8
   !
   real(kind=RP) :: wGL_1(0:1)
   data wGL_1 / &
         1.00000000000000000e+00_RP, &
         1.00000000000000000e+00_RP /

   real(kind=RP) :: wGL_2(0:2)
   data wGL_2 / &
         3.33333333333333315e-01_RP, &
         1.33333333333333326e+00_RP, &
         3.33333333333333315e-01_RP /

   real(kind=RP) :: wGL_3(0:3)
   data wGL_3 / &
         1.66666666666666657e-01_RP, &
         8.33333333333333370e-01_RP, &
         8.33333333333333370e-01_RP, &
         1.66666666666666657e-01_RP /

   real(kind=RP) :: wGL_4(0:4)
   data wGL_4 / &
         1.00000000000000006e-01_RP, &
         5.44444444444444398e-01_RP, &
         7.11111111111111138e-01_RP, &
         5.44444444444444398e-01_RP, &
         1.00000000000000006e-01_RP /

   real(kind=RP) :: wGL_5(0:5)
   data wGL_5 / &
         6.66666666666666657e-02_RP, &
         3.78474956297846943e-01_RP, &
         5.54858377035486461e-01_RP, &
         5.54858377035486461e-01_RP, &
         3.78474956297846943e-01_RP, &
         6.66666666666666657e-02_RP /

   real(kind=RP) :: wGL_6(0:6)
   data wGL_6 / &
         4.76190476190476164e-02_RP, &
         2.76826047361565741e-01_RP, &
         4.31745381209862611e-01_RP, &
         4.87619047619047619e-01_RP, &
         4.31745381209862611e-01_RP, &
         2.76826047361565741e-01_RP, &
         4.76190476190476164e-02_RP /

   real(kind=RP) :: wGL_7(0:7)
   data wGL_7 / &
         3.57142857142857123e-02_RP, &
         2.10704227143506007e-01_RP, &
         3.41122692483504408e-01_RP, &
         4.12458794658703720e-01_RP, &
         4.12458794658703720e-01_RP, &
         3.41122692483504408e-01_RP, &
         2.10704227143506007e-01_RP, &
         3.57142857142857123e-02_RP /

   real(kind=RP) :: wGL_8(0:8)
   data wGL_8 / &
         2.77777777777777762e-02_RP, &
         1.65495361560805576e-01_RP, &
         2.74538712500161652e-01_RP, &
         3.46428510973046389e-01_RP, &
         3.71519274376417241e-01_RP, &
         3.46428510973046389e-01_RP, &
         2.74538712500161652e-01_RP, &
         1.65495361560805576e-01_RP, &
         2.77777777777777762e-02_RP /
   
   !
   ! Gauss quadrature weights, orders 1..5
   !
   real(kind=RP) :: wG_1(0:1)
   data wG_1 / &
         1.00000000000000000e+00_RP, &
         1.00000000000000000e+00_RP /

   real(kind=RP) :: wG_2(0:2)
   data wG_2 / &
         0.55555555555555569_RP, &
         0.88888888888888884_RP, &
         0.55555555555555569_RP /

   real(kind=RP) :: wG_3(0:3)
   data wG_3 / &
         0.34785484513745385_RP, &
         0.65214515486254632_RP, &
         0.65214515486254632_RP, &
         0.34785484513745385_RP /

   real(kind=RP) :: wG_4(0:4)
   data wG_4 / &
         0.23692688505618911_RP, &
         0.47862867049936669_RP, &
         0.56888888888888889_RP, &
         0.47862867049936669_RP, &
         0.23692688505618911_RP /

   real(kind=RP) :: wG_5(0:5)
   data wG_5 / &
         0.17132449237917019_RP, &
         0.36076157304813833_RP, &
         0.46791393457269093_RP, &
         0.46791393457269093_RP, &
         0.36076157304813833_RP, &
         0.17132449237917019_RP /

   real(kind=RP) :: wG_6(0:6)
   data wG_6 / &
         0.12948496616886979_RP, &
         0.27970539148927670_RP, &
         0.38183005050511903_RP, &
         0.41795918367346940_RP, &
         0.38183005050511903_RP, &
         0.27970539148927670_RP, &
         0.12948496616886979_RP /

   real(kind=RP) :: wG_7(0:7)
   data wG_7 / &
         0.10122853629037630_RP, &
         0.22238103445337470_RP, &
         0.31370664587788727_RP, &
         0.36268378337836193_RP, &
         0.36268378337836193_RP, &
         0.31370664587788727_RP, &
         0.22238103445337470_RP, &
         0.10122853629037630_RP /

   real(kind=RP) :: wG_8(0:8)
   data wG_8 / &
         8.1274388361574565E-002_RP, &
         0.18064816069485751_RP, &
         0.26061069640293549_RP, &
         0.31234707704000264_RP, &
         0.33023935500125978_RP, &
         0.31234707704000264_RP, &
         0.26061069640293549_RP, &
         0.18064816069485751_RP, &
         8.1274388361574565E-002_RP /

   contains
      pure subroutine evalManufacturedSolution(x, y, z, t, gamma, Mach, res)
         implicit none
         real(kind=RP), intent(in)  :: x, y, z, t, gamma, Mach
         real(kind=RP), intent(out) :: res(5)

         real(kind=RP) :: s, rho_ms, u_ms, v_ms, w_ms, e_ms

         s = kx*x + ky*y + kz*z
         rho_ms = 3.0_RP + sin(s)*exp(t)
         u_ms = cos(s)*exp(-t)
         v_ms = sin(s)*exp(-t)
         w_ms = cos(s)*exp(-t)
         ! p = rho/(gamma*Mach**2), i.e. T = 1
         e_ms = 1.0_RP / (gamma*(gamma - 1.0_RP)*Mach**2) + (u_ms**2+v_ms**2+w_ms**2) / 2.0_RP

         res(1) = rho_ms
         res(2) = rho_ms * u_ms
         res(3) = rho_ms * v_ms
         res(4) = rho_ms * w_ms
         res(5) = rho_ms * e_ms
      end subroutine evalManufacturedSolution

      subroutine chooseIntegrationWeights(Nx, Ny, Nz, i, j, k, wi, wj, wk)
         implicit none
         integer, intent(in) :: Nx, Ny, Nz
         integer, intent(in) :: i, j, k
         REAL(KIND=RP), intent(out) :: wi, wj, wk

         if (.not. useGaussLobatto) then
            select case(Nx)
            case(1) ; wi = wG_1(i)
            case(2) ; wi = wG_2(i)
            case(3) ; wi = wG_3(i)
            case(4) ; wi = wG_4(i)
            case(5) ; wi = wG_5(i)
            case(6) ; wi = wG_6(i)
            case(7) ; wi = wG_7(i)
            case(8) ; wi = wG_8(i)
            case default
               print *, "MMS: order", Nx, "not tabulated (max=8)" ; error stop
            end select
            select case(Ny)
            case(1) ; wj = wG_1(j)
            case(2) ; wj = wG_2(j)
            case(3) ; wj = wG_3(j)
            case(4) ; wj = wG_4(j)
            case(5) ; wj = wG_5(j)
            case(6) ; wj = wG_6(j)
            case(7) ; wj = wG_7(j)
            case(8) ; wj = wG_8(j)
            case default
               print *, "MMS: order", Ny, "not tabulated (max=8)" ; error stop
            end select
            select case(Nz)
            case(1) ; wk = wG_1(k)
            case(2) ; wk = wG_2(k)
            case(3) ; wk = wG_3(k)
            case(4) ; wk = wG_4(k)
            case(5) ; wk = wG_5(k)
            case(6) ; wk = wG_6(k)
            case(7) ; wk = wG_7(k)
            case(8) ; wk = wG_8(k)
            case default
               print *, "MMS: order", Nz, "not tabulated (max=8)" ; error stop
            end select
         else ! Gauss - Lobatto
            select case(Nx)
            case(1) ; wi = wGL_1(i)
            case(2) ; wi = wGL_2(i)
            case(3) ; wi = wGL_3(i)
            case(4) ; wi = wGL_4(i)
            case(5) ; wi = wGL_5(i)
            case(6) ; wi = wGL_6(i)
            case(7) ; wi = wGL_7(i)
            case(8) ; wi = wGL_8(i)
            case default
               print *, "MMS: order", Nx, "not tabulated (max=8)" ; error stop
            end select
            select case(Ny)
            case(1) ; wj = wGL_1(j)
            case(2) ; wj = wGL_2(j)
            case(3) ; wj = wGL_3(j)
            case(4) ; wj = wGL_4(j)
            case(5) ; wj = wGL_5(j)
            case(6) ; wj = wGL_6(j)
            case(7) ; wj = wGL_7(j)
            case(8) ; wj = wGL_8(j)
            case default
               print *, "MMS: order", Ny, "not tabulated (max=8)" ; error stop
            end select
            select case(Nz)
            case(1) ; wk = wGL_1(k)
            case(2) ; wk = wGL_2(k)
            case(3) ; wk = wGL_3(k)
            case(4) ; wk = wGL_4(k)
            case(5) ; wk = wGL_5(k)
            case(6) ; wk = wGL_6(k)
            case(7) ; wk = wGL_7(k)
            case(8) ; wk = wGL_8(k)
            case default
               print *, "MMS: order", Nz, "not tabulated (max=8)" ; error stop
            end select
         end if
      end subroutine chooseIntegrationWeights
end module MMS
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
         SUBROUTINE UserDefinedStartup
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            IMPLICIT NONE  
         END SUBROUTINE UserDefinedStartup
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
!           ------------------------------------------------
!
            use smconstants
            use physicsstorage
            use hexmeshclass
            use fluiddata
            use MMS
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
            integer       :: i, j, k, eID
            real(kind=RP) :: x(NDIM)

#if defined(NAVIERSTOKES)
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          Ny => mesh % elements(eID) % Nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  x = mesh % elements(eID) % geom % x(:,i,j,k)
                  call evalManufacturedSolution(x(1), x(2), x(3), 0.0_RP, thermodynamics_ % gamma, dimensionless_ % Mach, &
                                                mesh % elements(eID) % storage % Q(:,i,j,k))
               end do;        end do;        end do
               end associate
            end do
#endif
         end subroutine UserDefinedInitialCondition

         subroutine UserDefinedState1(x, t, nHat, Q & 
#ifdef FLOW
            , thermodynamics_ &
            , dimensionless_  &
            , refValues_ & 
#endif
#ifdef CAHNHILLIARD
            , multiphase_ &
#endif
            )
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
         end subroutine UserDefinedState1

#ifdef FLOW

         subroutine UserDefinedGradVars1(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            use VariableConversion, only: GetGradientValues_f
            implicit none
            real(kind=RP), intent(in)          :: x(NDIM)
            real(kind=RP), intent(in)          :: t
            real(kind=RP), intent(in)          :: nHat(NDIM)
            real(kind=RP), intent(in)          :: Q(NCONS)
            real(kind=RP), intent(inout)       :: U(NGRAD)
            procedure(GetGradientValues_f)     :: GetGradients
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
            real(kind=RP), intent(in)    :: U_x(NGRAD)
            real(kind=RP), intent(in)    :: U_y(NGRAD)
            real(kind=RP), intent(in)    :: U_z(NGRAD)
            real(kind=RP), intent(inout) :: flux(NCONS)
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedNeumann1
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors &
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
!           ----------------------------------------------------------
!           Called before every time-step to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use MonitorsClass
            use fluiddata
            use physicsstorage
            IMPLICIT NONE
            CLASS(HexMesh)               :: mesh
            REAL(KIND=RP)                :: time
            REAL(KIND=RP)                :: dt
            type(Monitor_t), intent(in) :: monitors
#ifdef FLOW
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)    :: multiphase_
#endif
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
#ifdef FLOW
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
, multiphase_ &
#endif
#ifdef ACOUSTIC
, Qbase, Lambbase, Lamb &
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
            use MMS, only: kx, ky, kz
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
#ifdef ACOUSTIC
            real(kind=RP),             intent(in)  :: Qbase(NCONSB)
            real(kind=RP),             intent(in)  :: Lambbase(NDIM)
            real(kind=RP),             intent(in)  :: Lamb(NDIM)
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            real(kind=RP) :: sn, cs, et, mu, gamma, M

#if defined(NAVIERSTOKES)
            sn = sin(kx*x(1) + ky*x(2) + kz*x(3))
            cs = cos(kx*x(1) + ky*x(2) + kz*x(3))
            et = exp(time)
            ! Dimensionless viscosity (1/Re). No heat flux: the temperature is uniform
            mu    = dimensionless_ % mu
            gamma = thermodynamics_ % gamma
            M     = dimensionless_ % Mach

            S(IRHO) = (cs**2*et*kx + cs**2*et*kz + 2.0_RP*cs*et*ky*sn + 3.0_RP*cs*ky + et**2*sn - et*kx*sn**2 - &
                et*kz*sn**2 - 3.0_RP*kx*sn - 3.0_RP*kz*sn)/et
            S(IRHOU) = mu*(4.0_RP*cs*kx**2 + cs*kx*kz + 3.0_RP*cs*ky**2 + 3.0_RP*cs*kz**2 + &
                kx*ky*sn)/(3.0_RP*et) + (cs**3*et*kx + cs**3*et*kz + 2.0_RP*cs**2*et*ky*sn + 3.0_RP*cs**2*ky - &
                2.0_RP*cs*et*kx*sn**2 - 2.0_RP*cs*et*kz*sn**2 - 3.0_RP*cs*et - 6.0_RP*cs*kx*sn - 6.0_RP*cs*kz*sn &
                - et*ky*sn**3 - 3.0_RP*ky*sn**2)/et**2 + cs*et*kx/(M**2*gamma)
            S(IRHOV) = mu*(cs*kx*ky + cs*ky*kz + 3.0_RP*kx**2*sn + 4.0_RP*ky**2*sn + &
                3.0_RP*kz**2*sn)/(3.0_RP*et) + (2.0_RP*cs**2*et*kx*sn + 2.0_RP*cs**2*et*kz*sn + 3.0_RP*cs**2*kx &
                + 3.0_RP*cs**2*kz + 3.0_RP*cs*et*ky*sn**2 + 6.0_RP*cs*ky*sn - et*kx*sn**3 - et*kz*sn**3 - &
                3.0_RP*et*sn - 3.0_RP*kx*sn**2 - 3.0_RP*kz*sn**2)/et**2 + cs*et*ky/(M**2*gamma)
            S(IRHOW) = mu*(3.0_RP*cs*kx**2 + cs*kx*kz + 3.0_RP*cs*ky**2 + 4.0_RP*cs*kz**2 + &
                ky*kz*sn)/(3.0_RP*et) + (cs**3*et*kx + cs**3*et*kz + 2.0_RP*cs**2*et*ky*sn + 3.0_RP*cs**2*ky - &
                2.0_RP*cs*et*kx*sn**2 - 2.0_RP*cs*et*kz*sn**2 - 3.0_RP*cs*et - 6.0_RP*cs*kx*sn - 6.0_RP*cs*kz*sn &
                - et*ky*sn**3 - 3.0_RP*ky*sn**2)/et**2 + cs*et*kz/(M**2*gamma)
            S(IRHOE) = 2.0_RP*mu*(2.0_RP*cs**2*kx**2 + cs**2*kx*kz + cs**2*ky**2 + 2.0_RP*cs**2*kz**2 + &
                2.0_RP*cs*kx*ky*sn + 2.0_RP*cs*ky*kz*sn - 2.0_RP*kx**2*sn**2 - kx*kz*sn**2 - ky**2*sn**2 - &
                2.0_RP*kz**2*sn**2)/(3.0_RP*et**2) + (2.0_RP*cs**4*et*kx + 2.0_RP*cs**4*et*kz + &
                4.0_RP*cs**3*et*ky*sn + 6.0_RP*cs**3*ky - 2.0_RP*cs**2*et**2*sn - 3.0_RP*cs**2*et*kx*sn**2 - &
                3.0_RP*cs**2*et*kz*sn**2 - 12.0_RP*cs**2*et - 12.0_RP*cs**2*kx*sn - 12.0_RP*cs**2*kz*sn - &
                3.0_RP*cs*ky*sn**2 - et**2*sn**3 - et*kx*sn**4 - et*kz*sn**4 - 6.0_RP*et*sn**2 - 3.0_RP*kx*sn**3 &
                - 3.0_RP*kz*sn**3)/(2.0_RP*et**3) + (cs**2*et*gamma*kx + cs**2*et*gamma*kz + &
                2.0_RP*cs*et*gamma*ky*sn + 3.0_RP*cs*gamma*ky + et**2*sn - et*gamma*kx*sn**2 - et*gamma*kz*sn**2 &
                - 3.0_RP*gamma*kx*sn - 3.0_RP*gamma*kz*sn)/(M**2*et*gamma*(gamma - 1.0_RP))
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
            use FTAssertions
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use MonitorsClass
            use MMS
#ifdef _HAS_MPI_
            use mpi
#endif
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
!
!           ---------------
!           Local variables
!           ---------------
!
#if defined(NAVIERSTOKES)
            CHARACTER(LEN=29)                  :: testName = "NS MMS h-mortars"
            INTEGER                            :: eID, i, j, k
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            REAL(KIND=RP)                      :: u(NCONS), x(NDIM), wi, wj, wk
            REAL(KIND=RP)                      :: error_mesh, expectedError
#ifdef _HAS_MPI_
            REAL(KIND=RP)                      :: localErrorMesh
            integer                            :: ierr
#endif

            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            ! L2 norm of the error in the conservative variables
            error_mesh = 0.0_RP
            DO eID = 1, mesh % no_of_elements
               associate( e => mesh % elements(eID) )
               DO k = 0, e % Nxyz(3) ; DO j = 0, e % Nxyz(2) ; DO i = 0, e % Nxyz(1)
                  x = e % geom % x(:,i,j,k)
                  call evalManufacturedSolution(x(1), x(2), x(3), time, thermodynamics_ % gamma, dimensionless_ % Mach, u)
                  call chooseIntegrationWeights(e % Nxyz(1), e % Nxyz(2), e % Nxyz(3), i, j, k, wi, wj, wk)
                  error_mesh = error_mesh + wi*wj*wk * e % geom % jacobian(i,j,k) &
                                          * sum((u - e % storage % Q(:,i,j,k))**2)
               END DO ; END DO ; END DO
               end associate
            END DO

#ifdef _HAS_MPI_
            if (mesh % no_of_elements .ne. mesh % no_of_allElements) then ! This is done because MPI_Process is giving compilation error
               localErrorMesh = error_mesh
               call mpi_allreduce(localErrorMesh, error_mesh, 1, MPI_DOUBLE, MPI_SUM, &
                                 MPI_COMM_WORLD, ierr)
            end if
#endif
            error_mesh = sqrt(error_mesh)
            print *, "Error in L2 norm: ", error_mesh

            ! The reference depends on the control file, identified by the polynomial order
            ! TODO: fill in the reference values and uncomment the assertion
            if ( all(mesh % elements(1) % Nxyz == [1,1,1]) ) then      ! hMortars_MMS_p1.control
               ! expectedError = 
            else if ( all(mesh % elements(1) % Nxyz == [3,3,3]) ) then ! hMortars_MMS_p3.control
               expectedError = 3.6179333676541878E-004_rp 
            else if ( all(mesh % elements(1) % Nxyz == [3,3,1]) ) then ! hMortars_MMS_p331.control
               ! expectedError = 
            end if

            ! Relative check: FTAssertEqual falls back to an absolute tolerance when expectedValue <= tol
            ! CALL FTAssert(test = abs(error_mesh - expectedError) <= 1.0e-8_RP*abs(expectedError), &
            !               msg  = "L2 norm (relative error)")

            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)

            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
            ELSE
               WRITE(6,*) testName, " ... Failed"
               error stop 99
            END IF
            WRITE(6,*)

            CALL finalizeSharedAssertionsManager
            CALL detachSharedAssertionsManager
#endif
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
