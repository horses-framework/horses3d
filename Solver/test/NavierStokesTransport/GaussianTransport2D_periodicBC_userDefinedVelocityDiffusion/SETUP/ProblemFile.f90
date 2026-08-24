#include "Includes.h"
module AJRConcentration
   USE SMConstants
   implicit none
   public

   real(kind=rp), parameter :: D0 = 0.05_rp
   
   contains
   pure function analyticalSolutionConcentration(xyz, t) result(c)
      USE SMConstants
      USE PhysicsStorage
      use FluidData
      IMPLICIT NONE  
!
!           ---------
!           Arguments
!           ---------
!
      REAL(KIND=RP), intent(in) :: xyz(3)
      REAL(KIND=RP), intent(in) :: t
!
!           ---------------
!           Local variables
!           ---------------
!
      REAL(KIND=RP)                         :: dim, x, y, z, c
      
      dim = 2.0_rp
      x = xyz(1)
      y = xyz(2)
      
      c = 2.0_rp + sin(PI * x) * cos(PI * y) * exp(- t)
      
   END function analyticalSolutionConcentration

   pure subroutine prescribedVelocity(xyz, t, vel)
      use SMConstants
      implicit none
      REAL(KIND=RP), intent(in) :: xyz(3)
      REAL(KIND=RP), intent(in) :: t
      REAL(KIND=RP), intent(out) :: vel(1:NDIM)
      vel(IX) = 1.0_rp + 0.5_rp * sin(PI * xyz(1)) * cos(PI * xyz(2))
      vel(IY) = 1.0_rp + 0.5_rp * cos(PI * xyz(1)) * sin(PI * xyz(2))
      vel(IZ) = 0.0_rp
   end subroutine prescribedVelocity

   pure subroutine prescribedVelocityDivergence(xyz, t, div)
      use SMConstants
      implicit none
      REAL(KIND=RP), intent(in) :: xyz(3)
      REAL(KIND=RP), intent(in) :: t
      REAL(KIND=RP), intent(out) :: div
      div = PI * cos(PI * xyz(1)) * cos(PI * xyz(2))
   end subroutine prescribedVelocityDivergence

   pure subroutine prescribedDiffusion(xyz, t, transportD)
      use SMConstants
      implicit none
      REAL(KIND=RP), intent(in) :: xyz(3)
      REAL(KIND=RP), intent(in) :: t
      REAL(KIND=RP), intent(out) :: transportD(1:NDIM, 1:NDIM)
      transportD = 0.0_rp
      transportD(IX,IX) = D0 * (1.0_rp + 0.5_rp * sin(PI * xyz(1)) * sin(PI * xyz(2)))
      transportD(IY,IY) = D0 * (1.0_rp + 0.5_rp * sin(PI * xyz(1)) * sin(PI * xyz(2)))
   end subroutine prescribedDiffusion
end module AJRConcentration

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
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
                                        )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            class(HexMesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
#if defined(TRANSPORT)
            integer :: eID, fe, fa, Nx, Ny, Nz

            ! Allocate the transport velocity, transport velocity divergence, and transportD
            do eID = 1, mesh % no_of_elements
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)
               ! Element
               allocate( mesh % elements(eID) % storage % transportVelocity(1:NDIM,0:Nx,0:Ny,0:Nz) )
               allocate( mesh % elements(eID) % storage % transportVelocityDivergence(0:Nx,0:Ny,0:Nz) )
               allocate( mesh % elements(eID) % storage % transportD(1:NDIM,1:NDIM,0:Nx,0:Ny,0:Nz) )
               ! Face
               do fe = 1, size(mesh % elements(eID) % faceIDs)
                  fa = mesh % elements(eID) % faceIDs(fe)
                  Nx = mesh % faces(fa) % Nf(1)
                  Ny = mesh % faces(fa) % Nf(2)
                   ! It might be already allocated from the other element
                  if (.not. allocated( mesh % faces(fa) % storage(1) % transportVelocity) ) then
                     allocate( mesh % faces(fa) % storage(1) % transportVelocity(1:NDIM,0:Nx,0:Ny) )
                     allocate( mesh % faces(fa) % storage(1) % transportD(1:NDIM,1:NDIM,0:Nx,0:Ny) )
                     allocate( mesh % faces(fa) % storage(2) % transportVelocity(1:NDIM,0:Nx,0:Ny) )
                     allocate( mesh % faces(fa) % storage(2) % transportD(1:NDIM,1:NDIM,0:Nx,0:Ny) )
                  end if
               end do
            end do
#endif
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedInitialCondition(mesh &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
                                        )
!
!           ------------------------------------------------
!           called to set the initial condition for the flow
!              - by default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            use smconstants
            use physicsstorage
            use hexmeshclass
            use fluiddata
            use AJRConcentration
            implicit none
            class(hexmesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            REAL(KIND=RP) :: x(3)
            INTEGER       :: i, j, k, eID
            integer       :: Nx, Ny, Nz

#if defined(TRANSPORT)
            ! x(1) = x0
            ! x(2) = y0
            ! x(3) = 0.0_rp
            ! maxc = analyticalSolutionConcentration(x, 0.0_rp) ! The maximum is at (x0+ut,y0+vt)
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx 
                        x = mesh % elements(eID) % geom % x(:,i,j,k)
                        mesh % elements(eID) % storage % Q(IRHO,i,j,k) = 1.0_rp
                        mesh % elements(eID) % storage % Q(IRHOU:IRHOE,i,j,k) = 0.0_rp
                        mesh % elements(eID) % storage % Q(IC,i,j,k) = analyticalSolutionConcentration(x, 0.0_rp)

                        call prescribedVelocity(x, 0.0_rp, mesh % elements(eID) % storage % transportVelocity(:,i,j,k))
                        call prescribedVelocityDivergence(x, 0.0_rp, mesh % elements(eID) % storage % transportVelocityDivergence(i,j,k))
                        call prescribedDiffusion(x, 0.0_rp, mesh % elements(eID) % storage % transportD(:,:,i,j,k))
                     END DO
                  END DO
               END DO 

               ! print *, "initial c: ", minval(mesh % elements(eID) % storage % Q(IC,:,:,:)), maxval(mesh % elements(eID) % storage % Q(IC,:,:,:))
               
            END DO
#endif            

         end subroutine UserDefinedInitialCondition
#if defined(NAVIERSTOKES)
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
!
!           -------------------------------------------------
!           Used to define an user defined boundary condition
!           -------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
         end subroutine UserDefinedState1

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

         subroutine UserDefinedNeumann1(x, t, nHat, U_x, U_y, U_z)
!
!           --------------------------------------------------------
!           Used to define a Neumann user defined boundary condition
!           --------------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: U_x(NGRAD)
            real(kind=RP), intent(inout)  :: U_y(NGRAD)
            real(kind=RP), intent(inout)  :: U_z(NGRAD)
         end subroutine UserDefinedNeumann1
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
#if defined(NAVIERSTOKES)
            use MonitorsClass
#endif
#if defined(TRANSPORT)
            use PhysicsStorage_NSTPT
            use AJRConcentration
#endif
            IMPLICIT NONE
            class(HexMesh)               :: mesh
            real(kind=RP)                :: time
            real(kind=RP)                :: dt
#if defined(NAVIERSTOKES)
            type(Monitor_t), intent(in) :: monitors
#else
            logical, intent(in) :: monitors
#endif
            REAL(KIND=RP) :: x(3)
            INTEGER       :: i, j, k, eID
            integer       :: Nx, Ny, Nz

#if defined(TRANSPORT)
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx 
                        x = mesh % elements(eID) % geom % x(:,i,j,k)
                        ! This is done to ensure that Sunderland's law is satisfied.
                        mesh % elements(eID) % storage % Q(IRHO,i,j,k) = 1.0_rp
                        mesh % elements(eID) % storage % Q(IRHOU:IRHOE,i,j,k) = 0.0_rp
                        ! Prescribe the transport data
                        call prescribedVelocity(x, time, mesh % elements(eID) % storage % transportVelocity(:,i,j,k))
                        call prescribedVelocityDivergence(x, time, mesh % elements(eID) % storage % transportVelocityDivergence(i,j,k))
                        call prescribedDiffusion(x, time, mesh % elements(eID) % storage % transportD(:,:,i,j,k))
                     END DO
                  END DO
               END DO 
               
            END DO
#endif 
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
#if defined(NAVIERSTOKES)
         subroutine UserDefinedSourceTermNS(xyz, Q, t, S, thermodynamics_, dimensionless_, refValues_)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
#if defined(TRANSPORT)
            use AJRConcentration
#endif
            IMPLICIT NONE
            real(kind=RP),             intent(in)  :: xyz(NDIM)
            real(kind=RP),             intent(in)  :: Q(NCONS)
            real(kind=RP),             intent(in)  :: t
            real(kind=RP),             intent(inout) :: S(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_

            real(kind=RP) :: x,y,z
            real(kind=RP) :: Pr, M, gamma, lambda, mu, Re, kappa
            x = xyz(1)
            y = xyz(2)
            z = xyz(3)
            Pr = dimensionless_ % Pr
            M = dimensionless_ % Mach
            gamma = thermodynamics_ % gamma
            lambda = thermodynamics_ % lambda
            mu = refValues_ % mu
            Re = dimensionless_ % Re
            kappa = mu * dimensionless_ % mu_to_kappa
!
!           Usage example
!           -------------
!           S(:) = x(1) + x(2) + x(3) + time
#if defined(TRANSPORT)
            S = 0.0_rp
            S(IC) = (pi**2*D0*sin(2.0_rp*pi*y)/4.0_rp + pi**2*D0*sin(pi*(x - y)) + pi**2*D0*sin(pi*(x + y)) + pi**2*D0*sin(pi*(2.0_rp*x - 2.0_rp*y))/4.0_rp - pi**2*D0*sin(pi*(2.0_rp*x + 2.0_rp*y))/4.0_rp - sin(pi*(x - y))/2.0_rp - sin(pi*(x + y))/2.0_rp + pi*sin(pi*(2.0_rp*x - 2.0_rp*y))/8.0_rp + pi*sin(pi*(2.0_rp*x + 2.0_rp*y))/8.0_rp + pi*cos(pi*(x + y)))*exp(-t)
#endif
         end subroutine UserDefinedSourceTermNS
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual &
#if defined(NAVIERSTOKES)
                                                    , thermodynamics_ &
                                                    , dimensionless_  &
                                                    , refValues_ & 
#endif   
#if defined(CAHNHILLIARD)
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
            use AJRConcentration
#ifdef _HAS_MPI_
            use mpi
#endif
            IMPLICIT NONE
            class(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#if defined(CAHNHILLIARD)
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
            CHARACTER(LEN=29)                  :: testName           = "Gaussian diffusion 2D"
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, Nx, Ny, Nz
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            integer                            :: rank
            integer                            :: fid
            REAL(KIND=RP) :: sol, solh, x(NDIM), wi, wj, wk
            REAL(KIND=RP) :: error_mesh, error_elem
#ifdef _HAS_MPI_
            REAL(KIND=RP) :: localErrorMesh
            integer :: ierr
#endif
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
            


            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            ! Compute point-wise L2 norm
            error_mesh = 0.0_rp
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)


               error_elem = 0.0_rp
               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx
                        x = mesh % elements(eID) % geom % x(:,i,j,k)
                        sol = analyticalSolutionConcentration(x, time)
                        solh = mesh % elements(eID) % storage % Q(IC,i,j,k)
                        
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
                        error_elem = error_elem + wi*wj*wk &
                                 * mesh % elements(eID) % geom % jacobian(i,j,k) &
                                 * abs(sol - solh)**2

                        ! error_elem = error_elem + norm2(u - uh)
                     END DO
                  END DO
               END DO
               error_mesh = error_mesh + error_elem
            END DO

            ! p=1, error: 9.5547610295434798E-004
            ! p=2, error: 6.4970907833013173E-007
            ! p=3, error: 1.3997998316342829E-007
            ! p=4, error: 2.6778186532990858E-011

#ifdef _HAS_MPI_
            if (mesh % no_of_elements .ne. mesh % no_of_allElements) then ! This is done because MPI_Process is giving compilation error
               localErrorMesh = error_mesh
               call mpi_allreduce(localErrorMesh, error_mesh, 1, MPI_DOUBLE, MPI_SUM, &
                                 MPI_COMM_WORLD, ierr)
            end if
#endif
            error_mesh = sqrt(error_mesh / mesh % no_of_allElements)
            print *, "Error in L2 norm: ", error_mesh

            CALL FTAssertEqual(expectedValue = 1.0_rp + 6.4970907833013173E-007, &
                               actualValue   = 1.0_rp + error_mesh, &
                               tol           = 1.0e-7_RP, &
                               msg           = "L2 norm")

            
            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
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
