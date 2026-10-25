#include "Includes.h"
module AnalyticalSolution
   use SMConstants
   implicit none
   
   real(kind=rp), parameter :: inner_radius = 0.8_rp

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

   pure function evalSolution(x, y, t) result(res)
      implicit none
      real(kind=rp), intent(in) :: x, y, t
      real(kind=rp) :: res

      res = exp(-1.0_rp * (POW2(x-t) + POW2(y))) + 1.0_rp

   end function evalSolution

end module AnalyticalSolution

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
!              - by default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            use smconstants
            use physicsstorage
            use hexmeshclass
            use fluiddata
            use AnalyticalSolution
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
!           local variables
!           ---------------
!
            integer        :: eid, i, j, k
            real(kind=RP)  :: qq, u, v, w, p
#if defined(NAVIERSTOKES)
            real(kind=RP)  :: Q(NCONS), phi, theta
#endif

!
!           ---------------------------------------
!           Navier-Stokes default initial condition
!           ---------------------------------------
!
#if defined(NAVIERSTOKES)
            associate ( gammaM2 => dimensionless_ % gammaM2, &
                        gamma => thermodynamics_ % gamma )
            theta = refvalues_ % AOAtheta*(pi/180.0_RP)
            phi   = refvalues_ % AOAphi*(pi/180.0_RP)
      
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  u  = 1.0_rp
                  v  = 0.0_rp
                  w  = 0.0_rp
      
                  ! q(1) = exp(-1.0_rp * (POW2(mesh % elements(eID) % geom % x(1,i,j,k)) + POW2(mesh % elements(eID) % geom % x(2,i,j,k)))) + 1.0_rp
                  q(1) = evalSolution(mesh % elements(eID) % geom % x(1,i,j,k), mesh % elements(eID) % geom % x(2,i,j,k), -0.5_rp)
                  p    = 1.0_RP
                  q(2) = q(1)*u
                  q(3) = q(1)*v
                  q(4) = q(1)*w
                  q(5) = p/(gamma - 1._RP) + 0.5_RP*q(1)*(u**2 + v**2 + w**2)

                  mesh % elements(eID) % storage % q(:,i,j,k) = q 
               end do;        end do;        end do
               end associate
            end do

            end associate
#endif
         end subroutine UserDefinedInitialCondition
#ifdef FLOW
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
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

         subroutine UserDefinedGradVars1(x, t, nHat, Q, U, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)          :: x(NDIM)
            real(kind=RP), intent(in)          :: t
            real(kind=RP), intent(in)          :: nHat(NDIM)
            real(kind=RP), intent(in)          :: Q(NCONS)
            real(kind=RP), intent(inout)       :: U(NGRAD)
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
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors)
!
!           ----------------------------------------------------------
!           Called before every time-step to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)               :: mesh
            REAL(KIND=RP)                :: time
            REAL(KIND=RP)                :: dt
            type(Monitor_t), intent(in) :: monitors
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
#ifdef FLOW
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
, multiphase_ &
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
!
!           ---------------
!           Local variables
!           ---------------
!
            integer  :: i, j, k, eID
!
!           Usage example
!           -------------
!           S(:) = x(1) + x(2) + x(3) + time
   
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
            use AnalyticalSolution
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
            CHARACTER(LEN=29)                  :: testName           = "Density pulse h-mortars"
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, Nx, Ny, Nz
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
            REAL(KIND=RP) :: sol, solh, x(NDIM), wi, wj, wk
            REAL(KIND=RP) :: error_mesh, error_elem
#ifdef _HAS_MPI_
            REAL(KIND=RP) :: localErrorMesh
            integer :: ierr
#endif
!
!           ----------------------------------------------------------------
!           Expected solution: L2 error of the density with respect to the
!           advected pulse after 100 iterations (CFL 0.3, Roe, Mach 0.5,
!           Gauss-Lobatto, periodic boundaries)
!           ----------------------------------------------------------------
!
#if defined(NAVIERSTOKES)

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
                        sol = evalSolution(x(1), x(2), -0.5_rp + time)
                        solh = mesh % elements(eID) % storage % Q(1,i,j,k)

                        call chooseIntegrationWeights(Nx, Ny, Nz, i, j, k, wi, wj, wk)

                        error_elem = error_elem + wi*wj*wk &
                                 * mesh % elements(eID) % geom % jacobian(i,j,k) &
                                 * (sol - solh)**2
                     END DO
                  END DO
               END DO
               error_mesh = error_mesh + error_elem
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

            ! p=2 --> 2.7153096812333025E-005
            ! p=3 --> 3.8495148016535691E-007
            ! p=4 --> 2.7707680226187883E-009

            CALL FTAssertEqual(expectedValue = 2.7153096812333025E-005_rp, &
                               actualValue   = error_mesh, &
                               tol           = 1.0e-7_RP, &
                               msg           = "L2 norm")


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
      
