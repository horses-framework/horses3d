module MMS
   use SMConstants
   use SolutionFile
   implicit none
   public evalManufacturedSolution, compareCubeFile

   contains
      pure subroutine evalManufacturedSolution(x, y, z, t, res)
         use SMConstants
         implicit none
         real(rp), intent(in)  :: x, y, z, t
         real(rp), intent(out) :: res(5)

         real(rp) :: rho_ms, u_ms, v_ms, w_ms, e_ms

         rho_ms = 3.0_rp + sin(x+y+z)*exp(t)
         u_ms = cos(x+y+z)*exp(-t)
         v_ms = sin(x+y+z)*exp(-t)
         w_ms = cos(x+y+z)*exp(-t)
         e_ms = rho_ms * (u_ms**2+v_ms**2+w_ms**2) / 2.0_rp + 1.0_rp

         res(1) = rho_ms
         res(2) = rho_ms * u_ms
         res(3) = rho_ms * v_ms
         res(4) = rho_ms * w_ms
         res(5) = rho_ms * e_ms
      end subroutine

      pure subroutine evalLambVector(x, y, z, t, res)
         use SMConstants
         implicit none
         real(rp), intent(in)  :: x, y, z, t
         real(rp), intent(out) :: res(NDIM)

         res(IX) = (exp(-t)*sin(x+y+z)+exp(-t)*cos(x+y+z))*exp(-t)*sin(x+y+z)
         res(IY) = (-2.0_rp*exp(-t)*sin(x+y+z)-2.0_rp*exp(-t)*cos(x+y+z))*exp(-t)*cos(x+y+z)
         res(IZ) = (exp(-t)*sin(x+y+z)+exp(-t)*cos(x+y+z))*exp(-t)*sin(x+y+z)
      end subroutine

      pure function evalLambVectorDivergence(x, y, z, t) result(res)
         use SMConstants
         implicit none
         real(rp), intent(in)  :: x, y, z, t
         real(rp)              :: res
!
!        The Lamb vector depends only on x+y+z, and Lx+Ly+Lz = -2 exp(-2t) cos(2(x+y+z))
!
         res = 4.0_rp*exp(-2.0_rp*t)*sin(2.0_rp*(x+y+z))
      end function
!
!     Reads a file written with "lamb vector cube" record by record. Its records must be
!     the elements inside the cube, in order, with the analytical values (Lamb vector if
!     nVars = NDIM, its divergence if nVars = 1).
!
      subroutine compareCubeFile(fileName, nVars, elements, inside, t, tol, success)
         use SMConstants
         use ElementClass
         implicit none
         character(len=*), intent(in)    :: fileName
         integer,          intent(in)    :: nVars
         type(Element),    intent(in)    :: elements(:)
         integer,          intent(in)    :: inside(:)
         real(rp),         intent(in)    :: t, tol
         logical,          intent(inout) :: success

         integer               :: fid, n, i, j, k, rank, dims(4)
         real(rp)              :: x(NDIM), expected(NDIM), maxError
         real(rp), allocatable :: u_h(:,:,:,:)

         if ( getSolutionFileNoOfElements(fileName) /= size(inside) ) then
            print *, trim(fileName), ": number of elements ", getSolutionFileNoOfElements(fileName), &
                     ", expected ", size(inside)
            success = .false.
            return
         end if

         maxError = 0.0_rp
         open(newunit=fid, file=fileName, status="old", action="read", form="unformatted", access="stream")
         read(fid, pos=POS_INIT_DATA)
         do n = 1, size(inside)
            associate( e => elements(inside(n)) )
            read(fid) rank, dims
            if ( rank /= 4 .or. any(dims /= [nVars, e % Nxyz + 1]) ) then
               print *, trim(fileName), ": record ", n, " has dimensions ", dims
               success = .false.
               close(fid)
               return
            end if
            allocate(u_h(nVars, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
            read(fid) u_h
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               x = e % geom % x(:,i,j,k)
               if ( nVars == NDIM ) then
                  call evalLambVector(x(1), x(2), x(3), t, expected)
               else
                  expected(1) = evalLambVectorDivergence(x(1), x(2), x(3), t)
               end if
               maxError = max(maxError, maxval(abs(u_h(:,i,j,k) - expected(1:nVars))))
            end do                ; end do                ; end do
            deallocate(u_h)
            end associate
         end do
         close(fid)

         print *, trim(fileName), ": max. pointwise error ", maxError
         if ( maxError > tol ) success = .false.

      end subroutine compareCubeFile
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
            use MMS
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
            REAL(KIND=RP) :: rho , u , v , w , p
            REAL(KIND=RP) :: L, u_0, rho_0, p_0
            integer       :: Nx, Ny, Nz
#if defined(NAVIERSTOKES)
            
            L     = 1.0_RP
            u_0   = 1.0_RP
            rho_0 = 1.0_RP 
            p_0   = 100.0_RP

            associate( gamma => thermodynamics_ % gamma ) 
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx 
                        x = mesh % elements(eID) % geom % x(:,i,j,k)
                        call evalManufacturedSolution(x(1), x(2), x(3), 0.0_rp, mesh % elements(eID) % storage % Q(:,i,j,k))
                     END DO
                  END DO
               END DO 
               
            END DO 
            end associate
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
            IMPLICIT NONE
            class(HexMesh)               :: mesh
            real(kind=RP)                :: time
            real(kind=RP)                :: dt
#if defined(NAVIERSTOKES)
            type(Monitor_t), intent(in) :: monitors
#else
            logical, intent(in) :: monitors
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
            S    = 0.0_RP
            S(IRHO) = ((exp(t)*sin(x + y + z) - 2.0_rp*sin(x + y + z)**2 + sin(x + y + z)*cos(x + y + z) - 6.0_rp*exp(-t)*sin(x + y + z) + 3.0_rp*exp(-t)*cos(x + y + z))*exp(t) + exp(t)*sin(x + y + z)*cos(x + y + z) + 2.0_rp*exp(t)*cos(x + y + z)**2)*exp(-t)
            S(IRHOU) = (Re*(gamma*exp(t)*sin(x + y + z)**2*cos(x + y + z) + 2.0_rp*gamma*exp(t)*cos(x + y + z)**3.0_rp - exp(t)*sin(x + y + z)**2*cos(x + y + z) + exp(t)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)*cos(x + y + z)/2.0_rp) + Re*(gamma*exp(t)*sin(x + y + z)**3*cos(x + y + z) + 2.0_rp*gamma*exp(t)*sin(x + y + z)*cos(x + y + z)**3.0_rp + 6.0_rp*gamma*exp(-t)*sin(x + y + z)*cos(x + y + z) - exp(t)*sin(x + y + z)**3*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*cos(x + y + z)**3.0_rp + 2.0_rp*exp(t)*sin(x + y + z)*cos(x + y + z) - 6.0_rp*exp(-t)*sin(x + y + z)**2.0_rp - 30.0_rp*exp(-t)*sin(x + y + z)*cos(x + y + z) + 6.0_rp*exp(-t)*cos(x + y + z)**2.0_rp)*exp(t)/2.0_rp + Re*(gamma*((cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + 3)*exp(t)*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 12.0_rp*sin(x + y + z)*cos(x + y + z))*exp(-t)*sin(x + y + z)/4.0_rp + 3.0_rp*gamma*((cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + 3)*exp(t)*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 12.0_rp*sin(x + y + z)*cos(x + y + z))*exp(-2.0_rp*t)/4.0_rp + gamma*exp(t)*cos(x + y + z) - gamma*exp(-t)*sin(x + y + z)**2*cos(x + y + z) + gamma*exp(-t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - ((cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + 3)*exp(t)*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 12.0_rp*sin(x + y + z)*cos(x + y + z))*exp(-t)*sin(x + y + z)/4.0_rp - 3.0_rp*((cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + 3)*exp(t)*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 12.0_rp*sin(x + y + z)*cos(x + y + z))*exp(-2.0_rp*t)/4.0_rp - exp(t)*cos(x + y + z) - sin(x + y + z)*cos(x + y + z) - exp(-t)*sin(x + y + z)**3.0_rp + exp(-t)*sin(x + y + z)**2*cos(x + y + z) - 3.0_rp*exp(-t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + exp(-t)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)*cos(x + y + z)/2.0_rp - 3.0_rp*exp(-t)*cos(x + y + z))*exp(2.0_rp*t) - mu*(-exp(-t)*sin(x + y + z) - 11.0_rp*exp(-t)*cos(x + y + z))*exp(2.0_rp*t)/3.0_rp)*exp(-2.0_rp*t)/Re
            S(IRHOV) = (Re*(gamma*exp(t)*sin(x + y + z)**2*cos(x + y + z) + 2.0_rp*gamma*exp(t)*cos(x + y + z)**3.0_rp + exp(t)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)*cos(x + y + z) - 2.0_rp*exp(t)*cos(x + y + z)**3.0_rp) + Re*(gamma*exp(t)*sin(x + y + z)**3*cos(x + y + z) + 2.0_rp*gamma*exp(t)*sin(x + y + z)*cos(x + y + z)**3.0_rp + 6.0_rp*gamma*exp(-t)*sin(x + y + z)*cos(x + y + z) - exp(t)*sin(x + y + z)**3*cos(x + y + z) + 2.0_rp*exp(t)*sin(x + y + z)**2.0_rp - 2.0_rp*exp(t)*sin(x + y + z)*cos(x + y + z)**3.0_rp - 12.0_rp*exp(-t)*sin(x + y + z)**2.0_rp + 6.0_rp*exp(-t)*sin(x + y + z)*cos(x + y + z) + 12.0_rp*exp(-t)*cos(x + y + z)**2.0_rp)*exp(t)/2.0_rp + Re*(gamma*((cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + 3)*exp(t)*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 12.0_rp*sin(x + y + z)*cos(x + y + z))*exp(-t)*sin(x + y + z)/4.0_rp + 3.0_rp*gamma*((cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + 3)*exp(t)*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 12.0_rp*sin(x + y + z)*cos(x + y + z))*exp(-2.0_rp*t)/4.0_rp + gamma*exp(t)*cos(x + y + z) - gamma*exp(-t)*sin(x + y + z)**2*cos(x + y + z) + gamma*exp(-t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - ((cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + 3)*exp(t)*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 12.0_rp*sin(x + y + z)*cos(x + y + z))*exp(-t)*sin(x + y + z)/4.0_rp - 3.0_rp*((cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + 3)*exp(t)*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 12.0_rp*sin(x + y + z)*cos(x + y + z))*exp(-2.0_rp*t)/4.0_rp - exp(t)*cos(x + y + z) - sin(x + y + z)**2.0_rp - 2.0_rp*exp(-t)*sin(x + y + z)**3.0_rp + 3.0_rp*exp(-t)*sin(x + y + z)**2*cos(x + y + z) - exp(-t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 3.0_rp*exp(-t)*sin(x + y + z) + exp(-t)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)*cos(x + y + z))*exp(2.0_rp*t) - mu*(-10.0_rp*exp(-t)*sin(x + y + z) - 2.0_rp*exp(-t)*cos(x + y + z))*exp(2.0_rp*t)/3.0_rp)*exp(-2.0_rp*t)/Re
            S(IRHOW) = (Re*(gamma*exp(t)*sin(x + y + z)**2*cos(x + y + z) + 2.0_rp*gamma*exp(t)*cos(x + y + z)**3.0_rp - exp(t)*sin(x + y + z)**2*cos(x + y + z) + exp(t)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)*cos(x + y + z)/2.0_rp) + Re*(gamma*exp(t)*sin(x + y + z)**3*cos(x + y + z) + 2.0_rp*gamma*exp(t)*sin(x + y + z)*cos(x + y + z)**3.0_rp + 6.0_rp*gamma*exp(-t)*sin(x + y + z)*cos(x + y + z) - exp(t)*sin(x + y + z)**3*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*cos(x + y + z)**3.0_rp + 2.0_rp*exp(t)*sin(x + y + z)*cos(x + y + z) - 6.0_rp*exp(-t)*sin(x + y + z)**2.0_rp - 30.0_rp*exp(-t)*sin(x + y + z)*cos(x + y + z) + 6.0_rp*exp(-t)*cos(x + y + z)**2.0_rp)*exp(t)/2.0_rp + Re*(gamma*((cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + 3)*exp(t)*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 12.0_rp*sin(x + y + z)*cos(x + y + z))*exp(-t)*sin(x + y + z)/4.0_rp + 3.0_rp*gamma*((cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + 3)*exp(t)*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 12.0_rp*sin(x + y + z)*cos(x + y + z))*exp(-2.0_rp*t)/4.0_rp + gamma*exp(t)*cos(x + y + z) - gamma*exp(-t)*sin(x + y + z)**2*cos(x + y + z) + gamma*exp(-t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - ((cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + 3)*exp(t)*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 12.0_rp*sin(x + y + z)*cos(x + y + z))*exp(-t)*sin(x + y + z)/4.0_rp - 3.0_rp*((cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + 3)*exp(t)*cos(x + y + z) - 2.0_rp*exp(t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 12.0_rp*sin(x + y + z)*cos(x + y + z))*exp(-2.0_rp*t)/4.0_rp - exp(t)*cos(x + y + z) - sin(x + y + z)*cos(x + y + z) - exp(-t)*sin(x + y + z)**3.0_rp + exp(-t)*sin(x + y + z)**2*cos(x + y + z) - 3.0_rp*exp(-t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) + exp(-t)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)*cos(x + y + z)/2.0_rp - 3.0_rp*exp(-t)*cos(x + y + z))*exp(2.0_rp*t) - mu*(-exp(-t)*sin(x + y + z) - 11.0_rp*exp(-t)*cos(x + y + z))*exp(2.0_rp*t)/3.0_rp)*exp(-2.0_rp*t)/Re
            S(IRHOE) = -3.0_rp*M**2*gamma**2*kappa*(-5.0_rp*sin(x + y + z)*cos(2.0_rp*(x + y + z))/4.0_rp - 3.0_rp*sin(x + y + z)/4.0_rp - sin(2.0_rp*(x + y + z))*cos(x + y + z) + 3.0_rp*exp(-t)*sin(x + y + z)**2.0_rp - 3.0_rp*exp(-t)*cos(x + y + z)**2.0_rp)*exp(-t)/Re + 3.0_rp*M**2*gamma**2*kappa*exp(-2.0_rp*t)*sin(x + y + z)**2.0_rp/Re - 3.0_rp*M**2*gamma**2*kappa*exp(-2.0_rp*t)*cos(x + y + z)**2.0_rp/Re + 3.0_rp*M**2*gamma*kappa*(-5.0_rp*sin(x + y + z)*cos(2.0_rp*(x + y + z))/4.0_rp - 3.0_rp*sin(x + y + z)/4.0_rp - sin(2.0_rp*(x + y + z))*cos(x + y + z) + 3.0_rp*exp(-t)*sin(x + y + z)**2.0_rp - 3.0_rp*exp(-t)*cos(x + y + z)**2.0_rp)*exp(-t)/Re - 3.0_rp*M**2*gamma*kappa*exp(-2.0_rp*t)*sin(x + y + z)**2.0_rp/Re + 3.0_rp*M**2*gamma*kappa*exp(-2.0_rp*t)*cos(x + y + z)**2.0_rp/Re + gamma*(5.0_rp*exp(-t)*cos(x + y + z)/8.0_rp + 3.0_rp*exp(-t)*cos(3.0_rp*x + 3.0_rp*y + 3.0_rp*z)/8.0_rp - 3.0_rp*exp(-2.0_rp*t)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)/2.0_rp)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - gamma*(5.0_rp*exp(-t)*cos(x + y + z)/8.0_rp + 3.0_rp*exp(-t)*cos(3.0_rp*x + 3.0_rp*y + 3.0_rp*z)/8.0_rp - 3.0_rp*exp(-2.0_rp*t)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)/2.0_rp)*cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)/2.0_rp + gamma*(5.0_rp*exp(-t)*cos(x + y + z)/8.0_rp + 3.0_rp*exp(-t)*cos(3.0_rp*x + 3.0_rp*y + 3.0_rp*z)/8.0_rp - 3.0_rp*exp(-2.0_rp*t)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)/2.0_rp)/2.0_rp + 3.0_rp*gamma*(-exp(-t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)/2.0_rp + exp(-t)*cos(x + y + z)*cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)/4.0_rp + 3.0_rp*exp(-t)*cos(x + y + z)/4.0_rp - 3.0_rp*exp(-2.0_rp*t)*sin(x + y + z)*cos(x + y + z))*exp(-t)*sin(x + y + z) + 6.0_rp*gamma*(-exp(-t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)/2.0_rp + exp(-t)*cos(x + y + z)*cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)/4.0_rp + 3.0_rp*exp(-t)*cos(x + y + z)/4.0_rp - 3.0_rp*exp(-2.0_rp*t)*sin(x + y + z)*cos(x + y + z))*exp(-t)*cos(x + y + z) - 2.0_rp*gamma*sin(x + y + z)**2.0_rp + 2.0_rp*gamma*sin(x + y + z)*cos(x + y + z) + 2.0_rp*gamma*cos(x + y + z)**2.0_rp - gamma*exp(-t)*sin(x + y + z)**5.0_rp + gamma*exp(-t)*sin(x + y + z)**4*cos(x + y + z) - gamma*exp(-t)*sin(x + y + z)**3*cos(x + y + z)**2.0_rp + 2.0_rp*gamma*exp(-t)*sin(x + y + z)**2*cos(x + y + z)**3.0_rp + 2.0_rp*gamma*exp(-t)*sin(x + y + z)*cos(x + y + z)**4.0_rp - 6.0_rp*gamma*exp(-t)*sin(x + y + z) + 3.0_rp*gamma*exp(-t)*cos(x + y + z) - 5.0_rp*gamma*exp(-2.0_rp*t)*sin(x + y + z)**4.0_rp + 9.0_rp*gamma*exp(-2.0_rp*t)*sin(x + y + z)**3*cos(x + y + z)/2.0_rp - 6.0_rp*gamma*exp(-2.0_rp*t)*sin(x + y + z)**2*cos(x + y + z)**2.0_rp + 7.0_rp*gamma*exp(-2.0_rp*t)*sin(x + y + z)*cos(x + y + z)**3.0_rp + 4.0_rp*gamma*exp(-2.0_rp*t)*cos(x + y + z)**4.0_rp - 6.0_rp*gamma*exp(-3.0_rp*t)*sin(x + y + z)**3.0_rp + 3.0_rp*gamma*exp(-3.0_rp*t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 3.0_rp*gamma*exp(-3.0_rp*t)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)*cos(x + y + z) + 6.0_rp*gamma*exp(-3.0_rp*t)*cos(x + y + z)**3.0_rp + (-exp(-t)*sin(x + y + z)*cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)/4.0_rp - 3.0_rp*exp(-t)*sin(x + y + z)/4.0_rp - 3.0_rp*exp(-2.0_rp*t)*sin(x + y + z)**2.0_rp - 6.0_rp*exp(-2.0_rp*t)*cos(x + y + z)**2.0_rp)*exp(t)*sin(x + y + z) + exp(t)*sin(x + y + z) + sin(x + y + z)**4.0_rp/2.0_rp + sin(x + y + z)**2*cos(x + y + z)**2.0_rp + 3.0_rp*exp(-t)*sin(x + y + z)**3.0_rp/2.0_rp + 3.0_rp*exp(-t)*sin(x + y + z)*cos(x + y + z)**2.0_rp - 3.0_rp*exp(-t)*sin(x + y + z)*cos(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)/4.0_rp - 9.0_rp*exp(-t)*sin(x + y + z)/4.0_rp - exp(-2.0_rp*t)*sin(x + y + z)**4.0_rp - 3.0_rp*exp(-2.0_rp*t)*sin(x + y + z)**2*cos(x + y + z)**2.0_rp - 9.0_rp*exp(-2.0_rp*t)*sin(x + y + z)**2.0_rp + 2.0_rp*exp(-2.0_rp*t)*sin(x + y + z)*cos(x + y + z)**3.0_rp + 2.0_rp*exp(-2.0_rp*t)*cos(x + y + z)**4.0_rp - 18.0_rp*exp(-2.0_rp*t)*cos(x + y + z)**2.0_rp - 3.0_rp*exp(-3.0_rp*t)*sin(x + y + z)**3.0_rp + 9.0_rp*exp(-3.0_rp*t)*sin(x + y + z)**2*cos(x + y + z)/2.0_rp - 3.0_rp*exp(-3.0_rp*t)*sin(x + y + z)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z) - 18.0_rp*exp(-3.0_rp*t)*sin(x + y + z)*cos(x + y + z)**2.0_rp + 3.0_rp*exp(-3.0_rp*t)*sin(2.0_rp*x + 2.0_rp*y + 2.0_rp*z)*cos(x + y + z) + 3.0_rp*exp(-3.0_rp*t)*cos(x + y + z)**3.0_rp - 4.0_rp*mu*exp(-2.0_rp*t)*sin(x + y + z)**2.0_rp/Re + 8.0_rp*mu*exp(-2.0_rp*t)*sin(x + y + z)*cos(x + y + z)/(3.0_rp*Re) + 4.0_rp*mu*exp(-2.0_rp*t)*cos(x + y + z)**2.0_rp/Re
 

 
 
   
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
            use SMConstants
            use FTAssertions
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use MonitorsClass
            use MMS
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
#if defined(NAVIERSTOKES)
            CHARACTER(LEN=29)                  :: testName = "Lamb vector cube"
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            ! Same cube as "lamb vector cube" in LambVectorCube.control
            real(kind=RP), parameter           :: cube(6) = [0.0_RP, 2.0_RP, -1.0_RP, 1.5_RP, 0.5_RP, 3.0_RP]
            integer, allocatable               :: inside(:)
            logical                            :: success

            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            call mesh % findElementsInsideCube(cube, inside)
            CALL FTAssertEqual(expectedValue = .true., &
                               actualValue   = size(inside) > 0 .and. size(inside) < size(mesh % elements), &
                               msg           = "The cube should contain some, but not all, elements.")

            success = .true.
            call compareCubeFile("RESULTS/LambVectorCube_0000000020.Lamb.hsol", NDIM, mesh % elements, inside, &
                                 0.1_RP, 0.01_RP, success)
            CALL FTAssertEqual(expectedValue = .true., &
                               actualValue   = success, &
                               msg           = "Error when comparing the Lamb vector saved in the cube.")

            success = .true.
            call compareCubeFile("RESULTS/LambVectorCube_0000000020.divLamb.hsol", 1, mesh % elements, inside, &
                                 0.1_RP, 0.05_RP, success)
            CALL FTAssertEqual(expectedValue = .true., &
                               actualValue   = success, &
                               msg           = "Error when comparing the Lamb vector divergence saved in the cube.")

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
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
