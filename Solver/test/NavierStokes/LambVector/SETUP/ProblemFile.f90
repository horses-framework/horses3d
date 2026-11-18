module MMS
   use SMConstants
   use SolutionFile
   implicit none
   public evalManufacturedSolution, compareLambVector, compareLambVectorDivergence

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

      pure subroutine evalLambVectorIntegral(x, y, z, t, res)
         use SMConstants
         implicit none
         real(rp), intent(in)  :: x, y, z, t
         real(rp), intent(out) :: res(NDIM)

         res(IX) = (-sin(x + y + z)**2 - sin(x + y + z)*cos(x + y + z))*exp(-2.0_rp*t)/2.0_rp + sin(x + y + z)**2/2.0_rp + sin(x + y + z)*cos(x + y + z)/2.0_rp
         res(IY) = (sin(x + y + z)*cos(x + y + z) + cos(x + y + z)**2)*exp(-2.0_rp*t) - sin(x + y + z)*cos(x + y + z) - cos(x + y + z)**2
         res(IZ) = (-sin(x + y + z)**2 - sin(x + y + z)*cos(x + y + z))*exp(-2.0_rp*t)/2.0_rp + sin(x + y + z)**2/2.0_rp + sin(x + y + z)*cos(x + y + z)/2.0_rp
      end subroutine

      subroutine compareLambVector(Nelem, elements, success)
         use SMConstants
         use ElementClass
         implicit none
         integer, intent(in) :: Nelem
         type(Element)  :: elements(Nelem)
         logical, intent(inout) :: success

         ! Local variables
         INTEGER                            :: eID
         INTEGER                            :: i, j, k, Nx, Ny, Nz
         integer                            :: fid, pos
         REAL(KIND=RP) :: t, u(NDIM), x(NDIM)
         real(rp), allocatable :: u_h(:,:,:,:)
         real(rp) :: Qerror

         t = 0.1_rp

         open(newunit=fid, file="RESULTS/TGV_0000000020.Lamb.hsol", status="old", action="read", &
            form="unformatted" , access="stream")

         do eID = 1, SIZE( elements)
               Nx =  elements(eID) % Nxyz(1)
               Ny =  elements(eID) % Nxyz(2)
               Nz =  elements(eID) % Nxyz(3)

               ! Compare the values obtained (saved in file) with the analytical ones
               allocate(u_h(1:NDIM,0:Nx,0:Ny,0:Nz))
               pos = POS_INIT_DATA + ( elements(eID) % globID)*5_AddrInt*SIZEOF_INT + 1_AddrInt*NDIM* elements(eID) % offsetIO*SIZEOF_RP
               read(fid, pos=pos) u_h

               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx
                        x =  elements(eID) % geom % x(:,i,j,k)
                        call evalLambVector(x(1),x(2),x(3),t, u)

                        Qerror = maxval(abs(u_h(:,i,j,k) - u))
                        if (Qerror > 0.01_rp) then
                           print *, "eID: ", eID
                           write(STD_OUT,'(3(A,I0))') "i: ", i, ", j: ", j, ", k: ", k
                           print *, "Lamb vector analytical: ", u
                           print *, "Lamb vector read:       ", u_h(:,i,j,k)
                           success = .false.
                           return
                        end if
                     end do
                  end do
               end do
               deallocate(u_h)
         end do

      end subroutine compareLambVector

      subroutine compareLambVectorStats(Nelem, elements, t, success)
         use SMConstants
         use ElementClass
         implicit none
         integer, intent(in) :: Nelem
         type(Element)  :: elements(Nelem)
         real(rp), intent(in) :: t
         logical, intent(out) :: success

         ! Local variables
         INTEGER                            :: eID
         INTEGER                            :: i, j, k, Nx, Ny, Nz
         integer                            :: fid, pos
         REAL(KIND=RP) :: u(NDIM), x(NDIM)
         real(rp), allocatable :: u_h(:,:,:,:)
         real(rp) :: Qerror

         open(newunit=fid, file="RESULTS/TGV.Lamb.stats.hsol", status="old", action="read", &
            form="unformatted" , access="stream")

         do eID = 1, SIZE(elements)
               Nx =  elements(eID) % Nxyz(1)
               Ny =  elements(eID) % Nxyz(2)
               Nz =  elements(eID) % Nxyz(3)

               ! Compare the values obtained (saved in file) with the analytical ones
               allocate(u_h(1:NDIM,0:Nx,0:Ny,0:Nz))
               pos = POS_INIT_DATA + ( elements(eID) % globID)*5_AddrInt*SIZEOF_INT + 1_AddrInt*NDIM* elements(eID) % offsetIO*SIZEOF_RP
               read(fid, pos=pos) u_h

               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx
                        x =  elements(eID) % geom % x(:,i,j,k)
                        call evalLambVectorIntegral(x(1),x(2),x(3),t, u)
                        u = u / t

                        Qerror = maxval(abs(u_h(:,i,j,k) - u))
                        if (Qerror > 0.05_rp) then
                           print *, "eID: ", eID
                           write(STD_OUT,'(3(A,I0))') "i: ", i, ", j: ", j, ", k: ", k
                           print *, "Lamb vector integral analytical: ", u
                           print *, "Lamb vector integral read:       ", u_h(:,i,j,k)
                           print *, "error: ", Qerror
                           success = .false.
                           return
                        end if
                     end do
                  end do
               end do
               deallocate(u_h)
         end do

      end subroutine compareLambVectorStats

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

      subroutine compareLambVectorDivergence(Nelem, elements, t, maxError)
         use SMConstants
         use ElementClass
         implicit none
         integer,       intent(in)  :: Nelem
         type(Element), intent(in)  :: elements(Nelem)
         real(rp),      intent(in)  :: t
         real(rp),      intent(out) :: maxError

         integer               :: eID, i, j, k, Nx, Ny, Nz, fid
         integer(AddrInt)      :: pos
         real(rp)              :: x(NDIM)
         real(rp), allocatable :: divL_h(:,:,:,:)

         maxError = 0.0_rp
         open(newunit=fid, file="RESULTS/TGV_0000000020.divLamb.hsol", status="old", action="read", &
            form="unformatted" , access="stream")

         do eID = 1, Nelem
            Nx = elements(eID) % Nxyz(1)
            Ny = elements(eID) % Nxyz(2)
            Nz = elements(eID) % Nxyz(3)

            allocate(divL_h(1,0:Nx,0:Ny,0:Nz))
            pos = POS_INIT_DATA + elements(eID) % globID*5_AddrInt*SIZEOF_INT + 1_AddrInt*elements(eID) % offsetIO*SIZEOF_RP
            read(fid, pos=pos) divL_h

            do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
               x = elements(eID) % geom % x(:,i,j,k)
               maxError = max(maxError, abs(divL_h(1,i,j,k) - evalLambVectorDivergence(x(1),x(2),x(3),t)))
            end do       ; end do       ; end do
            deallocate(divL_h)
         end do

         close(fid)

      end subroutine compareLambVectorDivergence
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
            CHARACTER(LEN=29)                  :: testName           = "MMS Lamb Vector"
            REAL(KIND=RP)                      :: Qerror
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, Nx, Ny, Nz
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: Qsuccess
            integer                            :: rank
            integer                            :: fid
            REAL(KIND=RP) :: u(NCONS), uh(NCONS), x(NDIM), wi, wj, wk
            REAL(KIND=RP) :: error_mesh, error_elem
            logical :: useGaussLobatto = .false.
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
            
            Qsuccess = .true.

            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            open(newunit=fid, file="MMS_Q_degree5", status="old", action="read", &
                      form="unformatted" , access="stream")

            ! Compute point-wise L2 norm
            error_mesh = 0.0_rp
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               ! Compare the values obtained with the expected (saved in file)
               allocate(QExpected(1:NCONS,0:Nx,0:Ny,0:Nz))
               read(fid) QExpected(:,:,:,:)
               Qerror = maxval(abs(QExpected - mesh % elements(eID) % storage % Q))
               if (Qerror > 1e-12) then
                 Qsuccess = .false.
               end if
               deallocate(QExpected)


               error_elem = 0.0_rp
               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx
                        x = mesh % elements(eID) % geom % x(:,i,j,k)
                        call evalManufacturedSolution(x(1),x(2),x(3),time,u)
                        uh = mesh % elements(eID) % storage % Q(:,i,j,k)
                        
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
                                 * norm2(u - uh)**2

                        ! error_elem = error_elem + norm2(u - uh)
                     END DO
                  END DO
               END DO
               error_mesh = error_mesh + error_elem
            END DO

            ! close(fid)
            
            error_mesh = sqrt(error_mesh / SIZE(mesh % elements))
            print *, "Error in L2 norm: ", error_mesh

            CALL FTAssertEqual(expectedValue = 3.4755642143319206E-005_rp, &
                               actualValue   = error_mesh, &
                               tol           = 1.0e-7_RP, &
                               msg           = "L2 norm")
            
            CALL FTAssertEqual(expectedValue = .true., &
                               actualValue   = Qsuccess, &
                               msg           = "Error when comparing obtained Q from the reference one in file.")
            
            !
            ! Lamb vector comparison
            !
            Qsuccess = .true.
            call compareLambVector(SIZE(mesh % elements), mesh % elements, Qsuccess)
            CALL FTAssertEqual(expectedValue = .true., &
                               actualValue   = Qsuccess, &
                               msg           = "Error when comparing the Lamb vector.")
            
            Qsuccess = .true.
            call compareLambVectorStats(SIZE(mesh % elements), mesh % elements, time, Qsuccess)
            CALL FTAssertEqual(expectedValue = .true., &
                               actualValue   = Qsuccess, &
                               msg           = "Error when comparing the Lamb vector stats.")

            call compareLambVectorDivergence(SIZE(mesh % elements), mesh % elements, 0.1_rp, Qerror)
            print *, "Max. pointwise error of div(Lamb vector): ", Qerror
            CALL FTAssertEqual(expectedValue = .true., &
                               actualValue   = Qerror < 0.05_RP, &
                               msg           = "Error when comparing the divergence of the Lamb vector.")

            
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

            ! Save solution to file for the test
            ! block
            !    integer :: fid
            !    open(newunit=fid, file="MMS_Q_degree5", status="replace", action="write", &
            !          form="unformatted" , access="stream")
            !    DO eID = 1, SIZE(mesh % elements)
            !       write(fid) mesh % elements(eID) % storage % Q
            !    end do
            !    close(fid)
            ! end block
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

