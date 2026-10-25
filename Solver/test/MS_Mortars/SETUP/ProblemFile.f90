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
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
            use ManufacturedSolutionsNS
#endif
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
            INTEGER       :: i, j, k, eID
            REAL(KIND=RP) :: rho , u , v , w , p, Q(NCONS)
            REAL(KIND=RP), DIMENSION(7) :: rC, uC, vC, wC, pC
            
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
            rC = (/1._RP, 0.15_RP    ,-0.10000_RP ,-0.12000_RP, 1.0_RP      , 0.50_RP    , 1.50_RP/)
            uC = (/1._RP, 0.0625_RP  ,-0.03750_RP ,-0.02250_RP, 1.5_RP      , 0.60_RP    , 0.50_RP/)
            vC = (/1._RP,-0.09375_RP , 0.05000_RP , 0.03750_RP, 0.5_RP      , 2._RP/3._RP, 1.25_RP/)
            wC = (/1._RP, 0.01875_RP ,-0.03125_RP , 0.04375_RP, 1.0_RP/3._RP, 1._RP/5._RP, 1.00_RP/)
            pC = (/1._RP, 0.2_RP     , 0.50000_RP ,-0.35000_RP, 2.0_RP      , 1.00_RP    , 1._RP/3._RP/)
            
            !rC = (/1._RP, 0.15_RP    ,-0.1_RP    , 0._RP, 1._RP      , 0.50_RP    , 0._RP/)
            !uC = (/1._RP, 0.0625_RP  ,-0.0375_RP , 0._RP, 1.5_RP     , 0.60_RP    , 0._RP/)
            !vC = (/1._RP,-0.09375_RP , 0.05_RP   , 0._RP, 0.5_RP     , 2._RP/3._RP, 0._RP/)
           ! wC = (/0._RP, 0._RP      , 0._RP     , 0._RP, 0._RP      , 0.00_RP    , 0._RP/)
           ! pC = (/1._RP, 0.2_RP     , 0.5_RP    , 0._RP, 2.0_RP     , 1.00_RP    , 0._RP/)
            associate( gamma => thermodynamics_ % gamma ) 

            call random_seed()

            do eID = 1, SIZE(mesh % elements)
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          Ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3)  )

               do k = 0, Nz; do j = 0, Ny; do i = 0, Nx 
                  associate( x => mesh % elements(eID) % geom % x(1,i,j,k), &
                             y => mesh % elements(eID) % geom % x(2,i,j,k), &
                             z => mesh % elements(eID) % geom % x(3,i,j,k) )
                 ! write(*,*)'/////////////////////'
                 ! write(*,*) 'eID', eID, 'x', x  
                 ! write(*,*) 'eID', eID, 'y', y  
                 ! write(*,*) 'eID', eID, 'z', z  
                 ! write(*,*)'/////////////////////'
                  rho  = rC(1) + rC(2)*Sin(pi*rC(5)*x) + rC(3)*Sin(pi*rC(6)*y) + rC(4)*Sin(pi*rC(7)*z) 
                  u    = uC(1) + uC(2)*Sin(pi*uC(5)*x) + uC(3)*Sin(pi*uC(6)*y) + uC(4)*Sin(pi*uC(7)*z) 
                  v    = vC(1) + vC(2)*Sin(pi*vC(5)*x) + vC(3)*Sin(pi*vC(6)*y) + vC(4)*Sin(pi*vC(7)*z) 
                  w    = wC(1) + wC(2)*Sin(pi*wC(5)*x) + wC(3)*Sin(pi*wC(6)*y) + wC(4)*Sin(pi*wC(7)*z) 
                  p    = pC(1) + pC(2)*Sin(pi*pC(5)*x) + pC(3)*Sin(pi*pC(6)*y) + pC(4)*Sin(pi*pC(7)*z) 
                        
                  Q(1) = rho
                  Q(2) = rho*u
                  Q(3) = rho*v
                  Q(4) = rho*w
                  Q(5) = p/(gamma - 1.0_RP) + 0.5_RP*rho*(u**2 + v**2 + w**2)

                  mesh % elements(eID) % storage % Q(:,i,j,k) = Q

                  end associate
               end do;       end do;       end do 
               end associate
            end do 
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

            REAL(KIND=RP) :: rho , u , v , w , p
            REAL(KIND=RP), DIMENSION(7) :: rC, uC, vC, wC, pC
            
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
            rC = (/1._RP, 0.15_RP    ,-0.10000_RP ,-0.12000_RP, 1.0_RP      , 0.50_RP    , 1.50_RP/)
            uC = (/1._RP, 0.0625_RP  ,-0.03750_RP ,-0.02250_RP, 1.5_RP      , 0.60_RP    , 0.50_RP/)
            vC = (/1._RP,-0.09375_RP , 0.05000_RP , 0.03750_RP, 0.5_RP      , 2._RP/3._RP, 1.25_RP/)
            wC = (/1._RP, 0.01875_RP ,-0.03125_RP , 0.04375_RP, 1.0_RP/3._RP, 1._RP/5._RP, 1.00_RP/)
            pC = (/1._RP, 0.2_RP     , 0.50000_RP ,-0.35000_RP, 2.0_RP      , 1.00_RP    , 1._RP/3._RP/)
            
            !rC = (/1._RP, 0.15_RP    ,-0.1_RP    , 0._RP, 1._RP      , 0.50_RP    , 0._RP/)
            !uC = (/1._RP, 0.0625_RP  ,-0.0375_RP , 0._RP, 1.5_RP     , 0.60_RP    , 0._RP/)
            !vC = (/1._RP,-0.09375_RP , 0.05_RP   , 0._RP, 0.5_RP     , 2._RP/3._RP, 0._RP/)
            !wC = (/0._RP, 0._RP      , 0._RP     , 0._RP, 0._RP      , 0.00_RP    , 0._RP/)
            !pC = (/1._RP, 0.2_RP     , 0.5_RP    , 0._RP, 2.0_RP     , 1.00_RP    , 0._RP/)
            associate( gamma => thermodynamics_ % gamma ) 

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
#endif            
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
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_)
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
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
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
            use NodalStorageClass, only: NodalStorage
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
            CHARACTER(LEN=29)                  :: testName           = "Taylor-Green vortex"
            REAL(KIND=RP)                      :: maxError, locErr(NCONS), error(NCONS), Q(NCONS)
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, N
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
            integer                            :: rank
            real(kind=RP), parameter           :: kinEn = 0.12499758737106952_RP
            real(kind=RP), parameter           :: kinEnRate = -4.2807169311969659E-004_RP
            real(kind=RP), parameter           :: enstrophy = 0.37499411028501956_RP 
            real(kind=RP), parameter           :: res(5) = [5.229366754747479E-005_RP, &
                                                            0.12783424260634596_RP, &
                                                            0.12783424273963268_RP, &
                                                            0.24980299744783380_RP, &
                                                            0.61006093083852786_RP ]
            real(kind=RP), parameter           :: wP1(2) = [ 1.0_RP, 1.0_RP ]                            
            real(kind=RP), parameter           :: wP2(3) = [ 0.5555555555555555555556_RP, &
                                                             0.8888888888888888888889_RP, &
                                                             0.555555555555555555556_RP   ]
            real(kind=RP), parameter           :: wP3(4) = [ 0.3478548451374538573731_RP, &
                                                             0.6521451548625461426269_RP, &
                                                             0.6521451548625461426269_RP, &
                                                             0.3478548451374538573731_RP  ]
            real(kind=RP), parameter           :: wP4(5) = [ 0.2369268850561890875143_RP, &
                                                             0.4786286704993664680413_RP, &
                                                             0.5688888888888888888889_RP, &
                                                             0.4786286704993664680413_RP, &
                                                             0.2369268850561890875143_RP  ]
            real(kind=RP), parameter           :: wP5(6) = [ 0.1713244923791703450403_RP, &
                                                             0.3607615730481386075698_RP, &
                                                             0.4679139345726910473899_RP, &
                                                             0.46791393457269104739_RP,   &
                                                             0.3607615730481386075698_RP, &
                                                             0.1713244923791703450403_RP  ]
            real(kind=RP), parameter           :: wP6(7) = [ 0.1294849661688696932706_RP, &
                                                             0.2797053914892766679015_RP, &
                                                             0.38183005050511894495_RP,   &
                                                             0.417959183673469387755_RP,  &
                                                             0.38183005050511894495_RP,   &
                                                             0.279705391489276667901_RP,  &
                                                             0.129484966168869693271_RP   ]
            real(kind=RP), allocatable         :: w(:)
            integer :: l
            !SAME P ORDER FOR ALL THE DIRECTIONS
            N = mesh% elements(1)% Nxyz(1)
            allocate( w(0:N) )
            if( N .eq. 1 ) then 
              w(0:N) = wP1
            elseif( N .eq. 2 ) then 
              w(0:N) = wP2
            elseif( N .eq. 3 ) then
              w(0:N) = wP3
            elseif( N .eq. 4 ) then
              w(0:N) = wP4
            elseif( N .eq. 5 ) then
              w(0:N) = wP5
            elseif( N .eq. 6 ) then
              w(0:N) = wP6
            else
              write(*,*) ' Weights are wrong for the selected polynomial order, error is not computed' 
              return
            end if 
            
            error = 0.0_RP

            do eID = 1, mesh % no_of_elements
              associate( e => mesh% elements(eID) )
              do k = 0, e% Nxyz(3); do j = 0, e% Nxyz(2); do i = 0, e% Nxyz(1)

                call UserDefinedState1(e% geom% x(:,i,j,k), time, (/0.0_RP, 0.0_RP, 0.0_RP/), Q, thermodynamics_, dimensionless_, refValues_)                            
                                                
                locErr = e% storage% Q(:,i,j,k) - Q
                                                
                error = error + e% geom% jacobian(i,j,k)*w(i)*w(j)*w(k)*(locErr*locErr)
                                                
              end do               ; end do               ; end do
             ! if (mesh%nonconforming) then 
               !write(*,*) 'volum eid', eID, e% geom % Volume  
                !write(*,*) 'jac eid', eID, e% geom % jacobian
                !write(*,*) 'eid :', eID, 'error :', locErr
                !if (eID==9) write(*,*) 'ratjacelem',  mesh% elements(eID)%geom % jacobian / mesh% elements(6)%geom % jacobian
                !if (eID==9) write(*,*) 'ratjGradXielem',  mesh% elements(eID)%geom % jGradXi / mesh% elements(6)%geom % jGradXi
               ! !if (eID==9) write(*,*) 'ratjGradEtaelem',  mesh% elements(eID)%geom % jGradEta / mesh% elements(6)%geom % jGradEta
               ! if (eID==9) write(*,*) 'ratjGradZetaelem',  mesh% elements(eID)%geom % jGradZeta / mesh% elements(6)%geom % jGradZeta
               ! if (eID==9) then 
                 ! do l=1,6 
                    !if (mesh%faces(e%faceIDs(l))%IsMortar==1) then 
                     ! write(*,*) 'fmasterjac', mesh%faces(e%faceIDs(l))%geom%jacobian 
                     ! write(*,*) 'fratjac', mesh%faces(e%faceIDs(l))%geom%jacobian / mesh%faces(e%faceIDs(l)+2)%geom%jacobian
                      !write(*,*) 'fratjac', mesh%faces(e%faceIDs(l))%geom%jacobian / mesh%faces(e%faceIDs(l)+3)%geom%jacobian
                     ! write(*,*) 'fratjac', mesh%faces(e%faceIDs(l))%geom%jacobian / mesh%faces(e%faceIDs(l)+4)%geom%jacobian
                     ! write(*,*) '//////////////////'
                     ! write(*,*) 'fmastersurf', mesh%faces(e%faceIDs(l))%geom%surface 
                     ! write(*,*) 'fratsurf', mesh%faces(e%faceIDs(l))%geom%surface / mesh%faces(e%faceIDs(l)+2)%geom%surface
                     !! write(*,*) 'fratsurf', mesh%faces(e%faceIDs(l))%geom%surface / mesh%faces(e%faceIDs(l)+3)%geom%surface
                     ! write(*,*) 'fratsurf', mesh%faces(e%faceIDs(l))%geom%surface / mesh%faces(e%faceIDs(l)+4)%geom%surface
                     ! write(*,*) '//////////////////'
                     ! write(*,*) 'frmasternormal', mesh%faces(e%faceIDs(l))%geom%normal 
                     ! write(*,*) 'fratnormal', mesh%faces(e%faceIDs(l))%geom%normal / mesh%faces(e%faceIDs(l)+2)%geom%normal
                     ! write(*,*) 'fratnormal', mesh%faces(e%faceIDs(l))%geom%normal / mesh%faces(e%faceIDs(l)+3)%geom%normal
                     ! write(*,*) 'fratnormal', mesh%faces(e%faceIDs(l))%geom%normal / mesh%faces(e%faceIDs(l)+4)%geom%normal
                     ! write(*,*) '//////////////////'
                      !write(*,*) 'fmastertang1', mesh%faces(e%faceIDs(l))%geom%t1 
                     ! write(*,*) 'fmastertang2', mesh%faces(e%faceIDs(l))%geom%t2
                    !  write(*,*) 'frattang', mesh%faces(e%faceIDs(l))%geom%t1 / mesh%faces(e%faceIDs(l)+2)%geom%t1
                     ! write(*,*) 'frattang', mesh%faces(e%faceIDs(l))%geom%t1 / mesh%faces(e%faceIDs(l)+3)%geom%t1
                      !write(*,*) 'frattang', mesh%faces(e%faceIDs(l))%geom%t1 / mesh%faces(e%faceIDs(l)+4)%geom%t1
                     ! write(*,*) '//////////////////'
                      !write(*,*) 'fmasterGradXi', mesh%faces(e%faceIDs(l))%geom%GradXi 
                     ! write(*,*) 'fratGradXi', mesh%faces(e%faceIDs(l))%geom%GradXi / mesh%faces(e%faceIDs(l)+2)%geom%GradXi
                     ! write(*,*) 'fratGradXi', mesh%faces(e%faceIDs(l))%geom%GradXi / mesh%faces(e%faceIDs(l)+3)%geom%GradXi
                    !  write(*,*) 'fratGradXi', mesh%faces(e%faceIDs(l))%geom%GradXi / mesh%faces(e%faceIDs(l)+4)%geom%GradXi
                     ! write(*,*) '////////////////////////////////////////////////////////////////////////'
                     ! write(*,*) 'fmorratjac', mesh%faces(e%faceIDs(l)+1)%geom%jacobian 
                     ! write(*,*) 'fmorratsurf', mesh%faces(e%faceIDs(l)+1)%geom%surface 
                     ! write(*,*) 'fmorratnormal', mesh%faces(e%faceIDs(l)+1)%geom%normal
                     ! write(*,*) 'fmorrattang1', mesh%faces(e%faceIDs(l)+1)%geom%t1 
                     ! write(*,*) 'fmorrattang2', mesh%faces(e%faceIDs(l)+1)%geom%t2
                     ! write(*,*) 'fmorratGradXi', mesh%faces(e%faceIDs(l)+1)%geom%GradXi 
                   ! end if 
                 ! end do 
              !end if 
            !else 
              
              !do l=1,6 
               ! if (mesh%faces(e%faceIDs(l))%FaceType==1) then 
                 !   write(*,*) 'faceGradXi', mesh%faces(e%faceIDs(l))%geom%GradXi 
                 !   write(*,*) 'facet1', mesh%faces(e%faceIDs(l))%geom%t1 
                 !   write(*,*) 'facet2', mesh%faces(e%faceIDs(l))%geom%t2 
                 !   write(*,*) 'facenormal', mesh%faces(e%faceIDs(l))%geom%normal
               ! end if
              !end do  
            !end if 
              end associate
            end do

            write(*,*) 'error rho =', error(IRHO)
            write(*,*) 'error rhoU =', error(IRHOU)
            write(*,*) 'error rhoV =', error(IRHOV)
            write(*,*) 'error rhoW =', error(IRHOW)
            write(*,*) 'error rhoE =', error(IRHOE)            

            deallocate( w )

            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
            CALL FTAssertEqual(expectedValue = res(1) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(1,1) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "continuity residual")

            CALL FTAssertEqual(expectedValue = res(2) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(2,1) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "x-momentum residual")

            CALL FTAssertEqual(expectedValue = res(3) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(3,1) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "y-momentum residual")

            CALL FTAssertEqual(expectedValue = res(4) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(4,1) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "z-momentum residual")

            CALL FTAssertEqual(expectedValue = res(5) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(5,1) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "energy residual")

            CALL FTAssertEqual(expectedValue = kinEn, &
                               actualValue   = monitors % volumeMonitors(1) % values(1,1), &
                               tol           = 1.0e-11_RP, &
                               msg           = "Kinetic Energy")

            CALL FTAssertEqual(expectedValue = kinEnRate + 1.0_RP, &
                               actualValue   = monitors % volumeMonitors(2) % values(1,1) + 1.0_RP, &
                               tol           = 1.0e-11_RP, &
                               msg           = "Kinetic Energy Rate")

            CALL FTAssertEqual(expectedValue = enstrophy, &
                               actualValue   = monitors % volumeMonitors(3) % values(1,1), &
                               tol           = 1.0e-11_RP, &
                               msg           = "Enstrophy")

            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
               WRITE(6,*) "This test case has no expected solution yet, only checks the residual after 100 iterations."
            ELSE
               WRITE(6,*) testName, " ... Failed"
               WRITE(6,*) "NOTE: Failure is expected when the max eigenvalue procedure is changed."
               WRITE(6,*) "      If that is done, re-compute the expected values and modify this procedure"
                STOP 99
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
      