#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      use HyperbolicDiscretizations
      use EllipticDiscretizations
      use DGIntegrals
      use ShockCapturing
      use MeshTypes
      use HexMeshClass
      use ElementClass
      use PhysicsStorage
      use Physics
      use MPI_Face_Class
      use MPI_Process_Info
      use DGSEMClass
      use ParticlesClass
      use FluidData
      use VariableConversion, only: NSGradientVariables_STATE, GetNSViscosity, NSGradientVariables_ENTROPY, &
                                    GetGradientValues_f, NSGradientVariables_ENERGY, get_laminar_mu_kappa, &
                                    set_getVelocityGradients, GetNSKinematicViscosity
      use ProblemFileFunctions, only: UserDefinedSourceTermNS_f
      use BoundaryConditions
      use SpallartAlmarasTurbulence
      use ManufacturedSolutionsNSSA
      use IBMClass
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public   ComputeTimeDerivative, ComputeTimeDerivativeIsolated, viscousDiscretizationKey
      public   Initialize_SpaceAndTimeMethods, Finalize_SpaceAndTimeMethods


      abstract interface
      SUBROUTINE computeElementInterfaceFluxF(f, masterFace )
            use FaceClass
            IMPLICIT NONE
            TYPE(Face)   , INTENT(inout) :: f
            type(Face), optional, intent(inout) :: masterFace
         end subroutine computeElementInterfaceFluxF

         SUBROUTINE computeMPIFaceFluxF(f)
            use FaceClass
            IMPLICIT NONE
            TYPE(Face)   , INTENT(inout) :: f
         end subroutine computeMPIFaceFluxF

         SUBROUTINE computeBoundaryFluxF(f, time, mesh)
            use SMConstants
            use FaceClass,  only: Face
            use HexMeshClass
            IMPLICIT NONE
            type(Face),    intent(inout) :: f
            REAL(KIND=RP)                :: time
            type(HexMesh), intent(in)    :: mesh
         end subroutine computeBoundaryFluxF
      end interface

      procedure(computeElementInterfaceFluxF), pointer :: computeElementInterfaceFlux
      procedure(computeMPIFaceFluxF),          pointer :: computeMPIFaceFlux
      procedure(computeBoundaryFluxF),         pointer :: computeBoundaryFlux

      procedure(GetGradientValues_f),           pointer :: GetGradients
      procedure(EllipticFlux_f),                pointer :: ViscousFlux

      character(len=LINE_LENGTH), parameter  :: viscousDiscretizationKey = "viscous discretization"
!
!     ========
      CONTAINS
!     ========
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_SpaceAndTimeMethods(controlVariables, sem)
         use FTValueDictionaryClass
         use Utilities, only: toLower
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         use WallFunctionConnectivity
         implicit none
         class(FTValueDictionary),  intent(in)  :: controlVariables
         class(DGSem)                           :: sem
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)       :: inviscidDiscretizationName
         character(len=LINE_LENGTH)       :: viscousDiscretizationName
         character(len=*), parameter      :: gradient_variables_key = "gradient variables"
         character(len=LINE_LENGTH)       :: gradient_variables
         real(RP)                         :: hnmin, hnmax

         if (.not. sem % mesh % child) then ! If this is a child mesh, all these constructs were already initialized for the parent mesh

            call hnRange(sem % mesh, hnmin, hnmax)

            if ( MPI_Process % isRoot ) then
               write(STD_OUT,'(/)')
               call Section_Header("Spatial discretization scheme")
               write(STD_OUT,'(/)')

               write(STD_OUT,'(30X,A,A30,1pG10.3)') "->", "Minimum h/N: ", hnmin
               write(STD_OUT,'(30X,A,A30,1pG10.3)') "->", "Maximum h/N: ", hnmax
               write(STD_OUT,'(/)')
            end if
   !
   !        Initialize inviscid discretization
   !        ----------------------------------
            inviscidDiscretizationName = controlVariables % stringValueForKey(inviscidDiscretizationKey,requestedLength = LINE_LENGTH)

            call toLower(inviscidDiscretizationName)

            select case ( trim(inviscidDiscretizationName) )
            case ( "standard" )
               if (.not. allocated(HyperbolicDiscretization)) allocate( StandardDG_t  :: HyperbolicDiscretization )

            case ( "split-form")
               if (.not. allocated(HyperbolicDiscretization)) allocate( SplitDG_t     :: HyperbolicDiscretization)

            case default
               write(STD_OUT,'(A,A,A)') 'Requested inviscid discretization "',trim(inviscidDiscretizationName),'" is not implemented.'
               write(STD_OUT,'(A)') "Implemented discretizations are:"
               write(STD_OUT,'(A)') "  * Standard"
               write(STD_OUT,'(A)') "  * Split-Form"
               errorMessage(STD_OUT)
               error stop

            end select

            call HyperbolicDiscretization % Initialize(controlVariables)
   !
   !        Initialize viscous discretization
   !        ---------------------------------
            if ( flowIsNavierStokes ) then
               if ( controlVariables % ContainsKey(gradient_variables_key) ) then
                  gradient_variables = controlVariables % StringValueForKey(gradient_variables_key, LINE_LENGTH)
                  call toLower(gradient_variables)

                  select case (trim(gradient_variables))
                  case ("state")
                     call SetGradientVariables(GRADVARS_STATE)
                     GetGradients => NSGradientVariables_STATE
                     ViscousFlux  => ViscousFlux_STATE
                     call set_getVelocityGradients(GRADVARS_STATE)

                  case ("entropy")
                     call SetGradientVariables(GRADVARS_ENTROPY)
                     GetGradients => NSGradientVariables_ENTROPY
                     ViscousFlux  => ViscousFlux_ENTROPY
                     call set_getVelocityGradients(GRADVARS_ENTROPY)

                  case ("energy")
                     call SetGradientVariables(GRADVARS_ENERGY)
                     GetGradients => NSGradientVariables_ENERGY
                     ViscousFlux  => ViscousFlux_ENERGY
                     call set_getVelocityGradients(GRADVARS_ENERGY)

                  case default
                     print*, 'Entropy variables "',trim(gradient_variables),'" are not currently implemented'
                     write(STD_OUT,'(A)') "Implemented options are:"
                     write(STD_OUT,'(A)') "  * State"
                     write(STD_OUT,'(A)') "  * Entropy"
                     write(STD_OUT,'(A)') "  * Energy"
                     errorMessage(STD_OUT)
                     error stop
                  end select

               else
!
!                 Set state as default option
!                 ---------------------------
                  call SetGradientVariables(GRADVARS_STATE)
                  GetGradients => NSGradientVariables_STATE
                  ViscousFlux  => ViscousFlux_STATE
                  call set_getVelocityGradients(GRADVARS_STATE)

               end if

               if ( .not. controlVariables % ContainsKey(viscousDiscretizationKey) ) then
                  print*, "Input file is missing entry for keyword: viscous discretization"
                  errorMessage(STD_OUT)
                  error stop
               end if

               viscousDiscretizationName = controlVariables % stringValueForKey(viscousDiscretizationKey, requestedLength = LINE_LENGTH)
               call toLower(viscousDiscretizationName)

               select case ( trim(viscousDiscretizationName) )
               case("br1")
                  allocate(BassiRebay1_t     :: ViscousDiscretization)

               case("br2")
                  allocate(BassiRebay2_t     :: ViscousDiscretization)

               case("ip")
                  allocate(InteriorPenalty_t :: ViscousDiscretization)

               case default
                  write(STD_OUT,'(A,A,A)') 'Requested viscous discretization "',trim(viscousDiscretizationName),'" is not implemented.'
                  write(STD_OUT,'(A)') "Implemented discretizations are:"
                  write(STD_OUT,'(A)') "  * BR1"
                  write(STD_OUT,'(A)') "  * BR2"
                  write(STD_OUT,'(A)') "  * IP"
                  errorMessage(STD_OUT)
                  error stop

               end select

               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_NSSA)
               call ViscousDiscretization % Describe

            else
               if (.not. allocated(ViscousDiscretization)) allocate(EllipticDiscretization_t :: ViscousDiscretization)

               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_NSSA)
!              Set state as default option
!              ---------------------------
               call SetGradientVariables(GRADVARS_STATE)
               GetGradients => NSGradientVariables_STATE
               ViscousFlux  => ViscousFlux_STATE
               call set_getVelocityGradients(GRADVARS_STATE)

            end if

   !
   !        Initialize models
   !        -----------------
            call InitializeTurbulenceModel(SAmodel,controlVariables)

            call Initialize_ShockCapturing(ShockCapturingDriver, controlVariables, sem, &
                                           ComputeTimeDerivative, ComputeTimeDerivativeIsolated)

            call ShockCapturingDriver % Describe

         end if
!        --------------
         if (.not. sem % mesh % child) then
               computeElementInterfaceFlux => computeElementInterfaceFlux_NSSA
               computeMPIFaceFlux          => computeMPIFaceFlux_NSSA
               computeBoundaryFlux         => computeBoundaryFlux_NSSA
         end if

      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_SpaceAndTimeMethods
         implicit none
         IF ( ALLOCATED(HyperbolicDiscretization) ) DEALLOCATE( HyperbolicDiscretization )
         if ( allocated(ShockCapturingDriver) )     deallocate( ShockCapturingDriver )

      end subroutine Finalize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeTimeDerivative( mesh, particles, time, mode, HO_Elements, element_mask, Level)
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target           :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)                   :: time
         integer, intent(in)             :: mode
         logical, intent(in), optional   :: HO_Elements
         logical, intent(in), optional   :: element_mask(:)
		 integer, intent(in), optional   :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k

         call SetBoundaryConditionsEqn(NS_BC)
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time)
         call mesh % ProlongSolutionToFaces(NCONS)

!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution(NCONS)
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         if ( computeGradients ) then
            call ViscousDiscretization % ComputeGradient( NCONS, NGRAD, mesh , time, GetGradients)
         end if

#ifdef _HAS_MPI_
!$omp single
         if ( flowIsNavierStokes ) then
            call mesh % UpdateMPIFacesGradients(NGRAD)
         end if
!$omp end single
#endif
!         call ComputeArtificialViscosity(mesh)
!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call TimeDerivative_ComputeQDot(mesh = mesh , &
                                         particles = particles, &
                                         t    = time)
!$omp end parallel
!
      END SUBROUTINE ComputeTimeDerivative
!
!////////////////////////////////////////////////////////////////////////
!
!     This routine computes the time derivative element by element, without considering the Riemann Solvers
!     This is useful for estimating the isolated truncation error
!
      SUBROUTINE ComputeTimeDerivativeIsolated( mesh, particles, time, mode, HO_Elements, element_mask, Level)
         use EllipticDiscretizationClass
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target           :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)                   :: time
         integer,             intent(in) :: mode
         logical, intent(in), optional   :: HO_Elements
         logical, intent(in), optional   :: element_mask(:)
		 integer, intent(in), optional   :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time)
         call mesh % ProlongSolutionToFaces(NCONS)
!
!        -----------------------------------------------------
!        Compute LOCAL gradients and prolong them to the faces
!        -----------------------------------------------------
!
         if ( computeGradients ) then
            CALL BaseClass_ComputeGradient( ViscousDiscretization, NCONS, NGRAD, mesh , time , GetGradients)
!
!           The prolongation is usually done in the viscous methods, but not in the BaseClass
!           ---------------------------------------------------------------------------------
            call mesh % ProlongGradientsToFaces(NGRAD)
         end if

!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call TimeDerivative_ComputeQDotIsolated(mesh = mesh , &
                                                 t    = time )
!$omp end parallel
!
      END SUBROUTINE ComputeTimeDerivativeIsolated

      subroutine TimeDerivative_ComputeQDot( mesh , particles, t)
         use WallFunctionConnectivity
         implicit none
         type(HexMesh)              :: mesh
         type(Particles_t)          :: particles
         real(kind=RP)              :: t
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID, iFace, iEl, iP, m 
         real(kind=RP)  :: mu_smag, delta, mu_t, eta, kinematic_viscocity, mu_dim, &
                           Source(NCONS), TurbulentSource(NCONS)
         logical     :: isfirst = .TRUE.
!
!        ***********************************************
!        Compute the viscosity at the elements and faces
!        ***********************************************
!
         if (flowIsNavierStokes) then
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, size(mesh % elements)
               associate(e => mesh % elements(eID))
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call get_laminar_mu_kappa(e % storage % Q(:,i,j,k), e % storage % mu_NS(1,i,j,k), e % storage % mu_NS(2,i,j,k))
               end do                ; end do                ; end do
               end associate
            end do
!$omp end do
         end if

!$omp do schedule(runtime) private(i,j,k,mu_t, mu_dim, kinematic_viscocity)
            do eID = 1, size(mesh % elements)
               associate(e => mesh % elements(eID))
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)

                  mu_dim = e % storage % mu_NS(1,i,j,k) / dimensionless % mu

                  call GetNSKinematicViscosity(mu_dim, e % storage % Q(IRHO,i,j,k), kinematic_viscocity )

                  call SAmodel % ComputeViscosity(e % storage % Q(IRHOTHETA,i,j,k), kinematic_viscocity,&
                                                  e % storage % Q(IRHO,i,j,k), mu_dim,&
                                                  mu_t, e % storage % mu_NS(3,i,j,k) , e % geom % x(:,i,j,k))

                  e % storage % mu_NS(1,i,j,k) = e % storage % mu_NS(1,i,j,k) + mu_t * dimensionless % mu
                  e % storage % mu_NS(2,i,j,k) = e % storage % mu_NS(2,i,j,k) + mu_t * dimensionless % mu * dimensionless % mut_to_kappa_SA

                  call SAmodel % ComputeSourceTerms(e % storage % Q(IRHOTHETA,i,j,k), kinematic_viscocity, e % storage % Q(IRHO,i,j,k), &
                                                    e % geom % dWall(i,j,k), e % storage % Q(:,i,j,k),   &
                                                    e % storage % U_x(:,i,j,k), &
                                                    e % storage % U_y(:,i,j,k), &
                                                    e % storage % U_z(:,i,j,k), e % storage % S_SA(:,i,j,k), e % geom % x(:,i,j,k))
               end do                ; end do                ; end do
               end associate
            end do
!$omp end do

!        Compute viscosity at interior and boundary faces
!        ------------------------------------------------
         call compute_viscosity_at_faces(size(mesh % faces_interior), 2, mesh % faces_interior, mesh)
         call compute_viscosity_at_faces(size(mesh % faces_boundary), 1, mesh % faces_boundary, mesh)
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime)
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_VolumetricContribution( mesh, mesh % elements(eID) , t)
         end do
!$omp end do
!
#if defined(_HAS_MPI_)
!$omp single
         if (ShockCapturingDriver % isActive) then
            call mesh % UpdateMPIFacesAviscflux(NCONS)
         end if
!$omp end single
#endif
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) private(fID, m)
         do iFace = 1, size(mesh % faces_interior)
            fID = mesh % faces_interior(iFace)
            if (mesh % faces(fID) % MortarType == MORTAR_BIG) then 
               associate(fstar=>mesh% faces(fID)%storage(1)%fStar)
                  fstar=0.0_RP
               end associate
               associate(fstar=>mesh% faces(fID)%storage(2)%fStar)
                  fstar=0.0_RP
               end associate
               do m=1,4
                  if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                     call computeElementInterfaceFlux(masterFace=mesh % faces(fID), f=mesh % faces(mesh % faces(fID)%Mortar(m)))
                  end if 
               end do 
         elseif (mesh % faces(fID) % MortarType == MORTAR_NONE) then 
            call computeElementInterfaceFlux(mesh % faces(fID))
            end if 
         end do
!$omp end do nowait



!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % faces_boundary)
            fID = mesh % faces_boundary(iFace)
            call computeBoundaryFlux(mesh % faces(fID), t, mesh)
         end do
!$omp end do
!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
!
!$omp do schedule(runtime) private(i,j,k,eID)
         do iEl = 1, size(mesh % elements_sequential)
            eID = mesh % elements_sequential(iEl)
            associate(e => mesh % elements(eID))
            call TimeDerivative_FacesContribution(e, t, mesh)

            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k)
            end do         ; end do          ; end do
            end associate
         end do
!$omp end do
!
!        ****************************
!       Wait until messages are sent
!        ****************************
!
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
!$omp single
            if ( flowIsNavierStokes ) then
               call mesh % GatherMPIFacesGradients(NGRAD)
            else
               call mesh % GatherMPIFacesSolution(NCONS)
            end if
!$omp end single
!
!           Compute viscosity at MPI faces
!           ------------------------------
            call compute_viscosity_at_faces(size(mesh % faces_mpi), 2, mesh % faces_mpi, mesh)
!
!$omp single
            if ( flowIsNavierStokes ) then
               if ( ShockCapturingDriver % isActive ) then
                  call mpi_barrier(MPI_COMM_WORLD, ierr)     ! TODO: This can't be the best way :(
                  call mesh % GatherMPIFacesAviscflux(NCONS)
               end if
            end if
!$omp end single
!
!           **************************************
!           Compute Riemann solver of shared faces
!           **************************************
!
!$omp do schedule(runtime) private(fID, m)
            do iFace = 1, size(mesh % faces_mpi)
               fID = mesh % faces_mpi(iFace)
               if (mesh% faces(fID)%MortarType == MORTAR_BIG) then 
                  associate(fstar=>mesh% faces(fID)%storage(1)%fStar)
                     fstar=0.0_RP
                  end associate
                  do m=1,4
                     if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                        call computeElementInterfaceFlux(masterFace=mesh % faces(fID), f=mesh % faces(mesh % faces(fID)%Mortar(m)))
                     end if 
                  end do 
               end if 
               call computeMPIFaceFlux(mesh % faces(fID))
            end do
!$omp end do

!$omp single
            if ( mesh % nonconforming ) then
               call mesh % UpdateMPIFacesMortarflux(NCONS)
            end if
      !$omp end single
      
      
      !$omp single
            if ( mesh % nonconforming ) then
               call mesh % GatherMPIFacesMortarFlux(NCONS)         
            end if
      !$omp end single
!
!           ***********************************************************
!           Surface integrals and scaling of elements with shared faces
!           ***********************************************************
!
!$omp do schedule(runtime) private(i,j,k,eID)
            do iEl = 1, size(mesh % elements_mpi)
               eID = mesh % elements_mpi(iEl)
               associate(e => mesh % elements(eID))
               call TimeDerivative_FacesContribution(e, t, mesh)

               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k)
               end do         ; end do          ; end do
               end associate
            end do
!$omp end do
!
!           Add an MPI Barrier
!           ------------------
!$omp single
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
         end if
#endif

!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               e % storage % QDot(6,i,j,k) =   e % storage % QDot(6,i,j,k) + e % storage % S_SA(6,i,j,k)
            end do                  ; end do                ; end do
            end associate
         end do
!$omp end do
!
!        *****************************************************************************************************************************
!        Compute contributions to source term
!        ATTENTION: This is deactivated for child multigrid meshes since they have specially tailored source terms (already computed).
!                   If you are going to add contributions to the source term, do it adding to e % storage % S_NS inside the condition!
!        *****************************************************************************************************************************
         if (.not. mesh % child) then
!
!           Add physical source term
!           ************************
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do
!
!           Add Particles source
!           ********************
            if ( particles % active ) then
!$omp do schedule(runtime)
               do eID = 1, mesh % no_of_elements
                  associate ( e => mesh % elements(eID) )
                     e % storage % S_NSP = 0.0_RP
                  end associate
               enddo
!$omp end do

!$omp do schedule(runtime)
               do i = 1, particles % injection % injected + 1
                  if (particles % particle(i) % active) then

                     associate ( eID => particles % particle(i) % eID )

                     call particles % AddSource(i, mesh % elements( eID ), &
                        t, thermodynamics, dimensionless, refValues)

                     ! If this is uncommented, private(j) should be added to openmp.
                        !this commented section permits the computation of source term in neighbor elements
                     !do j = 1, 6
                     !   if (particles % particle(i) % mesh%elements( eID )%NumberOfConnections(j) > 0) then
                     !      call particles % AddSource(i, &
                     !      mesh % elements( particles % particle(i) % mesh%elements( eID )%Connection(j)%ElementIDs(1)  ), &
                     !      t, thermodynamics, dimensionless, refValues)
                     !   else
                     !      !
                     !   end if
                     !end do
                     end associate
                  endif
               end do
!$omp end do
!$omp do schedule(runtime)
               do eID = 1, mesh % no_of_elements
                  associate ( e => mesh % elements(eID) )
                     e % storage % S_NS = e % storage % S_NS + e % storage % S_NSP
                  end associate
               enddo
!$omp end do

            end if

       end if !(.not. mesh % child)

!        ***********************
!        Now add the source term
!        ***********************
!               print *, "and here"

!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               e % storage % QDot(:,i,j,k) =   e % storage % QDot(:,i,j,k) + e % storage % S_NS(:,i,j,k)
            end do                  ; end do                ; end do
            end associate
         end do
!$omp end do


!
!        *********************
!        Add IBM source term
!        *********************
         if( mesh% IBM% active ) then
            if( .not. mesh% IBM% semiImplicit ) then 
!$omp do schedule(runtime) private(i,j,k,Source)
               do eID = 1, mesh % no_of_elements
                  associate ( e => mesh % elements(eID) )
                  do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                     if( e% isInsideBody(i,j,k) ) then
                        call mesh% IBM% SourceTerm( eID = eID, Q = e % storage % Q(:,i,j,k), Source = Source, wallfunction = .false. )
                        e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + Source
                     end if
                  end do                  ; end do                ; end do
                  end associate
               end do
!$omp end do      
               if( mesh% IBM% Wallfunction ) then
!$omp single
                  call mesh% IBM% GetBandRegionStates( mesh% elements )
!$omp end single
!$omp do schedule(runtime) private(i,j,k,TurbulentSource)
                  do iP = 1, mesh% IBM% NumOfForcingPoints
                     associate( e => mesh% elements(mesh% IBM% ImagePoints(iP)% element_index)  )
                     i = mesh% IBM% ImagePoints(iP)% local_position(1)
                     j = mesh% IBM% ImagePoints(iP)% local_position(2)
                     k = mesh% IBM% ImagePoints(iP)% local_position(3)
                     call mesh % IBM % SourceTermTurbulence( mesh% IBM% ImagePoints(iP), e% storage% Q(:,i,j,k), &
                                                             e% geom% normal(:,i,j,k), e% geom% dWall(i,j,k),    &
                                                             e% STL(i,j,k), TurbulentSource                      )          
                     e% storage% QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + TurbulentSource  
                     end associate
                  end do
!$omp end do
               end if 
            end if 
         end if

      end subroutine TimeDerivative_ComputeQDot

      subroutine compute_viscosity_at_faces(no_of_faces, no_of_sides, face_ids, mesh)
         implicit none
         integer, intent(in)           :: no_of_faces
         integer, intent(in)           :: no_of_sides
         integer, intent(in)           :: face_ids(no_of_faces)
         class(HexMesh), intent(inout) :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: iFace, i, j, side
         real(kind=RP) :: delta, mu_smag, mu_t, kinematic_viscocity, mu_dim

         if (flowIsNavierStokes) then
!$omp do schedule(runtime) private(i,j)
            do iFace = 1, no_of_faces
               associate(f => mesh % faces(face_ids(iFace)))
                  if (f % MortarType == MORTAR_BIG) cycle 
               do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
                  do side = 1, no_of_sides
                      call get_laminar_mu_kappa(f % storage(side) % Q(:,i,j), f % storage(side) % mu_NS(1,i,j), f % storage(side) % mu_NS(2,i,j))
                  end do
               end do              ; end do
               end associate
            end do
!$omp end do
         end if

!$omp do schedule(runtime) private(i,j,mu_t, kinematic_viscocity, mu_dim)
            do iFace = 1, no_of_faces
               associate(f => mesh % faces(face_ids(iFace)))
                  if (f % MortarType == MORTAR_BIG) cycle 
               do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
                  do side = 1, no_of_sides

                  mu_dim = f % storage(side) % mu_NS(1,i,j) / dimensionless % mu

                  call GetNSKinematicViscosity(mu_dim, f % storage(side) % Q(IRHO,i,j), kinematic_viscocity )
                  call SAmodel % ComputeViscosity(f % storage(side) % Q(IRHOTHETA,i,j), kinematic_viscocity, &
                                                  f % storage(side) % Q(IRHO,i,j), mu_dim, &
                                                  mu_t, f % storage(side) % mu_NS(3,i,j))

                  f % storage(side) % mu_NS(1,i,j) = f % storage(side) % mu_NS(1,i,j) + mu_t * dimensionless % mu
                  f % storage(side) % mu_NS(2,i,j) = f % storage(side) % mu_NS(2,i,j) + mu_t * dimensionless % mu * dimensionless % mut_to_kappa_SA

               end do                ; end do                ; end do
               end associate
            end do
!$omp end do


      end subroutine compute_viscosity_at_faces
!
!////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------------
!     This routine computes Qdot neglecting the interaction with neighboring elements
!     and boundaries. Therefore, the external states are not needed.
!     -------------------------------------------------------------------------------
      subroutine TimeDerivative_ComputeQDotIsolated( mesh , t )
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, fID, iP 
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
         real(kind=rp) :: Source(NCONS), TurbulentSource(NCONS)
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime)
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_StrongVolumetricContribution(mesh, mesh % elements(eID) , t)
         end do
!$omp end do
!
!        *******************
!        Scaling of elements
!        *******************
!
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, size(mesh % elements)
            associate(e => mesh % elements(eID))

            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k)
            end do         ; end do          ; end do
            end associate
         end do
!$omp end do
!
!        Add a source term
!        -----------------
         if (.not. mesh % child) then
!$omp do schedule(runtime) private(i,j, k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do
         end if

!$omp do schedule(runtime) private(i,j)
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S_NS(:,i,j,k)
            end do                  ; end do                ; end do
            end associate
         end do
!$omp end do

!
!        *********************
!        Add IBM source term
!        *********************
         if( mesh% IBM% active ) then
            if( .not. mesh% IBM% semiImplicit ) then 
!$omp do schedule(runtime) private(i,j,k,Source)
              do eID = 1, mesh % no_of_elements
                 associate ( e => mesh % elements(eID) )
                  do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                     if( e% isInsideBody(i,j,k) ) then
                        call mesh% IBM% SourceTerm( eID = eID, Q = e % storage % Q(:,i,j,k), Source = Source, wallfunction = .false. )
                        e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + Source
                     end if
                  end do                  ; end do                ; end do
                  end associate
               end do
!$omp end do      
               if( mesh% IBM% Wallfunction ) then
!$omp single
                  call mesh% IBM% GetBandRegionStates( mesh% elements )
!$omp end single 
!$omp do schedule(runtime) private(i,j,k,TurbulentSource)
                  do iP = 1, mesh% IBM% NumOfForcingPoints
                     associate( e    => mesh% elements(mesh% IBM% ImagePoints(iP)% element_index) )
                     i = mesh% IBM% ImagePoints(iP)% local_position(1)
                     j = mesh% IBM% ImagePoints(iP)% local_position(2)
                     k = mesh% IBM% ImagePoints(iP)% local_position(3)
                     call mesh % IBM % SourceTermTurbulence( mesh% IBM% ImagePoints(iP), e% storage% Q(:,i,j,k), &
                                                             e% geom% normal(:,i,j,k), e% geom% dWall(i,j,k),    &
                                                             e% STL(i,j,k), TurbulentSource                      )             
                     e% storage% QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + TurbulentSource  
                     end associate
                  end do
!$omp end do
               end if 
            end if
         end if

      end subroutine TimeDerivative_ComputeQDotIsolated
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_StrongVolumetricContribution(mesh, e, t)
         use HexMeshClass
         use ElementClass
         implicit none
         type(HexMesh)             :: mesh
         type(Element)             :: e
         real(kind=RP), intent(in) :: t

!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: inviscidContravariantFlux ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP) :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: viscousContravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP) :: AviscContravariantFlux   ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP) :: contravariantFlux         ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , EulerFlux, inviscidContravariantFlux )
!
!        Compute viscous contravariant flux
!        ----------------------------------
         if (flowIsNavierStokes) then

            call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NGRAD, ViscousFlux, GetNSViscosity, e, viscousContravariantFlux)
!
!           Compute the artificial dissipation
!           ----------------------------------
            if (ShockCapturingDriver % isActive) then
               call ShockCapturingDriver % ComputeViscosity(mesh, e, AviscContravariantFlux)
            else
               AviscContravariantFlux = 0.0_RP
            end if

         else

            viscousContravariantFlux = 0.0_RP
            AviscContravariantFlux   = 0.0_RP

         end if
!
!        ************************
!        Perform volume integrals
!        ************************
!
         select type ( HyperbolicDiscretization )
         type is (StandardDG_t)
!
!           Compute the total Navier-Stokes flux
!           ------------------------------------
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux - AviscContravariantFlux
!
!           Perform the Weak Volume Green integral
!           --------------------------------------
            e % storage % QDot = ScalarStrongIntegrals % StdVolumeGreen ( e , NCONS, contravariantFlux )

         type is (SplitDG_t)
            error stop ':: TimeDerivative_StrongVolumetricContribution not implemented for split form'
!~ !
!~ !           Compute sharp fluxes for skew-symmetric approximations
!~ !           ------------------------------------------------------
!~             call HyperbolicDiscretization % ComputeSplitFormFluxes(e, inviscidContravariantFlux, fSharp, gSharp, hSharp)
!~ !
!~ !           Perform the Weak volume green integral
!~ !           --------------------------------------
!~             viscousContravariantFlux = viscousContravariantFlux

!~             e % storage % QDot = -ScalarWeakIntegrals % SplitVolumeDivergence( e, fSharp, gSharp, hSharp, viscousContravariantFlux)

         end select

      end subroutine TimeDerivative_StrongVolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_VolumetricContribution(mesh, e, t)
         use HexMeshClass
         use ElementClass
         implicit none
         type(HexMesh)             :: mesh
         type(Element)             :: e
         real(kind=RP), intent(in) :: t

!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: inviscidContravariantFlux ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP) :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: viscousContravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP) :: AviscContravariantFlux    ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP) :: contravariantFlux         ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , EulerFlux, inviscidContravariantFlux )

!
!        Compute viscous contravariant flux
!        ----------------------------------
         if (flowIsNavierStokes) then

            call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NGRAD, ViscousFlux, GetNSViscosity, e , viscousContravariantFlux)
!
!           Compute the artificial dissipation
!           ----------------------------------
            if (ShockCapturingDriver % isActive) then
               call ShockCapturingDriver % ComputeViscosity(mesh, e, AviscContravariantFlux)
            else
               AviscContravariantFlux = 0.0_RP
            end if

         else

            viscousContravariantFlux = 0.0_RP
            AviscContravariantFlux   = 0.0_RP

         end if
!
!        ************************
!        Perform volume integrals
!        ************************
!

         select type ( HyperbolicDiscretization )
         type is (StandardDG_t)
!
!           Compute the total Navier-Stokes flux
!           ------------------------------------
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux - AviscContravariantFlux
!
!           Perform the Weak Volume Green integral
!           --------------------------------------
            e % storage % QDot = ScalarWeakIntegrals % StdVolumeGreen ( e, NCONS, contravariantFlux )

         type is (SplitDG_t)
!
!           Compute sharp fluxes for skew-symmetric approximations
!           ------------------------------------------------------
            call HyperbolicDiscretization % ComputeSplitFormFluxes(e, inviscidContravariantFlux, fSharp, gSharp, hSharp)
!
!           Perform the Weak volume green integral
!           --------------------------------------
            viscousContravariantFlux = viscousContravariantFlux + AviscContravariantFlux

            e % storage % QDot = -ScalarWeakIntegrals % SplitVolumeDivergence( e, fSharp, gSharp, hSharp, viscousContravariantFlux)

         end select

      end subroutine TimeDerivative_VolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_FacesContribution( e , t , mesh)
         use HexMeshClass
         implicit none
         type(Element)           :: e
         real(kind=RP)           :: t
         type(HexMesh)           :: mesh
         e % storage % QDot = e % storage % QDot - ScalarWeakIntegrals % StdFace( e, NCONS, &
                      mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % fStar, &
                      mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % fStar, &
                      mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % fStar, &
                      mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % fStar, &
                      mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % fStar, &
                      mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % fStar )

      end subroutine TimeDerivative_FacesContribution
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!        Riemann solver drivers
!        ----------------------
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE computeElementInterfaceFlux_NSSA(f, masterFace)
         use FaceClass
         use RiemannSolvers_NSSA
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f
        type(Face), optional, intent(inout) :: masterFace 

         integer       :: i, j
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: Avisc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu_left(3), mu_right(3)
        integer        :: Sidearray(2)

!
  !        ---------------------------
  !        Artificial viscosity fluxes
  !        ---------------------------
  !
         if ( ShockCapturingDriver % isActive ) then
            Avisc_flux = 0.5_RP * (f % storage(1) % AviscFlux + f % storage(2) % AviscFlux)
         else
            Avisc_flux = 0.0_RP
         end if
  !
  !        --------------
  !        Viscous fluxes
  !        --------------

         visc_flux = 0._RP

         if (flowIsNavierStokes) then
            DO j = 0, f % Nf(2)
               DO i = 0, f % Nf(1)

                  mu_left(1) = f % storage(1) % mu_NS(1,i,j)
                  mu_right(1) = f % storage(2) % mu_NS(1,i,j)

                  mu_left(2)  = f % storage(1) % mu_NS(3,i,j)
                  mu_right(2) = f % storage(2) % mu_NS(3,i,j)

                  mu_left(3) = f % storage(1) % mu_NS(2,i,j)
                  mu_right(3) = f % storage(2) % mu_NS(2,i,j)
  !
  !                 --------------
  !                 Viscous fluxes
  !                 --------------
  !
                  CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                     EllipticFlux = ViscousFlux, &
                                                     f = f, &
                                                     QLeft = f % storage(1) % Q(:,i,j), &
                                                     QRight = f % storage(2) % Q(:,i,j), &
                                                     U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                     U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                     U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                     U_xRight = f % storage(2) % U_x(:,i,j), &
                                                     U_yRight = f % storage(2) % U_y(:,i,j), &
                                                     U_zRight = f % storage(2) % U_z(:,i,j), &
                                                     mu_left = mu_left, mu_right = mu_right, &
                                                     nHat = f % geom % normal(:,i,j) , &
                                                     dWall = f % geom % dWall(i,j), &
                                                     flux  = visc_flux(:,i,j) )

               end do
            end do
         end if

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
  !
  !              --------------
  !              Invscid fluxes
  !              --------------
  !
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )


  !
  !              Multiply by the Jacobian
  !              ------------------------

               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j) - Avisc_flux(:,i,j)

            END DO
         END DO
  !
  !        ---------------------------
  !        Return the flux to elements
  !        ---------------------------
  !
      if (f % MortarType == MORTAR_NONE) then  
         Sidearray = (/1,2/)
         call f % ProjectFluxToElements(NCONS, flux, Sidearray)
         end if 
     if (f % MortarType == MORTAR_SMALL4 .and. present(masterFace)) then 
      Sidearray = (/1,0/)
      call masterFace % ProjectMortarFluxToElements(nEqn=NCONS, whichElements=Sidearray, &
      slaveFace=f, MortarFlux=flux)
      Sidearray = (/0,2/)
      call f % ProjectFluxToElements(NCONS, flux, Sidearray)
     end if 
     ! end if 

      END SUBROUTINE computeElementInterfaceFlux_NSSA

      SUBROUTINE computeMPIFaceFlux_NSSA(f)
         use FaceClass
         use RiemannSolvers_NSSA
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: Avisc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu_left(3), mu_right(3)
         integer       :: Sidearray(2)
!
!        ---------------------------
!        Artificial viscosity fluxes
!        ---------------------------
!
         if ( ShockCapturingDriver % isActive ) then
            Avisc_flux = 0.5_RP * (f % storage(1) % AviscFlux + f % storage(2) % AviscFlux)
         else
            Avisc_flux = 0.0_RP
         end if
!
!        --------------
!        Viscous fluxes
!        --------------

         visc_flux = 0._RP
         if (flowIsNavierStokes) then
            DO j = 0, f % Nf(2)
               DO i = 0, f % Nf(1)

                  mu_left(1) = f % storage(1) % mu_NS(1,i,j)
                  mu_left(3) = f % storage(1) % mu_NS(2,i,j)

                  mu_left(2)  = f % storage(1) % mu_NS(3,i,j)
                  mu_right(2) = f % storage(2) % mu_NS(3,i,j)

                  mu_right(1) = f % storage(2) % mu_NS(1,i,j)
                  mu_right(3) = f % storage(2) % mu_NS(2,i,j)
!
!                 --------------
!                 Viscous fluxes
!                 --------------
!
                  CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                     EllipticFlux = ViscousFlux, &
                                                     f = f, &
                                                     QLeft = f % storage(1) % Q(:,i,j), &
                                                     QRight = f % storage(2) % Q(:,i,j), &
                                                     U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                     U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                     U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                     U_xRight = f % storage(2) % U_x(:,i,j), &
                                                     U_yRight = f % storage(2) % U_y(:,i,j), &
                                                     U_zRight = f % storage(2) % U_z(:,i,j), &
                                                     mu_left  = mu_left, &
                                                     mu_right = mu_right, &
                                                     nHat = f % geom % normal(:,i,j) , &
                                                     dWall = f % geom % dWall(i,j), &
                                                     flux  = visc_flux(:,i,j) )

               end do
            end do

         end if

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!
!              --------------
!              Invscid fluxes
!              --------------
!
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )


!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j) - Avisc_flux(:,i,j)

            END DO
         END DO
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)

         Sidearray = (/thisSide, HMESH_NONE/)
         call f % ProjectFluxToElements(NCONS, flux, Sidearray )
         if (f % MortarType == MORTAR_SMALL4) then 
            call f% Interpolatesmall2big(NCONS, flux)

         end if 
      end subroutine ComputeMPIFaceFlux_NSSA

      SUBROUTINE computeBoundaryFlux_NSSA(f, time, mesh)
      USE ElementClass
      use FaceClass
      use RiemannSolvers_NSSA
      use WallFunctionBC
      use WallFunctionConnectivity
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
      type(HexMesh), intent(in)    :: mesh
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: inv_flux(NCONS)
      real(kind=RP)                   :: visc_flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: Avisc_flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: fStar(NCONS, 0:f % Nf(1), 0: f % Nf(2))
      real(kind=RP)                   :: mu, kappa, beta, delta
      real(kind=RP)                   :: fv_3d(NCONS,NDIM)
      integer                         :: Sidearray(2)
      logical                         :: useWallFuncFace
      real(kind=RP)                   :: wallFunV(NDIM, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: wallFunVavg(NDIM, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: wallFunRho(0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: wallFunMu(0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: wallFunY(0:f % Nf(1), 0:f % Nf(2))


      if ( ShockCapturingDriver % isActive ) then
         do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
            Avisc_flux(:,i,j) = f % storage(1) % Aviscflux(:,i,j) / f % geom % jacobian(i,j)
         end do              ; end do
      else
         Avisc_flux = 0.0_RP
      end if
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL BCs(f % zone) % bc % FlowState(f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j))
      end do               ; end do

      if ( flowIsNavierStokes ) then

          useWallFuncFace = .false.
          if (useWallFunc) then
              do i = 1, len(wallFunBCs)
                  if (trim(wallFunBCs(i)) .eq. trim(f % boundaryName)) then
                      useWallFuncFace = .true.
                      exit
                  end if
              end do
          end if
          if (useWallFuncFace) then
              call WallFunctionGatherFlowVariables(mesh, f, wallFunV, wallFunRho, wallFunMu, wallFunY, wallFunVavg)
          end if

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

               mu    = f % storage(1) % mu_ns(1,i,j)
               beta  = f % storage(1) % mu_ns(3,i,j)
               kappa = f % storage(1) % mu_ns(2,i,j)

               call ViscousFlux(NCONS,NGRAD,f % storage(1) % Q(:,i,j), &
                                            f % storage(1) % U_x(:,i,j), &
                                            f % storage(1) % U_y(:,i,j), &
                                            f % storage(1) % U_z(:,i,j), &
                                            mu, beta, kappa, fv_3d)

               visc_flux(:,i,j) =   fv_3d(:,IX)*f % geom % normal(IX,i,j) &
                                  + fv_3d(:,IY)*f % geom % normal(IY,i,j) &
                                  + fv_3d(:,IZ)*f % geom % normal(IZ,i,j)

               visc_flux(:,i,j) = visc_flux(:,i,j) + Avisc_flux(:,i,j)

               if (useWallFuncFace) then
                   call WallViscousFlux(wallFunV(:,i,j), wallFunY(i,j), f % geom % normal(:,i,j), &
                                        wallFunRho(i,j), wallFunMu(i,j), wallFunVavg(:,i,j), &
                                        visc_flux(:,i,j), f % storage(1) % u_tau_NS(i,j))
               end if 

               CALL BCs(f % zone) % bc % FlowNeumann(&
                                              f % geom % x(:,i,j), &
                                              time, &
                                              f % geom % normal(:,i,j), &
                                              f % storage(1) % Q(:,i,j), &
                                              f % storage(1) % U_x(:,i,j), &
                                              f % storage(1) % U_y(:,i,j), &
                                              f % storage(1) % U_z(:,i,j), &
                                              visc_flux(:,i,j))


            end do
         end do
      else
         visc_flux = 0.0_RP

      end if

      DO j = 0, f % Nf(2)
         DO i = 0, f % Nf(1)
!
!           Hyperbolic part
!           -------------
            CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                               QRight = f % storage(2) % Q(:,i,j), &
                               nHat   = f % geom % normal(:,i,j), &
                               t1     = f % geom % t1(:,i,j), &
                               t2     = f % geom % t2(:,i,j), &
                               flux   = inv_flux)

            fStar(:,i,j) = (inv_flux - visc_flux(:,i,j) ) * f % geom % jacobian(i,j)
         END DO
      END DO

      Sidearray = (/1, HMESH_NONE/)
      call f % ProjectFluxToElements(NCONS, fStar, Sidearray)

      END SUBROUTINE computeBoundaryFlux_NSSA

end module SpatialDiscretization