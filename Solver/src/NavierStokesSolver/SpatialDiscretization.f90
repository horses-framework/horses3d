#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      use HyperbolicDiscretizations
      use EllipticDiscretizations
      use LESModels
      use ShockCapturing
      use DGIntegrals
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
                                    set_getVelocityGradients
      use ProblemFileFunctions, only: UserDefinedSourceTermNS_f
      use BoundaryConditions
      use IBMClass
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public   ComputeTimeDerivative, ComputeTimeDerivativeIsolated, viscousDiscretizationKey
      public   Initialize_SpaceAndTimeMethods, Finalize_SpaceAndTimeMethods

      procedure(GetGradientValues_f), pointer :: GetGradients
      procedure(EllipticFlux_f),      pointer :: ViscousFlux

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
               stop

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
                     stop
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
                  stop
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
                  stop

               end select

               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_NS)
               call ViscousDiscretization % Describe

            else
               if (.not. allocated(ViscousDiscretization)) allocate(EllipticDiscretization_t :: ViscousDiscretization)
               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_NS)
!
!              Set state as default option
!              ---------------------------
               call SetGradientVariables(GRADVARS_STATE)
               GetGradients => NSGradientVariables_STATE
               ViscousFlux  => ViscousFlux_STATE
               call set_getVelocityGradients(GRADVARS_STATE)

            end if

!
!           Initialize models
!           -----------------
            call InitializeLESModel(LESModel, controlVariables)
!
!           Initialize Shock-Capturing
!           --------------------------
            call Initialize_ShockCapturing(ShockCapturingDriver, controlVariables, sem, &
                                           ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
            call ShockCapturingDriver % Describe

         end if

      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_SpaceAndTimeMethods
         implicit none


         if ( allocated(HyperbolicDiscretization) ) deallocate( HyperbolicDiscretization )
         if ( allocated(LESModel) )                 deallocate( LESModel )
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
         logical :: HOElements

         if (present(HO_Elements)) then
            HOElements = HO_Elements
         else
            HOElements = .false.
         end if

         call SetBoundaryConditionsEqn(NS_BC)
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time, Level)
		 if (present(Level)) then
			call mesh % ProlongSolutionToFaces(NCONS, HO_Elements=HO_Elements, Level=Level)
		 else 
			call mesh % ProlongSolutionToFaces(NCONS, HO_Elements)
		 end if 
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
			 if (present(Level)) then
				call ViscousDiscretization % ComputeGradient( NCONS, NGRAD, mesh , time, GetGradients, HO_Elements=HO_Elements, Level=Level)
			 else 
				call ViscousDiscretization % ComputeGradient( NCONS, NGRAD, mesh , time, GetGradients, HO_Elements)
			 end if 
            
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
         if (HOElements) then
            call TimeDerivative_ComputeQDotHO(mesh = mesh , &
                                          particles = particles, &
                                          t    = time)
         else 
			 if (present(Level)) then
				call TimeDerivative_ComputeQDotMLRK(mesh = mesh , &
                                          particles = particles, &
                                          t    = time, Level=Level)
			 else 
				call TimeDerivative_ComputeQDot(mesh = mesh , &
                                          particles = particles, &
                                          t    = time)
			 end if
         end if
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
         TYPE(HexMesh), target            :: mesh
         type(Particles_t)                :: particles
         REAL(KIND=RP)                    :: time
         integer,             intent(in)  :: mode
         logical,   intent(in), optional  :: HO_Elements
         logical, intent(in), optional    :: element_mask(:)
		 integer, intent(in), optional    :: Level
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
         use TripForceClass, only: randomTrip
         use ActuatorLine, only: farm, ForcesFarm
         use SpongeClass, only: sponge, addSourceSponge
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
         integer     :: eID , i, j, k, ierr, fID, iFace, iEl, iP, STLNum, n, m  , ii, jj 
         real(kind=RP)  :: mu_smag, delta, Source(NCONS), TurbulentSource(NCONS), Q_target(NCONS)
         real(kind=RP), allocatable :: Source_HO(:,:,:,:)
         integer,       allocatable :: i_(:), j_(:), k_(:)
         real(kind=RP) :: fStarAux(1:NCONS,0:mesh%faces(1)%NfRight(1),0:mesh%faces(1)%NfRight(1))

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


         if ( LESModel % active) then
!$omp do schedule(runtime) private(i,j,k,delta,mu_smag)
            do eID = 1, size(mesh % elements)
               associate(e => mesh % elements(eID))
               delta = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call LESModel % ComputeViscosity(delta, e % geom % dWall(i,j,k), e % storage % Q(:,i,j,k),   &
                                                                                   e % storage % U_x(:,i,j,k), &
                                                                                   e % storage % U_y(:,i,j,k), &
                                                                                   e % storage % U_z(:,i,j,k), &
                                                                                   e % storage % mu_turb_NS(i,j,k) )
                                                                                   ! mu_smag)
                  ! e % storage % mu_NS(1,i,j,k) = e % storage % mu_NS(1,i,j,k) + mu_smag
                  ! e % storage % mu_NS(2,i,j,k) = e % storage % mu_NS(2,i,j,k) + mu_smag * dimensionless % mu_to_kappa
                  e % storage % mu_NS(1,i,j,k) = e % storage % mu_NS(1,i,j,k) + e % storage % mu_turb_NS(i,j,k)
                  e % storage % mu_NS(2,i,j,k) = e % storage % mu_NS(2,i,j,k) + e % storage % mu_turb_NS(i,j,k) * dimensionless % mu_to_kappa
               end do                ; end do                ; end do
               end associate
            end do
!$omp end do
      end if
!
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
!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % faces_interior)
            fID = mesh % faces_interior(iFace)
            if (mesh % faces(fID) % IsMortar==3) then 
               associate(fstar=>mesh% faces(fID)%storage(1)%fStar)
                  fstar=0.0_RP
               end associate
               associate(fstar=>mesh% faces(fID)%storage(2)%fStar)
                  fstar=0.0_RP
               end associate
            end if 
            if (mesh % faces(fID) % IsMortar==1) then 
               associate(fstar=>mesh% faces(fID)%storage(1)%fStar)
                  fstar=0.0_RP
               end associate
               associate(fstar=>mesh% faces(fID)%storage(2)%fStar)
                  fstar=0.0_RP
               end associate
               do m=1,4
                  if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                     call computeElementInterfaceFlux(fma=mesh % faces(fID), f=mesh % faces(mesh % faces(fID)%Mortar(m)), m=m)
                  end if 
               end do 
            elseif (mesh % faces(fID) % IsMortar==0) then 
               call computeElementInterfaceFlux(f=mesh % faces(fID))
   
            end if 
         end do
!$omp end do nowait

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (mesh % SlidingMesh % active) then 
!$omp do schedule(runtime) private(fID)
            do iFace = 1, size(mesh % mortar_faces)
               fID = mesh % mortar_faces(iFace)%ID
               associate(fstar=>mesh % faces(mesh % mortar_faces(fID)%Mortar(1))%storage(1)%fStar)
                  fstar=0.0_RP
               end associate
               associate(fstar=>mesh % faces(mesh % mortar_faces(fID)%Mortar(1))%storage(2)%fStar)
                  fstar=0.0_RP
               end associate
               associate(fstar=>mesh % faces(mesh % mortar_faces(fID)%Mortar(2))%storage(1)%fStar)
                  fstar=0.0_RP
               end associate
               associate(fstar=>mesh % faces(mesh % mortar_faces(fID)%Mortar(2))%storage(2)%fStar)
                  fstar=0.0_RP
               end associate
            end do  
!$omp end do nowait
         end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (mesh % SlidingMesh % active) then 
!$omp do schedule(runtime) private(fID)
            do iFace = 1, size(mesh % mortar_faces)
               fID = mesh % mortar_faces(iFace)%ID
               if (fID .ne. iFace) write(*,*)'line497 spatialdisc fID ne Iface'
               call computeElementInterfaceFlux(fma=mesh % faces(mesh % mortar_faces(fID)%Mortar(1)), fmb=mesh % faces(mesh % mortar_faces(fID)%Mortar(2)), &
               f=mesh % mortar_faces(iFace), m=m,sliding=.true.)
            end do 
!$omp end do nowait
         end if 
!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % faces_boundary)
            fID = mesh % faces_boundary(iFace)
            call computeBoundaryFlux(mesh % faces(fID), t, mesh)
         end do
!$omp end do
  !  if (mesh%sliding) then 
  !     fStarAux=0.0_RP
  !     do iFace=1, size(mesh%elements)
  !        if (mesh%elements(iFace)%sliding_newnodes) then 
  !           associate(fStar=>mesh%faces(mesh%elements(iFace)%faceIDs(5))%storage(2)%fStar)
  !             fStarAux=fStar
   !            fstar=0.0_RP
  !!           do j = 0, mesh%faces(mesh%elements(iFace)%faceIDs(5)) % NfRight(2)   ; do i = 0, mesh%faces(mesh%elements(iFace)%faceIDs(5)) % NfRight(1)   
   !             call leftIndexes2Right(i,j,mesh%faces(mesh%elements(iFace)%faceIDs(5)) % NfRight(1), mesh%faces(mesh%elements(iFace)%faceIDs(5)) % NfRight(2), &
  !!              mesh%faces(mesh%elements(iFace)%faceIDs(5)) % rotation, ii, jj)
   !             fstar(1:NCONS,ii,jj) = fStarAux(1:NCONS,i,j) 
   !          end do                        ; end do
   !          fStar=-fStar
   !          end associate


   !          associate(fStar=>mesh%faces(mesh%elements(iFace)%faceIDs(5))%storage(1)%fStar)
   !            fStarAux=fStar
   !            fstar=0.0_RP
   !          do j = 0, mesh%faces(mesh%elements(iFace)%faceIDs(5)) % NfRight(2)   ; do i = 0, mesh%faces(mesh%elements(iFace)%faceIDs(5)) % NfRight(1)   
   !             call leftIndexes2Right(i,j,mesh%faces(mesh%elements(iFace)%faceIDs(5)) % NfRight(1), mesh%faces(mesh%elements(iFace)%faceIDs(5)) % NfRight(2), &
   !             mesh%faces(mesh%elements(iFace)%faceIDs(5)) % rotation, ii, jj)
   !             fstar(1:NCONS,ii,jj) = fStarAux(1:NCONS,i,j) 
   !          end do                        ; end do
   !          fStar=-fStar
   !          end associate
   !       end if 
   !    end do
   ! end if 

         !call mesh % faces(1)%TestMortar(mesh%faces(mesh%elements(387)%faceIDs(5)), mesh%faces(mesh%elements(440)%faceIDs(5)),NCONS)

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
!        Wait until messages are sent
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
!$omp do schedule(runtime) private(fID)
            do iFace = 1, size(mesh % faces_mpi)
               fID = mesh % faces_mpi(iFace)
               if (mesh% faces(fID)%IsMortar==1) then 
                  !write(*,*) 'big mortar face mpi'
                  associate(fstar=>mesh% faces(fID)%storage(1)%fStar)
                     fstar=0.0_RP
                  end associate
                  do m=1,4
                     if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                        !write(*,*) mesh % faces(fID)%Mortar(m), mesh % faces(mesh % faces(fID)%Mortar(m))%IsMortar
                        call computeElementInterfaceFlux(fma=mesh % faces(fID), f=mesh % faces(mesh % faces(fID)%Mortar(m)))
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
! !
! !           Add an MPI Barrier
! !           ------------------
! !$omp single
!             call mpi_barrier(MPI_COMM_WORLD, ierr)
! !$omp end single
         end if
#endif
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
               ! All terms are calculated indepentenly and overwritten in case one gauss point has more than one contribution
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
                  call randomTrip % getTripSource( e % geom % x(:,i,j,k), e % storage % S_NS(:,i,j,k) )
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do
            ! for the sponge, loops are in the internal subroutine as values are precalculated
            call addSourceSponge(sponge,mesh)
            call ForcesFarm(farm, mesh, t)
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
!
!        ***********************
!        Now add the source term
!        ***********************
!$omp do schedule(runtime) private(i,j,k)
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
!$omp do schedule(runtime) private(i,j,k,Source,Q_target)
                  do eID = 1, mesh % no_of_elements  
                     associate ( e => mesh % elements(eID) ) 
                     do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                        if( e% isInsideBody(i,j,k) ) then
                           if( mesh% IBM% stl(e% STL(i,j,k))% move ) then 
                              Q_target = mesh% IBM% MaskVelocity( e% storage% Q(:,i,j,k), NCONS, e% STL(i,j,k), e% geom% x(:,i,j,k), t )
                              call mesh% IBM% SourceTerm( eID = eID, Q = e % storage % Q(:,i,j,k), Q_target = Q_target, Source = Source, wallfunction = .false. )
                           else 
                              call mesh% IBM% SourceTerm( eID = eID, Q = e % storage % Q(:,i,j,k), Source = Source, wallfunction = .false. )
                           end if 
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
                     associate( e    => mesh% elements(mesh% IBM% ImagePoints(iP)% element_index), &
                                e_in => mesh% elements(mesh% IBM% ImagePoints(iP)% element_in)     )
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        ---------------
!     Compute Qdot for MLRK Scheme - This split is  necessary due to problem on AnisFAS time integration
!        ---------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine TimeDerivative_ComputeQDotMLRK( mesh , particles, t, Level)
      use WallFunctionConnectivity
         use TripForceClass, only: randomTrip
         use ActuatorLine, only: farm, ForcesFarm
         use SpongeClass, only: sponge, addSourceSponge
         implicit none
         type(HexMesh)              :: mesh
         type(Particles_t)          :: particles
         real(kind=RP)              :: t
		 integer, intent(in), optional   :: Level
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID, iFace, iEl, iP, STLNum, n , locLevel, lID
         real(kind=RP)  :: mu_smag, delta, Source(NCONS), TurbulentSource(NCONS), Q_target(NCONS)
         real(kind=RP), allocatable :: Source_HO(:,:,:,:)
         integer,       allocatable :: i_(:), j_(:), k_(:)
		 
		 if (present(Level)) then
            locLevel = Level
         else
            locLevel = 1
         end if
		 
		 associate ( MLRK => mesh % MLRK)
	     associate ( MLIter_eID      => MLRK % MLIter_eID, &
                     MLIter          => MLRK % MLIter,      &
                     MLIter_eID_Seq  => MLRK % MLIter_eID_Seq, &
                     MLIter_eID_MPI  => MLRK % MLIter_eID_MPI,  &
					 MLIter_fID      => MLRK % MLIter_fID,  &
					 MLIter_fID_Interior => MLRK %  MLIter_fID_Interior, &
					 MLIter_fID_Boundary => MLRK %  MLIter_fID_Boundary, & 
					 MLIter_fID_MPI      => MLRK %  MLIter_fID_MPI )
!
!        ***********************************************
!        Compute the viscosity at the elements and faces
!        ***********************************************
!
         if (flowIsNavierStokes) then
!$omp do schedule(runtime) private(i,j,k,eID)
            do lID = 1, MLIter(locLevel,1)
			   eID = MLIter_eID(lID)
               associate(e => mesh % elements(eID))
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call get_laminar_mu_kappa(e % storage % Q(:,i,j,k), e % storage % mu_NS(1,i,j,k), e % storage % mu_NS(2,i,j,k))
               end do                ; end do                ; end do
               end associate
            end do
!$omp end do
         end if


         if ( LESModel % active) then
!$omp do schedule(runtime) private(i,j,k,delta,mu_smag,eID)
            do lID = 1, MLIter(locLevel,1)
			   eID = MLIter_eID(lID)
               associate(e => mesh % elements(eID))
               delta = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call LESModel % ComputeViscosity(delta, e % geom % dWall(i,j,k), e % storage % Q(:,i,j,k),   &
                                                                                   e % storage % U_x(:,i,j,k), &
                                                                                   e % storage % U_y(:,i,j,k), &
                                                                                   e % storage % U_z(:,i,j,k), &
                                                                                   e % storage % mu_turb_NS(i,j,k) )
                                                                                   ! mu_smag)
                  ! e % storage % mu_NS(1,i,j,k) = e % storage % mu_NS(1,i,j,k) + mu_smag
                  ! e % storage % mu_NS(2,i,j,k) = e % storage % mu_NS(2,i,j,k) + mu_smag * dimensionless % mu_to_kappa
                  e % storage % mu_NS(1,i,j,k) = e % storage % mu_NS(1,i,j,k) + e % storage % mu_turb_NS(i,j,k)
                  e % storage % mu_NS(2,i,j,k) = e % storage % mu_NS(2,i,j,k) + e % storage % mu_turb_NS(i,j,k) * dimensionless % mu_to_kappa
               end do                ; end do                ; end do
               end associate
            end do
!$omp end do
      end if
!
!        Compute viscosity at interior and boundary faces
!        ------------------------------------------------
         call compute_viscosity_at_faces(MLIter(locLevel,3), 2, MLIter_fID_Interior(1:MLIter(locLevel,3)), mesh)
         call compute_viscosity_at_faces(MLIter(locLevel,4), 1, MLIter_fID_Boundary(1:MLIter(locLevel,4)), mesh)
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) private(eID)
         do lID = 1, MLIter(locLevel,1)
			eID = MLIter_eID(lID)
            call TimeDerivative_VolumetricContribution( mesh, mesh % elements(eID) , t)
         end do
!$omp end do

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
!$omp do schedule(runtime) private(fID)
         do iFace = 1, MLIter(locLevel,3)
            fID = MLIter_fID_Interior(iFace)
            call computeElementInterfaceFlux(mesh % faces(fID))
         end do
!$omp end do nowait

!$omp do schedule(runtime) private(fID)
         do iFace = 1, MLIter(locLevel,4)
            fID = MLIter_fID_Boundary(iFace)
            call computeBoundaryFlux(mesh % faces(fID), t, mesh)
         end do
!$omp end do
!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
!
!$omp do schedule(runtime) private(i,j,k,eID)
         do iEl = 1, MLIter(locLevel,5)
            eID = MLIter_eID_Seq(iEl)
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
!        Wait until messages are sent
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
            call compute_viscosity_at_faces(MLIter(locLevel,7), 2, MLIter_fID_MPI(1:MLIter(locLevel,7)), mesh)

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
!$omp do schedule(runtime) private(fID)
            do iFace = 1, MLIter(locLevel,7)
               fID = MLIter_fID_MPI(iFace)
               call computeMPIFaceFlux(mesh % faces(fID))
            end do
!$omp end do
!
!           ***********************************************************
!           Surface integrals and scaling of elements with shared faces
!           ***********************************************************
!
!$omp do schedule(runtime) private(i,j,k,eID)
            do iEl = 1, MLIter(locLevel,6)
               eID = MLIter_eID_MPI(iEl)
               associate(e => mesh % elements(eID))
               call TimeDerivative_FacesContribution(e, t, mesh)

               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k)
               end do         ; end do          ; end do
               end associate
            end do
!$omp end do
! !
! !           Add an MPI Barrier
! !           ------------------
! !$omp single
!             call mpi_barrier(MPI_COMM_WORLD, ierr)
! !$omp end single
         end if
#endif
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
!$omp do schedule(runtime) private(i,j,k,eID)
            do lID = 1, MLIter(locLevel,1)
			   eID = MLIter_eID(lID)
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               ! All terms are calculated indepentenly and overwritten in case one gauss point has more than one contribution
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
                  call randomTrip % getTripSource( e % geom % x(:,i,j,k), e % storage % S_NS(:,i,j,k) )
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do
            ! for the sponge, loops are in the internal subroutine as values are precalculated
            call addSourceSponge(sponge,mesh)
            call ForcesFarm(farm, mesh, t, Level=locLevel)
!
!           Add Particles source
!           ********************
            if ( particles % active ) then
!$omp do schedule(runtime) private(eID)
               do lID = 1, MLIter(locLevel,1)
				  eID = MLIter_eID(lID)
                  associate ( e => mesh % elements(eID) )
                     e % storage % S_NSP = 0.0_RP
                  end associate
               end do
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
!$omp do schedule(runtime) private(eID)
               do lID = 1, MLIter(locLevel,1)
				  eID = MLIter_eID(lID)
                  associate ( e => mesh % elements(eID) )
                     e % storage % S_NS = e % storage % S_NS + e % storage % S_NSP
                  end associate
               end do
!$omp end do
            end if
         end if !(.not. mesh % child)
!
!        ***********************
!        Now add the source term
!        ***********************
!$omp do schedule(runtime) private(i,j,k,eID)
		 do lID = 1, MLIter(locLevel,1)
		    eID = MLIter_eID(lID)
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
!$omp do schedule(runtime) private(i,j,k,Source,Q_target,eID)
                  do lID = 1, MLIter(locLevel,1)
					 eID = MLIter_eID(lID)
                     associate ( e => mesh % elements(eID) ) 
                     do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                        if( e% isInsideBody(i,j,k) ) then
                           if( mesh% IBM% stl(e% STL(i,j,k))% move ) then 
                              Q_target = mesh% IBM% MaskVelocity( e% storage% Q(:,i,j,k), NCONS, e% STL(i,j,k), e% geom% x(:,i,j,k), t )
                              call mesh% IBM% SourceTerm( eID = eID, Q = e % storage % Q(:,i,j,k), Q_target = Q_target, Source = Source, wallfunction = .false. )
                           else 
                              call mesh% IBM% SourceTerm( eID = eID, Q = e % storage % Q(:,i,j,k), Source = Source, wallfunction = .false. )
                           end if 
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
                     associate( e    => mesh% elements(mesh% IBM% ImagePoints(iP)% element_index), &
                                e_in => mesh% elements(mesh% IBM% ImagePoints(iP)% element_in)     )
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
		 
		 end associate
		 end associate

      end subroutine TimeDerivative_ComputeQDotMLRK

!*******************************************************************************
!     Computes the Q-dot only for the High-Order elements (p>1)
!     Required for the hybrid Euler RK3 temporal scheme
!*******************************************************************************
   subroutine TimeDerivative_ComputeQDotHO( mesh , particles, t)
      use WallFunctionConnectivity
         use TripForceClass, only: randomTrip
         use ActuatorLine, only: farm, ForcesFarm
         use SpongeClass, only: sponge, addSourceSponge
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
         integer     :: eID , i, j, k, ierr, fID, iFace, iEl, iP, STLNum, n 
         real(kind=RP)  :: mu_smag, delta, Source(NCONS), TurbulentSource(NCONS), Q_target(NCONS)
         real(kind=RP), allocatable :: Source_HO(:,:,:,:)
         integer,       allocatable :: i_(:), j_(:), k_(:)
!
!        ***********************************************
!        Compute the viscosity at the elements and faces
!        ***********************************************
!
         write(*,*)'HO element spatial discretization line 774'
         if (flowIsNavierStokes) then
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, size(mesh % HO_Elements)
               associate(e => mesh % elements(mesh % HO_Elements(eID)))
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call get_laminar_mu_kappa(e % storage % Q(:,i,j,k), e % storage % mu_NS(1,i,j,k), e % storage % mu_NS(2,i,j,k))
               end do                ; end do                ; end do
               end associate
            end do
!$omp end do
         end if


         if ( LESModel % active) then
!$omp do schedule(runtime) private(i,j,k,delta,mu_smag)
            do eID = 1, size(mesh % HO_Elements)
               associate(e => mesh % elements(mesh % HO_Elements(eID)))
               delta = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call LESModel % ComputeViscosity(delta, e % geom % dWall(i,j,k), e % storage % Q(:,i,j,k),   &
                                                                                    e % storage % U_x(:,i,j,k), &
                                                                                    e % storage % U_y(:,i,j,k), &
                                                                                    e % storage % U_z(:,i,j,k), &
                                                                                    e % storage % mu_turb_NS(i,j,k) )
                                                                                    ! mu_smag)
                  ! e % storage % mu_NS(1,i,j,k) = e % storage % mu_NS(1,i,j,k) + mu_smag
                  ! e % storage % mu_NS(2,i,j,k) = e % storage % mu_NS(2,i,j,k) + mu_smag * dimensionless % mu_to_kappa
                  e % storage % mu_NS(1,i,j,k) = e % storage % mu_NS(1,i,j,k) + e % storage % mu_turb_NS(i,j,k)
                  e % storage % mu_NS(2,i,j,k) = e % storage % mu_NS(2,i,j,k) + e % storage % mu_turb_NS(i,j,k) * dimensionless % mu_to_kappa
               end do                ; end do                ; end do
               end associate
            end do
!$omp end do
      end if
!
!        Compute viscosity at interior and boundary faces
!        ------------------------------------------------
      call compute_viscosity_at_faces(size(mesh % HO_FacesInterior), 2, mesh % HO_FacesInterior, mesh)
      call compute_viscosity_at_faces(size(mesh % faces_boundary), 1, mesh % faces_boundary, mesh)

!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) private(eID)
      do iEl = 1 , size(mesh % HO_Elements)
         eID = mesh % HO_Elements(iEl)
         call TimeDerivative_VolumetricContribution( mesh, mesh % elements(eID) , t)
      end do
!$omp end do

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
!$omp do schedule(runtime) private(fID)
      do iFace = 1, size(mesh % HO_FacesInterior)
         fID = mesh % HO_FacesInterior(iFace)
         write(*,*)'HOOOOOO L842'
         call computeElementInterfaceFlux(mesh % faces(fID))
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
      do iEl = 1, size(mesh % HO_ElementsSequential)
         eID = mesh % HO_ElementsSequential(iEl)
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
!        Wait until messages are sent
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
!$omp do schedule(runtime) private(fID)
            do iFace = 1, size(mesh % faces_mpi)
               fID = mesh % faces_mpi(iFace)
               call computeMPIFaceFlux(mesh % faces(fID))
            end do
!$omp end do
!
!           ***********************************************************
!           Surface integrals and scaling of elements with shared faces
!           ***********************************************************
!
!$omp do schedule(runtime) private(i,j,k,eID)
            do iEl = 1, size(mesh % HO_ElementsMPI)
               eID = mesh % HO_ElementsMPI(iEl)
               associate(e => mesh % elements(eID))
               call TimeDerivative_FacesContribution(e, t, mesh)

               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k)
               end do         ; end do          ; end do
               end associate
            end do
!$omp end do
! !
! !           Add an MPI Barrier
! !           ------------------
! !$omp single
!             call mpi_barrier(MPI_COMM_WORLD, ierr)
! !$omp end single
         end if
#endif
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
            do eID = 1, size(mesh % HO_Elements)
               associate ( e => mesh % elements(mesh % HO_Elements(eID)) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
                  call randomTrip % getTripSource( e % geom % x(:,i,j,k), e % storage % S_NS(:,i,j,k) )
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do

            call addSourceSponge(sponge,mesh)
            call ForcesFarm(farm, mesh, t)
!
!           Add Particles source
!           ********************
            if ( particles % active ) then
!$omp do schedule(runtime)
               do eID = 1, size(mesh % HO_Elements)
                  associate ( e => mesh % elements(mesh % HO_Elements(eID)) )
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
               do eID = 1, size(mesh % HO_Elements)
                  associate ( e => mesh % elements(mesh % HO_Elements(eID)) )
                     e % storage % S_NS = e % storage % S_NS + e % storage % S_NSP
                  end associate
               enddo
!$omp end do
            end if ! end particles
         end if !(.not. mesh % child)
!
!        ***********************
!        Now add the source term
!        ***********************
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, size(mesh % HO_Elements)
            associate ( e => mesh % elements(mesh % HO_Elements(eID)) )
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
!$omp do schedule(runtime) private(i,j,k,Source,Q_target)
               do eID = 1, size(mesh % HO_Elements)
                  
                  associate ( e => mesh % elements(mesh % HO_Elements(eID)) )

                  do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                     if( e% isInsideBody(i,j,k) ) then
                        if( mesh% IBM% stl(e% STL(i,j,k))% move ) then 
                           Q_target = mesh% IBM% MaskVelocity( e% storage% Q(:,i,j,k), NCONS, e% STL(i,j,k), e% geom% x(:,i,j,k), t )
                           call mesh% IBM% SourceTerm( eID = eID, Q = e % storage % Q(:,i,j,k), Q_target = Q_target, Source = Source, wallfunction = .false. )
                        else 
                           call mesh% IBM% SourceTerm( eID = eID, Q = e % storage % Q(:,i,j,k), Source = Source, wallfunction = .false. )
                        end if 
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
                     associate( e    => mesh% elements(mesh% IBM% ImagePoints(iP)% element_index), &
                                 e_in => mesh% elements(mesh% IBM% ImagePoints(iP)% element_in)     )
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

      end subroutine TimeDerivative_ComputeQDotHO

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
         real(kind=RP) :: delta, mu_smag

         if (flowIsNavierStokes) then
!$omp do schedule(runtime) private(i,j)
            do iFace = 1, no_of_faces
               associate(f => mesh % faces(face_ids(iFace)))
                  if (f % IsMortar==1 .OR. f % IsMortar==3) cycle  
               do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
                  do side = 1, no_of_sides
                     !write(*,*) 'Q(1)',f % storage(side) % Q(1,i,j)
                      call get_laminar_mu_kappa(f % storage(side) % Q(:,i,j), f % storage(side) % mu_NS(1,i,j), f % storage(side) % mu_NS(2,i,j))
                  end do
               end do              ; end do
               end associate
            end do
!$omp end do
         end if

         if (mesh % SlidingMesh % active) then 
!$omp do schedule(runtime) private(i,j)
            do iFace = 1, size(mesh%mortar_faces)
               associate(f => mesh % mortar_faces(iFace))  
               do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
                  do side = 1, no_of_sides
                     !write(*,*) 'Q1 mortar_faces', f % storage(side) % Q(1,i,j)
                     !write(*,*) 'Q2 mortar_faces', f % storage(side) % Q(2,i,j)
                     !write(*,*) 'Q3 mortar_faces', f % storage(side) % Q(3,i,j)
                     !write(*,*) 'Q4 mortar_faces', f % storage(side) % Q(4,i,j)
                     !write(*,*) 'Q5 mortar_faces', f % storage(side) % Q(5,i,j)
                      call get_laminar_mu_kappa(f % storage(side) % Q(:,i,j), f % storage(side) % mu_NS(1,i,j), f % storage(side) % mu_NS(2,i,j))
                  end do
               end do              ; end do
               end associate
            end do
!$omp end do
         end if 
         if ( LESModel % Active ) then
!$omp do schedule(runtime) private(i,j,delta,mu_smag)
            do iFace = 1, no_of_faces
               associate(f => mesh % faces(face_ids(iFace)))

               delta = sqrt(f % geom % surface / product(f % Nf + 1))
               if (f % IsMortar==1 .OR. f % IsMortar==3) cycle 
               do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
                  do side = 1, no_of_sides
                     call LESModel % ComputeViscosity(delta, f % geom % dWall(i,j), f % storage(side) % Q(:,i,j),   &
                                                                                    f % storage(side) % U_x(:,i,j), &
                                                                                    f % storage(side) % U_y(:,i,j), &
                                                                                    f % storage(side) % U_z(:,i,j), &
                                                                                    mu_smag)
                     f % storage(side) % mu_NS(1,i,j) = f % storage(side) % mu_NS(1,i,j) + mu_smag
                     f % storage(side) % mu_NS(2,i,j) = f % storage(side) % mu_NS(2,i,j) + mu_smag * dimensionless % mu_to_kappa
                  end do
               end do              ; end do
               end associate
            end do
!$omp end do
         end if

         if ( LESModel % Active .and. mesh % SlidingMesh % active) then
!$omp do schedule(runtime) private(i,j,delta,mu_smag)
                        do iFace = 1, size(mesh%mortar_faces)
                           associate(f => mesh % mortar_faces(iFace))
            
                           delta = sqrt(f % geom % surface / product(f % Nf + 1))
                           do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
                              do side = 1, no_of_sides
                                 call LESModel % ComputeViscosity(delta, f % geom % dWall(i,j), f % storage(side) % Q(:,i,j),   &
                                                                                                f % storage(side) % U_x(:,i,j), &
                                                                                                f % storage(side) % U_y(:,i,j), &
                                                                                                f % storage(side) % U_z(:,i,j), &
                                                                                                mu_smag)
                                 f % storage(side) % mu_NS(1,i,j) = f % storage(side) % mu_NS(1,i,j) + mu_smag
                                 f % storage(side) % mu_NS(2,i,j) = f % storage(side) % mu_NS(2,i,j) + mu_smag * dimensionless % mu_to_kappa
                              end do
                           end do              ; end do
                           end associate
                        end do
!$omp end do
                     end if


      end subroutine compute_viscosity_at_faces
!
!////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------------
!     This routine computes Qdot neglecting the interaction with neighboring elements
!     and boundaries. Therefore, the external states are not needed.
!     -------------------------------------------------------------------------------
      subroutine TimeDerivative_ComputeQDotIsolated( mesh , t )
         use TripForceClass, only: randomTrip
         use ActuatorLine, only: farm, ForcesFarm
         use SpongeClass, only: sponge, addSourceSponge
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
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
                  call randomTrip % getTripSource( e % geom % x(:,i,j,k), e % storage % S_NS(:,i,j,k) )
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do
            call ForcesFarm(farm, mesh, t)
            call addSourceSponge(sponge,mesh)
         end if

!$omp do schedule(runtime) private(i,j,k)
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
         endif

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
         real(kind=RP) :: viscousContravariantFlux ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP) :: AviscContravariantFlux   ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP) :: contravariantFlux        ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )

!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
          if (mesh%slidingflux .and. e%sliding) then 
!         call SlidingMeshFluxCalculation(e, SlidingMeshFlux)
          call contravariantSMFlux ( e , inviscidContravariantFlux )
            ! print*, SlidingMeshFlux/inviscidContravariantFlux 
!            inviscidContravariantFlux = inviscidContravariantFlux - SlidingMeshFlux
         else 
            call HyperbolicDiscretization % ComputeInnerFluxes ( e , EulerFlux, inviscidContravariantFlux )
         end if 
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
!
!           Compute sharp fluxes for skew-symmetric approximations
!           ------------------------------------------------------
            call HyperbolicDiscretization % ComputeSplitFormFluxes(e, inviscidContravariantFlux, fSharp, gSharp, hSharp)
!
!           Perform the Weak volume green integral
!           --------------------------------------
            viscousContravariantFlux = viscousContravariantFlux + AviscContravariantFlux

            e % storage % QDot = -ScalarStrongIntegrals % SplitVolumeDivergence( e, fSharp, gSharp, hSharp, viscousContravariantFlux)

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
         !real(kind=RP) :: SlidingMeshFlux           ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         integer :: i,j,k

!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
!         call HyperbolicDiscretization % ComputeInnerFluxes ( e , EulerFlux, inviscidContravariantFlux )
         if (mesh%slidingflux .and. e%sliding) then 
!         call SlidingMeshFluxCalculation(e, SlidingMeshFlux)
          call contravariantSMFlux ( e , inviscidContravariantFlux )
            ! print*, SlidingMeshFlux/inviscidContravariantFlux 
!            inviscidContravariantFlux = inviscidContravariantFlux - SlidingMeshFlux
         else 
            call HyperbolicDiscretization % ComputeInnerFluxes ( e , EulerFlux, inviscidContravariantFlux )
         end if 
        
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
      subroutine computeElementInterfaceFlux(f, fma, fmb, m, sliding )
        use FaceClass
        use RiemannSolvers_NS
        implicit none
        type(Face)   , intent(inout) :: f
        type(Face), optional, intent(inout) :: fma 
        type(Face), optional, intent(inout) :: fmb
        integer, optional, intent(in) :: m 
        logical, optional , intent(in) :: sliding
        !type(Face), optional, intent(inout) :: fmb 
        !type(Face), optional, intent(inout) :: fmc 
        !type(Face), optional, intent(inout) :: fmd 

        integer       :: i, j
        real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
        real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
        real(kind=RP) :: Avisc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
        real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
        real(kind=RP) :: mu_left(3), mu_right(3)
        integer       :: Sidearray(2)



        !if (f % IsMortar == 0 .OR. f % IsMortar==2) then 
      !if (f % IsMortar == 0 ) then 
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
  !
          ! write(*,*)'face',f%ID,'h in line 1511', f%geom%h,'is mortar',f%ismortar
           if (flowIsNavierStokes) then
              do j = 0, f % Nf(2)
                 do i = 0, f % Nf(1)

                    mu_left(1) = f % storage(1) % mu_NS(1,i,j)
                    mu_left(2) = 0.0_RP
                    mu_left(3) = f % storage(1) % mu_NS(2,i,j)

                    mu_right(1) = f % storage(2) % mu_NS(1,i,j)
                    mu_right(2) = 0.0_RP
                    mu_right(3) = f % storage(2) % mu_NS(2,i,j)

                    call ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
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
           else
              visc_flux = 0.0_RP
           end if

           do j = 0, f % Nf(2)
              do i = 0, f % Nf(1)
  !
  !              --------------
  !              Invscid fluxes
  !              --------------
  !   
               
                 call RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                   QRight = f % storage(2) % Q(:,i,j), &
                                   nHat   = f % geom % normal(:,i,j), &
                                   t1     = f % geom % t1(:,i,j), &
                                   t2     = f % geom % t2(:,i,j), &
                                   flux   = inv_flux(:,i,j) )
  !
  !              Multiply by the Jacobian
  !              ------------------------
                 flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j) - Avisc_flux(:,i,j)

              end do
           end do
  !
  !        ---------------------------
  !        Return the flux to elements
  !        ---------------------------
  !
      if (.not.present(sliding)) then     
         if (f % IsMortar==0) then 
            Sidearray = (/1,2/)
            call f % ProjectFluxToElements(NCONS, flux, Sidearray)
         end if 
         if (f % IsMortar==2 .and. present(fma)) then 
            Sidearray = (/1,0/)
            call fma % ProjectMortarFluxToElements(nEqn=NCONS, whichElements=Sidearray, &
               fma=f, MortarFlux=flux)
               Sidearray = (/0,2/)
               call f % ProjectFluxToElements(NCONS, flux, Sidearray)
         end if 
      else 
         !write(*,*) 'projecting flux of mortr',f%ID,'to faces fma',fma%ID,'and fmb', fmb%ID
         !write(*,*)'element of face',fma%ID,'=',fma%elementIDs(1)
         !write(*,*)'element of face',fmb%ID,'=',fmb%elementIDs(1)
         Sidearray = (/1,0/)
         call fma % ProjectMortarFluxToElements(nEqn=NCONS, whichElements=Sidearray, &
         fma=f, MortarFlux=flux, sliding= .true.) 
         Sidearray = (/2,0/)
         call fmb % ProjectMortarFluxToElements(nEqn=NCONS, whichElements=Sidearray, &
         fma=f, MortarFlux=flux,  sliding=.true.) 
     end if 

     end subroutine computeElementInterfaceFlux

      subroutine computeMPIFaceFlux(f)
         use FaceClass
         use RiemannSolvers_NS
         implicit none
         type(Face)   , intent(inout) :: f
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
!
         if (flowIsNavierStokes) then
            do j = 0, f % Nf(2)
               do i = 0, f % Nf(1)

                  mu_left(1) = f % storage(1) % mu_NS(1,i,j)
                  mu_left(2) = 0.0_RP
                  mu_left(3) = f % storage(1) % mu_NS(2,i,j)

                  mu_right(1) = f % storage(2) % mu_NS(1,i,j)
                  mu_right(2) = 0.0_RP
                  mu_right(3) = f % storage(2) % mu_NS(2,i,j)

                  call ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
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
         else
            visc_flux = 0.0_RP
         end if

         do j = 0, f % Nf(2)
            do i = 0, f % Nf(1)
!
!              --------------
!              Invscid fluxes
!              --------------
!
               call RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j) - Avisc_flux(:,i,j)

            end do
         end do
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)

         Sidearray = (/thisSide, HMESH_NONE/)
         call f % ProjectFluxToElements(NCONS, flux, Sidearray )
         if (f % IsMortar==2) then 
            !write(*,*) 'this side', thisSide
            call f% Interpolatesmall2big(NCONS, flux)

         end if 
      end subroutine ComputeMPIFaceFlux

      SUBROUTINE computeBoundaryFlux(f, time, mesh)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_NS
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

      if (f % IsMortar .ne. 0) then 
         write(*,*) 'bface problem mortar...'
      end if 
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

         do j = 0, f % Nf(2)
            do i = 0, f % Nf(1)
               mu    = f % storage(1) % mu_ns(1,i,j)
               beta  = 0.0_RP
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

      do j = 0, f % Nf(2)
         do i = 0, f % Nf(1)
!
!           Hyperbolic part
!           -------------
            !write(*,*) 'riemann solver for boudnary face', f% ID, 'element', f % elementIDs(1), f % elementIDs(2)
           !! write(*,*) 'qleft=',f % storage(1) % Q(:,i,j)
            !write(*,*) 'qright=',f % storage(2) % Q(:,i,j)
            !write(*,*) '///////////////////////////////////'
            call RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                               QRight = f % storage(2) % Q(:,i,j), &
                               nHat   = f % geom % normal(:,i,j), &
                               t1     = f % geom % t1(:,i,j), &
                               t2     = f % geom % t2(:,i,j), &
                               flux   = inv_flux)

            fStar(:,i,j) = (inv_flux - visc_flux(:,i,j) ) * f % geom % jacobian(i,j)
         end do
      end do

      Sidearray = (/1, HMESH_NONE/)
      call f % ProjectFluxToElements(NCONS, fStar, Sidearray)

      end subroutine computeBoundaryFlux
      subroutine contravariantSMFlux(e, contravariantFlux)
         implicit none
         type(element), intent(in)  :: e
         real(kind=RP), intent(out) :: contravariantFlux(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP)              :: omega(IX:IZ), u_mesh, v_mesh, w_mesh, u, v, w, p
         real(kind=RP)              :: F(1:NCONS,IX:IZ)
         
         integer :: i,j,k
            ! define angular velocity
            omega(1) = 0.0_RP
            omega(2) = PI/20.0_RP
      !     omega(2) = 0.0_RP
            omega(3) = 0.0_RP
            ! loop inside the element 
         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2)    ; do i = 0, e%Nxyz(1)
         ! compute the velocity of the mesh 
            u_mesh = (e % geom % x(IZ,i,j,k)*omega(2)-omega(3)*e % geom % x(IY,i,j,k))
            v_mesh = (e % geom % x(IX,i,j,k)*omega(3)-omega(1)*e % geom % x(IZ,i,j,k))
            w_mesh = (e % geom % x(IY,i,j,k)*omega(1)-omega(2)*e % geom % x(IX,i,j,k))

            u = e % storage % Q(IRHOU,i,j,k) / e % storage % Q(IRHO,i,j,k)
            v = e % storage % Q(IRHOV,i,j,k) / e % storage % Q(IRHO,i,j,k)
            w = e % storage % Q(IRHOW,i,j,k) / e % storage % Q(IRHO,i,j,k)
            p = thermodynamics % gammaMinus1 * (e % storage % Q(IRHOE,i,j,k) - 0.5_RP * ( e % storage % Q(IRHOU,i,j,k) * u + e % storage % Q(IRHOV,i,j,k) * v + e % storage % Q(IRHOW,i,j,k) * w ) )
   !
   !        X-Flux
   !        ------         
            F(IRHO , IX ) = e % storage % Q(IRHOU,i,j,k) - (e % storage % Q(IRHO,i,j,k)  * u_mesh)
            F(IRHOU, IX ) = e % storage % Q(IRHOU,i,j,k) * u + p - (e % storage % Q(IRHOU,i,j,k) * u_mesh)
            F(IRHOV, IX ) = e % storage % Q(IRHOU,i,j,k) * v - (e % storage % Q(IRHOV,i,j,k) * u_mesh)
            F(IRHOW, IX ) = e % storage % Q(IRHOU,i,j,k) * w - (e % storage % Q(IRHOW,i,j,k) * u_mesh)
            F(IRHOE, IX ) = ( e % storage % Q(IRHOE,i,j,k) + p ) * u - (e % storage % Q(IRHOE,i,j,k) * u_mesh)
   !
   !        Y-Flux
   !        ------
            F(IRHO , IY ) = e % storage % Q(IRHOV,i,j,k) - (e % storage % Q(IRHO,i,j,k)  * v_mesh)
            F(IRHOU ,IY ) = e % storage % Q(IRHOV,i,j,k) * u - (e % storage % Q(IRHOU,i,j,k) * v_mesh)
            F(IRHOV ,IY ) = e % storage % Q(IRHOV,i,j,k) * v + p - (e % storage % Q(IRHOV,i,j,k) * v_mesh)
            F(IRHOW ,IY ) = e % storage % Q(IRHOV,i,j,k) * w - (e % storage % Q(IRHOW,i,j,k) * v_mesh)
            F(IRHOE ,IY ) = ( e % storage % Q(IRHOE,i,j,k) + p ) * v - (e % storage % Q(IRHOE,i,j,k) * v_mesh)
   !
   !        Z-Flux
   !        ------
            F(IRHO ,IZ) = e % storage % Q(IRHOW,i,j,k) - (e % storage % Q(IRHO,i,j,k)  * w_mesh)
            F(IRHOU,IZ) = e % storage % Q(IRHOW,i,j,k) * u - (e % storage % Q(IRHOU,i,j,k) * w_mesh)
            F(IRHOV,IZ) = e % storage % Q(IRHOW,i,j,k) * v - (e % storage % Q(IRHOV,i,j,k) * w_mesh)
            F(IRHOW,IZ) = e % storage % Q(IRHOW,i,j,k) * w + p - (e % storage % Q(IRHOW,i,j,k) * w_mesh)
            F(IRHOE,IZ) = ( e % storage % Q(IRHOE,i,j,k) + p ) * w - (e % storage % Q(IRHOE,i,j,k) * w_mesh)

            
            contravariantFlux(:,i,j,k,IX) =    F(:,IX) * e % geom % jGradXi(IX,i,j,k)   &
                                             + F(:,IY) * e % geom % jGradXi(IY,i,j,k)   &
                                             + F(:,IZ) * e % geom % jGradXi(IZ,i,j,k)
   
            contravariantFlux(:,i,j,k,IY) =    F(:,IX) * e % geom % jGradEta(IX,i,j,k)  &
                                             + F(:,IY) * e % geom % jGradEta(IY,i,j,k)  &
                                             + F(:,IZ) * e % geom % jGradEta(IZ,i,j,k)
   
            contravariantFlux(:,i,j,k,IZ) =    F(:,IX) * e % geom % jGradZeta(IX,i,j,k) &
                                             + F(:,IY) * e % geom % jGradZeta(IY,i,j,k) &
                                             + F(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)

         enddo;enddo;enddo

        end subroutine contravariantSMFlux

        subroutine SlidingMeshFluxCalculation(e,  contravariantFlux)
         type(element), intent(in)  :: e
         real(kind=RP), intent(out) :: contravariantFlux(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP)              :: u_g, v_g, w_g, omega(3), F(NCONS,NDIM)
         integer :: i,j,k
         ! define angular velocity
         omega(1) = 0.0_RP
         omega(2) = PI/20.0_RP
    !     omega(2) = 0.0_RP
         omega(3) = 0.0_RP
         ! loop inside the element 
         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2)    ; do i = 0, e%Nxyz(1)
         ! compute the velocity of the grid (adimensional)
            u_g = (e % geom % x(IZ,i,j,k)*omega(2)-omega(3)*e % geom % x(IY,i,j,k))!/refValues%V
            v_g = (e % geom % x(IX,i,j,k)*omega(3)-omega(1)*e % geom % x(IZ,i,j,k))!/refValues%V
            w_g = (e % geom % x(IY,i,j,k)*omega(1)-omega(2)*e % geom % x(IX,i,j,k))!/refValues%V
         ! compute cartesian fluxes due to mesh movement
         ! print*, i,j,k,e % geom % x(:,i,j,k) , omega , u_g,v_g,w_g
         !  X-Flux
         !  ------   
            F(IRHO , IX) = e % storage % Q(IRHO,i,j,k)  * u_g
            F(IRHOU, IX) = e % storage % Q(IRHOU,i,j,k) * u_g 
            F(IRHOV, IX) = e % storage % Q(IRHOV,i,j,k) * u_g
            F(IRHOW, IX) = e % storage % Q(IRHOW,i,j,k) * u_g
            F(IRHOE, IX) = e % storage % Q(IRHOE,i,j,k) * u_g
         
         !  Y-Flux
         !  ------
            F(IRHO , IY) = e % storage % Q(IRHO,i,j,k)  * v_g
            F(IRHOU, IY) = e % storage % Q(IRHOU,i,j,k) * v_g
            F(IRHOV, IY) = e % storage % Q(IRHOV,i,j,k) * v_g
            F(IRHOW, IY) = e % storage % Q(IRHOW,i,j,k) * v_g
            F(IRHOE, IY) = e % storage % Q(IRHOE,i,j,k) * v_g
         
         !  Z-Flux
         !  ------
            F(IRHO , IZ) = e % storage % Q(IRHO,i,j,k)  * w_g
            F(IRHOU, IZ) = e % storage % Q(IRHOU,i,j,k) * w_g
            F(IRHOV, IZ) = e % storage % Q(IRHOV,i,j,k) * w_g
            F(IRHOW, IZ) = e % storage % Q(IRHOW,i,j,k) * w_g
            F(IRHOE, IZ) = e % storage % Q(IRHOE,i,j,k) * w_g
   !        
         !  Contravariant flux computation 
   
            contravariantFlux(:,i,j,k,IX) =    F(:,IX) * e % geom % jGradXi(IX,i,j,k)   &
                                             + F(:,IY) * e % geom % jGradXi(IY,i,j,k)   &
                                             + F(:,IZ) * e % geom % jGradXi(IZ,i,j,k)
   
            contravariantFlux(:,i,j,k,IY) =    F(:,IX) * e % geom % jGradEta(IX,i,j,k)  &
                                             + F(:,IY) * e % geom % jGradEta(IY,i,j,k)  &
                                             + F(:,IZ) * e % geom % jGradEta(IZ,i,j,k)
   
            contravariantFlux(:,i,j,k,IZ) =    F(:,IX) * e % geom % jGradZeta(IX,i,j,k) &
                                             + F(:,IY) * e % geom % jGradZeta(IY,i,j,k) &
                                             + F(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)
   
      enddo;enddo;enddo
      end subroutine SlidingMeshFluxCalculation

      


end module SpatialDiscretization
