#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      use HyperbolicDiscretizations
      use EllipticDiscretizations
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
      use VariableConversion, only: iNSGradientVariables, GetiNSOneFluidViscosity, GetiNSTwoFluidsViscosity
      use ProblemFileFunctions
      use BoundaryConditions, only: BCs
      use ProblemFileFunctions, only: UserDefinedSourceTermNS_f
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public   ComputeTimeDerivative, ComputeTimeDerivativeIsolated, viscousDiscretizationKey
      public   Initialize_SpaceAndTimeMethods, Finalize_SpaceAndTimeMethods

      abstract interface
      SUBROUTINE computeElementInterfaceFluxF(f, masterFace)
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

      SUBROUTINE computeBoundaryFluxF(f, time)
         use SMConstants
         use FaceClass,  only: Face
         IMPLICIT NONE
         type(Face),    intent(inout) :: f
         REAL(KIND=RP)                :: time
      end subroutine computeBoundaryFluxF
      end interface
      
      character(len=LINE_LENGTH), parameter  :: viscousDiscretizationKey = "viscous discretization"
      procedure(GetViscosity_f), pointer, protected :: GetViscosity 
!
!     ========      
      CONTAINS 
!     ========      
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_SpaceAndTimeMethods(controlVariables, mesh)
         use FTValueDictionaryClass
         use Utilities, only: toLower
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         implicit none
         class(FTValueDictionary),  intent(in)  :: controlVariables
         class(HexMesh)                         :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)       :: inviscidDiscretizationName
         character(len=LINE_LENGTH)       :: viscousDiscretizationName
         
         if (.not. mesh % child) then ! If this is a child mesh, all these constructs were already initialized for the parent mesh
         
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,'(/)')
               call Section_Header("Spatial discretization scheme")
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

               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_iNS)
               call ViscousDiscretization % Describe

               select case (thermodynamics % number_of_fluids)
               case(1)
                  GetViscosity => GetiNSOneFluidViscosity
               case(2)
                  GetViscosity => GetiNSTwoFluidsViscosity
               end select

!
!        Compute wall distances
!        ----------------------
         call mesh % ComputeWallDistances
         
         end if

      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_SpaceAndTimeMethods
         implicit none
         IF ( ALLOCATED(HyperbolicDiscretization) ) DEALLOCATE( HyperbolicDiscretization )
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
         integer,             intent(in) :: mode
         logical, intent(in), optional   :: HO_Elements
         logical, intent(in), optional   :: element_mask(:)
		 integer, intent(in), optional   :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k, eID
!
!        *******************************************************************
!        Construct the auxiliary state for the fluxes with density positivity
!        *******************************************************************
!
!$omp do schedule(runtime)
         do eID = 1, size(mesh % elements)
            mesh % elements(eID) % storage % rho = mesh % elements(eID) % storage % Q(INSRHO,:,:,:)
            mesh % elements(eID) % storage % Q(INSRHO,:,:,:) = min(max(mesh % elements(eID) % storage % Q(INSRHO,:,:,:), thermodynamics % rho_min), &
                                                                   thermodynamics % rho_max)
         end do
!$omp end do nowait
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time)
         call mesh % ProlongSolutionToFaces(NCONS)
!
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
            CALL DGSpatial_ComputeGradient(mesh , time)
         end if

#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesGradients(NCONS)
!$omp end single
#endif
!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call ComputeNSTimeDerivative(mesh = mesh , &
                                         particles = particles, &
                                         t    = time)
!
!        ***************************************
!        Return the density to its default value
!        ***************************************
!
!$omp do schedule(runtime)
         do eID = 1, size(mesh % elements)
             mesh % elements(eID) % storage % Q(INSRHO,:,:,:) = mesh % elements(eID) % storage % rho 
         end do
!$omp end do

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
            CALL BaseClass_ComputeGradient( ViscousDiscretization, NCONS, NCONS, mesh , time, iNSGradientVariables)
!
!           The prolongation is usually done in the viscous methods, but not in the BaseClass
!           ---------------------------------------------------------------------------------
            call mesh % ProlongGradientsToFaces(NCONS)
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
!
!////////////////////////////////////////////////////////////////////////////////////
!
!           Navier--Stokes procedures
!           -------------------------
!
!////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ComputeNSTimeDerivative( mesh , particles, t )
         use SpongeClass, only: sponge, addSourceSponge
         use ActuatorLine, only: farm, ForcesFarm
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
         integer     :: eID , i, j, k, ierr, fID, m
         real(kind=RP)  :: mu_smag, delta
         real(kind=RP), dimension(NCONS)  :: Source
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) 
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_VolumetricContribution( mesh % elements(eID) , t)
         end do
!$omp end do nowait
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) 
         do fID = 1, size(mesh % faces) 
            associate( f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_INTERIOR) 
               if (f % IsMortar==1) then 
               associate(fstar=>mesh% faces(fID)%storage(1)%fStar)
                  fstar=0.0_RP
               end associate
               associate(fstar=>mesh% faces(fID)%storage(2)%fStar)
                  fstar=0.0_RP
               end associate
               do m=1,4
                  if (f%Mortar(m) .ne. 0) then 
            CALL computeElementInterfaceFlux_iNS(masterFace=f, f=mesh % faces(mesh % faces(fID)%Mortar(m))) 
                  end if 
               end do 
            elseif  (f % IsMortar==0) then
                  CALL computeElementInterfaceFlux_iNS( f ) 
               end if 
            case (HMESH_BOUNDARY) 
               CALL computeBoundaryFlux_iNS(f, t) 
 
            end select 
            end associate 
         end do 
!$omp end do 
!
!        **************************************************************
!        Surface integrals and scaling of elements without shared faces
!        **************************************************************
! 
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            if ( e % hasSharedFaces ) cycle
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
            call mesh % GatherMPIFacesGradients(NCONS)
!$omp end single
!
!           **************************************
!           Compute Riemann solver of shared faces
!           **************************************
!
!$omp do schedule(runtime) 
            do fID = 1, size(mesh % faces) 
               associate( f => mesh % faces(fID)) 
               select case (f % faceType) 
               case (HMESH_MPI)
                  if (mesh% faces(fID)%IsMortar==1) then 
                     !write(*,*) 'big mortar face mpi'
                     associate(fstar=>mesh% faces(fID)%storage(1)%fStar)
                        fstar=0.0_RP
                     end associate
                     do m=1,4
                        if (f%Mortar(m) .ne. 0) then 
                           CALL computeElementInterfaceFlux_iNS(masterFace=f, f=mesh % faces(mesh % faces(fID)%Mortar(m)))
                        end if 
                     end do 
                  end if  
                  CALL computeMPIFaceFlux_iNS( f ) 
               end select 
               end associate 
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
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, size(mesh % elements) 
               associate(e => mesh % elements(eID)) 
               if ( .not. e % hasSharedFaces ) cycle
               call TimeDerivative_FacesContribution(e, t, mesh) 
   
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k) 
               end do         ; end do          ; end do 
               end associate 
            end do
!$omp end do
!
!           Add a MPI Barrier
!           -----------------
!$omp single
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
         end if
#endif
!
!        ***********
!        Add gravity
!        ***********
!
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, size(mesh % elements)
               associate(e => mesh % elements(eID))
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  e % storage % QDot(INSRHOU:INSRHOW,i,j,k) = e % storage % QDot(INSRHOU:INSRHOW,i,j,k) + &
                                                        e % storage % Q(INSRHO,i,j,k) * &
                                    dimensionless % invFr2 * dimensionless % gravity_dir

               end do                ; end do                ; end do
               end associate
            end do
!$omp end do
!
!           ***************
!           Add source term
!           ***************
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               e % storage % S_NS = 0.0_RP
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do

! for the sponge, loops are in the internal subroutine as values are precalculated
            call addSourceSponge(sponge,mesh)
            call ForcesFarm(farm, mesh, t)

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
!           ********************
!           Add Particles source
!           ********************
            !if (.not. mesh % child) then
               !if ( particles % active ) then             
!!$omp do schedule(runtime)
                  !do eID = 1, size(mesh % elements)
                  !   call particles % AddSource(mesh % elements(eID), t, thermodynamics, dimensionless, refValues)
                  !end do
!!$omp end do
               !endif 
            !end if
!
!        *********************
!        Add IBM source term
!        *********************
! no wall function for INCNS
         if( mesh% IBM% active ) then
            if( .not. mesh% IBM% semiImplicit ) then 
!$omp do schedule(runtime) private(i,j,k,Source)
                  do eID = 1, mesh % no_of_elements  
                     associate ( e => mesh % elements(eID) ) 
                     do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                        if( e% isInsideBody(i,j,k) ) then
                           ! only without moving for now in INCNS
                           if( .not. mesh% IBM% stl(e% STL(i,j,k))% move ) then 
                              call mesh% IBM% SourceTerm( eID = eID, Q = e % storage % Q(:,i,j,k), Source = Source, wallfunction = .false. )
                           end if 
                           e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + Source
                        end if
                     end do                  ; end do                ; end do
                     end associate
                  end do
!$omp end do       
            end if 
         end if

      end subroutine ComputeNSTimeDerivative
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
         integer     :: eID , i, j, k, fID
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) 
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_StrongVolumetricContribution( mesh % elements(eID) , t)
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
         
      end subroutine TimeDerivative_ComputeQDotIsolated
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_StrongVolumetricContribution( e , t )
         use HexMeshClass
         use ElementClass
         implicit none
         type(Element)      :: e
         real(kind=RP)      :: t

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
         real(kind=RP) :: contravariantFlux         ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , iEulerFlux, inviscidContravariantFlux ) 
!
!        Compute viscous contravariant flux
!        ----------------------------------
         call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NCONS, iViscousFlux, GetViscosity, e , viscousContravariantFlux) 
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
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux 
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
!~             viscousContravariantFlux = viscousContravariantFlux + SVVContravariantFlux

!~             e % storage % QDot = -ScalarWeakIntegrals % SplitVolumeDivergence( e, fSharp, gSharp, hSharp, viscousContravariantFlux)

         end select

      end subroutine TimeDerivative_StrongVolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_VolumetricContribution( e , t )
         use HexMeshClass
         use ElementClass
         implicit none
         type(Element)      :: e
         real(kind=RP)      :: t

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
         real(kind=RP) :: contravariantFlux         ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , iEulerFlux, inviscidContravariantFlux ) 
!
!        Compute viscous contravariant flux
!        ----------------------------------
         call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NCONS, iViscousFlux, GetViscosity, e , viscousContravariantFlux) 
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
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux 
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
            viscousContravariantFlux = viscousContravariantFlux 

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
      SUBROUTINE computeElementInterfaceFlux_iNS(f, masterFace)
         use FaceClass
         use RiemannSolvers_iNS
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f  
         type(Face), optional, intent(inout) :: masterFace 

         integer       :: i, j
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: muL, muR, mu

         !if (f % IsMortar==0 .OR. f % IsMortar==2) then 
            DO j = 0, f % Nf(2)
               DO i = 0, f % Nf(1)

                  call GetViscosity(f % storage(1) % Q(INSRHO,i,j), muL)
                  call GetViscosity(f % storage(2) % Q(INSRHO,i,j), muR)
                  mu = 0.5_RP * (muL + muR)
   !      
   !              --------------
   !              Viscous fluxes
   !              --------------
   !      
                  CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NCONS, &
                                                   EllipticFlux = iViscousFlux, &
                                                   f = f, &
                                                   QLeft = f % storage(1) % Q(:,i,j), &
                                                   QRight = f % storage(2) % Q(:,i,j), &
                                                   U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                   U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                   U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                   U_xRight = f % storage(2) % U_x(:,i,j), &
                                                   U_yRight = f % storage(2) % U_y(:,i,j), &
                                                   U_zRight = f % storage(2) % U_z(:,i,j), &
                                                   mu_left  = [mu, 0.0_RP, 0.0_RP], &
                                                   mu_right = [mu, 0.0_RP, 0.0_RP], &
                                                   nHat = f % geom % normal(:,i,j) , &
                                                   dWall = f % geom % dWall(i,j), &
                                                   flux  = visc_flux(:,i,j) )

               end do
            end do

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
                  flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
                  
               END DO   
            END DO  
   !
   !        ---------------------------
   !        Return the flux to elements
   !        ---------------------------
   !
      if (f % IsMortar==0) then 
            call f % ProjectFluxToElements(NCONS, flux, (/1,2/))
       end if 
       if (f % IsMortar==2 .and. present(masterFace)) then 
         call masterFace % ProjectMortarFluxToElements(nEqn=NCONS, whichElements=(/1,0/), &
         slaveFace=f, MortarFlux=flux)
         call f % ProjectFluxToElements(NCONS, flux, (/0,2/))
      end if 
 !end if 

      END SUBROUTINE computeElementInterfaceFlux_iNS

      SUBROUTINE computeMPIFaceFlux_iNS(f)
         use FaceClass
         use RiemannSolvers_iNS
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu
!
!        --------------
!        Invscid fluxes
!        --------------
!
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               call GetViscosity(f % storage(1) % Q(INSRHO,i,j), mu)

               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NCONS, &
                                                  EllipticFlux = iViscousFlux, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu_left  = [mu, 0.0_RP, 0.0_RP], &
                                                  mu_right = [mu, 0.0_RP, 0.0_RP], &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

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
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectFluxToElements(NCONS, flux, (/thisSide, HMESH_NONE/))
         if (f % IsMortar==2) then 
            !write(*,*) 'this side', thisSide
            call f% Interpolatesmall2big(NCONS, flux)
           
         end if 
      end subroutine ComputeMPIFaceFlux_iNS

      SUBROUTINE computeBoundaryFlux_iNS(f, time)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_iNS
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: inv_flux(NCONS), fv_3d(NCONS,NDIM)
      real(kind=RP)                   :: visc_flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: fStar(NCONS, 0:f % Nf(1), 0: f % Nf(2))
      real(kind=RP)                   :: mu
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL BCs(f % zone) % bc % FlowState( &
                                      f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j))

      end do               ; end do

      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         call GetViscosity(f % storage(1) % Q(INSRHO,i,j), mu)

         call iViscousFlux(NCONS,NGRAD,f % storage(1) % Q(:,i,j), &
                                       f % storage(1) % U_x(:,i,j), &
                                       f % storage(1) % U_y(:,i,j), &
                                       f % storage(1) % U_z(:,i,j), &
                                       mu, 0.0_RP, 0.0_RP, fv_3d)

         visc_flux(:,i,j) =   fv_3d(:,IX)*f % geom % normal(IX,i,j) &
                            + fv_3d(:,IY)*f % geom % normal(IY,i,j) &
                            + fv_3d(:,IZ)*f % geom % normal(IZ,i,j) 

         CALL BCs(f % zone) % bc % FlowNeumann(&
                  f % geom % x(:,i,j),         & 
                  time,                        & 
                  f % geom % normal(:,i,j),    & 
                  f % storage(1) % Q(:,i,j),   & 
                  f % storage(1) % U_x(:,i,j), & 
                  f % storage(1) % U_y(:,i,j), & 
                  f % storage(1) % U_z(:,i,j), & 
                  visc_flux(:,i,j)              )
!
!        Hyperbolic part
!        -------------
         CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                            QRight = f % storage(2) % Q(:,i,j), &   
                            nHat   = f % geom % normal(:,i,j), &
                            t1     = f % geom % t1(:,i,j), &
                            t2     = f % geom % t2(:,i,j), &
                            flux   = inv_flux)

         fStar(:,i,j) = (inv_flux - visc_flux(:,i,j) ) * f % geom % jacobian(i,j)

      end do               ; end do

      call f % ProjectFluxToElements(NCONS, fStar, (/1, HMESH_NONE/))

      END SUBROUTINE computeBoundaryFlux_iNS
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              GRADIENT PROCEDURES
!              -------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_ComputeGradient( mesh , time)
         use HexMeshClass
         use PhysicsStorage, only: NCONS
         implicit none
         type(HexMesh)                  :: mesh
         real(kind=RP),      intent(in) :: time

         call ViscousDiscretization % ComputeGradient( NCONS, NCONS, mesh , time , iNSGradientVariables)

      end subroutine DGSpatial_ComputeGradient

      subroutine DensityLimiter(N,Q)
         implicit none
         integer,       intent(in)    :: N(3)
         real(kind=RP), intent(inout) :: Q(1:NCONS,0:N(1),0:N(2),0:N(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k 
         real(kind=RP) :: rhoIn01, p, rhomin, rhomax

         rhomin = thermodynamics % rho_min / refValues % rho
         rhomax = thermodynamics % rho_max / refValues % rho

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            rhoIn01 = (Q(INSRHO,i,j,k)-rhomin)/(rhomax-rhomin)

            if ( rhoIn01 .ge. 1.0_RP ) then
               Q(INSRHO,i,j,k) = rhomax

            elseif ( rhoIn01 .le. 0.0_RP ) then
               Q(INSRHO,i,j,k) = rhomin

            else
               !p = POW3(rhoIn01)*(6.0_RP*POW2(rhoIn01)-15.0_RP*rhoIn01+10.0_RP)
               p = rhoIn01

               Q(INSRHO,i,j,k) = (rhomax-rhomin)*p + rhomin
         
            end if

         end do         ; end do         ; end do

      end subroutine DensityLimiter
!
!////////////////////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization
