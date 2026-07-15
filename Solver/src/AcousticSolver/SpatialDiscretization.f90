#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      use HyperbolicDiscretizations
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
      use VariableConversion, only: NSGradientVariables_STATE, GetGradientValues_f
      use ProblemFileFunctions, only: UserDefinedSourceTermNS_f
      use BoundaryConditions
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public   ComputeTimeDerivative, ComputeTimeDerivativeIsolated, viscousDiscretizationKey
      public   Initialize_SpaceAndTimeMethods, Finalize_SpaceAndTimeMethods

      procedure(GetGradientValues_f), pointer :: getGradients_BaseFlow

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
         implicit none
         class(FTValueDictionary),  intent(in)  :: controlVariables
         class(DGSem)                           :: sem
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)       :: inviscidDiscretizationName
         ! character(len=LINE_LENGTH)       :: viscousDiscretizationName
         character(len=*), parameter      :: gradient_variables_key = "gradient variables"
         character(len=LINE_LENGTH)       :: gradient_variables
         real(RP)                         :: hnmin, hnmax

         ! call hnRange(sem % mesh, hnmin, hnmax)

         if ( MPI_Process % isRoot ) then
            write(STD_OUT,'(/)')
            call Section_Header("Spatial discretization scheme")
            ! write(STD_OUT,'(/)')
            ! write(STD_OUT,'(30X,A,A30,1pG10.3)') "->", "Minimum h/N: ", hnmin
            ! write(STD_OUT,'(30X,A,A30,1pG10.3)') "->", "Maximum h/N: ", hnmax
            write(STD_OUT,'(/)')
         end if
   !
   !     Initialize inviscid discretization
   !     ----------------------------------
         inviscidDiscretizationName = controlVariables % stringValueForKey(inviscidDiscretizationKey,requestedLength = LINE_LENGTH)

         call toLower(inviscidDiscretizationName)

         select case ( trim(inviscidDiscretizationName) )
         case ( "standard" )
            if (.not. allocated(HyperbolicDiscretization)) allocate( StandardDG_t  :: HyperbolicDiscretization )

         ! case ( "split-form")
            ! if (.not. allocated(HyperbolicDiscretization)) allocate( SplitDG_t     :: HyperbolicDiscretization)

         case default
            write(STD_OUT,'(A,A,A)') 'Requested inviscid discretization "',trim(inviscidDiscretizationName),'" is not implemented.'
            write(STD_OUT,'(A)') "Implemented discretizations are:"
            write(STD_OUT,'(A)') "  * Standard"
            ! write(STD_OUT,'(A)') "  * Split-Form"
            errorMessage(STD_OUT)
            error stop

         end select

         call HyperbolicDiscretization % Initialize(controlVariables)
   !
   !     Initialize gradients of base flow
   !     ---------------------------------
!
!        Set state as default option
!        ---------------------------
         ! call SetGradientVariables(GRADVARS_STATE)
         ! getGradients_BaseFlow => NSGradientVariables_STATE
         ! call set_getVelocityGradients(GRADVARS_STATE)

      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_SpaceAndTimeMethods
         implicit none

         if ( allocated(HyperbolicDiscretization) ) deallocate( HyperbolicDiscretization )

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
         ! use ActuatorLine, only: farm
         use SpongeClass, only: sponge, addSourceSponge
         use APESourceClass, only: apeSource, addAPESource
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
         integer     :: eID , i, j, k, ierr, fID, iFace, iEl, iP, STLNum, n, m
         integer,       allocatable :: i_(:), j_(:), k_(:)
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
                     call computeElementInterfaceFlux(masterFace1=mesh % faces(fID), f=mesh % faces(mesh % faces(fID)%Mortar(m)), m=m)
                  end if
               end do
            elseif (mesh % faces(fID) % MortarType == MORTAR_NONE) then
               call computeElementInterfaceFlux(f=mesh % faces(fID))
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
!        Wait until messages are sent
!        ****************************
!
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
!$omp single
           call mesh % GatherMPIFacesSolution(NCONS)
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
                        call computeElementInterfaceFlux(masterFace1=mesh % faces(fID), f=mesh % faces(mesh % faces(fID)%Mortar(m)))
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
!
!        *****************************************************************************************************************************
!        Compute contributions to source term
!        *****************************************************************************************************************************
!           Add physical source term
!           ************************
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               ! the source term is reset to 0 each time Qdot is calculated to enable the possibility to add source terms to
               ! different contributions and not accumulate each call
               e % storage % S_NS = 0.0_RP
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues, &
                     e % storage % Qbase(:,i,j,k), e % storage % Lambbase(:,i,j,k), e % storage % Lamb(:,i,j,k))
                  ! call farm % ForcesFarm(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), e % storage % S_NS(:,i,j,k), t)
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do
            ! for the sponge, loops are in the internal subroutine as values are precalculated
            call addSourceSponge(sponge,mesh)

            ! add source term for the acoustics
            call addAPESource(apeSource, mesh)
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

      end subroutine TimeDerivative_ComputeQDot
!
!////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------------
!     This routine computes Qdot neglecting the interaction with neighboring elements
!     and boundaries. Therefore, the external states are not needed.
!     -------------------------------------------------------------------------------
      subroutine TimeDerivative_ComputeQDotIsolated( mesh , t )
         ! use ActuatorLine, only: farm
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
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               e % storage % S_NS = 0.0_RP
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues, &
                        e % storage % Qbase(:,i,j,k), e % storage % Lambbase(:,i,j,k), e % storage % Lamb(:,i,j,k))
                  ! call farm % ForcesFarm(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), e % storage % S_NS(:,i,j,k), t)
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do

!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S_NS(:,i,j,k)
            end do                  ; end do                ; end do
            end associate
         end do
!$omp end do

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
         real(kind=RP) :: viscousContravariantFlux ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP) :: contravariantFlux        ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )

!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , APEFlux, inviscidContravariantFlux, useBaseFlow_ = .true. )
!
!        Compute viscous contravariant flux
!        ----------------------------------
         viscousContravariantFlux = 0.0_RP
!
!        ************************
!        Perform volume integrals
!        ************************
!        Compute the total Navier-Stokes flux
!        ------------------------------------
         contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux
!
!        Perform the Weak Volume Green integral
!        --------------------------------------
         e % storage % QDot = ScalarStrongIntegrals % StdVolumeGreen ( e , NCONS, contravariantFlux )
!
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
         real(kind=RP) :: viscousContravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP) :: contravariantFlux         ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )

!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , APEFlux, inviscidContravariantFlux, useBaseFlow_ = .true. )
!
!        Compute viscous contravariant flux
!        ----------------------------------
         viscousContravariantFlux = 0.0_RP
!
!        ************************
!        Perform volume integrals
!        ************************
!
!        Compute the total Navier-Stokes flux
!        ------------------------------------
         contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux
! 
!        Perform the Weak Volume Green integral
!        --------------------------------------
         e % storage % QDot = ScalarWeakIntegrals % StdVolumeGreen ( e, NCONS, contravariantFlux )

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
      subroutine computeElementInterfaceFlux(f, masterFace1, m)
         use FaceClass
         use RiemannSolvers_CAA
         implicit none
         type(Face)   , intent(inout) :: f
         type(Face), optional, intent(inout) :: masterFace1
         integer, optional, intent(in) :: m
         integer       :: i, j
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS)
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         integer       :: Sidearray(2)
!
!        --------------
!        Viscous fluxes
!        --------------
!
         visc_flux = 0.0_RP

         do j = 0, f % Nf(2)
            do i = 0, f % Nf(1)
!
!              --------------
!              Invscid fluxes
!              --------------
!
               call RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  QbaseL = f % storage(1) % Qbase(:,i,j), &
                                  QbaseR = f % storage(2) % Qbase(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux) * f % geom % jacobian(i,j)

            end do
         end do
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         if (f % MortarType == MORTAR_NONE) then
            Sidearray = (/1,2/)
            call f % ProjectFluxToElements(NCONS, flux, Sidearray)
         end if
         if (f % MortarType == MORTAR_SMALL4 .and. present(masterFace1)) then
            Sidearray = (/1,0/)
            call masterFace1 % ProjectMortarFluxToElements(nEqn=NCONS, whichElements=Sidearray, &
               slaveFace=f, MortarFlux=flux)
            Sidearray = (/0,2/)
            call f % ProjectFluxToElements(NCONS, flux, Sidearray)
         end if

      end subroutine computeElementInterfaceFlux

      subroutine computeMPIFaceFlux(f)
         use FaceClass
         use RiemannSolvers_CAA
         implicit none
         type(Face)   , intent(inout) :: f
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS)
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         integer       :: Sidearray(2)
!
!        ---------------------------
!        Artificial viscosity fluxes
!        ---------------------------
!
!        --------------
!        Viscous fluxes
!        --------------
!
         visc_flux = 0.0_RP

         do j = 0, f % Nf(2)
            do i = 0, f % Nf(1)
!
!              --------------
!              Invscid fluxes
!              --------------
!
               call RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  QbaseL = f % storage(1) % Qbase(:,i,j), &
                                  QbaseR = f % storage(2) % Qbase(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux) * f % geom % jacobian(i,j) 

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
         if (f % MortarType == MORTAR_SMALL4) then
            call f % Interpolatesmall2big(NCONS, flux)
         end if

      end subroutine ComputeMPIFaceFlux

      SUBROUTINE computeBoundaryFlux(f, time, mesh)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_CAA
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
INTEGER                         :: i, j
!
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: inv_flux(NCONS)
      real(kind=RP)                   :: visc_flux(NCONS)
      real(kind=RP)                   :: fStar(NCONS, 0:f % Nf(1), 0: f % Nf(2))
      real(kind=RP)                   :: fv_3d(NCONS,NDIM)
      integer                         :: Sidearray(2)

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

      visc_flux = 0.0_RP

      do j = 0, f % Nf(2)
         do i = 0, f % Nf(1)
!
!           Hyperbolic part
!           -------------
            call RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                               QRight = f % storage(2) % Q(:,i,j), &
                               QbaseL = f % storage(1) % Qbase(:,i,j), &
                               QbaseR = f % storage(1) % Qbase(:,i,j), &
                               nHat   = f % geom % normal(:,i,j), &
                               t1     = f % geom % t1(:,i,j), &
                               t2     = f % geom % t2(:,i,j), &
                               flux   = inv_flux)

            fStar(:,i,j) = (inv_flux - visc_flux ) * f % geom % jacobian(i,j)
         end do
      end do

      Sidearray = (/1, HMESH_NONE/)
      call f % ProjectFluxToElements(NCONS, fStar, Sidearray)

      end subroutine computeBoundaryFlux

end module SpatialDiscretization