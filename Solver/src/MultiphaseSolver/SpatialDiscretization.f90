#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      use HyperbolicDiscretizations
      use EllipticDiscretizations
      use DGIntegrals
      use MeshTypes
      use LESModels
      use HexMeshClass
      use ElementClass
      use PhysicsStorage
      use Physics
      use MPI_Face_Class
      use MPI_Process_Info
      use DGSEMClass
      use FluidData
      use VariableConversion, only: mGradientVariables, GetmOneFluidViscosity,&
                                    GetmTwoFluidsViscosity, chGradientVariables,&
                                    GetCHViscosity
      use BoundaryConditions, only: BCs, SetBoundaryConditionsEqn, NS_BC, C_BC, MU_BC
      use ProblemFileFunctions, only: UserDefinedSourceTermNS_f
      use ParticlesClass
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
      character(len=LINE_LENGTH), parameter     :: CHDiscretizationKey = "cahn-hilliard discretization"
      character(len=LINE_LENGTH), parameter  :: FLUID1_COMPRESSIBILITY_KEY = "fluid 1 sound speed square (m/s)"


      real(kind=RP), protected :: IMEX_S0 = 0.0_RP 
      real(kind=RP), protected :: IMEX_K0 = 1.0_RP
      logical :: use_non_constant_speed_of_sound = .false.
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
         character(len=LINE_LENGTH) :: inviscidDiscretizationName
         character(len=LINE_LENGTH) :: viscousDiscretizationName
         character(len=LINE_LENGTH) :: CHDiscretizationName
         
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
               print*, "There are no split-forms available for the Multiphase Solver"
               errorMessage(STD_OUT)
               stop
            case default
               write(STD_OUT,'(A,A,A)') 'Requested inviscid discretization "',trim(inviscidDiscretizationName),'" is not implemented.'
               write(STD_OUT,'(A)') "Implemented discretizations are:"
               write(STD_OUT,'(A)') "  * Standard"
               errorMessage(STD_OUT)
               stop 

            end select
               
            call HyperbolicDiscretization % Initialize(controlVariables)
   !
   !        Initialize viscous discretization
   !        ---------------------------------         
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

               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_MU)
               call ViscousDiscretization % Describe
!
!           Compute wall distances
!           ----------------------
            call mesh % ComputeWallDistances

!           Initialize models
!           -----------------
            call InitializeLESModel(LESModel, controlVariables)
!
!           Initialize Cahn--Hilliard discretization
!           ----------------------------------------         
            if ( .not. controlVariables % ContainsKey(CHDiscretizationKey) ) then
               print*, "Input file is missing entry for keyword: Cahn-Hilliard discretization"
               errorMessage(STD_OUT)
               stop
            end if
   
            CHDiscretizationName = controlVariables % stringValueForKey(CHDiscretizationKey, requestedLength = LINE_LENGTH)
            call toLower(CHDiscretizationName)
            
            select case ( trim(CHDiscretizationName) )
            case("br1")
               allocate(BassiRebay1_t     :: CHDiscretization)
   
            case("br2")
               allocate(BassiRebay2_t     :: CHDiscretization)
   
            case("ip")
               allocate(InteriorPenalty_t :: CHDiscretization)
   
            case default
               write(STD_OUT,'(A,A,A)') 'Requested viscous discretization "',trim(CHDiscretizationName),'" is not implemented.'
               write(STD_OUT,'(A)') "Implemented discretizations are:"
               write(STD_OUT,'(A)') "  * BR1"
               write(STD_OUT,'(A)') "  * BR2"
               write(STD_OUT,'(A)') "  * IP"
               errorMessage(STD_OUT)
               stop 
   
            end select
   
            use_non_constant_speed_of_sound = controlVariables % ContainsKey(FLUID1_COMPRESSIBILITY_KEY)

            call CHDiscretization % Construct(controlVariables, ELLIPTIC_CH)
            call CHDiscretization % Describe
			
			if ( .not. MPI_Process % isRoot ) return
			
            if(use_non_constant_speed_of_sound) then
               write(STD_OUT,'(A)') "  Implementing artificial compressibility with a non-constant speed of sound in each phase"
            else
               write(STD_OUT,'(A)') "  Implementing artificial compressibility with a constant ACM factor in each phase"
            endif
         
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
         integer, intent(in)             :: mode
         logical, intent(in), optional   :: HO_Elements
         logical, intent(in), optional   :: element_mask(:)
		 integer, intent(in), optional   :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                 :: k, eID, fID, i, j, ierr, locLevel, lID
         real(kind=RP)           :: sqrtRho, invMa2
         logical                 :: compute_element
		 logical, allocatable    :: face_mask(:)
         real(kind=RP)           :: mu_smag, delta

		 if (present(Level)) then
            locLevel = Level
         else
            locLevel = 1
         end if
		 
		 associate ( MLIter_eID => mesh % MLRK % MLIter_eID, MLIter => mesh % MLRK % MLIter , MLIter_fID => mesh % MLRK % MLIter_fID)																										   


        if (present(element_mask)) then
            call CreateFaceMask(mesh, element_mask, face_mask)
        endif

!$omp parallel shared(mesh, time, locLevel) private(compute_element)
!
!///////////////////////////////////////////////////
!        1st step: Get chemical potential
!///////////////////////////////////////////////////
!
!        ------------------------------------------
!        Update concentration with the state vector
!        ------------------------------------------
!
         select case (mode)
         case (CTD_IGNORE_MODE,CTD_IMEX_EXPLICIT)
!$omp do schedule(runtime) private(eID)
            do lID = 1, MLIter(locLevel,1)
			   eID = MLIter_eID(lID)
               compute_element = .true.
               if (present(element_mask)) compute_element = element_mask(eID)
               
               if (compute_element) then
                   mesh % elements(eID) % storage % c(1,:,:,:) = mesh % elements(eID) % storage % QNS(IMC,:,:,:)
               endif 
            end do
!$omp end do         
         end select
!
!        -------------------------------
!        Set memory to Cahn-Hilliard (C)
!        -------------------------------
!
!$omp single
         call mesh % SetStorageToEqn(C_BC)
         select case (mode)
         case (CTD_IGNORE_MODE,CTD_IMEX_EXPLICIT)
            call SetBoundaryConditionsEqn(C_BC)

         case (CTD_IMEX_IMPLICIT,CTD_LAPLACIAN)
            call SetBoundaryConditionsEqn(MU_BC)

         end select
!$omp end single
!
!        --------------------------------------------
!        Prolong Cahn-Hilliard concentration to faces
!        --------------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP, element_mask=element_mask, Level=locLevel)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution(NCOMP)
!$omp end single
#endif
!
!        ------------------------------------------------------------
!        Get concentration (lifted) gradients (also prolong to faces)
!        ------------------------------------------------------------
!
         call CHDiscretization % ComputeGradient(NCOMP, NCOMP, mesh, time, chGradientVariables, element_mask=element_mask, Level = locLevel)
!
!        --------------------
!        Update MPI Gradients
!        --------------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesGradients(NCOMP)
!$omp end single
#endif
!
!        ----------------------
!        Get chemical potential
!        ----------------------
!
!        Get the concentration Laplacian (into QDot => cDot)

         call ComputeLaplacian(mesh, time, element_mask, locLevel)

         select case (mode)
         case (CTD_IGNORE_MODE, CTD_IMEX_EXPLICIT)
!$omp do schedule(runtime) private(eID)
            do lID = 1, MLIter(locLevel,1)
			   eID = MLIter_eID(lID)
               compute_element = .true.
               if (present(element_mask)) compute_element = element_mask(eID)
               
               if (compute_element) then
!
!              + Linear part
                 !mesh % elements(eID) % storage % mu = - POW2(multiphase % eps)* mesh % elements(eID) % storage % QDot
                 mesh % elements(eID) % storage % mu = - 1.5_RP * multiphase % eps * multiphase % sigma * mesh % elements(eID) % storage % QDot
!  
!              + NonLinear part
                 !call AddQuarticDWPDerivative(mesh % elements(eID) % storage % c, mesh % elements(eID) % storage % mu)
                 call Multiphase_AddChemFEDerivative(mesh % elements(eID) % storage % c, mesh % elements(eID) % storage % mu)
               endif
            end do
!$omp end do         
         case (CTD_IMEX_IMPLICIT)
!$omp do schedule(runtime) private(eID)
            do lID = 1, MLIter(locLevel,1)
			   eID = MLIter_eID(lID)
               compute_element = .true.
               if (present(element_mask)) compute_element = element_mask(eID)
               
               if (compute_element) then
!  
!              + Linear part
                 !mesh % elements(eID) % storage % mu = - IMEX_K0 * POW2(multiphase % eps) * mesh % elements(eID) % storage % QDot &
                 mesh % elements(eID) % storage % mu = - 1.5_RP * IMEX_K0 * multiphase % eps * multiphase % sigma * mesh % elements(eID) % storage % QDot &
                                                       + IMEX_S0 * mesh % elements(eID) % storage % c
!              + Multiply by mobility
                 mesh % elements(eID) % storage % mu = multiphase % M0 * mesh % elements(eID) % storage % mu
               endif
            end do
!$omp end do         
         end select
!
!        -----------------------------------
!        Prolong chemical potential to faces
!        -----------------------------------
!
         select case(mode)
         case(CTD_LAPLACIAN)
         case default
!$omp single
         call mesh % SetStorageToEqn(MU_BC)
!$omp end single
         call mesh % ProlongSolutionToFaces(NCOMP, element_mask=element_mask, Level=locLevel)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution(NCOMP)
         call mesh % GatherMPIFacesSolution(NCOMP)
!$omp end single
#endif
         end select
!
!/////////////////////////////////////////////////////////////////////////////////
!        2nd step: If IMEX_IMPLCIIT, get the chemical potential laplacian and exit
!/////////////////////////////////////////////////////////////////////////////////
!
         select case (mode)
         case (CTD_IMEX_IMPLICIT)
!
!           ------------------------------------------------------------
!           Get concentration (lifted) gradients (also prolong to faces)
!           ------------------------------------------------------------
!
            call CHDiscretization % ComputeGradient(NCOMP, NCOMP, mesh, time, chGradientVariables, element_mask=element_mask)
!
!           --------------------
!           Update MPI Gradients
!           --------------------
!
#ifdef _HAS_MPI_
!$omp single
            call mesh % UpdateMPIFacesGradients(NCOMP)
            call mesh % GatherMPIFacesGradients(NCOMP)
!$omp end single
#endif
!
!           ----------------------
!           Get chemical potential
!           ----------------------
!
!           Get the concentration Laplacian (into QDot => cDot)

            call ComputeLaplacian(mesh, time, element_mask)
!
!           ------------------------------------------
!           *** WARNING! The storage leaves set to CH!
!           ------------------------------------------
!
!$omp single
            call mesh % SetStorageToEqn(C_BC)
            call SetBoundaryConditionsEqn(C_BC)
!$omp end single
         end select
!
!///////////////////////////////////////////////
!        3rd step: Navier-Stokes time derivative
!///////////////////////////////////////////////
!
         select case (mode)
         case (CTD_IGNORE_MODE, CTD_IMEX_EXPLICIT)
!$omp single         
         call mesh % SetStorageToEqn(NS_BC)
         call SetBoundaryConditionsEqn(NS_BC)
!$omp end single
!
!        -------------------------
!        Prolong solution to faces        
!        -------------------------
!
         call mesh % ProlongSolutionToFaces(NCONS, element_mask=element_mask, Level=locLevel)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution(NCONS)
         call mesh % GatherMPIFacesSolution(NCONS)
!$omp end single
#endif
!
!        -------------------------------------
!        Get the density and invMa2 in faces and elements
!        -------------------------------------
!
!$omp do schedule(runtime) private(eID)
         do lID = 1, MLIter(locLevel,1)
		    eID = MLIter_eID(lID)
            compute_element = .true.
            if (present(element_mask)) compute_element = element_mask(eID)
            
            if (compute_element) then

               mesh % elements(eID) % storage % rho = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2)) * mesh % elements(eID) % storage % Q(IMC,:,:,:)
               mesh % elements(eID) % storage % rho = min(max(mesh % elements(eID) % storage % rho, dimensionless % rho_min),dimensionless % rho_max)

               if (use_non_constant_speed_of_sound ) then
                  mesh % elements(eID) % storage % invMa2 = (sqrt(dimensionless % invMa2(1)/dimensionless % rho(1)) * min(max(mesh % elements(eID) % storage % Q(IMC,:,:,:),0.0_RP),1.0_RP) + sqrt(dimensionless % invMa2(2)/dimensionless % rho(2)) * (1.0_RP - min(max(mesh % elements(eID) % storage % Q(IMC,:,:,:),0.0_RP),1.0_RP)))**2 
                  mesh % elements(eID) % storage % invMa2 = mesh % elements(eID) % storage % invMa2*mesh % elements(eID) % storage % rho 
               else
                  mesh % elements(eID) % storage % invMa2 = dimensionless % invMa2(1)
               endif

            endif 

         end do
!$omp end do nowait

!$omp do schedule(runtime) private(fID)
         do lID = 1, MLIter(locLevel,2)
			fID = MLIter_fID(lID)
            compute_element = .true.
            if (present(element_mask)) compute_element = face_mask(fID)
            
            if (compute_element) then
			   
               mesh % faces(fID) % storage(1) % rho = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2))* mesh % faces(fID) % storage(1) % Q(IMC,:,:)
               mesh % faces(fID) % storage(2) % rho = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2))* mesh % faces(fID) % storage(2) % Q(IMC,:,:)

               mesh % faces(fID) % storage(1) % rho = min(max(mesh % faces(fID) % storage(1) % rho, dimensionless % rho_min),dimensionless % rho_max)
               mesh % faces(fID) % storage(2) % rho = min(max(mesh % faces(fID) % storage(2) % rho, dimensionless % rho_min),dimensionless % rho_max)


               if (use_non_constant_speed_of_sound ) then
                  mesh % faces(fID) % storage(1) % invMa2 = (sqrt(dimensionless % invMa2(1)/dimensionless % rho(1)) * min(max(mesh % faces(fID) % storage(1) % Q(IMC,:,:),0.0_RP),1.0_RP) + sqrt(dimensionless % invMa2(2)/dimensionless % rho(2)) * (1.0_RP - min(max(mesh % faces(fID) % storage(1) % Q(IMC,:,:),0.0_RP),1.0_RP)))**2
                  mesh % faces(fID) % storage(2) % invMa2 = (sqrt(dimensionless % invMa2(1)/dimensionless % rho(1)) * min(max(mesh % faces(fID) % storage(2) % Q(IMC,:,:),0.0_RP),1.0_RP) + sqrt(dimensionless % invMa2(2)/dimensionless % rho(2)) * (1.0_RP - min(max(mesh % faces(fID) % storage(2) % Q(IMC,:,:),0.0_RP),1.0_RP)))**2 

                  mesh % faces(fID) % storage(1) % invMa2 = mesh % faces(fID) % storage(1) % invMa2*mesh % faces(fID) % storage(1) % rho
                  mesh % faces(fID) % storage(2) % invMa2 = mesh % faces(fID) % storage(2) % invMa2*mesh % faces(fID) % storage(2) % rho

               else
                  mesh % faces(fID) % storage(1) % invMa2 = dimensionless % invMa2(1)
                  mesh % faces(fID) % storage(2) % invMa2 = dimensionless % invMa2(2)
               endif

            endif
         end do
!$omp end do
!
!        ----------------------------------------
!        Compute local entropy variables gradient
!        ----------------------------------------
!
         call ViscousDiscretization % ComputeLocalGradients( NCONS, NCONS, mesh , time , mGradientVariables, element_mask=element_mask,  Level = locLevel)
!
!        --------------------
!        Update MPI Gradients
!        --------------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesGradients(NCONS)
!$omp end single
#endif
!
!        -------------------------------------
!        Add the Non-Conservative term to QDot
!        -------------------------------------
!
!$omp do schedule(runtime) private(i,j,k,sqrtRho,invMa2,eID)
		 do lID = 1, MLIter(locLevel,1)  
		    eID = MLIter_eID(lID)
            compute_element = .true.
            if (present(element_mask)) compute_element = element_mask(eID)
            
            if (compute_element) then

               associate(e => mesh % elements(eID))
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  sqrtRho = sqrt(e % storage % rho(i,j,k))

                  invMa2 = e % storage % invMa2(i,j,k)

                  e % storage % QDot(IMC,i,j,k)      = 0.0_RP
                  e % storage % QDot(IMSQRHOU,i,j,k) = -0.5_RP*sqrtRho*(  e % storage % Q(IMSQRHOU,i,j,k)*e % storage % U_x(IGU,i,j,k) & 
                                                                        + e % storage % Q(IMSQRHOV,i,j,k)*e % storage % U_y(IGU,i,j,k) &   
                                                                        + e % storage % Q(IMSQRHOW,i,j,k)*e % storage % U_z(IGU,i,j,k) ) &
                                                       - e % storage % Q(IMC,i,j,k)*e % storage % U_x(IGMU,i,j,k)

                  e % storage % QDot(IMSQRHOV,i,j,k) = -0.5_RP*sqrtRho*(  e % storage % Q(IMSQRHOU,i,j,k)*e % storage % U_x(IGV,i,j,k) & 
                                                                        + e % storage % Q(IMSQRHOV,i,j,k)*e % storage % U_y(IGV,i,j,k) &   
                                                                        + e % storage % Q(IMSQRHOW,i,j,k)*e % storage % U_z(IGV,i,j,k) ) &
                                                       - e % storage % Q(IMC,i,j,k)*e % storage % U_y(IGMU,i,j,k)

                  e % storage % QDot(IMSQRHOW,i,j,k) = -0.5_RP*sqrtRho*(  e % storage % Q(IMSQRHOU,i,j,k)*e % storage % U_x(IGW,i,j,k) & 
                                                                        + e % storage % Q(IMSQRHOV,i,j,k)*e % storage % U_y(IGW,i,j,k) &   
                                                                        + e % storage % Q(IMSQRHOW,i,j,k)*e % storage % U_z(IGW,i,j,k) ) &
                                                       - e % storage % Q(IMC,i,j,k)*e % storage % U_z(IGMU,i,j,k)
               
                  ! e % storage % QDot(IMP,i,j,k) = -dimensionless % invMa2*(  e % storage % U_x(IGU,i,j,k) + e % storage % U_y(IGV,i,j,k) &
                  !                                                            + e % storage % U_z(IGW,i,j,k))    
                   e % storage % QDot(IMP,i,j,k) = - invMa2*(  e % storage % U_x(IGU,i,j,k) + e % storage % U_y(IGV,i,j,k) &
                                                                              + e % storage % U_z(IGW,i,j,k))                                                                                                                
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) * e % geom % jacobian(i,j,k)
               end do                ; end do                ; end do
               end associate
         endif 
         end do
!$omp end do

         call ViscousDiscretization % LiftGradients( NCONS, NCONS, mesh , time , mGradientVariables, element_mask=element_mask)

#ifdef _HAS_MPI_
!$omp single
         ! Not sure about the position of this w.r.t the MPI directly above
         call mesh % UpdateMPIFacesGradients(NCONS)
         call mesh % GatherMPIFacesGradients(NCONS)
!$omp end single
#endif   
!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         select case (mode)
         case(CTD_IMEX_EXPLICIT)
            call multiphase % SetStarMobility(0.0_RP)
         case(CTD_IGNORE_MODE)
            call multiphase % SetStarMobility(multiphase % M0)
         end select

         call ComputeNSTimeDerivative(mesh, time, element_mask, Level = locLevel)

         call multiphase % SetStarMobility(multiphase % M0)

         end select
!
!        -------------------------------------------------------------------------------
!        If IMEX_Explicit, compute cDot with the explicit part of the chemical potential
!        -------------------------------------------------------------------------------
!
         select case (mode)
         case(CTD_IMEX_EXPLICIT)
!$omp do schedule(runtime)
            do eID = 1, size(mesh % elements)
               compute_element = .true.
               if (present(element_mask)) compute_element = element_mask(eID)
               
               if (compute_element) then
!
!              + Linear part
                 mesh % elements(eID) % storage % mu = - IMEX_S0 * mesh % elements(eID) % storage % c &
                                                       - 1.5_RP*(1.0_RP - IMEX_K0)*multiphase % sigma*multiphase % eps*mesh % elements(eID) % storage % cDot
                 !mesh % elements(eID) % storage % mu = - IMEX_S0 * mesh % elements(eID) % storage % c &
                 !                                      - (1.0_RP - IMEX_K0)* POW2(multiphase % eps)*mesh % elements(eID) % storage % cDot
!  
!              + NonLinear part
                 !call AddQuarticDWPDerivative(mesh % elements(eID) % storage % c, mesh % elements(eID) % storage % mu)
                  call Multiphase_AddChemFEDerivative(mesh % elements(eID) % storage % c, mesh % elements(eID) % storage % mu)
               endif
            end do
!$omp end do         
!
!           -----------------------------------
!           Prolong chemical potential to faces
!           -----------------------------------
!
!$omp single
            call mesh % SetStorageToEqn(MU_BC)
            call SetBoundaryConditionsEqn(MU_BC)
!$omp end single
         call mesh % ProlongSolutionToFaces(NCOMP, .false., element_mask)
!
!           ------------------------------------------------------------
!           Get concentration (lifted) gradients (also prolong to faces)
!           ------------------------------------------------------------
!
            call CHDiscretization % ComputeGradient(NCOMP, NCOMP, mesh, time, chGradientVariables, .false., element_mask = element_mask)
!
!           --------------------------------
!           Get chemical potential laplacian
!           --------------------------------
!
!           Get the concentration Laplacian (into QDot => cDot)

            call ComputeLaplacian(mesh, time, element_mask)

!$omp single
            call mesh % SetStorageToEqn(NS_BC)
            call SetBoundaryConditionsEqn(NS_BC)
!$omp end single
!
!           -----------------------------------------
!           Add the Chemical potential to the NS QDot
!           -----------------------------------------
!
!$omp do schedule(runtime)
            do eID = 1, size(mesh % elements)
               compute_element = .true.
               if (present(element_mask)) compute_element = element_mask(eID)
               
               if (compute_element) then
                  mesh % elements(eID) % storage % QDot(IMC,:,:,:) =   mesh % elements(eID) % storage % QDot(IMC,:,:,:) &
                                                                     + multiphase % M0*mesh % elements(eID) % storage % cDot(1,:,:,:)
               endif
            end do
!$omp end do
         end select
!$omp end parallel
	  end associate
         if (present(element_mask)) deallocate(face_mask)
!
      END SUBROUTINE ComputeTimeDerivative
!
!////////////////////////////////////////////////////////////////////////////////////
!
!           Navier--Stokes procedures
!           -------------------------
!
!////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ComputeNSTimeDerivative( mesh , t, element_mask, Level )
         use SpongeClass, only: sponge, addSourceSponge
         use ActuatorLine, only: farm, ForcesFarm
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
		 integer, intent(in), optional   :: Level
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
         logical, intent(in), optional   :: element_mask(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID, m, iFace, iEl, locLevel,lID
         real(kind=RP) :: sqrtRho, invSqrtRho
         real(kind=RP)  :: mu_smag, delta
         real(kind=RP), dimension(NCONS)  :: Source
         logical        :: compute_element
         logical, allocatable :: face_mask(:)
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

         if (present(element_mask)) then
            call CreateFaceMask(mesh, element_mask, face_mask)
         endif     

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
									e % storage % mu_NS(1,i,j,k) = e % storage % mu_turb_NS(i,j,k)

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
            compute_element = .true.
            if (present(element_mask)) compute_element = element_mask(eID)
            
            if (compute_element) then
               call TimeDerivative_VolumetricContribution( mesh % elements(eID) , t)
            endif
         end do
!$omp end do nowait


!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) private(fID)
         do iFace = 1, MLIter(locLevel,3)
		    fID = MLIter_fID_Interior(iFace)
		    compute_element = .true.
			if (present(element_mask)) compute_element = face_mask(fID)
     
            if (compute_element) then
               if (mesh% faces(fID) % MortarType == MORTAR_BIG) then
                  associate(fstar=>mesh% faces(fID)%storage(1)%fStar)
                     fstar=0.0_RP
                  end associate
                  associate(fstar=>mesh% faces(fID)%storage(2)%fStar)
                     fstar=0.0_RP
                  end associate
                  do m=1,4
                     if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                        CALL computeElementInterfaceFlux_MU(masterFace=mesh% faces(fID), f=mesh % faces(mesh % faces(fID)%Mortar(m))) 
                     end if
                  end do 
               elseif (mesh % faces(fID) % MortarType == MORTAR_NONE) then
				      call computeElementInterfaceFlux_MU(mesh % faces(fID))
               end if
			end if 
         end do
!$omp end do nowait

!$omp do schedule(runtime) private(fID)
         do iFace = 1, MLIter(locLevel,4)
		    fID = MLIter_fID_Boundary(iFace)
		    compute_element = .true.
			if (present(element_mask)) compute_element = face_mask(fID)
     
            if (compute_element) then
				call computeBoundaryFlux_MU(mesh % faces(fID), t)	
            end if 				
         end do
!$omp end do
!
!        *************************************************************************************
!        Element without shared faces: Surface integrals, scaling of elements with Jacobian, 
!                                      sqrt(rho)
!        *************************************************************************************
! 
!$omp do schedule(runtime) private(i,j,k,eID,sqrtRho,invSqrtRho)
         do iEl = 1, MLIter(locLevel,5)
            eID = MLIter_eID_Seq(iEl)
            compute_element = .true.
            if (present(element_mask)) compute_element = element_mask(eID)
            
            if (compute_element) then

               associate(e => mesh % elements(eID)) 
			   
               call TimeDerivative_FacesContribution(e, t, mesh) 
               
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
                  sqrtRho = sqrt(e % storage % rho(i,j,k))
                  invSqrtRho = 1.0_RP / sqrtRho
!  
!               + Scale with Jacobian and sqrt(Rho)
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) * e % geom % InvJacobian(i,j,k)
                  e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) = e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) * invSqrtRho
!  
!               + Add gravity
                  e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) =   e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) & 
                                                                + sqrtRho * dimensionless % invFr2 * dimensionless % gravity_dir

               end do         ; end do          ; end do 

               end associate 
            endif
         end do
!$omp end do

!
!        ******************************************
!        Do the same for elements with shared faces
!        ******************************************
!
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
!$omp single
            call mesh % GatherMPIFacesGradients(NCONS)
!$omp end single

!
!           Compute viscosity at MPI faces
!           ------------------------------
            call compute_viscosity_at_faces(MLIter(locLevel,7), 2, MLIter_fID_MPI(1:MLIter(locLevel,7)), mesh)
!
!           **************************************
!           Compute Riemann solver of shared faces
!           **************************************
!
!$omp do schedule(runtime) private(fID)
            do iFace = 1, MLIter(locLevel,7)
               fID = MLIter_fID_MPI(iFace)
               compute_element = .true.
               if (present(element_mask)) compute_element = face_mask(fID)
               
               if (compute_element) then
                  if (mesh% faces(fID)%MortarType == MORTAR_BIG) then 
                     associate(fstar=>mesh% faces(fID)%storage(1)%fStar)
                        fstar=0.0_RP
                     end associate
                     do m=1,4
                        if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                           CALL computeElementInterfaceFlux_MU(masterFace=mesh% faces(fID), f=mesh % faces(mesh % faces(fID)%Mortar(m))) 
                        end if 
                     end do 
                  end if 
                  ! AJRTODO: Should it be here: ifelse MortarType == 0 ? Before you do it
                  CALL computeMPIFaceFlux_MU( mesh % faces(fID) )
               endif
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
!$omp do schedule(runtime) private(i,j,k, eID, sqrtRho, invSqrtRho)
            do iEl = 1, MLIter(locLevel,6)
               eID = MLIter_eID_MPI(iEl)
               compute_element = .true.
               if (present(element_mask)) compute_element = element_mask(eID)
               
               if (compute_element) then
                  associate(e => mesh % elements(eID))
                  call TimeDerivative_FacesContribution(e, t, mesh)
 
                  do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
                     sqrtRho = sqrt(e % storage % rho(i,j,k))
                     invSqrtRho = 1.0_RP / sqrtRho
!   
!                  + Scale with Jacobian and sqrt(Rho)
                     e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) * e % geom % InvJacobian(i,j,k)
                     e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) = e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) * invSqrtRho
!   
!                  + Add gravity
                     e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) =   e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) & 
                                                                   + sqrtRho * dimensionless % invFr2 * dimensionless % gravity_dir
   
                  end do         ; end do          ; end do 
               ! do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               !    e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k)
               ! end do         ; end do          ; end do
                  end associate
               endif
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

!           ***************
!           Add source term
!           ***************
!$omp do schedule(runtime) private(i,j,k, InvSqrtRho, eID)
            do lID = 1, MLIter(locLevel,1)
			   eID = MLIter_eID(lID)
               compute_element = .true.
               if (present(element_mask)) compute_element = element_mask(eID)
               
               if (compute_element) then

                  associate ( e => mesh % elements(eID) )
                  e % storage % S_NS = 0.0_RP
                  do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                     InvSqrtRho = 1.0_RP / sqrt(e % storage % rho(i,j,k))
                     call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues, multiphase)
                     ! scale UserDefinedSourceTerm momentum with sqrtRho
                     e % storage % S_NS(:,i,j,k) = e % storage % S_NS(:,i,j,k) * [1.0_RP,InvSqrtRho,InvSqrtRho,InvSqrtRho,1.0_RP]
                  end do                  ; end do                ; end do
                  end associate
               endif
            end do
!$omp end do

!for the sponge, loops are in the internal subroutine as values are precalculated
!The scale with sqrtRho is done in the subroutines, not done againg here
         call addSourceSponge(sponge,mesh)
         call ForcesFarm(farm, mesh, t, Level=locLevel)

! Add all the source terms
!$omp do schedule(runtime) private(i,j,k, Source, eID)
         do lID = 1, MLIter(locLevel,1)
		    eID = MLIter_eID(lID)
            compute_element = .true.
            if (present(element_mask)) compute_element = element_mask(eID)
            
            if (compute_element) then

               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S_NS(:,i,j,k)
               end do                  ; end do                ; end do
               end associate
            endif
         end do
!$omp end do  
!
!        *********************
!        Add IBM source term
!        *********************
! no wall function for MULTIPHASE
         if( mesh% IBM% active ) then
            if( .not. mesh% IBM% semiImplicit ) then 
!$omp do schedule(runtime) private(i,j,k,Source, eID)
                  do lID = 1, MLIter(locLevel,1)
					 eID = MLIter_eID(lID)
                     compute_element = .true.
                     if (present(element_mask)) compute_element = element_mask(eID)
                     
                     if (compute_element) then

                        associate ( e => mesh % elements(eID) ) 
                        do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                           if( e% isInsideBody(i,j,k) ) then
                              ! only without moving for now in MULTIPHASE
                              if( .not. mesh% IBM% stl(e% STL(i,j,k))% move ) then 
                                 call mesh% IBM% SourceTerm( eID = eID, Q = e % storage % Q(:,i,j,k), Source = Source, wallfunction = .false. )
                              end if 
                              e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + Source
                           end if
                        end do                  ; end do                ; end do
                        end associate
                     endif
                  end do
!$omp end do       
            end if 
         end if

	    end associate
	    end associate

         if (present(element_mask)) deallocate(face_mask)

      end subroutine ComputeNSTimeDerivative
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
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , mEulerFlux, inviscidContravariantFlux ) 
!
!        Compute viscous contravariant flux
!        ----------------------------------
         call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NCONS, mViscousFlux, GetmTwoFluidsViscosity, e , viscousContravariantFlux) 
!
!        ************************
!        Perform volume integrals
!        ************************
!
!
!        Compute the total Navier-Stokes flux
!        ------------------------------------
         contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux 
!
!        Perform the Weak Volume Green integral
!        --------------------------------------
         e % storage % QDot = e % storage % QDot + ScalarWeakIntegrals % StdVolumeGreen ( e, NCONS, contravariantFlux ) 


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
      SUBROUTINE computeElementInterfaceFlux_MU(f, masterFace)
         use FaceClass
         use RiemannSolvers_MU
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         type(Face), optional, intent(inout) :: masterFace 

         integer       :: i, j
         real(kind=RP) :: inv_fluxL(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: inv_fluxR(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: fluxL(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: fluxR(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: muL, muR
         real(kind=RP) :: UxL(1:NGRAD), UyL(1:NGRAD), UzL(1:NGRAD)
         real(kind=RP) :: UxR(1:NGRAD), UyR(1:NGRAD), UzR(1:NGRAD)

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

               call GetmTwoFluidsViscosity(f % storage(1) % Q(IMC,i,j), muL)
               call GetmTwoFluidsViscosity(f % storage(2) % Q(IMC,i,j), muR)

			   muL = muL +f % storage(1) % mu_NS(1,i,j)   ! Add subgrid viscosity
			   muR = muR +f % storage(2) % mu_NS(1,i,j)   ! Add subgrid viscosity	
!
!            - Premultiply velocity gradients by the viscosity
!              -----------------------------------------------
               UxL = [1.0_RP,muL,muL,muL,1.0_RP]*f % storage(1) % U_x(:,i,j) 
               UyL = [1.0_RP,muL,muL,muL,1.0_RP]*f % storage(1) % U_y(:,i,j) 
               UzL = [1.0_RP,muL,muL,muL,1.0_RP]*f % storage(1) % U_z(:,i,j) 

               UxR = [1.0_RP,muR,muR,muR,1.0_RP]*f % storage(2) % U_x(:,i,j) 
               UyR = [1.0_RP,muR,muR,muR,1.0_RP]*f % storage(2) % U_y(:,i,j) 
               UzR = [1.0_RP,muR,muR,muR,1.0_RP]*f % storage(2) % U_z(:,i,j) 
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NCONS, &
                                                  EllipticFlux = mViscousFlux, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = UxL, &
                                                  U_yLeft = UyL, &
                                                  U_zLeft = UzL, &
                                                  U_xRight = UxR, &
                                                  U_yRight = UyR, &
                                                  U_zRight = UzR, &
                                                  mu_left  = [1.0_RP, multiphase % M0_star, 0.0_RP], &
                                                  mu_right = [1.0_RP, multiphase % M0_star, 0.0_RP], &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  sigma = [multiphase % M0_star, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP], &
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
                                  rhoL   = f % storage(1) % rho(i,j), &
                                  rhoR   = f % storage(2) % rho(i,j), &
                                  muL    = f % storage(1) % mu(1,i,j), &
                                  muR    = f % storage(2) % mu(1,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  fL     = inv_fluxL(:,i,j), &
                                  fR     = inv_fluxR(:,i,j), &
                                  invMa2L= f % storage(1) % invMa2(i,j), &
                                  invMa2R= f % storage(2) % invMa2(i,j))

               
!
!              Multiply by the Jacobian
!              ------------------------
               fluxL(:,i,j) = ( inv_fluxL(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               fluxR(:,i,j) = ( inv_fluxR(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         if (f % MortarType == MORTAR_NONE) then 
         call f % ProjectFluxToElements(NCONS, fluxL, (/1, HMESH_NONE/))
         call f % ProjectFluxToElements(NCONS, fluxR, (/2, HMESH_NONE/)) 
      end if 
      if (f % MortarType == MORTAR_SMALL4 .and. present(masterFace)) then 

         call masterFace % ProjectMortarFluxToElements(nEqn=NCONS, whichElements=(/1,0/), slaveFace=f, MortarFlux=fluxL)
         call f % ProjectFluxToElements(NCONS, fluxR, (/2,0/))
      end if 
      END SUBROUTINE computeElementInterfaceFlux_MU

      SUBROUTINE computeMPIFaceFlux_MU(f)
         use FaceClass
         use RiemannSolvers_MU
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_fluxL(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: inv_fluxR(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: fluxL(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: fluxR(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2),2)
         real(kind=RP) :: muL, muR
         real(kind=RP) :: UxL(1:NGRAD), UyL(1:NGRAD), UzL(1:NGRAD)
         real(kind=RP) :: UxR(1:NGRAD), UyR(1:NGRAD), UzR(1:NGRAD)

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

               call GetmTwoFluidsViscosity(f % storage(1) % Q(IMC,i,j), muL)
               call GetmTwoFluidsViscosity(f % storage(2) % Q(IMC,i,j), muR)

			   muL = muL +f % storage(1) % mu_NS(1,i,j)   ! Add subgrid viscosity
			   muR = muR +f % storage(2) % mu_NS(1,i,j)   ! Add subgrid viscosity		
!
!            - Premultiply velocity gradients by the viscosity
!              -----------------------------------------------
               UxL = [1.0_RP,muL,muL,muL,1.0_RP]*f % storage(1) % U_x(:,i,j) 
               UyL = [1.0_RP,muL,muL,muL,1.0_RP]*f % storage(1) % U_y(:,i,j) 
               UzL = [1.0_RP,muL,muL,muL,1.0_RP]*f % storage(1) % U_z(:,i,j) 

               UxR = [1.0_RP,muR,muR,muR,1.0_RP]*f % storage(2) % U_x(:,i,j) 
               UyR = [1.0_RP,muR,muR,muR,1.0_RP]*f % storage(2) % U_y(:,i,j) 
               UzR = [1.0_RP,muR,muR,muR,1.0_RP]*f % storage(2) % U_z(:,i,j) 
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NCONS, &
                                                  EllipticFlux = mViscousFlux, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = UxL, &
                                                  U_yLeft = UyL, &
                                                  U_zLeft = UzL, &
                                                  U_xRight = UxR, &
                                                  U_yRight = UyR, &
                                                  U_zRight = UzR, &
                                                  mu_left  = [1.0_RP, multiphase % M0_star, 0.0_RP], &
                                                  mu_right = [1.0_RP, multiphase % M0_star, 0.0_RP], &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  sigma = [multiphase % M0_star, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP], &
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
                                  rhoL   = f % storage(1) % rho(i,j), &
                                  rhoR   = f % storage(2) % rho(i,j), &
                                  muL    = f % storage(1) % mu(1,i,j), &
                                  muR    = f % storage(2) % mu(1,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  fL     = inv_fluxL(:,i,j), &
                                  fR     = inv_fluxR(:,i,j),&
                                  invMa2L= f % storage(1) % invMa2(i,j), &
                                  invMa2R= f % storage(2) % invMa2(i,j)) 

               
!
!              Multiply by the Jacobian
!              ------------------------
               fluxL(:,i,j) = ( inv_fluxL(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               fluxR(:,i,j) = ( inv_fluxR(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         flux(:,:,:,1) = fluxL
         flux(:,:,:,2) = fluxR

         call f % ProjectFluxToElements(NCONS, flux(:,:,:,thisSide), (/thisSide, HMESH_NONE/))
         if (f % MortarType == MORTAR_SMALL4) then 
            call f% Interpolatesmall2big(NCONS, fluxL)
         
         end if 
      end subroutine ComputeMPIFaceFlux_MU

      SUBROUTINE computeBoundaryFlux_MU(f, time)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_MU
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
      REAL(KIND=RP)                   :: inv_fluxL(NCONS), inv_fluxR(NCONS), fv_3d(NCONS,NDIM)
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
         f % storage(2) % mu(1,i,j) = f % storage(1) % mu(1,i,j)
         CALL BCs(f % zone) % bc % FlowState( &
                                      f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j))

         f % storage(2) % rho(i,j) = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2))*f % storage(2) % Q(IMC,i,j)
         f % storage(2) % rho(i,j) = min(max(f % storage(2) % rho(i,j), dimensionless % rho_min),dimensionless % rho_max)
!   
!        --------------
!        Viscous fluxes
!        --------------
!   
         call GetmTwoFluidsViscosity(f % storage(1) % Q(IMC,i,j), mu)
		 mu = mu + f % storage(1) % mu_NS(1,i,j)   ! Add subgrid viscosity	

         call mViscousFlux(NCONS, NCONS, f % storage(1) % Q(:,i,j), &
                                         f % storage(1) % U_x(:,i,j), &
                                         f % storage(1) % U_y(:,i,j), &
                                         f % storage(1) % U_z(:,i,j), &
                                         mu, multiphase % M0_star, 0.0_RP, fv_3d)

         visc_flux(:,i,j) =   fv_3d(:,IX)*f % geom % normal(IX,i,j) &
                            + fv_3d(:,IY)*f % geom % normal(IY,i,j) &
                            + fv_3d(:,IZ)*f % geom % normal(IZ,i,j) 


         CALL BCs(f % zone) % bc % FlowNeumann(&
                                           f % geom % x(:,i,j), &
                                           time, &
                                           f % geom % normal(:,i,j), &
                                           f % storage(1) % Q(:,i,j), &
                                           f % storage(1) % U_x(:,i,j), &
                                           f % storage(1) % U_y(:,i,j), &
                                           f % storage(1) % U_z(:,i,j), visc_flux(:,i,j))

!
!           Hyperbolic part
!           -------------
            CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                               QRight = f % storage(2) % Q(:,i,j), &
                               rhoL   = f % storage(1) % rho(i,j), &
                               rhoR   = f % storage(2) % rho(i,j), &
                               muL    = f % storage(1) % mu(1,i,j), &
                               muR    = f % storage(2) % mu(1,i,j), &
                               nHat   = f % geom % normal(:,i,j), &
                               t1     = f % geom % t1(:,i,j), &
                               t2     = f % geom % t2(:,i,j), &
                               fL     = inv_fluxL, &
                               fR     = inv_fluxR,&
                               invMa2L= f % storage(1) % invMa2(i,j), &
                               invMa2R= f % storage(2) % invMa2(i,j)) 

            fStar(:,i,j) = (inv_fluxL - visc_flux(:,i,j) ) * f % geom % jacobian(i,j)
         end do   
      end do   

      call f % ProjectFluxToElements(NCONS, fStar, (/1, HMESH_NONE/))

      END SUBROUTINE computeBoundaryFlux_MU
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!        Laplacian procedures
!        --------------------
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ComputeLaplacian( mesh , t, element_mask, Level)
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
         logical, intent(in), optional   :: element_mask(:)
		 integer, intent(in), optional   :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID, m, iFace, iEl, lID, locLevel
         logical     :: compute_element
         logical, allocatable :: face_mask(:)

         if (present(element_mask)) then
            call CreateFaceMask(mesh, element_mask, face_mask)
         endif 

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
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) private(eID)
		 do lID = 1, MLIter(locLevel,1)
		    eID = MLIter_eID(lID)
            compute_element = .true.
            if (present(element_mask)) compute_element = element_mask(eID)
            
            if (compute_element) then
               call Laplacian_VolumetricContribution( mesh % elements(eID) , t)
            endif
         end do
!$omp end do nowait
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) private(fID)
         do iFace = 1, MLIter(locLevel,3)
            fID = MLIter_fID_Interior(iFace)
			compute_element = .true.
            if (present(element_mask)) compute_element = face_mask(fID)
            
            if (compute_element) then
               if(mesh% faces(fID) % MortarType == MORTAR_BIG) then 
                  associate(fstar=>mesh% faces(fID)%storage(1)%fStar)
                     fstar=0.0_RP
                  end associate
                  associate(fstar=>mesh% faces(fID)%storage(2)%fStar)
                     fstar=0.0_RP
                  end associate
                  do m=1,4
                     if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                        CALL Laplacian_computeElementInterfaceFlux(masterFace=mesh % faces(fID), f=mesh % faces(mesh % faces(fID)%Mortar(m)) ) 
                     end if 
                  end do 
               elseif (mesh % faces(fID) % MortarType == MORTAR_NONE) then
                  call Laplacian_computeElementInterfaceFlux(mesh % faces(fID))
               end if
			end if 
         end do
!$omp end do nowait

!$omp do schedule(runtime) private(fID)
         do iFace = 1, MLIter(locLevel,4)
            fID = MLIter_fID_Boundary(iFace)
			compute_element = .true.
            if (present(element_mask)) compute_element = face_mask(fID)
            
            if (compute_element) then
				call Laplacian_computeBoundaryFlux(mesh % faces(fID), t)
			end if 
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
            compute_element = .true.
            if (present(element_mask)) compute_element = element_mask(eID)
            
            if (compute_element) then

               associate(e => mesh % elements(eID)) 
               call Laplacian_FacesContribution(e, t, mesh) 

               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k)
               end do         ; end do          ; end do 
               end associate 
         endif
         end do
!$omp end do
!
!        ***********************************************************
!        Surface integrals and scaling of elements with shared faces
!        ***********************************************************
! 
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
!$omp single
            call mesh % GatherMPIFacesGradients(NCOMP)
!$omp end single
!
!           **************************************
!           Compute Riemann solver of shared faces
!           **************************************
!
!$omp do schedule(runtime) private(fID)
            do iFace = 1, MLIter(locLevel,7)
               fID = MLIter_fID_MPI(iFace)
               compute_element = .true.
               if (present(element_mask)) compute_element = face_mask(fID)
               
               if (compute_element) then
                  if (mesh % faces(fID)%MortarType == MORTAR_BIG) then 
                     associate(fstar=>mesh% faces(fID)%storage(1)%fStar)
                        fstar=0.0_RP
                     end associate
                     do m=1,4
                        if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                           CALL Laplacian_computeElementInterfaceFlux(masterFace=mesh % faces(fID), f=mesh % faces(mesh % faces(fID)%Mortar(m)) ) 
                        end if 
                     end do 
                  end if
                  ! AJRTODO: ifelse MortarType==0 as before ??
                  call Laplacian_computeMPIFaceFlux(mesh % faces(fID))
               end if
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
            do iEl = 1, MLIter(locLevel,6)
               eID = MLIter_eID_MPI(iEl)
               compute_element = .true.
               if (present(element_mask)) compute_element = element_mask(eID)
               
               if (compute_element) then
                  associate(e => mesh % elements(eID))
				  
                  call Laplacian_FacesContribution(e, t, mesh)

                  do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                     e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k)
                  end do         ; end do          ; end do
                  end associate
               endif
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
      end associate
	  end associate		
	  
      if (present(element_mask)) deallocate(face_mask)

      end subroutine ComputeLaplacian

      subroutine ComputeLaplacianNeumannBCs( mesh , t)
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID
!
!        **************************
!        Reset QDot and face fluxes
!        **************************
!
         do eID = 1, mesh % no_of_elements
            mesh % elements(eID) % storage % QDot = 0.0_RP
         end do
   
         do fID = 1, size(mesh % faces)
            mesh % faces(fID) % storage(1) % genericInterfaceFluxMemory = 0.0_RP
            mesh % faces(fID) % storage(2) % genericInterfaceFluxMemory = 0.0_RP
         end do
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) 
         do fID = 1, size(mesh % faces) 
            associate( f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_BOUNDARY) 
               CALL Laplacian_computeBoundaryFlux(f, t) 
            end select 
            end associate 
         end do 
!$omp end do 
!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
! 
!$omp do schedule(runtime) private(i, j, k)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            if ( e % hasSharedFaces ) cycle
            call Laplacian_FacesContribution(e, t, mesh) 
 
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k) 
            end do         ; end do          ; end do 
            end associate 
         end do
!$omp end do

      end subroutine ComputeLaplacianNeumannBCs
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Laplacian_VolumetricContribution( e , t )
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
         real(kind=RP) :: contravariantFlux  ( 1:NCOMP, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute contravariant flux
!        --------------------------
         e % storage % mu_NS(1,:,:,:) = 0.0_RP ! Need to zero out for Laplacian ComputeInnerFluxes
         call CHDiscretization  % ComputeInnerFluxes (NCOMP, NCOMP, CHDivergenceFlux, GetCHViscosity, e , contravariantFlux  ) 
!
!        ************************
!        Perform volume integrals
!        ************************
!
         e % storage % QDot = - ScalarWeakIntegrals % StdVolumeGreen ( e , NCOMP, contravariantFlux ) 

      end subroutine Laplacian_VolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Laplacian_FacesContribution( e , t , mesh)
         use HexMeshClass
         use PhysicsStorage
         implicit none
         type(Element)           :: e
         real(kind=RP)           :: t
         type(HexMesh)           :: mesh

         e % storage % QDot = e % storage % QDot + ScalarWeakIntegrals % StdFace(e, NCOMP, &
                      mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % fStar, &
                      mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % fStar, &
                      mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % fStar, &
                      mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % fStar, &
                      mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % fStar, &
                      mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % fStar )

      end subroutine Laplacian_FacesContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
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


         if ( LESModel % Active ) then
!$omp do schedule(runtime) private(i,j,delta,mu_smag,side)
            do iFace = 1, no_of_faces
               associate(f => mesh % faces(face_ids(iFace)))

               delta = sqrt(f % geom % surface / product(f % Nf + 1))
			  
               do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
                  do side = 1, no_of_sides
                     call LESModel % ComputeViscosity(delta, f % geom % dWall(i,j), f % storage(side) % Q(:,i,j),   &
                                                                                    f % storage(side) % U_x(:,i,j), &
                                                                                    f % storage(side) % U_y(:,i,j), &
                                                                                    f % storage(side) % U_z(:,i,j), &
                                                                                    mu_smag)
                     f % storage(side) % mu_NS(1,i,j) =  mu_smag
																											
                  end do
               end do              ; end do
               end associate
            end do
!$omp end do
         end if

      end subroutine compute_viscosity_at_faces
!
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
!        Riemann solver drivers 
!        ---------------------- 
! 
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
      subroutine Laplacian_computeElementInterfaceFlux(f, masterFace)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f  
         type(Face), optional, intent(inout) :: masterFace 

         integer       :: i, j, m
         real(kind=RP) :: flux(1:NCOMP,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               call GetCHViscosity(0.0_RP, mu)
               CALL CHDiscretization % RiemannSolver(nEqn = NCOMP, nGradEqn = NCOMP, &
                                                  EllipticFlux = CHDivergenceFlux, &
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
                                                  sigma = [1.0_RP], &
                                                  flux  = flux(:,i,j) )

               flux(:,i,j) = flux(:,i,j) * f % geom % jacobian(i,j)

            end do
         end do
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         if (f % MortarType == MORTAR_NONE) then 
         call f % ProjectFluxToElements(NCOMP, flux, (/1,2/))
      end if 
      if (f % MortarType == MORTAR_SMALL4 .and. present(masterFace)) then 
         call masterFace % ProjectMortarFluxToElements(nEqn=NCONS, whichElements=(/1,0/), slaveFace=f, MortarFlux=flux)
         call f % ProjectFluxToElements(NCONS, flux, (/0,2/))
        end if 
      end subroutine Laplacian_computeElementInterfaceFlux

      subroutine Laplacian_computeMPIFaceFlux(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: flux(1:NCOMP,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               call GetCHViscosity(0.0_RP, mu)
               CALL CHDiscretization % RiemannSolver(nEqn = NCOMP, nGradEqn = NCOMP, &
                                                  EllipticFlux = CHDivergenceFlux, &
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
                                                  sigma = [1.0_RP], &
                                                  flux  = flux(:,i,j) )

               flux(:,i,j) = flux(:,i,j) * f % geom % jacobian(i,j)

            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectFluxToElements(NCOMP, flux, (/thisSide, HMESH_NONE/))

         if (f % MortarType == MORTAR_SMALL4) then 
            call f% Interpolatesmall2big(NCONS, flux)
            
         end if 
      end subroutine Laplacian_ComputeMPIFaceFlux

      subroutine Laplacian_computeBoundaryFlux(f, time)
      USE ElementClass
      use FaceClass
      USE EllipticDiscretizations
      USE Physics
      use PhysicsStorage
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
      real(kind=RP)                   :: flux(NCOMP, 0:f % Nf(1), 0:f % Nf(2)), fv_3d(NCOMP,NDIM)
      real(kind=RP)                   :: mu
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
!
!        --------------
!        Viscous fluxes
!        --------------
!   
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         call GetCHViscosity(0.0_RP, mu)
         call CHDivergenceFlux(NCOMP, NCOMP, f % storage(1) % Q(:,i,j), &
                                             f % storage(1) % U_x(:,i,j), &
                                             f % storage(1) % U_y(:,i,j), &
                                             f % storage(1) % U_z(:,i,j), &
                                             mu, 0.0_RP, 0.0_RP, fv_3d)

         flux(:,i,j) =   fv_3d(:,IX)*f % geom % normal(IX,i,j) &
                       + fv_3d(:,IY)*f % geom % normal(IY,i,j) &
                       + fv_3d(:,IZ)*f % geom % normal(IZ,i,j) 

         CALL BCs(f % zone) % bc % NeumannForEqn(NCOMP, NCOMP, &
                                           f % geom % x(:,i,j), &
                                           time, &
                                           f % geom % normal(:,i,j), &
                                           f % storage(1) % Q(:,i,j), &
                                           f % storage(1) % U_x(:,i,j), &
                                           f % storage(1) % U_y(:,i,j), &
                                           f % storage(1) % U_z(:,i,j), flux(:,i,j))

         flux(:,i,j) = flux(:,i,j) * f % geom % jacobian(i,j)

      end do               ; end do

      call f % ProjectFluxToElements(NCOMP, flux, (/1, HMESH_NONE/))

      end subroutine Laplacian_computeBoundaryFlux

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

      end subroutine ComputeTimeDerivativeIsolated

      subroutine CreateFaceMask(mesh, element_mask, face_mask)
         use omp_lib
         implicit none
         type(HexMesh), intent(in)         :: mesh
         logical, intent(in)               :: element_mask(:)
         logical, allocatable, intent(out) :: face_mask(:)
     
         integer :: fID
         integer :: e1, e2
     
         allocate(face_mask(size(mesh % faces)))
         face_mask = .false.
     
         !$omp parallel do private(fID, e1, e2) schedule(runtime)
         do fID = 1, size(mesh % faces)
             associate(f => mesh % faces(fID))
                 e1 = f % elementIDs(1)
                 e2 = f % elementIDs(2)
     
                 if (e1 > 0 .and. e2 > 0) then
                     face_mask(fID) = element_mask(e1) .or. element_mask(e2)
                 else if (e1 > 0) then
                     face_mask(fID) = element_mask(e1)
                 else if (e2 > 0) then
                     face_mask(fID) = element_mask(e2)
                 end if
             end associate
         end do
         !$omp end parallel do
     
     end subroutine CreateFaceMask
     
end module SpatialDiscretization
