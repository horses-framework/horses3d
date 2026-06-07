!
!////////////////////////////////////////////////////////////////////////
!
!      Basic class for the discontinuous Galerkin spectral element
!      solution of conservation laws.
!
!////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
Module DGSEMClass
   use SMConstants
   USE NodalStorageClass         , only: CurrentNodes, NodalStorage, NodalStorage_t
   use MeshTypes                 , only: HOPRMESH
   use ElementClass
   USE HexMeshClass
   USE SlidingMeshProcedures
   USE PhysicsStorage
   use FileReadingUtilities      , only: getFileName
   use MPI_Process_Info          , only: MPI_Process
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
   use ManufacturedSolutionsNS
   use FWHGeneralClass
#elif defined(SPALARTALMARAS)
   use ManufacturedSolutionsNSSA
   use SpallartAlmarasTurbulence , only: Spalart_Almaras_t
#endif
   use MonitorsClass
   use ParticlesClass
   use Physics
   use FluidData
   use ProblemFileFunctions, only: UserDefinedInitialCondition_f
   use MPI_Utilities,        only: MPI_MinMax
#ifdef _HAS_MPI_
   use mpi
#endif

   IMPLICIT NONE

   private
   public   ComputeTimeDerivative_f, DGSem, ConstructDGSem
   public   DestructDGSEM, MaxTimeStep, ComputeMaxResiduals
   public   hnRange

   TYPE DGSem
      REAL(KIND=RP)                                           :: maxResidual
      integer                                                 :: nodes                 ! Either GAUSS or GAUSLOBATTO
      INTEGER                                                 :: numberOfTimeSteps
      INTEGER                                                 :: NDOF                        ! Number of degrees of freedom in this partition
      integer                                                 :: totalNDOF                   ! Number of degrees of freedom in the whole domain
      TYPE(HexMesh)                                           :: mesh
      LOGICAL                                                 :: ManufacturedSol = .FALSE.   ! Use manufactured solutions? default .FALSE.
      type(Monitor_t)                                         :: monitors
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
      type(FWHClass)                                          :: fwh
#endif
#if defined(SPALARTALMARAS)
      type(Spalart_Almaras_t), private    :: SAModel
#endif
#ifdef FLOW
      type(Particles_t)                                       :: particles
#else
      logical                                                 :: particles
#endif
      contains
         procedure :: construct => ConstructDGSem
         procedure :: destruct  => DestructDGSem
         procedure :: SaveSolutionForRestart
         procedure :: SetInitialCondition => DGSEM_SetInitialCondition
         procedure :: copy                => DGSEM_Assign
         generic   :: assignment(=)       => copy
   END TYPE DGSem

   abstract interface
      SUBROUTINE ComputeTimeDerivative_f( mesh, particles, time, mode, HO_Elements)
         use SMConstants
         use HexMeshClass
         use ParticlesClass
         IMPLICIT NONE
         type(HexMesh), target           :: mesh
#ifdef FLOW
         type(Particles_t)               :: particles
#else
         logical                         :: particles
#endif
         REAL(KIND=RP)                   :: time
         integer,             intent(in) :: mode
         logical, intent(in), optional   :: HO_Elements
      end subroutine ComputeTimeDerivative_f
   END INTERFACE

   CONTAINS
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructDGSem( self, meshFileName_, controlVariables, &
                                 polynomialOrder, Nx_, Ny_, Nz_, success, ChildSem )
      use ReadMeshFile
      use FTValueDictionaryClass
      use mainKeywordsModule
      use StopwatchClass
      use MPI_Process_Info
      use PartitionedMeshClass
      use MeshPartitioning
      use SurfaceMesh, only: surfacesMesh

      IMPLICIT NONE
!
!     --------------------------
!     Constructor for the class.
!     --------------------------
!
      CLASS(DGSem)                       :: self                               !<> Class to be constructed
      character(len=*),         optional :: meshFileName_
      class(FTValueDictionary)           :: controlVariables                   !<  Name of mesh file
      INTEGER, OPTIONAL                  :: polynomialOrder(3)                 !<  Uniform polynomial order
      INTEGER, OPTIONAL, TARGET          :: Nx_(:), Ny_(:), Nz_(:)             !<  Non-uniform polynomial order
      LOGICAL, OPTIONAL                  :: success                            !>  Construction finalized correctly?
      logical, optional                  :: ChildSem                           !<  Is this a (multigrid) child sem?
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                     :: i,j,k,el,bcset                     ! Counters
      INTEGER, POINTER            :: Nx(:), Ny(:), Nz(:)                ! Orders of every element in mesh (used as pointer to use less space)
      integer                     :: NelL(2), NelR(2)
      INTEGER                     :: nTotalElem                              ! Number of elements in mesh
      INTEGER                     :: fUnit
      integer                     :: dir2D
      integer                     :: ierr
      logical                     :: MeshInnerCurves                    ! The inner survaces of the mesh have curves?
      logical                     :: useRelaxPeriodic                   ! The periodic construction in direction z use a relative tolerance
      logical                     :: useWeightsPartition                ! Partitioning mesh using DOF of elements as weights
      real(kind=RP)               :: QbaseUniform(1:NCONS)
      character(len=*), parameter :: TWOD_OFFSET_DIR_KEY = "2d mesh offset direction"
#if (!defined(NAVIERSTOKES))
      logical, parameter          :: computeGradients = .true.
#endif

      if ( present(ChildSem) ) then
         if ( ChildSem ) self % mesh % child = .TRUE.
      end if

!
!     Measure preprocessing time
!     --------------------------
      if (.not. self % mesh % child) then
         call Stopwatch % CreateNewEvent("Preprocessing")
         call Stopwatch % Start("Preprocessing")
      end if

      if ( present( meshFileName_ ) ) then
!
!        Mesh file set up by input argument
!        ----------------------------------
         self % mesh % meshFileName = trim(meshFileName_)

      else
!
!        Mesh file set up by controlVariables
!        ------------------------------------
         self % mesh % meshFileName = controlVariables % stringValueForKey(meshFileNameKey, requestedLength = LINE_LENGTH)

      end if
!
!     ------------------------------
!     Discretization nodes selection
!     ------------------------------
!
      self % nodes = CurrentNodes
!
!     ---------------------------------------
!     Get polynomial orders for every element
!     ---------------------------------------
!
      IF (PRESENT(Nx_) .AND. PRESENT(Ny_) .AND. PRESENT(Nz_)) THEN
         Nx => Nx_
         Ny => Ny_
         Nz => Nz_
         nTotalElem = SIZE(Nx)
      ELSEIF (PRESENT(polynomialOrder)) THEN
         nTotalElem = NumOfElemsFromMeshFile( self % mesh % meshfileName )

         ALLOCATE (Nx(nTotalElem),Ny(nTotalElem),Nz(nTotalElem))
         Nx = polynomialOrder(1)
         Ny = polynomialOrder(2)
         Nz = polynomialOrder(3)
      ELSE
         error stop 'ConstructDGSEM: Polynomial order not specified'
      END IF

      if ( max(maxval(Nx),maxval(Ny),maxval(Nz)) /= min(minval(Nx),minval(Ny),minval(Nz)) ) self % mesh % anisotropic = .TRUE.

!
!     -------------------------------------------------------------
!     Construct the polynomial storage for the elements in the mesh
!     -------------------------------------------------------------
!
      call NodalStorage(0) % Construct(CurrentNodes, 0)   ! Always construct orders 0
      call NodalStorage(1) % Construct(CurrentNodes, 1)   ! and 1

      DO k=1, nTotalElem
         call NodalStorage(Nx(k)) % construct( CurrentNodes, Nx(k) )
         call NodalStorage(Ny(k)) % construct( CurrentNodes, Ny(k) )
         call NodalStorage(Nz(k)) % construct( CurrentNodes, Nz(k) )
      END DO
!
!     ------------------
!     Construct the mesh
!     ------------------
!
      if ( controlVariables % containsKey(TWOD_OFFSET_DIR_KEY) ) then
         select case ( controlVariables % stringValueForKey(TWOD_OFFSET_DIR_KEY,1))
         case("x")
            dir2D = 1
         case("y")
            dir2D = 2
         case("z")
            dir2D = 3
         case("3d","3D")
            dir2D = 0
         case default
            print*, "Unrecognized 2D mesh offset direction"
            error stop
            errorMessage(STD_OUT)
         end select

      else
         dir2D = 0

      end if

      if (controlVariables % containsKey("mesh inner curves")) then
         MeshInnerCurves = controlVariables % logicalValueForKey("mesh inner curves")
      else
         MeshInnerCurves = .true.
      end if


      useRelaxPeriodic = controlVariables % logicalValueForKey("periodic relative tolerance")
      useWeightsPartition = controlVariables % getValueOrDefault("partitioning with weights", .true.)
!
!     **********************************************************
!     *                  MPI PREPROCESSING                     *
!     **********************************************************
!
!     Initialization
!     --------------
      if (.not. self % mesh % child) then
         call Initialize_MPI_Partitions ( trim(controlVariables % stringValueForKey('partitioning', requestedLength = LINE_LENGTH)) )
!
!        Prepare the processes to receive the partitions
!        -----------------------------------------------
         if ( MPI_Process % doMPIAction ) then
            call RecvPartitionMPI( MeshFileType(self % mesh % meshFileName) == HOPRMESH )
         end if
!
!        Read the mesh by the root rank to perform the partitioning
!        ----------------------------------------------------------
         if ( MPI_Process % doMPIRootAction ) then
!
!           Construct the "full" mesh
!           -------------------------
            call constructMeshFromFile( self % mesh, self % mesh % meshFileName, CurrentNodes, Nx, Ny, Nz, MeshInnerCurves , dir2D, useRelaxPeriodic, success )
!
!           Perform the partitioning
!           ------------------------
            call PerformMeshPartitioning  (self % mesh, MPI_Process % nProcs, mpi_allPartitions, useWeightsPartition)
!
!           Send the partitions
!           -------------------
            call SendPartitionsMPI( MeshFileType(self % mesh % meshFileName) == HOPRMESH )
!
!           Destruct the full mesh
!           ----------------------
            call self % mesh % Destruct()

         end if

      end if
!
!     **********************************************************
!     *              MESH CONSTRUCTION                         *
!     **********************************************************
!
      if (MPI_Process % isRoot) write(STD_OUT,'(/,5X,A)') "Reading mesh..."
      CALL constructMeshFromFile( self % mesh, self % mesh % meshFileName, CurrentNodes, Nx, Ny, Nz, MeshInnerCurves , dir2D, useRelaxPeriodic, success )
      if (.not. self % mesh % child) call mpi_partition % ConstructGeneralInfo (self % mesh % no_of_allElements)   

!
!     **********************************************************
!     *              SLIDING MESH CONSTRUCTION            *
!     **********************************************************

!     
!     Sliding Mesh parameter
!     -----------------------------------
      call self% mesh% SlidingMesh% read_info( controlVariables )

      if (self% mesh% SlidingMesh% isConfigured) then 

         call AdvanceSlidingMesh(self% mesh, self% mesh% SlidingMesh% radius, self% mesh% SlidingMesh%center &
         , self% mesh% numBFacePoints, self% mesh% nodeType, .FALSE.)

      end if 

!
!     
!     Immersed boundary method parameter
!     -----------------------------------

      call self% mesh% IBM% read_info( controlVariables )

!     Compute wall distances
!     ----------------------
#if defined(NAVIERSTOKES)
      call self % mesh % ComputeWallDistances
#endif
      IF(.NOT. success) RETURN
!
!     construct surfaces mesh
!     -----------------------
      call surfacesMesh % construct(controlVariables, self % mesh)
!
!     ----------------------------
!     Get the final number of DOFS
!     ----------------------------
!
      self % NDOF = 0
      DO k=1, self % mesh % no_of_elements
         associate(e => self % mesh % elements(k))
         self % NDOF = self % NDOF + (e % Nxyz(1) + 1) * (e % Nxyz(2) + 1) * (e % Nxyz(3) + 1)
         end associate
      END DO
!
!     ----------------------------------
!     Get the final total number of DOFs
!     ----------------------------------
!
      if ( MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
         call mpi_allreduce(self % NDOF, self % totalNDOF, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      else
         self % totalNDOF = self % NDOF
      end if

!
!     **********************************************************
!     *              IMMERSED BOUNDARY CONSTRUCTION            *
!     **********************************************************
!
      if( self% mesh% IBM% active ) then
         if( .not. self % mesh % child ) then
            call self% mesh% IBM% GetDomainExtreme( self% mesh% elements )
            call self% mesh% IBM% construct( controlVariables )
         end if
!
!        ------------------------------------------------
!        building the IBM mask and the IBM band region
!        ------------------------------------------------
!
         call self% mesh% IBM% build( self% mesh% elements, self% mesh% no_of_elements, self% mesh% NDOF, self% mesh% child )
      end if

!
!
!     ------------------------
!     Allocate and zero memory
!     ------------------------
!
      call self % mesh % AllocateStorage(self % NDOF, controlVariables,computeGradients)
!
!     --------------------
!     Initialize Base Flow
!     --------------------
!
#if defined(ACOUSTIC)
      ! start by default with no flow conditions
      QbaseUniform = [1.0_RP,0.0_RP,0.0_RP,0.0_RP,1.0_RP/(dimensionless % gammaM2)]
      call self % mesh % SetUniformBaseFlow(QbaseUniform)
      call self % mesh % ProlongBaseSolutionToFaces(NCONS)
#ifdef _HAS_MPI_
!$omp single
      call self % mesh % UpdateMPIFacesBaseSolution(NCONS)
      ! not efficient, but only done once
      ! we can wait for the communication with more computation in between, but will need to be in a different subroutine
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call self % mesh % GatherMPIFacesBaseSolution(NCONS)
!$omp end single
#endif
!
#endif
!
!     ----------------------------------------------------
!     Get manufactured solution source term (if requested)
!     ----------------------------------------------------
!
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
      IF (self % ManufacturedSol) THEN
      IF (flowIsNavierStokes) THEN
         DO el = 1, SIZE(self % mesh % elements)
            DO k=0, Nz(el)
               DO j=0, Ny(el)
                  DO i=0, Nx(el)                     
                        CALL ManufacturedSolutionSourceNS(self % mesh % elements(el) % geom % x(:,i,j,k), &
                                                          0._RP, &
                                                          self % mesh % elements(el) % storage % S_NS (:,i,j,k)  )
                  END DO
               END DO
            END DO
         END DO
      ELSE
         DO el = 1, SIZE(self % mesh % elements)
            DO k=0, Nz(el)
               DO j=0, Ny(el)
                  DO i=0, Nx(el)
                         CALL ManufacturedSolutionSourceEuler(self % mesh % elements(el) % geom % x(:,i,j,k), &
                                                             0._RP, &
                                                          self % mesh % elements(el) % storage % S_NS (:,i,j,k)  )
                  END DO
               END DO
            END DO
         END DO
      END IF
      END IF

#elif defined(SPALARTALMARAS)
      IF (self % ManufacturedSol) THEN
         DO el = 1, SIZE(self % mesh % elements)
           DO k=0, Nz(el)
               DO j=0, Ny(el)
                  DO i=0, Nx(el)
                        CALL ManufacturedSolutionSourceNSSA(self % mesh % elements(el) % geom % x(:,i,j,k), &
                                                            self % mesh % elements(el) % geom % dwall(i,j,k), 0._RP, &
                                                            self % mesh % elements(el) % storage % S_NS (:,i,j,k)  )
                  END DO
               END DO
            END DO
         END DO
      END IF
#endif
!
!     ------------------
!     Build the monitors
!     ------------------
!
      call self % monitors % construct (self % mesh, controlVariables)
!
!     ------------------
!     Build the FWH general class
!     ------------------
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
      IF (flowIsNavierStokes) call self % fwh % construct(self % mesh, controlVariables)
#endif

! #if defined(NAVIERSTOKES)
! !
! !     -------------------
! !     Build the particles
! !     -------------------
! !

!       self % particles % active = controlVariables % logicalValueForKey("lagrangian particles")
!       if ( self % particles % active ) then
!             call self % particles % construct(self % mesh, controlVariables)
!       endif
! #endif

      NULLIFY(Nx,Ny,Nz)
!
!     Stop measuring preprocessing time
!     ----------------------------------
      if (.not. self % mesh % child) call Stopwatch % Pause("Preprocessing")

      END SUBROUTINE ConstructDGSem

!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructDGSem( self )
      use SurfaceMesh, only: surfacesMesh
      IMPLICIT NONE
      CLASS(DGSem) :: self
      INTEGER      :: k      !Counter

      CALL self % mesh % destruct

      call self % monitors % destruct

#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
      IF (flowIsNavierStokes) call self % fwh % destruct
#endif

      call surfacesMesh % destruct

      END SUBROUTINE DestructDGSem
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SaveSolutionForRestart( self, fUnit )
         IMPLICIT NONE
         CLASS(DGSem)     :: self
         INTEGER          :: fUnit
         INTEGER          :: k

         DO k = 1, SIZE(self % mesh % elements)
            WRITE(fUnit) self % mesh % elements(k) % storage % Q
         END DO

      END SUBROUTINE SaveSolutionForRestart

      subroutine DGSEM_SetInitialCondition( self, controlVariables, initial_iteration, initial_time )
         use FTValueDictionaryClass
         USE mainKeywordsModule
         use SurfaceMesh, only: surfacesMesh
         implicit none
         class(DGSEM)   :: self
         class(FTValueDictionary), intent(in)   :: controlVariables
         integer                                :: restartUnit
         integer,       intent(out)             :: initial_iteration
         real(kind=RP), intent(out)             :: initial_time
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)             :: solutionName
         logical                                :: saveGradients, loadFromNSSA, withSensor, saveLES
         procedure(UserDefinedInitialCondition_f) :: UserDefinedInitialCondition

         solutionName = controlVariables % stringValueForKey(solutionFileNameKey, requestedLength = LINE_LENGTH)
         solutionName = trim(getFileName(solutionName))

         withSensor = controlVariables % logicalValueForKey(saveSensorToSolutionKey)

         IF ( controlVariables % logicalValueForKey(restartKey) )     THEN
            loadFromNSSA = controlVariables % logicalValueForKey("ns from nssa")
            if (loadFromNSSA .and. MPI_Process % isRoot) write(STD_OUT,'(/,5X,A)') "NS restarting from RANS"
            CALL self % mesh % LoadSolutionForRestart(controlVariables, initial_iteration, initial_time, loadFromNSSA)
         ELSE

            call UserDefinedInitialCondition(self % mesh, FLUID_DATA_VARS)

            initial_time = 0.0_RP
            initial_iteration = 0
!
!           Save the initial condition
!           --------------------------
!
            saveGradients = controlVariables % logicalValueForKey(saveGradientsToSolutionKey)
            saveLES = controlVariables % logicalValueForKey(saveLESToSolutionKey)
            IF(controlVariables % stringValueForKey(solutionFileNameKey,LINE_LENGTH) /= "none")     THEN           
               write(solutionName,'(A,A,I10.10,A)') trim(solutionName), "_", initial_iteration, ".hsol"
               call self % mesh % SaveSolution(initial_iteration, initial_time, solutionName, saveGradients, withSensor, saveLES)
               !TDG: ADD PARTICLES WRITE WITH IFDEF
            END IF 
         END IF

         IF(controlVariables % stringValueForKey(solutionFileNameKey,LINE_LENGTH) /= "none")     THEN
            write(solutionName,'(A,A,I10.10)') trim(solutionName), "_", initial_iteration
            call self % mesh % Export( trim(solutionName) )

            call surfacesMesh % saveAllMesh(self % mesh, initial_iteration, controlVariables)
         END IF 

      end subroutine DGSEM_SetInitialCondition
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------------------------
!     Specifies how to copy a DGSem object:
!     -> Trivial but must be defined since self % monitors has a user-defined assignment
!     -> Must be impure because the assignment procedures of derived times are not pure
!        (they perform better being impure!!)
!     -----------------------------------------------------------------------------------
      subroutine DGSEM_Assign (to, from)
         implicit none
         class(DGSem), intent(inout) :: to
         type (DGSem), intent(in)    :: from

         to % maxResidual        = from % maxResidual
         to % nodes              = from % nodes
         to % numberOfTimeSteps  = from % numberOfTimeSteps
         to % NDOF               = from % NDOF
         to % totalNDOF          = from % totalNDOF

         to % mesh               = from % mesh
         to % ManufacturedSol    = from % ManufacturedSol

         to % monitors  = from % monitors
         to % particles = from % particles

      end subroutine DGSEM_Assign
!
!////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------
!  Compute maximum residual L_inf norm
!  -----------------------------------
!
   FUNCTION ComputeMaxResiduals(mesh) RESULT(maxResidual)
      use MPI_Process_Info
      IMPLICIT NONE
      CLASS(HexMesh), intent(in)  :: mesh
      REAL(KIND=RP), dimension(NCONS) :: maxResidual
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER       :: id , eq, ierr
      REAL(KIND=RP) :: localMaxResidual(NCONS)
      real(kind=RP) :: localR1, localR2, localR3, localR4, localR5, localR6, localc
      real(kind=RP) :: R1, R2, R3, R4, R5, R6, c

      maxResidual = 0.0_RP
      R1 = 0.0_RP
      R2 = 0.0_RP
      R3 = 0.0_RP
      R4 = 0.0_RP
      R5 = 0.0_RP
      R6 = 0.0_RP
      c    = 0.0_RP

!$omp parallel shared(maxResidual, R1, R2, R3, R4, R5, R6, c, mesh) default(private)
!$omp do reduction(max:R1,R2,R3,R4,R5, R6, c) schedule(runtime)
      DO id = 1, SIZE( mesh % elements )
#if defined FLOW && !(SPALARTALMARAS)
         localR1 = maxval(abs(mesh % elements(id) % storage % QDot(1,:,:,:)))
         localR2 = maxval(abs(mesh % elements(id) % storage % QDot(2,:,:,:)))
         localR3 = maxval(abs(mesh % elements(id) % storage % QDot(3,:,:,:)))
         localR4 = maxval(abs(mesh % elements(id) % storage % QDot(4,:,:,:)))
         localR5 = maxval(abs(mesh % elements(id) % storage % QDot(5,:,:,:)))
#else
         localR1 = maxval(abs(mesh % elements(id) % storage % QDot(1,:,:,:)))
         localR2 = maxval(abs(mesh % elements(id) % storage % QDot(2,:,:,:)))
         localR3 = maxval(abs(mesh % elements(id) % storage % QDot(3,:,:,:)))
         localR4 = maxval(abs(mesh % elements(id) % storage % QDot(4,:,:,:)))
         localR5 = maxval(abs(mesh % elements(id) % storage % QDot(5,:,:,:)))
         localR6 = maxval(abs(mesh % elements(id) % storage % QDot(6,:,:,:)))
#endif

#ifdef CAHNHILLIARD
         localc    = maxval(abs(mesh % elements(id) % storage % cDot(:,:,:,:)))
#endif

#if defined FLOW && !(SPALARTALMARAS)
         R1 = max(R1,localR1)
         R2 = max(R2,localR2)
         R3 = max(R3,localR3)
         R4 = max(R4,localR4)
         R5 = max(R5,localR5)
#elif defined(SPALARTALMARAS)
         R1 = max(R1,localR1)
         R2 = max(R2,localR2)
         R3 = max(R3,localR3)
         R4 = max(R4,localR4)
         R5 = max(R5,localR5)
         R6 = max(R6,localR6)
#endif

#ifdef CAHNHILLIARD
         c    = max(c, localc)
#endif
      END DO
!$omp end do
!$omp end parallel

#if defined FLOW && (!(SPALARTALMARAS))
      maxResidual(1:NCONS) = [R1, R2, R3, R4, R5]
#elif defined(SPALARTALMARAS)
      maxResidual(1:NCONS) = [R1, R2, R3, R4, R5, R6]
#endif

#if  defined(CAHNHILLIARD) && (!defined(FLOW))
      maxResidual(NCONS) = c
#endif

#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
         localMaxResidual = maxResidual
         call mpi_allreduce(localMaxResidual, maxResidual, NCONS, MPI_DOUBLE, MPI_MAX, &
                            MPI_COMM_WORLD, ierr)
      end if
#endif

   END FUNCTION ComputeMaxResiduals
!
!////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------
!  Estimate the maximum time-step of the system. This
!  routine computes a heuristic based on the smallest mesh
!  spacing (which goes as 1/N^2) AND the eigenvalues of the
!  particular hyperbolic system being solved (Euler or Navier Stokes). This is not
!  the only way to estimate the time-step, but it works in practice.
!  Other people use variations on this and we make no assertions that
!  it is the only or best way. Other variations look only at the smallest
!  mesh values, other account for differences across the element.
!     -------------------------------------------------------------------
   subroutine MaxTimeStep( self, cfl, dcfl, MaxDt , MaxDtVec)
      use VariableConversion
      use MPI_Process_Info
      implicit none
      !------------------------------------------------
      type(DGSem)                :: self
      real(kind=RP), intent(in)  :: cfl      !<  Advective cfl number
      real(kind=RP), optional, intent(in)  :: dcfl     !<  Diffusive cfl number
      real(kind=RP), intent(inout)  :: MaxDt
      real(kind=RP), allocatable, dimension(:), intent(inout), optional :: MaxDtVec
#ifdef FLOW
      !------------------------------------------------
      integer                       :: i, j, k, eID                     ! Coordinate and element counters
      integer                       :: N(3)                             ! Polynomial order in the three reference directions
      integer                       :: ierr                             ! Error for MPI calls
      real(kind=RP)                 :: eValues(3)                       ! Advective eigenvalues
      real(kind=RP)                 :: dcsi, deta, dzet                 ! Smallest mesh spacing
      real(kind=RP)                 :: dcsi2, deta2, dzet2              ! Smallest mesh spacing squared
      real(kind=RP)                 :: lamcsi_a, lamzet_a, lameta_a     ! Advective eigenvalues in the three reference directions
      real(kind=RP)                 :: lamcsi_v, lamzet_v, lameta_v     ! Diffusive eigenvalues in the three reference directions
      real(kind=RP)                 :: jac, mu, T                       ! Mapping Jacobian, viscosity and temperature
      real(kind=RP)                 :: kinematicviscocity, musa, etasa
      real(kind=RP)                 :: Q(NCONS)                           ! The solution in a node
      real(kind=RP)                 :: TimeStep_Conv, TimeStep_Visc     ! Time-step for convective and diffusive terms
      real(kind=RP)                 :: localMax_dt_v, localMax_dt_a     ! Time step to perform MPI reduction
      type(NodalStorage_t), pointer :: spAxi_p, spAeta_p, spAzeta_p     ! Pointers to the nodal storage in every direction
      external                      :: ComputeEigenvaluesForState       ! Advective eigenvalues
#if defined(INCNS) || defined(MULTIPHASE)
      logical :: flowIsNavierStokes = .true.
#endif
#if defined(SPALARTALMARAS)
      external                            :: ComputeEigenvaluesForStateSA
#elif defined(ACOUSTIC)
      external                            :: ComputeEigenvaluesForStateCAA
#endif
      !--------------------------------------------------------
!     Initializations
!     ---------------
      TimeStep_Conv = huge(1._RP)
      TimeStep_Visc = huge(1._RP)
      if (present(MaxDtVec)) MaxDtVec = huge(1._RP)
!$omp parallel shared(self,TimeStep_Conv,TimeStep_Visc,NodalStorage,cfl,dcfl,flowIsNavierStokes,MaxDtVec) default(private)
!$omp do reduction(min:TimeStep_Conv,TimeStep_Visc) schedule(runtime)
      do eID = 1, SIZE(self % mesh % elements)
         N = self % mesh % elements(eID) % Nxyz
         spAxi_p => NodalStorage(N(1))
         spAeta_p => NodalStorage(N(2))
         spAzeta_p => NodalStorage(N(3))

         if ( N(1) .ne. 0 ) then
            dcsi = 1.0_RP / abs( spAxi_p   % x(1) - spAxi_p   % x (0) )

         else
            dcsi = 0.0_RP

         end if

         if ( N(2) .ne. 0 ) then
            deta = 1.0_RP / abs( spAeta_p  % x(1) - spAeta_p  % x (0) )

         else
            deta = 0.0_RP

         end if

         if ( N(3) .ne. 0 ) then
            dzet = 1.0_RP / abs( spAzeta_p % x(1) - spAzeta_p % x (0) )

         else
            dzet = 0.0_RP

         end if

         if (flowIsNavierStokes) then
            dcsi2 = dcsi*dcsi
            deta2 = deta*deta
            dzet2 = dzet*dzet
         end if

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
!
!           ------------------------------------------------------------
!           The maximum eigenvalues for a particular state is determined
!           by the physics.
!           ------------------------------------------------------------
!
            Q(1:NCONS) = self % mesh % elements(eID) % storage % Q(1:NCONS,i,j,k)

#if defined(SPALARTALMARAS)
            CALL ComputeEigenvaluesForStateSA( Q , eValues )
#elif defined(ACOUSTIC)
            CALL ComputeEigenvaluesForStateCAA( Q , self % mesh % elements(eID) % storage % Qbase(:,i,j,k), eValues )
#else
            CALL ComputeEigenvaluesForState( Q , eValues )
#endif
            jac      = self % mesh % elements(eID) % geom % jacobian(i,j,k)
!
!           ----------------------------
!           Compute contravariant values
!           ----------------------------
!
            lamcsi_a =abs( self % mesh % elements(eID) % geom % jGradXi(IX,i,j,k)   * eValues(IX) + &
                           self % mesh % elements(eID) % geom % jGradXi(IY,i,j,k)   * eValues(IY) + &
                           self % mesh % elements(eID) % geom % jGradXi(IZ,i,j,k)   * eValues(IZ) ) * dcsi

            lameta_a =abs( self % mesh % elements(eID) % geom % jGradEta(IX,i,j,k)  * eValues(IX) + &
                           self % mesh % elements(eID) % geom % jGradEta(IY,i,j,k)  * eValues(IY) + &
                           self % mesh % elements(eID) % geom % jGradEta(IZ,i,j,k)  * eValues(IZ) ) * deta

            lamzet_a =abs( self % mesh % elements(eID) % geom % jGradZeta(IX,i,j,k) * eValues(IX) + &
                           self % mesh % elements(eID) % geom % jGradZeta(IY,i,j,k) * eValues(IY) + &
                           self % mesh % elements(eID) % geom % jGradZeta(IZ,i,j,k) * eValues(IZ) ) * dzet

            TimeStep_Conv = min( TimeStep_Conv, cfl*abs(jac)/(lamcsi_a+lameta_a+lamzet_a) )
            if (present(MaxDtVec)) MaxDtVec(eID) = min( MaxDtVec(eID), cfl*abs(jac)/(lamcsi_a+lameta_a+lamzet_a) )

#if defined(NAVIERSTOKES)
            if (flowIsNavierStokes) then
               T        = Temperature(Q)
               mu       = SutherlandsLaw(T)

#if defined(SPALARTALMARAS)

              call GetNSKinematicViscosity(mu, self % mesh % elements(eID) % storage % Q(IRHO,i,j,k), kinematicviscocity )
              call self % SAModel % ComputeViscosity( self % mesh % elements(eID) % storage % Q(IRHOTHETA,i,j,k), kinematicviscocity, &
                                               self % mesh % elements(eID) % storage % Q(IRHO,i,j,k), mu, musa, etasa)
              mu = mu + musa
#endif
               lamcsi_v = mu * dcsi2 * abs(sum(self % mesh % elements(eID) % geom % jGradXi  (:,i,j,k)))
               lameta_v = mu * deta2 * abs(sum(self % mesh % elements(eID) % geom % jGradEta (:,i,j,k)))
               lamzet_v = mu * dzet2 * abs(sum(self % mesh % elements(eID) % geom % jGradZeta(:,i,j,k)))
               TimeStep_Visc = min( TimeStep_Visc, dcfl*abs(jac)/(lamcsi_v+lameta_v+lamzet_v) )
               if (present(MaxDtVec)) MaxDtVec(eID) = min( MaxDtVec(eID), &
                                                      dcfl*abs(jac)/(lamcsi_v+lameta_v+lamzet_v)  )
            end if
#else
            TimeStep_Visc = huge(1.0_RP)
#endif

         end do ; end do ; end do
      end do
!$omp end do
!$omp end parallel

#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
         localMax_dt_v = TimeStep_Visc
         localMax_dt_a = TimeStep_Conv
         call mpi_allreduce(localMax_dt_a, TimeStep_Conv, 1, MPI_DOUBLE, MPI_MIN, &
                            MPI_COMM_WORLD, ierr)
         call mpi_allreduce(localMax_dt_v, TimeStep_Visc, 1, MPI_DOUBLE, MPI_MIN, &
                            MPI_COMM_WORLD, ierr)
      end if
#endif

      if (TimeStep_Conv  < TimeStep_Visc) then
         self % mesh % dt_restriction = DT_CONV
         MaxDt  = TimeStep_Conv
      else
         self % mesh % dt_restriction = DT_DIFF
         MaxDt  = TimeStep_Visc
      end if
#endif
   end subroutine MaxTimeStep
!
!////////////////////////////////////////////////////////////////////////
!
   subroutine hnRange(mesh, hnmin, hnmax)
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(HexMesh), intent(in)  :: mesh
      real(RP),      intent(out) :: hnmin
      real(RP),      intent(out) :: hnmax
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: eID, ierr
      real(RP) :: hn


      hnmin = huge(1.0_RP)
      hnmax = -huge(1.0_RP)
      do eID = 1, mesh % no_of_elements
         hn = mesh % elements(eID) % hn
         hnmin = min(hn, hnmin)
         hnmax = max(hn, hnmax)
      end do

      call MPI_MinMax(hnmin, hnmax)

   end subroutine hnRange
!
end module DGSEMClass
