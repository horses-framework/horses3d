!
!//////////////////////////////////////////////////////
!
!      FAS Multigrid Class
!        Provides the routines for solving a time step with nonlinear multigrid procedures.
!        Available smoothers are:
!           -> RK3:         Explicit smoother based on the Williamson's low-storage 3rd order Runge-Kutta. Only valid for STEADY_STATE cases since the time-stepping is advanced in every level.
!           -> BlockJacobi: Implicit smoother based on the matrix-free version of the classical Block-Jacobi method. STEADY_STATE or TIME_ACCURATE.
!           -> GMRES:       Implicit smoother based on the matrix-free GMRES method. STEADY_STATE or TIME_ACCURATE.
!
!        **NOTE: The implicit smoothers are currently using BDF time-integration. This can be changed...
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "Includes.h"
module FASMultigridClass
   use SMConstants
   use ExplicitMethods
   use DGSEMClass
   use PhysicsStorage
   use Physics
   use InterpolationMatrices
   use MultigridTypes
   use TimeIntegratorDefinitions
   use BDFTimeIntegrator
   use FileReadingUtilities      , only: getFileName
   use MPI_Process_Info          , only: MPI_Process
   use FileReadingUtilities      , only: getIntArrayFromString
   use MatrixClass
   use CSRMatrixClass         , only: csrMat_t
   use AnalyticalJacobian     , only: AnJacobian_t
   use NumericalJacobian      , only: NumJacobian_t
   use JacobianComputerClass  , only: JacobianComputer_t, GetJacobianFlag
   use DenseMatUtilities
   use MPI_Utilities          , only: infNorm, L2Norm! , MPI_SumAll
   use mkl_spblas
   use IBMClass
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS)) && (!(TRANSPORT))
   use ManufacturedSolutionsNS
#elif defined(SPALARTALMARAS)
   use ManufacturedSolutionsNSSA
#endif

   implicit none

   private
   public FASMultigrid_t

   type :: LUpivots_t
   !-----Variables-----------------------------------------------------------
         integer      , dimension(:)  , allocatable :: v   ! LU pivots
   end type LUpivots_t

!
!  Multigrid class
!  ---------------
   type :: FASMultigrid_t
      type(DGSem)              , pointer      :: p_sem                 ! Pointer to DGSem class variable of current system
      type(FASMultigrid_t)     , pointer      :: Child                 ! Next coarser multigrid solver
      type(FASMultigrid_t)     , pointer      :: Parent                ! Next finer multigrid solver
      type(MGSolStorage_t)     , allocatable  :: MGStorage(:)          ! Storage
      integer                                 :: MGlevel               ! Current Multigrid level
      real(kind=RP),             allocatable  :: lts_dt(:)             ! dt array for LTS

      ! variables for implicit time integration
      class(JacobianComputer_t), allocatable       :: Jacobian           ! Jacobian
      class(Matrix_t), allocatable                 :: A                  ! Jacobian matrix
      class(LUpivots_t), allocatable, dimension(:) :: LUpivots
      real(kind=RP), dimension(:)  , allocatable   :: dQ
      real(kind=RP), dimension(:)  , allocatable   :: dQ0
      real(kind=RP), dimension(:)  , allocatable   :: SGS_RHS
      integer                                      :: JacobianComputation = NUMERICAL_JACOBIAN
      logical                                      :: computeA              !< Compute A in this level?
      integer                                      :: DimPrb                ! problem size
      integer                                      :: GlobalDimPrb                ! global problem size
#ifdef HAS_MKL
      type(matrix_descr)                           :: descrA
      type(sparse_matrix_t)                        :: csrA
#endif

      contains
         procedure :: construct
         procedure :: solve
         procedure :: Smooth
         procedure :: destruct
         procedure :: SetPreviousSolution => FAS_SetPreviousSolution    ! For implicit smoothing, it's necessary to store the previous solution(s) in all levels
         procedure :: TakePseudoStep ! solve for Dual Time Stepping 

   end type FASMultigrid_t
!
!  ----------------
!  Module variables
!  ----------------
!
   procedure(SmoothIt_t), pointer :: SmoothIt

   !! Parameters
   integer, parameter :: MAX_SWEEPS_DEFAULT = 10000

   ! Other variables
   integer        :: MaxN           ! Maximum polynomial order of the mesh
   integer        :: NMIN           ! Minimum polynomial order allowed
   integer        :: MGlevels       ! Total number of multigrid levels
   integer        :: deltaN         !
   integer        :: nelem          ! Number of elements
   integer        :: num_of_allElems
   integer        :: Smoother       ! Current smoother being used
   integer        :: MaxSweeps      ! Maximum number of sweeps in a smoothing process
   logical        :: MGOutput       ! Display output?
   logical        :: FMG = .FALSE.  ! Use Full Multigrid algorithm?
   logical        :: PostFCycle,PostSmooth ! Post smoothing options
   logical        :: SmoothFine     !
   logical        :: ManSol         ! Does this case have manufactured solutions?
   logical        :: Compute_dt
   logical        :: SaveFMGFile =.FALSE.
   character(len=LINE_LENGTH)              :: FMGSolutionFile
   logical                                 :: saveGradients
   logical                                 :: saveSensor
   real(kind=RP)  :: SmoothFineFrac ! Fraction that must be smoothed in fine before going to coarser level
   real(kind=RP), target  :: cfl            ! Advective cfl number
   real(kind=RP), target  :: dcfl           ! Diffusive cfl number
   real(kind=RP), target  :: dt             ! dt
   integer, allocatable :: MGSweepsPre(:) ! Number of pre- and post-smoothings operations on each level
   integer, allocatable :: MGSweepsPost(:) ! Number of post- and post-smoothings operations on each level
   integer        :: Preconditioner       ! Current smoother being used
   integer        :: CurrentMGCycle
!-----CFL-ramping-variables-----------------------------------------------------------
   character(len=LINE_LENGTH) :: CFLboost  = "none"
   character(len=LINE_LENGTH) :: DCFLboost = "none"
   real(kind=RP)  :: cfl_max                    ! Max. advective cfl number (for CFL boost)
   real(kind=RP)  :: dcfl_max                   ! Max. diffusive cfl number (for CFL boost)
   real(kind=RP)  :: cflboost_rate              ! Max. diffusive cfl number (for CFL boost)
   integer        :: erk_order = 5              ! Steady-state OptERK type
!-----DTS-variables-------------------------------------------------------------------
   logical        :: DualTimeStepping = .false.
   logical        :: Compute_Global_dt = .true.
   logical        :: PseudoConvergenceMonitor = .false.
   real(kind=RP)  :: conv_tolerance = 1e-6_RP
   real(kind=RP), target  :: p_cfl            ! Pseudo advective cfl number
   real(kind=RP), target  :: p_dcfl           ! Pseudo diffusive cfl number
   real(kind=RP), target  :: p_dt             ! Pseudo dt
!-----Implicit-relaxation-variables---------------------------------------------------
   integer        :: MatrixType = JACOBIAN_MATRIX_NONE
   integer        :: StepsForJac = 1e8
   integer        :: StepsSinceJac
   integer        :: kSGS
!-----Initilization-------------------------------------------------------------------
   integer              :: ini_Preconditioner(2)
   integer              :: ini_Smoother(2)
   real(kind=RP)        :: ini_res
   real(kind=RP)        :: ini_cfl(4)
   logical              :: mg_Initialization = .false.
   logical              :: mg_Initialization_Present = .false.
!========
 contains
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine construct(this,controlVariables,sem)
      use FTValueDictionaryClass
      use StopwatchClass
      implicit none
      !-----------------------------------------------------------
      class(FASMultigrid_t)  , intent(inout), target :: this               !<> Anisotropic FAS multigrid solver to be constructed
      type(FTValueDictionary), intent(in)            :: controlVariables   !<  Input variables
      type(DGSem)            , intent(in)   , target :: sem                !<  Fine sem class
      !-----------------------------------------------------------
      character(len=LINE_LENGTH)                     :: PostSmoothOptions
      integer                                        :: zoneID                ! Zone counter
      logical                                        :: conformingBoundaries  ! Is the mesh conforming on all boundaries?
      character(len=LINE_LENGTH)                     :: tmpc
      integer                                        :: i
      !-----------------------------------------------------------

      call Stopwatch % Pause("Solver")
      call Stopwatch % Start("Preprocessing")

!
!     ----------------------------------
!     Read important variables from file
!     ----------------------------------
!
!
!     Read cfl and dcfl numbers
!     -------------------------
      if (controlVariables % containsKey("cfl")) then
#if defined(NAVIERSTOKES)
         Compute_dt = .TRUE.
         cfl = controlVariables % doublePrecisionValueForKey("cfl")
         ! if there is no dcfl, use the same number as for cfl
         if (flowIsNavierStokes) then
            if (controlVariables % containsKey("dcfl")) then
               dcfl       = controlVariables % doublePrecisionValueForKey("dcfl")
            else
               dcfl = cfl
            end if
         end if
         ! CFL boosting
         if (controlVariables % containsKey("cfl boost")) then
            CFLboost = controlVariables % StringValueForKey("cfl boost",LINE_LENGTH)
         end if
         ! DCFL boosting
         if (controlVariables % containsKey("dcfl boost")) then
            DCFLboost = controlVariables % StringValueForKey("dcfl boost",LINE_LENGTH)
         end if
         ! boost rate
         if (controlVariables % containsKey("cfl boost rate")) then
            CFLboost_rate = controlVariables % doublePrecisionValueForKey("cfl boost rate")
         else
            CFLboost_rate = 0.1_RP
         end if

         ! max cfl
         if (controlVariables % containsKey("cfl max")) then
            cfl_max = controlVariables % doublePrecisionValueForKey("cfl max")
         else
            cfl_max = 1.0_RP
         end if
         if (controlVariables % containsKey("dcfl max")) then
            dcfl_max = controlVariables % doublePrecisionValueForKey("dcfl max")
         else
            dcfl_max = 1.0_RP
         end if

#elif defined(CAHNHILLIARD)
         print*, "Error, use fixed time step to solve Cahn-Hilliard equations"
         errorMessage(STD_OUT)
         error stop
#endif
      elseif (controlVariables % containsKey("dt")) then
         Compute_dt = .false.
         Compute_Global_dt = .false.
         dt = controlVariables % doublePrecisionValueForKey("dt")
      else
         error stop '"cfl" (and "dcfl" if Navier-Stokes) or "dt" keywords must be specified for the FAS integrator'
      end if

!
!     Pseudo time stepping variables
!     -------------------------
      if ( trim(controlVariables % StringValueForKey("simulation type",LINE_LENGTH)) == "time-accurate" ) then
         ! is the simulation type is time-accurate and we use FAS
         DualTimeStepping = .true.
         Compute_dt = .true.

         if (controlVariables % containsKey("pseudo convergence monitor")) then
            PseudoConvergenceMonitor = controlVariables % logicalValueForKey("pseudo convergence monitor")
         end if

         if (controlVariables % containsKey("convergence tolerance")) &
            conv_tolerance = controlVariables % doublePrecisionValueForKey("convergence tolerance")

            if (controlVariables % containsKey("bdf order")) call BDF_SetOrder( controlVariables % integerValueForKey("bdf order") )
            call BDFInitialiseQ(sem % mesh)

         if (controlVariables % containsKey("pseudo cfl")) then
#if defined(NAVIERSTOKES)
            p_cfl = controlVariables % doublePrecisionValueForKey("pseudo cfl")
            if (flowIsNavierStokes) then
               if (controlVariables % containsKey("pseudo dcfl")) then
                  p_dcfl       = controlVariables % doublePrecisionValueForKey("pseudo dcfl")
               else
                  p_dcfl = p_cfl
               end if
            end if
#elif defined(CAHNHILLIARD)
            print*, "Error, use fixed time step to solve Cahn-Hilliard equations"
            errorMessage(STD_OUT)
            error stop
#endif
         elseif (controlVariables % containsKey("pseudo dt")) then
            p_dt = controlVariables % doublePrecisionValueForKey("pseudo dt")
            Compute_dt = .false.
         else
            error stop '"pseudo cfl" or "pseudo dt" keywords must be specified for the time-accurate FAS integrator'
         end if
      end if ! time-accurate

!
!     Read multigrid variables
!     -------------------------
      if (.NOT. controlVariables % containsKey("multigrid levels")) then
         print*, 'Fatal error: "multigrid levels" keyword is needed by the FASMultigrid solver'
         error stop
      end if

      MGlevels  = controlVariables % IntegerValueForKey("multigrid levels")

      if (controlVariables % containsKey("delta n")) then
         deltaN = controlVariables % IntegerValueForKey("delta n")
      else
         deltaN = 1
      end if

!
!     Number of sweeps
!     -------------------

      if (.not. allocated(MGSweepsPre)) allocate(MGSweepsPre(MGlevels))
      if (.not. allocated(MGSweepsPost)) allocate(MGSweepsPost(MGlevels))
      MGsweepsPre  = 1
      MGsweepsPost = 1

      if (controlVariables % containsKey("mg sweeps pre" ) .AND. &
          controlVariables % containsKey("mg sweeps post") ) then
        MGsweepsPre = controlVariables % IntegerValueForKey("mg sweeps pre")
        MGsweepsPost =controlVariables % IntegerValueForKey("mg sweeps post")
      else if (controlVariables % containsKey("mg sweeps" ) ) then
        MGsweepsPre = controlVariables % IntegerValueForKey("mg sweeps")
        MGsweepsPost = controlVariables % IntegerValueForKey("mg sweeps")
      end if
      if (controlVariables % containsKey("mg sweeps coarsest")) then
        MGsweepsPre(1) = controlVariables % IntegerValueForKey("mg sweeps coarsest")
        MGsweepsPost(1) =controlVariables % IntegerValueForKey("mg sweeps coarsest")
      end if

      if ( controlVariables % containsKey("mg sweeps pre exact") .and. controlVariables % containsKey("mg sweeps post exact") ) then
         tmpc = controlVariables % StringValueForKey("mg sweeps pre exact",LINE_LENGTH)
         MGSweepsPre = getIntArrayFromString(tmpc)
         tmpc = controlVariables % StringValueForKey("mg sweeps post exact",LINE_LENGTH)
         MGSweepsPost = getIntArrayFromString(tmpc)
      else if (controlVariables % containsKey("mg sweeps exact") ) then
         tmpc = controlVariables % StringValueForKey("mg sweeps exact",LINE_LENGTH)
         MGSweepsPre = getIntArrayFromString(tmpc)
         MGSweepsPost = getIntArrayFromString(tmpc)
      end if
!
!     Select the smoother
!     -------------------
      select case (controlVariables % StringValueForKey("mg smoother",LINE_LENGTH))
         case('Euler')
            Smoother = Euler_SMOOTHER
         case('RK3')
            Smoother = RK3_SMOOTHER
         case('RK5')
            Smoother = RK5_SMOOTHER
         case('RKOpt')
            Smoother = RKOpt_SMOOTHER
            ! Select order of ERK scheme
            if ( controlVariables % containsKey("rk order") ) then
               erk_order = controlVariables % IntegerValueForKey("rk order")
               if (erk_order .gt. 7) then
                  erk_order = 7
                  print *, "FASMultigrid :: ERK Order too high, switching to 7."
               else if (erk_order .lt. 2) then
                  erk_order = 2
                  print *, "FASMultigrid :: ERK Order too low, switching to 2."
               end if
            end if
         case('IRK')
            Smoother = IRK_SMOOTHER
            MatrixType = JACOBIAN_MATRIX_CSR
         case('SGS')
            Smoother = SGS_SMOOTHER
            MatrixType = JACOBIAN_MATRIX_CSR

            if (controlVariables % containsKey("k gauss seidel")) then
               kSGS = controlVariables % IntegerValueForKey("k gauss seidel")
            else
               kSGS = 1
            end if

         case('ILU')
            Smoother = ILU_SMOOTHER
            MatrixType = JACOBIAN_MATRIX_CSR
         case('BIRK5')
            Smoother = BIRK5_SMOOTHER
            MatrixType = JACOBIAN_MATRIX_DENSE
         case default
            if (MPI_Process % isRoot) write(STD_OUT,*) '"mg smoother" not recognized. Defaulting to RK3.'
            Smoother = RK3_SMOOTHER
      end select

!
!     Additional options for implicit smoothers
!     -------------------
      if (Smoother .ge. IMPLICIT_SMOOTHER_IDX) then
         if (controlVariables % containsKey("compute jacobian every")) StepsForJac = controlVariables % integerValueForKey("compute jacobian every")
      end if

!
!     Select the preconditioner for smoothing
!     ---------------------------------------
      if (controlVariables % containsKey("mg preconditioner")) then
         select case (controlVariables % StringValueForKey("mg preconditioner",LINE_LENGTH))
         case('LTS')
            Preconditioner = PRECONDIIONER_LTS
         case default
            Preconditioner = PRECONDIIONER_NONE
         end select
      else
         Preconditioner = PRECONDIIONER_NONE
      end if

!
!     More control parameters for mg cycle
!     -------------------------------
      PostSmoothOptions = controlVariables % StringValueForKey("postsmooth option",LINE_LENGTH)
      if (trim(PostSmoothOptions) == 'f-cycle') then
         PostFCycle = .true.
      elseif (trim(PostSmoothOptions) == 'smooth') then
         PostSmooth = .true.
      end if

      if (controlVariables % containsKey("smooth fine")) then
         SmoothFine = .TRUE.
         SmoothFineFrac = controlVariables % doublePrecisionValueForKey("smooth fine")
      else
         SmoothFine = .FALSE.
      end if

      if (controlVariables % containsKey("max mg sweeps")) then
         MaxSweeps = controlVariables % IntegerValueForKey("max mg sweeps")
      else
         MaxSweeps = MAX_SWEEPS_DEFAULT
      end if

      if (controlVariables % logicalValueForKey("fasfmg save solutions") ) then
         SaveFMGFile = .TRUE.
         saveGradients = controlVariables % logicalValueForKey("save gradients with solution")
         saveSensor    = controlVariables % logicalValueForKey("save sensor with solution")
         FMGSolutionFile = trim(getFileName(controlVariables % stringValueForKey("solution file name", requestedLength = LINE_LENGTH)))
      end if

!
!     Read variables for the initial solution cycle
!     -------------------------
      if (controlVariables % containsKey("mg initialization")) then
         mg_Initialization = controlVariables % logicalValueForKey("mg initialization")
         if (mg_Initialization) then
            mg_Initialization_Present = .true.

            if (controlVariables % containsKey("initial residual")) then
               ini_res = controlVariables % doublePrecisionValueForKey("initial residual")
            else
               ini_res = 1.0d0
            end if

            if (controlVariables % containsKey("initial preconditioner")) then
               select case (controlVariables % StringValueForKey("initial preconditioner",LINE_LENGTH))
               case('LTS')
                  ini_Preconditioner(1) = PRECONDIIONER_LTS
               case default
                  ini_Preconditioner(1) = PRECONDIIONER_NONE
               end select
            else
               ini_Preconditioner(1) = PRECONDIIONER_LTS
            end if
            ini_Preconditioner(2) = Preconditioner ! desired preconditioner
            ini_Smoother(1) = RK5_SMOOTHER ! smoother for initialization
            ini_Smoother(2) = Smoother ! desired smoother

            ! cfl/dcfl for initialization
            if (controlVariables % containsKey("initial cfl")) then
               ini_cfl(1) = controlVariables % doublePrecisionValueForKey("initial cfl")
               ini_cfl(2) = controlVariables % doublePrecisionValueForKey("initial cfl")
            else
               ini_cfl(1) = 0.5d0
               ini_cfl(2) = 0.5d0
            end if
            ! desired cfl/dcfl
            if (controlVariables % containsKey("cfl")) then
               ini_cfl(3) = cfl
               ini_cfl(4) = dcfl
            else if (controlVariables % containsKey("pseudo cfl")) then
               ini_cfl(3) = p_cfl
               ini_cfl(4) = p_dcfl
            else
               error stop "FASMultigridClass :: "
            end if

         end if
      end if
!
!     ------------------------------------------
!     Get the minimum multigrid polynomial order
!     ------------------------------------------
!
      NMIN = 1
!
!     -----------------------
!     Update module variables
!     -----------------------
!
      MGOutput       = controlVariables % logicalValueForKey("multigrid output")
      plotInterval   = controlVariables % integerValueForKey("output interval")
      ManSol         = sem % ManufacturedSol
      MaxN           = MAX(MAXVAL(sem % mesh % Nx),MAXVAL(sem % mesh % Ny),MAXVAL(sem % mesh % Nz))
      MGlevels       = MIN (MGlevels,MaxN - NMIN + 1)

      if (MPI_Process % isRoot) then
         write(STD_OUT,*) 'Constructing FAS Multigrid'
         write(STD_OUT,*) 'Number of levels:', MGlevels
      end if

      this % p_sem => sem

      nelem = SIZE(sem % mesh % elements)
      num_of_allElems = sem % mesh % no_of_allElements

!
!     --------------------------
!     Create linked solvers list
!     --------------------------
!
      call RecursiveConstructor(this, sem % mesh % Nx, sem % mesh % Ny, sem % mesh % Nz, MGlevels, controlVariables)

      call Stopwatch % Pause("Preprocessing")
      call Stopwatch % Start("Solver")

   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine RecursiveConstructor(Solver, N1x, N1y, N1z, lvl, controlVariables)
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS)) && (!(TRANSPORT))
   use ManufacturedSolutionsNS
#elif defined(SPALARTALMARAS)
   use ManufacturedSolutionsNSSA
#endif
      use FTValueDictionaryClass
      implicit none
      type(FASMultigrid_t), target  :: Solver
      integer, dimension(:)         :: N1x,N1y,N1z      !<  Order of approximation for every element in current solver
      integer                       :: lvl              !<  Current multigrid level
      type(FTValueDictionary)       :: controlVariables !< Control variables (for the construction of coarse sems
      !----------------------------------------------
      integer, dimension(nelem)           :: N2x,N2y,N2z            !   Order of approximation for every element in child solver
      integer, dimension(num_of_allElems) :: N2xAll,N2yAll,N2zAll   !   Order of approximation for every element in child solver
      integer                             :: i,j,k, iEl             !   Counter
      logical                             :: success                ! Did the creation of sem succeed?
      type(FASMultigrid_t), pointer       :: Child_p                ! Pointer to Child
      integer                             :: Q1,Q2,Q3,Q4            ! Sizes of vector Q (conserved solution) used for allocation. In this version the argument MOLD of ALLOCATE is not used since several versions of gfortran don't support it yet...
      !----------------------------------------------
      !
      integer :: Nxyz(3), fd, l
      integer :: N1(3), N2(3)
      integer, dimension(:), allocatable :: nnz_perblock

      Solver % MGlevel = lvl
      Solver % DimPrb = Solver % p_sem % NDOF * NCONS
      Solver % globalDimPrb = Solver % p_sem % totalNDOF * NCONS
!
!     --------------------------
!     Allocate Multigrid storage
!     --------------------------
!
      allocate (Solver % MGStorage(nelem))
!$omp parallel do private(Q1,Q2,Q3,Q4) schedule(runtime)
      DO k = 1, nelem
         Q1 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,1)
         Q2 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,2) - 1
         Q3 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,3) - 1
         Q4 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,4) - 1
         allocate(Solver % MGStorage(k) % Q    (Q1,0:Q2,0:Q3,0:Q4))
         allocate(Solver % MGStorage(k) % E    (Q1,0:Q2,0:Q3,0:Q4))
         allocate(Solver % MGStorage(k) % S    (Q1,0:Q2,0:Q3,0:Q4))
         allocate(Solver % MGStorage(k) % Scase(Q1,0:Q2,0:Q3,0:Q4))

         if (DualTimeStepping) then
            allocate(Solver % MGStorage(k) % R    (Q1,0:Q2,0:Q3,0:Q4))
            allocate(Solver % MGStorage(k) % Q0   (Q1,0:Q2,0:Q3,0:Q4))
         end if

         Solver % MGStorage(k) % Scase = 0._RP
      end DO
!$omp end parallel do

      ! allocate storage for implicit relaxation
      if (Smoother .ge. IMPLICIT_SMOOTHER_IDX) then
         if (.not. (allocated(Solver % dQ)) ) allocate( Solver % dQ(Solver % DimPrb) )
         Solver % dQ = 0._RP

         select case (Smoother)
         case (ILU_SMOOTHER)
            if (.not. (allocated(Solver % dQ0)) ) allocate( Solver % dQ0(Solver % DimPrb) )
            Solver % dQ0 = 0._RP
         case (SGS_SMOOTHER)
            if (.not. (allocated(Solver % dQ0)) ) allocate( Solver % dQ0(Solver % DimPrb) )
            Solver % dQ0 = 0._RP
            if (.not. (allocated(Solver % SGS_RHS)) ) allocate( Solver % SGS_RHS(Solver % DimPrb) )
            Solver % SGS_RHS = 0._RP
         end select

      end if

      ! allocate array for LTS
      if (Preconditioner .eq. PRECONDIIONER_LTS .or. mg_Initialization) allocate( Solver % lts_dt(nelem))
!
!     --------------------------------------------------------------
!     Fill MGStorage(iEl) % Scase if required (manufactured solutions)
!        (only for lower meshes)
!     --------------------------------------------------------------
!
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS)) && (!(TRANSPORT))
      if (ManSol) then
         DO iEl = 1, nelem

            DO k=0, Solver % p_sem % mesh % Nz(iEl)
               DO j=0, Solver % p_sem % mesh % Ny(iEl)
                  DO i=0, Solver % p_sem % mesh % Nx(iEl)
                     if (flowIsNavierStokes) then
                        call ManufacturedSolutionSourceNS(Solver % p_sem % mesh % elements(iEl) % geom % x(:,i,j,k), &
                                                          0._RP, &
                                                          Solver % MGStorage(iEl) % Scase (:,i,j,k)  )
                     else
                        call ManufacturedSolutionSourceEuler(Solver % p_sem % mesh % elements(iEl) % geom % x(:,i,j,k), &
                                                             0._RP, &
                                                             Solver % MGStorage(iEl) % Scase (:,i,j,k)  )
                     end if
                  end DO
               end DO
            end DO
         end DO
      end if
#elif defined(SPALARTALMARAS)
      if (ManSol) then
         DO iEl = 1, nelem
            DO k=0, Solver % p_sem % mesh % Nz(iEl)
               DO j=0, Solver % p_sem % mesh % Ny(iEl)
                  DO i=0, Solver % p_sem % mesh % Nx(iEl)
                        CALL ManufacturedSolutionSourceNSSA(Solver % p_sem % mesh % elements(iEl) % geom % x(:,i,j,k), &
                                                            Solver % p_sem % mesh % elements(iEl) % geom % dwall(i,j,k), 0._RP, &
                                                            Solver % MGStorage(iEl) % Scase (:,i,j,k)  )
                  END DO
               END DO
            END DO
         END DO
      END IF
#endif
!
!     -------------------------------------------
!     Assemble Jacobian for implicit smoothing
!     -------------------------------------------

      if ( Smoother .ge. IMPLICIT_SMOOTHER_IDX ) then

         Solver % JacobianComputation = GetJacobianFlag()
         select case (Solver % JacobianComputation)
            case (NOTDEF_JACOBIAN )    ; allocate(Solver % Jacobian)
            case (NUMERICAL_JACOBIAN ) ; allocate(NumJacobian_t :: Solver % Jacobian)
            case (ANALYTICAL_JACOBIAN) ; allocate(AnJacobian_t  :: Solver % Jacobian)
            case default
               error stop 'Invalid jacobian type'
         end select
         call Solver % Jacobian % construct(Solver % p_sem % mesh, NCONS, controlVariables)
!
!        Construct Jacobian
!        ---------------------
            select case ( MatrixType )
               case (JACOBIAN_MATRIX_NONE)
               case (JACOBIAN_MATRIX_DENSE)
               !
               ! Construct blocks only
               ! ---------------------
                  allocate(DenseBlockDiagMatrix_t :: Solver % A)
                  call Solver % A % construct (num_of_Blocks = Solver % p_sem % mesh % no_of_elements)

                  allocate(nnz_perblock(nelem))

                  do k=1,nelem
                     nnz_perblock(k) = NCONS*(Solver % p_sem % mesh % elements(k) % Nxyz(1)+1)*&
                           (Solver % p_sem % mesh % elements(k) % Nxyz(2)+1)*&
                           (Solver % p_sem % mesh % elements(k) % Nxyz(3)+1)
                  end do
                  call Solver % A % PreAllocate (nnzs=nnz_perblock)
                  call Solver % Jacobian % Configure (Solver % p_sem % mesh, NCONS, Solver % A)

                  ! allocate vectors for pivots
                  allocate (Solver % LUpivots(nelem))
                  do k = 1,nelem
                     allocate ( Solver % LUpivots(k) % v( nnz_perblock(k)) )
                  end do

                  deallocate(nnz_perblock)
               case (JACOBIAN_MATRIX_CSR)
               !
               ! Construct full Jacobian matrix
               ! ---------------------
                  allocate(csrMat_t :: Solver % A)
                  call Solver % A % Construct(num_of_Rows = Solver % DimPrb, num_of_TotalRows = Solver % globalDimPrb)
                  call Solver % Jacobian % Configure (Solver % p_sem % mesh, NCONS, Solver % A)
               case default

            end select
         Solver % computeA = .TRUE.
      end if

      if (lvl > 1) then
         ALLOCATE  (Solver % Child)
         Child_p => Solver % Child
         Solver % Child % Parent => Solver

!
!        ---------------------------------------------
!        Create restriction and prolongation operators
!        ---------------------------------------------
!
         DO k=1, nelem
            call CreateInterpolationOperators(N1x(k), N2x(k), MaxN, MGlevels, lvl-1, DeltaN, Solver % p_sem % nodes)
            call CreateInterpolationOperators(N1y(k), N2y(k), MaxN, MGlevels, lvl-1, DeltaN, Solver % p_sem % nodes)
            call CreateInterpolationOperators(N1z(k), N2z(k), MaxN, MGlevels, lvl-1, DeltaN, Solver % p_sem % nodes)
         end DO

         ! Create DGSEM class for child
         ALLOCATE (Child_p % p_sem)

         !<old
         ! Get the polynomial orders for constructing the child sem
         N2xAll = 0
         N2yAll = 0
         N2zAll = 0
         do k=1, nelem
            N2xAll( Solver % p_sem % mesh % elements(k) % globID ) = N2x(k)
            N2yAll( Solver % p_sem % mesh % elements(k) % globID ) = N2y(k)
            N2zAll( Solver % p_sem % mesh % elements(k) % globID ) = N2z(k)
         end do
         
         ! setting the IBM level & saving the KDtree

         if( Solver% p_sem% mesh% IBM% active ) call Child_p% p_sem% mesh% IBM% copy( Solver% p_sem% mesh% IBM, lvl )

         call Child_p % p_sem % construct (controlVariables = controlVariables,                                          &
                                           Nx_ = N2xAll,    Ny_ = N2yAll,    Nz_ = N2zAll,                               &
                                           success = success,                                                            &
                                           ChildSem = .TRUE.  )
 
         if (.NOT. success) error stop "Multigrid: Problem creating coarse solver."

         if (DualTimeStepping) then
!$omp do private(N1,N2) schedule(runtime)
            DO k = 1, nelem
               N1 = Solver  % p_sem % mesh % elements (k) % Nxyz
               N2 = Child_p % p_sem % mesh % elements (k) % Nxyz
               call Interp3DArrays(NCONS, N1, Solver  % p_sem % mesh % elements(k) % storage % Q, &
                                          N2, Child_p % p_sem % mesh % elements(k) % storage % Q )
            end DO
!$omp end do
            call BDFInitialiseQ(Child_p % p_sem % mesh)
         end if

         call RecursiveConstructor(Solver % Child, N2x, N2y, N2z, lvl - 1, controlVariables)
      end if

   end subroutine RecursiveConstructor
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Driver of the FAS multigrid solving procedure
!  ---------------------------------------------
   subroutine solve(this, timestep, t, ComputeTimeDerivative, ComputeTimeDerivativeIsolated, FullMG, tol)
      implicit none
      !-------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this
      integer                              :: timestep
      real(kind=RP)        , intent(in)    :: t
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivativeIsolated
      logical           , OPTIONAL         :: FullMG
      real(kind=RP)     , OPTIONAL         :: tol        !<  Tolerance for full multigrid
      !-------------------------------------------------
      integer :: maxVcycles = 40, i
      real(kind=RP) :: rnorm, xnorm
      integer :: firstIdx, lastIdx, eID
      real(kind=RP), pointer :: fassolve_dt, fassolve_cfl, fassolve_dcfl

      ThisTimeStep = timestep

      if (PRESENT(FullMG)) then
         if (FullMG) then
            if (.NOT. PRESENT(tol)) error stop 'FASFMG needs tolerance'
            FMG = .TRUE.
         else
            FMG = .FALSE.
         end if
      else
         FMG = .FALSE.
      end if

      if (.not. DualTimeStepping) then
         if (this % computeA) then
            StepsSinceJac = 0
         else
            StepsSinceJac = StepsSinceJac + 1
            if (StepsSinceJac .eq. StepsForJac) then
               call computeA_AllLevels(this,MGlevels)
               StepsSinceJac = 0
            end if
         end if
      end if
!
!     -----------------------
!     Perform multigrid cycle
!     -----------------------
!
      if (DualTimeStepping) then
         fassolve_dt  => p_dt
         fassolve_cfl  => p_cfl
         fassolve_dcfl => p_dcfl
      else
         fassolve_dt  => dt
         fassolve_cfl  => cfl
         fassolve_dcfl => dcfl
      end if

      if (FMG) then
         call FASFMGCycle(this,t,tol,MGlevels, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
      else
         call FASVCycle(this,t,MGlevels,MGlevels, ComputetimeDerivative, ComputeTimeDerivativeIsolated)
      end if

      call CFLRamp(cfl_max,fassolve_cfl,cflboost_rate,CFLboost)
      call CFLRamp(dcfl_max,fassolve_dcfl,cflboost_rate,DCFLboost)

   end subroutine solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Driver of the Dual time-stepping procedure.
!  Q_{m+1} = Q_m + d\tau <(> (Q_m - Q_n)/dt + R(Q_m) <)> (BDF1 example)
!  ---------------------------------------------
   subroutine TakePseudoStep(this, timestep, t, ComputeTimeDerivative, ComputeTimeDerivativeIsolated, FullMG, tol)
      implicit none
      !-------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this
      integer                              :: timestep
      real(kind=RP)        , intent(in)    :: t
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivativeIsolated
      logical           , OPTIONAL         :: FullMG
      real(kind=RP)     , OPTIONAL         :: tol        !<  Tolerance for full multigrid
      !-------------------------------------------------
      integer :: i, id
      integer :: tau_maxit = 10000
      real(kind=RP) :: dQdtau_norm, Qdot_norm
      real(kind=RP) :: tk
!
!     -----------------------
!     Solve local steady-state problem
!     -----------------------
!
      if (this % computeA) then
         StepsSinceJac = 0
      else
         StepsSinceJac = StepsSinceJac + 1
         if (StepsSinceJac .eq. StepsForJac) then
            call computeA_AllLevels(this,MGlevels)
            StepsSinceJac = 0
         end if
      end if

      tk = t
      if (Compute_Global_dt) call MaxTimeStep( self=this % p_sem, cfl=cfl, dcfl=dcfl , MaxDt=dt)

      call ComputeTimeDerivative( this % p_sem % mesh, this % p_sem % particles, tk, CTD_IGNORE_MODE)
      Qdot_norm = MAXVAL(ComputeMaxResiduals(this % p_sem % mesh))
      call ComputePseudoTimeDerivative(this % p_sem % mesh, tk, dt)
      dQdtau_norm = MAXVAL(ComputeMaxResiduals(this % p_sem % mesh))

      ! set previous solution
      do i= 1, bdf_order
         if (i .eq. bdf_order) then
!$omp parallel do schedule(runtime)
            do id = 1, SIZE(this % p_sem % mesh % elements )
               this % p_sem % mesh % elements(id) % storage % prevQ(i) % Q = this % p_sem % mesh % elements(id) % storage % Q
            end do ! id
!$omp end parallel do
         else
!$omp parallel do schedule(runtime)
            do id = 1, SIZE(this % p_sem % mesh % elements )
               this % p_sem % mesh % elements(id) % storage % prevQ(i) % Q = this % p_sem % mesh % elements(id) % storage % prevQ(i+1) % Q
            end do ! id
!$omp end parallel do
         end if
      end do

      if( this% p_sem% mesh% IBM% active ) then
          if( any(this% p_sem% mesh% IBM%  stl(:)% move) .and. MGlevels-1 > 1 ) call FAS_movingIBM( this% child, dt, MGlevels-1 )
      end if

      do i = 1, tau_maxit

         call this % solve(i, tk, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
         dQdtau_norm = MAXVAL(ComputeMaxResiduals(this % p_sem % mesh))
         if (PseudoConvergenceMonitor) then
            if (MPI_Process % isRoot ) write(STD_OUT,'(30X,A,I4,A,ES10.3)') "Pseudo Iter= ", i, ", Res= ", dQdtau_norm
         end if
         if (dQdtau_norm .le. conv_tolerance) exit

      end do

      dQdtau_norm = MAXVAL(ComputeMaxResiduals(this % p_sem % mesh))
      if (MPI_Process % isRoot ) write(STD_OUT,'(20X,A,I4,A,ES10.3)') "--- Pseudo time step converged in ", i, " iterations to Res= ", dQdtau_norm

      tk = tk + dt

!$omp parallel do schedule(runtime)
      do id = 1, SIZE(this % p_sem % mesh % elements )
         this % p_sem % mesh % elements(id) % storage % prevQ(1) % Q = this % p_sem % mesh % elements(id) % storage % Q
      end do ! id
!$omp end parallel do

      call ComputeTimeDerivative( this % p_sem % mesh, this % p_sem % particles, tk, CTD_IGNORE_MODE)

      if (mg_Initialization_Present) mg_Initialization = .true. ! set back explicit initialization for the next step

   end subroutine TakePseudoStep
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------
!  Signals computeA=.TRUE. in all MG levels
!  -----------------------------------------
   recursive subroutine computeA_AllLevels(this,lvl)
      implicit none
      class(FASMultigrid_t), intent(inout) :: this     !<  Current level solver
      integer :: lvl

      this % computeA = .TRUE.

      if (lvl>1) call computeA_AllLevels(this % Child,lvl-1)

   end subroutine computeA_AllLevels
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------
!  Recursive subroutine to perform a v-cycle
!  -----------------------------------------
   recursive subroutine FASVCycle(this,t,lvl,MGlevels, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
      implicit none
      !----------------------------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this     !<  Current level solver
      real(kind=RP)        , intent(in)    :: t        !<  Simulation time
      integer              , intent(in)    :: lvl      !<  Current multigrid level
      integer              , intent(in)    :: MGlevels !<  Number of finest multigrid level
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivativeIsolated
      !----------------------------------------------------------------------------
      integer                       :: iEl,iEQ              !Element/equation counter
      type(FASMultigrid_t), pointer :: Child_p              !Pointer to child
      integer                       :: N1(3), N2(3)
      real(kind=RP)                 :: maxResidual(NCONS)
      real(kind=RP)                 :: PrevRes
      real(kind=RP)                 :: NewRes
      integer                       :: sweepcount           ! Number of sweeps done in a point in time
      real(kind=RP)                 :: invdt
      integer                       :: k, stat
      real(kind=RP), pointer :: fasvcycle_dt, fasvcycle_cfl, fasvcycle_dcfl
      !----------------------------------------------------------------------------
#if defined(NAVIERSTOKES)
!
!     -----------------------------------------
!     ============ Update FAS info ============
!     -----------------------------------------
!
!     Update explicit initialization
!     ------------------------------
      call FAS_UpdateInitialization(this,lvl)

!     Check if we solve local or global problem
!     -----------------------------------------
      if (DualTimeStepping) then
         fasvcycle_dt  => p_dt
         fasvcycle_cfl  => p_cfl
         fasvcycle_dcfl => p_dcfl
      else
         fasvcycle_dt  => dt
         fasvcycle_cfl  => cfl
         fasvcycle_dcfl => dcfl
      end if

!     Compute dt and 1/dt
!     -------------------
      if (Compute_dt) call MaxTimeStep(self=this % p_sem, cfl=fasvcycle_cfl, dcfl=fasvcycle_dcfl, MaxDt=fasvcycle_dt )
      invdt = -1._RP/fasvcycle_dt

!     Taking care of Jacobian (for implicit residual relaxation)
!     ----------------------------------------------------------
      call FAS_ComputeAndFactorizeJacobian(this, lvl, t, invdt, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)

!
!     -------------------------------------------------------
!     ============ Recursive V-Cycle starts here ============
!     -------------------------------------------------------
!
!     Pre-smoothing procedure
!     -----------------------
      sweepcount = 0
      do
         call this % Smooth(MGSweepsPre(lvl),t, ComputeTimeDerivative)
         sweepcount = sweepcount + MGSweepsPre(lvl)

         if (MGOutput) call PlotResiduals( lvl , sweepcount,this % p_sem % mesh)

         if (SmoothFine .AND. lvl > 1) then ! .AND. .not. FMG
            if (FMG .and. MAXVAL(ComputeMaxResiduals(this % p_sem % mesh)) < 0.1_RP) exit
            call MGRestrictToChild(this,lvl-1,t, ComputeTimeDerivative)
            call ComputeTimeDerivative(this % Child % p_sem % mesh,this % Child % p_sem % particles, t, CTD_IGNORE_MODE)
            if (DualTimeStepping) call ComputePseudoTimeDerivative(this % Child % p_sem % mesh, t, dt)

            if (MAXVAL(ComputeMaxResiduals(this % p_sem % mesh)) < SmoothFineFrac * MAXVAL(ComputeMaxResiduals(this % Child % p_sem % mesh))) exit
         else
            exit
         end if

         if (sweepcount .ge. MaxSweeps) exit
      end do

      PrevRes = MAXVAL(ComputeMaxResiduals(this % p_sem % mesh))


      if (lvl > 1) then
         if (.not. SmoothFine) call MGRestrictToChild(this,lvl-1,t, ComputeTimeDerivative)
!
!        --------------------
!        Perform V-Cycle here
!        --------------------
!
         call FASVCycle(this % Child, t, lvl-1,MGlevels, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)

         Child_p => this % Child
!
!        -------------------------------------------
!        Interpolate coarse-grid error to this level
!        -------------------------------------------
!
!$omp parallel
!$omp do private(N1,N2) schedule(runtime)
         DO iEl = 1, nelem
            N1 = Child_p % p_sem % mesh % elements (iEl) % Nxyz
            N2 = this    % p_sem % mesh % elements (iEl) % Nxyz
            call Interp3DArrays(NCONS, N1, Child_p % MGStorage(iEl) % E, N2, this % MGStorage(iEl) % E )
         end DO
!$omp end do
!
!        -----------------------------------------------
!        Correct solution with coarse-grid approximation
!        -----------------------------------------------
!
!$omp barrier
!$omp do schedule(runtime)
         DO iEl = 1, nelem
            this % p_sem % mesh % elements(iEl) % storage % Q = &
                              this % p_sem % mesh % elements(iEl) % storage % Q + this % MGStorage(iEl) % E
         end DO
!$omp end do
!$omp end parallel

      end if

!
!     ------------------------
!     Post-smoothing procedure
!     ------------------------
!

      sweepcount = 0
      DO
         call this % Smooth(MGSweepsPost(lvl), t, ComputeTimeDerivative)
         sweepcount = sweepcount + MGSweepsPost(lvl)

         NewRes = MAXVAL(ComputeMaxResiduals(this % p_sem % mesh))

         if (MGOutput) call PlotResiduals( lvl, sweepcount , this % p_sem % mesh)

         if (sweepcount .ge. MaxSweeps) exit

         if (lvl > 1 .and. PostFCycle) then
            if (NewRes > PrevRes) then
               call MGRestrictToChild(this,lvl-1,t, ComputeTimeDerivative)
               call FASVCycle(this,t,lvl-1,lvl, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
            else
               exit
            end if
         elseif (PostSmooth .or. PostFCycle) then
            !if (FMG .and. MAXVAL(ComputeMaxResiduals(this % p_sem % mesh)) < 0.1_RP) exit
            if (NewRes < PrevRes) exit
         else
            exit
         end if

      end DO

!
!     -------------------------
!     Compute coarse-grid error
!     -------------------------
!
      if (lvl < MGlevels) then
!$omp parallel do schedule(runtime)
         DO iEl = 1, nelem
            this % MGStorage(iEl) % E = this % p_sem % mesh % elements(iEl) % storage % Q - this % MGStorage(iEl) % Q
         end DO
!$omp end parallel do
      end if

#endif
   end subroutine FASVCycle
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------
!  Recursive subroutine to perform a full multigrid cycle
!  ------------------------------------------------------
   recursive subroutine FASFMGCycle(this,t,tol,lvl, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
      implicit none
      !----------------------------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this    !<> Current level solver
      real(kind=RP)        , intent(in)    :: t       !<  Simulation time
      real(kind=RP)        , intent(in)    :: tol     !<  Convergence tolerance
      integer              , intent(in)    :: lvl     !<  Current multigrid level
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivativeIsolated
      !----------------------------------------------------------------------------
      integer        :: iEl, iEQ             ! Element and equation counters
      integer        :: N1(3), N2(3)
      real(kind=RP)  :: maxResidual(NCONS)   ! Maximum residual in each equation
      integer        :: counter              ! Iteration counter
      character(len=LINE_LENGTH) :: FMGFile
      !----------------------------------------------------------------------------
#if defined(NAVIERSTOKES)
!
!     ------------------------------------------
!     At the beginning, go to the coarsest level
!        (the initial condition must be passed)
!     ------------------------------------------
!
      if (lvl > 1) then
!$omp parallel do private(N1,N2) schedule(runtime)
         DO iEl = 1, nelem
            N1 = this         % p_sem % mesh % elements (iEl) % Nxyz
            N2 = this % Child % p_sem % mesh % elements (iEl) % Nxyz
            call Interp3DArrays(NCONS, N1, this % p_sem % mesh % elements(iEl) % storage % Q, &
                                       N2, this % Child % p_sem % mesh % elements(iEl) % storage % Q )
         end DO
!$omp end parallel do

         call FASFMGCycle(this % Child,t,tol,lvl-1, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
      end if
!
!     ------------------------------
!     Save FMG solution if requested
!     ------------------------------
!

      if (SaveFMGFile) then
         write(FMGFile,'(A,A,I2.2,A)')  trim( FMGSolutionFile ), '_FMG_', lvl, '.hsol'
         call this % p_sem % mesh % SaveSolution(0,0._RP,trim(FMGFile),saveGradients,saveSensor)
      end if

!
!     ----------------------
!     Perform a V-Cycle here
!     ----------------------
!
      counter = 0
      if (lvl > 1 ) then
         DO
            counter = counter + 1
            call FASVCycle(this,t,lvl,lvl, ComputeTimeDerivative, ComputeTimeDerivativeIsolated) ! FMG is still for STEADY_STATE. TODO: Change that
            maxResidual = ComputeMaxResiduals(this % p_sem % mesh)
            if (maxval(maxResidual) <= tol) exit
         end DO
      else
         DO
            counter = counter + 1
            call this % Smooth(1,t,ComputeTimeDerivative)

            maxResidual = ComputeMaxResiduals(this % p_sem % mesh)

            if (MOD(counter,100)==0) call PlotResiduals( lvl ,counter, this % p_sem % mesh)
            if (maxval(maxResidual) <= tol) exit
         end DO
      end if
      call PlotResiduals( lvl ,counter, this % p_sem % mesh,.TRUE.)
!
!     --------------------------------------------------
!     If not on finest, Interpolate to next (finer) grid
!     --------------------------------------------------
!
      if (lvl < MGlevels) then
!$omp parallel do private(N1,N2) schedule(runtime)
         DO iEl = 1, nelem
            N1 = this          % p_sem % mesh % elements (iEl) % Nxyz
            N2 = this % Parent % p_sem % mesh % elements (iEl) % Nxyz
            call Interp3DArrays(NCONS, N1, this % p_sem % mesh % elements(iEl) % storage % Q, &
                                       N2, this % Parent % p_sem % mesh % elements(iEl) % storage % Q )
         end DO
!$omp end parallel do
      end if
#endif
   end subroutine FASFMGCycle

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------
!  Subroutine for restricting the solution and residual to the child solver
!  ------------------------------------------
   subroutine MGRestrictToChild(this,lvl,t, ComputeTimeDerivative)
      implicit none
      !-------------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this     !<  Current level solver
      integer              , intent(IN)    :: lvl
      real(kind=RP)        , intent(IN)    :: t
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      !-------------------------------------------------------------
      class(FASMultigrid_t), pointer       :: Child_p  ! The child
      integer  :: iEl
      integer  :: iEQ
      integer  :: N1(3), N2(3)
      !-------------------------------------------------------------
#if defined(NAVIERSTOKES)

      Child_p => this % Child

!$omp parallel
!$omp do private(N1,N2) schedule(runtime)
      DO iEl = 1, nelem
         N1 = this    % p_sem % mesh % elements (iEl) % Nxyz
         N2 = Child_p % p_sem % mesh % elements (iEl) % Nxyz

!           Restrict solution
!           -----------------
         call Interp3DArrays(NCONS, N1, this % p_sem % mesh % elements(iEl) % storage % Q, &
                                    N2, Child_p % p_sem % mesh % elements(iEl) % storage % Q )

!           Restrict residual
!           -----------------
         call Interp3DArrays(NCONS, N1, this % p_sem % mesh % elements(iEl) % storage % Qdot, &
                                    N2, Child_p % MGStorage(iEl) % S)
      end DO
!$omp end do

!
!     **********************************************************************
!     **********************************************************************
!              Now arrange all the storage in the child solver
!     **********************************************************************
!     **********************************************************************
!

!     ------------------------------------
!     Copy solution from fine grid to MGStorage
!        ... and clear source term
!     ------------------------------------
!
!$omp barrier
!$omp do schedule(runtime)
      DO iEl = 1, nelem
         Child_p % MGStorage(iEl) % Q = Child_p % p_sem % mesh % elements(iEl) % storage % Q
         Child_p % p_sem % mesh % elements(iEl) % storage % S_NS =  0.0_RP
      end DO
!$omp end do
!$omp end parallel

      if( Child_p% p_sem% mesh% IBM% TimePenal ) Child_p% p_sem% mesh% IBM% penalization = dt

!
!     -------------------------------------------
!     If not on finest level, correct source term
!     -------------------------------------------
!
      call ComputeTimeDerivative(Child_p % p_sem % mesh,Child_p % p_sem % particles, t, CTD_IGNORE_MODE)
      if (DualTimeStepping) call ComputePseudoTimeDerivative(Child_p % p_sem % mesh, t, dt)

!$omp parallel do schedule(runtime)
      DO iEl = 1, nelem
         Child_p % p_sem % mesh % elements(iEl) % storage % S_NS = Child_p % MGStorage(iEl) % S - &
                                                                Child_p % p_sem % mesh % elements(iEl) % storage % Qdot !- Child_p % MGStorage(iEl) % Scase
      end DO
!$omp end parallel do

#endif
   end subroutine MGRestrictToChild
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------
!  Subroutine that destructs a FAS integrator
!  ------------------------------------------
   subroutine destruct(this)
      implicit none
      !-----------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this
      !-----------------------------------------------------------

      call RecursiveDestructor(this,MGlevels)

   end subroutine destruct

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine RecursiveDestructor(Solver,lvl)
      implicit none
      !-----------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: Solver
      integer                              :: lvl
      integer                              :: k
      !-----------------------------------------------------------

      ! First go to coarsest level
      if (lvl > 1) call RecursiveDestructor(Solver % Child,lvl-1)

      ! Deallocate LTS
      if (allocated(Solver % lts_dt)) deallocate(Solver % lts_dt)

      !Destruct Multigrid storage
      deallocate (Solver % MGStorage) ! allocatable components are automatically deallocated

      ! Destruct linear solver (if present)
      if (Smoother .ge. IMPLICIT_SMOOTHER_IDX) then
         ! call Solver % linsolver % destroy FINDME
         ! deallocate ( Solver % linsolver )
      end if

      if (lvl < MGlevels) then
         call Solver % p_sem % destruct()
         deallocate (Solver % p_sem)
         nullify    (Solver % Parent)
      else
         nullify    (Solver % p_sem)
      end if

      if (lvl > 1) then
         deallocate (Solver % Child)
      end if

      if ( allocated(Solver % LUpivots) ) then
         do k = 1,nelem
            deallocate ( Solver % LUpivots(k) % v )
         end do
         deallocate ( Solver % LUpivots )
      end if

   end subroutine RecursiveDestructor
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Smooth(this,SmoothSweeps,t, ComputeTimeDerivative)
      implicit none
      !-------------------------------------------------------------
      class(FASMultigrid_t)  , intent(inout), target :: this     !<> Anisotropic FAS multigrid
      integer                , intent(in)            :: SmoothSweeps
      real(kind=RP)          , intent(in)            :: t
      procedure(ComputeTimeDerivative_f)                     :: ComputeTimeDerivative
      !-------------------------------------------------------------
      real(kind=RP), pointer :: smoother_dt, smoother_cfl, smoother_dcfl
      integer :: sweep
      !-------------------------------------------------------------

      if (DualTimeStepping) then
         smoother_dt   => p_dt
         smoother_cfl  => p_cfl
         smoother_dcfl => p_dcfl
      else
         smoother_dt   => dt
         smoother_cfl  => cfl
         smoother_dcfl => dcfl
      end if

      select case (Preconditioner)
      case (PRECONDIIONER_LTS)
         ! compute LTS
         if (Compute_dt) then
            call MaxTimeStep( self=this % p_sem, cfl=smoother_cfl, dcfl=smoother_dcfl , MaxDt=smoother_dt, MaxDtVec=this % lts_dt)
         else
            error stop "FASMultigrid :: LTS needs cfd & dcfl."
         end if

         if( this% p_sem% mesh% IBM% TimePenal ) this% p_sem% mesh% IBM% penalization = smoother_dt

         select case (Smoother)
            ! Euler Smoother
            case (Euler_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )

                  call TakeExplicitEulerStep ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dt_vec=this % lts_dt, dts=DualTimeStepping, global_dt=dt )

                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )
               end do
            ! RK3 smoother
            case (RK3_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )

                  call TakeRK3Step ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dt_vec=this % lts_dt, dts=DualTimeStepping, global_dt=dt )

                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )
               end do
            ! RK5 smoother
            case (RK5_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )

                  call TakeRK5Step ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dt_vec=this % lts_dt, dts=DualTimeStepping, global_dt=dt )

                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )
               end do
            ! RK5 smoother opt for Steady State
            case (RKOpt_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )

                  call TakeRKOptStep ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, N_STAGES=erk_order, dt_vec=this % lts_dt, dts=DualTimeStepping, global_dt=dt )

                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )
               end do
            case default
               error stop "FASMultigrid :: No smoother defined for the multigrid."
         end select ! Smoother
      case (PRECONDIIONER_NONE)
         select case (Smoother)
            ! Euler Smoother
            case (Euler_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  if (Compute_dt) call MaxTimeStep(self=this % p_sem, cfl=smoother_cfl, dcfl=smoother_dcfl, MaxDt=smoother_dt )

                  if ( this% p_sem% mesh% IBM% TimePenal ) this% p_sem% mesh% IBM% penalization = smoother_dt
                  if( this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )

                  call TakeExplicitEulerStep ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dts=DualTimeStepping, global_dt=dt )

                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )
               end do
            ! RK3 smoother
            case (RK3_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  if (Compute_dt) call MaxTimeStep(self=this % p_sem, cfl=smoother_cfl, dcfl=smoother_dcfl, MaxDt=smoother_dt )

                  if ( this% p_sem% mesh% IBM% TimePenal ) this% p_sem% mesh% IBM% penalization = smoother_dt
                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )

                  call TakeRK3Step ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dts=DualTimeStepping, global_dt=dt )

                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )
               end do
            ! RK5 smoother
            case (RK5_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  if (Compute_dt) call MaxTimeStep(self=this % p_sem, cfl=smoother_cfl, dcfl=smoother_dcfl, MaxDt=smoother_dt )

                  if ( this% p_sem% mesh% IBM% TimePenal ) this% p_sem% mesh% IBM% penalization = smoother_dt
                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )

                  call TakeRK5Step ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dts=DualTimeStepping, global_dt=dt )

                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )
               end do
            ! RK Opt smoother
            case (RKOpt_SMOOTHER)
               do sweep = 1, SmoothSweeps
                  if (Compute_dt) call MaxTimeStep(self=this % p_sem, cfl=smoother_cfl, dcfl=smoother_dcfl, MaxDt=smoother_dt )

                  if ( this% p_sem% mesh% IBM% TimePenal ) this% p_sem% mesh% IBM% penalization = smoother_dt
                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )

                  call TakeRKOptStep ( mesh=this % p_sem % mesh, particles=this % p_sem % particles, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, N_STAGES=erk_order, dts=DualTimeStepping, global_dt=dt )

                  if(  this% p_sem% mesh% IBM% active ) call this% p_sem% mesh% IBM% SemiImplicitCorrection( this% p_sem% mesh% elements, t, smoother_dt )
               end do
!
!           Implicit smoothers
!           ------------------
            case (IRK_SMOOTHER)
               error stop "FASMultigrid :: IRK Smoother not ready."
            case (SGS_SMOOTHER)

               call this % p_sem % mesh % storage % local2globalq (this % p_sem % mesh % storage % NDOF)

               do sweep = 1, SmoothSweeps
                  if (Compute_dt) call MaxTimeStep(self=this % p_sem, cfl=smoother_cfl, dcfl=smoother_dcfl, MaxDt=smoother_dt )
                  if( this% p_sem% mesh% IBM% TimePenal ) this% p_sem% mesh% IBM% penalization = smoother_dt
                  call TakeSGSStep ( this=this, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dts=DualTimeStepping, global_dt=dt )
               end do

            case (ILU_SMOOTHER)

               call this % p_sem % mesh % storage % local2globalq (this % p_sem % mesh % storage % NDOF)

               do sweep = 1, SmoothSweeps
                  if (Compute_dt) call MaxTimeStep(self=this % p_sem, cfl=smoother_cfl, dcfl=smoother_dcfl, MaxDt=smoother_dt )
                  if( this% p_sem% mesh% IBM% TimePenal ) this% p_sem% mesh% IBM% penalization = smoother_dt
                  call TakeILUStep ( this=this, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dts=DualTimeStepping, global_dt=dt )
               end do

            case (BIRK5_SMOOTHER)

               call this % p_sem % mesh % storage % local2globalq (this % p_sem % mesh % storage % NDOF)

               do sweep = 1, SmoothSweeps
                  if (Compute_dt) call MaxTimeStep(self=this % p_sem, cfl=smoother_cfl, dcfl=smoother_dcfl, MaxDt=smoother_dt )
                  if( this% p_sem% mesh% IBM% TimePenal ) this% p_sem% mesh% IBM% penalization = smoother_dt
                  call TakeBIRK5Step ( this=this, t=t, deltaT=smoother_dt, &
                     ComputeTimeDerivative=ComputeTimeDerivative, dts=DualTimeStepping, global_dt=dt )
               end do

            case default
               error stop "FASMultigrid :: Smoother not specified."
         end select ! Smoother
      end select ! Preconditioner

   end subroutine Smooth
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine TakeBIRK5Step( this, t, deltaT, ComputeTimeDerivative, dts, global_dt )
!
!     ----------------------------------
!     5th order diagonally implicit Runge-Kutta scheme from Bassi 2009
!     ----------------------------------
!
      implicit none
!
!     -----------------
!     Input parameters:
!     -----------------
!
      class(FASMultigrid_t)  ,intent(inout), target :: this
      real(KIND=RP)   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)    :: ComputeTimeDerivative
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
!
!     ---------------
!     Local variables
!     ---------------
!
      real(kind=rp), dimension(5) :: a = (/0.2_RP, 0.25_RP, 0.333333_RP, 0.5_RP, 1.0_RP/)
      integer :: k, id, i
      real(kind=RP), allocatable :: x_loc(:) ! Local x
!-------------------------------------------------------------

      this % dQ = 0._RP
      do k = 1,5

         tk = t + a(k)*deltaT

         call ComputeTimeDerivative( this % p_sem % mesh, this % p_sem % particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(this % p_sem % mesh, t, global_dt)
         end if
         call this % p_sem % mesh % storage % local2globalqdot (this % p_sem % mesh % storage % NDOF)


         select type (Adense => this % A)
            type is (DenseBlockDiagMatrix_t)
!$omp parallel do private(x_loc) schedule(runtime)
            do id = 1, SIZE( this % p_sem % mesh % elements )
               allocate( x_loc(Adense % BlockSizes(id)) )
               call SolveLU(ALU      = Adense % Blocks(id) % Matrix, &
                           LUpivots = this % LUpivots(id) % v, &
                           x = x_loc, &
                           b = this % p_sem % mesh % storage % Qdot(Adense % BlockIdx(id):Adense % BlockIdx(id+1)-1) * a(k) )
!$omp critical
               this % dQ(Adense % BlockIdx(id):Adense % BlockIdx(id+1)-1) = x_loc
!$omp end critical
               deallocate(x_loc)
            end do ! id
!$omp end parallel do
         end select

         this % p_sem % mesh % storage % Q = this % p_sem % mesh % storage % Q - this % dQ
         call this % p_sem % mesh % storage % global2localq
      end do ! k

!$omp parallel do schedule(runtime)
      do k=1, this % p_sem % mesh % no_of_elements
         if ( any(isnan(this % p_sem % mesh % elements(k) % storage % Q))) then
            print*, "Numerical divergence obtained in solver."
            call exit(99)
         endif
      end do
!$omp end parallel do

   end subroutine TakeBIRK5Step
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine TakeSGSStep( this, t, deltaT, ComputeTimeDerivative, dts, global_dt )
!
!     ----------------------------------
!     5th order diagonally implicit Runge-Kutta scheme from Bassi 2009
!     ----------------------------------
!
      implicit none
!
!     -----------------
!     Input parameters:
!     -----------------
!
      class(FASMultigrid_t)  ,intent(inout), target :: this
      real(KIND=RP)   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)    :: ComputeTimeDerivative
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
!
!     ---------------
!     Local variables
!     ---------------
!
      integer :: k, id, i
      integer :: stat
      real(kind=rp), dimension(5) :: a = (/0.2_RP, 0.25_RP, 0.333333_RP, 0.5_RP, 1.0_RP/)
!-------------------------------------------------------------

      this % dQ  = 0._RP
      this % dQ0 = 0._RP
      this % SGS_RHS = 0.0_RP

      do k = 1,5

         tk = t + a(k)*deltaT

         call ComputeTimeDerivative( this % p_sem % mesh, this % p_sem % particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(this % p_sem % mesh, t, global_dt)
         end if
         call this % p_sem % mesh % storage % local2globalqdot (this % p_sem % mesh % storage % NDOF)


         select type (Acsr => this % A)
            type is (csrMat_t)

#ifdef HAS_MKL

            do i=1,kSGS

               Acsr % mkl_options % trans = SPARSE_OPERATION_NON_TRANSPOSE
               Acsr % mkl_options % descrA % type = SPARSE_MATRIX_TYPE_TRIANGULAR
               Acsr % mkl_options % descrA % mode = SPARSE_FILL_MODE_LOWER
               Acsr % mkl_options % descrA % diag = SPARSE_DIAG_NON_UNIT
               this % SGS_RHS = Acsr % LMatVecMul( this % dQ,  .false., 1 )
               this % SGS_RHS = this % SGS_RHS + this % p_sem % mesh % storage % Qdot * a(k)
               call Acsr % ForwSub(this % SGS_RHS, this % dQ0, this % DimPrb, 1.0_RP)

               Acsr % mkl_options % descrA % mode = SPARSE_FILL_MODE_UPPER
               this % SGS_RHS = Acsr % UMatVecMul( this % dQ0,  .false., 1 )
               this % SGS_RHS = this % SGS_RHS + this % p_sem % mesh % storage % Qdot * a(k)
               call Acsr % BackSub(this % SGS_RHS, this % dQ, this % DimPrb, 1.0_RP)

            end do


#else
            error stop "SGS smoother needs MKL."
#endif

         end select

         this % p_sem % mesh % storage % Q = this % p_sem % mesh % storage % Q - this % dQ
         call this % p_sem % mesh % storage % global2localq

      end do

!$omp parallel do schedule(runtime)
      do k=1, this % p_sem % mesh % no_of_elements
         if ( any(isnan(this % p_sem % mesh % elements(k) % storage % Q))) then
            print*, "Numerical divergence obtained in solver."
            call exit(99)
         endif
      end do
!$omp end parallel do

   end subroutine TakeSGSStep
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine TakeILUStep( this, t, deltaT, ComputeTimeDerivative, dts, global_dt )
!
!     ----------------------------
!     Semi-implicit ILU0 smoother.
!     ----------------------------
!
      implicit none
!
!     -----------------
!     Input parameters:
!     -----------------
!
      class(FASMultigrid_t)  ,intent(inout), target :: this
      real(KIND=RP)   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)    :: ComputeTimeDerivative
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
!
!     ---------------
!     Local variables
!     ---------------
!
      integer :: k, id, i, stat
      real(kind=rp), dimension(5) :: a = (/0.2_RP, 0.25_RP, 0.333333_RP, 0.5_RP, 1.0_RP/)
!-------------------------------------------------------------

      this % dQ  = 0._RP
      this % dQ0 = 0._RP

      do k = 1,5

         tk = t + a(k)*deltaT

         call ComputeTimeDerivative( this % p_sem % mesh, this % p_sem % particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(this % p_sem % mesh, t, global_dt)
         end if
         call this % p_sem % mesh % storage % local2globalqdot (this % p_sem % mesh % storage % NDOF)

         select type (Acsr => this % A)
            type is (csrMat_t)

#ifdef HAS_MKL

            Acsr % mkl_options % trans = SPARSE_OPERATION_NON_TRANSPOSE
            Acsr % mkl_options % descrA % type = SPARSE_MATRIX_TYPE_TRIANGULAR
            Acsr % mkl_options % descrA % mode = SPARSE_FILL_MODE_LOWER
            Acsr % mkl_options % descrA % diag = SPARSE_DIAG_UNIT
            call Acsr % ForwSub(this % p_sem % mesh % storage % Qdot * a(k), this % dQ0, this % DimPrb, 1.0_RP)

            Acsr % mkl_options % descrA % mode = SPARSE_FILL_MODE_UPPER
            Acsr % mkl_options % descrA % diag = SPARSE_DIAG_NON_UNIT
            call Acsr % BackSub(this % dQ0, this % dQ, this % DimPrb, 1.0_RP)

#else
            error stop " FASMultigridClass :: ILU smoother needs MKL."
#endif

         end select

         this % p_sem % mesh % storage % Q = this % p_sem % mesh % storage % Q - this % dQ
         call this % p_sem % mesh % storage % global2localq

      end do

!$omp parallel do schedule(runtime)
      do k=1, this % p_sem % mesh % no_of_elements
         if ( any(isnan(this % p_sem % mesh % elements(k) % storage % Q))) then
            print*, "Numerical divergence obtained in solver."
            call exit(99)
         endif
      end do
!$omp end parallel do
   end subroutine TakeILUStep
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine FAS_SetPreviousSolution(this,lvl)
      implicit none
      !-------------------------------------------------------------
      class(FASMultigrid_t), target, intent(inout) :: this     !<> Anisotropic FAS multigrid
      integer                      , intent(in)    :: lvl      !<  Current level
      !-------------------------------------------------------------
      integer :: N1(3), N2(3), eID
      !-------------------------------------------------------------

!
!     Set the previous solution in this level
!     ---------------------------------------

      if (lvl == MGlevels) then
         call BDF_SetPreviousSolution(this % p_sem % mesh % storage)
      else
         call BDF_SetPreviousSolution(this % p_sem % mesh % storage, NotANewStep = .TRUE. )
      end if
!
!     Send the solution to the next (coarser) level
!     ---------------------------------------------

      if (lvl > 1) then
!$omp parallel do private(N1,N2) schedule(runtime)
         do eID = 1, nelem
            N1 = this % p_sem % mesh % elements (eID) % Nxyz
            N2 = this % Child % p_sem % mesh % elements (eID) % Nxyz

!           Restrict solution
!           -----------------
            call Interp3DArrays(NCONS, N1, this % p_sem % mesh % elements(eID) % storage % Q, &
                                       N2, this % Child % p_sem % mesh % elements(eID) % storage % Q )
         end do
!$omp end parallel do

         call FAS_SetPreviousSolution(this % Child,lvl-1)
      end if
   end subroutine FAS_SetPreviousSolution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine FAS_UpdateInitialization(this, lvl)
      implicit none
      !-------------------------------------------------------------
      class(FASMultigrid_t), target, intent(inout) :: this     !<> Anisotropic FAS multigrid
      integer              , intent(in)    :: lvl      !<  Current multigrid level
      !-------------------------------------------------------------
      !-------------------------------------------------------------

      if (mg_Initialization) then
         Preconditioner = ini_Preconditioner(1)
         Smoother = ini_Smoother(1)
         this % computeA = .false.

         if (DualTimeStepping) then
            p_cfl = ini_cfl(1)
            p_dcfl = ini_cfl(2)
         else
            cfl = ini_cfl(1)
            dcfl = ini_cfl(2)
         end if

         if (MAXVAL(ComputeMaxResiduals(this % p_sem % mesh)) .le. ini_res) then
            mg_Initialization = .false.
            call computeA_AllLevels(this,lvl)
            Preconditioner = ini_Preconditioner(2)
            Smoother = ini_Smoother(2)
            if (DualTimeStepping) then
               p_cfl = ini_cfl(3)
               p_dcfl = ini_cfl(4)
            else
               cfl = ini_cfl(3)
               dcfl = ini_cfl(4)
            end if
         end if
      end if

   end subroutine FAS_UpdateInitialization
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine FAS_ComputeAndFactorizeJacobian(this, lvl, t, invdt, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
      implicit none
      !-------------------------------------------------------------
      class(FASMultigrid_t), target, intent(inout) :: this     !<> FAS multigrid
      integer                      , intent(in)    :: lvl      !<  Current level
      real(kind=RP)        , intent(in)    :: t                !<  Simulation time
      real(kind=RP)        , intent(in)    :: invdt            !<  Simulation time
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivativeIsolated
      !-------------------------------------------------------------
      integer :: k
      !-------------------------------------------------------------

      if (this % computeA) then

         select type (Amat => this % A)

!        Compute full Jacobian in CSR format
!        -----------------------------------
         type is (csrMat_t)
            call this % Jacobian % Compute (this % p_sem, NCONS, t, this % A, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
            call Amat % shift(invdt)
            ! call Amat % Visualize('A.dat')

!           Factorization (depends on the smoother used)
!           --------------------------------------------
            select case (SMOOTHER)
            case (SGS_SMOOTHER)
            case (ILU_SMOOTHER)
               call Amat % ILU0Factorization()
            end select

#ifdef HAS_MKL
            call Amat % CreateMKL()
#else
            error stop "Full Jacobian smoothers need MKL."
#endif

!        Compute Jacobian blocks in local dense format
!        ---------------------------------------------
         type is (DenseBlockDiagMatrix_t)
            call this % Jacobian % Compute (sem=this % p_sem, nEqn=NCONS, time=t, matrix=this % A, TimeDerivative=ComputeTimeDerivative, &
            TimeDerivativeIsolated=ComputeTimeDerivativeIsolated, BlockDiagonalized=.true.)
            call Amat % shift(invdt)

!           Local LU factorization
!           ----------------------
!$omp parallel do schedule(runtime)
            do k=1, nelem
               call ComputeLUandOverwrite (A       = Amat % Blocks(k) % Matrix, &
                                          LUpivots = this % LUpivots(k) % v)
            end do
!$omp end parallel do
         end select

         this % computeA = .false.

      end if

   end subroutine FAS_ComputeAndFactorizeJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   
   recursive subroutine FAS_MovingIBM( Solver, dt, lvl )
   
      implicit none
      !-arguments------------------------------------------
      type(FASMultigrid_t), intent(inout) :: Solver
      real(kind=RP),        intent(in)    :: dt
      integer,              intent(in)    :: lvl

      call  Solver% p_sem% mesh% IBM% MoveBody( Solver% p_sem% mesh% elements,       &
                                                Solver% p_sem% mesh% no_of_elements, &
                                                Solver% p_sem% mesh% NDOF,           &
                                                .true., dt                           )

      if( lvl > 1 ) call FAS_movingIBM( Solver% Child, dt, lvl-1  )
   
   end subroutine FAS_MovingIBM
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module FASMultigridClass