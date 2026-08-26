!//////////////////////////////////////////////////////
!
!  Class for solving a linear system obtained from implicit time discretization. 
!
!  Usage:
!
!  Variables for the control file:  
!  no_levels :: number of MG levels, IF NOT (no_levels = 2)        
!  define levels x :: hard define N on each level, IF NOT \Delta N_{x} = 1
!  define levels y ...
!  define levels z ...
!  ----------------
!  TODO:
!  ----------------
!  1. Deallacotion: 
!     a. DGSem on each level
!     b. Jac on each level
!     c. Prol/Rest operators 
!  2. Generalize MG to operate on nonconforming p-mesh (elements with different pol. orders)
!  3. Add more smoothers - so far only diag-jacobi implemented.
!  ----------------
!
!//////////////////////////////////////////////////////
module LinearMultigridSolverClass
   use GenericLinSolverClass
   use CSRMatrixClass
   use SMConstants
   use PhysicsStorage
   use StorageClass
   use PolynomialInterpAndDerivsModule
   use GaussQuadrature
   use DGSEMClass
   use TimeIntegratorDefinitions
   use MatrixClass
   use NumericalJacobian      , only: NumJacobian_t
   use AnalyticalJacobian     , only: AnJacobian_t
   use JacobianComputerClass  , only: JacobianComputer_t, GetJacobianFlag
   use StopWatchClass
   use DenseMatUtilities
   ! use GenericSmoother
#include "Includes.h"

   implicit none

   
   private
   public LinearMultigridSolver_t, TemporaryElementStorage_t
   public SOLVER_PMG_NONE, SOLVER_GMRES_PMG
   public JACOBIANCOMP_MF, JACOBIANCOMP_MB
   public KSP_PRECONDITIONER_NONE, KSP_PRECONDITIONER_BJ, KSP_PRECONDITIONER_PMG 

   integer, parameter ::  SOLVER_PMG_NONE = 0
   integer, parameter ::  SOLVER_GMRES_PMG = 1

   integer, parameter ::  JACOBIANCOMP_MB = 0
   integer, parameter ::  JACOBIANCOMP_MF = 1

   integer, parameter ::  KSP_PRECONDITIONER_NONE = 0
   integer, parameter ::  KSP_PRECONDITIONER_BJ = 1
   integer, parameter ::  KSP_PRECONDITIONER_PMG = 2
   integer, parameter ::  KSP_PRECONDITIONER_PJ = 3
   integer, parameter ::  KSP_PRECONDITIONER_ILU = 4

   integer, parameter :: S_NOTDEF     = 0
   integer, parameter :: S_POINTJAC   = 1
   integer, parameter :: S_BLOCKJAC   = 2
   integer, parameter :: S_ILU        = 3
   integer, parameter :: S_BILU       = 4
!
!  ------------------------------------------------
!  Blocks For Smoother
!  ------------------------------------------------
!
   type :: BlockPrec_t
!-----Variables-----------------------------------------------------------
      real(KIND=RP), dimension(:,:), allocatable :: PLU        ! LU factorization of elemental preconditioner matrix
      integer      , dimension(:)  , allocatable :: LUpivots   ! LU pivots
   end type BlockPrec_t
!
!  ------------------------------------------------
!  ILU smoother
!  ------------------------------------------------
!
   type :: ILUSmooth_t
!-----Variables-----------------------------------------------------------
      type(csrMat_t) :: A ! ILU matrix
   contains
!-----Subroutines-----------------------------------------------------------
      procedure                                  :: Construct => MGS_ConstructILU
      procedure                                  :: Destruct  => MGS_DestructILU
   end type ILUSmooth_t
!
!  ------------------------------------------------
!  Block-Jacobi smoother
!  ------------------------------------------------
!
   type :: BJSmooth_t
!-----Variables-----------------------------------------------------------
      class(DenseBlockDiagMatrix_t), pointer :: A_p           ! pointer to Block-Jacobian matrix
      type(BlockPrec_t), allocatable         :: BlockPrec(:)
   contains
!-----Subroutines-----------------------------------------------------------
      procedure                                  :: Construct => MGS_ConstructBlockJacobi
      procedure                                  :: Destruct  => MGS_DestructBlockJacobi
   end type BJSmooth_t
!
!  ------------------------------------------------
!  Local temporary element storage.
!  ------------------------------------------------
!
   type :: TemporaryElementStorage_t
!-----Variables-----------------------------------------------------------
      real(kind=RP), dimension(:,:,:,:), allocatable :: X 
      real(kind=RP), dimension(:,:,:,:), allocatable :: R 
   end type TemporaryElementStorage_t
!
!  ------------------------------------------------
!  Local temporary element storage.
!  ------------------------------------------------
!
   type :: KSPForMG_t
!-----Variables-----------------------------------------------------------
      integer                    :: KrylovSpace    ! size of the Krylov Space
      integer                    :: Preconditioner ! preconditioner type 
      real(kind=RP), allocatable :: H(:,:)
      real(kind=RP), allocatable :: W(:)
      real(kind=RP), allocatable :: V(:,:)
      real(kind=RP), allocatable :: Z(:,:)
      real(kind=RP), allocatable :: Y(:)
      real(kind=RP), allocatable :: cc(:)
      real(kind=RP), allocatable :: ss(:)
      real(kind=RP), allocatable :: g(:)
   end type KSPForMG_t
!
!  ------------------------------------------------
!  Multigrid type (Linear Solver class extension)
!  ------------------------------------------------
!
   type, extends(GenericLinSolver_t) :: LinearMultigridSolver_t
!-----Variables-----------------------------------------------------------
      type(KSPForMG_t)                              :: KSP                   ! KSP 
      type(csrMat_t)                                :: A                     ! Matrix to solve
      type(BJSmooth_t)            ,     allocatable :: BJSmoother            ! BJ smoother
      type(ILUSmooth_t)           ,     allocatable :: ILUSmoother           ! ILU smoother
      real(kind=RP), dimension(:) ,     allocatable :: x                     ! Solution vector
      real(kind=RP), dimension(:) ,     allocatable :: b                     ! Right hand side
      real(kind=RP), dimension(:) ,     allocatable :: r                     ! Residual
      real(kind=RP)                                 :: rnorm                 ! L2 norm of residual
      real(kind=RP)                                 :: tol                   ! Tolerance
      integer                                       :: maxiter               ! Max. # iterations
      real(kind=RP)                                 :: Ashift                ! Jacobian recompute
      integer                                       :: MGlevel               ! Current level
      integer                                       :: Nx                    ! Polynomial in X
      integer                                       :: Ny                    ! Polynomial in Y
      integer                                       :: Nz                    ! Polynomial in Z
      type(LinearMultigridSolver_t), pointer        :: Child                 ! Coarser level: MGlevel-1
      type(LinearMultigridSolver_t), pointer        :: Parent                ! Finer level: MGlevel+1
      type(TemporaryElementStorage_t) , allocatable :: LocalStorage(:)       ! Storage
      real(kind=RP)                                 :: timesolve             ! Time at the solution
      real(kind=RP)                                 :: dt                    ! dt for the solution
      real(kind=RP), dimension(:) ,     allocatable :: resvec                ! array for the residual
      ! 1D prolongation/restriction operators in each direction 
      real(kind=RP), dimension(:,:) ,   allocatable :: RestX(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: ProlX(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: RestY(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: ProlY(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: RestZ(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: ProlZ(:,:)             
      ! 3D prolongation/restriction operators (1 el, 1 eq)
      real(kind=RP), dimension(:,:) ,   allocatable :: Rest3D(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: Prol3D(:,:)             
      ! 3D prolongation/restriction operators (1 el, nEqn eq)
      real(kind=RP), dimension(:,:) ,   allocatable :: nRest3D(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: nProl3D(:,:)             
      ! 3D mat prolongation/restriction operators size of 1 block-element 
      type(csrMat_t) :: Rest1elCSR            
      type(csrMat_t) :: Prol1elCSR             
      ! 3D full mat prolongation/restriction operators 
      type(csrMat_t) :: RestCSR             
      type(csrMat_t) :: ProlCSR             
      ! Matrix_free operators
      real(kind=RP), allocatable :: F_Ur(:)          ! Qdot at the beginning of solving procedure
      real(kind=RP), allocatable :: Ur(:)            ! Q at the beginning of solving procedure
      ! full Jacobian
      real(kind=RP), dimension(:,:) ,   allocatable :: Afull(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: Pfull(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: Rfull(:,:)             
      ! test for the Jacobian 
      real(kind=RP), dimension(:,:) ,   allocatable :: RestJacX(:,:)
      real(kind=RP), dimension(:,:) ,   allocatable :: RestJacY(:,:)
      real(kind=RP), dimension(:,:) ,   allocatable :: RestJacZ(:,:)
      
   contains
!-----Subroutines-----------------------------------------------------------
      ! Main routines
      procedure :: Construct            => MG_Construct 
      procedure :: Destroy              => MG_Destruct 
      procedure :: SetRHS               => MG_SetRHS
      procedure :: GetX                 => MG_GetX
      procedure :: Solve                => MG_Solve 
      procedure :: Getrnorm             => MG_Getrnorm
      procedure :: Getxnorm             => MG_Getxnorm
      procedure :: SetRHSValue          => MG_SetRHSValue 
      procedure :: SetRHSValues         => MG_SetRHSValues 
      procedure :: SetOperatorDt        => MG_SetOperatorDt
      procedure :: ReSetOperatorDt      => MG_ReSetOperatorDt

      ! TBD
      procedure :: MG_JacVec
      procedure :: MG_PrecVec
      procedure :: MG_BJsmooth
      procedure :: MG_ILUsmooth
      procedure :: MG_PJsmooth
      procedure :: MG_smooth

      ! Interpolation routines
      procedure :: MG_CreateProlongationOperator
      procedure :: MG_CreateRestrictionOperator
      !procedure :: MG_Create3DProlongationMatrix
      !procedure :: MG_Create3DRestrictionMatrix
      !procedure :: prolong              => MG_1DProlongation
      !procedure :: restrict             => MG_1DRestriction
      procedure :: MG_1DProlongation
      procedure :: MG_1DRestriction
      procedure :: MG_3DProlongation
      procedure :: MG_3DRestriction
      procedure :: MG_JacRestriction
      procedure :: MG_Create1elCSRInterpolationMats

      ! Subsolver/preconditioner routines
      procedure :: MG_KSP_Construct
      procedure :: MG_KSP_Destruct 
      procedure :: MG_Solve_MB_PMG_NONE
      procedure :: MG_Solve_MB_GMRES_PMG

   end type LinearMultigridSolver_t
!
!  ----------------
!  Module variables
!  ----------------
!
   integer                            :: no_levels        ! Total number of multigrid levels
   integer, dimension(:), allocatable :: MG_levels_x      ! multigrid levels
   integer, dimension(:), allocatable :: MG_levels_y      ! multigrid levels
   integer, dimension(:), allocatable :: MG_levels_z      ! multigrid levels
   integer, dimension(:), allocatable :: pre_smooths      ! pre-smoothing operations
   integer, dimension(:), allocatable :: pos_smooths      ! post-smoothing operations
   integer                            :: nelem            ! Number of elements
   integer                            :: S_SMOOTHER       ! Smoother type
   integer                            :: S_SOLVER         ! Solver type
   integer                            :: S_PRECONDITIONER ! Preconditioner type
   integer                            :: S_MATCOMP        ! matrix free/matrix based
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Construct(this,DimPrb,globalDimPrb,nEqn,controlVariables,sem,MatrixShiftFunc)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t),  intent(inout), target :: this
      integer                  , intent(in)            :: DimPrb
      integer                  , intent(in)            :: globalDimPrb        
      integer                  , intent(in)            :: nEqn
      type(FTValueDictionary)  , intent(in), optional  :: controlVariables
      type(DGSem), target                  , optional  :: sem
      procedure(MatrixShift_FCN)                       :: MatrixShiftFunc     ! TODO: Make this optional
      procedure(ComputeTimeDerivative_f)               :: ComputeTimeDerivative
!-----Local-Variables-----------------------------------------------------
      character(len=LINE_LENGTH) :: pc
      integer                    :: i
      real(kind=RP)              :: tmp_1
!  -----------------------------------------------------------------------

      ! Check dims, MPI and allocate Jacobian
      call this % GenericLinSolver_t % construct(DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
      MatrixShift => MatrixShiftFunc ! FIXME: ask whether we need this 

!     Set variables from controlVariables
!     --------------------------------------------------------------------
      if (.not. present(controlVariables)) stop 'Fatal error: MultigridSolver needs controlVariables.'
      if (present(controlVariables)) then 

        ! Read # of multigrid coarse grids
!       ------------------------------------------------------------------
        if ( controlVariables % containsKey("multigrid levels") ) then
            no_levels = controlVariables % integerValueForKey("multigrid levels")
        else
            no_levels = 2
        end if

        allocate(MG_levels_x(no_levels))
        allocate(MG_levels_y(no_levels))
        allocate(MG_levels_z(no_levels))
        allocate(pre_smooths(no_levels))
        allocate(pos_smooths(no_levels))
!       ------------------------------------------------------------------

        ! Read user-specified polynomial order on each level
!       ------------------------------------------------------------------
        if ( controlVariables % containsKey("multigrid levels") ) then
            ! Levels in X
            pc = controlVariables % StringValueForKey("define levels x",LINE_LENGTH)
            do i = 1, len_trim(pc)
              read(pc(i:i),'(i1)') MG_levels_x(i)
            end do
            ! Levels in Y
            pc = controlVariables % StringValueForKey("define levels y",LINE_LENGTH)
            do i = 1, len_trim(pc)
              read(pc(i:i),'(i1)') MG_levels_y(i)
            end do
            ! Levels in Z
            pc = controlVariables % StringValueForKey("define levels z",LINE_LENGTH)
            do i = 1, len_trim(pc)
              read(pc(i:i),'(i1)') MG_levels_z(i)
            end do
        else
            ERROR stop ':: FIXME: Default multigrid levels NOT defined.' ! TODO
        end if
!       ------------------------------------------------------------------

        ! Read type of solver and preconditioner 
!       ------------------------------------------------------------------
        if ( controlVariables % containsKey("multigrid type") ) then
        select case ( trim( controlVariables % StringValueForKey("multigrid type",LINE_LENGTH) ) )
           case ('gmres-none')
             S_SOLVER = SOLVER_GMRES_PMG
             S_PRECONDITIONER = KSP_PRECONDITIONER_NONE
           case ('gmres-pmg')
             S_SOLVER = SOLVER_GMRES_PMG
             S_PRECONDITIONER = KSP_PRECONDITIONER_PMG
           case ('gmres-bj')
             S_SOLVER = SOLVER_GMRES_PMG
             S_SMOOTHER = S_BLOCKJAC 
             S_PRECONDITIONER = KSP_PRECONDITIONER_BJ
           case ('gmres-ilu')
             S_SOLVER = SOLVER_GMRES_PMG
             S_SMOOTHER = S_ILU 
             S_PRECONDITIONER = KSP_PRECONDITIONER_ILU
           case ('gmres-pj')
             S_SOLVER = SOLVER_GMRES_PMG
             S_SMOOTHER = S_POINTJAC 
             S_PRECONDITIONER = KSP_PRECONDITIONER_PJ
           case ('pmg-none')
             S_SOLVER = SOLVER_PMG_NONE
             !S_PRECONDITIONER = KSP_PRECONDITIONER_NONE
           case default
             ERROR Stop "LinearMultigridSolver :: Wrong solver "
        end select
        end if
!       ------------------------------------------------------------------

        ! Specify the version of the solver (matrix-based or matrix-free) 
!       ------------------------------------------------------------------
        if ( controlVariables % containsKey("jacobian assembly") ) then
        select case ( trim( controlVariables % StringValueForKey("jacobian assembly",LINE_LENGTH) ) )
           case ('matrix-free')
             S_MATCOMP = JACOBIANCOMP_MF
           case ('matrix-based')
             S_MATCOMP = JACOBIANCOMP_MB
           case default
             ERROR Stop "MultigridSolver :: Wrong jacobian assemble key"
        end select
        end if
!       ------------------------------------------------------------------

        ! Smoother
!       ------------------------------------------------------------------
        if ( controlVariables % containsKey("smoother") ) then
        select case ( trim( controlVariables % StringValueForKey("smoother",LINE_LENGTH) ) )
           case ('point-jacobi')
             S_SMOOTHER = S_POINTJAC
             print *, 'GenericSmoother :: WARNING p. Jacobi smoother is not suitable for DGSEM Jacobians'
           case ('block-jacobi')
             S_SMOOTHER = S_BLOCKJAC
           case ('ilu')
             S_SMOOTHER = S_ILU
           case default
             ERROR Stop "MultigridSolver :: Wrong smoother "
        end select
        end if
!       ------------------------------------------------------------------

        ! Specify number of pre- and post-smoothing operations
!       ------------------------------------------------------------------
        if ( controlVariables % containsKey("pre smooths") ) then
          pc = controlVariables % StringValueForKey("pre smooths",LINE_LENGTH)
          pre_smooths = getnosmooths(pc,no_levels)
        end if

        if ( controlVariables % containsKey("post smooths") ) then
          pc = controlVariables % StringValueForKey("post smooths",LINE_LENGTH)
          pos_smooths = getnosmooths(pc,no_levels)
        end if
!       ------------------------------------------------------------------
      end if 
!     --------------------------------------------------------------------

      nelem = size(sem % mesh % elements)
      this % p_sem => sem ! this is for generic linear solver class
      this % DimPrb = DimPrb
  
      ! Construct initial variables for the finest level and call recursive constructor
      this % MGlevel = no_levels
      this % Nx = MG_levels_x(no_levels)
      this % Ny = MG_levels_y(no_levels)
      this % Nz = MG_levels_z(no_levels)
      call MG_Levels_Construct(this,no_levels,controlVariables,nEqn)

   end subroutine MG_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine MG_Levels_Construct(Me,lvl,controlVariables,nEqn)
!  ---------------------------------------------------------
!  Recursive suplement to the multigrid constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      type(LinearMultigridSolver_t),  intent(inout), target :: Me
      integer,                  intent(in)            :: nEqn
      integer,                  intent(in)            :: lvl
      class(FTValueDictionary), intent(in)            :: controlVariables                
!-----Local-Variables-----------------------------------------------------
      type(LinearMultigridSolver_t), pointer :: Child_p          ! Pointer to Child
      integer                          :: i,j,k
      integer                          :: nnzs(nelem)
      character(len=1024)              :: filename
!-------------------------------------------------------------------------
      logical                           :: success            

      ! for coarse jacobian
      integer                                 :: JacobianComputation ! Get Jacobian type 

      allocate ( Me % b(Me % DimPrb) ) 
      allocate ( Me % x(Me % DimPrb) ) 
      allocate ( Me % r(Me % DimPrb) ) 

      allocate ( Me % Ur(Me % DimPrb) ) 
      allocate ( Me % F_Ur(Me % DimPrb) ) 

      ! Jacobian
      select case (S_MATCOMP)
         case (JACOBIANCOMP_MB) 
            call Me % A % Construct(num_of_Rows = Me % DimPrb, withMPI = .false.)
            call Me % Jacobian % Configure (Me % p_sem % mesh, nEqn, Me % A)
         case (JACOBIANCOMP_MF)
         case default
           ERROR Stop "MultigridSolver :: Select MATCOMOP."
      end select

      ! Construct smoother
      select case (S_SMOOTHER) 
      case (S_POINTJAC)
      case (S_BLOCKJAC)
         allocate ( Me % BJSmoother )
         call Me % BJSmoother % Construct ( Me % p_sem, Me % Nx, Me % Ny, Me % Nz, nEqn )
         call Me % Jacobian % Configure (Me % p_sem % mesh, nEqn, Me % BJSmoother % A_p)
      case (S_ILU)
         allocate ( Me % ILUsmoother )
         call Me % ILUSmoother % Construct (Me % A) 
      case default
         ERROR Stop "Error! No smoother set."
      end select 

      ALLOCATE ( Me % LocalStorage(nelem))
      do k = 1, nelem
         allocate(Me % LocalStorage(k) % X (nEqn,0:Me%Nx,0:Me%Ny,0:Me%Nz))
         allocate(Me % LocalStorage(k) % R (nEqn,0:Me%Nx,0:Me%Ny,0:Me%Nz))
         Me % LocalStorage(k) % X = 0._RP
         Me % LocalStorage(k) % R = 0._RP
      end do   

      ! Construct coarser level
      ! ------------------------------------------------------------------
      if (lvl > 1) then
        ! Allocate the child
        ! ------------------------------------------------------------------
        allocate(Me % Child) ! allocate coarser MG level (child)
        Child_p => Me % Child ! set local pointer to a child
        Me % Child % Parent => Me ! set child's parent pointer to this level (Me)
        Child_p % MGlevel = lvl - 1
        Child_p % Nx = MG_levels_x(lvl - 1)
        Child_p % Ny = MG_levels_y(lvl - 1)
        Child_p % Nz = MG_levels_z(lvl - 1)
        ! ------------------------------------------------------------------

        ! Create coarse DGSem
        ! ------------------------------------------------------------------
        allocate (Child_p % p_sem)
        call Child_p % p_sem % construct (  controlVariables  = controlVariables, &
        polynomialOrder = (/ Child_p%Nx, Child_p%Ny, Child_p%Nz /), success           = success)
        if (.not. success) ERROR STOP "LinearMultigrid: Problem creating coarse solver."

        Child_p % DimPrb = Child_p % p_sem % mesh % NDOF * nEqn 

        ! Prolongation restriction operators
        ! ------------------------------------------------------------------
        allocate (Me  % RestX(0:Child_p%Nx,0:Me%Nx)                 )
        allocate (Child_p % ProlX(0:Me%Nx,0:Child_p%Nx)             )
        allocate (Me  % RestY(0:Child_p%Ny,0:Me%Ny)                 )
        allocate (Child_p % ProlY(0:Me%Ny,0:Child_p%Ny)             )
        allocate (Me  % RestZ(0:Child_p%Nz,0:Me%Nz)                 )
        allocate (Child_p % ProlZ(0:Me%Nz,0:Child_p%Nz)             )
        ! Test for the Jacobian
        allocate (Me  % RestJacX(0:Child_p%Nx,0:Me%Nx)                 )
        allocate (Me  % RestJacY(0:Child_p%Nx,0:Me%Nx)                 )
        allocate (Me  % RestJacZ(0:Child_p%Nx,0:Me%Nx)                 )

        ! 3D Prolongation restriction operators
        ! ------------------------------------------------------------------
        allocate (Me      % Rest3D( 1:( (Child_p%Nx+1) * (Child_p%Ny+1) * (Child_p%Nz+1) ), 1:( (Me%Nx+1) * (Me%Ny+1) * (Me%Nz+1) ) )  )
        allocate (Child_p % Prol3D( 1:( (Me%Nx+1) * (Me%Ny+1) * (Me%Nz+1) ), 1:( (Child_p%Nx+1) * (Child_p%Ny+1) * (Child_p%Nz+1) ) )  )
        allocate (Me      % nRest3D( 1:( (Child_p%Nx+1) * (Child_p%Ny+1) * (Child_p%Nz+1) ) * nEqn, & 
                                        1:( (Me%Nx+1) * (Me%Ny+1) * (Me%Nz+1) )* nEqn )  )
        allocate (Child_p % nProl3D( 1:( (Me%Nx+1) * (Me%Ny+1) * (Me%Nz+1) ) * nEqn, & 
                                        1:( (Child_p%Nx+1) * (Child_p%Ny+1) * (Child_p%Nz+1) * nEqn ) )  )

        call MG_CreateRestrictionOperator  ( Me , .true.)
        ! call MG_CreateRestrictionOperator  ( Me , .true.)
        call MG_CreateProlongationOperator ( Child_p )

        ! ------------------------------------------------------------------
        ! coarse Jacobian
         JacobianComputation = GetJacobianFlag()
         select case (JacobianComputation)
           case (ANALYTICAL_JACOBIAN) ; allocate(AnJacobian_t  :: Child_p % Jacobian)
           case (NUMERICAL_JACOBIAN ) ; allocate(NumJacobian_t :: Child_p % Jacobian)
           case default
              ERROR stop 'Invalid jacobian type. FIXME: '
         end select
         call Child_p % Jacobian % construct(Child_p % p_sem % mesh, nEqn, controlVariables)

        call MG_Levels_Construct(Me % Child,lvl-1,controlVariables,nEqn)
      end if
   end subroutine MG_Levels_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Destruct(this)
!  ---------------------------------------------------------
!  Destructor. 
!  ---------------------------------------------------------
     implicit none
!-----Arguments-----------------------------------------------------------
     class(LinearMultigridSolver_t), intent(inout) :: this

!  -----------------------------------------------------------------------

     call MG_Levels_Destruct(this, no_levels)

      deallocate(MG_levels_x)
      deallocate(MG_levels_y)
      deallocate(MG_levels_z)
      deallocate(pre_smooths)
      deallocate(pos_smooths)

   end subroutine MG_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine MG_Levels_Destruct(Me, lvl)
      !  ---------------------------------------------------------
      !  Destructor. 
      !  ---------------------------------------------------------
      implicit none
      !-----Arguments-----------------------------------------------------------
      type(LinearMultigridSolver_t),  intent(inout), target :: Me
      integer,                  intent(in)            :: lvl
      !  -----------------------------------------------------------------------
      type(LinearMultigridSolver_t), pointer :: Child_p       
      !  -----------------------------------------------------------------------
      
      deallocate ( Me % b ) 
      deallocate ( Me % x ) 
      deallocate ( Me % r ) 
      deallocate ( Me % Ur ) 
      deallocate ( Me % F_Ur )
      
      ! Jacobian
      select case (S_MATCOMP)
         case (JACOBIANCOMP_MB) 
            call Me % A % destruct
         case (JACOBIANCOMP_MF)
         case default
      end select

      ! Construct smoother
      select case (S_SMOOTHER) 
      case (S_POINTJAC)
      case (S_BLOCKJAC)
         call Me % BJSmoother % Destruct ! ( )
         deallocate ( Me % BJSmoother )
      case (S_ILU)
         call Me % ILUSmoother % Destruct ! (Me % A) 
         deallocate ( Me % ILUsmoother )
      case default
      end select

      if (lvl > 1) then

         Child_p => Me % Child 
         Me % Child % Parent => Me

         deallocate (Me  % RestX)
         deallocate (Child_p % ProlX)
         deallocate (Me  % RestY)
         deallocate (Child_p % ProlY)
         deallocate (Me  % RestZ)
         deallocate (Child_p % ProlZ)
         deallocate (Me  % RestJacX)
         deallocate (Me  % RestJacY)
         deallocate (Me  % RestJacZ)
   
         deallocate (Me      % Rest3D)
         deallocate (Child_p % Prol3D)
         deallocate (Me      % nRest3D)
         deallocate (Child_p % nProl3D)

         call MG_Levels_Destruct(Me % Child,lvl-1)

         deallocate (Child_p % p_sem)
         deallocate(Me % Child)
      end if
   end subroutine MG_Levels_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine MG_UpdateInfo(this, lvl, dt, time)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      real(kind=rp),                    intent(in)    :: time
      real(kind=rp),                    intent(in)    :: dt
      integer,                          intent(in)    :: lvl
!  -----------------------------------------------------------------------

      this % dt = dt
      this % timesolve = time
      if (lvl > 1) then  
         call MG_UpdateInfo(this % child, lvl-1,dt,time)
      end if

   end subroutine
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Solve(this, nEqn, nGradEqn, ComputeTimeDerivative, tol, maxiter, time, dt, ComputeA)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target   , intent(inout) :: this
      type(LinearMultigridSolver_t) , pointer                  :: pMG     
      integer                            , intent(in)    :: nEqn, nGradEqn
      procedure(ComputeTimeDerivative_f)                 :: ComputeTimeDerivative
      real(kind=rp)           , optional                 :: tol
      integer                 , optional                 :: maxiter
      real(kind=rp)           , optional                 :: time
      real(kind=rp)           , optional                 :: dt
      logical                 , optional , intent(inout) :: ComputeA
!-----Local-Variables-----------------------------------------------------
      class(csrMat_t), pointer :: pA
      integer                  :: niter
      real(kind=rp)            :: tmpsize
      integer                  :: i, j, k
      character(len=1024)      :: filename
      logical                  :: file_exists
      integer                  :: solver_type             
!  -----------------------------------------------------------------------

      call MG_UpdateInfo(this,no_levels, dt, time)

      if ( present(ComputeA)) then
         if (ComputeA) then
            call MG_ComputeJacobians( this,no_levels,ComputeTimeDerivative,Time,dt,nEqn )
            ComputeA = .FALSE.
         end if
      else
         call MG_ComputeJacobians( this,no_levels,ComputeTimeDerivative,Time,dt,nEqn )
      end if

      call this % child % p_sem % mesh % storage % local2globalq (this % child % p_sem % mesh % storage % NDOF)

      this % niter = 0
      this % maxiter = maxiter
      this % converged = .false.
      this % tol = tol 

      allocate( this % resvec(maxiter))
      this % resvec = 0._RP

      solver_type = S_SOLVER
      this % x = 0.0_RP ! setting initial solution
      ! this % x = this % p_sem % mesh % storage % Q ! setting initial solution

      select case (solver_type)
         case (SOLVER_PMG_NONE)
            call this % MG_Solve_MB_PMG_NONE(nEqn, nGradEqn, ComputeTimeDerivative, this % tol, this % maxiter, time, dt, ComputeA)
         case (SOLVER_GMRES_PMG)
            call this % MG_KSP_Construct(this % DimPrb, this % maxiter)
            call this % MG_Solve_MB_GMRES_PMG(nEqn, nGradEqn, ComputeTimeDerivative, this % tol, this % maxiter, time, dt, ComputeA)
            call this % MG_KSP_Destruct(this % DimPrb)
      end select 

      print *, "No iterations: ", this % niter-2
      print *, "Residual history: "
      write(*,"(ES14.7)") this % resvec(1:this%niter-1)

      deallocate( this % resvec )
   end subroutine MG_Solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Solve_MB_PMG_NONE(this, nEqn, nGradEqn, ComputeTimeDerivative, tol, maxiter, time, dt, ComputeA)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      integer,       intent(in)                       :: nEqn, nGradEqn
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
      real(kind=rp), optional                         :: tol
      integer      , optional                         :: maxiter
      real(kind=rp), optional                         :: time
      real(kind=rp), optional                         :: dt
      logical      , optional , intent(inout)         :: ComputeA
!-----Local-Variables-----------------------------------------------------
      class(csrMat_t), pointer                        :: pA
      integer                                         :: niter
      real(kind=rp)                                   :: tmpsize
      real(kind=rp)                                   :: shift
      integer                                         :: i, j, k
      character(len=1024)                             :: filename
      logical                                         :: file_exists
!  -----------------------------------------------------------------------

      do while ( (.not. this % converged ) .and. (this % niter .lt. this % maxiter) )

         call MG_VCycle( this, no_levels, nEqn, ComputeTimeDerivative)

         call this % MG_JacVec(this % r,this % x,ComputeTimeDerivative,S_MATCOMP)

         do i = 1 , this % DimPrb
            this % r(i) = this % b(i) - this % r(i)
         end do 

         this % rnorm = norm2(this % r)
         this % resvec(this % niter) = this % rnorm
         this % niter = this % niter + 1
         if (this % rnorm < this % tol) this % converged = .true.

      end do
   end subroutine MG_Solve_MB_PMG_NONE
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Solve_MB_GMRES_PMG(this, nEqn, nGradEqn, ComputeTimeDerivative, tol, maxiter, time, dt, ComputeA)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      !type(LinearMultigridSolver_t) , pointer               :: pMG     
      integer,       intent(in)                       :: nEqn, nGradEqn
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
      real(kind=rp), optional                         :: tol
      integer      , optional                         :: maxiter
      real(kind=rp), optional                         :: time
      real(kind=rp), optional                         :: dt
      logical      , optional , intent(inout)         :: ComputeA
!-----Local-Variables-----------------------------------------------------
      real(kind=RP), dimension(:) ,     allocatable   :: x0,tmp_ilu 
      !class(csrMat_t), pointer                        :: pA
      integer                                         :: niter
      real(kind=rp)                                   :: tmpsize
      character(len=1024)                             :: filename
      logical                                         :: file_exists
      integer                                         :: i,j,k, l, ii,kk, m
      integer                                         :: idx1, idx2
      real(kind = RP)                                 :: tmp1, tmp2
      real(kind = RP)                                 :: rhsnorm 
!  -----------------------------------------------------------------------
     
      allocate (x0(this%DimPrb))
      allocate (tmp_ilu(this%DimPrb))

      associate ( V => this % KSP % V, H => this % KSP % H, &
                  W => this % KSP % W, Y => this % KSP % Y, &
                  cc => this % KSP % cc, ss => this % KSP % ss, &
                  g => this % KSP % g, Z => this % KSP % Z)

      x0 = this % x

      call this % MG_JacVec(V(:,1),this % x,ComputeTimeDerivative,S_MATCOMP)

      V(:,1) = this % b - V(:,1)   
      g(1) = norm2(V(:,1))
      V(:,1) = V(:,1) / g(1)
      rhsnorm = norm2(this % b)

      H = 0.0_RP
      m = this % KSP % KrylovSpace

      do j = 1,m ! Krylov loop
         select case (this % KSP % Preconditioner)
            case (KSP_PRECONDITIONER_NONE)
               call this % MG_JacVec(W,V(:,j),ComputeTimeDerivative,S_MATCOMP)
            case (KSP_PRECONDITIONER_PMG)
               ! z = M^{-1}v
               call this % SetRHS (V(:,j))
               this % x = 0.0_RP
               call MG_VCycle( this, no_levels, nEqn, ComputeTimeDerivative)
               ! call MG_FMGCycle( this, no_levels, nEqn, ComputeTimeDerivative)
               Z(:,j) = this % x
               ! w = Az
               call this % MG_JacVec(W,Z(:,j),ComputeTimeDerivative,S_MATCOMP)

            case (KSP_PRECONDITIONER_BJ)
               ! z = M^{-1}v
               do i=1, this%BJSmoother % A_p % num_of_Blocks
                  idx1 = this%BJSmoother % A_p % BlockIdx(i)
                  idx2 = this%BJSmoother % A_p % BlockIdx(i+1)-1
                  call SolveLU(ALU      = this%BJSmoother % BlockPrec(i) % PLU, &
                               LUpivots = this%BJSmoother % BlockPrec(i) % LUpivots, &
                               x = Z(idx1:idx2,j), &
                               b = V(idx1:idx2,j) )
               end do
               ! w = Az
               call this % MG_JacVec(W,Z(:,j),ComputeTimeDerivative,S_MATCOMP)
            case (KSP_PRECONDITIONER_ILU)
               ! z = M^{-1}v
               call this % ILUsmoother % A % ForwSub(V(:,j) , tmp_ilu, this%DimPrb)
               call this % ILUsmoother % A % BackSub(tmp_ilu, Z(:,j) , this%DimPrb)
               ! w = Az
               call this % MG_JacVec(W, Z(:,j), ComputeTimeDerivative, S_MATCOMP)
            case (KSP_PRECONDITIONER_PJ)
               select case (S_MATCOMP)
                  case (JACOBIANCOMP_MB)
                     do i=1, size(this % b, 1)
                        Z(i,j) = 2._RP/3._RP * V(i,j) / this % A % Values(this % A % Diag(i))
                     end do
                  case (JACOBIANCOMP_MF)
                     ERROR STOP "Matrix-free not implemented for Point-Jacobi prec."
                  case default  
               end select
               call this % MG_JacVec(W,Z(:,j),ComputeTimeDerivative,S_MATCOMP)
            case default ! PC_NONE
         end select
         
         do i = 1,j
            H(i,j) = dot_product(W,V(:,i))
            W = W - H(i,j) * V(:,i)
         end do
         
         H(j+1,j) = norm2(W)
         
         if ((ABS(H(j+1,j)) .LT. this%tol)) then
            this%CONVERGED = .TRUE.
            this%rnorm = this%tol
            this%niter = this%niter + 1                 
            m = j
            exit
         end if

         V(:,j+1) =  W / H(j+1,j) 
         do i = 1, j-1
            tmp1 = H(i,j)
            tmp2 = H(i+1,j)
            H(i,j) = cc(i) * tmp1 + ss(i) * tmp2
            H(i+1,j) = cc(i) * tmp2 - ss(i) * tmp1
         end do 
         
         tmp1 = SQRT(H(j,j)*H(j,j) + H(j+1,j)*H(j+1,j) )
         
         if (ABS(tmp1) .LT. 1e-15_RP) then
            ERROR stop "GMRES Loop has ~0 vec"
            RETURN
         end if

         cc(j) = H(j,j) / tmp1
         ss(j) = H(j+1,j) / tmp1
         g(j+1) = -ss(j) * g(j)
         g(j) = cc(j) * g(j)
         H(j,j) = cc(j) * H(j,j) + ss(j) * H(j+1,j) 
         this%rnorm = ABS(g(j+1)) / rhsnorm
         this % resvec(j) = this % rnorm
         this%niter = this%niter + 1

         if (this%rnorm .LT. this%tol) then
            this%CONVERGED = .TRUE.
            m = j
            exit
         end if
         if (this%niter .GE. this%maxiter) then
            m = j
            this%CONVERGED = .FALSE.
            exit
         end if
      end do ! End of Krylov loop
      
      if (m > 0) then
         y(m) = g(m) / H(m,m)
         do ii = 1, m-1
            kk = m - ii
            tmp1 = g(kk)
            l = kk+1
            do while (l .LE. m)
               tmp1 = tmp1 - H(kk,l) * y(l)
               y(kk) = tmp1 / H(kk,kk)
               l = l+1
            end do 
         end do
      end if ! m > 0
      
      select case (this % KSP % Preconditioner)
         case (KSP_PRECONDITIONER_NONE)
            this%x = x0 + MATMUL(V(:,1:m),y(1:m))
         case (KSP_PRECONDITIONER_BJ) 
            this%x = x0 + MATMUL(Z(:,1:m),y(1:m))
         case default 
            this%x = x0 + MATMUL(Z(:,1:m),y(1:m))
      end select

      end associate

      deallocate (x0)
      deallocate (tmp_ilu)

   end subroutine MG_Solve_MB_GMRES_PMG
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine MG_VCycle(Me,lvl,nEqn,ComputeTimeDerivative)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      type(LinearMultigridSolver_t), target    :: Me
      type(LinearMultigridSolver_t) , pointer  :: Child_p          ! Pointer to Child
      integer, intent(in)                :: lvl
      integer, intent(in)                :: nEqn
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
!-----Local-Variables-----------------------------------------------------
      integer                           :: i ! counter
!  -----------------------------------------------------------------------

      if (lvl == 1) then
         call MG_smooth(Me,Me % x,Me % b, Me % DimPrb , pre_smooths(lvl), ComputeTimeDerivative)
      else
         Child_p => Me % Child ! set local pointer to a child
         Me % Child % Parent => Me ! set child's parent pointer to this level (Me)

         call MG_smooth(Me,Me % x,Me % b, Me % DimPrb , pre_smooths(lvl), ComputeTimeDerivative)

         call Me % MG_JacVec(Me % r, Me % x,ComputeTimeDerivative,S_MATCOMP)

         do i = 1 , Me % DimPrb
            Me % r(i) = Me % b(i) - Me % r(i)
         end do 

         call MG_3DRestriction( Me, nEqn ,.true.,.true.,.false.) ! Restrict to a coarse grid

         Me % child % b = Me % child % r


         Me % child % x = 0.0_RP
         call MG_VCycle( Me % Child, lvl-1, nEqn ,ComputeTimeDerivative)

         call MG_3DProlongation( Me, nEqn )
         call MG_smooth(Me,Me % x,Me % b, Me % DimPrb , pos_smooths(lvl), ComputeTimeDerivative)
      end if 
   end subroutine MG_VCycle
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine MG_FMGCycle(Me,lvl,nEqn,ComputeTimeDerivative)
!  ---------------------------------------------------------
!  Recursive subroutine to perform a full multigrid cycle.
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      type(LinearMultigridSolver_t), target    :: Me
      type(LinearMultigridSolver_t), pointer   :: Child_p          ! Pointer to Child
      integer, intent(in)                :: lvl
      integer, intent(in)                :: nEqn
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
!-----Local-Variables-----------------------------------------------------
      integer        :: iEl, iEQ             ! Element and equation counters
      integer        :: N1(3), N2(3)
      real(kind=RP)  :: maxResidual(NCONS)   ! Maximum residual in each equation
      integer        :: counter              ! Iteration counter
      character(len=LINE_LENGTH) :: FMGFile
!  -----------------------------------------------------------------------

!
!     ------------------------------------------
!     At the beginning, go to the coarsest level
!        (the initial condition must be passed)
!     ------------------------------------------
!
      if (lvl > 1) then
         call MG_3DRestriction( Me, nEqn ,.true.,.true.,.true.) ! Restrict to a coarse grid
         call MG_FMGCycle( Me % Child, lvl-1, nEqn ,ComputeTimeDerivative)
      end if
!
!     ----------------------
!     Perform a V-Cycle here
!     ----------------------
!
      if (lvl > 1 ) then
         call MG_VCycle( Me , lvl, nEqn ,ComputeTimeDerivative)
      else
         call MG_smooth(Me,Me % x,Me % b, Me % DimPrb , pre_smooths(lvl), ComputeTimeDerivative)
      end if
!
!     --------------------------------------------------
!     If not on finest, Interpolate to next (finer) grid
!     --------------------------------------------------
! 
      if (lvl < no_levels) then
         call MG_3DProlongation( Me % Parent, nEqn )
      end if
   end subroutine MG_FMGCycle
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine MG_ComputeJacobians(Me,lvl,ComputeTimeDerivative,Time,dt,nEqn)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      use DenseMatUtilities
      implicit none
!-----Arguments-----------------------------------------------------------
      type(LinearMultigridSolver_t), target    :: Me
      type(LinearMultigridSolver_t) , pointer  :: Child_p          ! Pointer to Child
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      integer, intent(in)                :: nEqn
      integer, intent(in)                :: lvl
      real(kind=rp), intent(in)          :: dt
      real(kind=RP), intent(in)          :: Time
!-----Local-Variables-----------------------------------------------------
      integer             :: i,j
      character(len=1024) :: filename
      REAL(KIND=RP) :: A
!  -----------------------------------------------------------------------


      select case (S_MATCOMP)
      case (JACOBIANCOMP_MB)
         call Me % Jacobian % Compute (Me % p_sem, nEqn, time, Me % A, ComputeTimeDerivative)
         call Me % SetOperatorDt(dt)
      case (JACOBIANCOMP_MF)
         call Me % Jacobian % Compute (Me % p_sem, nEqn, time, Me % BJSmoother % A_p, ComputeTimeDerivative, BlockDiagonalized=.TRUE.)
         call Me % SetOperatorDt(dt)
         A = matrixshift(dt)
         call Me % BJSmoother % A_p % shift( A )
      end select 


      select case (S_SMOOTHER)
      case (S_POINTJAC)
         error stop "Point Jacobi not implemented."
      case (S_ILU)
         call Me % ILUSmoother % Construct (Me % A) 
         call ComputeILU(Me % ILUsmoother)
      case (S_BLOCKJAC)
         select case (S_MATCOMP)
         case (JACOBIANCOMP_MB)
            call getDBDfromCSR(Me % A,Me % BJSmoother % A_p) ! Get Jacobian diag-blocks
            call ComputeBlockPrec(Me % BJSmoother)
         case (JACOBIANCOMP_MF)
            call ComputeBlockPrec(Me % BJSmoother)
         end select 
      case default
         error stop "No smoother selected."
      end select

      if (lvl > 1) then
         Child_p => Me % Child ! set local pointer to a child
         Me % Child % Parent => Me ! set child's parent pointer to this level (Me)

         call MG_1DRestriction(Me, nEqn, 'soljac')
         call MG_ComputeJacobians(Me % Child,lvl-1,ComputeTimeDerivative, Time, dt, nEqn)

      end if

   end subroutine MG_ComputeJacobians
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_JacVec(this,xout,xin,ComputeTimeDerivative,MatFLAG)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      integer,                          intent(in)    :: MatFLAG 
      real(kind=RP),                    intent(in)    :: xin( this % DimPrb)
      real(kind=RP),                    intent(inout) :: xout( this % DimPrb)
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
!-----Local-Variables-----------------------------------------------------
      real(kind=RP) :: shift 
!  -----------------------------------------------------------------------

      select case (MatFLAG)
      case (JACOBIANCOMP_MB)
         xout = CSR_MatVecMul( this % A, xin ) 
      case (JACOBIANCOMP_MF)
         shift = MatrixShift(this % dt)
         call this % p_sem % mesh % storage % local2GlobalQ(this % p_sem % NDOF)
         this % Ur   = this % p_sem % mesh % storage % Q
         this % F_Ur = MF_p_F (this % p_sem, this % DimPrb, this % Ur, this % timesolve + this % dt,ComputeTimeDerivative) 
         call MF_JacVecMul(this % p_sem, this % DimPrb, this % Ur, this % F_Ur, xin, xout, this % dt, this % timesolve, shift, ComputeTimeDerivative)
      case default
         error stop "Wrong Jacobian computation."
      end select

   end subroutine MG_JacVec
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_PrecVec(this)
!  ---------------------------------------------------------
!  TODO. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
!  -----------------------------------------------------------------------

      error stop "TBD"
   end subroutine MG_precVec
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_smooth(this,x,b,n,SmoothIters,ComputeTimeDerivative)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      use DenseMatUtilities
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      real(kind=rp),    intent(inout), dimension(:)   :: x                     ! solution
      real(kind=rp),    intent(in),    dimension(:)   :: b                     ! Right hand side
      integer,          intent(in)                    :: n                     ! System siz
      integer,          intent(in)                    :: SmoothIters           ! # of iterations
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
!-----Local-Variables-----------------------------------------------------
      real(kind=rp)                           :: tmp(n)            ! tmp solution vector
      integer                                 :: i,j, idx1, idx2 ! Counters
!  -----------------------------------------------------------------------

      select case (S_SMOOTHER)
         case (S_NOTDEF)
            ERROR stop 'GenericSmoother :: Wrong smoother type'
         case (S_POINTJAC)
            call MG_PJsmooth(this,x,b,n, SmoothIters, ComputeTimeDerivative)
         case (S_BLOCKJAC)
            call MG_BJsmooth(this,x,b,n, SmoothIters, ComputeTimeDerivative)
         case (S_ILU)
            call MG_ILUsmooth(this,x,b,n, SmoothIters, ComputeTimeDerivative)
      end select 

   end subroutine MG_smooth
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_BJsmooth(this,x,b,n,SmoothIters,ComputeTimeDerivative)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      use DenseMatUtilities
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      real(kind=rp),    intent(inout), dimension(:)   :: x                     ! solution
      real(kind=rp),    intent(in),    dimension(:)   :: b                     ! Right hand side
      integer,          intent(in)                    :: n                     ! System siz
      integer,          intent(in)                    :: SmoothIters           ! # of iterations
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
!-----Local-Variables-----------------------------------------------------
      real(kind=rp)                           :: tmp(n)            ! tmp solution vector
      integer                                 :: i,j, idx1, idx2 ! Counters
!  -----------------------------------------------------------------------

      do j = 1, SmoothIters
         tmp = 0.0_RP

         call this % MG_JacVec(this % r, x, ComputeTimeDerivative, S_MATCOMP)

         do i=1, this%BJSmoother % A_p % num_of_Blocks
            idx1 = this%BJSmoother % A_p % BlockIdx(i)
            idx2 = this%BJSmoother % A_p % BlockIdx(i+1)-1
            this % r(idx1:idx2) = b(idx1:idx2) - this % r(idx1:idx2)
         end do

         do i=1, this%BJSmoother % A_p % num_of_Blocks
            idx1 = this%BJSmoother % A_p % BlockIdx(i)
            idx2 = this%BJSmoother % A_p % BlockIdx(i+1)-1
            call SolveLU(ALU      = this%BJSmoother % BlockPrec(i) % PLU, &
                         LUpivots = this%BJSmoother % BlockPrec(i) % LUpivots, &
                         x = tmp(idx1:idx2), &
                         b = this % r  (idx1:idx2) )
            x(idx1:idx2) = x(idx1:idx2) + tmp(idx1:idx2)
         end do
      end do

   end subroutine MG_BJsmooth
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_PJsmooth(this,x,b,n,SmoothIters,ComputeTimeDerivative)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      use DenseMatUtilities
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      real(kind=rp),    intent(inout), dimension(:)   :: x                     ! solution
      real(kind=rp),    intent(in),    dimension(:)   :: b                     ! Right hand side
      integer,          intent(in)                    :: n                     ! System siz
      integer,          intent(in)                    :: SmoothIters           ! # of iterations
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
!-----Local-Variables-----------------------------------------------------
      real(kind=rp)                           :: tmp(n)            ! tmp solution vector
      integer                                 :: i,j, idx1, idx2 ! Counters
      real(kind=rp), parameter                :: wpj = 2._RP/3._RP  ! Weight (optimal for FD laplacian... but DGSEM?)
!  -----------------------------------------------------------------------

      do j = 1, SmoothIters
         do i=1,n
            this % r(i) = b(i) - this % r(i)
            x(i) = x(i) + wpj * this % r(i) / this % A % Values(this % A % Diag(i))
         end do
      end do

   end subroutine MG_PJsmooth
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_ILUsmooth(this,x,b,n,SmoothIters,ComputeTimeDerivative)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      use DenseMatUtilities
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      real(kind=rp),    intent(inout), dimension(:)   :: x                     ! solution
      real(kind=rp),    intent(in),    dimension(:)   :: b                     ! Right hand side
      integer,          intent(in)                    :: n                     ! System siz
      integer,          intent(in)                    :: SmoothIters           ! # of iterations
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
!-----Local-Variables-----------------------------------------------------
      real(kind=rp), dimension(:), allocatable :: tmp1, tmp2 ! tmp solution vector
      integer                                  :: i,j, idx1, idx2  ! Counters
!  -----------------------------------------------------------------------

      allocate(tmp1(n))
      allocate(tmp2(n))
      do j = 1, SmoothIters
         tmp1 = 0.0_RP
         call this % MG_JacVec(this % r, x, ComputeTimeDerivative, S_MATCOMP)
         this % r = b - this % r

         call this % ILUsmoother % A % ForwSub(this % r, tmp2, n)
         call this % ILUsmoother % A % BackSub(tmp2    , tmp1, n)
         x = x + tmp1

      end do
      deallocate(tmp1)
      deallocate(tmp2)

   end subroutine MG_ILUsmooth
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_CreateProlongationOperator(this)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      use NodalStorageClass
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
!-----Local-Variables-----------------------------------------------------
      integer              :: x_Norigin ! Destination polynomial order X
      integer              :: x_Ndest   ! Destination polynomial order X
      type(NodalStorage_t) :: x_spAo    ! Origin nodal storage X
      type(NodalStorage_t) :: x_spAd    ! Destination nodal storage X
      integer              :: y_Norigin ! Destination polynomial order Y
      integer              :: y_Ndest   ! Destination polynomial order Y
      type(NodalStorage_t) :: y_spAo    ! Origin nodal storage Y
      type(NodalStorage_t) :: y_spAd    ! Destination nodal storage Y
      integer              :: z_Norigin ! Destination polynomial order Z
      integer              :: z_Ndest   ! Destination polynomial order Z
      type(NodalStorage_t) :: z_spAo    ! Origin nodal storage Z
      type(NodalStorage_t) :: z_spAd    ! Destination nodal storage Z
!  -----------------------------------------------------------------------

      x_Norigin  = this % Nx
      x_Ndest    = this % parent % Nx
      y_Norigin  = this % Ny
      y_Ndest    = this % parent % Ny
      z_Norigin  = this % Nz
      z_Ndest    = this % parent % Nz
      !-----------------------------------------------------
      call x_spAo % construct(this % p_sem % mesh % nodeType, x_Norigin)
      call x_spAd % construct(this % p_sem % mesh % nodeType, x_Ndest  )
      call y_spAo % construct(this % p_sem % mesh % nodeType, y_Norigin)
      call y_spAd % construct(this % p_sem % mesh % nodeType, y_Ndest  )
      call z_spAo % construct(this % p_sem % mesh % nodeType, z_Norigin)
      call z_spAd % construct(this % p_sem % mesh % nodeType, z_Ndest  )


      ! 1D X direction
      !-----------------------------------------------------
      call PolynomialInterpolationMatrix(x_Norigin, x_Ndest, x_spAo % x, x_spAo % wb, x_spAd % x, this % ProlX)
      !-----------------------------------------------------

      ! 1D Y direction
      !-----------------------------------------------------
      call PolynomialInterpolationMatrix(y_Norigin, y_Ndest, y_spAo % x, y_spAo % wb, y_spAd % x, this % ProlY)
      !-----------------------------------------------------

      ! 1D Z direction
      !-----------------------------------------------------
      call PolynomialInterpolationMatrix(z_Norigin, z_Ndest, z_spAo % x, z_spAo % wb, z_spAd % x, this % ProlZ)
      !-----------------------------------------------------

      ! 3D
      !-----------------------------------------------------
      call MG_Create3DProlongationMatrix(this % Prol3D, x_Norigin, y_Norigin, z_Norigin, x_Ndest, y_Ndest, z_Ndest, &
       x_spAo % x, y_spAo % x, z_spAo % x , x_spAd % x, y_spAd % x, z_spAd % x)
      !-----------------------------------------------------

      call x_spAo % destruct
      call x_spAd % destruct
      call y_spAo % destruct
      call y_spAd % destruct
      call z_spAo % destruct
      call z_spAd % destruct
   end subroutine MG_CreateProlongationOperator
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_CreateRestrictionOperator(this,Rweighted)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      use NodalStorageClass
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      logical, optional, intent(in) :: Rweighted    ! flag for weighted restriction operator 
!-----Local-Variables-----------------------------------------------------
      integer              :: x_Norigin ! x Destination polynomial order
      integer              :: x_Ndest   ! x Destination polynomial order
      type(NodalStorage_t) :: x_spAo    ! x Origin nodal storage
      type(NodalStorage_t) :: x_spAd    ! x Destination nodal storage
      integer              :: y_Norigin ! y Destination polynomial order
      integer              :: y_Ndest   ! y Destination polynomial order
      type(NodalStorage_t) :: y_spAo    ! y Origin nodal storage
      type(NodalStorage_t) :: y_spAd    ! y Destination nodal storage
      integer              :: z_Norigin ! z Destination polynomial order
      integer              :: z_Ndest   ! z Destination polynomial order
      type(NodalStorage_t) :: z_spAo    ! z Origin nodal storage
      type(NodalStorage_t) :: z_spAd    ! z Destination nodal storage
      ! only for restriction
      real(kind=RP), allocatable    :: Rtmp(:,:)
      integer                       :: i,j
!  -----------------------------------------------------------------------

      x_Norigin  = this % Nx
      x_Ndest    = this % child % Nx
      y_Norigin  = this % Ny
      y_Ndest    = this % child % Ny
      z_Norigin  = this % Nz
      z_Ndest    = this % child % Nz

      !-----------------------------------------------------
      call x_spAo % construct(this % p_sem % mesh % nodeType, x_Norigin)
      call x_spAd % construct(this % p_sem % mesh % nodeType, x_Ndest  )
      call y_spAo % construct(this % p_sem % mesh % nodeType, y_Norigin)
      call y_spAd % construct(this % p_sem % mesh % nodeType, y_Ndest  )
      call z_spAo % construct(this % p_sem % mesh % nodeType, z_Norigin)
      call z_spAd % construct(this % p_sem % mesh % nodeType, z_Ndest  )

      ! 1D X direction
      !-----------------------------------------------------
      allocate (Rtmp(0:x_Norigin,0:x_Ndest))
      call PolynomialInterpolationMatrix(x_Ndest, x_Norigin, x_spAd % x, x_spAd % wb, x_spAo % x, Rtmp)
      this % RestX = transpose(Rtmp)
      deallocate (Rtmp)
      this % RestJacX = this % RestX 
      if (present(Rweighted)) then
         if (Rweighted) then
            do j = 0, x_Norigin ; do i = 0, x_Ndest
               this % RestJacX(i,j) = this % RestJacX(i,j) * x_spAo % w(j) / x_spAd % w(i)
            end do            ; end do
         end if
      end if
      !-----------------------------------------------------

      ! 1D Y direction
      !-----------------------------------------------------
      allocate (Rtmp(0:y_Norigin,0:y_Ndest))
      call PolynomialInterpolationMatrix(y_Ndest, y_Norigin, y_spAd % x, y_spAd % wb, y_spAo % x, Rtmp)
      this % RestY = transpose(Rtmp)
      deallocate (Rtmp)
      this % RestJacY = this % RestY
      if (present(Rweighted)) then
         if (Rweighted) then
            do j = 0, y_Norigin ; do i = 0, y_Ndest
               this % RestJacY(i,j) = this % RestJacY(i,j) * y_spAo % w(j) / y_spAd % w(i)
            end do            ; end do
         end if
      end if
      !-----------------------------------------------------

      ! 1D Z direction
      !-----------------------------------------------------
      allocate (Rtmp(0:z_Norigin,0:z_Ndest))
      call PolynomialInterpolationMatrix(z_Ndest, z_Norigin, z_spAd % x, z_spAd % wb, z_spAo % x, Rtmp)
      this % RestZ = transpose(Rtmp)
      deallocate (Rtmp)
      this % RestJacZ = this % RestZ
      if (present(Rweighted)) then
         if (Rweighted) then
            do j = 0, z_Norigin ; do i = 0, z_Ndest
               this % RestJacZ(i,j) = this % RestJacZ(i,j) * z_spAo % w(j) / z_spAd % w(i)
            end do            ; end do
         end if
      end if
      !-----------------------------------------------------

      ! 3D
      !-----------------------------------------------------
      ! call MG_Create3DRestrictionMatrix(this % Rest3D, x_Norigin, y_Norigin, z_Norigin, x_Ndest, y_Ndest, z_Ndest, &
      !  x_spAo % x, y_spAo % x, z_spAo % x , x_spAd % x, y_spAd % x, z_spAd % x, &
      !  x_spAo % w, y_spAo % w, z_spAo % w , x_spAd % w, y_spAd % w, z_spAd % w, Rweighted)
      call MG_Create3DRestrictionMatrix(this % Rest3D, x_Norigin, y_Norigin, z_Norigin, x_Ndest, y_Ndest, z_Ndest, &
       x_spAo % x, y_spAo % x, z_spAo % x , x_spAd % x, y_spAd % x, z_spAd % x, &
       x_spAo % w, y_spAo % w, z_spAo % w , x_spAd % w, y_spAd % w, z_spAd % w, Rweighted)
      !-----------------------------------------------------

      call x_spAo % destruct
      call x_spAd % destruct
      call y_spAo % destruct
      call y_spAd % destruct
      call z_spAo % destruct
      call z_spAd % destruct

   end subroutine MG_CreateRestrictionOperator
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_1DProlongation(this,nEqn)
!  ---------------------------------------------------------
!  Prolong a 3D array using 1D interpolation operators.
!  NOT USED IN THE CURRENT VERSION. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      integer                         , intent(in)    :: nEqn
!-----Local-Variables-----------------------------------------------------
      class(LinearMultigridSolver_t), pointer               :: child_p
      integer                                         :: iEl
      integer, dimension(3)                           :: Norigin
      integer, dimension(3)                           :: Ndest
      integer, dimension(nelem)                       :: elsDOF
      !-----------------------------------------------------
      real(kind=RP), dimension(:), allocatable                                          :: x_org
      real(kind=RP), dimension(nEqn,0:this%Nx, 0:this%Ny, 0:this%Nz)                    :: xe_f
      real(kind=RP), dimension(nEqn,0:this%child%Nx, 0:this%child%Ny, 0:this%child%Nz)  :: xe_c
!  -----------------------------------------------------------------------

      child_p => this % child
      Norigin = (/ child_p % Nx, child_p % Ny, child_p % Nz /)
      Ndest = (/ this % Nx, this % Ny, this % Nz /)

      allocate ( x_org(this % DimPrb) ) 
      x_org = this % x 
      xe_f = 0.0_RP

      ! Vec to Elements
      !-----------------------------------------------------
      do iEl = 1, size(child_p % p_sem % mesh % elements)
         elsDOF(iEL) = child_p % p_sem % mesh  % elements(iEL) % storage % NDOF
      end do
      call MG_Vec2El ( child_p % x, child_p % LocalStorage , Norigin , child_p % DimPrb, nEqn, elsDOF , 'x' )
      !-----------------------------------------------------

      do iEl = 1, size(this % p_sem % mesh % elements)
         call MG_Interp3DArrays(nEqn, Norigin, child_p % LocalStorage(iEl) % x , &
                                      Ndest, this % LocalStorage(iEl) % x , &
                                      child_p % ProlX, child_p % ProlY, child_p % ProlZ )
      end do

      ! Elements to Vec
      !-----------------------------------------------------
      do iEl = 1, size(this % p_sem % mesh % elements)
         elsDOF(iEL) = this % p_sem % mesh  % elements(iEL) % storage % NDOF
      end do
      call MG_El2Vec ( this % LocalStorage , this % x, Ndest , this % DimPrb, nEqn, elsDOF , 'x' )
      !-----------------------------------------------------

      this % x = x_org + this % x
      deallocate( x_org )

   end subroutine MG_1DProlongation
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_1DRestriction(this, nEqn, r_vars)
!  ---------------------------------------------------------
!  Restrict a 3D array using 1D interpolation operators.
!  NOT USED IN THE CURRENT VERSION. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      integer                         , intent(in)    :: nEqn
      character(len=*)                , intent(in)    :: r_vars
!-----Local-Variables-----------------------------------------------------
      class(LinearMultigridSolver_t), pointer               :: child_p
      integer                                         :: iEl
      integer, dimension(3)                           :: Norigin
      integer, dimension(3)                           :: Ndest
      integer, dimension(nelem)                       :: elsDOF
      integer                                         :: i,j,k ! counters
      !-----------------------------------------------------
      real(kind=RP), dimension(nEqn,0:this%Nx, 0:this%Ny, 0:this%Nz)  :: re_f
      real(kind=RP), dimension(nEqn,0:this%Nx, 0:this%Ny, 0:this%Nz)  :: xe_f
      real(kind=RP), dimension(nEqn,0:this%child%Nx, 0:this%child%Ny, 0:this%child%Nz)  :: re_c
      real(kind=RP), dimension(nEqn,0:this%child%Nx, 0:this%child%Ny, 0:this%child%Nz)  :: xe_c
!  -----------------------------------------------------------------------

      child_p => this % child
      Norigin = (/ this % Nx, this % Ny, this % Nz /)
      Ndest = (/ child_p % Nx, child_p % Ny, child_p % Nz /)

      select case(r_vars)

         case('solres')

            ! Vec to Elements
            !-----------------------------------------------------
            do iEl = 1, size(this % p_sem % mesh % elements)
               elsDOF(iEL) = this % p_sem % mesh  % elements(iEL) % storage % NDOF
            end do
            call MG_Vec2El ( this % x, this % LocalStorage , Norigin , this % DimPrb, nEqn, elsDOF , 'x' )
            call MG_Vec2El ( this % r, this % LocalStorage , Norigin , this % DimPrb, nEqn, elsDOF , 'r' )
            !-----------------------------------------------------

            do iEl = 1, size(this % p_sem % mesh % elements)
               call MG_Interp3DArrays(nEqn, Norigin, this % LocalStorage(iEL) % x , &
                                            Ndest, child_p % LocalStorage(iEl) % x , &
                                            this % RestJacX, this % RestJacY, this % RestJacZ )
                                            !this % RestX, this % RestY, this % RestZ )

               call MG_Interp3DArrays(nEqn, Norigin, this % LocalStorage(iEl) % r , &
                                            Ndest, child_p % LocalStorage(iEl) % r, &
                                            this % RestJacX, this % RestJacY, this % RestJacZ )
                                            !this % RestX, this % RestY, this % RestZ )
            end do

            ! Elements to Vec
            !-----------------------------------------------------
            do iEl = 1, size(child_p % p_sem % mesh % elements)
               elsDOF(iEL) = child_p % p_sem % mesh  % elements(iEL) % storage % NDOF
            end do
            call MG_El2Vec ( child_p % LocalStorage , child_p % x, Ndest , this % DimPrb, nEqn, elsDOF , 'x' )
            call MG_El2Vec ( child_p % LocalStorage , child_p % r, Ndest , this % DimPrb, nEqn, elsDOF , 'r')
            !-----------------------------------------------------

         case('soljac')

            do iEl = 1, size(this % p_sem % mesh % elements)
               call MG_Interp3DArrays(nEqn, Norigin, this % p_sem % mesh % elements(iEl) % storage % Q, &
                                            Ndest, child_p % p_sem % mesh % elements(iEl) % storage % Q, &
                                            this % RestJacX, this % RestJacY, this % RestJacZ )
                                            !this % RestX, this % RestY, this % RestZ )
            end do

      end select

   end subroutine MG_1DRestriction
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_3DProlongation(this,nEqn)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      real(kind=RP), dimension(:,:), pointer          :: P_p
      integer                         , intent(in)    :: nEqn
!-----Local-Variables-----------------------------------------------------
      class(LinearMultigridSolver_t), pointer               :: child_p
      integer                                         :: iEl
      integer, dimension(3)                           :: Norigin
      integer, dimension(3)                           :: Ndest
      integer, dimension(nelem)                       :: elsDOF
      !-----------------------------------------------------
      integer                                         :: i,j,k,l
      integer                                         :: size_el_f, size_el_c ! size of the element excluding for Neqn=1
      real(kind=RP), dimension(:), allocatable        :: x_org
!  -----------------------------------------------------------------------

      child_p => this % child
      P_p => child_p % Prol3D 

      Norigin = (/ child_p % Nx, child_p % Ny, child_p % Nz /)
      Ndest = (/ this % Nx, this % Ny, this % Nz /)

      size_el_f = (Ndest(1)+1) * (Ndest(2)+1) * (Ndest(3)+1)
      size_el_c = (Norigin(1)+1) * (Norigin(2)+1) * (Norigin(3)+1)

      x_org = this % x 

      do iEl = 1, size(child_p % p_sem % mesh % elements) ! loop over the elements
         l = 1
         do i = 1, size_el_f*Neqn, Neqn ! loop over the fine grid vector size 
            j = i + size_el_f*Neqn * (iEl-1)
            do k = 0, (Neqn - 1)
               this % x (j+k) = this % x (j+k) + & 
                  dot_product(P_p(l,:),child_p % x ( (1 + (size_el_c*Neqn)*(iEl-1) + k) : ((size_el_c*Neqn)*iEl + k) : Neqn ))
            end do
            l = l + 1
         end do
      end do
   end subroutine MG_3DProlongation
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_3DRestriction(this, nEqn , do_x, do_r, do_b)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      integer                         , intent(in)    :: nEqn
      logical                         , intent(in)    :: do_x
      logical                         , intent(in)    :: do_r
      logical                         , intent(in)    :: do_b
!-----Local-Variables-----------------------------------------------------
      class(LinearMultigridSolver_t), pointer               :: child_p
      real(kind=RP), dimension(:,:), pointer          :: R_p
      integer                                         :: iEl
      integer, dimension(3)                           :: Norigin
      integer, dimension(3)                           :: Ndest
      integer, dimension(nelem)                       :: elsDOF
      integer                                         :: i,j,k,l ! counters
      !-----------------------------------------------------
      integer                                         :: size_el_f, size_el_c ! size of the element excluding for Neqn=1
      real(kind=RP), dimension(:), allocatable        :: x_org
!  -----------------------------------------------------------------------

      child_p => this % child
      R_p => this % Rest3D 
      Norigin = (/ this % Nx, this % Ny, this % Nz /)
      Ndest = (/ child_p % Nx, child_p % Ny, child_p % Nz /)

      size_el_c = (Ndest(1)+1) * (Ndest(2)+1) * (Ndest(3)+1)
      size_el_f = (Norigin(1)+1) * (Norigin(2)+1) * (Norigin(3)+1)

      if (do_x) then
         do iEl = 1, size(child_p % p_sem % mesh % elements) ! loop over the elements
            l = 1
            do i = 1, size_el_c*Neqn, Neqn ! loop over the fine grid vector size 
               j = i + size_el_c*Neqn * (iEl-1)
               do k = 0, (Neqn - 1)
                  child_p % x (j+k) = dot_product(R_p(l,:),this % x ( (1 + (size_el_f*Neqn)*(iEl-1) + k) : ((size_el_f*Neqn)*iEl + k) : Neqn ))
                  ! print *, j+k, l, 1 + (size_el_c*Neqn)*(iEl-1) + k
               end do
               l = l + 1
            end do
         end do
      end if

      if (do_r) then
         do iEl = 1, size(child_p % p_sem % mesh % elements) ! loop over the elements
            l = 1
            do i = 1, size_el_c*Neqn, Neqn ! loop over the fine grid vector size 
               j = i + size_el_c*Neqn * (iEl-1)
               do k = 0, (Neqn - 1)
                  child_p % r (j+k) = dot_product(R_p(l,:),this % r ( (1 + (size_el_f*Neqn)*(iEl-1) + k) : ((size_el_f*Neqn)*iEl + k) : Neqn ))
                  ! print *, j+k, l, 1 + (size_el_c*Neqn)*(iEl-1) + k
               end do
               l = l + 1
            end do
         end do
      end if

      if (do_b) then
         do iEl = 1, size(child_p % p_sem % mesh % elements) ! loop over the elements
            l = 1
            do i = 1, size_el_c*Neqn, Neqn ! loop over the fine grid vector size 
               j = i + size_el_c*Neqn * (iEl-1)
               do k = 0, (Neqn - 1)
                  child_p % b (j+k) = dot_product(R_p(l,:),this % b ( (1 + (size_el_f*Neqn)*(iEl-1) + k) : ((size_el_f*Neqn)*iEl + k) : Neqn ))
               end do
               l = l + 1
            end do
         end do
      end if

   end subroutine MG_3DRestriction
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_JacRestriction(this, nEqn )
!  ---------------------------------------------------------
!  Perform R * A * P operation, where:
!
!  A is a (NDOF x NDOF) CSR matrix,
!  R is a (Nc x Nf) full matrix,  
!  P is a (Nf x Nc) full matrix.
!  
!  Nc = (N_x_coarse + 1) * (N_y_coarse + 1) * (N_z_coarse + 1)
!  Nf = (N_x_fine + 1) * (N_y_fine + 1) * (N_z_fine + 1)
!  NDOF = Nel * Neqn * Nf
!  
!  N_x_coarse/fine is a pol. order in x direction.
!  Nel is a number of elements.
!  Neqn is a number of solved equations. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      integer                         , intent(in)    :: nEqn
!-----Local-Variables-----------------------------------------------------
      class(LinearMultigridSolver_t), pointer         :: child_p
      real(kind=RP), dimension(:,:), pointer          :: R_p
      real(kind=RP), dimension(:,:), pointer          :: P_p
      integer                                         :: iEl_row
      integer                                         :: iEl_col
      integer, dimension(3)                           :: Norigin
      integer, dimension(3)                           :: Ndest
      integer, dimension(nelem)                       :: elsDOF
      integer                                         :: i,j,k,l ! counters
      !-----------------------------------------------------
      integer                                         :: size_el_f, size_el_c 
      integer                                         :: size_el_eqn_f, size_el_eqn_c 
      integer                                         :: bStart_row, bEnd_row, bStart_col, bEnd_col
      integer                                         :: icol
      real(kind=RP), dimension(:), allocatable        :: tmp1
      real(kind=RP), dimension(:), allocatable        :: tmp2
      real(kind=RP), dimension(:), allocatable        :: tmp3
      real(kind=RP), dimension(:), allocatable        :: tmp4
      real(kind=RP), dimension(:,:), allocatable      :: tmp_mat
!  -----------------------------------------------------------------------


      ! set pointers for convenience
      child_p => this    % child
      R_p     => this    % Rest3D 
      P_p     => child_p % Prol3D

      ! array of coarse and fine pol orders
      Norigin = (/ this % Nx, this % Ny, this % Nz /)
      Ndest = (/ child_p % Nx, child_p % Ny, child_p % Nz /)

      ! size of an element for one equation
      size_el_c = (Ndest(1)+1) * (Ndest(2)+1) * (Ndest(3)+1)
      size_el_f = (Norigin(1)+1) * (Norigin(2)+1) * (Norigin(3)+1)

      ! total size of a block-element in the jacobian matrix A  
      size_el_eqn_c = size_el_c * Neqn
      size_el_eqn_f = size_el_f * Neqn

      ! allocate 
      allocate( tmp1(size_el_eqn_f) )
      allocate( tmp2(size_el_eqn_f) )
      allocate( tmp3(size_el_eqn_f) )
      allocate( tmp4(size_el_eqn_f) )
      allocate( tmp_mat(size_el_eqn_f,size_el_eqn_f) )

      tmp1 = 0.0_RP ! a single column from the CSR Jacobian matrix
      tmp_mat = 0.0_RP ! local dense matrix

      ! ( this%Cols(k) .gt. (i) * Adbd % BlockSizes(i)) ! check if we exceed the diagonal block
      ! ( this%Cols(k) .gt. (i-1) * Adbd % BlockSizes(i) )  ! check if we are in the right diagonal block

      do iEl_row = 1, size(child_p % p_sem % mesh % elements) ! loop over the elements (rows)
         do iEl_col = 1, size(child_p % p_sem % mesh % elements) ! loop over the elements (rows)

            ! find cols and rows range within the block
            bstart_row = 1 + size_el_eqn_f * (iEl_row - 1)
            bend_row = size_el_eqn_f * iEl_row

            bstart_col = 1 + size_el_eqn_f * (iEl_col - 1)
            bend_col = size_el_eqn_f * iEl_col

         end do
      end do

      ! allocate 
      deallocate( tmp1 )
      deallocate( tmp2 )
      deallocate( tmp3 )
      deallocate( tmp4 )
      deallocate( tmp_mat )

      error stop "TO BE FINISHED"

   end subroutine MG_JacRestriction
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   subroutine MG_Interp3DArrays(Nvars, Nin, inArray, Nout, outArray, InterpX, InterpY, InterpZ )
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      integer                                                        , intent(in)    :: Nvars
      integer      , dimension(3)                                    , intent(in)    :: Nin
      integer      , dimension(3)                                    , intent(in)    :: Nout
      real(kind=RP), dimension(Nvars,0:Nin (1), 0:Nin (2), 0:Nin (3)), intent(in)    :: inArray
      real(kind=RP), dimension(Nvars,0:Nout(1), 0:Nout(2), 0:Nout(3)), intent(out)   :: outArray
      real(kind=RP), dimension(0:Nout(1), 0:Nin(1))                  , intent(in)    :: InterpX
      real(kind=RP), dimension(0:Nout(2), 0:Nin(2))                  , intent(in)    :: InterpY
      real(kind=RP), dimension(0:Nout(3), 0:Nin(3))                  , intent(in)    :: InterpZ
!-----Local-Variables-----------------------------------------------------
      integer :: i,j,k,l,m,n
!  -----------------------------------------------------------------------
      
      outArray = 0.0_RP
      
      do n = 0, Nin(3)  ; do k = 0, Nout(3)
         do m = 0, Nin(2)  ; do j = 0, Nout(2)   
            do l = 0, Nin(1)  ; do i = 0, Nout(1)
               outArray(:,i,j,k) = outArray(:,i,j,k) +   InterpX (i,l) &
                                                       * InterpY (j,m) &
                                                       * InterpZ (k,n) &
                                                       * inArray(:,l,m,n)
            end do             ; end do
         end do             ; end do
      end do             ; end do

   end subroutine MG_Interp3DArrays
!
!////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Create3DProlongationMatrix(Mat,N1x,N1y,N1z,N2x,N2y,N2z,x1,y1,z1,x2,y2,z2)
!  ---------------------------------------------------------
!  Creates a 3D Lagrange interpolation matrix from a grid with 
!  coordinates x1, y1, z1 (origin) to a grid with coordinates
!  x2, y2, z2 (destination).
!  Author (original): David Kopriva 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      real(kind=rp)               ,intent(inout) :: Mat((N2x + 1) * (N2y + 1) * (N2z + 1),(N1x + 1) * (N1y + 1) * (N1z + 1)) !<>  
      integer                     ,intent(in)    :: N1x,N1y,N1z  !<  Origin order
      integer                     ,intent(in)    :: N2x,N2y,N2z  !<  Destination order
      real(kind=rp), dimension(:) ,intent(in)    :: x1(0:N1x),y1(0:N1y),z1(0:N1z)     !<  Nodes in origin
      real(kind=rp), dimension(:) ,intent(in)    :: x2(0:N2x),y2(0:N2y),z2(0:N2z)     !<  Nodes in destination
!-----Local-Variables-----------------------------------------------------
      integer :: i,j,k,l,m,n      ! Coordinate counters
      integer :: r,s              ! Matrix index counters
      ! integer :: NDOFEL1, NDOFEL2 ! Degrees of freedom in origin and destination
!  -----------------------------------------------------------------------

      do k=0, N1z
         do j=0, N1y
            do i=0, N1x
               r = i + j*(N1x + 1) + k*(N1x + 1)*(N1y + 1) + 1            ! Column index
               do n=0, N2z
                  do m=0, N2y
                     do l=0, N2x
                        s = l + m*(N2x + 1) + n*(N2x + 1)*(N2y + 1) + 1   ! Row index
                        Mat(s,r) =  MG_LagrangeInterpolationNoBar(x2(l),N1x,x1,i) * &
                                    MG_LagrangeInterpolationNoBar(y2(m),N1y,y1,j) * &
                                    MG_LagrangeInterpolationNoBar(z2(n),N1z,z1,k)
                     end do
                  end do
               end do
            end do
         end do
      end do

   end subroutine MG_Create3DProlongationMatrix
!
!////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Create3DRestrictionMatrix(Mat,N1x,N1y,N1z,N2x,N2y,N2z,x1,y1,z1,x2,y2,z2,w1x,w1y,w1z,w2x,w2y,w2z,Rweighted)
!  ---------------------------------------------------------
!  Creates an L2-3D Lagrange interpolation matrix from a grid  
!  with coordinates x1, y1, z1 (origin) to a grid with coordinates
!  x2, y2, z2 (destination).
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      ! real(kind=rp)               ,intent(inout) :: Mat((N1x + 1) * (N1y + 1) * (N1z + 1),(N2x + 1) * (N2y + 1) * (N2z + 1)) !<>
      real(kind=rp)               ,intent(inout) :: Mat((N2x + 1) * (N2y + 1) * (N2z + 1),(N1x + 1) * (N1y + 1) * (N1z + 1)) !<>  
      integer                     ,intent(in)    :: N1x,N1y,N1z  !<  Origin order
      integer                     ,intent(in)    :: N2x,N2y,N2z  !<  Destination order
      real(kind=rp), dimension(:) ,intent(in)    :: x1 (0:N1x),y1 (0:N1y),z1 (0:N1z)     !<  Nodes in origin
      real(kind=rp), dimension(:) ,intent(in)    :: x2 (0:N2x),y2 (0:N2y),z2 (0:N2z)     !<  Nodes in destination
      real(kind=rp), dimension(:) ,intent(in)    :: w1x(0:N1x),w1y(0:N1y),w1z(0:N1z)     !<  Weights in origin
      real(kind=rp), dimension(:) ,intent(in)    :: w2x(0:N2x),w2y(0:N2y),w2z(0:N2z)     !<  Weights in destination
      logical                     ,intent(in)    :: Rweighted !<  flag
!-----Local-Variables-----------------------------------------------------
      integer       :: i,j,k,l,m,n      ! Coordinate counters
      integer       :: r,s              ! Matrix index counters
      real(kind=rp) :: MASSterm         ! 
!  -----------------------------------------------------------------------
      
      do k=0, N1z
         do j=0, N1y
            do i=0, N1x
               r = i + j*(N1x + 1) + k*(N1x + 1)*(N1y + 1) + 1            ! Column index
               do n=0, N2z
                  do m=0, N2y
                     do l=0, N2x
                        s = l + m*(N2x + 1) + n*(N2x + 1)*(N2y + 1) + 1   ! Row index
                        
                        Mat(s,r) = MG_LagrangeInterpolationNoBar(x1(i),N2x,x2,l) * &
                                   MG_LagrangeInterpolationNoBar(y1(j),N2y,y2,m) * &
                                   MG_LagrangeInterpolationNoBar(z1(k),N2z,z2,n) 
                        if ( Rweighted ) Mat(s,r) = Mat(s,r) * w1x(i) * w1y(j) * w1z(k) 

                     end do
                  end do
               end do
            end do
         end do
      end do
      
      if ( Rweighted ) then
         do n=0, N2z
            do m=0, N2y
               do l=0, N2x
                  s = l + m*(N2x + 1) + n*(N2x + 1)*(N2y + 1) + 1   ! Row index
                  MASSterm = w2x(l) * w2y(m) * w2z(n)
                  Mat(s,:) = Mat(s,:) / MASSterm
               end do
            end do
         end do
      end if
      
   end subroutine MG_Create3DRestrictionMatrix
!
!////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Create1elCSRInterpolationMats(this, nEqn)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), target, intent(inout) :: this
      integer                         , intent(in)    :: nEqn
!-----Local-Variables-----------------------------------------------------
      class(LinearMultigridSolver_t), pointer               :: child_p
      real(kind=RP), dimension(:,:), pointer          :: R_p
      real(kind=RP), dimension(:,:), pointer          :: P_p
      !----------------------------------------------------------
      integer                                         :: size_el_c, size_el_f, size_el_eqn_c, size_el_eqn_f
      integer                                         :: i,j,k,l,m,n      ! Coordinate counters
      integer, dimension(3)                           :: Norigin
      integer, dimension(3)                           :: Ndest
      real(kind=RP), dimension(:,:), allocatable      :: P1el
      real(kind=RP), dimension(:,:), allocatable      :: R1el
      character(len=1024)                             :: filename
!  -----------------------------------------------------------------------

      child_p => this    % child

      ! array of coarse and fine pol orders
      Norigin = (/ this % Nx, this % Ny, this % Nz /)
      Ndest = (/ child_p % Nx, child_p % Ny, child_p % Nz /)

      ! size of an element for one equation
      size_el_c = (Ndest(1)+1) * (Ndest(2)+1) * (Ndest(3)+1)
      size_el_f = (Norigin(1)+1) * (Norigin(2)+1) * (Norigin(3)+1)

      ! total size of a block-element in the jacobian matrix A  
      size_el_eqn_c = size_el_c * Neqn
      size_el_eqn_f = size_el_f * Neqn


      ! first, construct 1el dense matrices 
      allocate( P1el(size_el_eqn_f,size_el_eqn_c) )
      allocate( R1el(size_el_eqn_c,size_el_eqn_f) )



      print *, size(child_p % Prol3D,1), size(child_p % Prol3D,2)
      print *, size(P1el,1), size(P1el,2)

      do i = 1, size_el_f
         do j = 0, Neqn-1
            P1el(1 + (i-1)*Neqn + j,1 + j*size_el_c:size_el_c + j*size_el_c) = child_p % Prol3D(i,:)
         end do
      end do

      do i = 1, size_el_c
         do j = 0, Neqn-1
            R1el(1 + (i-1)*Neqn + j,1 + j*size_el_f:size_el_f + j*size_el_f) = this % Rest3D(i,:)
         end do
      end do

      allocate( this % Rfull(nelem * size_el_eqn_c, nelem * size_el_eqn_f) )
      allocate( child_p % Pfull(nelem * size_el_eqn_f, nelem * size_el_eqn_c) )

      this % Rfull = 0.0_RP
      child_p % Pfull = 0.0_RP

      print *, nelem

      do i = 1, nelem
         this % Rfull ((1 + (i-1)*size_el_eqn_c):(i*size_el_eqn_c) , (1 + (i-1)*size_el_eqn_f):(i*size_el_eqn_f)) = R1el
         child_p % Pfull ((1 + (i-1)*size_el_eqn_f):(i*size_el_eqn_f) , (1 + (i-1)*size_el_eqn_c):(i*size_el_eqn_c)) = P1el
      end do 

      deallocate( P1el )
      deallocate( R1el )
   end subroutine MG_Create1elCSRInterpolationMats 
!
!/////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------
!     Compute the value of the interpolant WITHOUT using the barycentric formula.
!     --------------------------------------------------------------------------
!
   function MG_LagrangeInterpolationNoBar( x, N, nodes, j) RESULT(l)
!  ---------------------------------------------------------
!  Constructor. 
!  Author (original): David Kopriva 
!  ---------------------------------------------------------
      use Utilities, only: almostEqual
!-----Input---------------------------------------------------------------
      REAL(KIND=RP)                 :: x     !<  Point of evaluation of interpolant
      INTEGER                       :: N     !<  Polynomial order  
      REAL(KIND=RP), DIMENSION(0:N) :: nodes !<  Nodes of Lagrange interpolation
      INTEGER                       :: j     !<  Index of polynomial to be found
!-----Output--------------------------------------------------------------
      REAL(KIND=RP)                 :: l     !>  Lagrange interpolant
!-----Local-Variables-----------------------------------------------------
      INTEGER                       :: i
      REAL(KIND=RP)                 :: numerator, denominator
      REAL(KIND=RP), DIMENSION(0:N) :: values
!  -----------------------------------------------------------------------
      
      values      = 0.0_RP
      values(j)   = 1.0_RP
      numerator   = 1.0_RP
      denominator = 1.0_RP

      DO i = 0, N
         IF( AlmostEqual( x, nodes(i) ) )    THEN
            l = values(i)
            RETURN 
         ELSE IF (j.ne.i) THEN
         numerator   = numerator*(x - nodes(i))    
         denominator = denominator*(nodes(j) - nodes(i))
         END IF 
      END DO
      l = numerator/denominator

   END FUNCTION MG_LagrangeInterpolationNoBar      
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Vec2El(vec,els,N,vecsize,nEqn,eldof,var2interp)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      integer,       dimension(3),                           intent(in)    :: N
      integer,                                               intent(in)    :: vecsize
      integer,                                               intent(in)    :: nEqn
      integer,       dimension(nelem),                       intent(in)    :: eldof
      real(kind=RP), dimension(vecsize),                     intent(in)    :: vec
      ! real(kind=RP), dimension(nEqn,0:N(1), 0:N(2), 0:N(3)), intent(inout) :: els
      type(TemporaryElementStorage_t) , dimension(nelem)   , intent(inout) :: els
      character(len=*)                                     , intent(in)    :: var2interp
!-----Local-Variables-----------------------------------------------------
      integer :: eID, firstIdx, lastIdx ! counters
!  -----------------------------------------------------------------------
      
      firstIdx = 1
      do eID = 1, nelem
         lastIdx = firstIdx + eldof(eID) * nEqn
         select case (var2interp)
         case ('x')
            els(eID) % x(1:nEqn,0:N(1),0:N(2),0:N(3)) = reshape (vec (firstIdx:lastIdx-1),(/ nEqn, N(1)+1, N(2)+1, N(3)+1/))
         case ('r')
            els(eID) % r(1:nEqn,0:N(1),0:N(2),0:N(3)) = reshape (vec (firstIdx:lastIdx-1),(/ nEqn, N(1)+1, N(2)+1, N(3)+1/))
         end select
         firstIdx = lastIdx
      end do

   end subroutine MG_Vec2El
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_KSP_Construct(this,DimPrb,KrylovSpace)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), intent(inout) :: this
      integer                 , intent(in)    :: DimPrb
      integer                 , intent(in)    :: KrylovSpace
      integer                    :: KrylovSpace2=20
!  -----------------------------------------------------------------------

      !allocate(this % KSPForMG_t)

      this % KSP % KrylovSpace = KrylovSpace2
      this % KSP % Preconditioner = S_PRECONDITIONER
      allocate(this % KSP % Z (DimPrb,KrylovSpace2+1))
      allocate(this % KSP % V (DimPrb,KrylovSpace2+1))
      allocate(this % KSP % H (KrylovSpace2+1,KrylovSpace2))
      allocate(this % KSP % W (DimPrb))
      allocate(this % KSP % y (KrylovSpace2))
      allocate(this % KSP % cc(KrylovSpace2+1))
      allocate(this % KSP % ss(KrylovSpace2+1))
      allocate(this % KSP % g (KrylovSpace2+1))
   end subroutine MG_KSP_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_KSP_Destruct(this,DimPrb)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), intent(inout) :: this
      integer                 , intent(in)    :: DimPrb
!  -----------------------------------------------------------------------


      deallocate(this % KSP % Z )
      deallocate(this % KSP % V )
      deallocate(this % KSP % H )
      deallocate(this % KSP % W )
      deallocate(this % KSP % y )
      deallocate(this % KSP % cc)
      deallocate(this % KSP % ss)
      deallocate(this % KSP % g )

      !deallocate(this % KSP)
   end subroutine MG_KSP_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_El2Vec(els,vec,N,vecsize,nEqn,eldof,var2interp)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      integer, dimension(3), intent(in)                                    :: N
      integer, intent(in)                                                  :: vecsize
      integer, intent(in)                                                  :: nEqn
      integer, dimension(nelem), intent(in)                                :: eldof
      real(kind=RP), dimension(vecsize), intent(inout)                     :: vec
      ! real(kind=RP), dimension(nEqn,0:N(1), 0:N(2), 0:N(3)), intent(in)     :: els
      type(TemporaryElementStorage_t) , dimension(nelem)   , intent(inout) :: els
      character(len=*)                                     , intent(in)    :: var2interp
!-----Local-Variables-----------------------------------------------------
      integer           :: eID, firstIdx, lastIdx ! counters
!  -----------------------------------------------------------------------
      
      firstIdx = 1
      do eID = 1, nelem
         lastIdx = firstIdx + eldof(eID) * nEqn
         select case (var2interp)
         case ('x')
            vec (firstIdx : lastIdx - 1) = reshape ( els(eID) % x(1:nEqn,0:N(1),0:N(2),0:N(3)) , (/ eldof(eID) * nEqn /) )
         case ('r')
            vec (firstIdx : lastIdx - 1) = reshape ( els(eID) % r(1:nEqn,0:N(1),0:N(2),0:N(3)) , (/ eldof(eID) * nEqn /) )
         end select
         firstIdx = lastIdx
      end do

   end subroutine MG_El2Vec
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_SetOperatorDt(this, dt)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t),     intent(inout)     :: this
      real(kind=RP),                     intent(in)        :: dt
!-----Local-Variables-----------------------------------------------------
      real(kind=RP)                                        :: shift
      real(kind=RP)                                        :: eps = 1e-10
!  -----------------------------------------------------------------------
      
      this % dt = dt
      shift = MatrixShift(dt) !
      if (ABS(shift) .GT. eps) then                  
         call this % A % shift(shift) ! A = A + shift * I
         this % Ashift = shift
      end if
   end subroutine MG_SetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_ReSetOperatorDt(this, dt)
!  ---------------------------------------------------------
!  Removes previous shift in order to insert new one 
!              (important when Jacobian is reused)
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t),     intent(inout)     :: this
      real(kind=RP),                     intent(in)        :: dt
!-----Local-Variables-----------------------------------------------------
      real(kind=RP)                                        :: shift
      real(kind=RP)                                        :: eps = 1e-10
!  -----------------------------------------------------------------------
         
      this % dt = dt
      shift = MatrixShift(dt) !
      if (ABS(shift) .GT. eps) then
         call this % A % Reshift (shift) ! A = A + shift * I
         this % Ashift = shift
      end if
   end subroutine MG_ReSetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_SetRHS(this, RHS)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), intent(inout) :: this
      real(kind=RP)           , intent(in)    :: RHS(this % DimPrb)
!  -----------------------------------------------------------------------
       
      this % b = RHS
       
   end subroutine MG_SetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function MG_GetX(this) result(x)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), intent(inout) :: this
      real(kind=RP)                           :: x(this % DimPrb)
!  -----------------------------------------------------------------------

      x = this % x

   end function MG_GetX
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function MG_Getxnorm(this,TypeOfNorm) result(xnorm)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), intent(inout)  :: this
      character(len=*)                         :: TypeOfNorm
      real(kind=RP)                            :: xnorm
!  -----------------------------------------------------------------------
      
      select case (TypeOfNorm)
         case ('infinity')
            xnorm = MAXVAL(ABS(this % x))
         case ('l2')
            xnorm = NORM2(this % x)
         case default
            stop 'LinearMultigridSolverClass ERROR: Norm not implemented yet'
      end select

   end function MG_Getxnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function MG_Getrnorm(this) result(rnorm)
!  ---------------------------------------------------------
!  Infinity norm. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t), intent(inout) :: this
      real(kind=RP)                            :: rnorm
      !-----------------------------------------------------------
      real(kind=RP)                            :: residual(this % DimPrb)
!  -----------------------------------------------------------------------
      
      rnorm = this % rnorm
      
      
   end function MG_Getrnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_SetRHSValue(this, irow, value)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t) , intent(inout) :: this
      integer                  , intent(in)    :: irow
      real(kind=RP)            , intent(in)    :: value
!  -----------------------------------------------------------------------
      
      this % b (irow) = value
      
   end subroutine MG_SetRHSValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 !
   subroutine MG_SetRHSValues(this, nvalues, irow, values)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(LinearMultigridSolver_t)    , intent(inout)     :: this
      integer                     , intent(in)        :: nvalues
      integer      , dimension(1:), intent(in)        :: irow
      real(kind=RP), dimension(1:), intent(in)        :: values
      integer                                         :: i
!  -----------------------------------------------------------------------

      do i=1, nvalues
         if (irow(i)<0) cycle
         this % b(irow(i)) = values(i)
      end do
      
   end subroutine MG_SetRHSValues
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function getnosmooths( inputline, arrdim ) result(n)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      integer               :: arrdim, cstart, cend, i, j
      integer, dimension(arrdim) :: n
      character ( len = * ) :: inputline
!-----Local-Variables-----------------------------------------------------
!  -----------------------------------------------------------------------
!
      ! cstart = index(inputline,'[')
      ! cend   = index(inputline, ']', .true. )
      do i = 1, arrdim
         j = 1 + (i-1)*3
         read( inputline( j:j+2 ), * ) n(i)
      end do
!
   end function getnosmooths
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------
!  Convert the DBD matrix to a CSR matrix
!  --------------------------------------
   subroutine getCSRfromDBD(this,Acsr)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      ! class(Matrix_t), intent(in)    :: this          !<  This matrix
      class(DenseBlockDiagMatrix_t), intent(in)    :: this          !<  This matrix
      class(csrMat_t)              , intent(inout) :: Acsr      !<  Facorized matrix
!-----Local-Variables-----------------------------------------------------
      integer :: ii, jj
      integer :: bID
!  -----------------------------------------------------------------------
      
      if (this % num_of_Rows /= Acsr % num_of_Rows) then
         print*, 'getCSRfromDBD :: ERROR: Matrix dimensions mismatch:', this % num_of_Rows, Acsr % num_of_Rows
         stop
      end if
      
      call Acsr % PreAllocate()
      call Acsr % Reset
      
      call Acsr % SpecifyBlockInfo(this % BlockIdx, this % BlockSizes)
      
!     Fill the Matrix
!     ---------------
      do bID=1, this % num_of_Blocks

         do jj=1, this % BlockSizes(bID)
            do ii=1, this % BlockSizes(bID)
                  call Acsr % SetBlockEntry(bID,bID,ii,jj, this % Blocks (bID) % Matrix(ii,jj))
            end do
         end do

      end do
      
      call Acsr % assembly()
      
   end subroutine getCSRfromDBD
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------
!  Convert the DBD matrix to a CSR matrix
!  Author: Wojciech Laskowski (wj.laskowski@upm.es) 
!  --------------------------------------
   subroutine getDBDfromCSR(this,Adbd)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments-----------------------------------------------------------
      class(csrMat_t),               intent(in)    :: this          
      class(DenseBlockDiagMatrix_t), intent(inout) :: Adbd    
!-----Local-Variables-----------------------------------------------------
      integer                     :: i,j,k ! counters
      real(kind=rp), allocatable  :: Mat(:,:) ! local dense matrix (one block/element)
!  -----------------------------------------------------------------------
      
!     Dimension check
!     ---------------
      if (this % num_of_Rows /= Adbd % num_of_Rows) then
         print *, 'getDBDfromCSR :: ERROR: Matrix dimensions mismatch:', this % num_of_Rows, Adbd % num_of_Rows
         stop
      end if

!     Fill the DBD Matrix
!     ---------------
      if (size(Adbd%BlockSizes,1) .lt. 1) then
         ERROR Stop "BDB Matrix not allocated correctly."
      end if 
      allocate( Mat(Adbd%BlockSizes(1),Adbd%BlockSizes(1))) ! allocate local matrix, this one if for const. pol. order for all elements
      do i = 1, Adbd % num_of_blocks
         Mat = 0.d0
         do j = 1, Adbd % BlockSizes(i) 
            do k = this % Rows(j + (i-1) * Adbd % BlockSizes(i) ), this % Rows(j + (i-1) * Adbd % BlockSizes(i) + 1)-1
               if ( this%Cols(k) .gt. (i) * Adbd % BlockSizes(i)) exit ! check if we exceed the diagonal block
               if ( this%Cols(k) .gt. (i-1) * Adbd % BlockSizes(i) ) then ! check if we are in the right diagonal block
                  Mat( j , this % Cols(k) - (i-1) * Adbd % BlockSizes(i) ) = this % Values (k)
               end if
            end do
         end do
         Adbd % Blocks(i) % Matrix = Mat
      end do
      deallocate(Mat) ! deallocate 
      
   end subroutine getDBDfromCSR
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ComputeBlockPrec(this)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      ! use DenseMatUtilities
      implicit none
!-----Arguments-----------------------------------------------------------
      class(BJSmooth_t), target, intent(inout) :: this            !  Iterative solver class
!-----Local-Variables-----------------------------------------------------
      integer :: k      ! Counter
!  -----------------------------------------------------------------------

      do k=1, this % A_p % num_of_Blocks
         ! this % A_p % Blocks(k) % Matrix = inverse( this % A_p % Blocks(k) % Matrix ) ! direct factorziation and save at the same place
         call ComputeLU (A        = this % A_p % Blocks(k) % Matrix, &
                         ALU      = this % BlockPrec(k) % PLU, &
                         LUpivots = this % BlockPrec(k) % LUpivots)
      end do
      
   end subroutine ComputeBlockPrec
!   
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MGS_ConstructBlockJacobi(this,p_sem,Nx,Ny,Nz,nEqn)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(BJSmooth_t)            , intent(inout)      :: this                     ! Matrix to solve
      type(DGSem)                  , intent(in)         :: p_sem
      integer                      , intent(in)         :: Nx              
      integer                      , intent(in)         :: Ny              
      integer                      , intent(in)         :: Nz              
      integer                      , intent(in)         :: nEqn              
!-----Local-Variables-----------------------------------------------------
      integer :: k,j            ! Counters
      integer :: ndofelm        ! Dummies
      integer :: nnzs(nelem)
!  -----------------------------------------------------------------------

      allocate (this % A_p)
      call this % A_p % construct (num_of_Blocks = p_sem % mesh % no_of_elements)
      do j=1,p_sem % mesh % no_of_elements
          nnzs(j) = nEqn*(Nx+1)*(Ny+1)*(Nz+1) 
      end do
      call this % A_p % PreAllocate (nnzs=nnzs)

      allocate (this % BlockPrec(this % A_p % num_of_blocks))
      do k = 1, this % A_p % num_of_blocks
         ndofelm = this % A_p % BlockSizes(k)
         allocate (this % BlockPrec(k) % PLU(ndofelm,ndofelm) )
         allocate (this % BlockPrec(k) % LUpivots   (ndofelm) )
      end do
   end subroutine MGS_ConstructBlockJacobi
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MGS_DestructBlockJacobi(this)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(BJSmooth_t), intent(inout)             :: this                     ! Matrix to solve
!-----Local-Variables---------------------------------------------
      integer                                 :: i,j              ! Counters
!  -----------------------------------------------------------------------
   end subroutine MGS_DestructBlockJacobi
!   
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MGS_ConstructILU(this,A)
!  ---------------------------------------------------------
!  Constructor. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(ILUSmooth_t)            , intent(inout)   :: this ! Matrix to solve
      type(csrMat_t)                , intent(in)      :: A
!-----Local-Variables---------------------------------------------
!  -----------------------------------------------------------------------
#ifdef HAS_MKL
      call this % A % constructWithCSRArrays  (A % Rows, A % Cols, A % Values)
#else
      error stop "ILU smoother needs MKL."
#endif
   end subroutine MGS_ConstructILU
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MGS_DestructILU(this)
!  ---------------------------------------------------------
!  FIXME TBD. 
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(ILUSmooth_t), intent(inout)             :: this                     ! Matrix to solve
!-----Local-Variables---------------------------------------------
      integer                                 :: i,j              ! Counters
!  -----------------------------------------------------------------------
   end subroutine MGS_DestructILU
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ComputeILU(this)
!  ---------------------------------------------------------
!  Incomplete factorisation using only MKL routines. The method and the specific paremteres for the
!  input (ipar, dpar) can be found here: 
!  https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/sparse-solver-routines/preconditioners-based-on-incomplete-lu-factorization-technique/dcsrilu0.html 
!  ---------------------------------------------------------
      implicit none
!-----Arguments---------------------------------------------------
      class(ILUSmooth_t)            , intent(inout)   :: this                     ! Matrix to solve
      ! type(csrMat_t)                , intent(in)      :: A
!-----Local-Variables---------------------------------------------
      real(kind=rp), allocatable  :: bilu0(:) ! tmp array containing values of factorised matrix
      real(kind=rp)               :: dpar(128)
      integer                     :: ipar(128) 
      integer                     :: ierr 
      !--------------------------------------------
      integer, dimension(:), allocatable  :: rows,cols
      real(kind=RP), dimension(:), allocatable  :: vals,bilu_test

      type(csrMat_t) :: Atmp
!  -----------------------------------------------------------------------

#ifdef HAS_MKL
      ! initialisation
      ipar(2)  = 6
      ipar(6)  = 1
      ipar(31) = 0

      dpar(31) = 1.e-16
      dpar(32) = 1.e-10

      allocate(bilu0(size(this % A % Values,1)))
      call dcsrilu0  ( this % A % num_of_Rows, this % A % Values, this % A % Rows, this % A % Cols, bilu0 , ipar , dpar , ierr )
      if (ierr .ne. 0) then
         print *, "Error in dscrilu0, ierr: ", ierr
      endif

      ! call this % A % Visualize('Ajac_prefac.dat') ! visualisation before factorisation
      this % A % Values = bilu0
      ! call this % A % Visualize('Ajac_posfac.dat') ! visualisation after factorisation
      deallocate(bilu0)
#else
      error stop "ILU smoother needs MKL."
#endif
   end subroutine ComputeILU
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module LinearMultigridSolverClass
