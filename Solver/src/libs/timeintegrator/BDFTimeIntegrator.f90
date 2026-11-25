!
!   Module for integrating in time using the Backward Differentiation Formulas (BDF)
!
!////////////////////////////////////////////////////////////////////////
MODULE BDFTimeIntegrator
   use SMConstants
   USE PhysicsStorage
   use HexMeshClass
   USE LinearSolverClass
   USE CSRMatrixClass
   USE FTValueDictionaryClass
   use TimeIntegratorDefinitions
   use MatrixClass
   use DGSEMClass
   use StorageClass              , only: SolutionStorage_t
   use StopwatchClass            , only: Stopwatch
   use Utilities                 , only: AlmostEqual, toLower
   use MPI_Process_Info          , only: MPI_Process
   implicit none
   
   PRIVATE                          
   PUBLIC BDFIntegrator_t, TakeBDFStep, ComputeRHS, UpdateNewtonSol, BDF_SetPreviousSolution, bdf_order, BDF_MatrixShift, BDF_SetOrder
   public BDFCoeff, BDFInitialiseQ
   
!
!  ********************
!  BDF integrator class
!  ********************
   type BDFIntegrator_t
      
      class(GenericLinSolver_t), allocatable :: linsolver     !  Linear solver
      integer                                :: maxLinSolIter
      real(kind=RP)                          :: firstNorm
      real(kind=RP)                          :: NewtonFactor
      real(kind=RP)                          :: LinSolTolFactor
      integer                                :: StepsForJac   !· Maximum number of steps that should be taken for computing a new Jacobian matrix
      integer                                :: StepsSinceJac !  
      integer                                :: MaxNewtonIter
      logical                                :: JacByConv     !· .TRUE. if the Jacobian must be computed only when the convergence is bad
      logical                                :: TimeAccurate  !· .TRUE. if this is a time-accurate simulation
      logical                                :: UserNewtonTol !· .TRUE. if the newton tolerance is specified by the user
      real(kind=RP)                          :: NewtonTol     !  Specified Newton tolerance
      real(kind=RP)                          :: inner_dt   
      
      contains
         procedure :: construct => ConstructBDFIntegrator
         procedure :: destruct  => DestructBDFIntegrator
         procedure :: TakeStep  => TakeBDFStep
   end type BDFIntegrator_t
   
!
!  Module variables
!  ----------------
   
   logical       :: computeA = .TRUE.  ! Compute Jacobian? (only valid if it is meant to be computed according to the convergence)
   logical       :: Adaptive_dt = .TRUE.
   ! integer       :: bdf_order       ! BDF order specified by user
   integer       :: order           ! BDF order to be used
   integer       :: StepsTaken = 0

!
!  BDF coefficients for constant time-step
!  ---------------------------------------
   real(kind=RP), parameter :: BDFCoeff(6,5) = &
!                    a_1             a_2     a_3           a_4             a_5          a_6
         reshape( (/ 1.0_RP        , -1._RP, 0._RP       , 0._RP         , 0._RP      , 0._RP        ,  &   ! BDF1
                     1.5_RP        , -2._RP, 0.5_RP      , 0._RP         , 0._RP      , 0._RP        ,  &   ! BDF2
                     11._RP/6_RP   , -3._RP, 3._RP/2._RP , -1._RP/3._RP  , 0._RP      , 0._RP        ,  &   ! BDF3
                     25._RP/12_RP  , -4._RP, 3._RP       , -4._RP/3._RP  , 1._RP/4._RP, 0._RP        ,  &   ! BDF4
                     137._RP/60_RP , -5._RP, 5._RP       , -10._RP/3._RP , 5._RP/4._RP, -1._RP/5._RP /) &   ! BDF5
                                                                                                      , (/6,5/) )
   integer, parameter :: MAX_ORDER_CONS_DT = 5
   
!
!  Default parameters for Newton iterative procedure
!  -------------------------------------------------
   real(kind=RP), parameter   :: NEWTON_MIN_CONVRATE = 0.1_RP     ! Minimum convergence rate for Newton method... If newton loop convergence rate passes this value, inner_dt is decreased
   real(kind=RP), parameter   :: NEWTON_MAX_CONVRATE = 1.7_RP     ! Maximum convergence rate for Newton method... If newton loop convergence rate passes this value, inner_dt is increased
   real(kind=RP), parameter   :: NEWTON_TOL_DEFAULT  = 1.e-6_RP   ! Default convergence tolerance
   integer      , parameter   :: MAX_NEWTON_ITER = 30          ! If newton iter reaches this limit, this iteration is marked as  not converged 
   integer      , parameter   :: LIM_NEWTON_ITER = 12          ! If Newton converges but this limit is reached, jacobian matrix will be recomputed
   logical                    :: PRINT_NEWTON_INFO
   
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine ConstructBDFIntegrator(this,controlVariables,sem)
      implicit none
      !--------------------------------------------------------
      class(BDFIntegrator_t) , intent(inout) :: this
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(DGSem)            , intent(in)    :: sem
      !--------------------------------------------------------
      integer :: DimPrb, globalDimPrb
      character(len=LINE_LENGTH) :: simulationType
      !--------------------------------------------------------
      
!
!     Get general definitions
!     -----------------------

      Adaptive_dt = controlVariables % logicalValueForKey("implicit adaptive dt")
      this % JacByConv = controlVariables % LogicalValueForKey("jacobian by convergence")
      simulationType = controlVariables % StringValueForKey("simulation type",LINE_LENGTH)
      call toLower(simulationType)
      if (simulationType == 'time-accurate') then
         this % TimeAccurate = .TRUE.
      else
         this % TimeAccurate = .FALSE.
      end if
      
      if (controlVariables % containsKey("newton max iter")) then
         this % MaxNewtonIter = controlVariables % integerValueForKey("newton max iter")
      else
         this % MaxNewtonIter = MAX_NEWTON_ITER
      end if
      
      if (controlVariables % containsKey("linsolver max iter")) then
         this % maxLinSolIter = controlVariables % integerValueForKey("linsolver max iter")
      else
         this % maxLinSolIter = 500
      end if

      if (controlVariables % containsKey("newton first norm")) then
         this % firstNorm = controlVariables % doublePrecisionValueForKey("newton first norm")
      else
         this % firstNorm = 2.d-1
      end if

      if (controlVariables % containsKey("newton factor")) then
         this % NewtonFactor = controlVariables % doublePrecisionValueForKey("newton factor")
      else
         this % NewtonFactor = 1e-3_RP
      end if

      if (controlVariables % containsKey("linsolver tol factor")) then
         this % LinSolTolFactor = controlVariables % doublePrecisionValueForKey("linsolver tol factor")
      else
         this % LinSolTolFactor = 0.5_RP
      end if

      PRINT_NEWTON_INFO = controlVariables % logicalValueForKey("print newton info") .and. MPI_Process % isRoot

      if (controlVariables % containsKey("newton tolerance")) then
         this % UserNewtonTol = .TRUE.
         this % NewtonTol = controlVariables % doublePrecisionValueForKey("newton tolerance")
      else
         this % UserNewtonTol = .FALSE.
         this % NewtonTol = NEWTON_TOL_DEFAULT
      end if
      
      this % StepsForJac = controlVariables % integerValueForKey("compute jacobian every")
!
!     Setup linear solver
!     -------------------
      DimPrb = sem % NDOF * NCONS
      globalDimPrb = sem % totalNDOF * NCONS
      
      select case ( trim(controlVariables % StringValueForKey("linear solver",LINE_LENGTH)) )
         case('petsc')
            allocate (PetscKspLinearSolver_t :: this % linsolver)
         case('pardiso')
            allocate (MKLPardisoSolver_t     :: this % linsolver)
         case('matrix-free smooth')
            allocate (MatFreeSmooth_t        :: this % linsolver)
         case('matrix-free gmres')
            allocate (MatFreeGMRES_t         :: this % linsolver)
         case('smooth')
            allocate (IterativeSolver_t      :: this % linsolver)
         case('multigrid')
            allocate (LinearMultigridSolver_t      :: this % linsolver)
         case('static-condensation')
            allocate (StaticCondSolver_t     :: this % linsolver)
         case default
            write(STD_OUT,*) "Keyword 'linear solver' missing... Using PETSc as default"
            allocate (PetscKspLinearSolver_t :: this % linsolver)
      end select
      
      call this % linsolver % construct (DimPrb,globalDimPrb,NCONS,controlVariables,sem,BDF_MatrixShift)             !Constructs linear solver 
      
!
!     Setup BDF methods
!     -----------------
      call BDF_SetOrder( controlVariables % integerValueForKey("bdf order") )
      
      ! Check that the BDF order is consistent
      if (bdf_order > 1) then
         if ( (.not. controlVariables % containsKey("dt") ) .or. Adaptive_dt) then
            error stop ':: "bdf order">1 is only valid with fixed time-step sizes'
         end if
      end if
      
!     Create Stopwatch for solving
!     ----------------------------
      call Stopwatch % CreateNewEvent("BDF Newton-Solve")
      
   end subroutine ConstructBDFIntegrator
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
   subroutine DestructBDFIntegrator(this)
      implicit none
      !--------------------------------------------------------
      class(BDFIntegrator_t), intent(inout) :: this
      !--------------------------------------------------------
      
      call this % linsolver % destroy
      deallocate (this % linsolver)
      
   end subroutine DestructBDFIntegrator
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
   subroutine TakeBDFStep (this, sem, t , dt, ComputeTimeDerivative)
      implicit none
      !----------------------------------------------------------------------
      class(BDFIntegrator_t), intent(inout) :: this
      type(DGSem),                  INTENT(inout)           :: sem                  !<>DGSem class with solution storage 
      real(kind=RP),                INTENT(in)              :: t                    !< Time at the beginning of time step
      real(kind=RP),                INTENT(in)              :: dt                   !< Initial (outer) time step (the subroutine can use a smaller one depending on convergence)
      procedure(ComputeTimeDerivative_f)                            :: ComputeTimeDerivative
      !----------------------------------------------------------------------
      
      real(kind=RP) :: time               ! Time at the beginning of each inner(!) time step
      integer                                               :: k, newtonit
      
      real(kind=RP)                                         :: ConvRate
      logical                                               :: CONVERGED
      !----------------------------------------------------------------------
      
      if ((.not. this % TimeAccurate) .and. (.not. this % UserNewtonTol)) then
         this % NewtonTol = sem % MaxResidual* this % NewtonFactor
      end if
      
      !**************************
      ! If the Jacobian must only be computed sometimes
      if (this % JacByConv) then
         if ( (.not. computeA) .and. (.not. AlmostEqual (this % inner_dt,dt) ) ) then
            call this % linsolver % ReSetOperatorDt(this % inner_dt)
         end if
      end if
      ! 
      !**************************
      
      this % inner_dt = dt            ! first inner_dt is the outer step dt 
      time = t
      
      call sem % mesh % storage % local2GlobalQ(sem % NDOF)
!
!     ********************
!     Sub-time-step solver
!     ********************
      do
         
!        Set previous solution for inner time-step
!        -----------------------------------------
         
         call BDF_SetPreviousSolution(sem % mesh % storage)
         
!        Perform Newton iterative procedure
!        ----------------------------------
         
         if (computeA) then
            this % StepsSinceJac = 0
         else
            this % StepsSinceJac = this % StepsSinceJac + 1
            if (this % StepsSinceJac == this % StepsForJac) then
               computeA = .TRUE.
               this % StepsSinceJac = 0
            end if
         end if
         call NewtonSolve(sem, time+this % inner_dt, this % inner_dt, this % linsolver, this % NewtonTol, this % MaxNewtonIter, this % maxLinSolIter, this % firstNorm, this % LinSolTolFactor, &
                          this % JacByConv,ConvRate, newtonit,CONVERGED, ComputeTimeDerivative)
         
!        Actions if Newton converged
!        ***************************
         if (CONVERGED) then
            time = time + this % inner_dt
            
!           Check convergence to know if the Jacobian must be computed
!           ----------------------------------------------------------
            if (this % JacByConv .and. (newtonit > LIM_NEWTON_ITER) ) then   !Recomputes jacobian Matrix if convergence rate is poor
               if (PRINT_NEWTON_INFO) then
                  write(STD_OUT,*) "Convergence rate is poor, Jacobian matrix will be computed in next iteration..."
               end if
               computeA = .TRUE.                                        
            end if
            
!           Check if the sub time-stepping is done
!           --------------------------------------
            if (ABS((time)-(t+dt)) < 10 * EPSILON(1._RP)) then       ! If outer t+dt is reached, the time integration is done
               exit                                            
            end if
            
!           Increase dt if good convergence in previous step
!           ------------------------------------------------
            if (Adaptive_dt .and. ConvRate > NEWTON_MAX_CONVRATE) then
               this % inner_dt = this % inner_dt * 2.0_RP
               if (this % JacByConv)  call this % linsolver % ReSetOperatorDt(this % inner_dt)    ! Resets the operator with the new dt
               
               if (PRINT_NEWTON_INFO) write(*,*) "Increasing  dt  = ", this % inner_dt
            end if
            
!           Adjust dt to prevent sub time-stepping to be be greater than outer Dt 
!           ---------------------------------------------------------------------
            if ( time+this % inner_dt > t + dt) then  ! Adjusts inner dt to achieve exact outer Dt in the last substep
               this % inner_dt = t + dt - time
               if (this % JacByConv)  call this % linsolver % ReSetOperatorDt(this % inner_dt)    ! Resets the operator with the new dt
               
               if (PRINT_NEWTON_INFO) write(*,*) "Adjusting dt = ", this % inner_dt
            end if
         
!        Actions if Newton did not converge
!        **********************************
         else
            
!           Reduce dt is allowed
!           --------------------
            if (Adaptive_dt) then
               this % inner_dt = this % inner_dt / 2._RP
               if (this % JacByConv)  call this % linsolver % ReSetOperatorDt(this % inner_dt)    ! Resets the operator with the new dt
               
               sem % mesh % storage % Q = sem % mesh % storage % PrevQ(:, sem % mesh % storage % prevSol_index(1))  ! restores Q in sem to begin a new newton iteration
                 
               if (PRINT_NEWTON_INFO) write(*,*) "Newton loop did not converge, trying a smaller dt = ", this % inner_dt
               
!           Warn if dt cannot be changed
!           ----------------------------
            else
               if (this % TimeAccurate) then
                  error stop 'Newton loop did not converge. Consider using a smaller dt or "implicit adaptive dt = .TRUE."'
               else
                  if (MPI_Process % isRoot) write(STD_OUT,*) 'WARNING: Newton loop did not converge.'
                  exit
               end if
            end if
         end if
      
      end do
 
      if (PRINT_NEWTON_INFO) write(*,'(A10,f5.2)') "ConvRate: ", ConvRate
      
      !**************************
      ! for computing sometimes
      if (this % JacByConv .AND. ConvRate <0.65_RP ) then
         computeA = .TRUE.
      end if
      ! for computing sometimes
      !**************************
      
!~       if (MAXVAL(maxResidual) > sem % maxResidual) computeA = .TRUE.
      
      call sem % mesh % storage % global2LocalQ
      call sem % mesh % storage % global2LocalQdot
      
   end subroutine TakeBDFStep
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------
!  Routine for performing a nonlinear Newton iterative procedure
!  -> This can be taken out of the BDFTimeIntegrator if needed
!        (but careful with Adaptive_dt and the Newton vars)
!  -------------------------------------------------------------
   subroutine NewtonSolve(sem, t, dt, linsolver, NEWTON_TOLERANCE, MaxNewtonIter, maxLinSolIter, firstNorm, LinSolTolFactor, JacByConv,ConvRate, niter,CONVERGED, ComputeTimeDerivative)
      implicit none
      !----------------------------------------------------------------------
      type(DGSem),                  intent(inout)           :: sem
      real(kind=RP),                intent(in)              :: t
      real(kind=RP),                intent(in)              :: dt              !< Inner dt
      class(GenericLinSolver_t),    intent(inout)           :: linsolver       !Linear operator is calculate outside this subroutine
      real(kind=RP),                intent(in)              :: NEWTON_TOLERANCE
      integer,                      intent(in)              :: MaxNewtonIter
      integer,                      intent(in)              :: maxLinSolIter
      real(kind=RP),                intent(in)              :: firstNorm
      real(kind=RP),                intent(in)              :: LinSolTolFactor
      logical,                      intent(in)              :: JacByConv         !< Must the Jacobian be computed for bad convergence? if .false., the Jacobian is computed at the beginning of every newton it
      real(kind=RP),                intent(out)             :: ConvRate
      integer,                      intent(out)             :: niter
      logical,                      intent(out)             :: CONVERGED   
      procedure(ComputeTimeDerivative_f)                    :: ComputeTimeDerivative
      !----------------------------------------------------------------------           
      integer              :: newtonit
      real(kind=RP)        :: norm, norm_old, rel_tol
      real(kind=RP)        :: linsolver_tol
      logical      , save  :: isfirst = .TRUE.
      real(kind=RP), save  :: norm1
      !----------------------------------------------------------------------
      
!
!     Initializations
!     ---------------
      
      if (isfirst) then
         norm = firstNorm  ! A value to define the initial linsolver_tol
         isfirst = .FALSE.
      else
         norm = norm1
      end if

      norm_old = -1.0_RP  !Must be initialized to -1 to avoid bad things in the first newton iter
      ConvRate = 1.0_RP
   
      if (PRINT_NEWTON_INFO) then
         write(*, "(A9,1X,A18,1X,A18,1X,A15,1X,A12,1X,A18)") "Newton it", "Newton abs_err", "Newton rel_err", "LinSolverErr", "# ksp iter", "Iter wall time (s)"
      end if
!
!     Newton loop
!     -----------
      DO newtonit = 1, MaxNewtonIter
         
         linsolver_tol = norm * ( LinSolTolFactor**(newtonit) )   ! Use another expression? 0.25?                   ! Nastase approach ("High-order discontinuous Galerkin methods using an hp-multigrid approach")
         ! linsolver_tol = 1e-12
         
         call ComputeRHS(sem, t, dt, linsolver, ComputeTimeDerivative )               ! Computes b (RHS) and stores it into linsolver
         
         call Stopwatch % Start("BDF Newton-Solve")
         call linsolver%solve( nEqn=NCONS, nGradEqn=NGRAD, tol = linsolver_tol, maxiter=maxLinSolIter, time= t, dt=dt, &
                              ComputeTimeDerivative = ComputeTimeDerivative, computeA = computeA)        ! Solve (J-I/dt)·x = (Q_r- U_n)/dt - Qdot_r
         call Stopwatch % Pause("BDF Newton-Solve")

         ! print *, " Solver time: ", Stopwatch % lastelapsedtime("BDF Newton-Solve")
         ! error stop "DONE"
         
         if (.NOT. linsolver%converged .and. Adaptive_dt) then                           ! If linsolver did not converge, return converged=false
            converged = .FALSE.
            return
         end if
         call UpdateNewtonSol(sem, linsolver)                    ! Q_r+1 = Q_r + x
         
         norm = linsolver%Getxnorm('infinity')
         
!        Sometimes, iterative methods take 0 iterations because of a high initial linsolver_tol
!        -> In such cases, here norm = 0 even though that's not true. As a workaround we take the residual norm
!        ------------------------------------------------------------------------------------------------------
         if ( AlmostEqual(norm,0._RP) ) then
            norm = linsolver%Getrnorm()
         end if
         
         if (norm_old > 0._RP) then
            ConvRate = ConvRate + (LOG10(norm_old/norm)-ConvRate)/newtonit 
         end if
         norm_old = norm
         niter = newtonit
         if (newtonit == 1) then
            norm1 = norm
            rel_tol = norm1 * NEWTON_TOLERANCE
         end if
         if (PRINT_NEWTON_INFO) then
            write(*, "(I9,1X,ES18.3,1X,ES18.3,1X,ES15.3,1X,I12,1X,F18.5)")newtonit, norm, norm/norm1, linsolver%Getrnorm(),&
                                                      linsolver%niter, Stopwatch % lastElapsedTime("BDF Newton-Solve")
         end if
         
         if (ConvRate < NEWTON_MIN_CONVRATE .OR. newtonit == MaxNewtonIter .OR. ISNAN(norm)) then
            if (PRINT_NEWTON_INFO) write(STD_OUT,*) 'ConvRate: ', ConvRate
            converged = .FALSE.
            return
         end if
        
         if (norm < max(rel_tol,NEWTON_TOLERANCE)) then ! Careful: this may not be appropriate for unsteady simulations
            converged = .TRUE. 
            return
         end if
         
      end do
   
   end subroutine NewtonSolve
!  
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ComputeRHS(sem, t, dt, linsolver, ComputeTimeDerivative )
      implicit none
      !----------------------------------------------------------------
      type(DGSem),                intent(inout)    :: sem
      real(kind=RP),              intent(in)       :: t
      real(kind=RP),              intent(in)       :: dt
      class(GenericLinSolver_t),  intent (inout)   :: linsolver
      procedure(ComputeTimeDerivative_f)                   :: ComputeTimeDerivative
      !----------------------------------------------------------------
      real(kind=RP)                                :: value
      real(kind=RP)  :: RHS(NCONS*sem % NDOF)
      !----------------------------------------------------------------
      
      call sem % mesh % storage % global2LocalQ
      call ComputeTimeDerivative( sem % mesh, sem % particles, t, CTD_IGNORE_MODE)
      call sem % mesh % storage % local2GlobalQdot(sem % NDOF)
      
      RHS = BDF_GetRHS(sem % mesh % storage, dt)
      
      call linsolver % SetRHS(RHS)
      
      call linsolver % AssemblyRHS     ! b must be assembled before using
   end subroutine ComputeRHS
!  
!/////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine UpdateNewtonSol(sem, linsolver)

      type(DGSem),                     intent(inout)    :: sem
      class(GenericLinSolver_t),       intent(inout)    :: linsolver
      
      sem % mesh % storage % Q = sem % mesh % storage % Q  + linsolver % GetX()
      
   end subroutine UpdateNewtonSol
!
!////////////////////////////////////////////////////////////////////////////////////////////
!  TODO: Move from here....
   subroutine WriteEigenFiles(Mat,sem,FileName)
      IMPLICIT NONE
!
!     -----------------------------------------------------------
!     Writes files for performing eigenvalue analysis using TAUev
!        This only works for isotropic order meshes.........................TODO: Change that
!     -----------------------------------------------------------
!
      type(csrMat_t)    :: Mat      !< Jacobian matrix
      type(DGSem)       :: sem      !< DGSem class containing mesh
      character(len=*)  :: FileName !< ...
!     -----------------------------------------------------------
      integer           :: fd
!     -----------------------------------------------------------
      
      ! .frm file
      OPEN(newunit=fd, file=TRIM(FileName)//'.frm', action='write')
         write(fd,*)
         write(fd,*) SIZE(Mat % Values), SIZE(Mat % Rows)-1, 1, NCONS, 1
         write(fd,*) sem % mesh % elements(1) % Nxyz(1), SIZE(sem % mesh % elements)
      CLOSE (fd)
      
      ! .amg file
      call Mat % Visualize(TRIM(FileName)//'.amg',FirstRow=.FALSE.)
      
      ! .coo file
      call sem % mesh % WriteCoordFile(NCONS, TRIM(FileName)//'.coo')
      
      
   end subroutine WriteEigenFiles
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine BDF_SetOrder(order)
      implicit none
      !------------------------------------------------------
      integer, intent(in) :: order
      !------------------------------------------------------
      
      if (order > MAX_ORDER_CONS_DT) then
         write(STD_OUT,*) 'WARNING :: Maximum BDF order for constant time-step is 5. Using 1 by default.'
         bdf_order = 1
      else
         bdf_order = order
      end if
      
   end subroutine BDF_SetOrder
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine BDF_SetPreviousSolution(storage,NotANewStep)
      implicit none
      !------------------------------------------------------
      type(SolutionStorage_t), intent(inout) :: storage
      logical, optional :: NotANewStep
      !------------------------------------------------------
      integer :: i      ! Counter
      !------------------------------------------------------
      
      if (present(NotANewStep)) then
         if (.not. NotANewStep) StepsTaken = StepsTaken + 1
      else
         StepsTaken = StepsTaken + 1
      end if
      
      order = min(StepsTaken, bdf_order)
      
      call storage % SetGlobalPrevQ(storage % Q)
   end subroutine BDF_SetPreviousSolution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function BDF_MatrixShift(dt) result(Ashift)
      implicit none
      !------------------------------------------------------
      real(kind=RP), intent(in) :: dt
      real(kind=RP)             :: Ashift
      !------------------------------------------------------
      
      Ashift = -BDFCoeff(1,order)/dt
      
   end function BDF_MatrixShift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function BDF_GetRHS(storage,dt) result(RHS)
      implicit none
      !------------------------------------------------------
      type(SolutionStorage_t), intent(in) :: storage
      real(kind=RP)          , intent(in) :: dt
      real(kind=RP)                       :: RHS(size(storage % Q))
      !------------------------------------------------------
      integer :: k
      real(kind=RP) :: invdt
      !------------------------------------------------------
      
      invdt = 1._RP/dt
      
      RHS = storage % Q * BDFCoeff(1,order)*invdt - storage % Qdot
      
      do k=1, order
         RHS = RHS + BDFCoeff(k+1,order) * storage % PrevQ(:,storage % prevSol_index(k)) * invdt
      end do
      
   end function BDF_GetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine BDFInitialiseQ(mesh)
      implicit none
!-----Arguments-----------------------------------------------------------
      type(HexMesh)                                        :: mesh
!-----Local-Variables-----------------------------------------------------
      integer :: i, id
!-------------------------------------------------------------------------
      
      do i = 1, bdf_order
!$omp parallel do schedule(runtime)
         do id = 1, SIZE( mesh % elements )
            mesh % elements(id) % storage % prevQ(i) % Q = mesh % elements(id) % storage % Q
         end do ! id
!$omp end parallel do
      end do
      
   end subroutine BDFInitialiseQ
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE BDFTimeIntegrator