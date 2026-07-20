!//////////////////////////////////////////////////////
!
!  This module provides the routines for computing the analytical Jacobian matrix
!  -> Only for p-conforming representations (TODO: make general)
!  -> The Jacobian of the BCs is temporarily computed numerically
!
!//////////////////////////////////////////////////////
#include "Includes.h"
module AnalyticalJacobian
   use SMConstants
   use ElementClass
   use HexMeshClass
   use NodalStorageClass
   use PhysicsStorage
   use Physics
   use JacobianDefinitions             , only: JACEPS
   use JacobianComputerClass           , only: JacobianComputer_t
   use MatrixClass
   use DGSEMClass                      , only: DGSem, ComputeTimeDerivative_f
   use StopWatchClass
   use MeshTypes
#ifndef ACOUSTIC
   use EllipticDiscretizations
#endif
   use MPI_Process_Info                , only: MPI_Process
   use ElementConnectivityDefinitions  , only: axisMap, normalAxis
   use BoundaryConditions              , only: BCs
   use FaceClass                       , only: Face
   use Utilities                       , only: dot_product
   use ConnectivityClass               , only: Connectivity
   use FTValueDictionaryClass
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
   use RiemannSolvers_NS               , only: RiemannSolver_dFdQ, RiemannSolver
   use HyperbolicDiscretizations       , only: HyperbolicDiscretization
   use VariableConversion              , only: NSGradientVariables_STATE
   use FluidData_NS, only: dimensionless
#elif defined(SPALARTALMARAS)
   use RiemannSolvers_NSSA               , only: RiemannSolver_dFdQ, RiemannSolver
   use HyperbolicDiscretizations       , only: HyperbolicDiscretization
   use VariableConversion              , only: NSGradientVariables_STATE
   use FluidData_NSSA, only: dimensionless
#endif
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none
   
   private
   public AnJacobian_t
   
   integer, parameter :: other(2) = [2, 1]
   
!
!  **************************************************
!  Main type for the analytical Jacobian computations
!  **************************************************
   type, extends(JacobianComputer_t) :: AnJacobian_t
      
      contains
         procedure :: Construct => AnJacobian_Construct
         procedure :: Compute   => AnJacobian_Compute
   end type AnJacobian_t
   
   
#if defined(NAVIERSTOKES)
   real(kind=RP) :: Identity(NCONS,NCONS) ! identity matrix. TODO: Define only once in the constructor (When this is a class!!)
#endif
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------
!  AnJacobian_Construct:
!  Main class constructor
!  -------------------------------------------------------
   subroutine AnJacobian_Construct(this, mesh, nEqn, controlVariables)
      implicit none
      !-arguments-----------------------------------------
      class(AnJacobian_t)  , intent(inout) :: this
      type(HexMesh)        , intent(inout) :: mesh
      integer              , intent(in)    :: nEqn
      type(FTValueDictionary)  , intent(in)    :: controlVariables
      !-local-variables-----------------------------------
      integer :: fID, eID
      !---------------------------------------------------
      
!
!     Construct parent
!     ----------------
      call this % JacobianComputer_t % construct (mesh, nEqn, controlVariables)
!
!     Create stopwatch event
!     ----------------------
      call Stopwatch % CreateNewEvent("Analytical Jacobian construction")
!
!     Construct specific storage
!     --------------------------
      ! Global storage
      mesh % storage % anJacobian = .TRUE.
      ! Elements:
!$omp parallel do schedule(runtime)
      do eID=1, mesh % no_of_elements
         call mesh % storage % elements(eID) % ConstructAnJac ()
      end do
!$omp end parallel do
      
      ! Faces:
!$omp parallel do schedule(runtime)
      do fID=1, size(mesh % faces)
         call mesh % faces(fID) % storage % ConstructAnJac (NDIM)
      end do
!$omp end parallel do
      
      !TODO: Add conformity check
      
   end subroutine AnJacobian_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------
!  Subroutine for computing the analytical Jacobian matrix
!  -------------------------------------------------------
   subroutine AnJacobian_Compute(this, sem, nEqn, time, Matrix, TimeDerivative, TimeDerivativeIsolated, eps_in, BlockDiagonalized, mode)
      implicit none
      !-arguments----------------------------------
      class(AnJacobian_t)      , intent(inout)     :: this
      type(DGSem)              , intent(inout)     :: sem
      integer,                   intent(in)        :: nEqn
      real(kind=RP)            , intent(in)        :: time
      class(Matrix_t)          , intent(inout)     :: Matrix
      procedure(ComputeTimeDerivative_f), optional :: TimeDerivative    ! Not needed here...
      procedure(ComputeTimeDerivative_f), optional :: TimeDerivativeIsolated
      real(kind=RP)  , optional, intent(in)        :: eps_in            ! Not needed here...
      logical        , optional, intent(in)        :: BlockDiagonalized !<? Construct only the block diagonal?
      integer        , optional, intent(in)        :: mode
#if defined(NAVIERSTOKES) 
      !--------------------------------------------
      integer :: nnz
      integer :: nelem, ierr
      logical :: BlockDiagonal
      integer :: eID, i
      !--------------------------------------------
      
      
      call Stopwatch % Start("Analytical Jacobian construction")
!
!     Initializations
!     ---------------
      if ( present(BlockDiagonalized) ) then
         BlockDiagonal = BlockDiagonalized
      else
         BlockDiagonal = .FALSE.
      end if
      
      nelem = size(sem % mesh % elements)

      !TODO: Check if there was p-adaptation and reconstruct if necessary
      
!
!     *************************************************
!     If Jacobian matrix was not preallocated, allocate
!     *************************************************
!
      if (.not. this % preAllocate) then
         select type(Matrix_p => Matrix)
!
!           If block-diagonal matrix, construct with size of blocks
!           -------------------------------------------------------
            type is(DenseBlockDiagMatrix_t)
               call Matrix_p % Preallocate(nnzs=this % ndofelm_l)
            type is(SparseBlockDiagMatrix_t)
               call Matrix_p % Preallocate(nnzs=this % ndofelm_l)
               
!
!           If matrix is CSR, standard preallocate with LinkedListMatrix
!           ------------------------------------------------------------
            type is(csrMat_t)
               call Matrix_p % Preallocate()
!
!           Otherwise, construct with nonzeros in each row
!           ----------------------------------------------
            class default 
               if (BlockDiagonal) then
                  nnz = maxval(this % ndofelm_l)
               else
                  nnz = 7*maxval(this % ndofelm_l) ! 20180201: hard-coded to 7, since the analytical Jacobian can only compute off-diagonal blocks for compact schemes (neighbors' effect)
               end if
               call Matrix_p % Preallocate(nnz)
         end select
         
         call Matrix % SpecifyBlockInfo(this % firstIdx,this % ndofelm)
      end if
      
      
      call Matrix % Reset (ForceDiagonal = .TRUE.)
      
!$omp parallel
!     ------------------
!     Prolong Q to faces
!     ------------------
!
      call sem % mesh % ProlongSolutionToFaces(NCONS)
!
!     ----------------
!     Update MPI Faces
!     ----------------
!
#ifdef _HAS_MPI_
!$omp single
      call sem % mesh % UpdateMPIFacesSolution(NCONS)
!$omp end single
#endif
!
!     ******************************************************************
!     If the physics has an elliptic component, compute the DG gradients
!     -> This routine also projects the appropriate grads to the faces
!     ******************************************************************
!
      if (flowIsNavierStokes) then
         call ViscousDiscretization % ComputeGradient (nEqn, nEqn, sem % mesh, time, NSGradientVariables_STATE)
#ifdef _HAS_MPI_
!$omp single
         call sem % mesh % UpdateMPIFacesGradients(NGRAD)
!$omp end single
#endif
      end if
!
!     ************************************************************************
!     Compute the Jacobian of the Numerical Fluxes \hat{f}^a and \hat{f}^{\nu}
!     ************************************************************************
!
      call ComputeNumericalFluxJacobian(sem % mesh,nEqn,time)
!
!     ***************
!     Diagonal blocks
!     ***************
!
      call AnalyticalJacobian_DiagonalBlocks(sem % mesh, nEqn,time, Matrix)
!
!     *******************
!     Off-Diagonal blocks
!     *******************
!
      if (.not. BlockDiagonal) call AnalyticalJacobian_OffDiagonalBlocks(sem % mesh,Matrix)
!$omp end parallel
!
!     ************************
!     Finish assembling matrix
!     ************************
!
!     Add an MPI Barrier
!     ------------------
#ifdef _HAS_MPI_
      if (MPI_Process % doMPIAction) then
         call mpi_barrier(MPI_COMM_WORLD, ierr)
      end if
#endif
      call Matrix % Assembly()
!
!     *********
!     Finish up
!     *********
!
      call Stopwatch % Pause("Analytical Jacobian construction")
      
      if ( this % verbose ) then
         write(STD_OUT,'(A,ES10.3,A)') "Analytical Jacobian construction: ", Stopwatch % lastElapsedtime("Analytical Jacobian construction"), ' seconds'
      end if
#else
      error stop ':: Analytical Jacobian only for NS'
#endif

   ! call Matrix % Visualize('Jacobian.txt')
   end subroutine AnJacobian_Compute
#if defined(NAVIERSTOKES) 
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------------------------
!  Subroutine for adding the faces' contribution to the diagonal blocks of the Jacobian matrix
!  -------------------------------------------------------------------------------------------
   subroutine AnalyticalJacobian_DiagonalBlocks(mesh,nEqn,time,Matrix)
      implicit none
      !--------------------------------------------
      type(HexMesh), target    , intent(inout) :: mesh
      integer,                   intent(in)    :: nEqn
      real(kind=RP)            , intent(in)    :: time
      class(Matrix_t)          , intent(inout) :: Matrix
      !--------------------------------------------
      integer :: eID, fID
      type(Element), pointer :: e   
      !--------------------------------------------
!
!     Project flux Jacobian to corresponding element
!     ----------------------------------------------
!$omp do schedule(runtime)
      do fID = 1, size(mesh % faces)
         if (mesh % faces(fID)% MortarType == MORTAR_NONE .OR. mesh % faces(fID)% MortarType == MORTAR_SMALL4) then 

            if (mesh % faces(fID)% MortarType == MORTAR_NONE) call mesh % faces(fID) % ProjectFluxJacobianToElements(nEqn,LEFT ,LEFT )   ! dF/dQL to the left element 
            if (.not. (mesh % faces(fID) % faceType == HMESH_BOUNDARY )) call mesh % faces(fID) % ProjectFluxJacobianToElements(nEqn,RIGHT,RIGHT)   ! dF/dQR to the right element

         elseif(mesh % faces(fID)% MortarType == MORTAR_BIG) then 

            call mesh % faces(fID) % ProjectMortarFluxJacobianToElements(nEqn,LEFT ,LEFT,  mesh % faces(fID+1),mesh % faces(fID+2), &
            mesh % faces(fID+3), mesh % faces(fID+4))   ! dF/dQL to the left element 

         end if 
      end do
!$omp end do
!
!     Project flux Jacobian with respect to gradients to corresponding elements
!     -------------------------------------------------------------------------
      if (flowIsNavierStokes) then
!$omp do schedule(runtime)
         do fID = 1, size(mesh % faces)
            if (mesh % faces(fID)% MortarType == MORTAR_NONE .OR. mesh % faces(fID)% MortarType == MORTAR_SMALL4) then 

               if (mesh % faces(fID)% MortarType == MORTAR_NONE) call mesh % faces(fID) % ProjectGradJacobianToElements(LEFT, LEFT)   ! dF/dQL to the left element 
               if (.not. (mesh % faces(fID) % faceType == HMESH_BOUNDARY)) call mesh % faces(fID) % ProjectGradJacobianToElements(RIGHT,RIGHT)   ! dF/dQR to the right element

            elseif(mesh % faces(fID)% MortarType == MORTAR_BIG) then 

               call mesh % faces(fID) % ProjectMortarGradJacobianToElements(LEFT ,LEFT,  mesh % faces(fID+1),mesh % faces(fID+2), &
               mesh % faces(fID+3), mesh % faces(fID+4))   ! dF/dQL to the left element 

            end if 
         end do
!$omp end do
      end if
!
!     Compute each element's diagonal block
!     -------------------------------------
!$omp do schedule(runtime) private(e)
      do eID = 1, size(mesh % elements)
         e => mesh % elements(eID)
         call Local_SetDiagonalBlock( e, &
                                      mesh % faces( e % faceIDs(EFRONT ) ), &
                                      mesh % faces( e % faceIDs(EBACK  ) ), &
                                      mesh % faces( e % faceIDs(EBOTTOM) ), &
                                      mesh % faces( e % faceIDs(ERIGHT ) ), &
                                      mesh % faces( e % faceIDs(ETOP   ) ), &
                                      mesh % faces( e % faceIDs(ELEFT  ) ), &
                                      Matrix,                               &
                                      mesh% IBM                             )
      end do
!$omp end do
      nullify (e)
   end subroutine AnalyticalJacobian_DiagonalBlocks   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  Subroutine for adding the faces' contribution to the off-diagonal blocks of the Jacobian matrix
!  -> Note: Only the interior faces contribute to the off-diagonal blocks
!  -> only for p-conforming meshes
!  -----------------------------------------------------------------------------------------------
   subroutine AnalyticalJacobian_OffDiagonalBlocks(mesh,Matrix)
      implicit none
      !--------------------------------------------
      type(HexMesh), target    , intent(inout) :: mesh
      class(Matrix_t)          , intent(inout) :: Matrix
      !--------------------------------------------
      integer :: eID, fID, elSide, side
      type(Element), pointer :: e_plus
      type(Face)   , pointer :: f
      !--------------------------------------------
      
!
!     Project flux Jacobian to opposed elements (RIGHT to LEFT and LEFT to RIGHT)
!     ---------------------------------------------------------------------------
!$omp do schedule(runtime)
      do fID = 1, size(mesh % faces)
         if (mesh % faces(fID) % faceType /= HMESH_BOUNDARY) then
            if (mesh % faces(fID) % MortarType == MORTAR_NONE) then 

               call mesh % faces(fID) % ProjectFluxJacobianToElements(NCONS, LEFT ,RIGHT)   ! dF/dQR to the left element
               call mesh % faces(fID) % ProjectFluxJacobianToElements(NCONS, RIGHT,LEFT )   ! dF/dQL to the right element 

            elseif (mesh % faces(fID) % MortarType == MORTAR_SMALL4) then 

               call mesh % faces(fID) % ProjectFluxJacobianToElements(NCONS, RIGHT,LEFT )   ! dF/dQL to the right element 

            elseif (mesh % faces(fID) % MortarType == MORTAR_BIG) then 

               call mesh % faces(fID) % ProjectMortarFluxJacobianToElements(NCONS, LEFT ,RIGHT, mesh % faces(fID+1), mesh % faces(fID+2), &
               mesh % faces(fID+3), mesh % faces(fID+4))   ! dF/dQR to the left element

            end if 
         end if
      end do
!$omp end do

!
!     Project flux Jacobian with respect to gradients to opposed elements (RIGHT to LEFT and LEFT to RIGHT)
!     -----------------------------------------------------------------------------------------------------
      if (flowIsNavierStokes) then
!$omp do schedule(runtime)
         do fID = 1, size(mesh % faces)
            if (mesh % faces(fID) % faceType /= HMESH_BOUNDARY) then
               if (mesh % faces(fID) % MortarType == MORTAR_NONE) then 

                  call mesh % faces(fID) % ProjectGradJacobianToElements(LEFT ,RIGHT)   ! dF/dGradQR to the left element
                  call mesh % faces(fID) % ProjectGradJacobianToElements(RIGHT,LEFT )   ! dF/dGradQL to the right element 

               elseif (mesh % faces(fID) % MortarType == MORTAR_SMALL4) then 

                  call mesh % faces(fID) % ProjectGradJacobianToElements(RIGHT,LEFT ) 

               elseif (mesh % faces(fID) % MortarType == MORTAR_BIG) then 

                  call mesh % faces(fID) % ProjectMortarGradJacobianToElements(LEFT ,RIGHT, mesh % faces(fID+1), mesh % faces(fID+2), &
                  mesh % faces(fID+3), mesh % faces(fID+4))
                  
               end if 
            end if
         end do
!$omp end do
      end if
      
!
!     Compute the off-diagonal blocks for each element's equations
!     ------------------------------------------------------------
!$omp do schedule(runtime) private(e_plus,elSide,fID,side,f)
      do eID = 1, mesh % no_of_elements
         e_plus => mesh % elements(eID)
!
!        One block for every neighbor element
!        ------------------------------------
         do elSide = 1, 6
            if (e_plus % NumberOfConnections(elSide) == 0) cycle
            
            fID  = e_plus % faceIDs(elSide)
            side = e_plus % faceSide(elSide)
            
            f => mesh % faces(fID)
            
            call Local_GetOffDiagonalBlock(f,e_plus,e_plus % Connection(elSide),side,Matrix)
         end do
      end do
!$omp end do
      nullify (f, e_plus)
   end subroutine AnalyticalJacobian_OffDiagonalBlocks
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  
!  -----------------------------------------------------------------------------------------------
   subroutine ComputeNumericalFluxJacobian(mesh,nEqn,time)
      implicit none
      !--------------------------------------------
      type(HexMesh), intent(inout)    :: mesh
      integer,       intent(in)       :: nEqn
      real(kind=RP), intent(in)       :: time
      !--------------------------------------------
      integer :: fID, ierr
      !--------------------------------------------
!
!     Compute the non-shared faces
!     ----------------------------
!
!$omp do schedule(runtime)
      do fID = 1, size(mesh % faces)
         select case (mesh % faces(fID) % faceType)
            case (HMESH_INTERIOR)
               if (mesh % faces(fID)% MortarType == MORTAR_NONE) then 
                  call ComputeInterfaceFluxJacobian(mesh % faces(fID))
               elseif (mesh % faces(fID)%MortarType == MORTAR_BIG) then 
                  call  ComputeInterfaceMortarFluxJacobian(mesh % faces(fID), mesh % faces(fID+1), mesh % faces(fID+2), &
                  mesh % faces(fID+3), mesh % faces(fID+4))
               end if 
            case (HMESH_BOUNDARY)
               call ComputeBoundaryFluxJacobian(mesh % faces(fID),time)
         end select
      end do
!$omp end do
!
!     Compute the MPI faces
!     ---------------------
!
#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
!
!        Wait until the MPI messages have been received
!        ----------------------------------------------
!$omp single
         if ( flowIsNavierStokes ) then 
            call mesh % GatherMPIFacesGradients(NGRAD)
         else  
            call mesh % GatherMPIFacesSolution(NCONS)
         end if
!$omp end single
!$omp do schedule(runtime)
         do fID = 1, size(mesh % faces)
            select case (mesh % faces(fID) % faceType)
               case (HMESH_MPI)
                  call ComputeInterfaceFluxJacobian(mesh % faces(fID))
            end select
         end do
!$omp end do
!
!        Add an MPI Barrier
!        ------------------
!$omp single
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
      end if
#endif
   end subroutine ComputeNumericalFluxJacobian

   subroutine ComputeInterfaceMortarFluxJacobian(f, fma, fmb, fmc, fmd ) 
      implicit none
      !--------------------------------------------
      type(Face), intent(inout) :: f
      type(Face), intent(inout) :: fma
      type(Face), intent(inout) :: fmb 
      type(Face), intent(inout) :: fmc 
      type(Face), intent(inout) :: fmd 
      !--------------------------------------------
      integer :: i,j, lm 
      integer :: Nfm(4,2) 

      Nfm(1,:)=fma % Nf
      Nfm(2,:)=fmb % Nf
      Nfm(3,:)=fmc % Nf
      Nfm(4,:)=fmd % Nf
      !--------------------------------------------
      
     do lm=1,4 
      do j = 0, Nfm(lm,2) ; do i = 0, Nfm(lm,1)
!
!        Get numerical flux jacobian on the face point (i,j)
!        ---------------------------------------------------
        select case(lm)
        case(1)
         
         call RiemannSolver_dFdQ(ql   = fma % storage(LEFT)  % Q(:,i,j), &
                                 qr   = fma % storage(RIGHT) % Q(:,i,j), &
                                 nHat = fma % geom % normal (:,i,j)    , &
                                 dfdq_num = fma % storage(LEFT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqL
                                 side = LEFT)
         call RiemannSolver_dFdQ(ql   = fma % storage(LEFT)  % Q(:,i,j), &
                                 qr   = fma % storage(RIGHT) % Q(:,i,j), &
                                 nHat = fma % geom % normal (:,i,j)    , &
                                 dfdq_num = fma % storage(RIGHT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqR
                                 side = RIGHT)
!
!        Scale with the mapping Jacobian
!        -------------------------------
         fma % storage(LEFT ) % dFStar_dqF (:,:,i,j) = fma % storage(LEFT  ) % dFStar_dqF (:,:,i,j) * fma % geom % jacobian(i,j)
         fma % storage(RIGHT) % dFStar_dqF (:,:,i,j) = fma % storage(RIGHT ) % dFStar_dqF (:,:,i,j) * fma % geom % jacobian(i,j)

        case(2) 

        call RiemannSolver_dFdQ(ql   = fmb % storage(LEFT)  % Q(:,i,j), &
                                 qr   = fmb % storage(RIGHT) % Q(:,i,j), &
                                 nHat = fmb % geom % normal (:,i,j)    , &
                                 dfdq_num = fmb % storage(LEFT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqL
                                 side = LEFT)
        call RiemannSolver_dFdQ(ql   = fmb % storage(LEFT)  % Q(:,i,j), &
                                 qr   = fmb % storage(RIGHT) % Q(:,i,j), &
                                 nHat = fmb % geom % normal (:,i,j)    , &
                                 dfdq_num = fmb % storage(RIGHT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqR
                                 side = RIGHT)
!
!        Scale with the mapping Jacobian
!        -------------------------------
         fmb % storage(LEFT ) % dFStar_dqF (:,:,i,j) = fmb % storage(LEFT  ) % dFStar_dqF (:,:,i,j) * fmb % geom % jacobian(i,j)
         fmb % storage(RIGHT) % dFStar_dqF (:,:,i,j) = fmb % storage(RIGHT ) % dFStar_dqF (:,:,i,j) * fmb % geom % jacobian(i,j)

        case(3)
        
        call RiemannSolver_dFdQ(ql   = fmc % storage(LEFT)  % Q(:,i,j), &
                                 qr   = fmc % storage(RIGHT) % Q(:,i,j), &
                                 nHat = fmc % geom % normal (:,i,j)    , &
                                 dfdq_num = fmc % storage(LEFT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqL
                                 side = LEFT)
        call RiemannSolver_dFdQ(ql   = fmc % storage(LEFT)  % Q(:,i,j), &
                                 qr   = fmc % storage(RIGHT) % Q(:,i,j), &
                                 nHat = fmc % geom % normal (:,i,j)    , &
                                 dfdq_num = fmc % storage(RIGHT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqR
                                 side = RIGHT)
!
!        Scale with the mapping Jacobian
!        -------------------------------
         fmc % storage(LEFT ) % dFStar_dqF (:,:,i,j) = fmc % storage(LEFT  ) % dFStar_dqF (:,:,i,j) * fmc % geom % jacobian(i,j)
         fmc % storage(RIGHT) % dFStar_dqF (:,:,i,j) = fmc % storage(RIGHT ) % dFStar_dqF (:,:,i,j) * fmc % geom % jacobian(i,j)

        case(4)

        call RiemannSolver_dFdQ(ql   = fmd % storage(LEFT)  % Q(:,i,j), &
                                 qr   = fmd % storage(RIGHT) % Q(:,i,j), &
                                 nHat = fmd % geom % normal (:,i,j)    , &
                                 dfdq_num = fmd % storage(LEFT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqL
                                 side = LEFT)
        call RiemannSolver_dFdQ(ql   = fmd % storage(LEFT)  % Q(:,i,j), &
                                 qr   = fmd % storage(RIGHT) % Q(:,i,j), &
                                 nHat = fmd % geom % normal (:,i,j)    , &
                                 dfdq_num = fmd % storage(RIGHT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqR
                                 side = RIGHT)
!
!        Scale with the mapping Jacobian
!        -------------------------------
         fmd % storage(LEFT ) % dFStar_dqF (:,:,i,j) = fmd % storage(LEFT  ) % dFStar_dqF (:,:,i,j) * fmd % geom % jacobian(i,j)
         fmd % storage(RIGHT) % dFStar_dqF (:,:,i,j) = fmd % storage(RIGHT ) % dFStar_dqF (:,:,i,j) * fmd % geom % jacobian(i,j)
        end select 
      end do             ; end do
#ifndef SPALARTALMARAS
      if (flowIsNavierStokes)  then 
      if (lm==1) call ViscousDiscretization % RiemannSolver_Jacobians(fma)
      if (lm==2) call ViscousDiscretization % RiemannSolver_Jacobians(fmb)
      if (lm==3) call ViscousDiscretization % RiemannSolver_Jacobians(fmc)
      if (lm==4) call ViscousDiscretization % RiemannSolver_Jacobians(fmd)
      end if 
#endif

    end do 

   end subroutine ComputeInterfaceMortarFluxJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  Computes the Jacobian of the numerical flux on a boundary as
!     dF*    dF*(Q⁺,Q⁻,n^)     dF*(Q⁺,Q⁻,n^)     dQ⁻
!    ---- = --------------- + --------------- * ----,
!     dQ⁺      dQ⁺               dQ⁻             dQ⁺
!  where the last term is the Jacobian of the boundary condition.... +: internal state
!                                                                    -: external state
!  -----------------------------------------------------------------------------------------------
   subroutine ComputeBoundaryFluxJacobian(f,time)
      implicit none
      !--------------------------------------------
      type(Face), intent(inout) :: f
      real(kind=RP), intent(in) :: time
      !--------------------------------------------
      integer :: i,j, n, m
      real(kind=RP) :: BCjac(NCONS,NCONS,0:f % Nf(1),0:f % Nf(2))
      real(kind=RP) :: dF_dQ     (NCONS,NCONS)
      real(kind=RP) :: dF_dgradQ (NCONS,NCONS,NDIM)
      !--------------------------------------------
      
!
!     *********************
!     Inviscid contribution
!     *********************
!
      do j = 0, f % Nf(2) ; do i = 0, f % Nf(1) 
!
!        Get external state
!        ------------------
         
         f % storage(RIGHT) % Q(:,i,j) = f % storage(LEFT) % Q(:,i,j)
         CALL BCs(f % zone) % bc % StateForEqn( NCONS, &
                                      f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(RIGHT) % Q(:,i,j))
!
!        Get numerical flux jacobian on the face point (i,j)
!        ---------------------------------------------------

         call RiemannSolver_dFdQ(ql   = f % storage(LEFT)  % Q(:,i,j), &
                                 qr   = f % storage(RIGHT) % Q(:,i,j), &
                                 nHat = f % geom % normal (:,i,j)    , &
                                 dfdq_num = f % storage(LEFT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqL
                                 side = LEFT)
         
         call RiemannSolver_dFdQ(ql   = f % storage(LEFT)  % Q(:,i,j), &
                                 qr   = f % storage(RIGHT) % Q(:,i,j), &
                                 nHat = f % geom % normal (:,i,j)    , &
                                 dfdq_num = f % storage(RIGHT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqR
                                 side = RIGHT)
!
!        Scale with the mapping Jacobian
!        -------------------------------
         f % storage(LEFT ) % dFStar_dqF (:,:,i,j) = f % storage(LEFT ) % dFStar_dqF (:,:,i,j) * f % geom % jacobian(i,j)
         f % storage(RIGHT) % dFStar_dqF (:,:,i,j) = f % storage(RIGHT) % dFStar_dqF (:,:,i,j) * f % geom % jacobian(i,j)
         
!        Correct dFstar/dQL with the Jacobian of the boundary condition
!        --------------------------------------------------------------
         call ExternalStateJacobian( f % geom % x(:,i,j), &
                                     time, &
                                     f % geom % normal(:,i,j), &
                                     f % storage(1) % Q(:,i,j),&
                                     f % storage(2) % Q(:,i,j),&
                                     f % zone, &
                                     BCjac(:,:,i,j) )
         
         f % storage(LEFT ) % dFStar_dqF (:,:,i,j) = f % storage(LEFT ) % dFStar_dqF (:,:,i,j) &
                                            + matmul(f % storage(RIGHT) % dFStar_dqF (:,:,i,j),BCjac(:,:,i,j))
      end do             ; end do
      
!
!     ********************
!     Viscous contribution - Temporarily done fully numerically. TODO: This can be improved
!     ********************
!
      if (flowIsNavierStokes) then
      
         call f % ProjectBCJacobianToElements(NCONS,BCjac)
      
         do j = 0, f % Nf(2) ; do i = 0, f % Nf(1) 
            
            call ViscousExternalStateJacobian(f, &
                                              f % geom % x(:,i,j), &
                                              time, &
                                              f % geom % normal(:,i,j), &
                                              f % geom % dWall(i,j),&
                                              f % geom % t1(:,i,j), &
                                              f % geom % t2(:,i,j), &
                                              f % storage(1) % Q(:,i,j),&
                                              f % storage(1) % U_x(:,i,j),&
                                              f % storage(1) % U_y(:,i,j),&
                                              f % storage(1) % U_z(:,i,j),&
                                              f % zone, &
                                              dF_dQ, &
                                              dF_dgradQ )
                                        
            f % storage(LEFT ) % dFStar_dqF (:,:,i,j)     = f % storage(LEFT ) % dFStar_dqF (:,:,i,j)  + dF_dQ     * f % geom % jacobian(i,j)
            
            associate ( dF_dGradQ_out => f % storage(LEFT ) % dFv_dGradQF(:,:,:,i,j), &
                        nHat => f % geom % normal(:,i,j))
         
            dF_dGradQ_out = 0._RP
            do n = 1, NDIM
               dF_dGradQ_out(:,:,1) = dF_dGradQ_out(:,:,1) + dF_dgradQ(:,:,n) * f % geom % GradXi  (n,i,j)
               dF_dGradQ_out(:,:,2) = dF_dGradQ_out(:,:,2) + dF_dgradQ(:,:,n) * f % geom % GradEta (n,i,j)
               dF_dGradQ_out(:,:,3) = dF_dGradQ_out(:,:,3) + dF_dgradQ(:,:,n) * f % geom % GradZeta(n,i,j)
            end do
            
            dF_dGradQ_out = dF_dGradQ_out * f % geom % jacobian(i,j)
            
            end associate
            
         end do             ; end do
      end if
      
   end subroutine ComputeBoundaryFluxJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  Computes the Jacobians of the numerical flux for an inner face as
!                         dF*(Q⁺,Q⁻,n^)                        dF*(Q⁺,Q⁻,n^) 
!           dFStar/dqL = ---------------         dFStar/dqR = ---------------
!                          dQ⁺                                  dQ⁻         
!  -----------------------------------------------------------------------------------------------
   subroutine ComputeInterfaceFluxJacobian(f)
      implicit none
      !--------------------------------------------
      type(Face), intent(inout) :: f
      !--------------------------------------------
      integer :: i,j
      !--------------------------------------------
      
      do j = 0, f % Nf(2) ; do i = 0, f % Nf(1) 
!
!        Get numerical flux jacobian on the face point (i,j)
!        ---------------------------------------------------
         
         call RiemannSolver_dFdQ(ql   = f % storage(LEFT)  % Q(:,i,j), &
                                 qr   = f % storage(RIGHT) % Q(:,i,j), &
                                 nHat = f % geom % normal (:,i,j)    , &
                                 dfdq_num = f % storage(LEFT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqL
                                 side = LEFT)
         call RiemannSolver_dFdQ(ql   = f % storage(LEFT)  % Q(:,i,j), &
                                 qr   = f % storage(RIGHT) % Q(:,i,j), &
                                 nHat = f % geom % normal (:,i,j)    , &
                                 dfdq_num = f % storage(RIGHT) % dFStar_dqF (:,:,i,j), & ! this is dFStar/dqR
                                 side = RIGHT)
!
!        Scale with the mapping Jacobian
!        -------------------------------
         f % storage(LEFT ) % dFStar_dqF (:,:,i,j) = f % storage(LEFT  ) % dFStar_dqF (:,:,i,j) * f % geom % jacobian(i,j)
         f % storage(RIGHT) % dFStar_dqF (:,:,i,j) = f % storage(RIGHT ) % dFStar_dqF (:,:,i,j) * f % geom % jacobian(i,j)
         
      end do             ; end do
#ifndef SPALARTALMARAS
      if (flowIsNavierStokes) call ViscousDiscretization % RiemannSolver_Jacobians(f)
#endif

   end subroutine ComputeInterfaceFluxJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  This routine obtains the Jacobian of the boundary condition numerically.
!  This can be optimized introducing the analytical Jacobian of every single implemented BC...
!  -----------------------------------------------------------------------------------------------
   subroutine ExternalStateJacobian(x,time,nHat,Qin,Qex,zone,BCjac)
      implicit none
      !--------------------------------------------
      real(kind=RP), intent(in)       :: x(3)
      real(kind=RP), intent(in)       :: time
      real(kind=RP), intent(in)       :: nHat(3)
      real(kind=RP), intent(in)       :: Qin(NCONS)
      real(kind=RP), intent(in)       :: Qex(NCONS)
      integer,       intent(in)       :: zone
      real(kind=RP), intent(out)      :: BCjac(NCONS,NCONS)
      !--------------------------------------------
      real(kind=RP) :: newQext (NCONS)
      real(kind=RP) :: q(NCONS), buffer
      real(kind=RP),parameter :: eps = 1.e-8_RP
      integer :: i
      !--------------------------------------------
      
      q = Qin
      
      do i = 1, NCONS
         buffer = q(i)
         q(i) = q(i) + eps
         newQext = q
         CALL BCs(zone) % bc % StateForEqn( NCONS, x, time, nHat, newQext)
         
         BCjac(:,i) = (newQext-Qex)/eps
         
         q(i) = buffer
      end do
      
   end subroutine ExternalStateJacobian
!~!
!~!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!~!
!~!  -----------------------------------------------------------------------------------------------
!~!  This routine obtains the Jacobian of the Neumann boundary condition numerically.
!~!  This can be optimized introducing the analytical Jacobian of every single implemented BC...
!~!  -----------------------------------------------------------------------------------------------
!~   subroutine ExternalGradientJacobian(x,time,nHat,Qin,Qex,zone,BCjac)
!~      implicit none
!~      !--------------------------------------------
!~      real(kind=RP), intent(in)       :: x(3)
!~      real(kind=RP), intent(in)       :: time
!~      real(kind=RP), intent(in)       :: nHat(3)
!~      real(kind=RP), intent(in)       :: Qin(NCONS)
!~      real(kind=RP), intent(in)       :: Qex(NCONS)
!~      integer,       intent(in)       :: zone
!~      real(kind=RP), intent(out)      :: BCjac(NCONS,NCONS)
!~      !--------------------------------------------
!~      real(kind=RP) :: newQext (NCONS)
!~      real(kind=RP) :: q(NCONS), buffer
!~      real(kind=RP),parameter :: eps = 1.e-8_RP
!~      integer :: i
!~      !--------------------------------------------
      
!~      q = Qin
!~      !! THIS IS NOT FINISHED YET!
      
!~      do dir=1, NDIM
!~         do i = 1, NCONS
            
!~            gradQR(:,1) = Q_xL
!~            gradQR(:,2) = Q_yL
!~            gradQR(:,3) = Q_zL
            
!~            gradQR(i,dir) = gradQR(i,dir) + eps
            
!~            CALL BCs(f % zone) % bc % FlowNeumann(&
!~                                              x, &
!~                                              time, &
!~                                              nHat, &
!~                                              qr, &
!~                                              gradQR(:,1), &
!~                                              gradQR(:,2), &
!~                                              gradQR(:,3) )
            
            
            
            
!~            call ViscousDiscretization % RiemannSolver ( NCONS, &
!~                                                         NCONS, &
!~                                                         f, &
!~                                                         q , qr , &
!~                                                         Q_xL , Q_yL , Q_zL , &
!~                                                         gradQR(:,1) , gradQR(:,2) , gradQR(:,3) , &
!~                                                         mu, beta, kappa, &
!~                                                         nHat , dWall, newFlux )
            
!~            dF_dgradQ(:,i,dir) = (newFlux-flux)/eps
            
!~         end do
!~      end do
      
!~   end subroutine ExternalGradientJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  This routine obtains the Jacobian of the boundary viscous flux numerically.
!  -----------------------------------------------------------------------------------------------
   subroutine ViscousExternalStateJacobian(f, x,time,nHat, dWall, t1, t2,QL, Q_xL, Q_yL, Q_zL ,zone,dF_dQ,dF_dgradQ)
      implicit none
      !--------------------------------------------
      type(Face)   , intent(in)       :: f
      real(kind=RP), intent(in)       :: x(3)
      real(kind=RP), intent(in)       :: time
      real(kind=RP), intent(in)       :: nHat(3)
      real(kind=RP), intent(in)       :: dWall
      real(kind=RP), intent(in)       :: t1(3)
      real(kind=RP), intent(in)       :: t2(3)
      real(kind=RP), intent(in)       :: QL(NCONS)
      real(kind=RP), intent(in)       :: Q_xL(NCONS)
      real(kind=RP), intent(in)       :: Q_yL(NCONS)
      real(kind=RP), intent(in)       :: Q_zL(NCONS)
      integer,       intent(in)       :: zone
      real(kind=RP), intent(out)      :: dF_dQ(NCONS,NCONS)
      real(kind=RP), intent(out)      :: dF_dgradQ(NCONS,NCONS,NDIM)
      !--------------------------------------------
      real(kind=RP) :: newFlux (NCONS), flux(NCONS), fv_3d(NCONS,NDIM)
      real(kind=RP) :: q(NCONS), buffer, qr(NCONS)
      real(kind=RP) :: mu, beta, kappa
      real(kind=RP) :: gradQL(NCONS,NDIM)
      real(kind=RP),parameter :: eps = 1.e-8_RP
      integer :: i, dir
      !--------------------------------------------
      
      dF_dQ     = 0._RP
      dF_dgradQ = 0._RP
      
      mu    = dimensionless % mu !+ f % storage(1) % mu_art(1,i,j)
      beta  = 0._RP !f % storage(1) % mu_art(2,i,j)
      kappa = dimensionless % kappa !+ f % storage(1) % mu_art(3,i,j)
      
      q = QL
!
!     Get base flux
!     -------------
      gradQL(:,1) = Q_xL
      gradQL(:,2) = Q_yL
      gradQL(:,3) = Q_zL
      
!
!     Compute the viscous flux
!     ------------------------
      call ViscousFlux_STATE(NCONS,NCONS,q,Q_xL,Q_yL,Q_zL,mu,beta,kappa,fv_3d)
      flux = fv_3d(:,IX)*nHat(IX) + fv_3d(:,IY)*nHat(IY) + fv_3d(:,IZ)*nHat(IZ)
      CALL BCs(f % zone) % bc % FlowNeumann( & 
                                    x,       & 
                                    time,    & 
                                    nHat,    & 
                                    q,       & 
                                    Q_xL,    & 
                                    Q_yL,    & 
                                    Q_zL,    & 
                                    flux )
      
!
!     Get contribution to dF_dq
!     -------------------------
      do i = 1, NCONS
         buffer = q(i)
         q(i) = q(i) + eps
         
         call ViscousFlux_STATE(NCONS,NCONS,q,Q_xL,Q_yL,Q_zL,mu,beta,kappa,fv_3d)
         newflux = fv_3d(:,IX)*nHat(IX) + fv_3d(:,IY)*nHat(IY) + fv_3d(:,IZ)*nHat(IZ)
         CALL BCs(f % zone) % bc % FlowNeumann( & 
                                       x,       & 
                                       time,    & 
                                       nHat,    & 
                                       q,       & 
                                       Q_xL,    & 
                                       Q_yL,    & 
                                       Q_zL,    & 
                                       newflux )
         
         dF_dQ(:,i) =  dF_dQ(:,i) - (newFlux-flux)/eps
         
         q(i) = buffer
      end do
         
!
!     Get contribution to dF_dgradQ
!     -----------------------------
      do dir=1, NDIM
         do i = 1, NCONS
            buffer = gradQL(i,dir)
            gradQL(i,dir) = gradQL(i,dir) + eps

            call ViscousFlux_STATE(NCONS,NCONS,q,gradQL(:,1),gradQL(:,2),gradQL(:,3),mu,beta,kappa,fv_3d)
            newflux = fv_3d(:,IX)*nHat(IX) + fv_3d(:,IY)*nHat(IY) + fv_3d(:,IZ)*nHat(IZ)
            CALL BCs(f % zone) % bc % FlowNeumann( & 
                                       x,       & 
                                       time,    & 
                                       nHat,    & 
                                       q,       & 
                                       gradQL(:,1), & 
                                       gradQL(:,2), & 
                                       gradQL(:,3), & 
                                       newflux )

            dF_dgradQ(:,i,dir) = (newFlux-flux)/eps 
            gradQL(i,dir) = buffer
            
         end do
      end do
         
   end subroutine ViscousExternalStateJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  RiemannSolverJacobian:
!     This routine obtains the Jacobian of any Riemann Solver numerically.
!     Currently not used, but can be activated if one wants to use a Riemann Solver that has no
!     analytical Jacobian implemented
!  -----------------------------------------------------------------------------------------------
   subroutine RiemannSolverJacobian(QLeft, QRight, nHat, t1, t2, side, BCjac)
      implicit none
      !--------------------------------------------
      real(kind=RP), intent(in)       :: QLeft (NCONS)
      real(kind=RP), intent(in)       :: QRight(NCONS)
      real(kind=RP), intent(in)       :: nHat(3)
      real(kind=RP), intent(in)       :: t1(3)
      real(kind=RP), intent(in)       :: t2(3)
      integer      , intent(in)       :: side
      real(kind=RP), intent(out)      :: BCjac(NCONS,NCONS)
      !--------------------------------------------
      real(kind=RP) :: flux (NCONS)
      real(kind=RP) :: newQext (NCONS)
      real(kind=RP) :: q(NCONS), buffer
      real(kind=RP),parameter :: eps = 1.e-8_RP
      integer :: i
      !--------------------------------------------
      
      call RiemannSolver (Qleft, QRight, nHat, t1, t2, flux)
      
      select case (side)
      case(LEFT)
         q = QLeft
         do i = 1, NCONS
            buffer = q(i)
            q(i) = q(i) + eps
            
            call RiemannSolver (q, QRight, nHat, t1, t2, newQext)
            
            BCjac(:,i) = (newQext-flux)/eps
            
            q(i) = buffer
         end do
         
      case(RIGHT)
         q = QRight
         do i = 1, NCONS
            buffer = q(i)
            q(i) = q(i) + eps
            
            call RiemannSolver (QLeft, q, nHat, t1, t2, newQext)
            
            BCjac(:,i) = (newQext-flux)/eps
            
            q(i) = buffer
         end do
      end select
   end subroutine RiemannSolverJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------
!  Local_SetDiagonalBlock:
!     Subroutine to set a diagonal block in the Jacobian
!  -----------------------------------------------------
   subroutine Local_SetDiagonalBlock(e, fF, fB, fO, fR, fT, fL, Matrix, IBM)
      use IBMClass
      implicit none
      !-------------------------------------------
      type(Element)     , intent(inout) :: e
      type(Face), target, intent(in)    :: fF, fB, fO, fR, fT, fL !< The six faces of the element
      class(Matrix_t)   , intent(inout) :: Matrix
      type(IBM_type),     intent(inout) :: IBM
      !-------------------------------------------
      real(kind=RP) :: MatEntries(NCONS,NCONS)
      real(kind=RP) :: dFdQ      (NCONS,NCONS,NDIM,0:e%Nxyz(1),0:e%Nxyz(2),0:e%Nxyz(3))
      integer :: i, j             ! Matrix indices
      integer :: r                ! Additional index for counting
      integer :: i1, j1, k1, eq1  ! variable counters for row
      integer :: i2, j2, k2, eq2  ! variable counters for column
      integer :: i12, j12, k12    ! variable counters when (row index) = (column index)
      integer :: baseRow, baseCol ! Position of NCONS by NCONS miniblock of Jacobian
      integer :: nXi, nEta        ! Number of nodes in every direction
      integer :: EtaSpa, ZetaSpa  ! Spacing for these two coordinate directions
      integer :: Deltas           ! A variable to know if enough deltas are zero, in which case this is a zero entry of the matrix..
      integer :: sideF
      integer :: sideB
      integer :: sideO
      integer :: sideR
      integer :: sideT
      integer :: sideL
      real(kind=RP), dimension(:,:,:,:)      , pointer :: dfdq_fr, dfdq_ba, dfdq_bo, dfdq_to, dfdq_le, dfdq_ri
      real(kind=RP), dimension(:,:,:,:,:)    , pointer :: dfdGradQ_fr, dfdGradQ_ba, dfdGradQ_bo, dfdGradQ_to, dfdGradQ_ri, dfdGradQ_le
      real(kind=RP), dimension(:,:,:,:,:,:,:), pointer :: dF_dgradQ
      real(kind=RP) :: nF( NDIM, 0:e % Nxyz(1), 0:e % Nxyz(3) )   ! Normal vecs for FRONT  face
      real(kind=RP) :: nB( NDIM, 0:e % Nxyz(1), 0:e % Nxyz(3) )   ! Normal vecs for BACK   face
      real(kind=RP) :: nO( NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2) )   ! Normal vecs for BOTTOM face
      real(kind=RP) :: nR( NDIM, 0:e % Nxyz(2), 0:e % Nxyz(3) )   ! Normal vecs for RIGHT  face
      real(kind=RP) :: nT( NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2) )   ! Normal vecs for TOP    face
      real(kind=RP) :: nL( NDIM, 0:e % Nxyz(2), 0:e % Nxyz(3) )   ! Normal vecs for LEFT   face
      real(kind=RP) :: xiAux  (NCONS,NCONS,4)
      real(kind=RP) :: etaAux (NCONS,NCONS,4)
      real(kind=RP) :: zetaAux(NCONS,NCONS,4)
      real(kind=RP), pointer :: Gvec_xi(:,:,:), Gvec_eta (:,:,:), Gvec_zeta(:,:,:)
      real(kind=RP) :: dqHat_dqp ! Derivative of solution numerical flux, qHat, with respect to q⁺ (viscous method)
      real(kind=RP) :: dqHat_dqm ! Derivative of solution numerical flux, qHat, with respect to q⁻ (viscous method)
      real(kind=RP) :: temp
      real(kind=RP) :: JacF( NCONS,NCONS, 0:e % Nxyz(1), 0:e % Nxyz(3) )   ! Jacobian for FRONT  face
      real(kind=RP) :: JacB( NCONS,NCONS, 0:e % Nxyz(1), 0:e % Nxyz(3) )   ! Jacobian for BACK   face
      real(kind=RP) :: JacO( NCONS,NCONS, 0:e % Nxyz(1), 0:e % Nxyz(2) )   ! Jacobian for BOTTOM face
      real(kind=RP) :: JacR( NCONS,NCONS, 0:e % Nxyz(2), 0:e % Nxyz(3) )   ! Jacobian for RIGHT  face
      real(kind=RP) :: JacT( NCONS,NCONS, 0:e % Nxyz(1), 0:e % Nxyz(2) )   ! Jacobian for TOP    face
      real(kind=RP) :: JacL( NCONS,NCONS, 0:e % Nxyz(2), 0:e % Nxyz(3) )   ! Jacobian for LEFT   face
      type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta
      
      optional :: IBM      
      
      !-------------------------------------------
      spAxi   => NodalStorage(e % Nxyz(1))
      spAeta  => NodalStorage(e % Nxyz(2))
      spAzeta => NodalStorage(e % Nxyz(3))
!
!     *******************
!     Initial definitions
!     *******************
!
      nXi   = e % Nxyz(1) + 1
      nEta  = e % Nxyz(2) + 1
      EtaSpa  = NCONS*nXi
      ZetaSpa = NCONS*nXi*nEta
      
!     TODO: Fill this only once in class constructor
!     ----------------------------------------------
      Identity = 0._RP
      do i=1, NCONS
         Identity(i,i) = 1._RP
      end do
      
      dF_dgradQ => e % storage % dF_dgradQ
      
      sideF = e % faceSide(EFRONT )
      sideB = e % faceSide(EBACK  )
      sideO = e % faceSide(EBOTTOM)
      sideR = e % faceSide(ERIGHT )
      sideT = e % faceSide(ETOP   )
      sideL = e % faceSide(ELEFT  )
      
!     Get normal vectors in element indexing
!     --------------------------------------      
      call AnJac_GetNormalVec(fF, sideF, nF)
      call AnJac_GetNormalVec(fB, sideB, nB)
      call AnJac_GetNormalVec(fO, sideO, nO)
      call AnJac_GetNormalVec(fR, sideR, nR)
      call AnJac_GetNormalVec(fT, sideT, nT)
      call AnJac_GetNormalVec(fL, sideL, nL)
      
!
!     *********************
!     Inviscid contribution
!     *********************
!
#ifndef SPALARTALMARAS     
      call HyperbolicDiscretization % ComputeInnerFluxJacobian( e, dFdQ) 
      if (flowIsNavierStokes) then
         call ViscousDiscretization % ComputeInnerFluxJacobian( e, dF_dgradQ, dFdQ)
         ! TODO: Read this from Viscous discretization
         dqHat_dqp = 0.5_RP
         dqHat_dqm = 0.5_RP
         
         call AnJac_GetFaceJac(fF, sideF, JacF, dqHat_dqp, dqHat_dqm)
         call AnJac_GetFaceJac(fB, sideB, JacB, dqHat_dqp, dqHat_dqm)
         call AnJac_GetFaceJac(fO, sideO, JacO, dqHat_dqp, dqHat_dqm)
         call AnJac_GetFaceJac(fR, sideR, JacR, dqHat_dqp, dqHat_dqm)
         call AnJac_GetFaceJac(fT, sideT, JacT, dqHat_dqp, dqHat_dqm)
         call AnJac_GetFaceJac(fL, sideL, JacL, dqHat_dqp, dqHat_dqm)
         
      end if
#endif    
!
!     Pointers to the flux Jacobians with respect to q on the faces
!     -------------------------------------------------------------
      dfdq_fr(1:,1:,0:,0:) => fF % storage(sideF) % dFStar_dqEl(1:,1:,0:,0:,sideF)
      dfdq_ba(1:,1:,0:,0:) => fB % storage(sideB) % dFStar_dqEl(1:,1:,0:,0:,sideB)
      dfdq_bo(1:,1:,0:,0:) => fO % storage(sideO) % dFStar_dqEl(1:,1:,0:,0:,sideO)
      dfdq_to(1:,1:,0:,0:) => fT % storage(sideT) % dFStar_dqEl(1:,1:,0:,0:,sideT)
      dfdq_ri(1:,1:,0:,0:) => fR % storage(sideR) % dFStar_dqEl(1:,1:,0:,0:,sideR)
      dfdq_le(1:,1:,0:,0:) => fL % storage(sideL) % dFStar_dqEl(1:,1:,0:,0:,sideL)
      
!     Xi contribution (dj*dk)
!     -----------------------
      do k12 = 0, e % Nxyz(3) ; do j12 = 0, e % Nxyz(2) ; do i2 = 0, e % Nxyz(1) ; do i1 = 0, e % Nxyz(1)
         
         MatEntries = 0._RP
         
         MatEntries = MatEntries + &
                           (       dFdQ(:,:,1,i2,j12,k12) * spAXi   % hatD(i1,i2) &                           ! Volumetric contribution
                            -   dfdq_ri(:,:,j12,k12) * spAXi  % b(i1,RIGHT ) * spAXi  % v(i2,RIGHT )    & ! Face 4 Right
                            -   dfdq_le(:,:,j12,k12) * spAXi  % b(i1,LEFT  ) * spAXi  % v(i2,LEFT  )  )   ! Face 6 Left
                            
         MatEntries = MatEntries * e % geom % invJacobian(i1,j12,k12) ! Scale with Jacobian from mass matrix
            
         baseRow = i1*NCONS + j12*EtaSpa + k12*ZetaSpa
         baseCol = i2*NCONS + j12*EtaSpa + k12*ZetaSpa
         
         do eq2 = 1, NCONS ; do eq1 = 1, NCONS 
            i = eq1 + baseRow  ! row index (1-based)
            j = eq2 + baseCol  ! column index (1-based)
            
            call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))  ! TODO: This can be improved by setting the whole matrix at once
            
         end do            ; end do
      end do                  ; end do                  ; end do                 ; end do
      
!     Eta contribution (di*dk)
!     ------------------------
      do k12 = 0, e % Nxyz(3) ; do j2 = 0, e % Nxyz(2) ; do j1 = 0, e % Nxyz(2) ; do i12 = 0, e % Nxyz(1)
         
         MatEntries = 0._RP
         
         MatEntries = MatEntries + &
                           (       dFdQ(:,:,2,i12,j2,k12) * spAEta  % hatD(j1,j2)                           & ! Volumetric contribution
                            -   dfdq_fr(:,:,i12,k12) * spAeta % b(j1,FRONT ) * spAeta % v(j2,FRONT )    & ! Face 1 Front
                            -   dfdq_ba(:,:,i12,k12) * spAeta % b(j1,BACK  ) * spAeta % v(j2,BACK  )  )   ! Face 2 Back
                            
         MatEntries = MatEntries * e % geom % invJacobian(i12,j1,k12) ! Scale with Jacobian from mass matrix
            
         baseRow = i12*NCONS + j1*EtaSpa + k12*ZetaSpa
         baseCol = i12*NCONS + j2*EtaSpa + k12*ZetaSpa
         
         do eq2 = 1, NCONS ; do eq1 = 1, NCONS 
            i = eq1 + baseRow  ! row index (1-based)
            j = eq2 + baseCol  ! column index (1-based)
            
            call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
            
         end do            ; end do
      end do                  ; end do                  ; end do                 ; end do
      
!     Zeta contribution (di*dj)
!     -------------------------
      do k2 = 0, e % Nxyz(3)  ; do k1 = 0, e % Nxyz(3)   ; do j12 = 0, e % Nxyz(2)  ; do i12 = 0, e % Nxyz(1)
         
         MatEntries = 0._RP
         
         MatEntries = MatEntries + &
                           (       dFdQ(:,:,3,i12,j12,k2) * spAZeta % hatD(k1,k2)                           & ! Volumetric contribution
                            -   dfdq_bo(:,:,i12,j12) * spAZeta% b(k1,BOTTOM) * spAZeta% v(k2,BOTTOM)    & ! Face 3 Bottom
                            -   dfdq_to(:,:,i12,j12) * spAZeta% b(k1,TOP   ) * spAZeta% v(k2,TOP   )  )   ! Face 5 Top
         
         MatEntries = MatEntries * e % geom % invJacobian(i12,j12,k1) ! Scale with Jacobian from mass matrix
            
         baseRow = i12*NCONS + j12*EtaSpa + k1*ZetaSpa
         baseCol = i12*NCONS + j12*EtaSpa + k2*ZetaSpa
         
         do eq2 = 1, NCONS ; do eq1 = 1, NCONS 
            i = eq1 + baseRow  ! row index (1-based)
            j = eq2 + baseCol  ! column index (1-based)
            
            call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
            
         end do            ; end do
      end do                  ; end do                   ; end do                   ; end do
      nullify (dfdq_fr, dfdq_ba, dfdq_bo, dfdq_to, dfdq_le, dfdq_ri)
      
!
!     ********************
!     Viscous contribution
!     ********************
!      
      if (flowIsNavierStokes) then
!
!        Pointers to the flux Jacobians with respect to grad(q) on the faces
!        -------------------------------------------------------------------
         dfdGradQ_fr(1:,1:,1:,0:,0:) => fF % Storage(sideF) % dFv_dGradQEl(1:,1:,1:,0:,0:,sideF)
         dfdGradQ_ba(1:,1:,1:,0:,0:) => fB % Storage(sideB) % dFv_dGradQEl(1:,1:,1:,0:,0:,sideB)
         dfdGradQ_bo(1:,1:,1:,0:,0:) => fO % Storage(sideO) % dFv_dGradQEl(1:,1:,1:,0:,0:,sideO)
         dfdGradQ_to(1:,1:,1:,0:,0:) => fT % Storage(sideT) % dFv_dGradQEl(1:,1:,1:,0:,0:,sideT)
         dfdGradQ_ri(1:,1:,1:,0:,0:) => fR % Storage(sideR) % dFv_dGradQEl(1:,1:,1:,0:,0:,sideR)
         dfdGradQ_le(1:,1:,1:,0:,0:) => fL % Storage(sideL) % dFv_dGradQEl(1:,1:,1:,0:,0:,sideL)

!
!        (dj*dk) contribution
!        ********************
         do k12 = 0, e % Nxyz(3) ; do j12 = 0, e % Nxyz(2) ; do i2 = 0, e % Nxyz(1) ; do i1 = 0, e % Nxyz(1)
            MatEntries = 0._RP
            
!           Get Xi auxiliary terms
!           ---------------------
            xiAux   = 0._RP
            
            do r=0, e % Nxyz(1)
               
               temp = spAxi   % hatD(i1,r) * e % geom % invJacobian(r,j12,k12)
               
               Gvec_xi   => dF_dgradQ(:,:,:,1,r,j12,k12)
               
               xiAux(:,:,4)   = xiAux(:,:,4)   + temp *  ( spAXi   % b(r,LEFT  ) * spAXi   % v(i2,LEFT  ) * matmul( dot_product( Gvec_xi  , nL(:,j12,k12) ) , JacL(:,:,j12,k12) ) &
                                                          +spAXi   % b(r,RIGHT ) * spAXi   % v(i2,RIGHT ) * matmul( dot_product( Gvec_xi  , nR(:,j12,k12) ) , JacR(:,:,j12,k12) ) )
               temp = temp * spAxi   % hatD(r,i2)
               
               xiAux  (:,:,1:3) = xiAux  (:,:,1:3) + temp * Gvec_xi
            end do
            
!
!           Compute contributions to the Xi-component of the flux
!           -----------------------------------------------------
            MatEntries =  MatEntries &
                  +  dot_product(xiAux  (:,:,1:3), e % geom % jGradXi  (:,i2,j12,k12) ) - xiAux(:,:,4)        & ! Volumetric contribution: xi-gradient to xi-component of the flux
                  +   dfdGradQ_ri(:,:,1,j12,k12) * spAXi   % b(i1,RIGHT ) * spAXi   % vd(i2,RIGHT ) & ! Right face contribution (outer surface integral)
                  +   dfdGradQ_le(:,:,1,j12,k12) * spAXi   % b(i1,LEFT  ) * spAXi   % vd(i2,LEFT  )   ! Left face contribution  (outer surface integral)
            
            MatEntries = MatEntries * e % geom % invJacobian(i1,j12,k12) ! Scale with Jacobian from mass matrix
            
            baseRow = i1*NCONS + j12*EtaSpa + k12*ZetaSpa
            baseCol = i2*NCONS + j12*EtaSpa + k12*ZetaSpa
            
            !-------------------------------------
            do eq2 = 1, NCONS ; do eq1 = 1, NCONS 
               i = eq1 + baseRow  ! row index (1-based)
               j = eq2 + baseCol  ! column index (1-based)
               
               call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
               
            end do            ; end do
            !-------------------------------------
         end do                  ; end do                  ; end do                 ; end do
!
!        (di*dk) contribution
!        ********************
         do k12 = 0, e % Nxyz(3) ; do j2 = 0, e % Nxyz(2) ; do j1 = 0, e % Nxyz(2) ; do i12 = 0, e % Nxyz(1)
            MatEntries = 0._RP
            
!           Get Eta auxiliary terms
!           ----------------------
            etaAux   = 0._RP
            
            do r=0, e % Nxyz(2)
               
               temp = spAEta  % hatD(j1,r) * e % geom % invJacobian(i12,r,k12)
               
               Gvec_eta  => dF_dgradQ(:,:,:,2,i12,r,k12)
               
               etaAux(:,:,4)  = etaAux(:,:,4)  + temp * ( spAEta  % b(r,FRONT ) * spAEta  % v(j2,FRONT ) * matmul( dot_product( Gvec_eta , nF(:,i12,k12) ), JacF(:,:,i12,k12)) &
                                                         +spAEta  % b(r,BACK  ) * spAEta  % v(j2,BACK  ) * matmul( dot_product( Gvec_eta , nB(:,i12,k12) ), JacB(:,:,i12,k12)) )
               temp = temp * spAEta  % hatD(r,j2)
               
               etaAux (:,:,1:3) = etaAux (:,:,1:3) + temp * Gvec_eta
            end do
!
!           Compute contributions to the Eta-component of the flux
!           ------------------------------------------------------
            MatEntries =  MatEntries &
                  + dot_product(etaAux (:,:,1:3), e % geom % jGradEta (:,i12,j2,k12) ) - etaAux(:,:,4)        & ! Volumetric contribution: eta-gradient to eta-component of the flux
                  +   dfdGradQ_fr(:,:,2,i12,k12) * spAEta  % b(j1,FRONT ) * spAEta  % vd(j2,FRONT ) & ! Front face contribution (outer surface integral)
                  +   dfdGradQ_ba(:,:,2,i12,k12) * spAEta  % b(j1,BACK  ) * spAEta  % vd(j2,BACK  )   ! Back face contribution (outer surface integral)
            
            MatEntries = MatEntries * e % geom % invJacobian(i12,j1,k12) ! Scale with Jacobian from mass matrix
            
            baseRow = i12*NCONS + j1*EtaSpa + k12*ZetaSpa
            baseCol = i12*NCONS + j2*EtaSpa + k12*ZetaSpa
            
            !-------------------------------------
            do eq2 = 1, NCONS ; do eq1 = 1, NCONS 
               i = eq1 + baseRow  ! row index (1-based)
               j = eq2 + baseCol  ! column index (1-based)
               
               call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
               
            end do            ; end do
            !-------------------------------------
         end do                  ; end do                  ; end do                 ; end do
!
!        (di*dj) contribution
!        ********************
         do k2 = 0, e % Nxyz(3)  ; do k1 = 0, e % Nxyz(3)   ; do j12 = 0, e % Nxyz(2)  ; do i12 = 0, e % Nxyz(1)
            MatEntries = 0._RP
            
!           Get Zeta auxiliary terms
!           -----------------------
            zetaAux   = 0._RP
                  
            do r=0, e % Nxyz(3)
               
               temp = spAZeta % hatD(k1,r) * e % geom % invJacobian(i12,j12,r)
               
               Gvec_zeta => dF_dgradQ(:,:,:,3,i12,j12,r)
               
               zetaAux(:,:,4) = zetaAux(:,:,4) + temp * ( spAZeta % b(r,BOTTOM) * spAZeta % v(k2,BOTTOM) * matmul( dot_product( Gvec_zeta, nO(:,i12,j12) ), JacO(:,:,i12,j12)) &
                                                         +spAZeta % b(r,TOP   ) * spAZeta % v(k2,TOP   ) * matmul( dot_product( Gvec_zeta, nT(:,i12,j12) ), JacT(:,:,i12,j12)) )
               temp = temp * spAZeta % hatD(r,k2)
               
               zetaAux(:,:,1:3) = zetaAux(:,:,1:3) + temp * Gvec_zeta
            end do
!
!           Compute contributions to the Zeta-component of the flux
!           -------------------------------------------------------
            MatEntries = MatEntries &
                  + dot_product(zetaAux(:,:,1:3), e % geom % jGradZeta(:,i12,j12,k2) ) - zetaAux(:,:,4)       & ! Volumetric contribution: zeta-gradient to zeta-component of the flux
                  +   dfdGradQ_bo(:,:,3,i12,j12) * spAZeta % b(k1,BOTTOM) * spAZeta % vd(k2,BOTTOM) & ! Bottom face contribution (outer surface integral) 
                  +   dfdGradQ_to(:,:,3,i12,j12) * spAZeta % b(k1,TOP   ) * spAZeta % vd(k2,TOP   )   ! Top  face contribution (outer surface integral)
            
            MatEntries = MatEntries * e % geom % invJacobian(i12,j12,k1) ! Scale with Jacobian from mass matrix
            
            baseRow = i12*NCONS + j12*EtaSpa + k1*ZetaSpa
            baseCol = i12*NCONS + j12*EtaSpa + k2*ZetaSpa
            
            !-------------------------------------
            do eq2 = 1, NCONS ; do eq1 = 1, NCONS 
               i = eq1 + baseRow  ! row index (1-based)
               j = eq2 + baseCol  ! column index (1-based)
               
               call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
               
            end do            ; end do
            !-------------------------------------
         end do                  ; end do                   ; end do                   ; end do
!
!        (di) contribution
!        *****************
         do k2 = 0, e % Nxyz(3) ; do k1 = 0, e % Nxyz(3)
            do j2 = 0, e % Nxyz(2) ; do j1 = 0, e % Nxyz(2)
               do i12 = 0, e % Nxyz(1)
                  MatEntries = 0._RP
                  
                  Gvec_eta  => dF_dgradQ(:,:,:,2,i12,j2,k1)
                  Gvec_zeta => dF_dgradQ(:,:,:,3,i12,j1,k2)
                  
!
!                 Contributions to the Eta-component of the flux
!                 **********************************************
                  MatEntries = MatEntries &
!
!                    Volumetric contribution: zeta-gradient to eta-component of the flux
!                    -------------------------------------------------------------------
                        + spAEta  % hatD(j1,j2) * e % geom % invJacobian(i12,j2,k1) *                                                  &
                           (  spAZeta % hatD(k1,k2) * dot_product( Gvec_eta , e % geom % jGradZeta(:,i12,j2,k2) )                      &
                            - (  spAZeta % b(k1,BOTTOM) * spAZeta % v(k2,BOTTOM) * matmul ( dot_product( Gvec_eta , nO(:,i12,j2) ), JacO(:,:,i12,j2) ) &
                               + spAZeta % b(k1,TOP   ) * spAZeta % v(k2,TOP   ) * matmul ( dot_product( Gvec_eta , nT(:,i12,j2) ), JacT(:,:,i12,j2) ) &
                              )                                                                                                  &
                           )                                                                                                              &
                        
                        +   dfdGradQ_bo(:,:,2,i12,j1) * spAZeta % b(k1,BOTTOM) * spAEta  % D(j1,j2) * spAZeta % v(k2,BOTTOM) & ! Bottom face outer surface integral
                        +   dfdGradQ_to(:,:,2,i12,j1) * spAZeta % b(k1,TOP   ) * spAeta  % D(j1,j2) * spAZeta % v(k2,TOP   )   ! Top face outer surface integral
!
!                 Contributions to the Zeta-component of the flux
!                 ***********************************************
                  MatEntries = MatEntries &
!
!                    Volumetric contribution: eta-gradient to zeta-component of the flux
!                    -------------------------------------------------------------------
                        + spAZeta % hatD(k1,k2) * e % geom % invJacobian(i12,j1,k2) *                                                  &
                           (  spAEta  % hatD(j1,j2) * dot_product( Gvec_zeta, e % geom % jGradEta (:,i12,j2,k2) )                      &
                            - (  spAEta  % b(j1,FRONT ) * spAEta  % v(j2,FRONT ) * matmul( dot_product( Gvec_zeta, nF(:,i12,k2) ), JacF(:,:,i12,k2) ) &
                               + spAEta  % b(j1,BACK  ) * spAEta  % v(j2,BACK  ) * matmul( dot_product( Gvec_zeta, nB(:,i12,k2) ), JacB(:,:,i12,k2) ) &
                              )                                                                                                  &
                           )                                                                                                              &
                        
                        +   dfdGradQ_fr(:,:,3,i12,k1) * spAEta  % b(j1,FRONT ) * spAZeta % D(k1,k2) * spAEta  % v(j2,FRONT ) & ! Front face outer surface integral
                        +   dfdGradQ_ba(:,:,3,i12,k1) * spAEta  % b(j1,BACK  ) * spAZeta % D(k1,k2) * spAEta  % v(j2,BACK  )   ! Back face outer surface integral
               
               
                  MatEntries = MatEntries * e % geom % invJacobian(i12,j1,k1) ! Scale with Jacobian from mass matrix
                  
                  baseRow = i12*NCONS + j1*EtaSpa + k1*ZetaSpa
                  baseCol = i12*NCONS + j2*EtaSpa + k2*ZetaSpa
                  
                  !-------------------------------------
                  do eq2 = 1, NCONS ; do eq1 = 1, NCONS 
                     i = eq1 + baseRow  ! row index (1-based)
                     j = eq2 + baseCol  ! column index (1-based)
                     
                     call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
                     
                  end do            ; end do
                  !-------------------------------------
               end do
            end do                 ; end do
         end do                 ; end do
!
!        (dj) contribution
!        *****************
         do k2 = 0, e % Nxyz(3) ; do k1 = 0, e % Nxyz(3)
            do j12 = 0, e % Nxyz(2)
               do i2 = 0, e % Nxyz(1) ; do i1 = 0, e % Nxyz(1)
                  MatEntries = 0._RP
                  
                  Gvec_xi   => dF_dgradQ(:,:,:,1,i2,j12,k1)
                  Gvec_zeta => dF_dgradQ(:,:,:,3,i1,j12,k2)
               
!
!                 Contributions to the Xi-component of the flux
!                 *********************************************
                  MatEntries = MatEntries &
!
!                    Volumetric contribution: zeta-gradient to xi-component of the flux
!                    ------------------------------------------------------------------
                        + spAXi   % hatD(i1,i2) * e % geom % invJacobian(i2,j12,k1) *                                                  &
                           (  spAZeta % hatD(k1,k2) * dot_product( Gvec_xi  , e % geom % jGradZeta(:,i2,j12,k2) )                      & 
                            - (  spAZeta % b(k1,BOTTOM) * spAZeta % v(k2,BOTTOM) * matmul( dot_product( Gvec_xi  , nO(:,i2,j12) ), JacO(:,:,i2,j12) ) &
                               + spAZeta % b(k1,TOP   ) * spAZeta % v(k2,TOP   ) * matmul( dot_product( Gvec_xi  , nT(:,i2,j12) ), JacT(:,:,i2,j12) ) &
                              )                                                                                                  &
                           )                                                                                                              &
                        
                        +   dfdGradQ_bo(:,:,1,i1,j12) * spAZeta % b(k1,BOTTOM) * spAXi   % D(i1,i2) * spAZeta % v(k2,BOTTOM) & ! Bottom face outer surface integral
                        +   dfdGradQ_to(:,:,1,i1,j12) * spAZeta % b(k1,TOP   ) * spAXi   % D(i1,i2) * spAZeta % v(k2,TOP   )   ! Top face outer surface integral
!
!                 Contributions to the Zeta-component of the flux
!                 ***********************************************
                  MatEntries = MatEntries &
!
!                    Volumetric contribution: xi-gradient to zeta-component of the flux
!                    ------------------------------------------------------------------
                        + spAZeta % hatD(k1,k2) * e % geom % invJacobian(i1,j12,k2) *                                                  &
                           (  spAXi   % hatD(i1,i2) * dot_product( Gvec_zeta, e % geom % jGradXi  (:,i2,j12,k2) )                      &
                            - (  spAXi   % b(i1,LEFT  ) * spAXi   % v(i2,LEFT  ) * matmul( dot_product( Gvec_zeta, nL(:,j12,k2) ), JacL(:,:,j12,k2) ) &
                               + spAXi   % b(i1,RIGHT ) * spAXi   % v(i2,RIGHT ) * matmul( dot_product( Gvec_zeta, nR(:,j12,k2) ), JacR(:,:,j12,k2) ) &
                              )                                                                                                  &
                           )                                                                                                              &
                        
                        +   dfdGradQ_ri(:,:,3,j12,k1) * spAXi   % b(i1,RIGHT ) * spAZeta % D(k1,k2) * spAXi   % v(i2,RIGHT ) & ! Right face outer surface integral
                        +   dfdGradQ_le(:,:,3,j12,k1) * spAXi   % b(i1,LEFT  ) * spAZeta % D(k1,k2) * spAXi   % v(i2,LEFT  )   ! Left  face outer surface integral
               
               
                  MatEntries = MatEntries * e % geom % invJacobian(i1,j12,k1) ! Scale with Jacobian from mass matrix
                  
                  baseRow = i1*NCONS + j12*EtaSpa + k1*ZetaSpa
                  baseCol = i2*NCONS + j12*EtaSpa + k2*ZetaSpa
                  
                  !-------------------------------------
                  do eq2 = 1, NCONS ; do eq1 = 1, NCONS 
                     i = eq1 + baseRow  ! row index (1-based)
                     j = eq2 + baseCol  ! column index (1-based)
                     
                     call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
                     
                  end do            ; end do
                  !-------------------------------------
               end do                 ; end do
            end do
         end do                 ; end do
!
!        (dk) contribution
!        *****************
         do k12 = 0, e % Nxyz(3)
            do j2 = 0, e % Nxyz(2) ; do j1 = 0, e % Nxyz(2)
               do i2 = 0, e % Nxyz(1) ; do i1 = 0, e % Nxyz(1)
                  MatEntries = 0._RP
                  
                  Gvec_xi   => dF_dgradQ(:,:,:,1,i2,j1,k12)
                  Gvec_eta  => dF_dgradQ(:,:,:,2,i1,j2,k12)
!
!                 Contributions to the Xi-component of the flux
!                 *********************************************
                  MatEntries = MatEntries &
!
!                    Volumetric contribution: eta-gradient to xi-component of the flux
!                    ------------------------------------------------------------------
                        + spAXi   % hatD(i1,i2) * e % geom % invJacobian(i2,j1,k12) *                                                  &
                           (  spAEta  % hatD(j1,j2) * dot_product( Gvec_xi  , e % geom % jGradEta (:,i2,j2,k12) )                      &
                            - (  spAEta  % b(j1,FRONT ) * spAEta  % v(j2,FRONT ) * matmul( dot_product( Gvec_xi  , nF(:,i2,k12) ), JacF(:,:,i2,k12) ) &
                               + spAEta  % b(j1,BACK  ) * spAEta  % v(j2,BACK  ) * matmul( dot_product( Gvec_xi  , nB(:,i2,k12) ), JacB(:,:,i2,k12) ) &
                              )                                                                                                  &
                           )                                                                                                              &
                           
                        +   dfdGradQ_fr(:,:,1,i1,k12) * spAEta  % b(j1,FRONT ) * spAXi   % D(i1,i2) * spAEta  % v(j2,FRONT ) & ! Front face outer surface integral
                        +   dfdGradQ_ba(:,:,1,i1,k12) * spAEta  % b(j1,BACK  ) * spAXi   % D(i1,i2) * spAEta  % v(j2,BACK  )   ! Back  face outer surface integral
!
!                 Contributions to the Eta-component of the flux
!                 **********************************************
                  MatEntries = MatEntries &
!
!                    Volumetric contribution: xi-gradient to eta-component of the flux
!                    -----------------------------------------------------------------
                        + spAEta  % hatD(j1,j2) * e % geom % invJacobian(i1,j2,k12) *                                                  &
                           (  spAXi   % hatD(i1,i2) * dot_product( Gvec_eta , e % geom % jGradXi  (:,i2,j2,k12) )                      &
                            - (  spAXi   % b(i1,LEFT  ) * spAXi   % v(i2,LEFT  ) * matmul( dot_product( Gvec_eta , nL(:,j2,k12) ), JacL(:,:,j2,k12) ) &
                               + spAXi   % b(i1,RIGHT ) * spAXi   % v(i2,RIGHT ) * matmul( dot_product( Gvec_eta , nR(:,j2,k12) ), JacL(:,:,j2,k12) ) &
                              )                                                                                                  &
                           )                                                                                                              &
                        
                        +   dfdGradQ_ri(:,:,2,j1,k12) * spAXi   % b(i1,RIGHT ) * spAeta  % D(j1,j2) * spAXi   % v(i2,RIGHT ) & ! Right face outer surface integral
                        +   dfdGradQ_le(:,:,2,j1,k12) * spAXi   % b(i1,LEFT  ) * spAeta  % D(j1,j2) * spAXi   % v(i2,LEFT  )   ! Left face outer surface integral
               
                  MatEntries = MatEntries * e % geom % invJacobian(i1,j1,k12) ! Scale with Jacobian from mass matrix
                  
                  baseRow = i1*NCONS + j1*EtaSpa + k12*ZetaSpa
                  baseCol = i2*NCONS + j2*EtaSpa + k12*ZetaSpa
                  
                  !-------------------------------------
                  do eq2 = 1, NCONS ; do eq1 = 1, NCONS 
                     i = eq1 + baseRow  ! row index (1-based)
                     j = eq2 + baseCol  ! column index (1-based)
                     
                     call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
                     
                  end do            ; end do
                  !-------------------------------------
               end do                 ; end do
            end do                 ; end do
         end do
            
         nullify (dfdGradQ_fr, dfdGradQ_ba, dfdGradQ_bo, dfdGradQ_to, dfdGradQ_ri, dfdGradQ_le )
      end if
      
!
!     adding IBM contributes
!     ---------------------- 
      if( present(IBM) .and. IBM% active ) then
         do k12 = 0, e% Nxyz(3); do j12 = 0, e% Nxyz(2) ; do i12 = 0, e% Nxyz(1)
            if( e% isInsideBody(i12,j12,k12) ) then
               associate( Q => e% storage% Q(:,i12,j12,k12) )
               call IBM% semiImplicitJacobian( e% GlobID, Q, MatEntries )

               baseRow = i12*NCONS + j12*EtaSpa + k12*ZetaSpa
               baseCol = i12*NCONS + j12*EtaSpa + k12*ZetaSpa
 
               ! IBM contributes affect only the diagonal of the block diag. matrix
               !-------------------------------------
               do eq2 = 1, NCONS ; do eq1 = 1, NCONS 
                  i = eq1 + baseRow  ! row index (1-based)
                  j = eq2 + baseCol  ! column index (1-based)
                  call Matrix % AddToBlockEntry (e % GlobID, e % GlobID, i, j, MatEntries(eq1,eq2))
               end do; end do

               end associate
            end if
         end do; end do; end do
      end if
      
      nullify(dF_dgradQ)
      nullify(spAxi,spAeta,spAzeta)
   end subroutine Local_SetDiagonalBlock
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------
!  Local_GetOffDiagonalBlock:
!     Routine to compute the off-diagonal block that results of connecting e_plus with e_minus for 
!     the equations that correspond to e_plus and inserting it in the Matrix
!
!     -> Currently, this routine only works for p-conforming representations
!           TODO: Add routine for p-nonconforming representations (block computation is more expensive and delta cycling is ruled out)
!  -----------------------------------------------------------------------------------------------
   subroutine Local_GetOffDiagonalBlock(f,e_plus,e_minus,side,Matrix)
      implicit none
      !-arguments----------------------------------------------------------------------
      type(Face), target, intent(in)    :: f       !<  Face connecting elements
      type(Element)     , intent(in)    :: e_plus  !<  The off-diagonal block is the contribution to this element's equations
      type(Connectivity), intent(in)    :: e_minus !<  Element that connects with e_plus through "f"
      integer           , intent(in)    :: side    !<  side of face where e_plus is
      class(Matrix_t)   , intent(inout) :: Matrix  !<> Jacobian matrix
      !-local-variables----------------------------------------------------------------
      integer :: i, j                     ! Matrix indexes
      integer :: i1, j1, k1, eq1          ! variable counters
      integer :: i2, j2, k2, eq2          ! variable counters
      integer :: baseRow, baseCol         ! Position of NCONS by NCONS miniblock of Jacobian
      integer :: r                        ! Additional index for counting in the normal direction (only needed for viscous fluxes)
      integer :: nXi1, nEta1              ! Number of nodes in every direction
      integer :: EtaSpa1, ZetaSpa1        ! Spacing for these two coordinate directions
      integer :: nXi2, nEta2              ! Number of nodes in every direction
      integer :: EtaSpa2, ZetaSpa2        ! Spacing for these two coordinate directions
      integer :: elSide_plus              ! Element side where f is on e⁻
      integer :: elSide_minus             ! Element side where f is on e⁺
      integer :: elInd_plus(3)            ! Element indexes on e⁺
      integer :: elInd_minus(3)           ! Element indexes on e⁻
      integer :: normAx_plus              ! Normal axis to f on e⁺
      integer :: normAx_minus             ! Normal axis to f on e⁻
      integer :: tanAx_plus (2)           ! Tangent axes to f on e⁺
      integer :: tanAx_minus(2)           ! Tangent axes to f on e⁻
      integer :: normAxSide_plus          ! Side of the normal axis that is in contact with f on e⁺
      integer :: normAxSide_minus         ! Side of the normal axis that is in contact with f on e⁻
      integer :: faceInd_plus (2)         ! Face indexes on e⁺
      integer :: faceInd_minus(2)         ! Face indexes on e⁻
      integer :: NxyFace_plus (2)         ! Polynomial orders of face on element e⁺
      integer :: NxyFace_minus(2)         ! Polynomial orders of face on element e⁻ (only needed for viscous fluxes)
      integer :: faceInd_plus2minus(2)    ! Face indexes on e⁺ passed to the reference frame of e⁻
      integer :: faceInd_minus2plus(2)    ! Face indexes on e⁻ passed to the reference frame of e⁺ (only needed for viscous fluxes)
      integer :: Deltas                   ! Number of Kronecker deltas /= 0
      integer :: elInd_plus_norm(3)       ! Element indexes on e⁺, but the index in the normal direction has been replaced by a specified value "r"
      integer :: elInd_plus_tan1(3)       ! Element indexes on e⁺, but the index in the first tangent direction has been replaced by the index on e⁻ (in the reference frame of e⁺)
      integer :: elInd_plus_tan2(3)       ! Element indexes on e⁺, but the index in the second tangent direction has been replaced by the index on e⁻ (in the reference frame of e⁺)
      real(kind=RP) :: MatEntries(NCONS,NCONS)           ! Values of the matrix entries
      type(NodalStorage_t), target  :: spA_plus (3)      ! Nodal storage in the different directions for e_plus  - local copy
      type(NodalStorage_t), target  :: spA_minus(3)      ! Nodal storage in the different directions for e_minus - local copy
      type(NodalStorage_t), pointer :: spAnorm_plus      ! Nodal storage in the direction that is normal to the face for e⁺
      type(NodalStorage_t), pointer :: spAtan1_plus      ! Nodal storage in the tangent direction "1" to the face for e⁺ (only needed for viscous fluxes)
      type(NodalStorage_t), pointer :: spAtan2_plus      ! Nodal storage in the tangent direction "2" to the face for e⁺ (only needed for viscous fluxes)
      type(NodalStorage_t), pointer :: spAnorm_minus     ! Nodal storage in the direction that is normal to the face for e⁻
      type(NodalStorage_t), pointer :: spAtan1_minus     ! Nodal storage in the tangent direction "1" to the face for e⁻ (only needed for viscous fluxes)
      type(NodalStorage_t), pointer :: spAtan2_minus     ! Nodal storage in the tangent direction "2" to the face for e⁻ (only needed for viscous fluxes)
      real(kind=RP)       , pointer :: dfdq(:,:,:,:)     ! 
      real(kind=RP)       , pointer :: df_dGradQ_f(:,:,:,:,:)     ! Pointer to the Jacobian with respect to gradQ on face
      real(kind=RP)       , pointer :: df_dGradQ_e(:,:,:,:,:,:,:) ! Pointer to the Jacobian with respect to gradQ on face
      real(kind=RP) :: dtan1, dtan2             ! Kronecker deltas in the tangent directions
      real(kind=RP) :: dtan1_minus, dtan2_minus ! Kronecker deltas in the tangent directions in the reference frame of e⁻ (only needed for viscous fluxes)
      real(kind=RP) :: a_minus
      real(kind=RP) :: normAux(NCONS,NCONS)
      real(kind=RP), pointer :: Gvec_norm(:,:,:)       ! Auxiliary vector containing values of dFv_dgradQ in the direction normal to the face
      real(kind=RP), pointer :: Gvec_tan1(:,:,:)       ! Auxiliary vector containing values of dFv_dgradQ in the first tangent direction to the face
      real(kind=RP), pointer :: Gvec_tan2(:,:,:)       ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
      real(kind=RP), allocatable :: nHat(:,:,:)
      !--------------------------------------------------------------------------------
!
!     ***********
!     Definitions
!     ***********
!
      a_minus = 0.5_RP  ! Temp... TODO: read from ViscousDiscretization
      
      ! Entry spacing for element e⁺
      nXi1     = e_plus % Nxyz(1) + 1
      nEta1    = e_plus % Nxyz(2) + 1
      EtaSpa1  = NCONS*nXi1
      ZetaSpa1 = NCONS*nXi1*nEta1
      
      ! Entry spacing for element e⁻
      nXi2     = e_minus % Nxyz(1) + 1
      nEta2    = e_minus % Nxyz(2) + 1
      EtaSpa2  = NCONS*nXi2
      ZetaSpa2 = NCONS*nXi2*nEta2
      
      ! Element sides
      elSide_plus  = f % elementSide(side)
      elSide_minus = f % elementSide(other(side))
      
      ! Normal and tangent axes
      normAx_plus  = normalAxis (elSide_plus)
      normAx_minus = normalAxis (elSide_minus)
      tanAx_plus   = axisMap(:,f % elementSide(    side     ) )
      tanAx_minus  = axisMap(:,f % elementSide( other(side) ) )
      
      ! Side of axis where f is
      if (normAx_plus < 0) then
         normAxSide_plus = LEFT
      else
         normAxSide_plus = RIGHT
      end if
      normAx_plus = abs(normAx_plus)
      if (normAx_minus < 0) then
         normAxSide_minus = LEFT
      else
         normAxSide_minus = RIGHT
      end if
      normAx_minus = abs(normAx_minus)
      
      ! Nodal storage
      ! --------------------------------------
      ! TODO: Why this doesn't work since ifort ver. 19.1?
      ! --------------------------------------
      ! spA_plus  = NodalStorage(e_plus  % Nxyz)
      ! spA_minus = NodalStorage(e_minus % Nxyz)
      
      spA_plus(1)  = NodalStorage(e_plus  % Nxyz(1))
      spA_plus(2)  = NodalStorage(e_plus  % Nxyz(2))
      spA_plus(3)  = NodalStorage(e_plus  % Nxyz(3))
      spA_minus(1)  = NodalStorage(e_minus  % Nxyz(1))
      spA_minus(2)  = NodalStorage(e_minus  % Nxyz(2))
      spA_minus(3)  = NodalStorage(e_minus  % Nxyz(3))

      spAnorm_plus  => spA_plus( normAx_plus   )
      spAtan1_plus  => spA_plus( tanAx_plus(1) )
      spAtan2_plus  => spA_plus( tanAx_plus(2) )
      
      spAnorm_minus => spA_minus( normAx_minus   )
      spAtan1_minus => spA_minus( tanAx_minus(1) )
      spAtan2_minus => spA_minus( tanAx_minus(2) )
      
      ! Polynomial orders
      NxyFace_plus  = e_plus  % Nxyz ( tanAx_plus  )
      NxyFace_minus = e_minus % Nxyz ( tanAx_minus )
      
!     Get normal vectors in right element indexing (e⁻)
!     -------------------------------------------------
      allocate ( nHat(3,0:NxyFace_minus(1),0:NxyFace_minus(2)) )
      call AnJac_GetNormalVec(f, other(side), nHat)
      
!
!     *********************
!     Inviscid contribution
!     *********************
!
      
!
!     Pointers to the flux Jacobians with respect to q on the faces
!     -------------------------------------------------------------
      dfdq(1:,1:,0:,0:) => f % storage(side) % dFStar_dqEl(1:,1:,0:,0:,other(side))
      
      do k2 = 0, e_minus % Nxyz(3) ; do j2 = 0, e_minus % Nxyz(2) ; do i2 = 0, e_minus % Nxyz(1)
         do k1 = 0, e_plus % Nxyz(3) ; do j1 = 0, e_plus % Nxyz(2) ; do i1 = 0, e_plus % Nxyz(1) 
            
            elInd_plus  = [i1, j1, k1]
            elInd_minus = [i2, j2, k2]
            
            faceInd_plus  = elInd_plus ( tanAx_plus  )
            faceInd_minus = elInd_minus( tanAx_minus )
            
            call indexesOnOtherFace(faceInd_plus(1),faceInd_plus(2), NxyFace_plus(1),NxyFace_plus(2), f % rotation, side, faceInd_plus2minus(1),faceInd_plus2minus(2))
            
            ! "delta" cycling
            if ( any(faceInd_plus2minus /= faceInd_minus) ) cycle
            
            MatEntries = -    dfdq(:,:,faceInd_plus(1),faceInd_plus(2)) &
                            * spAnorm_plus  % b(elInd_plus ( normAx_plus  ), normAxSide_plus ) &
                            * spAnorm_minus % v(elInd_minus( normAx_minus ), normAxSide_minus )  & 
                            * e_plus % geom % invJacobian(i1,j1,k1)
            
            baseRow = i1*NCONS + j1*EtaSpa1 + k1*ZetaSpa1
            baseCol = i2*NCONS + j2*EtaSpa2 + k2*ZetaSpa2
            
            do eq2 = 1, NCONS ; do eq1 = 1, NCONS
               i = eq1 + baseRow ! row index (1-based)
               j = eq2 + baseCol ! column index (1-based)
               
               call Matrix % AddToBlockEntry (e_plus % globID, e_minus % globID, i, j, MatEntries(eq1,eq2))
            end do            ; end do
            
         end do                ; end do                ; end do
      end do                ; end do                ; end do
      nullify(dfdq)
!
!     ********************
!     Viscous contribution
!     ********************
!
      if (flowIsNavierStokes) then
!
!        Pointers to the flux Jacobians with respect to grad(q) on the faces
!        -------------------------------------------
         df_dGradQ_f(1:,1:,1:,0:,0:)    => f % Storage(side) % dFv_dGradQEl (1:,1:,1:,0:,0:, other(side) )
         df_dGradQ_e(1:,1:,1:,1:,0:,0:,0:) => e_plus % storage  % dF_dGradQ    (1:,1:,1:,1:,0:,0:,0:)
         do k2 = 0, e_minus % Nxyz(3) ; do j2 = 0, e_minus % Nxyz(2) ; do i2 = 0, e_minus % Nxyz(1)
            do k1 = 0, e_plus % Nxyz(3) ; do j1 = 0, e_plus % Nxyz(2) ; do i1 = 0, e_plus % Nxyz(1)
               
               elInd_plus  = [i1, j1, k1]
               elInd_minus = [i2, j2, k2]
               
               faceInd_plus  = elInd_plus ( tanAx_plus  )
               faceInd_minus = elInd_minus( tanAx_minus )
               
               call indexesOnOtherFace(faceInd_plus (1),faceInd_plus (2), NxyFace_plus (1),NxyFace_plus (2), f % rotation, side       , faceInd_plus2minus(1),faceInd_plus2minus(2))
               call indexesOnOtherFace(faceInd_minus(1),faceInd_minus(2), NxyFace_minus(1),NxyFace_minus(2), f % rotation, other(side), faceInd_minus2plus(1),faceInd_minus2plus(2))
               
               elInd_plus_tan1 = elInd_plus
               elInd_plus_tan1(tanAx_plus(1)) = faceInd_minus2plus(1)
               
               elInd_plus_tan2 = elInd_plus
               elInd_plus_tan2(tanAx_plus(2)) = faceInd_minus2plus(2)
               
!              Delta computation on "plus" side (used for delta cycling)
!              ---------------------------------------------------------
               Deltas = 0
               if ( faceInd_plus(1) == faceInd_minus2plus(1) ) then
                  dtan1 = 1._RP
                  Deltas = Deltas + 1
               else
                  dtan1 = 0._RP
               end if
               if ( faceInd_plus(2) == faceInd_minus2plus(2) ) then
                  dtan2 = 1._RP
                  Deltas = Deltas + 1
               else
                  dtan2 = 0._RP
               end if
               
               ! delta cycling
               if (Deltas < 1) cycle
               
!              Delta computation on "minus" side
!              ---------------------------------
               if ( faceInd_minus(1) == faceInd_plus2minus(1) ) then
                  dtan1_minus = 1._RP
               else
                  dtan1_minus = 0._RP
               end if
               if ( faceInd_minus(2) == faceInd_plus2minus(2) ) then
                  dtan2_minus = 1._RP
               else
                  dtan2_minus = 0._RP
               end if
               
!              Computation of the matrix entry
!              -------------------------------
               
!              Get NORMAL auxiliary term
!              ------------------------
               normAux = 0._RP
               if (dtan1*dtan2 > 0.5_RP) then
                  do r=0, e_plus  % Nxyz ( normAx_plus  )
                     
                     elInd_plus_norm = elInd_plus
                     elInd_plus_norm(normAx_plus) = r
                     
                     Gvec_norm => df_dGradQ_e(:,:,:,normAx_plus,elInd_plus_norm(1),elInd_plus_norm(2),elInd_plus_norm(3))
                     
                     normAux  = normAux  + &
                                e_plus % geom % invJacobian(elInd_plus_norm(1),elInd_plus_norm(2),elInd_plus_norm(3)) &
                              * spAnorm_plus  % hatD(elInd_plus ( normAx_plus  ),r) &
                              * spAnorm_plus  % b(r,normAxSide_plus ) &
                              * spAnorm_minus % v(elInd_minus( normAx_minus ),normAxSide_minus) &
                              * dot_product( Gvec_norm , nHat(:,faceInd_minus(1),faceInd_minus(2)) )
                  end do
                  normAux = normAux * a_minus
               end if
               
               Gvec_tan1 => df_dGradQ_e(:,:,:,tanAx_plus(1),elInd_plus_tan1(1),elInd_plus_tan1(2),elInd_plus_tan1(3))
               Gvec_tan2 => df_dGradQ_e(:,:,:,tanAx_plus(2),elInd_plus_tan2(1),elInd_plus_tan2(2),elInd_plus_tan2(3))
               
               MatEntries = 0._RP
!
!              Volumetric contribution (inner fluxes)
!              **************************************
!                    Normal direction:
               MatEntries = MatEntries + dtan1 * dtan2 * normAux
!                    Tangent direction 1:
               MatEntries = MatEntries + dtan2 * a_minus                                                           &
                     * e_plus % geom % invJacobian(elInd_plus_tan1(1),elInd_plus_tan1(2),elInd_plus_tan1(3))   &
                     * spAtan1_plus  % hatD(faceInd_plus(1),faceInd_minus2plus(1))                             &
                     * spAnorm_plus  % b(elInd_plus ( normAx_plus  ),normAxSide_plus )                         &
                     * spAnorm_minus % v(elInd_minus( normAx_minus ),normAxSide_minus)                         &
                     * dot_product( Gvec_tan1 , nHat(:,faceInd_minus(1),faceInd_minus(2)) ) 
!                    Tangent direction 2:
               MatEntries = MatEntries + dtan1 * a_minus                                                           &
                     * e_plus % geom % invJacobian(elInd_plus_tan2(1),elInd_plus_tan2(2),elInd_plus_tan2(3))   &
                     * spAtan2_plus  % hatD(faceInd_plus(2),faceInd_minus2plus(2))                             &
                     * spAnorm_plus  % b(elInd_plus ( normAx_plus  ),normAxSide_plus )                         &
                     * spAnorm_minus % v(elInd_minus( normAx_minus ),normAxSide_minus)                         &
                     * dot_product( Gvec_tan2, nHat(:,faceInd_minus(1),faceInd_minus(2)) )
!
!              Faces contribution (surface integrals from the outer equation) - PENALTY TERM IS BEING CONSIDERED IN THE INVISCID PART - TODO: Reorganize storage to put it explicitly in another place (needed for purely viscous equations)
!                 The tangent directions here are taken in the reference frame of e⁻
!              *********************************************************************
!
!
               MatEntries = MatEntries &
                   +   df_dGradQ_f(:,:,tanAx_minus(1),faceInd_plus(1),faceInd_plus(2))                               &
                     * spAnorm_plus  % b   (elInd_plus ( normAx_plus  ), normAxSide_plus  )                             &
                     * spAtan1_minus % D   (faceInd_plus2minus(1),faceInd_minus(1))                                     &
                     * spAnorm_minus % v   (elInd_minus( normAx_minus ), normAxSide_minus ) * dtan2_minus               & ! Tangent direction 1
                     
                   +   df_dGradQ_f(:,:, normAx_minus ,faceInd_plus(1),faceInd_plus(2))                               &
                     * spAnorm_plus  % b   (elInd_plus ( normAx_plus  ), normAxSide_plus  )                             &
                     * spAnorm_minus % vd  (elInd_minus( normAx_minus ), normAxSide_minus ) * dtan1_minus * dtan2_minus & ! Normal direction
                     
                   +   df_dGradQ_f(:,:,tanAx_minus(2),faceInd_plus(1),faceInd_plus(2))                               &
                     * spAnorm_plus  % b   (elInd_plus ( normAx_plus  ), normAxSide_plus  )                             &
                     * spAtan2_minus % D   (faceInd_plus2minus(2),faceInd_minus(2))                                     &
                     * spAnorm_minus % v   (elInd_minus( normAx_minus ), normAxSide_minus ) * dtan1_minus                 ! Tangent direction 2 
               
               MatEntries = MatEntries * e_plus % geom % invJacobian(i1,j1,k1) ! Scale with Jacobian from mass matrix
               
               baseRow = i1*NCONS + j1*EtaSpa1 + k1*ZetaSpa1
               baseCol = i2*NCONS + j2*EtaSpa2 + k2*ZetaSpa2
               
               do eq2 = 1, NCONS ; do eq1 = 1, NCONS 
                  i = eq1 + baseRow ! row index (1-based)
                  j = eq2 + baseCol ! column index (1-based)
               
                  call Matrix % AddToBlockEntry (e_plus % globID, e_minus % globID, i, j, MatEntries(eq1,eq2))
               end do            ; end do
               
            end do                ; end do                ; end do
         end do                ; end do                ; end do
         nullify(df_dGradQ_f)
         nullify(df_dGradQ_e)
      end if
!
!     *********
!     Finish up
!     *********
!
      nullify(spAnorm_plus, spAtan1_plus, spAtan2_plus, spAnorm_minus, spAtan1_minus, spAtan2_minus)
      
   end subroutine Local_GetOffDiagonalBlock
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------
!  AnJac_GetNormalVec: 
!  Get normal vector for element indexing and scale with surface Jacobian
!     (Only for p-conforming meshes)
!  ----------------------------------------------------------------------
   subroutine AnJac_GetNormalVec(f, faceSide, nEl)
      implicit none
      !-arguments-----------------------------------------------------
      type(Face)     , intent(in)    :: f
      integer        , intent(in)    :: faceSide
      real(kind=RP)  , intent(inout) :: nEl (:,0:,0:)
      !-local-variables-----------------------------------------------
      integer :: i,j, ii, jj
      !---------------------------------------------------------------
      
      select case (faceSide)
         case (LEFT)
            do j = 0, f % Nf(2)   ; do i = 0, f % Nf(1)
               nEl(:,i,j) = f % geom % normal(:,i,j) * f % geom % jacobian(i,j)
            end do                        ; end do
         case (RIGHT)
            do j = 0, f % Nf(2)   ; do i = 0, f % Nf(1)
               call leftIndexes2Right(i,j,f % Nf(1), f % Nf(2), f % rotation, ii, jj)
               
               nEl(1:NDIM,ii,jj) = - f % geom % normal (1:NDIM,i,j) * f % geom % jacobian(i,j)
            end do                        ; end do
      end select
      
   end subroutine AnJac_GetNormalVec  
   
   !
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------
!  AnJac_GetFaceJac: 
!  Get face Jacobian... This routine is used for the surface integrals of the grad (inner) equation
!  Perform the matrix operation only on boundaries....
!  ----------------------------------------------------------------------
   subroutine AnJac_GetFaceJac(f, faceSide, faceJac, dqHat_dqp, dqHat_dqm)
      implicit none
      !-arguments-----------------------------------------------------
      type(Face)     , intent(in)    :: f
      integer        , intent(in)    :: faceSide
      real(kind=RP)  , intent(inout) :: faceJac (:,:,0:,0:)
      real(kind=RP)  , intent(in)    :: dqHat_dqp, dqHat_dqm
      !-local-variables-----------------------------------------------
      integer :: i,j, ii, jj
      !---------------------------------------------------------------
      
      select case (faceSide)
         case (LEFT)
            if (f % FaceType == HMESH_BOUNDARY) then
               do j = 0, f % NelLeft(2)   ; do i = 0, f % NelLeft(1)
                  faceJac(:,:,i,j) = Identity * dqHat_dqp + dqHat_dqm * f % storage(1) % BCJac(:,:,i,j)
               end do                        ; end do
            else
               do j = 0, f % NelLeft(2)   ; do i = 0, f % NelLeft(1)
                  faceJac(:,:,i,j) = Identity * dqHat_dqp
               end do                        ; end do
            end if
         case (RIGHT)
            do j = 0, f % NelRight(2)   ; do i = 0, f % NelRight(1)
               faceJac(:,:,i,j) = Identity * dqHat_dqp
            end do                        ; end do
      end select
      
   end subroutine AnJac_GetFaceJac  
   
#endif
end module AnalyticalJacobian
