!//////////////////////////////////////////////////////
!
!     TODO1: Store FaceStorage in SolutionStorage
!     TODO2: Remove physics-related pointers... Allocate one storage for each physics (AND REMOVE ALL POINTERS)
!     TODO3: Allocate implicit methods storage (global arrays) from implicit classes
#include "Includes.h"

module StorageClass
   use, intrinsic :: iso_c_binding
   use SMConstants
   use PhysicsStorage
   use InterpolationMatrices, only: Interp3DArrays, Tset
   use NodalStorageClass    , only: NodalStorage
   implicit none

   private
   public   ElementStorage_t, FaceStorage_t, SolutionStorage_t
   public   GetStorageEquations

   enum, bind(C)
      enumerator :: OFF = 0, NS, C, MU, NSSA, CAA
   end enum

   type Statistics_t
      real(kind=RP), dimension(:,:,:,:),  allocatable    :: data
      contains
         procedure   :: Construct => Statistics_Construct
         procedure   :: Destruct  => Statistics_Destruct
   end type Statistics_t

!
!  Class for pointing to previous solutions in an element
!  ******************************************************
   type ElementPrevSol_t
      real(kind=RP), dimension(:,:,:,:),  pointer, contiguous     :: Q ! Pointers to the appropriate storage (NS or CH)
#ifdef FLOW
      real(kind=RP), dimension(:,:,:,:),  allocatable     :: QNS
#endif
#ifdef CAHNHILLIARD
      real(kind=RP), dimension(:,:,:,:),  allocatable     :: c
#endif
   end type ElementPrevSol_t

   type RKStep_t
      real(kind=RP), dimension(:,:,:,:), allocatable :: K      ! Explicit Runge-Kutta coefficients
      real(kind=RP), dimension(:,:,:,:), allocatable :: hatK   ! Implicit Runge-Kutta coefficients
   end type RKStep_t
!
!  Class for storing variables element-wise
!  ****************************************
   type ElementStorage_t
      logical                                         :: computeGradients
      logical                                         :: anJacobian = .FALSE.    ! Does we need analytical Jacobian specific storage?
      integer                                         :: prevSol_num
      integer                                         :: RKSteps_num
      integer                                         :: currentlyLoaded
      integer                                         :: NDOF              ! Number of degrees of freedom of element
      integer                                         :: Nxyz(NDIM)
      real(kind=RP)                                   :: min_lcl_dst       ! Minimum local distance between nodal points for CFL calculation (physical space)
      integer                                         :: first_sensed      ! Time-steps since the element was first sensed
      real(kind=RP)                                   :: prev_sensor       ! Previous value of the sensor
      real(kind=RP)                                   :: sensor            ! Value of the sensor
      real(kind=RP), dimension(:,:,:,:),  pointer, contiguous     :: Q           ! Pointers to the appropriate storage (NS or CH)
      real(kind=RP), dimension(:,:,:,:),  pointer, contiguous     :: QDot        !
      real(kind=RP), dimension(:,:,:,:),  pointer, contiguous     :: U_x         !
      real(kind=RP), dimension(:,:,:,:),  pointer, contiguous     :: U_y         !
      real(kind=RP), dimension(:,:,:,:),  pointer, contiguous     :: U_z         !
      type(ElementPrevSol_t),  allocatable :: PrevQ(:)           ! Previous solution
      type(RKStep_t),          allocatable :: RKSteps(:)         ! Runge-Kutta stages
#ifdef FLOW
      real(kind=RP),           allocatable :: QNS(:,:,:,:)           ! NSE State vector
      real(kind=RP),           allocatable :: rho(:,:,:)             ! Temporal storage for the density
      real(kind=RP), private,  allocatable :: QDotNS(:,:,:,:)        ! NSE State vector time derivative
      real(kind=RP), private,  allocatable :: U_xNS(:,:,:,:)         ! NSE x-gradients
      real(kind=RP), private,  allocatable :: U_yNS(:,:,:,:)         ! NSE y-gradients
      real(kind=RP), private,  allocatable :: U_zNS(:,:,:,:)         ! NSE z-gradients
      real(kind=RP),           allocatable :: G_NS(:,:,:,:)          ! NSE auxiliary storage
      real(kind=RP),           allocatable :: S_NS(:,:,:,:)          ! NSE source term
      real(kind=RP),           allocatable :: S_NSP(:,:,:,:)         ! NSE Particles source term
      real(kind=RP),           allocatable :: QbaseSponge(:,:,:,:)   ! Base Flow State vector for sponges
      real(kind=RP),           allocatable :: intensitySponge(:,:,:) ! Intensity for sponges
      real(kind=RP),           allocatable :: QLowRK(:,:,:,:)        ! NSE State vector solved with a lower order RK method (required for adaptive time-stepping)
#ifndef ACOUSTIC
      real(kind=RP),           allocatable :: mu_NS(:,:,:,:)       ! (mu, beta, kappa) artificial
      real(kind=RP),           allocatable :: mu_turb_NS(:,:,:)    ! mu of LES
#endif
      real(kind=RP),           allocatable :: dF_dgradQ(:,:,:,:,:,:,:) ! NSE Jacobian with respect to gradQ
      type(Statistics_t)                   :: stats                ! NSE statistics
      real(kind=RP)                        :: artificialDiss
#endif
#ifdef SPALARTALMARAS
      real(kind=RP),           allocatable ::  S_SA(:,:,:,:)
      real(kind=RP),           allocatable :: mu_SA(:,:,:,:)         ! EddyViscocityVector, EddyetaVector
#endif
#ifdef ACOUSTIC
      real(kind=RP),           allocatable :: Qbase(:,:,:,:)         ! Base flow State vector
      real(kind=RP),           allocatable :: grada2base(:,:,:,:)    ! Gradient of the squared sound velocity
      real(kind=RP),           allocatable :: Lambbase(:,:,:,:)      ! Lamb vector base
      real(kind=RP),           allocatable :: Lamb(:,:,:,:)          ! Lamb vector
      real(kind=RP),           allocatable :: Lambread(:,:,:,:,:)    ! Lamb vector read at different time instants (NDIM,i,j,k,NInterpolation)
#endif
#ifdef CAHNHILLIARD
      real(kind=RP), dimension(:,:,:,:),   allocatable :: c     ! CHE concentration
      real(kind=RP), dimension(:,:,:,:),   allocatable :: cDot  ! CHE concentration time derivative
      real(kind=RP), dimension(:,:,:,:),   allocatable :: c_x   ! CHE concentration x-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: c_y   ! CHE concentration y-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: c_z   ! CHE concentration z-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: mu    ! CHE chemical potential
      real(kind=RP), dimension(:,:,:,:),   allocatable :: mu_x  ! CHE chemical potential x-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: mu_y  ! CHE chemical potential y-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: mu_z  ! CHE chemical potential z-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: v     ! CHE flow field velocity
      real(kind=RP), dimension(:,:,:,:),   allocatable :: G_CH  ! CHE auxiliary storage
#endif
#ifdef MULTIPHASE
      real(kind=RP),           allocatable :: invMa2(:,:,:)        !Storage for the density artificial compressibility factor
#endif
      contains
         procedure   :: Assign              => ElementStorage_Assign
         generic     :: assignment(=)       => Assign
         procedure   :: Construct           => ElementStorage_Construct
         procedure   :: Destruct            => ElementStorage_Destruct
         procedure   :: InterpolateSolution => ElementStorage_InterpolateSolution
         procedure   :: PointStorage        => ElementStorage_PointStorage
         procedure   :: constructAnJac      => ElementStorage_ConstructAnJac
#ifdef FLOW
         procedure   :: SetStorageToNS    => ElementStorage_SetStorageToNS
#endif
#ifdef CAHNHILLIARD
         procedure   :: SetStorageToCH_c  => ElementStorage_SetStorageToCH_c
         procedure   :: SetStorageToCH_mu => ElementStorage_SetStorageToCH_mu
#endif
   end type ElementStorage_t

!
!  Class for storing variables in the whole domain
!  ***********************************************
   type SolutionStorage_t
      integer                                    :: NDOF
      logical                                    :: AdaptedQ     = .FALSE.
      logical                                    :: AdaptedQdot  = .FALSE.
      logical                                    :: AdaptedPrevQ = .FALSE.
      logical                                    :: anJacobian   = .FALSE.
      integer                                    :: prevSol_num    = -1
      integer                      , allocatable :: prevSol_index(:)           ! Indexes for the previous solutions

      type(ElementStorage_t)       , allocatable :: elements(:)
      real(kind=RP),                 pointer     :: Q(:)
      real(kind=RP),                 pointer     :: QDot(:)
      real(kind=RP),                 pointer     :: PrevQ(:,:)
#ifdef FLOW
      real(kind=RP), dimension(:)  , allocatable :: QdotNS
      real(kind=RP), dimension(:)  , allocatable :: QNS
      real(kind=RP), dimension(:,:), allocatable :: PrevQNS ! Previous solution(s) in the whole domain
#endif
#ifdef CAHNHILLIARD
      real(kind=RP), dimension(:)  , allocatable :: cDot
      real(kind=RP), dimension(:)  , allocatable :: c
      real(kind=RP), dimension(:,:), allocatable :: Prevc(:,:)
#endif
      contains
         procedure :: construct        => SolutionStorage_Construct
         procedure :: local2GlobalQ    => SolutionStorage_local2GlobalQ
         procedure :: local2GlobalQdot => SolutionStorage_local2GlobalQdot
         procedure :: SetGlobalPrevQ   => SolutionStorage_SetGlobalPrevQ
         procedure :: global2LocalQ    => SolutionStorage_global2LocalQ
         procedure :: global2LocalQdot => SolutionStorage_global2LocalQdot
         procedure :: Destruct         => SolutionStorage_Destruct
         procedure :: SignalAdaptation => SolutionStorage_SignalAdaptation
         procedure :: PointStorage     => SolutionStorage_PointStorage
         procedure :: copy             => SolutionStorage_Assign
         generic   :: assignment(=)    => copy
   end type SolutionStorage_t
!
!  Class for storing variables in the faces
!  ****************************************
   type FaceStorage_t
      logical                                          :: constructed = .FALSE.
      logical                                          :: computeGradients
      logical                                          :: computeQdot
      logical                                          :: anJacobian =.FALSE.         ! Has Jacobian storage?
      integer                                          :: NDIM
      integer                                          :: currentlyLoaded
      integer                                          :: Nf(2), Nel(2)
      real(kind=RP), dimension(:,:,:),     pointer     :: Q
      real(kind=RP), dimension(:,:,:),     pointer     :: Qdot
      real(kind=RP), dimension(:,:,:),     pointer     :: U_x, U_y, U_z
      real(kind=RP), dimension(:,:,:),     pointer     :: FStar
      real(kind=RP), dimension(:,:,:),     allocatable :: AviscFlux
      real(kind=RP), dimension(:,:,:,:),   pointer     :: unStar
      real(kind=RP), dimension(:),         allocatable :: genericInterfaceFluxMemory ! unStar and fStar point to this memory simultaneously. This seems safe.
#ifdef FLOW
      real(kind=RP), dimension(:,:,:),     allocatable :: QNS
      real(kind=RP), dimension(:,:,:),     allocatable :: QdotNS
      real(kind=RP), dimension(:,:,:),     allocatable :: U_xNS, U_yNS, U_zNS
      real(kind=RP), dimension(:,:),       allocatable :: rho
      real(kind=RP), dimension(:,:,:),     allocatable :: mu_NS
      real(kind=RP), dimension(:,:),       allocatable :: u_tau_NS
      real(kind=RP), dimension(:,:),       allocatable :: wallNodeDistance ! for BC walls, distance to the first fluid node
#ifdef ACOUSTIC
      real(kind=RP), dimension(:,:,:),     allocatable :: Qbase ! Base flow State vector
      real(kind=RP), dimension(:,:,:),     allocatable :: grada2base ! Gradient of the squared sound velocity
      real(kind=RP), dimension(:,:,:),     allocatable :: Lambbase ! Lamb vector base
      real(kind=RP), dimension(:,:,:),     allocatable :: Lamb ! Lamb vector
      real(kind=RP), dimension(:,:,:,:),   allocatable :: Lambread ! Lamb vector read
#endif
#ifdef MULTIPHASE
      real(kind=RP), dimension(:,:),       allocatable :: invMa2
#endif
!
!     Inviscid Jacobians
!     ------------------
!     * On the face (mortar points):
      real(kind=RP), allocatable :: dFStar_dqF(:,:,:,:)
!                   storage(side) % dFStar_dqF(:,:,i,j)
!                           |                  |_| |_|
!                           |                   |   |
!                           |                   |   |__Coordinate indexes in face
!                           |                   |______Jacobian for this component
!                           |__________________________1 for dFStar/dqL and 2 for dFStar/dqR
!
!     * On the coordinates that match the element's (face-element points):
      real(kind=RP), allocatable :: dFStar_dqEl(:,:,:,:,:)
!                   storage(side) % dFStar_dqeL(:,:,i,j,:)
!                           |                   |_| |_| |
!                           |                    |   |  |_1 for dFStar/dqL and 2 for dFStar/dqR
!                           |                    |   |____Coordinate indexes in face
!                           |                    |________Jacobian for this component
!                           |_____________________________1 for element on the left, 2 for element on the right
!
!     Viscous Jacobians
!     -----------------
!     * On the face (mortar points):
      real(kind=RP), allocatable :: dFv_dGradQF(:,:,:,:,:)
!                   storage(side) % dFv_dGradQF(:,:,:,i,j)
!                           |                   |_| | |_|
!                           |                    |  |  |
!                           |                    |  |  |__Coordinate indexes in face
!                           |                    |  | ____1 for inner term, 2 for outer term
!                           |                    |  |_____∇q component: 1, 2, 3
!                           |                    |________Jacobian for this component
!                           |_____________________________1 for dFv*/d∇qL and 2 for dFv*/d∇qR
!
!     * On the coordinates that match the element's (face-element points):
      real(kind=RP), allocatable :: dFv_dGradQEl(:,:,:,:,:,:)
!                   storage(side) % dFv_dGradQEl(:,:,:,i,j,:)
!                           |                    |_| | |_| |
!                           |                     |  |  |  |_1 for dFv*/d∇qL and 2 for dFv*/d∇qR
!                           |                     |  |  |____Coordinate indexes in face
!                           |                     |  |_______∇q component: 1, 2, 3
!                           |                     |__________Jacobian for this component
!                           |________________________________1 for element on the left, 2 for element on the right
!
!     * Jacobian of the boundary condition (only needs to be stored on boundary faces for viscous physics to apply the BC to the grad equation):
      real(kind=RP), allocatable :: BCJac(:,:,:,:)
!                                         |_| |_|
!                                          |   |
!                                          |   |__Coordinate indexes in face
!                                          |______Jacobian for this component
!
#endif
!
#ifdef CAHNHILLIARD
      real(kind=RP), dimension(:,:,:),   allocatable :: c
      real(kind=RP), dimension(:,:,:),   allocatable :: c_x
      real(kind=RP), dimension(:,:,:),   allocatable :: c_y
      real(kind=RP), dimension(:,:,:),   allocatable :: c_z
      real(kind=RP), dimension(:,:,:),   allocatable :: mu
      real(kind=RP), dimension(:,:,:),   allocatable :: mu_x
      real(kind=RP), dimension(:,:,:),   allocatable :: mu_y
      real(kind=RP), dimension(:,:,:),   allocatable :: mu_z
      real(kind=RP), dimension(:,:,:),   allocatable :: v
#endif
      contains
         procedure   :: Construct      => FaceStorage_Construct
         procedure   :: Destruct       => FaceStorage_Destruct
         procedure   :: PointStorage   => FaceStorage_PointStorage
         procedure   :: ConstructAnJac => FaceStorage_ConstructAnJac
#ifdef FLOW
         procedure   :: SetStorageToNS => FaceStorage_SetStorageToNS
#endif
#ifdef CAHNHILLIARD
         procedure   :: SetStorageToCH_c  => FaceStorage_SetStorageToCH_c
         procedure   :: SetStorageToCH_mu => FaceStorage_SetStorageToCH_mu
#endif
         procedure :: copy             => FaceStorage_Assign
         generic   :: assignment(=)    => copy
         final     :: FaceStorage_Finalize
   end type FaceStorage_t
!
!  ========
   contains
!  ========
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           Global Storage procedures
!           --------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------
!     The global solution arrays are only allocated if needed
!     -------------------------------------------------------
      subroutine SolutionStorage_Construct(self, NDOF, Nx, Ny, Nz, computeGradients, analyticalJac, prevSol_num, RKSteps_num)
         implicit none
         !-arguments---------------------------------------------
         class(SolutionStorage_t), target, intent(inout) :: self
         integer                 , intent(in)    :: NDOF
         integer, dimension(:)   , intent(in)    :: Nx, Ny, Nz
         logical                 , intent(in)    :: computeGradients                   !<  Compute gradients?
         logical                 , intent(in)    :: analyticalJac                      !<  Create storage for analytical Jacobian?
         integer, optional       , intent(in)    :: prevSol_num
         integer, optional       , intent(in)    :: RKSteps_num
         !-local-variables---------------------------------------
         integer :: k, eID, num_of_elems
         !-------------------------------------------------------

         self % NDOF = NDOF
         self % anJacobian = analyticalJac

         if ( present(prevSol_num) ) then
            self % prevSol_num = prevSol_num
            if ( prevSol_num > 0 ) then
               allocate ( self % prevSol_index(prevSol_num) )
               self % prevSol_index = (/ (k, k=1, prevSol_num) /)
#ifdef FLOW
               allocate ( self % PrevQNS(NCONS*NDOF, prevSol_num) )
               self % PrevQ    => self % PrevQNS
#endif
#ifdef CAHNHILLIARD
               allocate ( self % Prevc  (NCOMP*NDOF, prevSol_num) )
               self % PrevQ    => self % Prevc
#endif
            end if
            if ( prevSol_num >= 0 ) then
#ifdef FLOW
               allocate ( self % QNS   (NCONS*NDOF) )
               allocate ( self % QdotNS(NCONS*NDOF) )
               self % Q    => self % QNS
               self % Qdot => self % QdotNS
#endif
#ifdef CAHNHILLIARD
               allocate ( self % c     (NCOMP*NDOF) )
               allocate ( self % cDot  (NCOMP*NDOF) )
               self % Q    => self % c
               self % Qdot => self % cDot
#endif
            end if
         end if

         num_of_elems = size(Nx)
         allocate (self % elements(num_of_elems) )

         if ( present(RKSteps_num) .and. present(prevSol_num)) then
!$omp parallel do schedule(runtime)
            do eID=1,  num_of_elems
               call self % elements(eID) % construct( Nx(eID), Ny(eID), Nz(eID), computeGradients, analyticalJac, prevSol_num, RKSteps_num)
            end do
!$omp end parallel do
         elseif ( present(prevSol_num) ) then
!$omp parallel do schedule(runtime)
            do eID=1,  num_of_elems
               call self % elements(eID) % construct( Nx(eID), Ny(eID), Nz(eID), computeGradients, analyticalJac, prevSol_num, 0)
            end do
!$omp end parallel do
         elseif ( present(RKSteps_num) ) then
!$omp parallel do schedule(runtime)
            do eID=1,  num_of_elems
               call self % elements(eID) % construct( Nx(eID), Ny(eID), Nz(eID), computeGradients, analyticalJac, -1, RKSteps_num)
            end do
!$omp end parallel do
         else
!$omp parallel do schedule(runtime)
            do eID=1,  num_of_elems
               call self % elements(eID) % construct( Nx(eID), Ny(eID), Nz(eID), computeGradients, analyticalJac,-1,0)
            end do
!$omp end parallel do


         end if

      end subroutine SolutionStorage_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------
!     Load the local solutions into a global array
!     --------------------------------------------
      pure subroutine SolutionStorage_local2GlobalQ(self, NDOF)
         implicit none
         !----------------------------------------------
         class(SolutionStorage_t), target, intent(inout) :: self
         integer                 , intent(in)    :: NDOF
         !----------------------------------------------
         integer :: firstIdx, lastIdx, eID
         !----------------------------------------------

!        Allocate storage
!        ****************
#ifdef FLOW
         if (self % AdaptedQ .or. (.not. allocated(self % QNS) ) ) then
            self % NDOF = NDOF
            safedeallocate (self % QNS   )
            allocate ( self % QNS   (NCONS*NDOF) )

            self % AdaptedQ = .FALSE.
         end if
#endif
#ifdef CAHNHILLIARD
         if (self % AdaptedQ .or. (.not. allocated(self % c) ) ) then
            self % NDOF = NDOF
            safedeallocate (self % c)
            allocate ( self % c(NCOMP*NDOF) )

            self % AdaptedQ = .FALSE.
         end if
#endif
!
!        Load solution
!        *************

#ifdef FLOW
         firstIdx = 1
         do eID=1, size(self % elements)
            lastIdx = firstIdx + self % elements(eID) % NDOF * NCONS
            self % QNS (firstIdx : lastIdx - 1) = reshape ( self % elements(eID) % QNS , (/ self % elements(eID) % NDOF *NCONS /) )
            firstIdx = lastIdx
         end do
#endif
#ifdef CAHNHILLIARD
         firstIdx = 1
         do eID=1, size(self % elements)
            lastIdx = firstIdx + self % elements(eID) % NDOF * NCOMP
            self % c (firstIdx : lastIdx - 1) = reshape ( self % elements(eID) % c , (/ self % elements(eID) % NDOF *NCOMP /) )
            firstIdx = lastIdx
         end do
#endif

      end subroutine SolutionStorage_local2GlobalQ
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------
!     Load the local Qdot into a global array
!     ---------------------------------------
      pure subroutine SolutionStorage_local2GlobalQdot(self, NDOF)
         implicit none
         !----------------------------------------------
         class(SolutionStorage_t), target, intent(inout) :: self
         integer                 , intent(in)    :: NDOF
         !----------------------------------------------
         integer :: firstIdx, lastIdx, eID
         !----------------------------------------------

         self % NDOF = NDOF

#ifdef FLOW
         if (self % AdaptedQdot .or. (.not. allocated(self % QdotNS) ) ) then
            safedeallocate (self % QdotNS)
            allocate ( self % QdotNS(NCONS*NDOF) )
         end if
#endif
#ifdef CAHNHILLIARD
         if (self % AdaptedQdot .or. (.not. allocated(self % cDot) ) ) then
            safedeallocate ( self % cDot )
            allocate ( self % cDot(NCOMP*NDOF) )
         end if
#endif
         self % AdaptedQdot = .FALSE.

!
!        Load solution
!        *************

#ifdef FLOW
         firstIdx = 1
         do eID=1, size(self % elements)
            lastIdx = firstIdx + self % elements(eID) % NDOF * NCONS
            self % QdotNS (firstIdx : lastIdx - 1) = reshape ( self % elements(eID) % QdotNS , (/ self % elements(eID) % NDOF * NCONS/) )
            firstIdx = lastIdx
         end do
#endif
#ifdef CAHNHILLIARD
         firstIdx = 1
         do eID=1, size(self % elements)
            lastIdx = firstIdx + self % elements(eID) % NDOF * NCOMP
            self % cDot (firstIdx : lastIdx - 1) = reshape ( self % elements(eID) % cDot , (/ self % elements(eID) % NDOF * NCOMP /) )
            firstIdx = lastIdx
         end do
#endif
      end subroutine SolutionStorage_local2GlobalQdot
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------
!     Load the local solutions into a global array
!     --------------------------------------------
      pure subroutine SolutionStorage_SetGlobalPrevQ(self, Q)
         implicit none
         !-arguments---------------------------------------------
         class(SolutionStorage_t), target, intent(inout) :: self
         real(kind=RP)                   , intent(in)    :: Q(:)
         !-local-variables---------------------------------------
         integer :: k, oldest_index
         !-------------------------------------------------------

         if (self % prevSol_num < 1) return

!        Allocate global storage
!        ***********************
#ifdef FLOW
         if (self % AdaptedPrevQ .or. (.not. allocated(self % PrevQNS) ) ) then
            safedeallocate (self % PrevQNS)
            allocate ( self % PrevQNS (NCONS * self % NDOF, self % prevSol_num) )

            self % AdaptedPrevQ = .FALSE.

            ! TODO: Adapt previous solutions...
         end if
#endif
#ifdef CAHNHILLIARD
         if (self % AdaptedPrevQ .or. (.not. allocated(self % PrevC) ) ) then
            safedeallocate (self % PrevC)
            allocate ( self % PrevC(NCOMP * self % NDOF, self % prevSol_num) )

            self % AdaptedPrevQ = .FALSE.

            ! TODO: Adapt previous solutions...
         end if
#endif

!
!        Load solutions
!        **************

         oldest_index = self % prevSol_index ( self % prevSol_num )
         do k=self % prevSol_num, 2, -1
            self % prevSol_index(k) = self % prevSol_index(k-1)
         end do
         self % prevSol_index(1) = oldest_index

         self % PrevQ (:,self % prevSol_index(1)) = Q

      end subroutine SolutionStorage_SetGlobalPrevQ
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine SolutionStorage_global2LocalQ(self)
         implicit none
         !----------------------------------------------
         class(SolutionStorage_t), target, intent(inout)    :: self
         !----------------------------------------------
         integer :: firstIdx, lastIdx, eID
         integer :: nEqn
         !----------------------------------------------

         ! Temporary only checking first element!
         select case (self % elements(1) % currentlyLoaded)
            case (NS,NSSA,CAA)
#ifdef FLOW
               nEqn = NCONS
#endif
            case (C,MU)
#ifdef CAHNHILLIARD
               nEqn = NCOMP
#endif
         end select

         firstIdx = 1
         do eID=1, size(self % elements)
            associate ( N => self % elements(eID) % Nxyz )
            lastIdx = firstIdx + self % elements(eID) % NDOF * nEqn
            self % elements(eID) % Q(1:nEqn,0:N(1),0:N(2),0:N(3)) = reshape ( self % Q (firstIdx:lastIdx-1) , (/ nEqn, N(1)+1, N(2)+1, N(3)+1/) )
            firstIdx = lastIdx
            end associate
         end do

      end subroutine SolutionStorage_global2LocalQ
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine SolutionStorage_global2LocalQdot(self)
         implicit none
         !----------------------------------------------
         class(SolutionStorage_t), target, intent(inout) :: self
         !----------------------------------------------
         integer :: firstIdx, lastIdx, eID
         integer :: nEqn
         !----------------------------------------------

         ! Temporary only checking first element!
         select case (self % elements(1) % currentlyLoaded)
            case (NS,NSSA,CAA)
#ifdef FLOW
               nEqn = NCONS
#endif
            case (C,MU)
#ifdef CAHNHILLIARD
               nEqn = NCOMP
#endif
         end select

         firstIdx = 1
         do eID=1, size(self % elements)
            associate ( N => self % elements(eID) % Nxyz )
            lastIdx = firstIdx + self % elements(eID) % NDOF * nEqn
            self % elements(eID) % Qdot(1:nEqn,0:N(1),0:N(2),0:N(3)) = reshape ( self % Qdot (firstIdx:lastIdx-1) , (/ nEqn, N(1)+1, N(2)+1, N(3)+1/) )
            firstIdx = lastIdx
            end associate
         end do

      end subroutine SolutionStorage_global2LocalQdot
!
!/////////////////////////////////////////////////
!
      pure subroutine SolutionStorage_SignalAdaptation(self)
         implicit none
         class(SolutionStorage_t), intent(inout) :: self

         self % AdaptedQ     = .TRUE.
         self % AdaptedQdot  = .TRUE.
         self % AdaptedPrevQ = .TRUE.
      end subroutine SolutionStorage_SignalAdaptation
!
!/////////////////////////////////////////////////
!
      pure subroutine SolutionStorage_PointStorage(self)
         implicit none
         class(SolutionStorage_t), intent(inout), target :: self
!
!        Point SolutionStorage
!        ---------------------

         select case ( self % elements(1) % currentlyLoaded ) ! Using element 1 (to be deprecated)
            case (OFF)
               self % Q     => NULL()
               self % Qdot  => NULL()
               self % PrevQ => NULL()
#ifdef FLOW
            case (NS,NSSA,CAA)
               self % Q     => self % QNS
               self % Qdot  => self % QdotNS
               self % PrevQ => self % PrevQNS
#endif
#ifdef CAHNHILLIARD
            case (C,MU)
               self % Q     => self % c
               self % Qdot  => self % cDot
               self % PrevQ => self % Prevc
#endif
         end select

!
!        Point elements' Storage
!        -----------------------
         call self % elements % PointStorage()

      end subroutine SolutionStorage_PointStorage
!
!/////////////////////////////////////////////////
!
      pure subroutine SolutionStorage_Destruct(self)
         implicit none
         class(SolutionStorage_t), intent(inout) :: self

         self % prevSol_num = 0
         self % AdaptedQ     = .FALSE.
         self % AdaptedQdot  = .FALSE.
         self % AdaptedPrevQ = .FALSE.
         self % anJacobian   = .FALSE.

         safedeallocate(self % prevSol_index)

#ifdef FLOW
         safedeallocate(self % QNS)
         safedeallocate(self % QdotNS)
         safedeallocate(self % PrevQNS)
#endif
#ifdef CAHNHILLIARD
         safedeallocate(self % c)
         safedeallocate(self % cDot)
         safedeallocate(self % PrevC)
#endif

         if ( allocated(self % elements) ) then
            call self % elements % destruct
            deallocate (self % elements)
         end if

      end subroutine SolutionStorage_Destruct
!
!/////////////////////////////////////////////////
!
!     (We need an special assign procedure)
!
      subroutine SolutionStorage_Assign(to, from)
         implicit none
         !-arguments------------------------------------
         class(SolutionStorage_t), intent(inout) :: to
         type(SolutionStorage_t),  intent(in)    :: from
         !-local-variables------------------------------
         integer :: num_of_elems
         integer :: eID
         !----------------------------------------------
!
!        Copy the storage
!        ----------------
         call to % destruct

         to % NDOF         =  from % NDOF
         to % AdaptedQ     =  from % AdaptedQ
         to % AdaptedQdot  =  from % AdaptedQdot
         to % AdaptedPrevQ =  from % AdaptedPrevQ
         to % prevSol_num  =  from % prevSol_num
         to % anJacobian   =  from % anJacobian



         num_of_elems = size(from % elements)
         allocate ( to % elements (num_of_elems) )
!$omp parallel do schedule(runtime)
         do eID=1, num_of_elems
            to % elements(eID) = from % elements(eID)
         end do
!$omp end parallel do

         if ( to % prevSol_num > 0 ) then
            allocate ( to % prevSol_index ( size(from % prevSol_index) ) )
            to % prevSol_index=  from % prevSol_index
#ifdef FLOW
            allocate ( to % PrevQNS ( size(from % PrevQNS,1),size(from % PrevQNS,2) ) )
            to % PrevQNS      =  from % PrevQNS
#endif
#ifdef CAHNHILLIARD
            allocate ( to % Prevc ( size(from % Prevc,1),size(from % Prevc,2) ) )
            to % Prevc        =  from % Prevc
#endif
         end if

         if ( to % prevSol_num >= 0 ) then
#ifdef FLOW
            allocate ( to % QdotNS ( size(from % QdotNS) ) )
            to % QdotNS       =  from % QdotNS

            allocate ( to % QNS ( size(from % QNS) ) )
            to % QNS          =  from % QNS
#endif
#ifdef CAHNHILLIARD
            allocate ( to % cDot ( size(from % cDot) ) )
            to % cDot         =  from % cDot

            allocate ( to % c ( size(from % c) ) )
            to % c            =  from % c
#endif
         end if

!
!        Point the storage
!        -----------------
         call to % PointStorage()

      end subroutine SolutionStorage_Assign
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!           Element Storage procedures
!           --------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      elemental subroutine ElementStorage_Construct(self, Nx, Ny, Nz, computeGradients, analyticalJac, prevSol_num, RKSteps_num)
         implicit none
         !------------------------------------------------------------
         class(ElementStorage_t), intent(inout) :: self                               !<> Storage to be constructed
         integer                , intent(in)    :: Nx, Ny, Nz                         !<  Polynomial orders in every direction
         logical                , intent(in)    :: computeGradients                   !<  Compute gradients?
         logical                , intent(in)    :: analyticalJac                      !<  Analytical Jacobian specific storage(?)
         integer                , intent(in)    :: prevSol_num
         integer                , intent(in)    :: RKSteps_num
         !------------------------------------------------------------
         integer :: k
         !------------------------------------------------------------

         self % computeGradients = computeGradients
         self % prevSol_num = prevSol_num
         self % RKSteps_num = RKSteps_num

!
!        --------------------------------
!        Get number of degrees of freedom
!        --------------------------------
!
         self % NDOF = (Nx + 1) * (Ny + 1) * (Nz + 1)
         self % Nxyz = [Nx, Ny, Nz]
!
!        ----------------
!        Volume variables
!        ----------------
!
#ifdef FLOW
         allocate ( self % QNS   (1:NCONS,0:Nx,0:Ny,0:Nz) )
         allocate ( self % QLowRK(1:NCONS,0:Nx,0:Ny,0:Nz) )
         allocate (self % QbaseSponge(1:NCONS,0:Nx,0:Ny,0:Nz) )
         allocate (self % intensitySponge(0:Nx,0:Ny,0:Nz) )
         allocate ( self % QdotNS(1:NCONS,0:Nx,0:Ny,0:Nz) )
         allocate ( self % rho   (0:Nx,0:Ny,0:Nz) )
         ! Previous solution
         if ( prevSol_num /= 0 ) then
            allocate ( self % PrevQ(prevSol_num) )
            do k=1, prevSol_num
               allocate ( self % PrevQ(k) % QNS(1:NCONS,0:Nx,0:Ny,0:Nz) )
            end do
         end if

         ALLOCATE( self % G_NS   (NCONS,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % S_NS   (NCONS,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % S_NSP  (NCONS,0:Nx,0:Ny,0:Nz) )
#if defined (SPALARTALMARAS)
         ALLOCATE( self % S_SA  (NCONS,0:Nx,0:Ny,0:Nz) )
#endif
#if defined (ACOUSTIC)
         ALLOCATE( self % Qbase  (NCONSB,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % grada2base  (NDIM,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % Lambbase  (NDIM,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % Lamb  (NDIM,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % Lambread  (NDIM,0:Nx,0:Ny,0:Nz,NLambInterpolation) )
#endif
         if (computeGradients) then
            ALLOCATE( self % U_xNS (NGRAD,0:Nx,0:Ny,0:Nz) )
            ALLOCATE( self % U_yNS (NGRAD,0:Nx,0:Ny,0:Nz) )
            ALLOCATE( self % U_zNS (NGRAD,0:Nx,0:Ny,0:Nz) )
         end if

#ifndef ACOUSTIC
         allocate( self % mu_NS(1:3,0:Nx,0:Ny,0:Nz) )
         allocate( self % mu_turb_NS(0:Nx,0:Ny,0:Nz) )
#endif

         if (analyticalJac) call self % constructAnJac      ! TODO: This is actually not specific for NS

!
!        Point to NS by default
!        ----------------------
         call self % SetStorageToNS
#endif
#ifdef CAHNHILLIARD
         allocate ( self % c   (1:NCOMP,0:Nx,0:Ny,0:Nz) )
         allocate ( self % cDot(1:NCOMP,0:Nx,0:Ny,0:Nz) )
         ! Previous solution
         if ( prevSol_num /= 0 ) then
            if ( .not. allocated(self % PrevQ)) then
               allocate ( self % PrevQ(prevSol_num) )
            end if
         end if

         do k=1, prevSol_num
            allocate ( self % PrevQ(k) % c(1:NCOMP,0:Nx,0:Ny,0:Nz) )
         end do

         allocate(self % c_x (NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % c_y (NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % c_z (NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % mu  (NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % mu_x(NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % mu_y(NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % mu_z(NCOMP, 0:Nx, 0:Ny, 0:Nz))
         ALLOCATE(self % G_CH(NCOMP,0:Nx,0:Ny,0:Nz) )
         allocate(self % v   (1:NDIM, 0:Nx, 0:Ny, 0:Nz))
#endif

#ifdef MULTIPHASE
         if ( RKSteps_num .gt. 0 ) then
            allocate(self % RKSteps(RKSteps_num))

            do k = 1, RKSteps_num
               allocate(self % RKSteps(k) % K(1:NCOMP,0:Nx, 0:Ny, 0:Nz))
               allocate(self % RKSteps(k) % hatK(1:NCONS,0:Nx, 0:Ny, 0:Nz))
            enddo
         end if

         allocate ( self % invMa2   (0:Nx,0:Ny,0:Nz) )
#endif
!
!        -----------------
!        Initialize memory
!        -----------------
!
#ifdef FLOW
         self % G_NS            = 0.0_RP
         self % S_NS            = 0.0_RP
         self % S_NSP           = 0.0_RP
         self % QNS             = 0.0_RP
         self % QLowRK          = 0.0_RP
         self % QbaseSponge     = 0.0_RP
         self % intensitySponge = 0.0_RP
         self % QDotNS          = 0.0_RP
         self % rho             = 0.0_RP
#ifndef ACOUSTIC  
         self % mu_NS  = 0.0_RP
         self % mu_turb_NS  = 0.0_RP
#endif
#if defined (SPALARTALMARAS)
         self % S_SA   = 0.0_RP
#endif
#if defined (ACOUSTIC)
         self % Qbase  = 0.0_RP
         self % grada2base = 0.0_RP
         self % Lambbase = 0.0_RP
         self % Lamb = 0.0_RP
         self % Lambread = 0.0_RP
#endif
         if (computeGradients) then
            self % U_xNS = 0.0_RP
            self % U_yNS = 0.0_RP
            self % U_zNS = 0.0_RP
         end if

         self % artificialDiss = 0.0_RP
#endif

#ifdef CAHNHILLIARD
         self % c     = 0.0_RP
         self % c_x   = 0.0_RP
         self % c_y   = 0.0_RP
         self % c_z   = 0.0_RP
         self % mu    = 0.0_RP
         self % mu_x  = 0.0_RP
         self % mu_y  = 0.0_RP
         self % mu_z  = 0.0_RP
         self % G_CH  = 0.0_RP
         self % v     = 0.0_RP
#endif

#ifdef MULTIPHASE
         self % invMa2     = 0.0_RP
#endif
         self % first_sensed = huge(1)
         self % prev_sensor = 1.0_RP
         self % sensor = 1.0_RP  ! Activate the sensor by default (first time-step when SC is on)

      end subroutine ElementStorage_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ElementStorage_ConstructAnJac:
!     Construct analytical Jacobian (specific) storage
!
      elemental subroutine ElementStorage_ConstructAnJac(self)
         implicit none
         !-arguments----------------------------------------
         class(ElementStorage_t), intent(inout) :: self
         !--------------------------------------------------
#ifdef FLOW
         self % anJacobian = .TRUE.

!        Allocation
!        ----------
         allocate( self % dF_dgradQ( NCONS, NCONS, NDIM, NDIM, 0:self % Nxyz(1), 0:self % Nxyz(2), 0:self % Nxyz(3) ) )

!        Zero storage
!        ------------
         self % dF_dgradQ = 0._RP
#endif
      end subroutine ElementStorage_ConstructAnJac
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      elemental subroutine ElementStorage_Assign(to, from)
!
!        **********************************
!        We need a special assign procedure
!        **********************************
!
         implicit none
         class(ElementStorage_t), intent(inout) :: to
         type(ElementStorage_t),  intent(in)    :: from

         call to % destruct
         call to % construct (from % Nxyz(1), &
                              from % Nxyz(2), &
                              from % Nxyz(3), &
                              from % computeGradients, & ! TODO: Fix this: it is not being used!!
                              from % anJacobian, &
                              from % prevSol_num, &
                              from % RKSteps_num )
!
!        Copy the storage
!        ----------------
         to % currentlyLoaded = from % currentlyLoaded
         to % NDOF            = from % NDOF
         to % Nxyz            = from % Nxyz
!
!        Copy the sensors
!        ----------------
         to % first_sensed = from % first_sensed
         to % prev_sensor  = from % prev_sensor
         to % sensor       = from % sensor

#ifdef FLOW
         to % QNS    = from % QNS
         to % QLowRK = from % QLowRK
         to % QbaseSponge = from % QbaseSponge
         to % intensitySponge = from % intensitySponge

         if (to % computeGradients) then
            to % U_xNS  = from % U_xNS
            to % U_yNS  = from % U_yNS
            to % U_zNS  = from % U_zNS
         end if
         to % QDotNS = from % QDotNS
         to % G_NS   = from % G_NS
         to % S_NS   = from % S_NS
         to % S_NSP  = from % S_NSP

#if defined (SPALARTALMARAS)
         to % S_SA   = from % S_SA
#endif
#if defined (ACOUSTIC)
         to % Qbase   = from % Qbase
         to % grada2base = from % grada2base
         to % Lambbase   = from % Lambbase
         to % Lamb   = from % Lamb
         to % Lambread   = from % Lambread
#endif

#ifndef ACOUSTIC
         to % mu_NS     = from % mu_NS
         to % mu_turb_NS     = from % mu_turb_NS
#endif
         to % stats     = from % stats

         if (to % anJacobian) then
            to % dF_dgradQ = from % dF_dgradQ
         end if
#endif
#ifdef CAHNHILLIARD
         to % c    = from % c
         to % c_x  = from % c_x
         to % c_y  = from % c_y
         to % c_z  = from % c_z
         to % mu   = from % mu
         to % mu_x = from % mu_x
         to % mu_y = from % mu_y
         to % mu_z = from % mu_z
         to % v    = from % v
         to % cDot = from % cDot
         to % G_CH = from % G_CH
#endif

         call to % PointStorage()

      end subroutine ElementStorage_Assign
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      elemental subroutine ElementStorage_PointStorage (self)
         implicit none
         class(ElementStorage_t), intent(inout) :: self

         select case ( self % currentlyLoaded )
            case (OFF)
               self % Q    => NULL()
               self % U_x  => NULL()
               self % U_y  => NULL()
               self % U_z  => NULL()
               self % QDot => NULL()
#ifdef FLOW
            case (NS,NSSA,CAA)
               call self % SetStorageToNS
#endif
#ifdef CAHNHILLIARD
            case (C)
               call self % SetStorageToCH_c

            case (MU)
               call self % SetStorageToCH_mu
#endif
         end select
      end subroutine ElementStorage_PointStorage
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      elemental subroutine ElementStorage_Destruct(self)
         implicit none
         class(ElementStorage_t), intent(inout) :: self
         integer                 :: num_prevSol, k

         self % currentlyLoaded = OFF
         self % NDOF = 0

         self % Q => NULL()
         self % QDot => NULL()
         self % U_x => NULL()
         self % U_y => NULL()
         self % U_z => NULL()

#ifdef FLOW
         safedeallocate(self % QNS)
         safedeallocate(self % QDotNS)
         safedeallocate(self % QbaseSponge)
         safedeallocate(self % intensitySponge)
         safedeallocate(self % QLowRK)

         if ( allocated(self % PrevQ) ) then
            num_prevSol = size(self % PrevQ)
            do k=1, num_prevSol
               safedeallocate( self % PrevQ(k) % QNS )
            end do
         end if

         safedeallocate(self % G_NS)
         safedeallocate(self % S_NS)
         safedeallocate(self % S_NSP)

#if defined (SPALARTALMARAS)
         safedeallocate(self % S_SA)
#endif
#if defined (ACOUSTIC)
         safedeallocate(self % Qbase)
         safedeallocate(self % grada2base)
         safedeallocate(self % Lambbase)
         safedeallocate(self % Lamb)
         safedeallocate(self % Lambread)
#endif

         if (self % computeGradients) then
            safedeallocate(self % U_xNS)
            safedeallocate(self % U_yNS)
            safedeallocate(self % U_zNS)
         end if
#ifndef ACOUSTIC
         safedeallocate(self % mu_NS)
         safedeallocate(self % mu_turb_NS)
#endif
         safedeallocate(self % rho)

         !if (self % anJacobian) then ! Not needed since there's only one variable (= one if)
            safedeallocate(self % dF_dgradQ)
         !end if

         ! Destruct statistics
         call self % stats % destruct()
#endif

#ifdef CAHNHILLIARD
         safedeallocate(self % c)
         safedeallocate(self % cDot)

         if ( allocated(self % PrevQ) ) then
            num_prevSol = size(self % PrevQ)
            do k=1, num_prevSol
               safedeallocate( self % PrevQ(k) % c )
            end do
         end if

         safedeallocate(self % c_x)
         safedeallocate(self % c_y)
         safedeallocate(self % c_z)
         safedeallocate(self % mu)
         safedeallocate(self % mu_x)
         safedeallocate(self % mu_y)
         safedeallocate(self % mu_z)
         safedeallocate(self % G_CH)
         safedeallocate(self % v)
#endif

#ifdef MULTIPHASE
         safedeallocate(self % invMa2)
#endif

         safedeallocate(self % PrevQ)

!        Deallocate RKSteps
!        ------------------
#ifdef MULTIPHASE
         if ( allocated(self % RKSteps) ) then
            do k = 1, size(self % RKSteps)
               safedeallocate(self % RKSteps(k) % K)
               safedeallocate(self % RKSteps(k) % hatK)
            end do
            deallocate(self % RKSteps)
         end if
#endif

      end subroutine ElementStorage_Destruct
#ifdef FLOW
      pure subroutine ElementStorage_SetStorageToNS(self)
!
!        *****************************************
!        This subroutine selects the Navier-Stokes
!        state vector as current storage.
!        *****************************************
!
         implicit none
         class(ElementStorage_t), target, intent(inout) :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: k

#ifndef SPALARTALMARAS
         self % currentlyLoaded = NS
#elif defined(ACOUSTIC)
         self % currentlyLoaded = CAA
#else
         self % currentlyLoaded = NSSA
#endif
         self % Q   (1:,0:,0:,0:) => self % QNS
         if (self % computeGradients) then
            self % U_x (1:,0:,0:,0:) => self % U_xNS
            self % U_y (1:,0:,0:,0:) => self % U_yNS
            self % U_z (1:,0:,0:,0:) => self % U_zNS
         end if
         self % QDot(1:,0:,0:,0:) => self % QDotNS

         do k = 1, size(self % prevQ)
            self % prevQ(k) % Q(1:,0:,0:,0:) => self % prevQ(k) % QNS
         end do

      end subroutine ElementStorage_SetStorageToNS
#endif
#ifdef CAHNHILLIARD
      pure subroutine ElementStorage_SetStorageToCH_c(self)
!
!        *********************************************
!        This subroutine selects the concentration as
!        current storage.
!        *********************************************
!
         implicit none
         class(ElementStorage_t), target, intent(inout) :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: k

         self % currentlyLoaded = C
!
!        Point to the one dimensional pointers with generic arrays
!        ---------------------------------------------------------
         self % Q   (1:,0:,0:,0:) => self % c
         self % U_x (1:,0:,0:,0:) => self % c_x
         self % U_y (1:,0:,0:,0:) => self % c_y
         self % U_z (1:,0:,0:,0:) => self % c_z
         self % QDot(1:,0:,0:,0:) => self % cDot

         do k = 1, size(self % prevQ)
            self % prevQ(k) % Q(1:,0:,0:,0:) => self % prevQ(k) % c
         end do


      end subroutine ElementStorage_SetStorageToCH_c

      pure subroutine ElementStorage_SetStorageToCH_mu(self)
!
!        *************************************************
!        This subroutine selects the chemical potential as
!        current storage, with the particularity that
!        selects also cDot as QDot.
!        *************************************************
!
         implicit none
         class(ElementStorage_t), target, intent(inout) :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: k

         self % currentlyLoaded = MU

         self % Q   (1:,0:,0:,0:) => self % mu
         self % U_x (1:,0:,0:,0:) => self % mu_x
         self % U_y (1:,0:,0:,0:) => self % mu_y
         self % U_z (1:,0:,0:,0:) => self % mu_z
         self % QDot(1:,0:,0:,0:) => self % cDot

         do k = 1, size(self % prevQ)
            self % prevQ(k) % Q(1:,0:,0:,0:) => self % prevQ(k) % c
         end do

      end subroutine ElementStorage_SetStorageToCH_mu
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------
!     Interpolate solution to another element storage
!            this % Q  ->  other % Q
!     -----------------------------------------------
      impure elemental subroutine ElementStorage_InterpolateSolution(this,other,nodes,with_gradients)
         implicit none
         !-arguments----------------------------------------------
         class(ElementStorage_t), intent(inout), target    :: this
         type(ElementStorage_t) , intent(inout), target    :: other
         integer                , intent(in)               :: nodes
         logical, optional      , intent(in)               :: with_gradients
         !-local-variables----------------------------------------
         logical                       :: gradients
         !--------------------------------------------------------
#if (defined(NAVIERSTOKES) || defined(MULTIPHASE))
         if ( present(with_gradients) ) then
            gradients = with_gradients
         else
            gradients = .FALSE.
         end if

         ! Copy the solution if the polynomial orders are the same, if not, interpolate
         if (all(this % Nxyz == other % Nxyz)) then
            other % Q = this % Q
            other % QbaseSponge = this % QbaseSponge
            other % QLowRK = this % QLowRK
         else
!$omp critical
            call NodalStorage(this  % Nxyz(1)) % construct(nodes,this  % Nxyz(1))
            call NodalStorage(this  % Nxyz(2)) % construct(nodes,this  % Nxyz(2))
            call NodalStorage(this  % Nxyz(3)) % construct(nodes,this  % Nxyz(3))
            call NodalStorage(other % Nxyz(1)) % construct(nodes,other % Nxyz(1))
            call NodalStorage(other % Nxyz(2)) % construct(nodes,other % Nxyz(2))
            call NodalStorage(other % Nxyz(3)) % construct(nodes,other % Nxyz(3))
            !------------------------------------------------------------------
            ! Construct the interpolation matrices in every direction if needed
            !------------------------------------------------------------------
            call Tset ( this % Nxyz(1), other % Nxyz(1) ) % construct ( this % Nxyz(1), other % Nxyz(1) )  ! Xi
            call Tset ( this % Nxyz(2), other % Nxyz(2) ) % construct ( this % Nxyz(2), other % Nxyz(2) )  ! Eta
            call Tset ( this % Nxyz(3), other % Nxyz(3) ) % construct ( this % Nxyz(3), other % Nxyz(3) )  ! Zeta
!$omp end critical

            !---------------------------------------------
            ! Interpolate solution to new solution storage
            !---------------------------------------------
            call Interp3DArrays  (Nvars      = NCONS   , &
                                  Nin        = this  % Nxyz , &
                                  inArray    = this  % Q    , &
                                  Nout       = other % Nxyz , &
                                  outArray   = other % Q    )

            if (gradients .and. this % computeGradients .and. other % computeGradients) then
               call Interp3DArrays  (Nvars      = NGRAD  , &
                                     Nin        = this % Nxyz  , &
                                     inArray    = this % U_x   , &
                                     Nout       = other % Nxyz , &
                                     outArray   = other % U_x  )

               call Interp3DArrays  (Nvars      = NGRAD  , &
                                     Nin        = this  % Nxyz , &
                                     inArray    = this  % U_y  , &
                                     Nout       = other % Nxyz , &
                                     outArray   = other % U_y  )

               call Interp3DArrays  (Nvars      = NGRAD  , &
                                     Nin        = this  % Nxyz , &
                                     inArray    = this  % U_z  , &
                                     Nout       = other % Nxyz , &
                                     outArray   = other % U_z  )
            end if

            this % Q(1:,0:,0:,0:) => this % QbaseSponge
            other % Q(1:,0:,0:,0:) => other % QbaseSponge

            call Interp3DArrays  (Nvars      = NCONS   , &
                                  Nin        = this  % Nxyz , &
                                  inArray    = this  % Q    , &
                                  Nout       = other % Nxyz , &
                                  outArray   = other % Q    )

            this % Q(1:,0:,0:,0:) => this % QLowRK
            other % Q(1:,0:,0:,0:) => other % QLowRK

            call Interp3DArrays  (Nvars      = NCONS   , &
                                  Nin        = this  % Nxyz , &
                                  inArray    = this  % Q    , &
                                  Nout       = other % Nxyz , &
                                  outArray   = other % Q    )

            call this % SetStorageToNS()
            call other % SetStorageToNS()

         end if
#endif
      end subroutine ElementStorage_InterpolateSolution
!
!////////////////////////////////////////////////////////////////////////////////////////////
!
!        Face storage procedures
!        -----------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine FaceStorage_Construct(self, NDIM, Nf, Nel, computeGradients, analyticalJac, computeQdot)
         implicit none
         class(FaceStorage_t), intent(inout) :: self
         integer             , intent(in)    :: NDIM
         integer             , intent(in)    :: Nf(2)              ! Face polynomial order
         integer             , intent(in)    :: Nel(2)             ! Element face polynomial order
         logical             , intent(in)    :: computeGradients
         logical             , intent(in)    :: analyticalJac      !<? Construct analytical Jacobian storage?
         logical             , intent(in)    :: computeQdot
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: interfaceFluxMemorySize

         self % Nf  = Nf
         self % Nel = Nel
         self % NDIM = NDIM
         self % computeGradients = computeGradients
         self % computeQdot = computeQdot

         interfaceFluxMemorySize = 0

#ifdef FLOW
         ALLOCATE( self % QNS   (NCONS,0:Nf(1),0:Nf(2)) )

         if (computeGradients) then
            ALLOCATE( self % U_xNS(NGRAD,0:Nf(1),0:Nf(2)) )
            ALLOCATE( self % U_yNS(NGRAD,0:Nf(1),0:Nf(2)) )
            ALLOCATE( self % U_zNS(NGRAD,0:Nf(1),0:Nf(2)) )
         end if
!
         if (computeQdot) then
             ALLOCATE( self % QdotNS   (NCONS,0:Nf(1),0:Nf(2)) )
         end if
#if defined (ACOUSTIC)
         ALLOCATE( self % Qbase (NCONSB,0:Nf(1),0:Nf(2)) )
         ALLOCATE( self % grada2base (NDIM,0:Nf(1),0:Nf(2)) )
         ALLOCATE( self % Lambbase (NDIM,0:Nf(1),0:Nf(2)) )
         ALLOCATE( self % Lamb (NDIM,0:Nf(1),0:Nf(2)) )
         ALLOCATE( self % Lambread (NDIM,0:Nf(1),0:Nf(2),NLambInterpolation) )
#endif
!        Biggest Interface flux memory size is u\vec{n}
!        ----------------------------------------------
         interfaceFluxMemorySize = NGRAD * nDIM * product(Nf + 1)

         allocate( self % rho       (0:Nf(1),0:Nf(2)) )
#ifndef ACOUSTIC
         allocate( self % mu_NS     (1:3,0:Nf(1),0:Nf(2)) )
         allocate( self % u_tau_NS  (0:Nf(1),0:Nf(2)) )
         allocate( self % wallNodeDistance  (0:Nf(1),0:Nf(2)) )
#endif
         
         if (analyticalJac) call self % ConstructAnJac(NDIM) ! This is actually not specific for NS
#endif
#ifdef CAHNHILLIARD
         allocate(self % c   (NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % c_x (NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % c_y (NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % c_z (NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % mu  (NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % mu_x(NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % mu_y(NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % mu_z(NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % v   (1:NDIM, 0:Nf(1), 0:Nf(2)))
!
!        CH will never be the biggest memory requirement unless NSE are disabled
!        -----------------------------------------------------------------------
         interfaceFluxMemorySize = max(interfaceFluxMemorySize, NCOMP*nDIM*product(Nf+1))
#endif

#ifdef MULTIPHASE
         allocate( self % invMa2       (0:Nf(1),0:Nf(2)) )
#endif
!
!        Reserve memory for the interface fluxes
!        ---------------------------------------
         allocate(self % genericInterfaceFluxMemory(interfaceFluxMemorySize))

#ifdef FLOW
!
!        Point to NS by default
!        ----------------------
         call self % SetStorageToNS
#endif
!
!        -----------------
!        Initialize memory
!        -----------------
!
#ifdef FLOW
         self % QNS    = 0.0_RP

         if (computeGradients) then
            self % U_xNS = 0.0_RP
            self % U_yNS = 0.0_RP
            self % U_zNS = 0.0_RP
         end if

         if (computeQdot) then
            self % QdotNS = 0.0_RP
         end if
#if defined (ACOUSTIC)
         self % Qbase    = 0.0_RP
         self % grada2base = 0.0_RP
         self % Lambbase    = 0.0_RP
         self % Lamb    = 0.0_RP
         self % Lambread    = 0.0_RP
#endif

         self % rho    = 0.0_RP
#ifndef ACOUSTIC
         self % mu_NS  = 0.0_RP
         self % u_tau_NS = 0.0_RP
         self % wallNodeDistance = 0.0_RP
#endif

#endif

#ifdef NAVIERSTOKES
         allocate(self % AviscFlux(NCONS,0:Nf(1), 0:Nf(2)))
         self % AviscFlux = 0.0_RP
#endif

#ifdef CAHNHILLIARD
         self % c     = 0.0_RP
         self % c_x   = 0.0_RP
         self % c_y   = 0.0_RP
         self % c_z   = 0.0_RP
         self % mu    = 0.0_RP
         self % mu_x  = 0.0_RP
         self % mu_y  = 0.0_RP
         self % mu_z  = 0.0_RP
         self % v     = 0.0_RP
#endif

#ifdef MULTIPHASE
         self % invMa2    = 0.0_RP
#endif 

         self % constructed = .TRUE.
      end subroutine FaceStorage_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     FaceStorage_ConstructAnJac:
!     Subroutine to construct the analytical Jacobian (specific) storage.
!     -> Must be called during or after face storage construction
!     -> Only for Navier-Stokes (meanwhile)
!
      elemental subroutine FaceStorage_ConstructAnJac(self,NDIM)
         implicit none
         class(FaceStorage_t) , intent(inout) :: self
         integer              , intent(in)    :: NDIM

#ifdef FLOW
         self % anJacobian = .TRUE.

!        Allocate memory
!        ---------------
         allocate( self % dFStar_dqF  (NCONS,NCONS, 0: self % Nf(1), 0: self % Nf(2)) )
         allocate( self % dFStar_dqEl (NCONS,NCONS, 0:self % Nel(1), 0:self % Nel(2),2) )

         if (self % computeGradients) then
            allocate( self % dFv_dGradQF (NCONS,NCONS,NDIM,0: self % Nf(1),0: self % Nf(2)) )
            allocate( self % dFv_dGradQEl(NCONS,NCONS,NDIM,0:self % Nel(1),0:self % Nel(2),2) )
         end if

!        TODO: AMR, if Boundary
         allocate( self % BCJac       (NCONS,NCONS,0:self % Nel(1),0:self % Nel(2)) )

!        Zero memory
!        -----------
         self % dFStar_dqF   = 0.0_RP
         self % dFStar_dqEl  = 0.0_RP
         if (self % computeGradients) then
            self % dFv_dGradQF  = 0.0_RP
            self % dFv_dGradQEl = 0.0_RP
         end if
         self % BCJac        = 0.0_RP
#endif
      end subroutine

      elemental subroutine FaceStorage_Destruct(self)
         implicit none
         class(FaceStorage_t), intent(inout) :: self
         call FaceStorage_Finalize(self)
      end subroutine FaceStorage_Destruct

      elemental subroutine FaceStorage_Finalize(self)
         implicit none
         type(FaceStorage_t), intent(inout) :: self

         self % constructed = .FALSE.
         self % currentlyLoaded = OFF

#ifdef FLOW
         safedeallocate(self % QNS)
         if (self % computeGradients) then
            safedeallocate(self % U_xNS)
            safedeallocate(self % U_yNS)
            safedeallocate(self % U_zNS)
         end if
         if (self % computeQdot) then
             safedeallocate(self % QdotNS)
         end if
         safedeallocate(self % mu_NS)
         safedeallocate(self % u_tau_NS)
         safedeallocate(self % wallNodeDistance)
         safedeallocate(self % rho )
#if defined (ACOUSTIC)
         safedeallocate(self % Qbase )
         safedeallocate(self % grada2base )
         safedeallocate(self % Lambbase )
         safedeallocate(self % Lamb )
         safedeallocate(self % Lambread )
#endif

         self % anJacobian      = .FALSE.

         if (self % anJacobian) then
            safedeallocate(self % dFStar_dqF)
            safedeallocate(self % dFStar_dqEl)
            if (self % computeGradients) then
               safedeallocate(self % dFv_dGradQF)
               safedeallocate(self % dFv_dGradQEl)
            end if
            safedeallocate(self % BCJac )
         end if
#endif

#ifdef CAHNHILLIARD
         safedeallocate(self % c)
         safedeallocate(self % c_x)
         safedeallocate(self % c_y)
         safedeallocate(self % c_z)
         safedeallocate(self % mu)
         safedeallocate(self % mu_x)
         safedeallocate(self % mu_y)
         safedeallocate(self % mu_z)
         safedeallocate(self % v)
#endif

#ifdef MULTIPHASE
         safedeallocate(self % invMa2 )
#endif

         safedeallocate(self % genericInterfaceFluxMemory)

         self % Q      => NULL()
         self % U_x    => NULL() ; self % U_y => NULL() ; self % U_z => NULL()
         self % Qdot   => NULL()
         self % unStar => NULL()
         self % fStar  => NULL()

         safedeallocate(self % AviscFlux)

      end subroutine FaceStorage_Finalize
#ifdef FLOW
      pure subroutine FaceStorage_SetStorageToNS(self)
         implicit none
         class(FaceStorage_t), intent(inout), target    :: self
#ifndef SPALARTALMARAS
         self % currentlyLoaded = NS
#elif defined(ACOUSTIC)
         self % currentlyLoaded = CAA
#else
         self % currentlyLoaded = NSSA
#endif
!        Get sizes
!        ---------
         self % Q   (1:,0:,0:)            => self % QNS
         self % fStar(1:NCONS, 0:self % Nel(1), 0:self % Nel(2)) => self % genericInterfaceFluxMemory

         self % genericInterfaceFluxMemory = 0.0_RP

         if (self % computeGradients) then
            self % U_x (1:,0:,0:) => self % U_xNS
            self % U_y (1:,0:,0:) => self % U_yNS
            self % U_z (1:,0:,0:) => self % U_zNS
            self % unStar(1:NGRAD, 1:NDIM, 0:self % Nel(1), 0:self % Nel(2)) => self % genericInterfaceFluxMemory
         end if

         if (self % computeQdot) then
             self % QDot(1:,0:,0:) => self % QDotNS
         end if

      end subroutine FaceStorage_SetStorageToNS
#endif
#ifdef CAHNHILLIARD
      pure subroutine FaceStorage_SetStorageToCH_c(self)
         implicit none
         class(FaceStorage_t), intent(inout), target  :: self

         self % currentlyLoaded = C
!
!        Get sizes
!        ---------
         self % Q(1:,0:,0:)   => self % c
         self % U_x(1:,0:,0:) => self % c_x
         self % U_y(1:,0:,0:) => self % c_y
         self % U_z(1:,0:,0:) => self % c_z

         self % fStar(1:NCOMP,0:self % Nel(1),0:self % Nel(2))            => self % genericInterfaceFluxMemory
         self % unStar(1:NCOMP, 1:NDIM, 0:self % Nel(1), 0:self % Nel(2)) => self % genericInterfaceFluxMemory

         self % genericInterfaceFluxMemory = 0.0_RP

      end subroutine FaceStorage_SetStorageToCH_c

      pure subroutine FaceStorage_SetStorageToCH_mu(self)
         implicit none
         class(FaceStorage_t), intent(inout), target  :: self

         self % currentlyLoaded = MU

         self % Q(1:,0:,0:)   => self % mu
         self % U_x(1:,0:,0:) => self % mu_x
         self % U_y(1:,0:,0:) => self % mu_y
         self % U_z(1:,0:,0:) => self % mu_z

         self % fStar(1:NCOMP,0:self % Nel(1),0:self % Nel(2))            => self % genericInterfaceFluxMemory
         self % unStar(1:NCOMP, 1:NDIM, 0:self % Nel(1), 0:self % Nel(2)) => self % genericInterfaceFluxMemory

         self % genericInterfaceFluxMemory = 0.0_RP

      end subroutine FaceStorage_SetStorageToCH_mu
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      elemental subroutine FaceStorage_PointStorage (self)
         implicit none
         class(FaceStorage_t), intent(inout)  :: self

         select case ( self % currentlyLoaded )
            case (OFF)
               self % Q      => NULL()
               self % U_x    => NULL()
               self % U_y    => NULL()
               self % U_z    => NULL()
               self % FStar  => NULL()
               self % unStar => NULL()
#ifdef FLOW
            case (NS,NSSA,CAA)
               call self % SetStorageToNS
#endif
#ifdef CAHNHILLIARD
            case (C)
               call self % SetStorageToCH_c

            case (MU)
               call self % SetStorageToCH_mu
#endif
         end select
      end subroutine FaceStorage_PointStorage
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!

      elemental subroutine FaceStorage_Assign(to, from)
         implicit none
         !-arguments--------------------------------------
         class(FaceStorage_t), intent(inout) :: to
         class(FaceStorage_t), intent(in)    :: from
         !------------------------------------------------

         call to % destruct
         if (.not. from % constructed) return
         call to % construct(from % NDIM, from % Nf, from % Nel, from % computeGradients, from % anJacobian, from % computeQdot)

         to % currentlyLoaded = from % currentlyLoaded
         to % Nf = from % Nf
         to % Nel = from % Nel
         to % genericInterfaceFluxMemory = from % genericInterfaceFluxMemory

#ifdef FLOW
         to % QNS = from % QNS
#if defined (ACOUSTIC)
         to % Qbase = from % Qbase
         to % grada2base = from % grada2base
         to % Lambbase = from % Lambbase
         to % Lamb = from % Lamb
         to % Lambread = from % Lambread
#endif
         if (to % computeGradients) then
            to % U_xNS = from % U_xNS
            to % U_yNS = from % U_yNS
            to % U_zNS = from % U_zNS
         end if
         if (to % computeQdot) to % QdotNS = from % QdotNS
         to % rho = from % rho
#ifndef ACOUSTIC
         to % mu_NS  = from % mu_NS
         to % u_tau_NS  = from % u_tau_NS
         to % wallNodeDistance  = from % wallNodeDistance
#endif

         if (to % anJacobian) then
            to % dFStar_dqF = from % dFStar_dqF
            to % dFStar_dqEl = from % dFStar_dqEl
            if (to % computeGradients) then
               to % dFv_dGradQF = from % dFv_dGradQF
               to % dFv_dGradQEl = from % dFv_dGradQEl
            end if
            to % BCJac = from % BCJac
         end if
#endif
#ifdef CAHNHILLIARD
         to % c = from % c
         to % c_x = from % c_x
         to % c_y = from % c_y
         to % c_z = from % c_z
         to % mu = from % mu
         to % mu_x = from % mu_x
         to % mu_y = from % mu_y
         to % mu_z = from % mu_z
         to % v = from % v
#endif
#ifdef MULTIPHASE
         to % invMa2 = from % invMa2
#endif
         call to % PointStorage
      end subroutine FaceStorage_Assign
!
!/////////////////////////////////////////////////////////////////////////////////////
!
!           Statistics procedures
!           ---------------------
!
!/////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Statistics_Construct(self, no_of_variables, N)
         implicit none
         class(Statistics_t)           :: self
         integer,          intent(in)  :: no_of_variables
         integer,          intent(in)  :: N(3)
!
!        Allocate and initialize
!        -----------------------
         allocate( self % data(no_of_variables, 0:N(1), 0:N(2), 0:N(3) ) )
         self % data = 0.0_RP

      end subroutine Statistics_Construct

      pure subroutine Statistics_Destruct(self)
         implicit none
         class(Statistics_t), intent(inout) :: self

         safedeallocate( self % data )

      end subroutine Statistics_Destruct

      subroutine GetStorageEquations(off_, ns_, c_, mu_, nssa_, caa_)
         implicit none
         integer, intent(out) :: off_, ns_, c_, mu_, nssa_, caa_

         off_ = OFF
         ns_  = NS
         c_   = C
         mu_  = MU
         nssa_= NSSA
         caa_ = CAA

      end subroutine GetStorageEquations

end module StorageClass
