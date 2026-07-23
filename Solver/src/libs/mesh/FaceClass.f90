!!
!////////////////////////////////////////////////////////////////////////
!
!      Modification history:
!           Modified to 3D             5/27/15, 11:13 AM: David A. Kopriva
!           Added isotropic mortars    4/26/17, 11:12 AM: Andrés Rueda
!           Added anisotropic mortars  5/16/17, 11:11 AM: Andrés Rueda
!           Embedded mortars          11/13/17, 05:37 PM: Juan Manzanero
!
!////////////////////////////////////////////////////////////////////////
!
      Module FaceClass
         use SMConstants
         use MeshTypes
         use PolynomialInterpAndDerivsModule
         use GaussQuadrature
         use MappedGeometryClass
         use StorageClass                    , only: FaceStorage_t
         use PhysicsStorage
         use NodalStorageClass
         use InterpolationMatrices           , only: Tset, TsetM
         IMPLICIT NONE 
   
         private
         public   Face
!
!     ************************************************************************************
!
!           Face derived type. The face connects two elements,
!           =================  specified in "elementIDs". 
!
!        The face handles two problems:
!
!           1) Rotation: faces are not oriented in the same direction.
!              --------  The direction adopted is always that of the left
!                        element.
!           2) Different polynomial order: elements at each sidecan have different
!              --------------------------  polynomial orders.
!
!        We have defined the following quantities:
!
!           -> NelLeft/Right: The face original order in the element, without
!              ~~~~~~~~~~~~~  rotation or projections.
!
!           -> NfLeft/Right: The face original order rotated to match each other,
!              ~~~~~~~~~~~~  without projection.
!                       ** Thus, always NfLeft = NelLeft, and
!                          NfRight = rotation(NelRight)                           
!
!           -> Nf: The face polynomial order. This is the maximum along 
!              ~~  each side of both directions.
!
!        The process from the element storage to face storage is:
!
!          (0:Nxyz)           (0:NelLeft)                (0:Nf)
!        Left element -- Interpolation to Boundary -- Projection to face
!
!          (0:Nxyz)           (0:NelRight)             (0:NfRight)       (0:Nf)
!        Right element -- Interpolation to Boundary --  Rotation -- Projection to face
!
!     ************************************************************************************
!



      type Face
      integer, allocatable            :: Mortar(:)                !fID of the slave mortar 
      integer                         :: MortarType               !0 = h-conforming or p-mortar, 1 = big hp-master mortar, 2 = small hp-slave, 3 = sliding h-mortar (the face is connected to a mesh % mortar_faces)
      integer                         :: Mortarpos                !Mortar index (only for slave faces, from 1 to 2 or 4; 0 if not slave mortar face)
      logical                         :: flat
      integer                         :: ID                       ! face ID
      integer                         :: FaceType                 ! Type of face: 1 = HMESH_INTERIOR, 2 = HMESH_BOUNDARY, 3 = 2 = HMESH_MPI
      integer                         :: zone                     ! In the case of HMESH_BOUNDARY, which zone it belongs to
      integer                         :: rotation                 ! Relative orientation between faces
      integer                         :: NelLeft(2)               ! Left element face polynomial order
      integer                         :: NelRight(2)              ! Right element face polynomial order
      integer                         :: NfLeft(2)                ! Left face polynomial order
      integer                         :: NfRight(2)               ! Right face polynomial order
      integer                         :: Nf(2)                    ! Face polynomial order
      integer                         :: nodeIDs(4)               ! Face nodes ID
      integer                         :: elementIDs(2)            ! Convention is: 1 = Left, 2 = Right
      integer                         :: elementSide(2)     
      integer                         :: projectionType(2)  
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryName
      type(MappedGeometryFace)        :: geom
      type(FaceStorage_t)             :: storage(2)
      integer                         :: n_mpi_mortar
      real(kind=RP)                   :: offset(2)                   
      real(kind=RP)                   :: s(2)

         contains
            procedure   :: Construct                     => ConstructFace
            procedure   :: Destruct                      => DestructFace
            procedure   :: Print                         => PrintFace
            procedure   :: LinkWithElements              => Face_LinkWithElements
            procedure   :: AdaptSolutionToFace           => Face_AdaptSolutionToFace
            procedure   :: AdaptSolutionToMortarFace     => Face_AdaptSolutionToMortarFace
            procedure   :: AdaptGradientsToFace          => Face_AdaptGradientsToFace
            procedure   :: AdaptGradientsToMortarFace    => Face_AdaptGradientsToMortarFace
            procedure   :: AdaptAviscFluxToFace          => Face_AdaptAviscFluxToFace
            procedure   :: AdaptAviscFluxToMortarFace    => Face_AdaptAviscFluxToMortarFace

            procedure   :: ProjectFluxToElements                => Face_ProjectFluxToElements
            procedure   :: ProjectMortarFluxToElements          => Face_ProjectMortarFluxToElements
            procedure   :: ProjectGradientFluxToElements        => Face_ProjectGradientFluxToElements
            procedure   :: ProjectMortarGradientFluxToElements  => Face_ProjectMortarGradientFluxToElements
            procedure   :: Interpolatebig2small                => Face_Interpolatebig2small
            procedure   :: Interpolatesmall2big                => Face_Interpolatesmall2big
            procedure   :: Interpolatesmall2biggrad            => Face_Interpolatesmall2biggrad
#if defined(NAVIERSTOKES)
            procedure   :: ProjectFluxJacobianToElements       => Face_ProjectFluxJacobianToElements
            procedure   :: ProjectMortarFluxJacobianToElements => Face_ProjectMortarFluxJacobianToElements
            procedure   :: ProjectGradJacobianToElements       => Face_ProjectGradJacobianToElements
            procedure   :: ProjectMortarGradJacobianToElements => Face_ProjectMortarGradJacobianToElements
            procedure   :: ProjectBCJacobianToElements         => Face_ProjectBCJacobianToElements
#endif
#if defined(ACOUSTIC)
            procedure   :: AdaptBaseSolutionToFace       => Face_AdaptBaseSolutionToFace
            procedure   :: AdaptBaseSolutionToMortarFace => Face_AdaptBaseSolutionToMortarFace
            procedure   :: Interpolatebig2smallacoustic  => Face_Interpolatebig2smallacoustic
#endif
            procedure   :: copy           => Face_Assign
            generic     :: assignment(=)  => copy
      end type Face
!
!     ========
      CONTAINS
!     ========
!
      SUBROUTINE ConstructFace( self, ID, nodeIDs, elementID, side )
!
!        *******************************************************
!
!           This is not the full face construction, but a 
!         variable initialization. Just nodes are introduced and
!         the left element.
!
!        *******************************************************
!
         IMPLICIT NONE 
         class(Face) :: self
         integer     :: ID, nodeIDs(4), elementID, side
!
!        --------------------------
!        Set nodes and left element 
!        --------------------------
!
         self % ID             = ID
         self % nodeIDS        = nodeIDs
         self % elementIDs(1)  = elementID
         self % elementIDs(2)  = -1
         self % elementSide(1) = side
!
!        ------------
!        Set defaults
!        ------------
!
         self % FaceType       = HMESH_UNDEFINED
         self % elementIDs(2)  = HMESH_NONE
         self % elementSide(2) = HMESH_NONE
         self % boundaryName   = emptyBCName
         self % rotation       = 0
         self % zone           = 0
         self % MortarType     = MORTAR_NONE
         self % n_mpi_mortar   = 0
      end SUBROUTINE ConstructFace
!
!////////////////////////////////////////////////////////////////////////
!
      elemental SUBROUTINE DestructFace( self )
         IMPLICIT NONE 
         class(Face), intent(inout) :: self
         
         self % ID = -1
         self % FaceType = HMESH_NONE
         self % rotation = 0
         self % NelLeft = -1
         self % NelRight = -1       
         self % NfLeft = -1       
         self % NfRight = -1       
         self % Nf = -1         
         self % nodeIDs = -1             
         self % elementIDs = -1
         self % elementSide = -1
         self % projectionType = -1
         self % boundaryName = ""
         
         call self % geom % Destruct
         call self % storage % Destruct

         !safedeallocate(self % Mortar)

      end SUBROUTINE DestructFace
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintFace( self ) 
      IMPLICIT NONE
      class(Face) :: self
      PRINT *, "Face type = "   , self % FaceType
      PRINT *, "Element IDs: "  , self % elementIDs
      PRINT *, "Element Sides: ", self % elementSide
      IF ( self % FaceType == HMESH_INTERIOR )     THEN
         PRINT *, "Neighbor rotation: ", self  %  rotation
      ELSE
         PRINT *, "Boundary name = ", self % boundaryName
      end IF
      PRINT *, "-----------------------------------"
      end SUBROUTINE PrintFace
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Face_LinkWithElements( self, NelLeft, NelRight, nodeType, offset, s)
         IMPLICIT NONE
         class(Face)        ,     intent(INOUT) :: self        ! Current face
         integer,                 intent(in)    :: NelLeft(2)  ! Left element face polynomial order
         integer,                 intent(in)    :: NelRight(2) ! Right element face polynomial order
         integer,                 intent(in)    :: nodeType    ! Either Gauss or Gauss-Lobatto
         real(kind=RP), optional, intent(in)    :: offset(2)   ! If present, compute hp-mortar projection for 4:1
         real(kind=RP), optional, intent(in)    :: s(2)        ! If present, compute h-mortar projection for sliding elements

#if (!defined(NAVIERSTOKES)) && (!defined(INCNS))
         logical  :: computeGradients = .true.
#endif

!     
!     -------------------------------------------------------------
!     First, get face elements polynomial orders (without rotation)
!     -------------------------------------------------------------
!
      self % NelLeft = NelLeft
      self % NelRight = NelRight
!
!     ---------------------------------------------------------
!     Second, get face polynomial orders (considering rotation)
!     ---------------------------------------------------------
!
      self % NfLeft = self % NelLeft     ! Left elements are always oriented.
      
      SELECT CASE ( self % rotation )
      CASE ( 0, 2, 5, 7 ) ! Local x and y axis are parallel or antiparallel
         self % NfRight = self % NelRight
      CASE ( 1, 3, 4, 6 ) ! Local x and y axis are perpendicular
         self % NfRight(2)  = self % NelRight(1)
         self % NfRight(1)  = self % NelRight(2)
      end SELECT
!
!     ------------------------------------------------------------------
!     Third, the face polynomial order (the maximum for both directions)
!     ------------------------------------------------------------------
!
      self % Nf(1) = max(self % NfLeft(1), self % NfRight(1))
      self % Nf(2) = max(self % NfLeft(2), self % NfRight(2))
!
!     ------------------------------
!     Construct needed nodal storage
!     ------------------------------
!
      call NodalStorage(self % Nf(1)) % Construct(nodeType, self % Nf(1))
      call NodalStorage(self % NfLeft(1)) % Construct(nodeType, self % NfLeft(1))
      call NodalStorage(self % NfRight(1)) % Construct(nodeType, self % NfRight(1))

      call NodalStorage(self % Nf(2)) % Construct(nodeType, self % Nf(2))
      call NodalStorage(self % NfLeft(2)) % Construct(nodeType, self % NfLeft(2))
      call NodalStorage(self % NfRight(2)) % Construct(nodeType, self % NfRight(2))
!
!     -----------------------------------------------------------------------
!     Construction of the projection matrices (simple Lagrange interpolation)
!     -----------------------------------------------------------------------
!
      if (.not.present(offset)) then 
         call Tset(self % NfLeft(1), self % Nf(1)) % construct(self % NfLeft(1), self % Nf(1))
         call Tset(self % Nf(1), self % NfLeft(1)) % construct(self % Nf(1), self % NfLeft(1))
         
         call Tset(self % NfLeft(2), self % Nf(2)) % construct(self % NfLeft(2), self % Nf(2))
         call Tset(self % Nf(2), self % NfLeft(2)) % construct(self % Nf(2), self % NfLeft(2))

         call Tset(self % NfRight(1), self % Nf(1)) % construct(self % NfRight(1), self % Nf(1))
         call Tset(self % Nf(1), self % NfRight(1)) % construct(self % Nf(1), self % NfRight(1))
         
         call Tset(self % NfRight(2), self % Nf(2)) % construct(self % NfRight(2), self % Nf(2))
         call Tset(self % Nf(2), self % NfRight(2)) % construct(self % Nf(2), self % NfRight(2))


      end if 
      if (present(offset) .and. (.not.present(s))) then !4:1
         call TsetM(self % NfLeft(1), self % Nf(1), 2, 1) % construct(self % NfLeft(1), self % Nf(1), 0.5_RP, 0.5_RP, big2small)  
         call TsetM(self % Nf(1), self % NfLeft(1), 2, 2) % construct(self % Nf(1), self % NfLeft(1), 0.5_RP, 0.5_RP, small2big)  

         call TsetM(self % NfLeft(1), self % Nf(1), 1, 1) % construct(self % NfLeft(1), self % Nf(1), -0.5_RP, 0.5_RP, big2small) 
         call TsetM(self % Nf(1), self % NfLeft(1), 1, 2) % construct(self % Nf(1), self % NfLeft(1), -0.5_RP, 0.5_RP, small2big) 
         
         call TsetM(self % NfLeft(2), self % Nf(2), 2, 1) % construct(self % NfLeft(2), self % Nf(2), 0.5_RP, 0.5_RP, big2small)  
         call TsetM(self % Nf(2), self % NfLeft(2), 2, 2) % construct(self % Nf(2), self % NfLeft(2), 0.5_RP, 0.5_RP, small2big)  

         call TsetM(self % NfLeft(2), self % Nf(2), 1, 1) % construct(self % NfLeft(2), self % Nf(2), -0.5_RP, 0.5_RP, big2small)  
         call TsetM(self % Nf(2), self % NfLeft(2), 1, 2) % construct(self % Nf(2), self % NfLeft(2), -0.5_RP, 0.5_RP, small2big)  
      end if 
      if (present(offset) .and. present(s)) then !Sliding
         if (self% Mortarpos==0) then
            call TsetM(self % NfLeft(1), self % Nf(1), 1, 1) % construct(self % NfLeft(1), self % Nf(1), offset(2), s(1), big2small)
            call TsetM(self % Nf(1), self % NfLeft(1), 1, 2) % construct(self % Nf(1), self % NfLeft(1), offset(2), s(1), small2big)

            call TsetM(self % NfLeft(1), self % Nf(1), 2, 1) % construct(self % NfLeft(1), self % Nf(1), offset(1), s(2), big2small)    
            call TsetM(self % Nf(1), self % NfLeft(1), 2, 2) % construct(self % Nf(1), self % NfLeft(1), offset(1), s(2), small2big)
         else if (self%Mortarpos==1)  then
            call TsetM(self % NfLeft(1), self % Nf(1), 3, 1) % construct(self % NfLeft(1), self % Nf(1), offset(2), s(1), big2small)       
            call TsetM(self % Nf(1), self % NfLeft(1), 3, 2) % construct(self % Nf(1), self % NfLeft(1), offset(2), s(1), small2big)

            call TsetM(self % NfLeft(1), self % Nf(1), 4, 1) % construct(self % NfLeft(1), self % Nf(1), offset(1), s(2), big2small)    
            call TsetM(self % Nf(1), self % NfLeft(1), 4, 2) % construct(self % Nf(1), self % NfLeft(1), offset(1), s(2), small2big)
         end if 
      end if 
      

   
!
!     -----------------------  0- no projection
!     Set the projection type: 1- x needs projection
!     -----------------------  2- y needs projection
!                              3- both x and y need projection
      if ( (self % NfLeft(1) .eq. self % Nf(1)) .and. (self % NfLeft(2) .eq. self % Nf(2)) ) then
         self % projectionType(1) = 0
      elseif ( (self % NfLeft(1) .ne. self % Nf(1)) .and. (self % NfLeft(2) .eq. self % Nf(2)) ) then
         self % projectionType(1) = 1
      elseif ( (self % NfLeft(1) .eq. self % Nf(1)) .and. (self % NfLeft(2) .ne. self % Nf(2)) ) then
         self % projectionType(1) = 2
      elseif ( (self % NfLeft(1) .ne. self % Nf(1)) .and. (self % NfLeft(2) .ne. self % Nf(2)) ) then
         self % projectionType(1) = 3
      end if

      if ( (self % NfRight(1) .eq. self % Nf(1)) .and. (self % NfRight(2) .eq. self % Nf(2)) ) then
         self % projectionType(2) = 0
      elseif ( (self % NfRight(1) .ne. self % Nf(1)) .and. (self % NfRight(2) .eq. self % Nf(2)) ) then
         self % projectionType(2) = 1
      elseif ( (self % NfRight(1) .eq. self % Nf(1)) .and. (self % NfRight(2) .ne. self % Nf(2)) ) then
         self % projectionType(2) = 2
      elseif ( (self % NfRight(1) .ne. self % Nf(1)) .and. (self % NfRight(2) .ne. self % Nf(2)) ) then
         self % projectionType(2) = 3
      end if

   end SUBROUTINE Face_LinkWithElements
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Face_AdaptSolutionToFace(self, nEqn, Nelx, Nely, Qe, side, QdotE, computeQdot)
      use MappedGeometryClass
      implicit none
      class(Face),   intent(inout)              :: self
      integer,       intent(in)                 :: nEqn
      integer,       intent(in)                 :: Nelx, Nely
      real(kind=RP), intent(in)                 :: Qe(1:nEqn, 0:Nelx, 0:Nely)
      integer,       intent(in)                 :: side
      real(kind=RP), intent(in), optional       :: QdotE(1:nEqn, 0:Nelx, 0:Nely)
      logical,       intent(in), optional       :: computeQdot
!
!     ---------------
!     Local variables
!     ---------------
!
      integer       :: i, j, k, l, m, ii, jj, a, b, lm, inb, jnb , p, q
      real(kind=RP) :: Qe_rot(1:nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
      real(kind=RP) :: QdotE_rot(1:nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
      logical :: prolongQdot

      ! prolongQdot = present(QdotE)
      if (present(computeQdot)) then
          prolongQdot = computeQdot
      else
          prolongQdot = .FALSE.
      end if
      
      select case (side)
      case(1)

         if (self % MortarType == MORTAR_SMALL4) then 
            error stop 'MortarType SMALL4 reached in subroutine AdaptSolutionToFace expecting non-SMALL4: check calling logic'
         end if

         associate(Qf => self % storage(1) % Q)
         select case ( self % projectionType(1) )
         case (0)
            Qf = Qe
            if (prolongQdot) self % storage(1) % Qdot = QdotE

         case (1)
            Qf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfLeft(1)   ; do i = 0, self % Nf(1)
               Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * Qe(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            Qf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * Qe(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            Qf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfLeft(1) ; do i = 0, self % Nf(1)
                  Qf(:,i,j) = Qf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                            * Qe(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select
         end associate
      case(2) 

         associate( Qf => self % storage(2) % Q )
         do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
            call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
            Qe_rot(:,i,j) = Qe(:,ii,jj) 
         end do                        ; end do

         select case ( self % projectionType(2) )
         case (0)
            Qf = Qe_rot
         case (1)
            Qf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfRight(1)   ; do i = 0, self % Nf(1)
               Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * Qe_rot(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            Qf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * Qe_rot(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            Qf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfRight(1) ; do i = 0, self % Nf(1)
                  Qf(:,i,j) = Qf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                            * Qe_rot(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select
         end associate
      end select


   end subroutine Face_AdaptSolutionToFace
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine Face_AdaptSolutionToMortarFace(self, nEqn, Nelx, Nely, Qe, side, QdotE, computeQdot, fma, sliding)
      use MappedGeometryClass
      implicit none
      class(Face),   intent(inout)              :: self
      integer,       intent(in)                 :: nEqn
      integer,       intent(in)                 :: Nelx, Nely
      real(kind=RP), intent(in)                 :: Qe(1:nEqn, 0:Nelx, 0:Nely)
      integer,       intent(in)                 :: side
      real(kind=RP), intent(in), optional       :: QdotE(1:nEqn, 0:Nelx, 0:Nely)
      logical,       intent(in), optional       :: computeQdot
      type(Face), intent(inout)                 :: fma
      logical,       intent(in), optional       :: sliding

      integer       :: i, j, l, m, ii, jj

      real(kind=RP) :: MIntXi(0:fma%Nf(1), 0:fma%NfLeft(1))
      real(kind=RP) :: MIntEta(0:fma%Nf(2), 0:fma%NfLeft(2))
      real(kind=RP) :: MIntSliding(0:fma%Nf(1), 0:fma%NfLeft(1), 1:2)

      real(kind=RP) :: Qe_rot(1:nEqn, 0:Nelx, 0:Nely)
      real(kind=RP) :: tmp(1:nEqn, 0:fma%Nf(1), 0:fma%NfLeft(2))

      Qe_rot = 0.0_RP

      if ((fma % MortarType .ne. MORTAR_SMALL4) .and. .not. present(sliding)) then
         error stop 'MortarType SMALL4 reached in subroutine AdaptSolutionToMortarFace expecting non-SMALL4: check calling logic'
      end if

      if (.not. present(sliding)) then
         call GetMortarMInt(fma, MIntXi, MIntEta)
      else
         call GetSlidingMInt(fma, MIntSliding)
      end if

      if (.not. present(sliding)) then
         associate(Qf => fma % storage(1) % Q)
            Qf = 0.0_RP

            ! Pass 1: contraction direction 1 (xi)  -> tmp(i,l)
            tmp = 0.0_RP
            do l = 0, fma % NfLeft(2) ; do i = 0, fma % Nf(1) ; do m = 0, fma % NfLeft(1)
               tmp(:,i,l) = tmp(:,i,l) + MIntXi(i,m) * Qe(:,m,l)
            end do ; end do ; end do

            ! Pass 2: contraction direction 2 (eta) -> Qf(i,j)
            do j = 0, fma % Nf(2) ; do i = 0, fma % Nf(1) ; do l = 0, fma % NfLeft(2)
               Qf(:,i,j) = Qf(:,i,j) + MIntEta(j,l) * tmp(:,i,l)
            end do ; end do ; end do

         end associate

      else !sliding mortars
         if (side==1) then
            associate(Qf => fma % storage(1) % Q)
               Qf=0.0_RP
               do l = 0, self % NfLeft(1)  ; do j = 0, self % Nf(1)   ; do i = 0, self % Nf(1)
               Qf(:,i,j) = Qf(:,i,j)  + MIntSliding(j,l,1) * Qe(:,i,l)
               end do                  ; end do                   ; end do
            end associate
         else
            do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
               call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), fma % rotation, ii, jj)
               Qe_rot(:,i,j) = Qe(:,ii,jj)
            end do                        ; end do
            associate(Qf => fma % storage(2) % Q)
            Qf=0.0_RP
               do l = 0, self % NfLeft(1)  ; do j = 0, self % Nf(1)   ; do i = 0, self % Nf(1)
                  Qf(:,i,j) = Qf(:,i,j)  + MIntSliding(j,l,2) *  Qe_rot(:,i,l)
               end do                  ; end do                   ; end do
            end associate
         end if
      end if
   end subroutine Face_AdaptSolutionToMortarFace
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Face_AdaptGradientsToFace(self, nEqn, Nelx, Nely, Uxe, Uye, Uze, side)
      use MappedGeometryClass
      implicit none
      class(Face),   intent(inout)  :: self
      integer,       intent(in)     :: nEqn 
      integer,       intent(in)     :: Nelx, Nely
      real(kind=RP), intent(in)     :: Uxe(nEqn, 0:Nelx, 0:Nely)
      real(kind=RP), intent(in)     :: Uye(nEqn, 0:Nelx, 0:Nely)
      real(kind=RP), intent(in)     :: Uze(nEqn, 0:Nelx, 0:Nely)
      integer,       intent(in)     :: side

!
!     ---------------
!     Local variables
!     ---------------
!
      integer       :: i, j, k, l, m, ii, jj, a, b, lm
      real(kind=RP) :: Uxe_rot(nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
      real(kind=RP) :: Uye_rot(nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
      real(kind=RP) :: Uze_rot(nEqn, 0:self % NfRight(1), 0:self % NfRight(2))

      
      select case (side)
      case(1)
         if (self % MortarType == MORTAR_SMALL4) then 
            error stop 'MortarType SMALL4 reached in subroutine AdaptGradientsToFace expecting non-SMALL4: check calling logic'
         end if

         associate(Uxf => self % storage(1) % U_x, &
                     Uyf => self % storage(1) % U_y, &
                     Uzf => self % storage(1) % U_z   )
         select case ( self % projectionType(1) )
         case (0)
            Uxf = Uxe
            Uyf = Uye
            Uzf = Uze
         case (1)
            Uxf = 0.0_RP
            Uyf = 0.0_RP
            Uzf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfLeft(1)   ; do i = 0, self % Nf(1)
               Uxf(:,i,j) = Uxf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * Uxe(:,l,j)
               Uyf(:,i,j) = Uyf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * Uye(:,l,j)
               Uzf(:,i,j) = Uzf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * Uze(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            Uxf = 0.0_RP
            Uyf = 0.0_RP
            Uzf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               Uxf(:,i,j) = Uxf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * Uxe(:,i,l)
               Uyf(:,i,j) = Uyf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * Uye(:,i,l)
               Uzf(:,i,j) = Uzf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * Uze(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            Uxf = 0.0_RP
            Uyf = 0.0_RP
            Uzf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfLeft(1) ; do i = 0, self % Nf(1)
                  Uxf(:,i,j) = Uxf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                             * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                             * Uxe(:,m,l)
                  Uyf(:,i,j) = Uyf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                             * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                             * Uye(:,m,l)
                  Uzf(:,i,j) = Uzf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                             * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                             * Uze(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select
         end associate
      case(2)
         associate(Uxf => self % storage(2) % U_x, &
                     Uyf => self % storage(2) % U_y, &
                     Uzf => self % storage(2) % U_z   )
         do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
            call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
            Uxe_rot(:,i,j) = Uxe(:,ii,jj) 
            Uye_rot(:,i,j) = Uye(:,ii,jj) 
            Uze_rot(:,i,j) = Uze(:,ii,jj) 
         end do                        ; end do

         select case ( self % projectionType(2) )
         case (0)
            Uxf = Uxe_rot
            Uyf = Uye_rot
            Uzf = Uze_rot
         case (1)
            Uxf = 0.0_RP
            Uyf = 0.0_RP
            Uzf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfRight(1)   ; do i = 0, self % Nf(1)
               Uxf(:,i,j) = Uxf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * Uxe_rot(:,l,j)
               Uyf(:,i,j) = Uyf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * Uye_rot(:,l,j)
               Uzf(:,i,j) = Uzf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * Uze_rot(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            Uxf = 0.0_RP
            Uyf = 0.0_RP
            Uzf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               Uxf(:,i,j) = Uxf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * Uxe_rot(:,i,l)
               Uyf(:,i,j) = Uyf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * Uye_rot(:,i,l)
               Uzf(:,i,j) = Uzf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * Uze_rot(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            Uxf = 0.0_RP
            Uyf = 0.0_RP
            Uzf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfRight(1) ; do i = 0, self % Nf(1)
                  Uxf(:,i,j) = Uxf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                             * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                             * Uxe_rot(:,m,l)
                  Uyf(:,i,j) = Uyf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                             * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                             * Uye_rot(:,m,l)
                  Uzf(:,i,j) = Uzf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                             * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                             * Uze_rot(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select
         end associate
      end select
   
   end subroutine Face_AdaptGradientsToFace
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine Face_AdaptGradientsToMortarFace(self, nEqn, Nelx, Nely, Uxe, Uye, Uze, side, fma, sliding)
      use MappedGeometryClass
      implicit none
      class(Face),   intent(inout)  :: self
      integer,       intent(in)     :: nEqn 
      integer,       intent(in)     :: Nelx, Nely
      real(kind=RP), intent(in)     :: Uxe(nEqn, 0:Nelx, 0:Nely)
      real(kind=RP), intent(in)     :: Uye(nEqn, 0:Nelx, 0:Nely)
      real(kind=RP), intent(in)     :: Uze(nEqn, 0:Nelx, 0:Nely)
      integer,       intent(in)     :: side
      type(Face), intent(inout)      ::fma
      logical,       intent(in), optional       :: sliding
!
!     ---------------
!     Local variables
!     ---------------
!
      integer       :: i, j, l, m, ii, jj 
      real(kind=RP) :: MIntXi(0:fma%Nf(1), 0:fma%NfLeft(1))
      real(kind=RP) :: MIntEta(0:fma%Nf(2), 0:fma%NfLeft(2))
      real(kind=RP) :: MIntSliding(0:fma%Nf(1), 0:fma%NfLeft(1), 1:2)

      real(kind=RP)     :: Uxe_rot(nEqn, 0:Nelx, 0:Nely)
      real(kind=RP)     :: Uye_rot(nEqn, 0:Nelx, 0:Nely)
      real(kind=RP)     :: Uze_rot(nEqn, 0:Nelx, 0:Nely)

      real(kind=RP) :: tmpx(1:nEqn, 0:fma%Nf(1), 0:fma%NfLeft(2))
      real(kind=RP) :: tmpy(1:nEqn, 0:fma%Nf(1), 0:fma%NfLeft(2))
      real(kind=RP) :: tmpz(1:nEqn, 0:fma%Nf(1), 0:fma%NfLeft(2))

      Uxe_rot=0.0_RP
      Uye_rot=0.0_RP
      Uze_rot=0.0_RP

      if ((fma % MortarType .ne. MORTAR_SMALL4) .and. .not.present(sliding) ) then 
         error stop 'MortarType SMALL4 reached in subroutine AdaptGradientsToMortarFace expecting non-SMALL4: check calling logic'
      end if 
      if (.not.present(sliding)) then 
         call GetMortarMInt(fma, MIntXi, MIntEta)
      else 
         call GetSlidingMInt(fma, MIntSliding)
      end if 
   
      if (.not.present(sliding)) then 

         associate(Uxf => fma % storage(1) % U_x, &
            Uyf => fma % storage(1) % U_y, &
            Uzf => fma % storage(1) % U_z   )
            Uxf=0.0_RP
            Uyf=0.0_RP
            Uzf=0.0_RP

            ! Pass 1: contraction direction 1 (xi) -> tmp(i,l)
            tmpx = 0.0_RP
            tmpy = 0.0_RP
            tmpz = 0.0_RP
            do l = 0, fma % NfLeft(2) ; do i = 0, fma % Nf(1) ; do m = 0, fma % NfLeft(1)
               tmpx(:,i,l) = tmpx(:,i,l) + MIntXi(i,m) * Uxe(:,m,l)
               tmpy(:,i,l) = tmpy(:,i,l) + MIntXi(i,m) * Uye(:,m,l)
               tmpz(:,i,l) = tmpz(:,i,l) + MIntXi(i,m) * Uze(:,m,l)
            end do ; end do ; end do

            ! Pass 2: contraction direction 2 (eta) -> Uxf/Uyf/Uzf(i,j)
            do j = 0, fma % Nf(2) ; do i = 0, fma % Nf(1) ; do l = 0, fma % NfLeft(2)
               Uxf(:,i,j) = Uxf(:,i,j) + MIntEta(j,l) * tmpx(:,i,l)
               Uyf(:,i,j) = Uyf(:,i,j) + MIntEta(j,l) * tmpy(:,i,l)
               Uzf(:,i,j) = Uzf(:,i,j) + MIntEta(j,l) * tmpz(:,i,l)
            end do ; end do ; end do

      end associate
      else  !sliding mortars

         if (side==1) then 
         associate(Uxf => fma % storage(side) % U_x, &
            Uyf => fma % storage(side) % U_y, &
            Uzf => fma % storage(side) % U_z   )
            Uxf=0.0_RP
            Uyf=0.0_RP
            Uzf=0.0_RP
            do l = 0, self % NfLeft(1)  ; do j = 0, self % Nf(1)   ; do i = 0, self % Nf(1)
               Uxf(:,i,j) = Uxf(:,i,j) + MIntSliding(j,l,side) * Uxe(:,i,l)
               Uyf(:,i,j) = Uyf(:,i,j) + MIntSliding(j,l,side) * Uye(:,i,l)
               Uzf(:,i,j) = Uzf(:,i,j) + MIntSliding(j,l,side) * Uze(:,i,l)
            end do                  ; end do                   ; end do
         end associate 
         else 
            do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
               call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), fma % rotation, ii, jj)
               Uxe_rot(:,i,j) = Uxe(:,ii,jj) 
               Uye_rot(:,i,j) = Uye(:,ii,jj) 
               Uze_rot(:,i,j) = Uze(:,ii,jj) 
            end do                        ; end do
            associate(Uxf => fma % storage(side) % U_x, &
               Uyf => fma % storage(side) % U_y, &
               Uzf => fma % storage(side) % U_z   )
               Uxf=0.0_RP
               Uyf=0.0_RP
               Uzf=0.0_RP
               do l = 0, self % NfLeft(1)  ; do j = 0, self % Nf(1)   ; do i = 0, self % Nf(1)
                  Uxf(:,i,j) = Uxf(:,i,j) + MIntSliding(j,l,side) * Uxe_rot(:,i,l)
                  Uyf(:,i,j) = Uyf(:,i,j) + MIntSliding(j,l,side) * Uye_rot(:,i,l)
                  Uzf(:,i,j) = Uzf(:,i,j) + MIntSliding(j,l,side) * Uze_rot(:,i,l)
               end do                  ; end do                   ; end do
            end associate 
         end if 
      end if 
   end subroutine Face_AdaptGradientsToMortarFace
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine Face_AdaptAviscFluxToFace(self, nEqn, Nelx, Nely, AVn_e, side)
      use MappedGeometryClass
      implicit none
      class(Face),   intent(inout)  :: self
      integer,       intent(in)     :: nEqn
      integer,       intent(in)     :: Nelx, Nely
      real(kind=RP), intent(in)     :: AVn_e(1:nEqn, 0:Nelx, 0:Nely)
      integer,       intent(in)     :: side
!
!     ---------------
!     Local variables
!     ---------------
!
      integer       :: i, j, k, l, m, ii, jj
      real(kind=RP) :: AVn_e_rot(1:nEqn, 0:self % NfRight(1), 0:self % NfRight(2))



      select case (side)
      case(1)
         if (self % MortarType == MORTAR_SMALL4) then 
            error stop 'MortarType SMALL4 reached in subroutine AdaptAviscFluxToFace expecting non-SMALL4: check calling logic'
         end if
         associate(AVf => self % storage(1) % AviscFlux)
         select case ( self % projectionType(1) )
         case (0)
            AVf = AVn_e
         case (1)
            AVf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfLeft(1)   ; do i = 0, self % Nf(1)
               AVf(:,i,j) = AVf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * AVn_e(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            AVf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               AVf(:,i,j) = AVf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * AVn_e(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            AVf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfLeft(1) ; do i = 0, self % Nf(1)
                  AVf(:,i,j) = AVf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                             * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                             * AVn_e(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select
         end associate
      case(2)
         associate( AVf => self % storage(2) % AviscFlux )
         do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
            call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
            AVn_e_rot(:,i,j) = AVn_e(:,ii,jj) 
         end do                        ; end do

         select case ( self % projectionType(2) )
         case (0)
            AVf = AVn_e_rot
         case (1)
            AVf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfRight(1)   ; do i = 0, self % Nf(1)
               AVf(:,i,j) = AVf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * AVn_e_rot(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            AVf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               AVf(:,i,j) = AVf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * AVn_e_rot(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            AVf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfRight(1) ; do i = 0, self % Nf(1)
                  AVf(:,i,j) = AVf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                             * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                             * AVn_e_rot(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select

         AVf = -AVf
         end associate
      end select

   end subroutine Face_AdaptAviscFluxToFace
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine Face_AdaptAviscFluxToMortarFace(self, nEqn, Nelx, Nely, AVn_e, side, fma, sliding)
      use MappedGeometryClass
      implicit none
      class(Face),   intent(inout)  :: self
      integer,       intent(in)     :: nEqn
      integer,       intent(in)     :: Nelx, Nely
      real(kind=RP), intent(in)     :: AVn_e(1:nEqn, 0:Nelx, 0:Nely)
      integer,       intent(in)     :: side
      type(Face), intent(inout)     :: fma
      logical,       intent(in), optional :: sliding
!
!     ---------------
!     Local variables
!     ---------------
!
      integer       :: i, j, l, m, ii, jj
      real(kind=RP) :: MIntXi(0:fma%Nf(1), 0:fma%NfLeft(1))
      real(kind=RP) :: MIntEta(0:fma%Nf(2), 0:fma%NfLeft(2))
      real(kind=RP) :: MIntSliding(0:fma%Nf(1), 0:fma%NfLeft(1), 1:2)

      real(kind=RP) :: AVn_e_rot(1:nEqn, 0:Nelx, 0:Nely)
      real(kind=RP) :: tmp(1:nEqn, 0:fma%Nf(1), 0:fma%NfLeft(2))

      AVn_e_rot = 0.0_RP

      if ((fma % MortarType .ne. MORTAR_SMALL4) .and. .not. present(sliding)) then
         error stop 'AdaptAviscFluxToMortarFace called on non-slave mortar without sliding: check calling logic'
      end if

      if (.not. present(sliding)) then
         call GetMortarMInt(fma, MIntXi, MIntEta)
      else
         call GetSlidingMInt(fma, MIntSliding)
      end if

      if (.not. present(sliding)) then

         associate(AVf => fma % storage(1) % AviscFlux)
            AVf = 0.0_RP

            ! Pass 1: contraction direction 1 (xi) -> tmp(i,l)
            tmp = 0.0_RP
            do l = 0, fma % NfLeft(2) ; do i = 0, fma % Nf(1) ; do m = 0, fma % NfLeft(1)
               tmp(:,i,l) = tmp(:,i,l) + MIntXi(i,m) * AVn_e(:,m,l)
            end do ; end do ; end do

            ! Pass 2: contraction direction 2 (eta) -> AVf(i,j)
            do j = 0, fma % Nf(2) ; do i = 0, fma % Nf(1) ; do l = 0, fma % NfLeft(2)
               AVf(:,i,j) = AVf(:,i,j) + MIntEta(j,l) * tmp(:,i,l)
            end do ; end do ; end do

         end associate

      else  !sliding mortars
         if (side==1) then
            associate(AVf => fma % storage(side) % AviscFlux)
               AVf = 0.0_RP
               do l = 0, self % NfLeft(1)  ; do j = 0, self % Nf(1)   ; do i = 0, self % Nf(1)
                  AVf(:,i,j) = AVf(:,i,j) + MIntSliding(j,l,side) * AVn_e(:,i,l)
               end do                  ; end do                   ; end do
            end associate
         else
            do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
               call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), fma % rotation, ii, jj)
               AVn_e_rot(:,i,j) = AVn_e(:,ii,jj)
            end do                        ; end do
            associate(AVf => fma % storage(side) % AviscFlux)
               AVf = 0.0_RP
               do l = 0, self % NfLeft(1)  ; do j = 0, self % Nf(1)   ; do i = 0, self % Nf(1)
                  AVf(:,i,j) = AVf(:,i,j) + MIntSliding(j,l,side) * AVn_e_rot(:,i,l)
               end do                  ; end do                   ; end do
            end associate
         end if
      end if

   end subroutine Face_AdaptAviscFluxToMortarFace
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine Face_ProjectFluxToElements(self, nEqn, flux, whichElements)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      class(Face)       :: self
      integer,       intent(in)  :: nEqn
      real(kind=RP), optional, intent(in)  :: flux(1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
      integer,       intent(in)  :: whichElements(2)
   
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
      integer           :: i, j, ii, jj, l, m, side, q, p, lm
      real(kind=RP)     :: fStarAux(nEqn, 0:self % NfRight(1), 0:self % NfRight(2))

      do side = 1, 2
         select case ( whichElements(side) )
         case (1)    ! Prolong to left element
            if (self % MortarType == MORTAR_NONE) then
               associate(fStar => self % storage(1) % Fstar)
               select case ( self % projectionType(1) )
               case (0)
                  fStar(1:nEqn,:,:) = flux
               case (1)
                  fStar(1:nEqn,:,:) = 0.0
                  do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
                     fStar(1:nEqn,i,j) = fStar(1:nEqn,i,j) + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * flux(:,l,j)
                  end do                  ; end do                   ; end do
                  
               case (2)
                  fStar(1:nEqn,:,:) = 0.0
                  do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
                     fStar(1:nEqn,i,j) = fStar(1:nEqn,i,j) + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * flux(:,i,l)
                  end do                  ; end do                   ; end do
         
               case (3)
                  fStar(1:nEqn,:,:) = 0.0
                  do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
                     do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                        fStar(1:nEqn,i,j) = fStar(1:nEqn,i,j) +   Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                               * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                               * flux(:,m,l)
                     end do                 ; end do
                  end do                  ; end do
               end select

               end associate
            end if
            
      
         case (2)    ! Prolong to right element
!      
!           *********
!           1st stage: Projection
!           *********
!           
         if (self % MortarType == MORTAR_NONE .OR. self % MortarType == MORTAR_SMALL4) then
            select case ( self % projectionType(2) )
            case (0)
               fStarAux(1:nEqn,:,:) = flux
            case (1)
               fStarAux(1:nEqn,:,:) = 0.0
               do j = 0, self % NfRight(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NfRight(1)
                  fStarAux(1:nEqn,i,j) = fStarAux(1:nEqn,i,j) + Tset(self % Nf(1), self % NfRight(1)) % T(i,l) * flux(:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               fStarAux(1:nEqn,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                  fStarAux(1:nEqn,i,j) = fStarAux(1:nEqn,i,j) + Tset(self % Nf(2), self % NfRight(2)) % T(j,l) * flux(:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               fStarAux(1:nEqn,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfRight(1)
                     fStarAux(1:nEqn,i,j) = fStarAux(1:nEqn,i,j) +   Tset(self % Nf(1), self % NfRight(1)) % T(i,m) &
                                                            * Tset(self % Nf(2), self % NfRight(2)) % T(j,l) &
                                                            * flux(:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
!      
!           *********
!           2nd stage: Rotation
!           *********
!      
            associate(fStar => self % storage(2) % Fstar)
            do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
               call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
               fStar(1:nEqn,ii,jj) = fStarAux(1:nEqn,i,j) 
            end do                        ; end do
!
!           *********
!           3rd stage: Inversion
!           *********
!
            
            fStar = -fStar
            end associate
         end if 

      end select
      end do

   end subroutine Face_ProjectFluxToElements
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine Face_ProjectMortarFluxToElements(self, nEqn, whichElements, slaveFace, MortarFlux, sliding)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      class(Face)       :: self !Master mortar face
      integer,       intent(in)  :: nEqn
      integer,       intent(in)  :: whichElements(2)
      type(Face), intent(inout)  :: slaveFace !Slave (small) mortar face
      real(kind=RP), intent(in)  :: MortarFlux(1:nEqn, 0:slaveFace % Nf(1), 0:slaveFace % Nf(2))
      logical, intent(in), optional     :: sliding 
   
      integer           :: i, j, ii, jj, l, m, side

      real(kind=RP) :: MoutXi(0:slaveFace%NfLeft(1), 0:slaveFace%Nf(1))
      real(kind=RP) :: MoutEta(0:slaveFace%NfLeft(2), 0:slaveFace%Nf(2))
      real(kind=RP) :: MoutSliding(0:slaveFace%NfLeft(1), 0:slaveFace%Nf(1), 1:2)

      real(kind=RP)     :: fStarAux(nEqn, 0:slaveFace % NfLeft(1), 0:slaveFace % NfLeft(2))
      real(kind=RP)     :: fStarAux2(nEqn, 0:slaveFace % NfRight(1), 0:slaveFace % NfRight(2))
      real(kind=RP)     :: tmp(nEqn, 0:slaveFace % NfLeft(1), 0:slaveFace % Nf(2))

      if (.not.present(sliding)) then 
         call GetMortarMout(slaveFace, MoutXi, MoutEta)
      else 
         call GetSlidingMout(slaveFace, MoutSliding)
      end if 

      if (.not.present(sliding)) then 
         fStarAux(1:nEqn,:,:) = 0.0_RP 

         ! Pass 1: contraction direction 1 (i <- m) -> tmp(i,l)
         tmp(1:nEqn,:,:) = 0.0_RP
         do l = 0, slaveFace % Nf(2) ; do i = 0, slaveFace % NfLeft(1) ; do m = 0, slaveFace % Nf(1)
            tmp(1:nEqn,i,l) = tmp(1:nEqn,i,l) + MoutXi(i,m) * MortarFlux(1:nEqn,m,l)
         end do ; end do ; end do

         ! Pass 2: contraction direction 2 (j <- l) -> fStarAux(i,j)
         do j = 0, slaveFace % NfLeft(2) ; do i = 0, slaveFace % NfLeft(1) ; do l = 0, slaveFace % Nf(2)
            fStarAux(1:nEqn,i,j) = fStarAux(1:nEqn,i,j) + MoutEta(j,l) * tmp(1:nEqn,i,l)
         end do ; end do ; end do

         associate(fStar => self % storage(1) % Fstar)
            fStar=fStar + fStarAux
         end associate 
         associate(fStar => self % storage(2) % Fstar)
            fStar=fStar + (-fStarAux)
         end associate 

      else !sliding mortar
         fStarAux(1:nEqn,:,:) = 0.0_RP 
         if (whichElements(1)==1) then 
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               fStarAux(1:nEqn,i,j)= fStarAux(1:nEqn,i,j) + MoutSliding(j,l,1) * MortarFlux(1:nEqn,i,l)
            end do                  ; end do                   ; end do
            associate(fStar => self % storage(1) % Fstar)
               fStar=fStar + slaveFace%s(1) * fStarAux
            end associate 
            associate(fStar => self % storage(2) % Fstar)
               fStar=fStar + slaveFace%s(1) * (-fStarAux)
            end associate 

         else 
            fStarAux(1:nEqn,:,:) = 0.0_RP 
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               fStarAux(1:nEqn,i,j)= fStarAux(1:nEqn,i,j) + MoutSliding(j,l,2) * MortarFlux(1:nEqn,i,l)
            end do                  ; end do                   ; end do

            fStarAux2=0.0_RP
            fStarAux=slaveFace%s(1) * (fStarAux)
            associate(fStar => self % storage(2) % Fstar)
               do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                  call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), slaveFace % rotation, ii, jj)
                  fStarAux2(1:nEqn,ii,jj) = fStarAux(1:nEqn,i,j) 
               end do                        ; end do
               fStar=fStar+(-fStarAux2)
            end associate
            fStarAux2=0.0_RP

            associate(fStar => self % storage(1) % Fstar)
               do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                  call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), slaveFace % rotation, ii, jj)
                  fStarAux2(1:nEqn,ii,jj) = fStarAux(1:nEqn,i,j) 
               end do                        ; end do
               fStar=fStar+(-fStarAux2)
            end associate

         end if 
      end if  
   end subroutine Face_ProjectMortarFluxToElements
#if defined(NAVIERSTOKES) 
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!  Face_ProjectBCJacobianToElements:
!  Routine to project the Jacobian of the boundary condition to the element
!  --
   subroutine Face_ProjectBCJacobianToElements(self,nEqn,BCjacF)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      class(Face), target        :: self
      integer,       intent(in)  :: nEqn
      real(kind=RP), intent(in)  :: BCjacF(nEqn,nEqn,0:self % Nf(1), 0:self % Nf(2))
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                :: i, j, l, m
!     Always project to left element (this is only for boundaries)
!     ------------------------------------------------------------
      associate(BCJac => self % storage(1) % BCJac)
      
      select case ( self % projectionType(1) )
         case (0)
            BCJac(1:nEqn,1:nEqn,:,:) = BCJacF
         case (1)
            BCJac(1:nEqn,1:nEqn,:,:) = 0._RP
            do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
               BCJac(1:nEqn,1:nEqn,i,j) = BCJac(1:nEqn,1:nEqn,i,j) &
                                                            + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * BCJacF(:,:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            BCJac(1:nEqn,1:nEqn,:,:) = 0._RP
            do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
               BCJac(1:nEqn,1:nEqn,i,j) = BCJac(1:nEqn,1:nEqn,i,j) &
                                                            + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * BCJacF(:,:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            BCJac(1:nEqn,1:nEqn,:,:) = 0._RP
            do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
               do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                  BCJac(1:nEqn,1:nEqn,i,j) = BCJac(1:nEqn,1:nEqn,i,j) &
                                                            +  Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                            * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                            * BCJacF(:,:,m,l)
               end do                 ; end do
            end do                  ; end do
      end select
      
      end associate
   end subroutine Face_ProjectBCJacobianToElements
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!            
   subroutine Face_ProjectFluxJacobianToElements(self, nEqn, whichElement,whichderiv)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      class(Face), target        :: self
      integer,       intent(in)  :: nEqn
      integer,       intent(in)  :: whichElement
      integer,       intent(in)  :: whichderiv           !<  One can either transfer the derivative with respect to qL (LEFT) or to qR (RIGHT)
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                :: i, j, ii, jj, l, m, side, lm, a, b
      real(kind=RP), pointer :: fluxDeriv(:,:,:,:)
      real(kind=RP)          :: fStarAux(nEqn,nEqn, 0:self % NfRight(1), 0:self % NfRight(2))

      fluxDeriv(1:,1:,0:,0:) => self % storage(whichderiv) % dFStar_dqF
      select case ( whichElement )
         case (LEFT)    ! Prolong to left element
            if (self % MortarType == MORTAR_NONE) then 
               associate(dFStar_dq => self % storage(1) % dFStar_dqEl)
               select case ( self % projectionType(1) )
               case (0)
                  dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) = fluxDeriv
               case (1)
                  dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) = 0._RP
                  do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
                     dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) = dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) &
                                                                  + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * fluxDeriv(:,:,l,j)
                  end do                  ; end do                   ; end do
                  
               case (2)
                  dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) = 0._RP
                  do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
                     dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) = dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) &
                                                                  + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * fluxDeriv(:,:,i,l)
                  end do                  ; end do                   ; end do
         
               case (3)
                  dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) = 0._RP
                  do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
                     do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                        dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) = dFStar_dq(1:nEqn,1:nEqn,i,j,whichderiv) &
                                                               +  Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                                  * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                                  * fluxDeriv(:,:,m,l)
                     end do                 ; end do
                  end do                  ; end do
               end select
               end associate
            end if 

         case (RIGHT)    ! Prolong to right element
!      
!           *********
!           1st stage: Projection
!           *********
!       
            if (self % MortarType == MORTAR_NONE .OR. self % MortarType == MORTAR_SMALL4) then
               select case ( self % projectionType(2) )
               case (0)
                  fStarAux(1:nEqn,1:nEqn,:,:) = fluxDeriv
               case (1)
                  fStarAux(1:nEqn,1:nEqn,:,:) = 0.0
                  do j = 0, self % NfRight(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NfRight(1)
                     fStarAux(1:nEqn,1:nEqn,i,j) = fStarAux(1:nEqn,1:nEqn,i,j) + Tset(self % Nf(1), self % NfRight(1)) % T(i,l) * fluxDeriv(:,:,l,j)
                  end do                  ; end do                   ; end do
                  
               case (2)
                  fStarAux(1:nEqn,1:nEqn,:,:) = 0.0
                  do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                     fStarAux(1:nEqn,1:nEqn,i,j) = fStarAux(1:nEqn,1:nEqn,i,j) + Tset(self % Nf(2), self % NfRight(2)) % T(j,l) * fluxDeriv(:,:,i,l)
                  end do                  ; end do                   ; end do
         
               case (3)
                  fStarAux(1:nEqn,1:nEqn,:,:) = 0.0
                  do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   
                     do m = 0, self % Nf(1) ; do i = 0, self % NfRight(1)
                        fStarAux(1:nEqn,1:nEqn,i,j) = fStarAux(1:nEqn,1:nEqn,i,j) +   Tset(self % Nf(1), self % NfRight(1)) % T(i,m) &
                                                               * Tset(self % Nf(2), self % NfRight(2)) % T(j,l) &
                                                               * fluxDeriv(:,:,m,l)
                     end do                 ; end do
                  end do                  ; end do
               end select
!      
!           *********
!           2nd stage: Rotation
!           *********
!      
               associate(dFStar_dq => self % storage(2) % dFStar_dqEl)
               do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                  call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
                  dFStar_dq(1:nEqn,1:nEqn,ii,jj,whichderiv) = fStarAux(1:nEqn,1:nEqn,i,j) 
               end do                        ; end do
!
!           *********
!           3rd stage: Inversion
!           *********
!
               dFStar_dq = -dFStar_dq
               end associate
            end if 
      end select
      

   end subroutine Face_ProjectFluxJacobianToElements
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine Face_ProjectMortarFluxJacobianToElements(self, nEqn, whichElement,whichderiv, slaveFace, sliding) 
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      class(Face), target        :: self
      integer,       intent(in)  :: nEqn
      integer,       intent(in)  :: whichElement
      integer,       intent(in)  :: whichderiv !<  One can either transfer the derivative with respect to qL (LEFT) or to qR (RIGHT)
      type(Face), intent(inout)      :: slaveFace
      logical, intent(in), optional  :: sliding 

!
!     ---------------
!     Local variables
!     ---------------
!
      integer                :: i, j, l, m
      real(kind=RP), pointer :: fluxDeriv(:,:,:,:)
      real(kind=RP)     :: dFStar_dqFAux(nEqn, nEqn, 0: slaveFace % NfLeft(1),0:slaveFace % NfLeft(2))
      real(kind=RP)     :: Flux_tmp(nEqn,nEqn,0: slaveFace % NfLeft(1),0:slaveFace % Nf(2))
      real(kind=RP) :: MoutXi(0:slaveFace%NfLeft(1), 0:slaveFace%Nf(1))
      real(kind=RP) :: MoutEta(0:slaveFace%NfLeft(2), 0:slaveFace%Nf(2))

      !if (.not.present(sliding)) then 
         call GetMortarMout(slaveFace, MoutXi, MoutEta)
      !else 
         !call GetSlidingMout(slaveFace, MoutSliding)
      !end if 

         fluxDeriv(1:,1:,0:,0:) => self % storage(whichderiv) % dFStar_dqF
         dFStar_dqFAux=0.0_RP 

        ! Pass 1: contraction direction 1 (i <- m) -> Flux_tmp(i,l)
         Flux_tmp(:,:,:,:)=0.0_RP
         do l = 0, slaveFace % Nf(2) ; do i = 0, slaveFace % NfLeft(1) ; do m = 0, slaveFace % Nf(1)
           Flux_tmp(1:nEqn,1:nEqn,i,l) = Flux_tmp(1:nEqn,1:nEqn,i,l)  + MoutXi(i,m) * fluxDeriv(1:nEqn,1:nEqn,m,l)
        end do ; end do ; end do

        ! Pass 2: contraction direction 2 (j <- l) -> fStarAux(i,j)
        do j = 0, slaveFace % NfLeft(2) ; do i = 0, slaveFace % NfLeft(1) ; do l = 0, slaveFace % Nf(2)
           dFStar_dqFAux(1:nEqn,1:nEqn,i,j) = dFStar_dqFAux(1:nEqn,1:nEqn,i,j) + MoutEta(j,l) * Flux_tmp(1:nEqn,1:nEqn,i,l) 
        end do ; end do ; end do

      select case (whichElement)
      case (1)
         associate(dFStar_dq => self % storage(whichElement) % dFStar_dqEl)
            dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv)=dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) + &
            dFStar_dqFAux(1:nEqn,1:nEqn,:,:)
        end associate 
      case (2)
         associate(dFStar_dq => self % storage(whichElement) % dFStar_dqEl)
            dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv)=dFStar_dq(1:nEqn,1:nEqn,:,:,whichderiv) + &
            (-dFStar_dqFAux(1:nEqn,1:nEqn,:,:))
         end associate 
      end select 

   end subroutine Face_ProjectMortarFluxJacobianToElements
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#endif
   subroutine Face_ProjectGradientFluxToElements(self, nEqn, Hflux, whichElements, factor)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      class(Face)       :: self
      integer,       intent(in)  :: nEqn
      real(kind=RP), intent(in), optional  :: Hflux(nEqn, NDIM, 0:self % Nf(1), 0:self % Nf(2))
      integer,       intent(in)  :: whichElements(2)
      integer,       intent(in)  :: factor               ! A factor that relates LEFT and RIGHT fluxes

!
!     ---------------
!     Local variables
!     ---------------
!
      integer           :: i, j, ii, jj, l, m, side,p,q,lm
      real(kind=RP)     :: hStarAux(nEqn, NDIM, 0:self % NfRight(1), 0:self % NfRight(2))
      real(kind=RP), allocatable     :: Flux_tmp(:,:,:,:)

   do side = 1, 2
      select case ( whichElements(side) )
      case (1)    ! Prolong from left element
         if (self % MortarType == MORTAR_NONE) then
         select case ( self % projectionType(1) )
         case (0)
            self % storage(1) % unStar(:,:,:,:) = Hflux
         case (1)
            self % storage(1) % unStar(:,:,:,:) = 0.0
            do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
               self % storage(1) % unStar(:,:,i,j) = self % storage(1) % unStar(:,:,i,j) + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * Hflux(:,:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            self % storage(1) % unStar(:,:,:,:) = 0.0
            do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
               self % storage(1) % unStar(:,:,i,j) = self % storage(1) % unStar(:,:,i,j) + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * Hflux(:,:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            self % storage(1) % unStar(:,:,:,:) = 0.0
            do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
               do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                  self % storage(1) % unStar(:,:,i,j) = self % storage(1) % unStar(:,:,i,j) +   Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                            * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                            * Hflux(:,:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select
         end if ! MortarType == 0

      case (2)    ! Prolong from right element
!      
!           *********
!           1st stage: Projection
!           *********
            if (self % MortarType == MORTAR_NONE .OR. self % MortarType == MORTAR_SMALL4) then 
            select case ( self % projectionType(2) )
            case (0)
               HstarAux = Hflux
            case (1)
               HstarAux = 0.0
               do j = 0, self % NfRight(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NfRight(1)
                  HstarAux(:,:,i,j) = HstarAux(:,:,i,j) + Tset(self % Nf(1), self % NfRight(1)) % T(i,l) * Hflux(:,:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               HstarAux = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                  HstarAux(:,:,i,j) = HstarAux(:,:,i,j) + Tset(self % Nf(2), self % NfRight(2)) % T(j,l) * Hflux(:,:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               HstarAux = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfRight(1)
                     HstarAux(:,:,i,j) = HstarAux(:,:,i,j) +   Tset(self % Nf(1), self % NfRight(1)) % T(i,m) &
                                                             * Tset(self % Nf(2), self % NfRight(2)) % T(j,l) &
                                                             * Hflux(:,:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
!      
!           *********
!           2nd stage: Rotation
!           *********
!      
            do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
               call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
               self % storage(2) % unStar(:,:,ii,jj) = HstarAux(:,:,i,j) 
            end do                       ; end do
!
!           *********
!           3rd stage: Multiplication by a factor (inversion usually)
!           *********
!
            self % storage(2) % unStar = factor * self % storage(2) % unStar

         end if 
         end select
      end do

   end subroutine Face_ProjectGradientFluxToElements
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine Face_ProjectMortarGradientFluxToElements(self, nEqn, slaveFace, Hflux, whichElements, factor, sliding)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      class(Face)       :: self
      integer,       intent(in)  :: nEqn
      type(Face), intent(inout)  :: slaveFace
      real(kind=RP), intent(in), optional  :: Hflux(nEqn, NDIM, 0:slaveFace % Nf(1), 0:slaveFace % Nf(2))
      integer,       intent(in)  :: whichElements(2)
      integer,       intent(in)  :: factor
      logical, optional, intent(in) :: sliding
   
      real(kind=RP)     :: hStarAux(nEqn, NDIM, 0:slaveFace % NfLeft(1), 0:slaveFace % NfLeft(2))
      real(kind=RP)     :: hStarAux2(nEqn, NDIM, 0:slaveFace % NfRight(1), 0:slaveFace % NfRight(2))
   
      real(kind=RP)     :: MoutXi(0:slaveFace%NfLeft(1), 0:slaveFace%Nf(1))
      real(kind=RP)     :: MoutEta(0:slaveFace%NfLeft(2), 0:slaveFace%Nf(2))
      real(kind=RP)     :: MoutSliding(0:slaveFace%NfLeft(1), 0:slaveFace%Nf(1), 1:2)
   
      real(kind=RP)     :: tmp(nEqn, NDIM, 0:slaveFace % NfLeft(1), 0:slaveFace % Nf(2))
      integer           :: i, j, l, m, ii, jj
   
      if (.not.present(sliding)) then
         call GetMortarMout(slaveFace, MoutXi, MoutEta)
      else
         call GetSlidingMout(slaveFace, MoutSliding)
      end if
   
      hStarAux = 0.0_RP
      if (.not.present(sliding)) then
   
         ! Pass 1: contraction direction 1 (i <- m) -> tmp(i,l)
         tmp = 0.0_RP
         do l = 0, slaveFace % Nf(2) ; do i = 0, slaveFace % NfLeft(1) ; do m = 0, slaveFace % Nf(1)
            tmp(:,:,i,l) = tmp(:,:,i,l) + MoutXi(i,m) * Hflux(:,:,m,l)
         end do ; end do ; end do
   
         ! Pass 2: contraction direction 2 (j <- l) -> hStarAux(i,j)
         do j = 0, slaveFace % NfLeft(2) ; do i = 0, slaveFace % NfLeft(1) ; do l = 0, slaveFace % Nf(2)
            hStarAux(:,:,i,j) = hStarAux(:,:,i,j) + MoutEta(j,l) * tmp(:,:,i,l)
         end do ; end do ; end do
   
         associate(unStar => self % storage(1) % unStar)
               unStar=unStar+hStarAux
         end associate
   
         associate(unStar => self % storage(2) % unStar)
               unStar=unStar+ factor * hStarAux
         end associate
      else  !sliding
      if (whichElements(1)==1) then
         hStarAux = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
            hStarAux(:,:,i,j)= hStarAux(:,:,i,j) + MoutSliding(j,l,whichElements(1)) * Hflux(:,:,i,l)
         end do                  ; end do                   ; end do
   
         hStarAux=slaveFace%s(1) * (hStarAux)
         associate(unStar => self % storage(1) % unStar)
   
         end associate
            end if
            if (whichElements(1)==2) then
               hStarAux = 0.0_RP
               do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               hStarAux(:,:,i,j)= hStarAux(:,:,i,j) + MoutSliding(j,l,whichElements(1)) * Hflux(:,:,i,l)
                  end do                  ; end do                   ; end do
                  hStarAux=slaveFace%s(1) * (hStarAux)
                  associate(unStar => self % storage(2) % unStar)
                     do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                        call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), slaveFace % rotation, ii, jj)
                        hStarAux2(:,:,ii,jj) = hStarAux(:,:,i,j)
                     end do                        ; end do
                     unStar=unStar+factor*(hStarAux2)
                  end associate
   
            end if
      end if
   
   end subroutine Face_ProjectMortarGradientFluxToElements
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine Face_Interpolatebig2small(self, nEqn, fma, grad)
      use MappedGeometryClass
      implicit none
      CLASS(Face),   intent(inout)              :: self
      integer,       intent(in)                 :: nEqn
      type(Face), intent(inout)                 ::fma
      integer, optional                         :: grad

#ifdef _HAS_MPI_

      integer       :: i, j, l, m
      real(kind=RP) :: Um(1:nEqn, 0:fma%Nf(1), 0:fma%Nf(2))
      real(kind=RP) :: Qfm(1:nEqn, 0:fma%Nf(1), 0:fma%Nf(2))

      real(kind=RP) :: MIntXi(0:fma%Nf(1), 0:fma%NfLeft(1))
      real(kind=RP) :: MIntEta(0:fma%Nf(2), 0:fma%NfLeft(2))

      real(kind=RP) :: tmp(1:nEqn, 0:fma%Nf(1), 0:fma%NfLeft(2))

      call GetMortarMInt(fma, MIntXi, MIntEta)
      
      if (present(grad)) then 
         select case(grad)
         case(1)
         associate(Uxf => fma % storage(1) % U_x)
            Um=Uxf
            Uxf=0.0_RP

            tmp = 0.0_RP
            do l = 0, fma % NfLeft(2) ; do i = 0, fma % Nf(1) ; do m = 0, fma % NfLeft(1)
               tmp(:,i,l) = tmp(:,i,l) + MIntXi(i,m) * Um(:,m,l)
            end do ; end do ; end do

            do j = 0, fma % Nf(2) ; do i = 0, fma % Nf(1) ; do l = 0, fma % NfLeft(2)
               Uxf(:,i,j) = Uxf(:,i,j) + MIntEta(j,l) * tmp(:,i,l)
            end do ; end do ; end do
            end associate 
         case(2)
         associate(Uyf => fma % storage(1) % U_y)
            Um=Uyf
            Uyf=0.0_RP

            tmp = 0.0_RP
            do l = 0, fma % NfLeft(2) ; do i = 0, fma % Nf(1) ; do m = 0, fma % NfLeft(1)
               tmp(:,i,l) = tmp(:,i,l) + MIntXi(i,m) * Um(:,m,l)
            end do ; end do ; end do

            do j = 0, fma % Nf(2) ; do i = 0, fma % Nf(1) ; do l = 0, fma % NfLeft(2)
               Uyf(:,i,j) = Uyf(:,i,j) + MIntEta(j,l) * tmp(:,i,l)
            end do ; end do ; end do
         end associate 
         case(3)
         associate(Uzf => fma % storage(1) % U_z)
            Um=Uzf
            Uzf=0.0_RP

            tmp = 0.0_RP
            do l = 0, fma % NfLeft(2) ; do i = 0, fma % Nf(1) ; do m = 0, fma % NfLeft(1)
               tmp(:,i,l) = tmp(:,i,l) + MIntXi(i,m) * Um(:,m,l)
            end do ; end do ; end do

            do j = 0, fma % Nf(2) ; do i = 0, fma % Nf(1) ; do l = 0, fma % NfLeft(2)
               Uzf(:,i,j) = Uzf(:,i,j) + MIntEta(j,l) * tmp(:,i,l)
            end do ; end do ; end do
         end associate 
         end select 
      else 
            associate(Qf => fma % storage(1) % Q)
               Qfm=Qf
               Qf=0.0_RP

               tmp = 0.0_RP
               do l = 0, fma % NfLeft(2) ; do i = 0, fma % Nf(1) ; do m = 0, fma % NfLeft(1)
                  tmp(:,i,l) = tmp(:,i,l) + MIntXi(i,m) * Qfm(:,m,l)
               end do ; end do ; end do

               do j = 0, fma % Nf(2) ; do i = 0, fma % Nf(1) ; do l = 0, fma % NfLeft(2)
                  Qf(:,i,j) = Qf(:,i,j) + MIntEta(j,l) * tmp(:,i,l)
               end do ; end do ; end do

            end associate 

      end if 

#endif
   end subroutine Face_Interpolatebig2small
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   SUBROUTINE Face_Interpolatesmall2big(self, nEqn, flux_M, anJacobian)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      class(Face)       :: self
      integer,       intent(in)  :: nEqn
      real(kind=RP), intent(inout), optional  :: flux_M(1:nEqn, 0:self % Nf(1), 0:self % Nf(2))
      logical, intent(in), optional :: anJacobian

#ifdef _HAS_MPI_
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
      integer           :: i, j, l, m
      real(kind=RP)     :: fStarAux(nEqn, 0:self % NfLeft(1), 0:self % NfLeft(2))
      real(kind=RP)     :: dFStar_dqFAux(nEqn, nEqn, 0:self % NfLeft(1), 0:self % NfLeft(2))

      real(kind=RP)     :: MoutXi(0: self % NfLeft(1), 0:self % Nf(1))
      real(kind=RP)     :: MoutEta(0: self % NfLeft(2), 0:self % Nf(2))
      real(kind=RP)     :: tmp(nEqn, 0:self % NfLeft(1), 0:self % Nf(2))
      real(kind=RP)     :: Flux_tmp(nEqn,nEqn,0: self % NfLeft(1),0:self % Nf(2))

      call GetMortarMout(self, MoutXi, MoutEta)

      if (.not.present(anJacobian)) then 
         fStarAux(1:nEqn,:,:) = 0.0_RP 

         ! Pass 1: contraction direction 1 (i <- m) -> tmp(i,l)
         tmp(1:nEqn,:,:) = 0.0_RP
         do l = 0, self % Nf(2) ; do i = 0, self % NfLeft(1) ; do m = 0, self % Nf(1)
            tmp(1:nEqn,i,l) = tmp(1:nEqn,i,l) + MoutXi(i,m) * flux_M(1:nEqn,m,l)
         end do ; end do ; end do

         ! Pass 2: contraction direction 2 (j <- l) -> fStarAux(i,j)
         do j = 0, self % NfLeft(2) ; do i = 0, self % NfLeft(1) ; do l = 0, self % Nf(2)
            fStarAux(1:nEqn,i,j) = fStarAux(1:nEqn,i,j) + MoutEta(j,l) * tmp(1:nEqn,i,l)
         end do ; end do ; end do

         associate(Mflux => self % storage(1) % MortarFlux)
            Mflux=fStarAux
         end associate 

      else 
#if defined(NAVIERSTOKES)
         dFStar_dqFAux=0.0_RP 

         ! Pass 1: contraction direction 1 (i <- m) -> Flux_tmp(i,l)
         Flux_tmp(:,:,:,:)=0.0_RP
         do l = 0, self % Nf(2) ; do i = 0, self % NfLeft(1) ; do m = 0, self % Nf(1)
            Flux_tmp(1:nEqn,1:nEqn,i,l) = Flux_tmp(1:nEqn,1:nEqn,i,l)  + MoutXi(i,m) * self % storage(1) % dFStar_dqF(1:nEqn,1:nEqn,m,l)
         end do ; end do ; end do

         ! Pass 2: contraction direction 2 (j <- l) -> fStarAux(i,j)
         do j = 0, self % NfLeft(2) ; do i = 0, self % NfLeft(1) ; do l = 0, self % Nf(2)
            dFStar_dqFAux(1:nEqn,1:nEqn,i,j) = dFStar_dqFAux(1:nEqn,1:nEqn,i,j) + MoutEta(j,l) * Flux_tmp(1:nEqn,1:nEqn,i,l) 
         end do ; end do ; end do

         associate(MfluxJ => self % storage(1) % MortarFlux_J)
            MfluxJ=dFStar_dqFAux
         end associate 
#endif
      end if 
#endif
   END SUBROUTINE Face_Interpolatesmall2big
   

   SUBROUTINE Face_Interpolatesmall2biggrad(self, nEqn, Hflux, anJacobian)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      class(Face)       :: self
      integer,       intent(in)  :: nEqn
      real(kind=RP), intent(in)     :: Hflux(nEqn, NDIM, 0:self % Nf(1), 0:self % Nf(2))
      logical, optional, intent(in) :: anJacobian
#ifdef _HAS_MPI_
   !
   !     ---------------
   !     Local variables
   !     ---------------
   !
      integer           :: i, j, l, m
      real(kind=RP)     :: hStarAux(nEqn, NDIM, 0:self % NfLeft(1), 0:self % NfLeft(2))
      real(kind=RP)     :: MoutXi(0: self % NfLeft(1), 0:self % Nf(1))
      real(kind=RP)     :: MoutEta(0: self % NfLeft(2), 0:self % Nf(2))
      real(kind=RP)     :: tmp(nEqn, NDIM, 0:self % NfLeft(1), 0:self % Nf(2))

      real(kind=RP)     :: Flux_tmp(NCONS,NCONS,1:NDIM,0: self % NfLeft(1),0:self % Nf(2))
      real(kind=RP)     :: dFv_dGradQElAux(NCONS,NCONS,1:NDIM, 0:self % NfLeft(1),0:self % NfLeft(2))

      call GetMortarMout(self, MoutXi, MoutEta)

      if (.not.present(anJacobian)) then 
         hStarAux = 0.0_RP 

         ! Pass 1: contraction direction 1 (i <- m) -> tmp(i,l)
         tmp = 0.0_RP
         do l = 0, self % Nf(2) ; do i = 0, self % NfLeft(1) ; do m = 0, self % Nf(1)
            tmp(:,:,i,l) = tmp(:,:,i,l) + MoutXi(i,m) * Hflux(:,:,m,l)
         end do ; end do ; end do

         ! Pass 2: contraction direction 2 (j <- l) -> hStarAux(i,j)
         do j = 0, self % NfLeft(2) ; do i = 0, self % NfLeft(1) ; do l = 0, self % Nf(2)
            hStarAux(:,:,i,j) = hStarAux(:,:,i,j) + MoutEta(j,l) * tmp(:,:,i,l)
         end do ; end do ; end do

         associate(MunStar => self % storage(1) %GradMortarFlux)
            MunStar=hStarAux
         end associate 

      else

#if defined(NAVIERSTOKES)

         dFv_dGradQElAux=0.0_RP 

         ! Pass 1: contraction direction 1 (i <- m) -> Flux_tmp(i,l)
         Flux_tmp(:,:,:,:,:)=0.0_RP
         do l = 0, self % Nf(2) ; do i = 0, self % NfLeft(1) ; do m = 0, self % Nf(1)
         Flux_tmp(1:NCONS,1:NCONS, 1:NDIM, i,l) = Flux_tmp(1:NCONS,1:NCONS, 1:NDIM, i,l)  + MoutXi(i,m) * self % storage(1) % dFv_dGradQF(1:NCONS,1:NCONS,1:NDIM,m,l)
         end do ; end do ; end do

         ! Pass 2: contraction direction 2 (j <- l) -> fStarAux(i,j)
         do j = 0, self % NfLeft(2) ; do i = 0, self % NfLeft(1) ; do l = 0, self % Nf(2)
            dFv_dGradQElAux(1:NCONS,1:NCONS,1:NDIM,i,j) = dFv_dGradQElAux(1:NCONS,1:NCONS,1:NDIM,i,j) + MoutEta(j,l) * Flux_tmp(1:NCONS,1:NCONS, 1:NDIM,i,l) 
         end do ; end do ; end do

         associate(MfluxJ => self % storage(1) % GradMortarFlux_J)
            MfluxJ=dFv_dGradQElAux
         end associate 
#endif
      end if 
#endif
   END SUBROUTINE Face_Interpolatesmall2biggrad
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES)
   subroutine Face_ProjectGradJacobianToElements(self, whichElement, whichderiv)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      !---------------------------------------------------------
      class(Face), target        :: self
      integer,       intent(in)  :: whichElement
      integer,       intent(in)  :: whichderiv
      !---------------------------------------------------------
      integer                :: i, j, ii, jj, l, m, side, lm, a, b
      real(kind=RP), pointer :: fluxDeriv(:,:,:,:,:)
      real(kind=RP)          :: fStarAux(NCONS,NCONS,1:NDIM, 0:self % NfRight(1), 0:self % NfRight(2))

      
      fluxDeriv(1:,1:,1:,0:,0:) => self % storage(whichderiv) % dFv_dGradQF
      select case ( whichElement )
         case (LEFT)    ! Prolong to left element
         if (self % MortarType == MORTAR_NONE) then 
            associate(dFv_dGradQEl => self % storage(1) % dFv_dGradQEl)
            
            select case ( self % projectionType(1) )
            case (0)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv) = fluxDeriv
            case (1)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv) = 0._RP
               do j = 0, self % NelLeft(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NelLeft(1)
                  dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) = dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) &
                                                                           + Tset(self % Nf(1), self % NfLeft(1)) % T(i,l) * fluxDeriv(:,:,:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv) = 0._RP
               do l = 0, self % Nf(2)  ; do j = 0, self % NelLeft(2)   ; do i = 0, self % NelLeft(1)
                  dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) = dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) &
                                                                           + Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) * fluxDeriv(:,:,:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv) = 0._RP
               do l = 0, self % Nf(2)  ; do j = 0, self % NfLeft(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfLeft(1)
                     dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) = dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,i,j,whichderiv) &
                                                                                    +  Tset(self % Nf(1), self % NfLeft(1)) % T(i,m) &
                                                                                    * Tset(self % Nf(2), self % NfLeft(2)) % T(j,l) &
                                                                                    * fluxDeriv(:,:,:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
            
            end associate
         end if 

         case (RIGHT)    ! Prolong to right element
         if (self % MortarType == MORTAR_NONE .OR. self % MortarType == MORTAR_SMALL4) then
!      
!           *********
!           1st stage: Projection
!           *********
!      
            select case ( self % projectionType(2) )
            case (0)
               fStarAux(1:NCONS,1:NCONS,1:NDIM,:,:) = fluxDeriv
            case (1)
               fStarAux(1:NCONS,1:NCONS,1:NDIM,:,:) = 0.0
               do j = 0, self % NfRight(2)  ; do l = 0, self % Nf(1)   ; do i = 0, self % NfRight(1)
                  fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) = fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) + Tset(self % Nf(1), self % NfRight(1)) % T(i,l) * fluxDeriv(:,:,:,l,j)
               end do                  ; end do                   ; end do
               
            case (2)
               fStarAux(1:NCONS,1:NCONS,1:NDIM,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
                  fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) = fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) + Tset(self % Nf(2), self % NfRight(2)) % T(j,l) * fluxDeriv(:,:,:,i,l)
               end do                  ; end do                   ; end do
      
            case (3)
               fStarAux(1:NCONS,1:NCONS,1:NDIM,:,:) = 0.0
               do l = 0, self % Nf(2)  ; do j = 0, self % NfRight(2)   
                  do m = 0, self % Nf(1) ; do i = 0, self % NfRight(1)
                     fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) = fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) +   Tset(self % Nf(1), self % NfRight(1)) % T(i,m) &
                                                               * Tset(self % Nf(2), self % NfRight(2)) % T(j,l) &
                                                               * fluxDeriv(:,:,:,m,l)
                  end do                 ; end do
               end do                  ; end do
            end select
!      
!           ********* 
!           2nd stage: Rotation
!           *********
!      
            associate(dFv_dGradQEl => self % storage(2) % dFv_dGradQEl)
            
            do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
               call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
               dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,ii,jj,whichderiv) = fStarAux(1:NCONS,1:NCONS,1:NDIM,i,j) 
            end do                        ; end do
!
!           *********
!           3rd stage: Inversion
!           *********
!
            dFv_dGradQEl = -dFv_dGradQEl
            
            end associate
         end if 
      end select
      
      
   end subroutine Face_ProjectGradJacobianToElements
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine Face_ProjectMortarGradJacobianToElements(self, whichElement, whichderiv, slaveFace, sliding)
      use MappedGeometryClass
      use PhysicsStorage
      implicit none
      !---------------------------------------------------------
      class(Face), target        :: self
      integer,       intent(in)  :: whichElement
      integer,       intent(in)  :: whichderiv !<  One can either transfer the derivative with respect to qL (LEFT) or to qR (RIGHT)
      type(Face), intent(inout)      :: slaveFace
      logical, intent(in), optional  :: sliding 

!
!     ---------------
!     Local variables
!     ---------------
!
      integer                :: i, j, l, m
      real(kind=RP), pointer :: fluxDeriv(:,:,:,:,:)

      real(kind=RP)     :: Flux_tmp(NCONS,NCONS,1:NDIM,0:slaveFace % NfLeft(1), 0:slaveFace % Nf(2))
      real(kind=RP)     :: dFv_dGradQElAux(NCONS,NCONS,1:NDIM, 0:slaveFace % NfLeft(1), 0:slaveFace % NfLeft(2))
      real(kind=RP)     :: MoutXi(0:slaveFace%NfLeft(1), 0:slaveFace%Nf(1))
      real(kind=RP)     :: MoutEta(0:slaveFace%NfLeft(2), 0:slaveFace%Nf(2))

      !if (.not.present(sliding)) then 
         call GetMortarMout(slaveFace, MoutXi, MoutEta)
      !else 
         !call GetSlidingMout(slaveFace, MoutSliding)
      !end if 

         fluxDeriv(1:,1:,1:,0:,0:) => self % storage(whichderiv) % dFv_dGradQF
         dFv_dGradQElAux=0.0_RP 

        ! Pass 1: contraction direction 1 (i <- m) -> Flux_tmp(i,l)
         Flux_tmp(:,:,:,:,:)=0.0_RP
         do l = 0, slaveFace % Nf(2) ; do i = 0, slaveFace % NfLeft(1) ; do m = 0, slaveFace % Nf(1)
           Flux_tmp(1:NCONS,1:NCONS, 1:NDIM, i,l) = Flux_tmp(1:NCONS,1:NCONS, 1:NDIM, i,l)  + MoutXi(i,m) * fluxDeriv(1:NCONS,1:NCONS,1:NDIM,m,l)
        end do ; end do ; end do

        ! Pass 2: contraction direction 2 (j <- l) -> fStarAux(i,j)
        do j = 0, slaveFace % NfLeft(2) ; do i = 0, slaveFace % NfLeft(1) ; do l = 0, slaveFace % Nf(2)
         dFv_dGradQElAux(1:NCONS,1:NCONS,1:NDIM,i,j) = dFv_dGradQElAux(1:NCONS,1:NCONS,1:NDIM,i,j) + MoutEta(j,l) * Flux_tmp(1:NCONS,1:NCONS, 1:NDIM,i,l) 
        end do ; end do ; end do

      select case (whichElement)
      case (1)
         associate(dFv_dGradQEl => self % storage(whichElement) % dFv_dGradQEl)
            dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv)=dFv_dGradQEl(1:NCONS,1:NCONS,1:NDIM,:,:,whichderiv) + &
            dFv_dGradQElAux(1:NCONS, 1:NCONS, 1:NDIM,:,:)
        end associate 
      case (2)
         associate(dFv_dGradQEl => self % storage(whichElement) % dFv_dGradQEl)
            dFv_dGradQEl(1:NCONS,1:NCONS, 1:NDIM,:,:,whichderiv)=dFv_dGradQEl(1:NCONS,1:NCONS, 1:NDIM,:,:,whichderiv) + &
            (-dFv_dGradQElAux(1:NCONS, 1:NCONS, 1:NDIM,:,:))
         end associate 
      end select 

   end subroutine Face_ProjectMortarGradJacobianToElements

#endif
#if defined(ACOUSTIC)
   subroutine Face_AdaptBaseSolutionToFace(self, nEqn, Nelx, Nely, Qe, side)
      use MappedGeometryClass
      implicit none
      class(Face),   intent(inout)              :: self
      integer,       intent(in)                 :: nEqn
      integer,       intent(in)                 :: Nelx, Nely
      real(kind=RP), intent(in)                 :: Qe(1:nEqn, 0:Nelx, 0:Nely)
      integer,       intent(in)                 :: side
!
!     ---------------
!     Local variables
!     ---------------
!
      integer       :: i, j, k, l, m, ii, jj
      real(kind=RP) :: Qe_rot(1:nEqn, 0:self % NfRight(1), 0:self % NfRight(2))
      select case (side)
      case(1)
         associate(Qf => self % storage(1) % Qbase)
         select case ( self % projectionType(1) )
         case (0)
            Qf = Qe

         case (1)
            Qf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfLeft(1)   ; do i = 0, self % Nf(1)
               Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfLeft(1), self % Nf(1)) % T(i,l) * Qe(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            Qf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) * Qe(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            Qf = 0.0_RP
            do l = 0, self % NfLeft(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfLeft(1) ; do i = 0, self % Nf(1)
                  Qf(:,i,j) = Qf(:,i,j) +   Tset(self % NfLeft(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfLeft(2), self % Nf(2)) % T(j,l) &
                                            * Qe(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select
         end associate
      case(2)
         associate( Qf => self % storage(2) % Qbase )
         do j = 0, self % NfRight(2)   ; do i = 0, self % NfRight(1)
            call leftIndexes2Right(i,j,self % NfRight(1), self % NfRight(2), self % rotation, ii, jj)
            Qe_rot(:,i,j) = Qe(:,ii,jj) 
         end do                        ; end do

         select case ( self % projectionType(2) )
         case (0)
            Qf = Qe_rot
         case (1)
            Qf = 0.0_RP
            do j = 0, self % Nf(2)  ; do l = 0, self % NfRight(1)   ; do i = 0, self % Nf(1)
               Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfRight(1), self % Nf(1)) % T(i,l) * Qe_rot(:,l,j)
            end do                  ; end do                   ; end do
            
         case (2)
            Qf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   ; do i = 0, self % Nf(1)
               Qf(:,i,j) = Qf(:,i,j) + Tset(self % NfRight(2), self % Nf(2)) % T(j,l) * Qe_rot(:,i,l)
            end do                  ; end do                   ; end do
   
         case (3)
            Qf = 0.0_RP
            do l = 0, self % NfRight(2)  ; do j = 0, self % Nf(2)   
               do m = 0, self % NfRight(1) ; do i = 0, self % Nf(1)
                  Qf(:,i,j) = Qf(:,i,j) +   Tset(self % NfRight(1), self % Nf(1)) % T(i,m) &
                                            * Tset(self % NfRight(2), self % Nf(2)) % T(j,l) &
                                            * Qe_rot(:,m,l)
               end do                 ; end do
            end do                  ; end do
         end select
         end associate
      end select

   end subroutine Face_AdaptBaseSolutionToFace

   !
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine Face_AdaptBaseSolutionToMortarFace(self, nEqn, Nelx, Nely, Qe, side, fma)
      use MappedGeometryClass
      implicit none
      class(Face),   intent(inout)              :: self
      integer,       intent(in)                 :: nEqn
      integer,       intent(in)                 :: Nelx, Nely
      real(kind=RP), intent(in)                 :: Qe(1:nEqn, 0:Nelx, 0:Nely)
      integer,       intent(in)                 :: side
      type(Face), intent(inout)                 :: fma

      integer       :: i, j, l, m

      real(kind=RP) :: MIntXi(0:fma%Nf(1), 0:fma%NfLeft(1))
      real(kind=RP) :: MIntEta(0:fma%Nf(2), 0:fma%NfLeft(2))

      real(kind=RP) :: tmp(1:nEqn, 0:fma%Nf(1), 0:fma%NfLeft(2))

      if (fma % MortarType .ne. MORTAR_SMALL4) then
         error stop 'MortarType SMALL4 reached in subroutine AdaptBaseSolutionToMortarFace expecting non-SMALL4: check calling logic'
      end if

      call GetMortarMInt(fma, MIntXi, MIntEta)

      associate(Qf => fma % storage(1) % Qbase)
         Qf = 0.0_RP

         ! Pass 1: contraction direction 1 (xi)  -> tmp(i,l)
         tmp = 0.0_RP
         do l = 0, fma % NfLeft(2) ; do i = 0, fma % Nf(1) ; do m = 0, fma % NfLeft(1)
            tmp(:,i,l) = tmp(:,i,l) + MIntXi(i,m) * Qe(:,m,l)
         end do ; end do ; end do

         ! Pass 2: contraction direction 2 (eta) -> Qf(i,j)
         do j = 0, fma % Nf(2) ; do i = 0, fma % Nf(1) ; do l = 0, fma % NfLeft(2)
            Qf(:,i,j) = Qf(:,i,j) + MIntEta(j,l) * tmp(:,i,l)
         end do ; end do ; end do

      end associate

   end subroutine Face_AdaptBaseSolutionToMortarFace

   subroutine Face_Interpolatebig2smallacoustic(self, nEqn, fma)
      use MappedGeometryClass
      implicit none
      CLASS(Face),   intent(inout)              :: self
      integer,       intent(in)                 :: nEqn
      type(Face), intent(inout)                 ::fma

#ifdef _HAS_MPI_

      integer       :: i, j, l, m
      real(kind=RP) :: Qfm(1:nEqn, 0:fma%Nf(1), 0:fma%Nf(2))

      real(kind=RP) :: MIntXi(0:fma%Nf(1), 0:fma%NfLeft(1))
      real(kind=RP) :: MIntEta(0:fma%Nf(2), 0:fma%NfLeft(2))

      real(kind=RP) :: tmp(1:nEqn, 0:fma%Nf(1), 0:fma%NfLeft(2))

      call GetMortarMInt(fma, MIntXi, MIntEta)

      associate(Qf => fma % storage(1) % Qbase)
         Qfm=Qf
         Qf=0.0_RP

         tmp = 0.0_RP
         do l = 0, fma % NfLeft(2) ; do i = 0, fma % Nf(1) ; do m = 0, fma % NfLeft(1)
            tmp(:,i,l) = tmp(:,i,l) + MIntXi(i,m) * Qfm(:,m,l)
         end do ; end do ; end do

         do j = 0, fma % Nf(2) ; do i = 0, fma % Nf(1) ; do l = 0, fma % NfLeft(2)
            Qf(:,i,j) = Qf(:,i,j) + MIntEta(j,l) * tmp(:,i,l)
         end do ; end do ; end do

      end associate 

#endif
   end subroutine Face_Interpolatebig2smallacoustic
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE InterpolateToBoundary( u, v, Nx, Ny, Nz, which_dim , bValue , NEQ)
      use SMConstants
      use PhysicsStorage
      IMPLICIT NONE
      integer                      , INTENT(IN)    :: Nx, Ny, Nz
      real(kind=RP)                , intent(in)    :: u(1:NEQ,0:Nx,0:Ny,0:Nz) , v(0:)
      integer                      , intent(in)    :: which_dim
      REAL(KIND=RP)                , INTENT(INOUT) :: bValue(1:,0:,0:)
      integer                      , intent(in)    :: NEQ
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                                    :: i , j , k , eq

      select case (which_dim)
      case (IX)

         do j = 0 , Nz ; do i = 0 , Ny ; do k = 0 , Nx
            bValue(:,i,j) = bValue(:,i,j) + u(:,k,i,j) * v(k)
         end do        ; end do        ; end do

      case (IY)

         do j = 0 , Nz ; do k = 0 , Ny ; do i = 0 , Nx
            bValue(:,i,j) = bValue(:,i,j) + u(:,i,k,j) * v(k)
         end do        ; end do        ; end do

      case (IZ)

         do k = 0 , Nz ; do j = 0 , Ny ; do i = 0 , Nx
            bValue(:,i,j) = bValue(:,i,j) + u(:,i,j,k) * v(k)
         end do        ; end do        ; end do

      end select
      
      end SUBROUTINE InterpolateToBoundary
!
!////////////////////////////////////////////////////////////////////////
!
   elemental subroutine Face_Assign(to,from)
      implicit none
      class(Face), intent(inout) :: to
      class(Face), intent(in)    :: from
      
      
      to % flat = from % flat
      to % ID = from % ID
      to % FaceType = from % FaceType
      to % zone = from % zone
      to % rotation = from % rotation
      to % NelLeft = from % NelLeft
      to % NelRight = from % NelRight
      to % NfLeft = from % NfLeft
      to % NfRight = from % NfRight
      to % Nf = from % Nf
      to % nodeIDs = from % nodeIDs
      to % elementIDs = from % elementIDs
      to % elementSide = from % elementSide
      to % projectionType = from % projectionType
      to % boundaryName = from % boundaryName
      to % geom = from % geom
      to % storage = from % storage
      to % MortarType = from % MortarType
      to % Mortarpos = from % Mortarpos

   end subroutine Face_Assign

!
!////////////////////////////////////////////////////////////////////////
!
!  Get the projection matrices for sliding mortars
!
   subroutine GetSlidingMInt(fma, MInt)
      implicit none
      type(Face),    intent(in)  :: fma
      real(kind=RP), intent(out) :: MInt(0:fma%Nf(1), 0:fma%NfLeft(1), 1:2)

      if (fma % Mortarpos == 0) then
         MInt(:,:,1) = TsetM(fma % NfLeft(1), fma % Nf(1), 1, 1) % T
         MInt(:,:,2) = TsetM(fma % NfLeft(1), fma % Nf(1), 2, 1) % T
      else if (fma % Mortarpos == 1) then
         MInt(:,:,1) = TsetM(fma % NfLeft(1), fma % Nf(1), 3, 1) % T
         MInt(:,:,2) = TsetM(fma % NfLeft(1), fma % Nf(1), 4, 1) % T
      end if

   end subroutine GetSlidingMInt
!
!////////////////////////////////////////////////////////////////////////
!
   subroutine GetSlidingMout(fma, Mout)
      implicit none
      type(Face),    intent(in)  :: fma
      real(kind=RP), intent(out) :: Mout(0:fma%NfLeft(1), 0:fma%Nf(1), 1:2)

      if (fma % Mortarpos == 0) then
         Mout(:,:,1) = TsetM(fma % Nf(1), fma % NfLeft(1), 1, 2) % T
         Mout(:,:,2) = TsetM(fma % Nf(1), fma % NfLeft(1), 2, 2) % T
      else if (fma % Mortarpos == 1) then
         Mout(:,:,1) = TsetM(fma % Nf(1), fma % NfLeft(1), 3, 2) % T
         Mout(:,:,2) = TsetM(fma % Nf(1), fma % NfLeft(1), 4, 2) % T
      end if

   end subroutine GetSlidingMout
!
!////////////////////////////////////////////////////////////////////////
!
!  Get the projection matrices for hp (4:1) adaptations
!
   subroutine GetMortarMInt(f, MIntXi, MIntEta)
      implicit none
      type(Face),    intent(in)  :: f
      real(kind=RP), intent(out) :: MIntXi(0:f%Nf(1), 0:f%NfLeft(1))
      real(kind=RP), intent(out) :: MIntEta(0:f%Nf(2), 0:f%NfLeft(2))

      select case (f % Mortarpos)
      case (1)
         MIntXi  = TsetM(f % NfLeft(1), f % Nf(1), 1, 1) % T
         MIntEta = TsetM(f % NfLeft(2), f % Nf(2), 1, 1) % T
      case (2)
         MIntXi  = TsetM(f % NfLeft(1), f % Nf(1), 2, 1) % T
         MIntEta = TsetM(f % NfLeft(2), f % Nf(2), 1, 1) % T
      case (3)
         MIntXi  = TsetM(f % NfLeft(1), f % Nf(1), 1, 1) % T
         MIntEta = TsetM(f % NfLeft(2), f % Nf(2), 2, 1) % T
      case (4)
         MIntXi  = TsetM(f % NfLeft(1), f % Nf(1), 2, 1) % T
         MIntEta = TsetM(f % NfLeft(2), f % Nf(2), 2, 1) % T
      end select

   end subroutine GetMortarMInt
!
!////////////////////////////////////////////////////////////////////////
!
   subroutine GetMortarMout(f, MoutXi, MoutEta)
      implicit none
      type(Face),    intent(in)  :: f
      real(kind=RP), intent(out) :: MoutXi(0:f%NfLeft(1), 0:f%Nf(1))
      real(kind=RP), intent(out) :: MoutEta(0:f%NfLeft(2), 0:f%Nf(2))

      select case (f % Mortarpos)
      case (1)
         MoutXi  = TsetM(f % Nf(1), f % NfLeft(1), 1, 2) % T
         MoutEta = TsetM(f % Nf(2), f % NfLeft(2), 1, 2) % T
      case (2)
         MoutXi  = TsetM(f % Nf(1), f % NfLeft(1), 2, 2) % T
         MoutEta = TsetM(f % Nf(2), f % NfLeft(2), 1, 2) % T
      case (3)
         MoutXi  = TsetM(f % Nf(1), f % NfLeft(1), 1, 2) % T
         MoutEta = TsetM(f % Nf(2), f % NfLeft(2), 2, 2) % T
      case (4)
         MoutXi  = TsetM(f % Nf(1), f % NfLeft(1), 2, 2) % T
         MoutEta = TsetM(f % Nf(2), f % NfLeft(2), 2, 2) % T
      end select

   end subroutine GetMortarMout

end Module FaceClass
   
