!
!//////////////////////////////////////////////////////
!
!      Implements Algorithms:
!         Algorithm: 124: ElementClass (QuadElementClass)
!
!       The Quad Element class, Alg. 124. See Sec. 8.2.1.2. self has
!       been modified to add the association of a boundary name to an element
!       edge so that different boundary conditions can be applied to different
!       elements. The names of the boundaries (not necessarily the names of the
!       *boundary conditions* to be applied) are of length BC_STRING_LENGTH.
!       One will associate boundary conditions to boundaries in the routine
!       "ExternalState".
!
!       Modified 2D Code to move solution into element class. 5/14/15, 5:36 PM
!
!////////////////////////////////////////////////////////////////////////
!
   Module ElementClass
      USE SMConstants
      USE PolynomialInterpAndDerivsModule
      USE GaussQuadrature
      USE TransfiniteMapClass
      USE MappedGeometryClass
      USE MeshTypes
      USE ElementConnectivityDefinitions
      USE ConnectivityClass
      use StorageClass
      USE NodalStorageClass               , only: NodalStorage, NodalStorage_t
      use PhysicsStorage
      use Physics
      use VariableConversion
      use FacePatchClass                  , only: FacePatch
      IMPLICIT NONE

      private
      public   Element
      public   PrintElement, SetElementBoundaryNames, SurfInfo_t

!
!     -----------------------------------------------------------------------------------------
!     Type containing the information of the surfaces of an element, as read from the mesh file
!     -> This is used only for constructung the mesh by the readers
!     -----------------------------------------------------------------------------------------
      type SurfInfo_t
         ! Variables to specify that the element is a hex8
         logical         :: IsHex8 = .FALSE.
         real(kind=RP)   :: corners(NDIM,NODES_PER_ELEMENT)
         ! Variables for general elements
         type(FacePatch) :: facePatches(6)
         contains
         procedure :: destruct  => SurfInfo_Destruct
      end type SurfInfo_t
!
!     ---------------------
!     Main hex-element type
!     ---------------------
      TYPE Element
         integer                         :: MortarFaces(6)     !0=Conforming; 1=Master; 2= Slave, 3=Master_sliding; 4= Slave_sliding
         real(kind=RP)                   :: Psvv
         logical                         :: hasSharedFaces
         integer                         :: dir2D
         integer                         :: globDir(3)        ! If the global coordinate (GLOBAL) is aligned with the local coordinate (REFERENCE): globDir(GLOBAL) = REFERENCE
         integer                         :: eID               ! ID of this element
         integer                         :: globID            ! globalID of the element
         integer                         :: offsetIO          ! Offset from the first element for IO
         INTEGER                         :: nodeIDs(8)
         integer                         :: faceIDs(6)
         integer                         :: faceSide(6)
         integer                         :: MLevel            ! RK Level
         integer                         :: MLevelwN          ! RK Level with Neighbour
         real(kind=RP)                   :: ML_CFL            ! CFL storage for Multi Level RK	
         real(kind=RP)                   :: ML_error_ratio    ! Ratio between temporal and spatial error relative to the global ratio
         INTEGER, DIMENSION(3)           :: Nxyz              ! Polynomial orders in every direction (Nx,Ny,Nz)
         real(kind=RP)                   :: hn                ! Ratio of size and polynomial order
         TYPE(MappedGeometry)            :: geom
         CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryName(6)
         INTEGER                         :: NumberOfConnections(6)
         TYPE(Connectivity)              :: Connection(6)
         type(ElementStorage_t), pointer :: storage
         type(SurfInfo_t)                :: SurfInfo          ! Information about the geometry of the neighboring faces, as in the mesh file
         type(TransfiniteHexMap)         :: hexMap            ! High-order mapper
         logical, dimension(:,:,:), allocatable :: isInsideBody, isForcingPoint ! Immersed boundaty term -> if InsideBody(i,j,k) = true, the point(i,j,k) is inside the body (IB)	
         integer, dimension(:,:,:), allocatable :: STL !STL file the DoFbelongs to if isInsideBody = .true. (IB)
         integer                                :: IP_index
         logical                                :: sliding
         logical                                :: sliding_newnodes 
         contains
            procedure   :: Construct               => HexElement_Construct
            procedure   :: Destruct                => HexElement_Destruct
            procedure   :: ConstructGeometry       => HexElement_ConstructGeometry
            procedure   :: FindPointInLinElement   => HexElement_FindPointInLinearizedElement
            procedure   :: FindPointWithCoords     => HexElement_FindPointWithCoords
            procedure   :: EvaluateSolutionAtPoint => HexElement_EvaluateSolutionAtPoint
            procedure   :: ProlongSolutionToFaces  => HexElement_ProlongSolutionToFaces
            procedure   :: ProlongGradientsToFaces => HexElement_ProlongGradientsToFaces
            procedure   :: ProlongAviscFluxToFaces => HexElement_ProlongAviscFluxToFaces
            procedure   :: ComputeLocalGradient    => HexElement_ComputeLocalGradient
            procedure   :: pAdapt                  => HexElement_pAdapt
            procedure   :: copy                    => HexElement_Assign
            procedure   :: EvaluateGradientAtPoint => HexElement_EvaluateGradientAtPoint
            procedure   :: ConstructIBM            => HexElement_ConstructIBM
#if defined(ACOUSTIC)
               procedure   :: ProlongBaseSolutionToFaces  => HexElement_ProlongBaseSolutionToFaces
#endif
            generic     :: assignment(=)           => copy
      END TYPE Element

      CONTAINS
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE HexElement_Construct( self, Nx, Ny, Nz, nodeIDs, eID, globID)
         IMPLICIT NONE
         class(Element)      :: self
         integer, intent(in) :: Nx, Ny, Nz  !<  Polynomial orders
         integer, intent(in) :: nodeIDs(8)
         integer, intent(in) :: eID, globID

         self % eID                 = eID
         self % dir2D               = 0
         self % globDir             = 0
         self % globID              = globID
         self % nodeIDs             = nodeIDs
         self % Nxyz(1)             = Nx
         self % Nxyz(2)             = Ny
         self % Nxyz(3)             = Nz
         self % hn                  = 0.0_RP
         self % boundaryName        = emptyBCName
         self % hasSharedFaces      = .false.
         self % NumberOfConnections = 0
         self % MLevel              = 1
         self % ML_error_ratio      = 1.0_RP
         self % MortarFaces         = 0 
         self % sliding             = .false.
         self % sliding_newnodes    = .false.
!
!        ----------------------------------------
!        Solution Storage is allocated separately
!        ----------------------------------------
!
      END SUBROUTINE HexElement_Construct
!
!////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------------
!     Constructs the mapped geometry of the element (metric terms)
!     ------------------------------------------------------------
      subroutine HexElement_ConstructGeometry( self, hexMap)
         implicit none
         !--------------------------------------
         class(Element)         , intent(inout) :: self
         TYPE(TransfiniteHexMap), intent(in)    :: hexMap
         !--------------------------------------

         self % hexMap = hexMap
         CALL self % geom % Construct( NodalStorage(Self % Nxyz(1)), NodalStorage(Self % Nxyz(2)), NodalStorage(Self % Nxyz(3)), hexMap )

      end subroutine HexElement_ConstructGeometry

!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetElementBoundaryNames( self, names )
         use Utilities, only: toLower
         IMPLICIT NONE
         TYPE(Element)                   :: self
         CHARACTER(LEN=BC_STRING_LENGTH) :: names(6)
         INTEGER                         :: j

         DO j = 1, 6
            CALL toLower(names(j))
            self % boundaryName(j) = names(j)
         END DO
      END SUBROUTINE SetElementBoundaryNames
!
!////////////////////////////////////////////////////////////////////////
!
      elemental SUBROUTINE HexElement_Destruct( self )
         IMPLICIT NONE
         class(Element), intent(inout) :: self

         CALL self % geom % Destruct
         call self % Connection % destruct
         call self % hexMap % destruct

         call self % SurfInfo % destruct
         
         if( allocated(self% isInsideBody) )   deallocate(self% isInsideBody)
         if( allocated(self% isForcingPoint) ) deallocate(self% isForcingPoint)
         if( allocated(self% STL) )            deallocate(self% STL)

      END SUBROUTINE HexElement_Destruct
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintElement( self, id )
         IMPLICIT NONE
         TYPE(Element) :: self
         INTEGER      :: id
         PRINT *, id, self % nodeIDs
         PRINT *, "   ",self % boundaryName
      END SUBROUTINE PrintElement
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SaveSolutionStorageToUnit( self, fUnit )
         IMPLICIT NONE
!
!        -----------------------
!        Save for a restart file
!        -----------------------
!
         TYPE(Element) :: self
         INTEGER       :: fUnit

         WRITE(funit) self % storage % Q

      END SUBROUTINE SaveSolutionStorageToUnit
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LoadSolutionFromUnit( self, fUnit )
         IMPLICIT NONE
!
!        -----------------------
!        Save for a restart file
!        -----------------------
!
         TYPE(Element) :: self
         INTEGER       :: fUnit

         READ(funit) self % storage % Q

      END SUBROUTINE LoadSolutionFromUnit
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexElement_ProlongSolutionToFaces(self, nEqn, fFR, fBK, fBOT, fR, fT, fL, computeQdot, faces)
         use FaceClass
         implicit none
         class(Element),   intent(in)  :: self
         integer,          intent(in)  :: nEqn
         class(Face),      intent(inout) :: fFR, fBK, fBOT, fR, fT, fL
         logical,optional, intent(in)  :: computeQdot
         type(Face),optional, intent(inout) :: faces(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k, l, N(3), m
         real(kind=RP), dimension(1:nEqn, 0:self % Nxyz(1), 0:self % Nxyz(3)) :: QFR, QBK, QdotFR, QdotBK
         real(kind=RP), dimension(1:nEqn, 0:self % Nxyz(1), 0:self % Nxyz(2)) :: QBOT, QT, QdotBOT, QdotT
         real(kind=RP), dimension(1:nEqn, 0:self % Nxyz(2), 0:self % Nxyz(3)) :: QL, QR, QdotL, QdotR
         type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta
         logical :: prolongQdot

         N = self % Nxyz
         spAxi   => NodalStorage(N(1))
         spAeta  => NodalStorage(N(2))
         spAzeta => NodalStorage(N(3))

         if (present(computeQdot)) then
               prolongQdot = computeQdot
         else
               prolongQdot = .FALSE.
         end if
!
!        *************************
!        Prolong solution to faces
!        *************************
!
         QL   = 0.0_RP     ; QR   = 0.0_RP
         QFR  = 0.0_RP     ; QBK  = 0.0_RP
         QBOT = 0.0_RP     ; QT   = 0.0_RP

         QdotL   = 0.0_RP     ; QdotR   = 0.0_RP
         QdotFR  = 0.0_RP     ; QdotBK  = 0.0_RP
         QdotBOT = 0.0_RP     ; QdotT   = 0.0_RP

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            QL  (:,j,k)= QL  (:,j,k)+ self % storage % Q(:,i,j,k)* spAxi % v  (i,LEFT  )
            QR  (:,j,k)= QR  (:,j,k)+ self % storage % Q(:,i,j,k)* spAxi % v  (i,RIGHT )
            QFR (:,i,k)= QFR (:,i,k)+ self % storage % Q(:,i,j,k)* spAeta % v (j,FRONT )
            QBK (:,i,k)= QBK (:,i,k)+ self % storage % Q(:,i,j,k)* spAeta % v (j,BACK  )
            QBOT(:,i,j)= QBOT(:,i,j)+ self % storage % Q(:,i,j,k)* spAzeta % v(k,BOTTOM)
            QT  (:,i,j)= QT  (:,i,j)+ self % storage % Q(:,i,j,k)* spAzeta % v(k,TOP   )
            if (prolongQdot) then
                  QdotL  (:,j,k)= QdotL  (:,j,k)+ self % storage % Qdot(:,i,j,k)* spAxi % v  (i,LEFT  )
                  QdotR  (:,j,k)= QdotR  (:,j,k)+ self % storage % Qdot(:,i,j,k)* spAxi % v  (i,RIGHT )
                  QdotFR (:,i,k)= QdotFR (:,i,k)+ self % storage % Qdot(:,i,j,k)* spAeta % v (j,FRONT )
                  QdotBK (:,i,k)= QdotBK (:,i,k)+ self % storage % Qdot(:,i,j,k)* spAeta % v (j,BACK  )
                  QdotBOT(:,i,j)= QdotBOT(:,i,j)+ self % storage % Qdot(:,i,j,k)* spAzeta % v(k,BOTTOM)
                  QdotT  (:,i,j)= QdotT  (:,i,j)+ self % storage % Qdot(:,i,j,k)* spAzeta % v(k,TOP   )
            end if
         end do                   ; end do                   ; end do
         nullify (spAxi, spAeta, spAzeta)
         if (present(faces)) then 
            if (fL % MortarType == MORTAR_BIG) then 
               if (fL % n_mpi_mortar .ne. 0) then 
                  call fL   % AdaptSolutionToFace(nEqn, N(2), N(3), QL   , self % faceSide(ELEFT  ), QdotL, computeQdot)
               end if 
               do m=1,4
                  if (fL % Mortar(m) .ne. 0) then 
                     call fL % AdaptSolutionToMortarFace(nEqn, N(2), N(3), QL   , self % faceSide(ELEFT  ), QdotL, computeQdot,&
                     faces(fL % Mortar(m)))
                  end if 
               end do 
            elseif (fL % MortarType == MORTAR_SLIDING) then 
               do m=1,size(fL%Mortar)
                  call fL % AdaptSolutionToMortarFace(nEqn, N(2), N(3), QL   , self % faceSide(ELEFT  ), QdotL, computeQdot,&
                     faces(fL % Mortar(m)), .true.)
               end do 
            else 
               call fL   % AdaptSolutionToFace(nEqn, N(2), N(3), QL   , self % faceSide(ELEFT  ), QdotL, computeQdot)
            end if 
            if (fR % MortarType == MORTAR_BIG) then 
               if (fR % n_mpi_mortar .ne. 0) then 
                  call fR   % AdaptSolutionToFace(nEqn, N(2), N(3), QR   , self % faceSide(ERIGHT ), QdotR, computeQdot)
               end if 
               do m=1,4
                  if (fR % Mortar(m) .ne. 0) then 
                     call fR % AdaptSolutionToMortarFace(nEqn, N(2), N(3), QR   , self % faceSide(ERIGHT  ), QdotR, computeQdot,&
                     faces(fR % Mortar(m)))
                  end if 
               end do 
            elseif (fR % MortarType == MORTAR_SLIDING) then 
                  do m=1,size(fR%Mortar)
                     call fR % AdaptSolutionToMortarFace(nEqn, N(2), N(3), QR   , self % faceSide(ERIGHT  ), QdotR, computeQdot,&
                        faces(fR % Mortar(m)), .true.)
                  end do 
            else 
               call fR   % AdaptSolutionToFace(nEqn, N(2), N(3), QR   , self % faceSide(ERIGHT ), QdotR, computeQdot)
            end if 
            if (fFR % MortarType == MORTAR_BIG) then 
               if (fFR % n_mpi_mortar .ne. 0) then 
                  call fFR   % AdaptSolutionToFace(nEqn, N(1), N(3), QFR   , self % faceSide(EFRONT ), QdotFR, computeQdot)
               end if 
               do m=1,4
                  if (fFR % Mortar(m) .ne. 0) then 
                     call fFR % AdaptSolutionToMortarFace(nEqn, N(1), N(3), QFR   , self % faceSide(EFRONT  ), QdotFR, computeQdot,&
                     faces(fFR % Mortar(m)))
                  end if 
               end do 
               elseif (fFR % MortarType == MORTAR_SLIDING) then 
                  do m=1,size(fFR%Mortar)
                     call fFR % AdaptSolutionToMortarFace(nEqn, N(1), N(3), QFR   , self % faceSide(EFRONT  ), QdotFR, computeQdot,&
                        faces(fFR % Mortar(m)), .true.)
                  end do 
            else 
               call fFR   % AdaptSolutionToFace(nEqn, N(1), N(3), QFR   , self % faceSide(EFRONT ), QdotFR, computeQdot)
            end if 
            if (fBK % MortarType == MORTAR_BIG) then 
               if (fBK % n_mpi_mortar .ne. 0) then 
                  call fBK   % AdaptSolutionToFace(nEqn, N(1), N(3), QBK   , self % faceSide(EBACK  ), QdotBK, computeQdot)
               end if 
               do m=1,4
                  if (fBK% Mortar(m) .ne. 0) then 
                     call fBK % AdaptSolutionToMortarFace(nEqn, N(1), N(3), QBK   , self % faceSide(EBACK  ), QdotBK, computeQdot,&
                     faces(fBK % Mortar(m)))
                  end if
               end do 
               elseif (fBK % MortarType == MORTAR_SLIDING) then 
                  do m=1,size(fBK%Mortar)
                     call fBK % AdaptSolutionToMortarFace(nEqn, N(1), N(3), QBK   , self % faceSide(EBACK  ), QdotBK, computeQdot,&
                        faces(fBK % Mortar(m)), .true.)
                  end do 
            else 
               call fBK   % AdaptSolutionToFace(nEqn, N(1), N(3), QBK   , self % faceSide(EBACK  ), QdotBK, computeQdot)
            end if 
            if (fBOT % MortarType == MORTAR_BIG) then 
               if (fBOT % n_mpi_mortar .ne. 0) then 
                  call fBOT   % AdaptSolutionToFace(nEqn, N(1), N(2), QBOT   , self % faceSide(EBOTTOM), QdotBOT, computeQdot)
               end if 
               do m=1,4
                  if (fBOT % Mortar(m) .ne. 0) then 
                     call fBOT % AdaptSolutionToMortarFace(nEqn, N(1), N(2), QBOT   , self % faceSide(EBOTTOM  ), QdotBOT, computeQdot,&
                     faces(fBOT % Mortar(m)))
                  end if 
               end do 
               elseif (fBOT % MortarType == MORTAR_SLIDING) then 
                  do m=1,size(fBOT%Mortar)
                     call fBOT % AdaptSolutionToMortarFace(nEqn, N(1), N(2), QBOT   , self % faceSide(EBOTTOM  ), QdotBOT, computeQdot,&
                        faces(fBOT % Mortar(m)), .true.)
                  end do 
            else 
               call fBOT   % AdaptSolutionToFace(nEqn, N(1), N(2), QBOT   , self % faceSide(EBOTTOM), QdotBOT, computeQdot)
            end if 
            if (fT % MortarType == MORTAR_BIG) then 
               if (fT % n_mpi_mortar .ne. 0) then 
                  call fT   % AdaptSolutionToFace(nEqn, N(1), N(2), QT   , self % faceSide(ETOP   ), QdotT, computeQdot)
               end if 
               do m=1,4
                  if (fT% Mortar(m) .ne. 0) then 
                     call fT % AdaptSolutionToMortarFace(nEqn, N(1), N(2), QT   , self % faceSide(ETOP  ), QdotT, computeQdot,&
                     faces(fT % Mortar(m)))
                  end if 
               end do 
               elseif (fT % MortarType == MORTAR_SLIDING) then 
                  do m=1,size(fT%Mortar)
                     call fT % AdaptSolutionToMortarFace(nEqn, N(1), N(2), QT   , self % faceSide(ETOP  ), QdotT, computeQdot,&
                        faces(fT % Mortar(m)), .true.)
                  end do 
            else 
               call fT   % AdaptSolutionToFace(nEqn, N(1), N(2), QT   , self % faceSide(ETOP   ), QdotT, computeQdot)
            end if 
         else
            call fL   % AdaptSolutionToFace(nEqn, N(2), N(3), QL   , self % faceSide(ELEFT  ), QdotL, computeQdot)
            call fR   % AdaptSolutionToFace(nEqn, N(2), N(3), QR   , self % faceSide(ERIGHT ), QdotR, computeQdot)
            call fFR  % AdaptSolutionToFace(nEqn, N(1), N(3), QFR  , self % faceSide(EFRONT ), QdotFR, computeQdot)
            call fBK  % AdaptSolutionToFace(nEqn, N(1), N(3), QBK  , self % faceSide(EBACK  ), QdotBK, computeQdot)
            call fBOT % AdaptSolutionToFace(nEqn, N(1), N(2), QBOT , self % faceSide(EBOTTOM), QdotBOT, computeQdot)
            call fT   % AdaptSolutionToFace(nEqn, N(1), N(2), QT   , self % faceSide(ETOP   ), QdotT, computeQdot)
         end if

      end subroutine HexElement_ProlongSolutionToFaces

      subroutine HexElement_ProlongGradientsToFaces(self, nGradEqn, fFR, fBK, fBOT, fR, fT, fL, faces)
         use FaceClass
         implicit none
         class(Element),   intent(in)  :: self
         integer,          intent(in)  :: nGradEqn
         class(Face),      intent(inout) :: fFR, fBK, fBOT, fR, fT, fL
         type(Face),optional, intent(inout) :: faces(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k, l, m, N(3)
         real(kind=RP), dimension(nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(3)) :: UxFR, UyFR, UzFR
         real(kind=RP), dimension(nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(3)) :: UxBK, UyBK, UzBK
         real(kind=RP), dimension(nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(2)) :: UxBT, UyBT, UzBT
         real(kind=RP), dimension(nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(2)) :: UxT, UyT, UzT
         real(kind=RP), dimension(nGradEqn, 0:self % Nxyz(2), 0:self % Nxyz(3)) :: UxL, UyL, UzL
         real(kind=RP), dimension(nGradEqn, 0:self % Nxyz(2), 0:self % Nxyz(3)) :: UxR, UyR, UzR
         type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta

         N = self % Nxyz
         spAxi   => NodalStorage(N(1))
         spAeta  => NodalStorage(N(2))
         spAzeta => NodalStorage(N(3))
!
!        *************************
!        Prolong solution to faces
!        *************************
!
         UxL  = 0.0_RP ; UyL  = 0.0_RP ; UzL  = 0.0_RP
         UxR  = 0.0_RP ; UyR  = 0.0_RP ; UzR  = 0.0_RP
         UxFR = 0.0_RP ; UyFR = 0.0_RP ; UzFR = 0.0_RP
         UxBK = 0.0_RP ; UyBK = 0.0_RP ; UzBK = 0.0_RP
         UxBT = 0.0_RP ; UyBT = 0.0_RP ; UzBT = 0.0_RP
         UxT  = 0.0_RP ; UyT  = 0.0_RP ; UzT  = 0.0_RP

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            UxL (:,j,k) = UxL (:,j,k) + self % storage % U_x(:,i,j,k)* spAxi   % v (i,LEFT  )
            UxR (:,j,k) = UxR (:,j,k) + self % storage % U_x(:,i,j,k)* spAxi   % v (i,RIGHT )
            UxFR(:,i,k) = UxFR(:,i,k) + self % storage % U_x(:,i,j,k)* spAeta  % v (j,FRONT )
            UxBK(:,i,k) = UxBK(:,i,k) + self % storage % U_x(:,i,j,k)* spAeta  % v (j,BACK  )
            UxBT(:,i,j) = UxBT(:,i,j) + self % storage % U_x(:,i,j,k)* spAzeta % v (k,BOTTOM)
            UxT (:,i,j) = UxT (:,i,j) + self % storage % U_x(:,i,j,k)* spAzeta % v (k,TOP   )

            UyL (:,j,k) = UyL (:,j,k) + self % storage % U_y(:,i,j,k)* spAxi   % v (i,LEFT  )
            UyR (:,j,k) = UyR (:,j,k) + self % storage % U_y(:,i,j,k)* spAxi   % v (i,RIGHT )
            UyFR(:,i,k) = UyFR(:,i,k) + self % storage % U_y(:,i,j,k)* spAeta  % v (j,FRONT )
            UyBK(:,i,k) = UyBK(:,i,k) + self % storage % U_y(:,i,j,k)* spAeta  % v (j,BACK  )
            UyBT(:,i,j) = UyBT(:,i,j) + self % storage % U_y(:,i,j,k)* spAzeta % v (k,BOTTOM)
            UyT (:,i,j) = UyT (:,i,j) + self % storage % U_y(:,i,j,k)* spAzeta % v (k,TOP   )

            UzL (:,j,k) = UzL (:,j,k) + self % storage % U_z(:,i,j,k)* spAxi   % v (i,LEFT  )
            UzR (:,j,k) = UzR (:,j,k) + self % storage % U_z(:,i,j,k)* spAxi   % v (i,RIGHT )
            UzFR(:,i,k) = UzFR(:,i,k) + self % storage % U_z(:,i,j,k)* spAeta  % v (j,FRONT )
            UzBK(:,i,k) = UzBK(:,i,k) + self % storage % U_z(:,i,j,k)* spAeta  % v (j,BACK  )
            UzBT(:,i,j) = UzBT(:,i,j) + self % storage % U_z(:,i,j,k)* spAzeta % v (k,BOTTOM)
            UzT (:,i,j) = UzT (:,i,j) + self % storage % U_z(:,i,j,k)* spAzeta % v (k,TOP   )

         end do                   ; end do                   ; end do

         nullify (spAxi, spAeta, spAzeta)
         if (present(faces)) then 
            if (fL % MortarType == MORTAR_BIG) then 
               if (fL % n_mpi_mortar .ne. 0) then 
                  call fL   % AdaptGradientsToFace(nGradEqn, N(2), N(3), UxL , UyL , UzL , self % faceSide(ELEFT  ))
               end if 
               do m=1,4
                  if (fL % Mortar(m) .ne. 0) then 
                     call fL % AdaptGradientsToMortarFace(nGradEqn, N(2), N(3), UxL , UyL , UzL , self % faceSide(ELEFT), faces(fL % Mortar(m)))
                  end if 
               end do 
            elseif(fL % MortarType == MORTAR_SLIDING) then 
               do m=1,size(fL%Mortar)
                  call fL % AdaptGradientsToMortarFace(nGradEqn, N(2), N(3), UxL , UyL , UzL , self % faceSide(ELEFT), faces(fL % Mortar(m)), .true.)
               end do 
            else 
               call fL   % AdaptGradientsToFace(nGradEqn, N(2), N(3), UxL , UyL , UzL , self % faceSide(ELEFT  ))
            end if 

            if (fR % MortarType == MORTAR_BIG) then 
               if (fR % n_mpi_mortar .ne. 0) then 
                  call fR   % AdaptGradientsToFace(nGradEqn, N(2), N(3), UxR , UyR , UzR , self % faceSide(ERIGHT  ))
               end if 
               do m=1,4
                  if (fR % Mortar(m) .ne. 0) then 
                     call fR % AdaptGradientsToMortarFace(nGradEqn, N(2), N(3), UxR , UyR , UzR , self % faceSide(ERIGHT), faces(fR % Mortar(m)))
                  end if 
               end do 
            elseif(fR % MortarType == MORTAR_SLIDING) then 
               do m=1,size(fR%Mortar)
                  call fR % AdaptGradientsToMortarFace(nGradEqn, N(2), N(3), UxR , UyR , UzR , self % faceSide(ERIGHT), faces(fR % Mortar(m)), .true.)
               end do 
            else 
               call fR   % AdaptGradientsToFace(nGradEqn, N(2), N(3), UxR , UyR , UzR , self % faceSide(ERIGHT  ))
            end if 

            if (fFR % MortarType == MORTAR_BIG) then 
               if (fFR % n_mpi_mortar .ne. 0) then 
                  call fFR   % AdaptGradientsToFace(nGradEqn, N(1), N(3), UxFR , UyFR , UzFR , self % faceSide(EFRONT  ))
               end if 
               do m=1,4
                  if (fFR % Mortar(m) .ne. 0) then 
                     call fFR % AdaptGradientsToMortarFace(nGradEqn, N(1), N(3), UxFR , UyFR , UzFR , self % faceSide(EFRONT), faces(fFR % Mortar(m)))
                  end if 
               end do 
            elseif(fFR % MortarType == MORTAR_SLIDING) then 
               do m=1,size(fFR%Mortar)
                  call fFR % AdaptGradientsToMortarFace(nGradEqn, N(1), N(3), UxFR , UyFR , UzFR , self % faceSide(EFRONT), faces(fFR % Mortar(m)), .true.)
               end do 
            else 
               call fFR   % AdaptGradientsToFace(nGradEqn, N(1), N(3), UxFR , UyFR , UzFR , self % faceSide(EFRONT  ))
            end if 

            if (fBK % MortarType == MORTAR_BIG) then 
               if (fBK % n_mpi_mortar .ne. 0) then 
                  call fBK   % AdaptGradientsToFace(nGradEqn, N(1), N(3), UxBK , UyBK , UzBK , self % faceSide(EBACK  ))
               end if 
               do m=1,4
                  if (fBK % Mortar(m) .ne. 0) then 
                     call fBK % AdaptGradientsToMortarFace(nGradEqn, N(1), N(3), UxBK , UyBK , UzBK , self % faceSide(EBACK), faces(fBK % Mortar(m)))
                  end if 
               end do 
            elseif(fBK % MortarType == MORTAR_SLIDING) then 
               do m=1,size(fBK%Mortar)
                  call fBK % AdaptGradientsToMortarFace(nGradEqn, N(1), N(3), UxBK , UyBK , UzBK , self % faceSide(EBACK), faces(fBK % Mortar(m)), .true.)
               end do 
            else 
               call fBK   % AdaptGradientsToFace(nGradEqn, N(1), N(3), UxBK , UyBK , UzBK , self % faceSide(EBACK  ))
            end if 

            if (fBOT % MortarType == MORTAR_BIG) then 
               if (fBOT % n_mpi_mortar .ne. 0) then 
                  call fBOT   % AdaptGradientsToFace(nGradEqn, N(1), N(2), UxBT , UyBT , UzBT , self % faceSide(EBOTTOM  ))
               end if 
               do m=1,4
                  if (fBOT % Mortar(m) .ne. 0) then 
                     call fBOT % AdaptGradientsToMortarFace(nGradEqn, N(1), N(2), UxBT , UyBT , UzBT , self % faceSide(EBOTTOM), faces(fBOT % Mortar(m)))
                  end if 
               end do 
            elseif(fBOT % MortarType == MORTAR_SLIDING) then 
               do m=1,size(fBOT%Mortar)
                  call fBOT % AdaptGradientsToMortarFace(nGradEqn, N(1), N(2), UxBT , UyBT , UzBT , self % faceSide(EBOTTOM), faces(fBOT % Mortar(m)), .true.)
               end do 
            else 
               call fBOT   % AdaptGradientsToFace(nGradEqn, N(1), N(2), UxBT , UyBT , UzBT , self % faceSide(EBOTTOM  ))
            end if 

            if (fT % MortarType == MORTAR_BIG) then 
               if (fT % n_mpi_mortar .ne. 0) then 
                  call fT   % AdaptGradientsToFace(nGradEqn, N(1), N(2), UxT , UyT , UzT , self % faceSide(ETOP  ))
               end if 
               do m=1,4
                  if (fT % Mortar(m) .ne. 0) then 
                     call fT % AdaptGradientsToMortarFace(nGradEqn, N(1), N(2), UxT , UyT , UzT , self % faceSide(ETOP), faces(fT % Mortar(m)))
                  end if 
               end do 
            elseif(fT % MortarType == MORTAR_SLIDING) then 
               do m=1,size(fT%Mortar)
                  call fT % AdaptGradientsToMortarFace(nGradEqn, N(1), N(2), UxT , UyT , UzT , self % faceSide(ETOP), faces(fT % Mortar(m)), .true.)
               end do 
            else 
               call fT   % AdaptGradientsToFace(nGradEqn, N(1), N(2), UxT , UyT , UzT , self % faceSide(ETOP  ))
            end if 
            
         else
            call fL   % AdaptGradientsToFace(nGradEqn, N(2), N(3), UxL , UyL , UzL , self % faceSide(ELEFT  ))
            call fR   % AdaptGradientsToFace(nGradEqn, N(2), N(3), UxR , UyR , UzR , self % faceSide(ERIGHT ))
            call fFR  % AdaptGradientsToFace(nGradEqn, N(1), N(3), UxFR, UyFR, UzFR, self % faceSide(EFRONT ))
            call fBK  % AdaptGradientsToFace(nGradEqn, N(1), N(3), UxBK, UyBK, UzBK, self % faceSide(EBACK  ))
            call fBOT % AdaptGradientsToFace(nGradEqn, N(1), N(2), UxBT, UyBT, UzBT, self % faceSide(EBOTTOM))
            call fT   % AdaptGradientsToFace(nGradEqn, N(1), N(2), UxT , UyT , UzT , self % faceSide(ETOP   ))
         end if

      end subroutine HexElement_ProlongGradientsToFaces

      subroutine HexElement_ProlongAviscFluxToFaces(self, nEqn, AVflux, fFR, fBK, fBOT, fR, fT, fL, faces)
         use FaceClass
         implicit none
         class(Element),   intent(in)    :: self
         integer,          intent(in)    :: nEqn
         real(kind=RP) ,   intent(in)    :: AVflux(1:NCONS, 0:self%Nxyz(1), 0:self%Nxyz(2), 0:self%Nxyz(3), 1:NDIM)
         class(Face),      intent(inout) :: fFR, fBK, fBOT, fR, fT, fL
         type(Face),optional, intent(inout) :: faces(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                                              :: i, j, k, l, m, N(3)
         real(kind=RP), dimension(1:nEqn, 0:self % Nxyz(1), 0:self % Nxyz(3)) :: VFR, VBK
         real(kind=RP), dimension(1:nEqn, 0:self % Nxyz(1), 0:self % Nxyz(2)) :: VBOT, VT
         real(kind=RP), dimension(1:nEqn, 0:self % Nxyz(2), 0:self % Nxyz(3)) :: VL, VR
         type(NodalStorage_t), pointer                                        :: spAxi, spAeta, spAzeta

         N = self % Nxyz
         spAxi   => NodalStorage(N(1))
         spAeta  => NodalStorage(N(2))
         spAzeta => NodalStorage(N(3))
!
!        *************************
!        Prolong solution to faces
!        *************************
!
         VL   = 0.0_RP     ; VR   = 0.0_RP
         VFR  = 0.0_RP     ; VBK  = 0.0_RP
         VBOT = 0.0_RP     ; VT   = 0.0_RP

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            VL  (:,j,k) = VL  (:,j,k) - AVflux(:,i,j,k,IX) * spAxi   % v(i,LEFT  )
            VR  (:,j,k) = VR  (:,j,k) + AVflux(:,i,j,k,IX) * spAxi   % v(i,RIGHT )
            VFR (:,i,k) = VFR (:,i,k) - AVflux(:,i,j,k,IY) * spAeta  % v(j,FRONT )
            VBK (:,i,k) = VBK (:,i,k) + AVflux(:,i,j,k,IY) * spAeta  % v(j,BACK  )
            VBOT(:,i,j) = VBOT(:,i,j) - AVflux(:,i,j,k,IZ) * spAzeta % v(k,BOTTOM)
            VT  (:,i,j) = VT  (:,i,j) + AVflux(:,i,j,k,IZ) * spAzeta % v(k,TOP   )
         end do                   ; end do                   ; end do
         nullify (spAxi, spAeta, spAzeta)
         if (present(faces)) then 

            if (fL % MortarType == MORTAR_BIG) then 
               if (fL % n_mpi_mortar .ne. 0) then 
                  call fL   % AdaptAviscFluxToFace(nEqn, N(2), N(3), VL, self % faceSide(ELEFT  ))
               end if 
               do m=1,4
                  if (fL % Mortar(m) .ne. 0) then 
                     call fL % AdaptAviscFluxToMortarFace(nEqn, N(2), N(3), VL, self % faceSide(ELEFT), faces(fL % Mortar(m)))
                  end if 
               end do 
            elseif (fL % MortarType == MORTAR_SLIDING) then 
               do m=1,size(fL%Mortar)
                  call fL % AdaptAviscFluxToMortarFace(nEqn, N(2), N(3), VL, self % faceSide(ELEFT), faces(fL % Mortar(m)), .true.)
               end do 
            else 
               call fL   % AdaptAviscFluxToFace(nEqn, N(2), N(3), VL, self % faceSide(ELEFT  ))
            end if 

            if (fR % MortarType == MORTAR_BIG) then 
               if (fR % n_mpi_mortar .ne. 0) then 
                  call fR   % AdaptAviscFluxToFace(nEqn, N(2), N(3), VR, self % faceSide(ERIGHT  ))
               end if 
               do m=1,4
                  if (fR % Mortar(m) .ne. 0) then 
                     call fR % AdaptAviscFluxToMortarFace(nEqn, N(2), N(3), VR, self % faceSide(ERIGHT), faces(fR % Mortar(m)))
                  end if 
               end do 
            elseif (fR % MortarType == MORTAR_SLIDING) then 
               do m=1,size(fR%Mortar)
                  call fR % AdaptAviscFluxToMortarFace(nEqn, N(2), N(3), VR, self % faceSide(ERIGHT), faces(fR % Mortar(m)), .true.)
               end do 
            else 
               call fR   % AdaptAviscFluxToFace(nEqn, N(2), N(3), VR, self % faceSide(ERIGHT  ))
            end if 

            if (fFR % MortarType == MORTAR_BIG) then 
               if (fFR % n_mpi_mortar .ne. 0) then 
                  call fFR   % AdaptAviscFluxToFace(nEqn, N(1), N(3), VFR, self % faceSide(EFRONT  ))
               end if 
               do m=1,4
                  if (fFR % Mortar(m) .ne. 0) then 
                     call fFR % AdaptAviscFluxToMortarFace(nEqn, N(1), N(3), VFR, self % faceSide(EFRONT), faces(fFR % Mortar(m)))
                  end if 
               end do 
            elseif (fFR % MortarType == MORTAR_SLIDING) then 
               do m=1,size(fFR%Mortar)
                  call fFR % AdaptAviscFluxToMortarFace(nEqn, N(1), N(3), VFR, self % faceSide(EFRONT), faces(fFR % Mortar(m)), .true.)
               end do 
            else 
               call fFR   % AdaptAviscFluxToFace(nEqn, N(1), N(3), VFR, self % faceSide(EFRONT  ))
            end if 

            if (fBK % MortarType == MORTAR_BIG) then 
               if (fBK % n_mpi_mortar .ne. 0) then 
                  call fBK   % AdaptAviscFluxToFace(nEqn, N(1), N(3), VBK, self % faceSide(EBACK  ))
               end if 
               do m=1,4
                  if (fBK % Mortar(m) .ne. 0) then 
                     call fBK % AdaptAviscFluxToMortarFace(nEqn, N(1), N(3), VBK, self % faceSide(EBACK), faces(fBK % Mortar(m)))
                  end if 
               end do 
            elseif (fBK % MortarType == MORTAR_SLIDING) then 
               do m=1,size(fBK%Mortar)
                  call fBK % AdaptAviscFluxToMortarFace(nEqn, N(1), N(3), VBK, self % faceSide(EBACK), faces(fBK % Mortar(m)), .true.)
               end do 
            else 
               call fBK   % AdaptAviscFluxToFace(nEqn, N(1), N(3), VBK, self % faceSide(EBACK  ))
            end if 

            if (fBOT % MortarType == MORTAR_BIG) then 
               if (fBOT % n_mpi_mortar .ne. 0) then 
                  call fBOT   % AdaptAviscFluxToFace(nEqn, N(1), N(2), VBOT, self % faceSide(EBOTTOM  ))
               end if 
               do m=1,4
                  if (fBOT % Mortar(m) .ne. 0) then 
                     call fBOT % AdaptAviscFluxToMortarFace(nEqn, N(1), N(2), VBOT, self % faceSide(EBOTTOM), faces(fBOT % Mortar(m)))
                  end if 
               end do 
            elseif (fBOT % MortarType == MORTAR_SLIDING) then 
               do m=1,size(fBOT%Mortar)
                  call fBOT % AdaptAviscFluxToMortarFace(nEqn, N(1), N(2), VBOT, self % faceSide(EBOTTOM), faces(fBOT % Mortar(m)), .true.)
               end do 
            else 
               call fBOT   % AdaptAviscFluxToFace(nEqn, N(1), N(2), VBOT, self % faceSide(EBOTTOM  ))
            end if 

            if (fT % MortarType == MORTAR_BIG) then 
               if (fT % n_mpi_mortar .ne. 0) then 
                  call fT   % AdaptAviscFluxToFace(nEqn, N(1), N(2), VT, self % faceSide(ETOP  ))
               end if 
               do m=1,4
                  if (fT % Mortar(m) .ne. 0) then 
                     call fT % AdaptAviscFluxToMortarFace(nEqn, N(1), N(2), VT, self % faceSide(ETOP), faces(fT % Mortar(m)))
                  end if 
               end do 
            elseif (fT % MortarType == MORTAR_SLIDING) then 
               do m=1,size(fT%Mortar)
                  call fT % AdaptAviscFluxToMortarFace(nEqn, N(1), N(2), VT, self % faceSide(ETOP), faces(fT % Mortar(m)), .true.)
               end do 
            else 
               call fT   % AdaptAviscFluxToFace(nEqn, N(1), N(2), VT, self % faceSide(ETOP  ))
            end if 

         else
            call fL   % AdaptAviscFluxToFace(nEqn, N(2), N(3), VL   , self % faceSide(ELEFT  ))
            call fR   % AdaptAviscFluxToFace(nEqn, N(2), N(3), VR   , self % faceSide(ERIGHT ))
            call fFR  % AdaptAviscFluxToFace(nEqn, N(1), N(3), VFR  , self % faceSide(EFRONT ))
            call fBK  % AdaptAviscFluxToFace(nEqn, N(1), N(3), VBK  , self % faceSide(EBACK  ))
            call fBOT % AdaptAviscFluxToFace(nEqn, N(1), N(2), VBOT , self % faceSide(EBOTTOM))
            call fT   % AdaptAviscFluxToFace(nEqn, N(1), N(2), VT   , self % faceSide(ETOP   ))
         end if

      end subroutine HexElement_ProlongAviscFluxToFaces
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexElement_ComputeLocalGradient(self, nEqn, nGradEqn, GetGradientValues, set_mu)
!
!        ****************************************************************
!           This subroutine computes local gradients as:
!
!              nabla U = (1/J) \sum_i Ja^i \cdot \partial(u)/ \partial \xi^i
!
!        ****************************************************************
!
         implicit none
         class(Element),   intent(inout)  :: self
         integer,          intent(in)     :: nEqn
         integer,          intent(in)     :: nGradEqn
         procedure(GetGradientValues_f)   :: GetGradientValues
         logical,          intent(in)     :: set_mu
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k, l
         real(kind=RP)  :: U(1:nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(2), 0:self % Nxyz(3))
         real(kind=RP)  :: U_xi(1:nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(2), 0:self % Nxyz(3))
         real(kind=RP)  :: U_eta(1:nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(2), 0:self % Nxyz(3))
         real(kind=RP)  :: U_zeta(1:nGradEqn, 0:self % Nxyz(1), 0:self % Nxyz(2), 0:self % Nxyz(3))
         type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta

         associate( N => self % Nxyz )
         spAxi   => NodalStorage(N(1))
         spAeta  => NodalStorage(N(2))
         spAzeta => NodalStorage(N(3))
!
!        **********************
!        Get gradient variables
!        **********************
!
#ifdef MULTIPHASE
            do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
               call GetGradientValues(nEqn, nGradEqn, self % storage % Q(:,i,j,k), U(:,i,j,k), self % storage % rho(i,j,k) )
            end do         ; end do         ; end do
#else
            do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
               call GetGradientValues(nEqn, nGradEqn, self % storage % Q(:,i,j,k), U(:,i,j,k))
            end do         ; end do         ; end do
#endif
   
#ifdef MULTIPHASE
!
!        The multiphase solver needs the Chemical potential as first entropy variable
!        ----------------------------------------------------------------------------
         if ( set_mu ) U(IGMU,:,:,:) = self % storage % mu(1,:,:,:)
#endif
!
!        ************
!        Compute U_xi
!        ************
!
         U_xi = 0.0_RP
         do k = 0, N(3)   ; do j = 0, N(2) ; do l = 0, N(1) ; do i = 0, N(1)
            U_xi(:,i,j,k) = U_xi(:,i,j,k) + U(:,l,j,k) * spAxi % D(i,l)
         end do           ; end do         ; end do         ; end do
!
!        *************
!        Compute U_eta
!        *************
!
         U_eta = 0.0_RP
         do k = 0, N(3)   ; do l = 0, N(2) ; do j = 0, N(2) ; do i = 0, N(1)
            U_eta(:,i,j,k) = U_eta(:,i,j,k) + U(:,i,l,k) * spAeta % D(j,l)
         end do           ; end do         ; end do         ; end do
!
!        **************
!        Compute U_zeta
!        **************
!
         U_zeta = 0.0_RP
         do l = 0, N(3)   ; do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            U_zeta(:,i,j,k) = U_zeta(:,i,j,k) + U(:,i,j,l) * spAzeta % D(k,l)
         end do           ; end do         ; end do         ; end do

         nullify (spAxi, spAeta, spAzeta)
!
!        ******************************
!        Project to the cartesian basis
!        ******************************
!
         do k = 0, N(3) ; do j = 0, N(2)  ; do i = 0, N(1)
            self % storage % U_x(:,i,j,k) = (   U_xi(:,i,j,k) * self % geom % jGradXi(1,i,j,k) &
                                                + U_eta(:,i,j,k) * self % geom % jGradEta(1,i,j,k) &
                                                + U_zeta(:,i,j,k) * self % geom % jGradZeta(1,i,j,k))&
                                             * self % geom % InvJacobian(i,j,k)

            self % storage % U_y(:,i,j,k) = (    U_xi(:,i,j,k) * self % geom % jGradXi(2,i,j,k) &
                                                + U_eta(:,i,j,k) * self % geom % jGradEta(2,i,j,k) &
                                                + U_zeta(:,i,j,k) * self % geom % jGradZeta(2,i,j,k))&
                                             * self % geom % InvJacobian(i,j,k)

            self % storage % U_z(:,i,j,k) = (    U_xi(:,i,j,k) * self % geom % jGradXi(3,i,j,k) &
                                                + U_eta(:,i,j,k) * self % geom % jGradEta(3,i,j,k) &
                                                + U_zeta(:,i,j,k) * self % geom % jGradZeta(3,i,j,k))&
                                             * self % geom % InvJacobian(i,j,k)
         end do         ; end do          ; end do
         end associate

      end subroutine HexElement_ComputeLocalGradient
   !
   !////////////////////////////////////////////////////////////////////////
   !
   !     ----------------------------------------------------------------------------
   !     Checks if a point is inside a linearized version of the element (flat faces)
   !     ----------------------------------------------------------------------------
      logical function HexElement_FindPointInLinearizedElement(self, x, nodes)
         use NodeClass, only: Node
         use Utilities, only: SolveThreeEquationLinearSystem
         implicit none
         !-arguments-------------------------------------------------
         class(Element), intent(in) :: self
         real(kind=RP) , intent(in) :: x(NDIM)
         type(Node)    , intent(in) :: nodes(:)
         !-local-variables-------------------------------------------
         real(kind=RP)  :: F(NDIM)
         real(kind=RP)  :: xi(NDIM)
         real(kind=RP)  :: dx(NDIM)
         real(kind=RP)  :: Jac(NDIM,NDIM)
         real(kind=RP)  :: corners(NDIM,NODES_PER_ELEMENT)
         integer        :: i
!        Newton iterative solver parameters
!        ----------------------------------
         integer,       parameter   :: N_MAX_ITER = 50
         real(kind=RP), parameter   :: TOL = 1.0e-12_RP
         real(kind=RP), parameter   :: STEP = 1.0_RP
         real(kind=RP), parameter   :: INSIDE_TOL = 1.0e-08_RP
         !-arguments-------------------------------------------------


         do i = 1, NODES_PER_ELEMENT
            corners(:,i) = nodes(self % nodeIDs(i)) % x
         end do

!
!        Initial seed
!        ------------
         xi = 0.0_RP

         do i = 1 , N_MAX_ITER

            call Hex8TransfiniteMap ( xi, F, corners )

            F = F - x
!
!           Stopping criteria: there are several
!           ------------------------------------
            if ( maxval(abs(F)) .lt. TOL ) exit
            if ( abs(xi(1)) .ge. 2.5_RP ) exit
            if ( abs(xi(2)) .ge. 2.5_RP ) exit
            if ( abs(xi(3)) .ge. 2.5_RP ) exit
!
!           Perform a step
!           --------------

            call GradHex8TransfiniteMap (xi, Jac, corners)

            dx = solveThreeEquationLinearSystem( Jac , -F )
            xi = xi + STEP * dx

         end do

         if ( (abs(xi(1)) .lt. 1.0_RP + INSIDE_TOL) .and. &
               (abs(xi(2)) .lt. 1.0_RP + INSIDE_TOL) .and. &
               (abs(xi(3)) .lt. 1.0_RP + INSIDE_TOL)          ) then
!
!           Solution is valid
!           -----------------
            HexElement_FindPointInLinearizedElement = .true.

         else
!
!           Solution is not valid
!           ---------------------
            HexElement_FindPointInLinearizedElement = .false.

         end if


      end function HexElement_FindPointInLinearizedElement

!
!////////////////////////////////////////////////////////////////////////
!
      logical function HexElement_FindPointWithCoords(self, x, dir2D, xi)
!
!        **********************************************************
!
!           This function finds whether a point is inside or not
!           of the element. This is done solving
!           the mapping non-linear system
!
!        **********************************************************
!
!
         use Utilities, only: SolveThreeEquationLinearSystem
         implicit none
         class(Element),      intent(in)  :: self
         real(kind=RP),       intent(in)  :: x(NDIM)
         integer,             intent(in)  :: dir2D
         real(kind=RP),       intent(out) :: xi(NDIM)
!
!        ----------------------------------
!        Newton iterative solver parameters
!        ----------------------------------
!
         integer,       parameter   :: N_MAX_ITER = 50
         real(kind=RP), parameter   :: TOL = 1.0e-12_RP
         integer,       parameter   :: STEP = 1.0_RP
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i, j, k, iter
         real(kind=RP), parameter      :: INSIDE_TOL = 1.0e-08_RP
         real(kind=RP)                 :: lxi   (0:self % Nxyz(1))
         real(kind=RP)                 :: leta  (0:self % Nxyz(2))
         real(kind=RP)                 :: lzeta (0:self % Nxyz(3))
         real(kind=RP)                 :: dlxi   (0:self % Nxyz(1))
         real(kind=RP)                 :: dleta  (0:self % Nxyz(2))
         real(kind=RP)                 :: dlzeta (0:self % Nxyz(3))
         real(kind=RP)                 :: F(NDIM)
         real(kind=RP)                 :: Jac(NDIM,NDIM)
         real(kind=RP)                 :: dx(NDIM)
         type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta

         spAxi   => NodalStorage(self % Nxyz(1))
         spAeta  => NodalStorage(self % Nxyz(2))
         spAzeta => NodalStorage(self % Nxyz(3))

!
!        Initial seed
!        ------------
         xi = 0.0_RP

         do iter = 1 , N_MAX_ITER
!
!           Get Lagrange polynomials and derivatives
!           ----------------------------------------
            lxi     = spAxi % lj   (xi(1))
            leta    = spAeta % lj  (xi(2))
            lzeta   = spAzeta % lj (xi(3))

            F = 0.0_RP
            do k = 0, spAzeta % N   ; do j = 0, spAeta % N ; do i = 0, spAxi % N
               F = F + self % geom % x(:,i,j,k) * lxi(i) * leta(j) * lzeta(k)
            end do               ; end do             ; end do

            F = F - x
!
!           Stopping criteria: there are several
!           ------------------------------------
            if ( dir2D .gt. 0 ) F(dir2D) = 0.0_RP
            if ( maxval(abs(F)) .lt. TOL ) exit
            if ( abs(xi(1)) .ge. 2.5_RP ) exit
            if ( abs(xi(2)) .ge. 2.5_RP ) exit
            if ( abs(xi(3)) .ge. 2.5_RP ) exit
!
!           Perform a step
!           --------------
            dlxi    = spAxi % dlj  (xi(1))
            dleta   = spAeta % dlj (xi(2))
            dlzeta  = spAzeta % dlj(xi(3))

            Jac = 0.0_RP
            do k = 0, spAzeta % N   ; do j = 0, spAeta % N ; do i = 0, spAxi % N
               Jac(:,1) = Jac(:,1) + self % geom % x(:,i,j,k) * dlxi(i) * leta(j) * lzeta(k)
               Jac(:,2) = Jac(:,2) + self % geom % x(:,i,j,k) * lxi(i) * dleta(j) * lzeta(k)
               Jac(:,3) = Jac(:,3) + self % geom % x(:,i,j,k) * lxi(i) * leta(j) * dlzeta(k)
            end do               ; end do             ; end do

            if ( (all(abs(Jac(:,1)) .lt. epsilon(1.0_RP))) .and. (spAxi % N .eq. 0) ) then
               Jac(dir2D,1) = 1.0_RP
            end if

            if ( all(abs(Jac(:,2)) .lt. epsilon(1.0_RP)) .and. (spAeta % N .eq. 0)) then
               Jac(dir2D,2) = 1.0_RP
            end if

            if ( all(abs(Jac(:,3)) .lt. epsilon(1.0_RP)) .and. (spAzeta % N .eq. 0)) then
               Jac(dir2D,3) = 1.0_RP
            end if

            dx = solveThreeEquationLinearSystem( Jac , -F )
            xi = xi + STEP * dx

         end do

         nullify (spAxi, spAeta, spAzeta)

         if ( (abs(xi(1)) .lt. 1.0_RP + INSIDE_TOL) .and. &
               (abs(xi(2)) .lt. 1.0_RP + INSIDE_TOL) .and. &
               (abs(xi(3)) .lt. 1.0_RP + INSIDE_TOL)          ) then
!
!           Solution is valid
!           -----------------
            HexElement_FindPointWithCoords = .true.

         else
!
!           Solution is not valid
!           ---------------------
            HexElement_FindPointWithCoords = .false.

         end if

      end function HexElement_FindPointWithCoords
   
      function HexElement_EvaluateSolutionAtPoint(self, nEqn, xi)
         implicit none
         class(Element),   intent(in)    :: self
         integer,          intent(in)    :: nEqn
         real(kind=RP),    intent(in)    :: xi(NDIM)
         real(kind=RP)                   :: HexElement_EvaluateSolutionAtPoint(nEqn)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, k
         real(kind=RP)  :: lxi(0:self % Nxyz(1))
         real(kind=RP)  :: leta(0:self % Nxyz(2))
         real(kind=RP)  :: lzeta(0:self % Nxyz(3))
         real(kind=RP)  :: Q(nEqn)
         type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta

         spAxi   => NodalStorage(self % Nxyz(1))
         spAeta  => NodalStorage(self % Nxyz(2))
         spAzeta => NodalStorage(self % Nxyz(3))
!
!        Compute Lagrange basis
!        ----------------------
         lxi   = spAxi % lj(xi(1))
         leta  = spAeta % lj(xi(2))
         lzeta = spAzeta % lj(xi(3))
!
!        Compute the tensor product
!        --------------------------
         Q = 0.0_RP

         do k = 0, spAzeta % N   ; do j = 0, spAeta % N ; do i = 0, spAxi % N
            Q = Q + self % storage % Q(:,i,j,k) * lxi(i) * leta(j) * lzeta(k)
         end do               ; end do             ; end do

         HexElement_EvaluateSolutionAtPoint = Q

         nullify (spAxi, spAeta, spAzeta)
      end function HexElement_EvaluateSolutionAtPoint

      function HexElement_EvaluateGradientAtPoint(self, nEqn, xi, dir)
         implicit none
         class(Element),   intent(in)    :: self
         integer,          intent(in)    :: nEqn, dir
         real(kind=RP),    intent(in)    :: xi(NDIM)
         real(kind=RP)                   :: HexElement_EvaluateGradientAtPoint(nEqn)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, k
         real(kind=RP)  :: lxi(0:self % Nxyz(1))
         real(kind=RP)  :: leta(0:self % Nxyz(2))
         real(kind=RP)  :: lzeta(0:self % Nxyz(3))
         real(kind=RP)  :: U(nEqn)
         type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta

         spAxi   => NodalStorage(self % Nxyz(1))
         spAeta  => NodalStorage(self % Nxyz(2))
         spAzeta => NodalStorage(self % Nxyz(3))
!
!        Compute Lagrange basis
!        ----------------------
         lxi   = spAxi % lj(xi(1))
         leta  = spAeta % lj(xi(2))
         lzeta = spAzeta % lj(xi(3))
!
!        Compute the tensor product
!        --------------------------
         U = 0.0_RP

         select case( dir )
         case(IX)
            do k = 0, spAzeta % N   ; do j = 0, spAeta % N ; do i = 0, spAxi % N
               U = U + self % storage % U_x(:,i,j,k) * lxi(i) * leta(j) * lzeta(k)
            end do               ; end do             ; end do
         case(IY)
            do k = 0, spAzeta % N   ; do j = 0, spAeta % N ; do i = 0, spAxi % N
               U = U + self % storage % U_y(:,i,j,k) * lxi(i) * leta(j) * lzeta(k)
            end do               ; end do             ; end do
         case(IZ)
            do k = 0, spAzeta % N   ; do j = 0, spAeta % N ; do i = 0, spAxi % N
               U = U + self % storage % U_z(:,i,j,k) * lxi(i) * leta(j) * lzeta(k)
            end do               ; end do             ; end do
         end select

         HexElement_EvaluateGradientAtPoint = U

         nullify (spAxi, spAeta, spAzeta)
         
      end function HexElement_EvaluateGradientAtPoint
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------------
!  Adapts an element to new polynomial orders NNew
!  -> TODO: Previous solutions are not implemented
!  --------------------------------------------------------
      subroutine HexElement_pAdapt (self, NNew, nodes, saveGradients, prevSol_num)
         implicit none
         !-arguments--------------------------------------------
         class(Element), intent(inout) :: self
         integer       , intent(in)    :: NNew(NDIM)
         integer       , intent(in)    :: nodes
         logical       , intent(in)    :: saveGradients
         integer       , intent(in)    :: prevSol_num
         !-arguments--------------------------------------------
         logical                       :: anJacobian
         type(ElementStorage_t)        :: tempStorage
#if (!defined(NAVIERSTOKES))
         logical, parameter            :: computeGradients = .true.
#endif
         !-----------------------------------------------------

!
!        Reconstruct storage
!        -------------------

         anJacobian = self % storage % anJacobian

         call tempStorage % construct (self % Nxyz(1), self % Nxyz(2), self % Nxyz(3), computeGradients, anJacobian, prevSol_num,0)
         tempStorage = self % storage

         self % Nxyz = NNew
         self % hn = (self % geom % Volume / product(self % Nxyz + 1)) ** (1.0_RP / 3.0_RP)

         call self % storage % destruct()
         call self % storage % construct ( NNew(1), NNew(2), NNew(3), computeGradients, anJacobian, prevSol_num,0)

         call tempStorage % InterpolateSolution (self % storage, nodes, saveGradients)

         if (prevSol_num > 0) then
            ! TODO : call InterpolatePrevSol
         end if
         call tempStorage % destruct()

      end subroutine HexElement_pAdapt

      pure subroutine SurfInfo_Destruct (self)
         implicit none
         class(SurfInfo_t), intent(inout) :: self

         call self % facePatches % destruct
      end subroutine SurfInfo_Destruct

      elemental subroutine HexElement_Assign(to, from)
         implicit none
         class(Element), intent(inout) :: to
         class(Element), intent(in)    :: from

         to % hasSharedFaces = from % hasSharedFaces
         to % dir2D = from % dir2D
         to % globDir = from % globDir
         to % eID = from % eID
         to % globID = from % globID
         to % offsetIO = from % offsetIO
         to % nodeIDs = from % nodeIDs
         to % faceIDs = from % faceIDs
         to % faceSide = from % faceSide
         to % Nxyz = from % Nxyz
         to % hn = from % hn
         to % geom = from % geom
         to % boundaryName = from % boundaryName
         to % NumberOfConnections = from % NumberOfConnections
         to % Connection = from % Connection
         to % hexMap = from % hexMap
         to % SurfInfo = from % SurfInfo
!~            IGNORE to % storage

      end subroutine HexElement_Assign
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexElement_ConstructIBM( self, Nx, Ny, Nz, NumOfSTL )
         implicit none
         class(Element), intent(inout) :: self
         integer,        intent(in)    :: Nx, Ny, Nz, NumOfSTL  !<  Polynomial orders, num of stl files

         allocate(self% isInsideBody(0:Nx,0:Ny,0:Nz))
         allocate(self% isForcingPoint(0:Nx,0:Ny,0:Nz))
         allocate(self% STL(0:Nx,0:Ny,0:Nz))
         
         self% isInsideBody   = .false.
         self% isForcingPoint = .false.
         self% STL            = 0
         self% IP_index       = 0
   
         end subroutine HexElement_ConstructIBM
#if defined(ACOUSTIC)
subroutine HexElement_ProlongBaseSolutionToFaces(self, nEqn, fFR, fBK, fBOT, fR, fT, fL, faces)
         use FaceClass
         implicit none
         class(Element),   intent(in)  :: self
         integer,          intent(in)  :: nEqn
         class(Face),      intent(inout) :: fFR, fBK, fBOT, fR, fT, fL
         type(Face),optional, intent(inout) :: faces(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k, N(3), m
         real(kind=RP), dimension(1:nEqn, 0:self % Nxyz(1), 0:self % Nxyz(3)) :: QFR, QBK
         real(kind=RP), dimension(1:nEqn, 0:self % Nxyz(1), 0:self % Nxyz(2)) :: QBOT, QT
         real(kind=RP), dimension(1:nEqn, 0:self % Nxyz(2), 0:self % Nxyz(3)) :: QL, QR
         type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta

         N = self % Nxyz
         spAxi   => NodalStorage(N(1))
         spAeta  => NodalStorage(N(2))
         spAzeta => NodalStorage(N(3))
!
!        *************************
!        Prolong solution to faces
!        *************************
!
         QL   = 0.0_RP     ; QR   = 0.0_RP
         QFR  = 0.0_RP     ; QBK  = 0.0_RP
         QBOT = 0.0_RP     ; QT   = 0.0_RP

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            QL  (:,j,k)= QL  (:,j,k)+ self % storage % Qbase(:,i,j,k)* spAxi % v  (i,LEFT  )
            QR  (:,j,k)= QR  (:,j,k)+ self % storage % Qbase(:,i,j,k)* spAxi % v  (i,RIGHT )
            QFR (:,i,k)= QFR (:,i,k)+ self % storage % Qbase(:,i,j,k)* spAeta % v (j,FRONT )
            QBK (:,i,k)= QBK (:,i,k)+ self % storage % Qbase(:,i,j,k)* spAeta % v (j,BACK  )
            QBOT(:,i,j)= QBOT(:,i,j)+ self % storage % Qbase(:,i,j,k)* spAzeta % v(k,BOTTOM)
            QT  (:,i,j)= QT  (:,i,j)+ self % storage % Qbase(:,i,j,k)* spAzeta % v(k,TOP   )
         end do                   ; end do                   ; end do
         nullify (spAxi, spAeta, spAzeta)

         if (present(faces)) then
            if (fL % MortarType == MORTAR_BIG) then
               if (fL % n_mpi_mortar .ne. 0) then
                  call fL   % AdaptBaseSolutionToFace(nEqn, N(2), N(3), QL, self % faceSide(ELEFT  ))
               end if
               do m=1,4
                  if (fL % Mortar(m) .ne. 0) then
                     call fL % AdaptBaseSolutionToMortarFace(nEqn, N(2), N(3), QL, self % faceSide(ELEFT  ), faces(fL % Mortar(m)))
                  end if
               end do
            elseif (fL % MortarType == MORTAR_SLIDING) then
               error stop 'ProlongBaseSolutionToFaces: sliding mortars not supported for base solution'
            else
               call fL   % AdaptBaseSolutionToFace(nEqn, N(2), N(3), QL, self % faceSide(ELEFT  ))
            end if
            if (fR % MortarType == MORTAR_BIG) then
               if (fR % n_mpi_mortar .ne. 0) then
                  call fR   % AdaptBaseSolutionToFace(nEqn, N(2), N(3), QR, self % faceSide(ERIGHT ))
               end if
               do m=1,4
                  if (fR % Mortar(m) .ne. 0) then
                     call fR % AdaptBaseSolutionToMortarFace(nEqn, N(2), N(3), QR, self % faceSide(ERIGHT ), faces(fR % Mortar(m)))
                  end if
               end do
            elseif (fR % MortarType == MORTAR_SLIDING) then
               error stop 'ProlongBaseSolutionToFaces: sliding mortars not supported for base solution'
            else
               call fR   % AdaptBaseSolutionToFace(nEqn, N(2), N(3), QR, self % faceSide(ERIGHT ))
            end if
            if (fFR % MortarType == MORTAR_BIG) then
               if (fFR % n_mpi_mortar .ne. 0) then
                  call fFR   % AdaptBaseSolutionToFace(nEqn, N(1), N(3), QFR, self % faceSide(EFRONT ))
               end if
               do m=1,4
                  if (fFR % Mortar(m) .ne. 0) then
                     call fFR % AdaptBaseSolutionToMortarFace(nEqn, N(1), N(3), QFR, self % faceSide(EFRONT ), faces(fFR % Mortar(m)))
                  end if
               end do
            elseif (fFR % MortarType == MORTAR_SLIDING) then
               error stop 'ProlongBaseSolutionToFaces: sliding mortars not supported for base solution'
            else
               call fFR   % AdaptBaseSolutionToFace(nEqn, N(1), N(3), QFR, self % faceSide(EFRONT ))
            end if
            if (fBK % MortarType == MORTAR_BIG) then
               if (fBK % n_mpi_mortar .ne. 0) then
                  call fBK   % AdaptBaseSolutionToFace(nEqn, N(1), N(3), QBK, self % faceSide(EBACK  ))
               end if
               do m=1,4
                  if (fBK % Mortar(m) .ne. 0) then
                     call fBK % AdaptBaseSolutionToMortarFace(nEqn, N(1), N(3), QBK, self % faceSide(EBACK  ), faces(fBK % Mortar(m)))
                  end if
               end do
            elseif (fBK % MortarType == MORTAR_SLIDING) then
               error stop 'ProlongBaseSolutionToFaces: sliding mortars not supported for base solution'
            else
               call fBK   % AdaptBaseSolutionToFace(nEqn, N(1), N(3), QBK, self % faceSide(EBACK  ))
            end if
            if (fBOT % MortarType == MORTAR_BIG) then
               if (fBOT % n_mpi_mortar .ne. 0) then
                  call fBOT   % AdaptBaseSolutionToFace(nEqn, N(1), N(2), QBOT, self % faceSide(EBOTTOM))
               end if
               do m=1,4
                  if (fBOT % Mortar(m) .ne. 0) then
                     call fBOT % AdaptBaseSolutionToMortarFace(nEqn, N(1), N(2), QBOT, self % faceSide(EBOTTOM), faces(fBOT % Mortar(m)))
                  end if
               end do
            elseif (fBOT % MortarType == MORTAR_SLIDING) then
               error stop 'ProlongBaseSolutionToFaces: sliding mortars not supported for base solution'
            else
               call fBOT   % AdaptBaseSolutionToFace(nEqn, N(1), N(2), QBOT, self % faceSide(EBOTTOM))
            end if
            if (fT % MortarType == MORTAR_BIG) then
               if (fT % n_mpi_mortar .ne. 0) then
                  call fT   % AdaptBaseSolutionToFace(nEqn, N(1), N(2), QT, self % faceSide(ETOP   ))
               end if
               do m=1,4
                  if (fT % Mortar(m) .ne. 0) then
                     call fT % AdaptBaseSolutionToMortarFace(nEqn, N(1), N(2), QT, self % faceSide(ETOP  ), faces(fT % Mortar(m)))
                  end if
               end do
            elseif (fT % MortarType == MORTAR_SLIDING) then
               error stop 'ProlongBaseSolutionToFaces: sliding mortars not supported for base solution'
            else
               call fT   % AdaptBaseSolutionToFace(nEqn, N(1), N(2), QT, self % faceSide(ETOP   ))
            end if
         else
            call fL   % AdaptBaseSolutionToFace(nEqn, N(2), N(3), QL  , self % faceSide(ELEFT  ))
            call fR   % AdaptBaseSolutionToFace(nEqn, N(2), N(3), QR  , self % faceSide(ERIGHT ))
            call fFR  % AdaptBaseSolutionToFace(nEqn, N(1), N(3), QFR , self % faceSide(EFRONT ))
            call fBK  % AdaptBaseSolutionToFace(nEqn, N(1), N(3), QBK , self % faceSide(EBACK  ))
            call fBOT % AdaptBaseSolutionToFace(nEqn, N(1), N(2), QBOT, self % faceSide(EBOTTOM))
            call fT   % AdaptBaseSolutionToFace(nEqn, N(1), N(2), QT  , self % faceSide(ETOP   ))
         end if

      end subroutine HexElement_ProlongBaseSolutionToFaces
#endif
! 
END Module ElementClass
   