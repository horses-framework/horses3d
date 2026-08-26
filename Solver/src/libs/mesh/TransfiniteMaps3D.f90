!
! /////////////////////////////////////////////////////////////////////
!
!!     This is a collection of transfinite maps and operations
!!     on those maps.
!
!      PUBLIC METHODS:
!          SUBROUTINE Hex8TransfiniteMap( xi, x, corners )
!          SUBROUTINE GradHex8TransfiniteMap( xi, grad_x, corners )
!          SUBROUTINE GeneralHexTransfiniteMap( xi, x, cornerPoints, faceData )
!          SUBROUTINE GradGeneralHexTransfiniteMap( xi, grad_x, cornerPoints,&
!         &                                        faceData )
!
!      OTHER METHODS:
!          SUBROUTINE ComputeHexTransfiniteMap( xi, x, face, edge, corners )
!          SUBROUTINE ComputeGradHexTransfiniteMap( xi, grad_x, face, faceDer, &
!          &                                       edge, edgeDer, corners )
!
!
! //////////////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!
!!     Hex8TransfiniteMap: Compute the transfinite mapping for a hex8 element
!!>
!!     corners = the 8 corners in standard FE orientation
!!     u       = computational space variable (Xi, Eta, Zeta) in [-1,1]
!!     x       = physical space variable (x, y, z)
!!<
!!
!-------------------------------------------------------------------------------
!
   MODULE TransfiniteMapClass
      USE SMConstants
      USE FacePatchClass
      IMPLICIT NONE

      private
      public TransfiniteHexMap
      
      public GeneralHexTransfiniteMap, GradGeneralHexTransfiniteMap
      public GradHex8TransfiniteMap, Hex8TransfiniteMap
      
      TYPE TransfiniteHexMap
         REAL(KIND=RP)  , DIMENSION(3,8)              :: corners
         TYPE(FacePatch), DIMENSION(:)  , ALLOCATABLE :: faces
!
!        ========
         CONTAINS
!        ========
!         
         PROCEDURE :: constructWithCorners
         PROCEDURE :: constructWithFaces
         PROCEDURE :: destruct => destructTransfiniteHexMap
         PROCEDURE :: transfiniteMapAt
         PROCEDURE :: metricDerivativesAt
         PROCEDURE :: isHex8
         PROCEDURE :: setCorners

      END TYPE TransfiniteHexMap
!
!     ========
      CONTAINS 
!     ========
!
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE constructWithCorners(self, corners)  
         IMPLICIT NONE
         CLASS(TransFiniteHexMap) :: self
         REAL(KIND=RP)           :: corners(3,8)
         
         CALL self % setCorners(corners)
         
      END SUBROUTINE constructWithCorners
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE constructWithFaces(self, faces)  
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(TransFiniteHexMap) :: self
         TYPE(FacePatch)          :: faces(6)
!
!        ---------------
!        Local variables
!        ---------------
!
         REAL(KIND=RP) :: corners(3,8)
!
!        -----------------------
!        Save the boundary faces
!        -----------------------
!
         ALLOCATE( self % faces(6))
         self % faces = faces
!
!        -----------------------------
!        Compute and save the corners 
!        -----------------------------
!
         CALL ComputeFacePoint(faces(3),u = [-1.0_RP,-1.0_RP],p = corners(:,1))
         CALL ComputeFacePoint(faces(3),u = [ 1.0_RP,-1.0_RP],p = corners(:,2))
         CALL ComputeFacePoint(faces(3),u = [ 1.0_RP, 1.0_RP],p = corners(:,3))
         CALL ComputeFacePoint(faces(3),u = [-1.0_RP, 1.0_RP],p = corners(:,4))
         
         CALL ComputeFacePoint(faces(5),u = [-1.0_RP,-1.0_RP],p = corners(:,5))
         CALL ComputeFacePoint(faces(5),u = [ 1.0_RP,-1.0_RP],p = corners(:,6))
         CALL ComputeFacePoint(faces(5),u = [ 1.0_RP, 1.0_RP],p = corners(:,7))
         CALL ComputeFacePoint(faces(5),u = [-1.0_RP, 1.0_RP],p = corners(:,8))
         
         CALL self % setCorners(corners)
         
      END SUBROUTINE constructWithFaces
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE setCorners(self, corners)  
         IMPLICIT NONE  
         CLASS(TransFiniteHexMap) :: self
         REAL(KIND=RP)            :: corners(3,8)

          self % corners = corners

      END SUBROUTINE setCorners
!
!//////////////////////////////////////////////////////////////////////// 
      FUNCTION isHex8(self) RESULT(ans)  
         IMPLICIT NONE  
         CLASS(TransfiniteHexMap) :: self
         LOGICAL :: ans
         
         ans = .TRUE.
         IF(ALLOCATED(self % faces)) ans = .FALSE.
         
      END FUNCTION isHex8
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION transfiniteMapAt(self,u) RESULT(x)  
         IMPLICIT NONE
          CLASS(TransFiniteHexMap) :: self
          REAL(KIND=RP)            :: u(3), x(3)
          
          IF(ALLOCATED(self % faces))     THEN
          
             CALL GeneralHexTransfiniteMap(u = u,x = x,                  &
                                           cornerPoints = self % corners,&
                                           faceData = self % faces)
          ELSE
             CALL Hex8TransfiniteMap(u = u,x = x,corners = self % corners)
          END IF 
      END FUNCTION transfiniteMapAt
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION metricDerivativesAt(self,u) RESULT(grad_x)  
         IMPLICIT NONE
          CLASS(TransFiniteHexMap) :: self
          REAL(KIND=RP)            :: u(3), grad_x(3,3)
          
          IF(ALLOCATED(self % faces))     THEN
             CALL GradGeneralHexTransfiniteMap(u = u,grad_x = grad_x,        &
                                               cornerPoints = self % corners,&
                                               faceData = self % faces)
          ELSE
             CALL GradHex8TransfiniteMap(u = u, grad_x = grad_x, corners = self % corners)
          END IF 
      END FUNCTION metricDerivativesAt
!
!//////////////////////////////////////////////////////////////////////// 
! 
      pure SUBROUTINE destructTransfiniteHexMap(self)
         IMPLICIT NONE  
          CLASS(TransFiniteHexMap), intent(inout) :: self
          integer :: i
          IF(ALLOCATED(self % faces)) then
                  DO i=1, SIZE(self % faces)
                        CALL self % faces(i) % Destruct()
                  END DO
                  DEALLOCATE(self % faces)
          end if
     END SUBROUTINE destructTransfiniteHexMap
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Hex8TransfiniteMap( u, x, corners )
!
      USE SMConstants
      IMPLICIT NONE
!
      REAL(KIND=RP), DIMENSION(3,8), INTENT(IN)  :: corners
      REAL(KIND=RP), DIMENSION(3)  , INTENT(IN)  :: u
      REAL(KIND=RP), DIMENSION(3)  , INTENT(OUT) :: x
!
      INTEGER                                     :: j
      REAL(KIND=RP), DIMENSION(3)                 :: xi
!
      xi = 0.5_RP*(u+1.0_RP)
      DO j = 1,3
         x(j) =   corners(j,1)*(1._RP - xi(1))*(1._RP - xi(2))*(1._RP - xi(3)) &
     &          + corners(j,2)* xi(1)         *(1._RP - xi(2))*(1._RP - xi(3)) &
     &          + corners(j,3)* xi(1)         * xi(2)         *(1._RP - xi(3)) &
     &          + corners(j,4)*(1._RP - xi(1))* xi(2)         *(1._RP - xi(3)) &
     &          + corners(j,5)*(1._RP - xi(1))*(1._RP - xi(2))* xi(3)          &
     &          + corners(j,6)* xi(1)         *(1._RP - xi(2))* xi(3)          &
     &          + corners(j,7)* xi(1)         * xi(2)         * xi(3)          &
     &          + corners(j,8)*(1._RP - xi(1))* xi(2)         * xi(3)
      END DO
!
      RETURN

      END SUBROUTINE Hex8TransfiniteMap
!
!     ///////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!
!!     GradHex8TransfiniteMap: Compute the gradient of the transfinite mapping
!!     for a hex8 element
!!>
!!     corners = the 8 corners in standard FE orientation
!!     u       = computational space variable (Xi, Eta, Zeta) in [-1,1]
!!     grad_x  = physical space gradient
!!               grad_x(1,:) = (x_Xi, x_Eta, x_Zeta), etc., i.e.
!!               grad_x(i,j) = dx_(i)/dXi_(j)
!!<
!!
!-------------------------------------------------------------------------------
!
      SUBROUTINE GradHex8TransfiniteMap( u, grad_x, corners )
!
      USE SMConstants
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(3, 8), INTENT(IN)  :: corners
      REAL(KIND=RP), DIMENSION(3)   , INTENT(IN)  :: u
      REAL(KIND=RP), DIMENSION(3,3) , INTENT(OUT) :: grad_x
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(3) :: xi
      INTEGER                     :: i
!
      xi = 0.5_RP*(u+1.0_RP)
      DO i = 1,3
         grad_x(i,1) =  - corners(i,1)*(1._RP - xi(2))*(1._RP - xi(3)) &
                        + corners(i,2)*(1._RP - xi(2))*(1._RP - xi(3)) &
                        + corners(i,3)* xi(2)             *(1._RP - xi(3)) &
                        - corners(i,4)* xi(2)             *(1._RP - xi(3)) &
                        - corners(i,5)*(1._RP - xi(2))* xi(3)              &
                        + corners(i,6)*(1._RP - xi(2))* xi(3)              &
                        + corners(i,7)* xi(2)             * xi(3)              &
                        - corners(i,8)* xi(2)             * xi(3)
!
         grad_x(i,2) =  - corners(i,1)*(1._RP - xi(1))*(1._RP - xi(3)) &
                        - corners(i,2)* xi(1)             *(1._RP - xi(3)) &
                        + corners(i,3)* xi(1)             *(1._RP - xi(3)) &
                        + corners(i,4)*(1._RP - xi(1))*(1._RP - xi(3)) &
                        - corners(i,5)*(1._RP - xi(1))* xi(3)              &
                        - corners(i,6)* xi(1)             * xi(3)              &
                        + corners(i,7)* xi(1)             * xi(3)              &
                        + corners(i,8)*(1._RP - xi(1))* xi(3)
!
         grad_x(i,3) =  - corners(i,1)*(1._RP - xi(1))*(1._RP - xi(2))&
                        - corners(i,2)* xi(1)             *(1._RP - xi(2)) &
                        - corners(i,3)* xi(1)             * xi(2)              &
                        - corners(i,4)*(1._RP - xi(1))* xi(2)              &
                        + corners(i,5)*(1._RP - xi(1))*(1._RP - xi(2)) &
                        + corners(i,6)* xi(1)             *(1._RP - xi(2)) &
                        + corners(i,7)* xi(1)             * xi(2)              &
                        + corners(i,8)*(1._RP - xi(1))* xi(2)
      END DO
      grad_x = 0.5_RP*grad_x
!
      END SUBROUTINE GradHex8TransfiniteMap
!
!///////////////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!
!!     GeneralHexTransfiniteMap: Given the six faces and four corners of a hex
!!     element, and 
!      a computational point xi, return the physical space location x.
!!>
!!     cornerPoints = the 8 corners in standard FE orientation
!!     faceData     = Interpolation data that defines a face.
!!     u            = Computational space variable (u, v, w) in [-1,1]
!!     x            = Physical space location
!!
!!<
!!
!-------------------------------------------------------------------------------
!
      SUBROUTINE GeneralHexTransfiniteMap( u, x, cornerPoints, faceData )
!
!     ......................................................................
!     Given the six faces and four corners of a hex element, and a
!     computational point xi, return the physical space location x
!     ......................................................................
!
      USE SMConstants

      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)   , DIMENSION(3)   , INTENT(IN)  :: u
      REAL(KIND=RP)   , DIMENSION(3,8) , INTENT(IN)  :: cornerPoints
      TYPE(FacePatch) , DIMENSION(6)   , INTENT(IN)  :: faceData   
      REAL(KIND=RP)   , DIMENSION(3)   , INTENT(OUT) :: x
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) , DIMENSION(3, 6) :: facePoint
      REAL(KIND=RP) , DIMENSION(3,12) :: edgePoint
      REAL(KIND=RP) , DIMENSION(3)    :: xi
!
!     -------------
!     Compute edges
!     -------------
!
      CALL ComputeFacePoint( faceData(1), [u(1), -1.0_RP] , edgePoint(:,1))
      CALL ComputeFacePoint( faceData(1), [1.0_RP, u(3)],   edgePoint(:,2))
      CALL ComputeFacePoint( faceData(1), [u(1), 1.0_RP],   edgePoint(:,3))
      CALL ComputeFacePoint( faceData(1), [-1.0_RP, u(3)],  edgePoint(:,4))
      CALL ComputeFacePoint( faceData(2), [u(1), -1.0_RP] , edgePoint(:,5))
      CALL ComputeFacePoint( faceData(2), [1.0_RP, u(3)],   edgePoint(:,6))
      CALL ComputeFacePoint( faceData(2), [u(1), 1.0_RP],   edgePoint(:,7))
      CALL ComputeFacePoint( faceData(2), [-1.0_RP, u(3)],  edgePoint(:,8))
      CALL ComputeFacePoint( faceData(4), [u(2), -1.0_RP],  edgePoint(:,10))
      CALL ComputeFacePoint( faceData(6), [u(2), -1.0_RP],  edgePoint(:,9))
      CALL ComputeFacePoint( faceData(4), [u(2), 1.0_RP],   edgePoint(:,11))
      CALL ComputeFacePoint( faceData(6), [u(2), 1.0_RP],   edgePoint(:,12))
!
!     -------------
!     Compute faces
!     -------------
!
      CALL ComputeFacePoint(faceData(1), [u(1), u(3)], facePoint(:,1))
      CALL ComputeFacePoint(faceData(2), [u(1), u(3)], facePoint(:,2))
      CALL ComputeFacePoint(faceData(3), [u(1), u(2)], facePoint(:,3))
      CALL ComputeFacePoint(faceData(4), [u(2), u(3)], facePoint(:,4))
      CALL ComputeFacePoint(faceData(5), [u(1), u(2)], facePoint(:,5))
      CALL ComputeFacePoint(faceData(6), [u(2), u(3)], facePoint(:,6))
!
!     -------------------
!     Compute the mapping
!     -------------------
!
      xi = 0.5_RP*(u+1.0_RP)
      CALL ComputeHexTransfiniteMap(xi, x, facePoint, edgePoint, cornerPoints)
!
      END SUBROUTINE GeneralHexTransfiniteMap
!
!///////////////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!
!!     ComputeHexTransfiniteMap: Given the six faces and four corners of a hex 
!!     element, and 
!      a computational point xi, return the physical space location x.
!!>
!!     Compute the transfinite mapping for a general hex element
!
!!     xi      = computational space variable (Xi, Zeta, Eta) in [0,1]
!!     x       = resultant physical space variable (x, y, z)
!!
!!     face    = face interpolant location for appropriate local face coordinate
!!     edge    = edge interpolant location for appropriate edge coordinate
!!     corners = The 8 corners in standard FE format
!!<
!!
!-------------------------------------------------------------------------------
!
      SUBROUTINE ComputeHexTransfiniteMap( xi, x, face, edge, corners )
!
      USE SMConstants
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(3)   , INTENT(IN) :: xi
      REAL(KIND=RP), DIMENSION(3,6) , INTENT(IN) :: face
      REAL(KIND=RP), DIMENSION(3,12), INTENT(IN) :: edge
      REAL(KIND=RP), DIMENSION(3,8) , INTENT(IN) :: corners
!
      REAL(KIND=RP), DIMENSION(3), INTENT(OUT)   :: x
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER :: j
!
      DO j = 1,3
!
!        ------------------
!        Face contributions
!        ------------------
!
         x(j) =  face(j,6)*(1._RP - xi(1)) + face(j,4)*xi(1) &
               + face(j,1)*(1._RP - xi(2)) + face(j,2)*xi(2) &
               + face(j,3)*(1._RP - xi(3)) + face(j,5)*xi(3)
!
!        ------------------
!        Edge contributions
!        ------------------
!
         x(j) = x(j) - edge(j,1) *(1._RP - xi(2))*(1._RP - xi(3)) &
                     - edge(j,3) *(1._RP - xi(2))*         xi(3)  &
                     - edge(j,5) *         xi(2) *(1._RP - xi(3)) &
                     - edge(j,7) *         xi(2) *         xi(3)  &
                     - edge(j,9) *(1._RP - xi(1))*(1._RP - xi(3)) &
                     - edge(j,12)*(1._RP - xi(1))*         xi(3)  &
                     - edge(j,10)*(1._RP - xi(3))*         xi(1)  &
                     - edge(j,11)*         xi(1) *         xi(3)  &
                     - edge(j,4) *(1._RP - xi(1))*(1._RP - xi(2)) &
                     - edge(j,8) *(1._RP - xi(1))*         xi(2)  &
                     - edge(j,2) *         xi(1) *(1._RP - xi(2)) &
                     - edge(j,6) *         xi(1) *         xi(2)  
!
!        ------------------
!        Corner contributions
!        ------------------
!
         x(j) = x(j) + corners(j,1)*(1._RP-xi(1))*(1._RP-xi(2))*(1._RP-xi(3)) &
                     + corners(j,5)*(1._RP-xi(1))*(1._RP-xi(2))*       xi(3)  &
                     + corners(j,4)*(1._RP-xi(1))*       xi(2) *(1._RP-xi(3)) &
                     + corners(j,8)*(1._RP-xi(1))*       xi(2) *       xi(3)  &
                     + corners(j,2)*       xi(1)*(1._RP- xi(2))*(1._RP-xi(3)) &
                     + corners(j,6)*       xi(1)*(1._RP- xi(2))*       xi(3)  &
                     + corners(j,3)*       xi(1)*        xi(2)*(1._RP- xi(3)) &
                     + corners(j,7)*       xi(1)*        xi(2)*        xi(3)
      END DO
!
      END SUBROUTINE ComputeHexTransfiniteMap
!
!     ///////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!
!!     GradGeneralHexTransfiniteMap: Given the six faces and four cornerss of a 
!!     hex element, 
!!     and a computational point xi, return the jacouan matrix 
!!     grad_x = d x_i/d Xi_j
!!
!!>
!!    u            = computational space variable (xi, eta, zeta) in [-1,1]
!!    grad_x = physical space gradient grad_x(:,1) = (x_Xi, x_Eta, x_Zeta), etc.
!!    cornerPoints = the 8 corners in standard fe format
!!    faceData     = interpolant data for the 6 faces
!!<
!!
!-------------------------------------------------------------------------------
!
      SUBROUTINE GradGeneralHexTransfiniteMap( u, grad_x, cornerPoints,&
     &                                        faceData )
!
      USE SMConstants

      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)  , DIMENSION(3)   :: u
      REAL(KIND=RP)  , DIMENSION(3,3) :: grad_x
      REAL(KIND=RP)  , DIMENSION(3,8) :: cornerPoints
      TYPE(FacePatch), DIMENSION(6)   :: faceData
!
!     ---------------
!     Local Variables
!     ---------------
!
      
      REAL(KIND=RP) , DIMENSION(3,2)    :: grad2D
      REAL(KIND=RP) , DIMENSION(3,6)    :: facePoint
      REAL(KIND=RP) , DIMENSION(3,2,6)  :: faceDer
      REAL(KIND=RP) , DIMENSION(3,12)   :: edgePoint
      REAL(KIND=RP) , DIMENSION(3,12)   :: edgeDer
      REAL(KIND=RP) , DIMENSION(3)      :: xi
      REAL(KIND=RP)                     :: eps
      INTEGER                           :: i,j
!
!     -------------
!     Compute edges
!     -------------
!
      CALL ComputeFacePoint(faceData(1),[u(1)  ,-1._RP] , edgePoint(:,1))
      CALL ComputeFacePoint(faceData(1),[1.0_RP ,u(3) ], edgePoint(:,2))
      CALL ComputeFacePoint(faceData(1),[u(1)  ,1.0_RP], edgePoint(:,3))
      CALL ComputeFacePoint(faceData(1),[-1.0_RP,u(3) ], edgePoint(:,4))
      CALL ComputeFacePoint(faceData(2),[u(1)  ,-1._RP] , edgePoint(:,5))
      CALL ComputeFacePoint(faceData(2),[1.0_RP ,u(3) ], edgePoint(:,6))
      CALL ComputeFacePoint(faceData(2),[u(1)  ,1.0_RP], edgePoint(:,7))
      CALL ComputeFacePoint(faceData(2),[-1.0_RP,u(3) ], edgePoint(:,8))
      CALL ComputeFacePoint(faceData(4),[u(2)  ,-1._RP] ,edgePoint(:,10))
      CALL ComputeFacePoint(faceData(6),[u(2)  ,-1._RP] , edgePoint(:,9))
      CALL ComputeFacePoint(faceData(4),[u(2)  ,1.0_RP],edgePoint(:,11))
      CALL ComputeFacePoint(faceData(6),[u(2)  ,1.0_RP],edgePoint(:,12))
!
!     -------------
!     Compute faces
!     -------------
!
      CALL ComputeFacePoint(faceData(1), [u(1), u(3)], facePoint(:,1))
      CALL ComputeFacePoint(faceData(2), [u(1), u(3)], facePoint(:,2))
      CALL ComputeFacePoint(faceData(3), [u(1), u(2)], facePoint(:,3))
      CALL ComputeFacePoint(faceData(4), [u(2), u(3)], facePoint(:,4))
      CALL ComputeFacePoint(faceData(5), [u(1), u(2)], facePoint(:,5))
      CALL ComputeFacePoint(faceData(6), [u(2), u(3)], facePoint(:,6))
!
!     ------------------------
!     Compute edge derivatives
!     ------------------------
!
      CALL ComputeFaceDerivative(faceData(1), [u(1), -1._RP] , grad2D)
      edgeDer(:,1) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(1), [1.0_RP, u(3)], grad2D)
      edgeDer(:,2) = grad2D(:,2)
      CALL ComputeFaceDerivative(faceData(1), [u(1), 1.0_RP], grad2D)
      edgeDer(:,3) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(1), [-1.0_RP, u(3)], grad2D)
      edgeDer(:, 4) = grad2D(:,2)
      CALL ComputeFaceDerivative(faceData(2), [u(1), -1._RP] , grad2D)
      edgeDer(:, 5) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(2), [1.0_RP, u(3)], grad2D)
      edgeDer(:, 6) = grad2D(:,2)
      CALL ComputeFaceDerivative(faceData(2), [u(1), 1.0_RP], grad2D)
      edgeDer(:, 7) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(2), [-1.0_RP, u(3)], grad2D)
      edgeDer(:, 8) = grad2D(:,2)
      CALL ComputeFaceDerivative(faceData(6), [u(2), -1._RP] , grad2D)
      edgeDer(:, 9) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(4), [u(2), -1._RP] , grad2D)
      edgeDer(:,10) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(4), [u(2), 1.0_RP], grad2D)
      edgeDer(:,11) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(6), [u(2), 1.0_RP], grad2D)
      edgeDer(:,12) = grad2D(:,1)
!
!     ------------------------
!     Compute face derivatives
!     ------------------------
!
      CALL ComputeFaceDerivative(faceData(1),[u(1),u(3)], faceDer(:,:,1))
      CALL ComputeFaceDerivative(faceData(2),[u(1),u(3)], faceDer(:,:,2))
      CALL ComputeFaceDerivative(faceData(3),[u(1),u(2)], faceDer(:,:,3))
      CALL ComputeFaceDerivative(faceData(4),[u(2),u(3)], faceDer(:,:,4))
      CALL ComputeFaceDerivative(faceData(5),[u(1),u(2)], faceDer(:,:,5))
      CALL ComputeFaceDerivative(faceData(6),[u(2),u(3)], faceDer(:,:,6))
!
!     -------------------
!     Compute the mapping
!     -------------------
!
      CALL ComputeGradHexTransfiniteMap(u, grad_x, facePoint, 2.0_RP * faceDer, &
     &                                  edgePoint, 2.0_RP * edgeDer, cornerPoints)
!
!     ---------------------------
!     Zero out rounded quantities
!     ---------------------------
!
      eps = 100._RP*EPSILON(eps)
      DO i = 1,3
         DO j = 1,3
            IF( ABS(grad_x(i, j)) <= eps )   grad_x(i,j) = 0.0_RP
         END DO
      END DO
!
!     -------------------------------------------
!     At this point should renormalize the matrix
!     -------------------------------------------
!

!
      RETURN
      END SUBROUTINE GradGeneralHexTransfiniteMap
!
!     ///////////////////////////////////////////////////////////////////////
!
!
!-------------------------------------------------------------------------------
!
!!     ComputeGradHexTransfiniteMap: Compute the gradient for a general hex
!!                                   element after the incoming face data has
!!                                   been processed.
!!
!!>
!!     xi          = computational space variable (Xi, Eta, Zeta) in [0,1]
!!     grad_x(1,:) = (x_Xi, x_Eta, x_Zeta), etc., i.e.
!!                    grad_x(i,j) = dx_(i)/dXi_(j)
!!
!!     face    = face interpolant location for appropriate local face coordinate
!!     faceDer = face interpolant derivatives for appropriate local face 
!!               coordinate
!!     edge    = edge interpolant location for appropriate edge coordinate
!!     edgeDer = edge interpolant derivative for appropriate edge coordinate
!!     corners = The 8 corners in standard fe format
!!<
!!
!-------------------------------------------------------------------------------
!
      SUBROUTINE ComputeGradHexTransfiniteMap( u, grad_x, face, faceDer, &
     &                                        edge, edgeDer, corners )
!
      USE SMConstants
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(3)      , INTENT(IN) :: u
      REAL(KIND=RP), DIMENSION(3,6)    , INTENT(IN) :: face
      REAL(KIND=RP), DIMENSION(3,2,6)  , INTENT(IN) :: faceDer
      REAL(KIND=RP), DIMENSION(3,12)   , INTENT(IN) :: edge
      REAL(KIND=RP), DIMENSION(3,12)   , INTENT(IN) :: edgeDer
      REAL(KIND=RP), DIMENSION(3, 8)   , INTENT(IN) :: corners
!
      REAL(KIND=RP), DIMENSION(3,3), INTENT(OUT)    :: grad_x
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER       :: j
      REAL(KIND=RP) :: xi(3)
!
      xi = 0.5_RP*(u + 1.0_RP)
      DO j = 1,3
!
!        ------------
!        Xi derivative
!        ------------
!
         grad_x(j,1) = - face(j,6) + face(j,4) &
                       + faceDer(j,1,1)*(1._RP - xi(2)) + faceDer(j,1,2)*xi(2) &
                       + faceDer(j,1,3)*(1._RP - xi(3)) + faceDer(j,1,5)*xi(3)
!
         grad_x(j,1) = grad_x(j,1)                                   &
        &               - edgeDer(j,1)*(1._RP - xi(2))*(1._RP - xi(3)) &
                        - edgeDer(j,3)*(1._RP - xi(2))*         xi(3)  &
                        - edgeDer(j,5)*         xi(2) *(1._RP - xi(3)) &
                        - edgeDer(j,7)*         xi(2) *         xi(3)  &
                        + edge(j,9)   *(1._RP - xi(3))                 &
                        + edge(j,12)  *         xi(3)                  &
                        - edge(j,10)  *(1._RP - xi(3))                 &
                        - edge(j,11)  *         xi(3)                  &
                        + edge(j,4)   *(1._RP - xi(2))                 &
                        + edge(j,8)   *         xi(2)                  &
                        - edge(j,2)   *(1._RP - xi(2))                 &
                        - edge(j,6)   *         xi(2)
!
         grad_x(j,1) =  grad_x(j,1) &
                        - corners(j,1)*(1._RP -  xi(2))*(1._RP -  xi(3)) &
                        - corners(j, 5)*(1._RP -  xi(2))*          xi(3)  &
                        - corners(j, 4)*          xi(2) *(1._RP -  xi(3)) &
                        - corners(j, 8)*          xi(2) *          xi(3)  &
                        + corners(j,2)*(1._RP -  xi(2))*(1._RP -  xi(3)) &
                        + corners(j, 6)*(1._RP -  xi(2))*          xi(3)  &
                        + corners(j,3)*          xi(2) *(1._RP -  xi(3)) &
                        + corners(j, 7)*          xi(2) *          xi(3)
!        --------------
!        Eta derivative
!        --------------
!
         grad_x(j,2) =   faceDer(j,1,6)*(1._RP - xi(1)) + faceDer(j,1,4)*xi(1) &
                       - face(j,1) + face(j,2)                               &
                       + faceDer(j,2,3)*(1._RP - xi(3)) + faceDer(j,2,5)*xi(3)
!
         grad_x(j,2) = grad_x(j,2) + edge(j,1)*(1._RP - xi(3))                 &
                                   + edge(j,3)*         xi(3)                  &
                                   - edge(j,5)*(1._RP - xi(3))                 &
                                   - edge(j,7)*         xi(3)                  &
                                   - edgeDer(j,9) *(1._RP-xi(1))*(1._RP-xi(3)) &
                                   - edgeDer(j,12)*(1._RP-xi(1))*       xi(3)  &
                                   - edgeDer(j,10)*(1._RP-xi(3))*       xi(1)  &
                                   - edgeDer(j,11)*       xi(1) *       xi(3)  &
                                   + edge(j,4)    *(1._RP-xi(1))               &
                                   - edge(j,8)    *(1._RP-xi(1))               &
                                   + edge(j,2)    *       xi(1)                &
                                   - edge(j,6)    *       xi(1)
!
         grad_x(j,2) = grad_x(j,2) - corners(j,1)*(1._RP- xi(1))*(1._RP-xi(3)) &
                                   - corners(j,5)*(1._RP- xi(1))*       xi(3)  &
                                   + corners(j,4)*(1._RP- xi(1))*(1._RP-xi(3)) &
                                   + corners(j,8)*(1._RP- xi(1))*       xi(3)  &
                                   - corners(j,2)*        xi(1) *(1._RP-xi(3)) &
                                   - corners(j,6)*        xi(1) *       xi(3)  &
                                   + corners(j,3)*        xi(1) *(1._RP-xi(3)) &
                                   + corners(j,7)*        xi(1) *       xi(3)
!        ---------------
!        Zeta derivative
!        ---------------
!
         grad_x(j,3) =  faceDer(j,2,6)*(1._RP - xi(1)) + faceDer(j,2,4)*xi(1) &
                      + faceDer(j,2,1)*(1._RP - xi(2)) + faceDer(j,2,2)*xi(2) &
                      - face(j,3) + face(j,5)
!
         grad_x(j,3) = grad_x(j,3) + edge(j,1)   *(1._RP- xi(2))               &
                                   - edge(j,3)   *(1._RP- xi(2))               &
                                   + edge(j,5)   *        xi(2)                &
                                   - edge(j,7)   *        xi(2)                &
                                   + edge(j,9)   *(1._RP- xi(1))               &
                                   - edge(j,12)  *(1._RP- xi(1))               &
                                   + edge(j,10)  *        xi(1)                &
                                   - edge(j,11)  *        xi(1)                &
                                   - edgeDer(j,4)*(1._RP- xi(1))*(1._RP-xi(2)) &
                                   - edgeDer(j,8) *(1._RP-xi(1))*       xi(2)  &
                                   - edgeDer(j,2) *       xi(1) *(1._RP-xi(2)) &
                                   - edgeDer(j,6) *       xi(1) *       xi(2)
!
         grad_x(j,3) = grad_x(j,3) - corners(j,1)*(1._RP- xi(1))*(1._RP-xi(2)) &
                                   + corners(j,5)*(1._RP- xi(1))*(1._RP-xi(2)) &
                                   - corners(j,4)*(1._RP- xi(1))*       xi(2)  &
                                   + corners(j,8)*(1._RP- xi(1))*       xi(2)  &
                                   - corners(j,2)*        xi(1) *(1._RP-xi(2)) &
                                   + corners(j,6)*        xi(1) *(1._RP-xi(2)) &
                                   - corners(j,3)*        xi(1) *       xi(2)  &
                                   + corners(j,7)*        xi(1) *       xi(2)
      END DO

      grad_x = 0.5_RP * grad_x
!
      END SUBROUTINE ComputeGradHexTransfiniteMap
!
!     ///////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!
!!     GeneralHexGradAndMap: Given the six faces and four cornerss of a 
!!     hex element, and a computational point xi, return physical space location
!!     and the jacouan matrix grad_x = d x_i/d Xi_j
!!
!!>
!!    xi           = computational space variable (X, Y, Z)
!!    grad_x = physical space gradient grad_x(:,1) = (x_Xi, x_Eta, x_Zeta), etc.
!!    cornerPoints = the 8 corners in standard fe format
!!    faceData     = interpolant data for the 6 faces
!!<
!!
!-------------------------------------------------------------------------------
!
      SUBROUTINE GeneralHexGradAndMap( u, x, grad_x, cornerPoints,faceData )
!
      USE SMConstants

      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)   , DIMENSION(3)  , INTENT(IN)  :: u
      REAL(KIND=RP)   , DIMENSION(3)  , INTENT(OUT) :: x
      REAL(KIND=RP)   , DIMENSION(3,3), INTENT(OUT) :: grad_x
      REAL(KIND=RP)   , DIMENSION(3,8), INTENT(IN)  :: cornerPoints
      TYPE(FacePatch) , DIMENSION(6)  , INTENT(IN)  :: faceData
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) , DIMENSION(3)      :: xi
      REAL(KIND=RP) , DIMENSION(3,2)    :: grad2D
      REAL(KIND=RP) , DIMENSION(3,6)    :: facePoint
      REAL(KIND=RP) , DIMENSION(3,2,6)  :: faceDer
      REAL(KIND=RP) , DIMENSION(3,12)   :: edgePoint
      REAL(KIND=RP) , DIMENSION(3,12)   :: edgeDer
      REAL(KIND=RP)                     :: eps
      INTEGER                           :: i,j
!
!     -------------
!     Compute edges
!     -------------
!
      CALL ComputeFacePoint(faceData(1),[u(1)  ,-1._RP], edgePoint(:,1))
      CALL ComputeFacePoint(faceData(1),[1.0_RP ,u(3) ], edgePoint(:,2))
      CALL ComputeFacePoint(faceData(1),[u(1)  ,1.0_RP], edgePoint(:,3))
      CALL ComputeFacePoint(faceData(1),[-1.0_RP,u(3) ], edgePoint(:,4))
      CALL ComputeFacePoint(faceData(2),[u(1)  ,-1._RP], edgePoint(:,5))
      CALL ComputeFacePoint(faceData(2),[1.0_RP ,u(3) ], edgePoint(:,6))
      CALL ComputeFacePoint(faceData(2),[u(1)  ,1.0_RP], edgePoint(:,7))
      CALL ComputeFacePoint(faceData(2),[-1.0_RP,u(3) ], edgePoint(:,8))
      CALL ComputeFacePoint(faceData(4),[u(2)  ,-1._RP], edgePoint(:,10))
      CALL ComputeFacePoint(faceData(6),[u(2)  ,-1._RP], edgePoint(:,9))
      CALL ComputeFacePoint(faceData(4),[u(2)  ,1.0_RP], edgePoint(:,11))
      CALL ComputeFacePoint(faceData(6),[u(2)  ,1.0_RP], edgePoint(:,12))
!
!     -------------
!     Compute faces
!     -------------
!
      CALL ComputeFacePoint(faceData(1), [u(1), u(3)], facePoint(:,1))
      CALL ComputeFacePoint(faceData(2), [u(1), u(3)], facePoint(:,2))
      CALL ComputeFacePoint(faceData(3), [u(1), u(2)], facePoint(:,3))
      CALL ComputeFacePoint(faceData(4), [u(2), u(3)], facePoint(:,4))
      CALL ComputeFacePoint(faceData(5), [u(1), u(2)], facePoint(:,5))
      CALL ComputeFacePoint(faceData(6), [u(2), u(3)], facePoint(:,6))
!
!     ------------------------
!     Compute edge derivatives
!     ------------------------
!
      CALL ComputeFaceDerivative(faceData(1), [u(1), -1._RP] , grad2D)
      edgeDer(:,1) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(1), [1.0_RP, u(3)], grad2D)
      edgeDer(:,2) = grad2D(:,2)
      CALL ComputeFaceDerivative(faceData(1), [u(1), 1.0_RP], grad2D)
      edgeDer(:,3) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(1), [-1.0_RP,u(3)], grad2D)
      edgeDer(:, 4) = grad2D(:,2)
      CALL ComputeFaceDerivative(faceData(2), [u(1), -1._RP] , grad2D)
      edgeDer(:, 5) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(2), [1.0_RP, u(3)], grad2D)
      edgeDer(:, 6) = grad2D(:,2)
      CALL ComputeFaceDerivative(faceData(2), [u(1), 1.0_RP], grad2D)
      edgeDer(:, 7) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(2), [-1.0_RP,u(3)], grad2D)
      edgeDer(:, 8) = grad2D(:,2)
      CALL ComputeFaceDerivative(faceData(6), [u(2), -1._RP] , grad2D)
      edgeDer(:, 9) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(4), [u(2), -1._RP] , grad2D)
      edgeDer(:,10) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(4), [u(2), 1.0_RP], grad2D)
      edgeDer(:,11) = grad2D(:,1)
      CALL ComputeFaceDerivative(faceData(6), [u(2), 1.0_RP], grad2D)
      edgeDer(:,12) = grad2D(:,1)
!
!     ------------------------
!     Compute face derivatives
!     ------------------------
!
      CALL ComputeFaceDerivative(faceData(1),[u(1),u(3)], faceDer(:,:,1))
      CALL ComputeFaceDerivative(faceData(2),[u(1),u(3)], faceDer(:,:,2))
      CALL ComputeFaceDerivative(faceData(3),[u(1),u(2)], faceDer(:,:,3))
      CALL ComputeFaceDerivative(faceData(4),[u(2),u(3)], faceDer(:,:,4))
      CALL ComputeFaceDerivative(faceData(5),[u(1),u(2)], faceDer(:,:,5))
      CALL ComputeFaceDerivative(faceData(6),[u(2),u(3)], faceDer(:,:,6))
!
!     --------------------
!     Compute the mappings
!     --------------------
!
      xi = 0.5_RP*(u + 1.0_RP)
      CALL ComputeHexTransfiniteMap(xi, x, facePoint, edgePoint, cornerPoints)
      CALL ComputeGradHexTransfiniteMap(xi, grad_x, facePoint, faceDer, &
     &                                  edgePoint, edgeDer, cornerPoints)
!
!     ---------------------------
!     Zero out rounded quantities
!     ---------------------------
!
      eps = 100._RP*EPSILON(eps)
      DO i = 1,3
         DO j = 1,3
            IF( ABS(grad_x(i, j)) <= eps )   grad_x(i,j) = 0.0_RP
         END DO
      END DO
!
!     -------------------------------------------
!     At this point should renormalize the matrix
!     -------------------------------------------
!

!
      END SUBROUTINE GeneralHexGradAndMap
      
      END MODULE TransfiniteMapClass