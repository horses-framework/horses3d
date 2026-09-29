!
! //////////////////////////////////////////////////////////////////////////////
!
!!
!!     self class stores data needed to define a 2D interpolant. In
!!     self context, self means an iterpolant that defines a surface.
!
!      TYPE FacePatch
!
!      PUBLIC METHODS:
!         SUBROUTINE ConstructFacePatch         ( self, points, uKnots, vKnots )
!         SUBROUTINE DestructFacePatch          ( self )
!         SUBROUTINE ComputeFacePoint           ( self, u, p )
!         SUBROUTINE ComputeFaceDerivative      ( self, u, grad )
!         SUBROUTINE PrintFacePatch             ( self )
!         LOGICAL FUNCTION FaceIs4CorneredQuad  ( self )
!
!!    
! //////////////////////////////////////////////////////////////////////////////
#include "Includes.h"
!
!  ******
   MODULE FacePatchClass
!  ******
!
     USE SMConstants
     use PolynomialInterpAndDerivsModule
     IMPLICIT NONE
     PRIVATE
!
!    ---------------
!    Type definition
!    ---------------
!
!!   Stores the data needed to specify a surface through a 2D
!!   interpolant. The surface is a function of (u,v) in [-1,1]x[-1,1].
!    Arrays are dimensioned 1:nKnots.
!
     TYPE FacePatch
         REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: points
         REAL(KIND=RP), DIMENSION(:)    , ALLOCATABLE :: uKnots,vKnots
         REAL(kind=RP), dimension(:),     allocatable :: wbu, wbv
         REAL(kind=RP), dimension(:,:),   allocatable :: Du, Dv
         INTEGER      , DIMENSION(2)                  :: noOfKnots

!
!        ========         
         CONTAINS 
!        ========         
!
         PROCEDURE :: construct => ConstructFacePatch
         PROCEDURE :: destruct  => DestructFacePatch
         PROCEDURE :: setFacePoints
     END TYPE FacePatch
     
     PUBLIC:: FacePatch
     PUBLIC:: ConstructFacePatch, DestructFacePatch
     PUBLIC:: ComputeFacePoint  , ComputeFaceDerivative
     PUBLIC:: ProjectFaceToNewPoints
     PUBLIC:: PrintFacePatch    , FaceIs4CorneredQuad
!
!    ========
     CONTAINS
!    ========
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    The constructor for a 2D interpolant takes an array of points and
!!    an array of knots to define a sursurface.
!     -----------------------------------------------------------------
!
      SUBROUTINE ConstructFacePatch( self, uKnots, vKnots, points )
!
!     ---------
!     Arguments
!     ---------
!
      CLASS(FacePatch)                          :: self
      REAL(KIND=RP), DIMENSION(:,:,:), OPTIONAL :: points
      REAL(KIND=RP), DIMENSION(:)               :: uKnots,vKnots
!
!     ---------------
!     Local Variables
!     ---------------
!

!
!     ---------------------------------------------------------
!     The dimensions of the interpolant are given by the number 
!     of points in each  direction. Allocate memory to store
!     self data in.
!     ---------------------------------------------------------
!
      self % noOfKnots(1) = SIZE(uKnots)
      self % noOfKnots(2) = SIZE(vKnots)
!      
      !write(*,*) 'constructing face patch (allocating points)'
      if (.NOT. ALLOCATED(self % points)) ALLOCATE( self % points(3, self % noOfKnots(1), self % noOfKnots(2)) )
      if (.NOT. ALLOCATED(self % uKnots)) ALLOCATE( self % uKnots(self % noOfKnots(1)) )
      if (.NOT. ALLOCATED(self % vKnots)) ALLOCATE( self % vKnots(self % noOfKnots(2)) )
      if (.NOT. ALLOCATED(self % wbu))    allocate( self % wbu(self % noOfKnots(1)) )
      if (.NOT. ALLOCATED(self % wbv))    allocate( self % wbv(self % noOfKnots(2)) )
      if (.NOT. ALLOCATED(self % Du))     allocate( self % Du(self % noOfKnots(1), self % noOfKnots(1)) )
      if (.NOT. ALLOCATED(self % Dv))     allocate( self % Dv(self % noOfKnots(2), self % noOfKnots(2)) )
!
!     -------------------------
!     Save the points and knots
!     -------------------------
!
      self % uKnots  = uKnots
      self % vKnots  = vKnots
!
!     Compute the barycentric weights
!     -------------------------------
      call BarycentricWeights(self % noOfKnots(1)-1, self % uKnots, self % wbu)
      call BarycentricWeights(self % noOfKnots(2)-1, self % vKnots, self % wbv)
!
!     Compute the derivative matrices
!     -------------------------------
      call PolynomialDerivativeMatrix(self % noOfKnots(1)-1, self % uKnots, self % Du)
      call PolynomialDerivativeMatrix(self % noOfKnots(2)-1, self % vKnots, self % Dv)

      IF(PRESENT(points)) CALL self % setFacePoints( points )
!
      END SUBROUTINE ConstructFacePatch
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    Deletes allocated memory
!     -----------------------------------------------------------------
!
      elemental SUBROUTINE DestructFacePatch(self)
!
!     ---------
!     Arguments
!     ---------
!
      CLASS(FacePatch), intent(inout) :: self
      
      IF ( ALLOCATED( self % points ) )   DEALLOCATE( self % points )
      IF ( ALLOCATED( self % uKnots ) )   DEALLOCATE( self % uKnots )
      IF ( ALLOCATED( self % vKnots ) )   DEALLOCATE( self % vKnots )
      if (allocated(self%wbu)) deallocate(self%wbu)
      if (allocated(self%wbv)) deallocate(self%wbv)
      if (allocated(self%Du))  deallocate(self%Du)
      if (allocated(self%Dv))  deallocate(self%Dv)
      
      self % noOfKnots = 0
!
      END SUBROUTINE DestructFacePatch
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE setFacePoints(self,points)  
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(FacePatch)                :: self
         REAL(KIND=RP), DIMENSION(:,:,:) :: points
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: n, k, j
!
         self % points  = points
         
      END SUBROUTINE setFacePoints
!
!     //////////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------
!
!!     ComputeFacePoint: Compute the nodes on a surface
!!     by interpolating
!!
!!>
!!     u       = computational space variable (u, Zeta)
!!     p        = resultant physical space variable (x, y, z)
!!     self     = Interpolation data that defines a surface
!!<
!!
!     --------------------------------------------------------------------------
!
      SUBROUTINE ComputeFacePoint(self, u, p)
      IMPLICIT NONE
!     ---------
!     Arguments
!     ---------
!
      TYPE(FacePatch), INTENT(IN)                 :: self
      REAL(KIND=RP)  , DIMENSION(2), INTENT(IN)   :: u
      REAL(KIND=RP)  , DIMENSION(3), INTENT(OUT)  :: p
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER :: j
!
!     ------------------------------------------------------------------
!     Use bi-linear mapping if self is a flat side (for speed) otherwise
!     do general Lagrange interpolation
!     ------------------------------------------------------------------
!
      IF( FaceIs4CorneredQuad(self))     THEN
         DO j = 1,3
            p(j)  = self % points(j,1,1)*(1._RP - u(1))*(1._RP - u(2)) &
                  + self % points(j,2,1)*(1._RP + u(1))*(1._RP - u(2)) &
                  + self % points(j,2,2)*(1._RP + u(1))*(1._RP + u(2))  &
                  + self % points(j,1,2)*(1._RP - u(1))*(1._RP + u(2))
         END DO
         p = 0.25_RP * p
      ELSE
         CALL ComputePoly2D(self, u, p)
      END IF
      
      RETURN
      END SUBROUTINE ComputeFacePoint
!
!///////////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------
!
!!     ComputePoly2D: Compute the lagrange interpolant through the points 
!!     p(i, j) at the location (x, y)
!!
!!>
!!     u   = computational space variable (u, Zeta)
!!     p    = resultant physical space variable (x, y, z)
!!     self = Interpolation data that defines a surface
!!<
!!
!     --------------------------------------------------------------------------
!
      SUBROUTINE ComputePoly2D(self, u, p)
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FacePatch)             :: self
      REAL(KIND=RP), DIMENSION(2) :: u
      REAL(KIND=RP), DIMENSION(3) :: p
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP)                              :: l_j(self % noOfKnots(2))
      real(kind=RP)                              :: l_i(self % noOfKnots(1))
      INTEGER                                    :: i, j

      call InterpolatingPolynomialVector(u(1), self % noOfKnots(1)-1, self % uKnots, self % wbu, l_i)
      call InterpolatingPolynomialVector(u(2), self % noOfKnots(2)-1, self % vKnots, self % wbv, l_j)

      p = 0.0_RP

      DO j = 1, self % noOfKnots(2)
         do i = 1, self % noOfKnots(1)
            p = p + self % points(:,i,j) * l_i(i) * l_j(j)
         end do
      end do

      END SUBROUTINE ComputePoly2D
!
!     ///////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------
!
!!     ProjectFaceToNewPoints
!!
!     --------------------------------------------------------------------------
!
      subroutine ProjectFaceToNewPoints(patch,Nx,xi,Ny,eta,facecoords)
         
         implicit none
         type(FacePatch),  intent(in)     :: patch
         integer,          intent(in)     :: Nx
         real(kind=RP),    intent(in)     :: xi(Nx+1)
         integer,          intent(in)     :: Ny
         real(kind=RP),    intent(in)     :: eta(Ny+1)
         real(kind=RP),    intent(out)    :: faceCoords(1:3,Nx+1,Ny+1)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j
         real(kind=RP)  :: localCoords(2)
               
         do j = 1, Ny+1 ; do i = 1, Nx+1
            localCoords = (/ xi(i), eta(j) /)
            call ComputeFacePoint(patch, localCoords, faceCoords(:,i,j) )
         end do       ; end do

      end subroutine ProjectFaceToNewPoints
!
!     ///////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------
!
!!     ComputeSurfaceDerivative: Compute the derivative in the two local  
!!     coordinate directions on a subdomain surface by interpolating
!!     the surface data
!!
!     --------------------------------------------------------------------------
!
      SUBROUTINE ComputeFaceDerivative(self, u, grad)
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FacePatch)       :: self
      REAL(KIND=RP) , DIMENSION(2)   :: u
      REAL(KIND=RP) , DIMENSION(3,2) :: grad
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER :: j
!
!     ------------------------------------------------------------------
!     Use bi-linear mapping of self is a flat side (for speed) otherwise
!     do general lagrange interpolation
!     ------------------------------------------------------------------
!
      IF( FaceIs4CorneredQuad(self))     THEN
         DO j = 1,3
            grad(j,1)  =   -self % points(j,1,1)*(1._RP - u(2)) &
                          + self % points(j,2,1)*(1._RP - u(2)) &
                          + self % points(j,2,2)*(1._RP + u(2)) &
                          - self % points(j,1,2)*(1._RP + u(2))

            grad(j,2)  =  - self % points(j,1,1)*(1._RP - u(1)) &
                          - self % points(j,2,1)*(1._RP + u(1))      &
                          + self % points(j,2,2)*(1._RP + u(1))      &
                          + self % points(j,1,2)*(1._RP - u(1))
         END DO
         grad = 0.25_RP * grad
      ELSE
         CALL Compute2DPolyDeriv(self, u, grad)
      END IF

      RETURN
      END SUBROUTINE ComputeFaceDerivative
!
!     ///////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------
!
!!     Compute2DPolyDeriv: Compute the gradient of the 
!!                         lagrange interpolant through the points 
!!                         p(i, j) at the location (u, v). The result, grad,
!!                         is (dX/du,dX/dv) where X = (x,y,z)
!!
!     --------------------------------------------------------------------------
!
      SUBROUTINE Compute2DPolyDeriv(self, u, grad)
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FacePatch)      , INTENT(IN)  :: self
      REAL(KIND=RP), DIMENSION(2)   , INTENT(IN)  :: u
      REAL(KIND=RP), DIMENSION(3,2) , INTENT(OUT) :: grad
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: l_i(self % noOfKnots(1))
      REAL(KIND=RP) :: l_j(self % noOfKnots(2))
      REAL(KIND=RP) :: dl_i(self % noOfKnots(1))
      REAL(KIND=RP) :: dl_j(self % noOfKnots(2))
      INTEGER                                  :: i, j, k
!
      grad = 0.0_RP
!
!     ----------------------------
!     Compute Lagrange polynomials
!     ----------------------------
!
      call InterpolatingPolynomialVector(u(1), self % noOfKnots(1)-1, self % uKnots, self % wbu, l_i)
      call InterpolatingPolynomialVector(u(2), self % noOfKnots(2)-1, self % vKnots, self % wbv, l_j)

      dl_i = 0.0_RP
      do i = 1, self % noOfKnots(1)
         do j = 1, self % noOfKnots(1)
            dl_i(i) = dl_i(i) + self % Du(j,i) * l_i(j)
         end do   
      end do

      dl_j = 0.0_RP
      do i = 1, self % noOfKnots(2)
         do j = 1, self % noOfKnots(2)
            dl_j(i) = dl_j(i) + self % Dv(j,i) * l_j(j)
         end do   
      end do
!
!     ------------------------
!     Evaluate the polynomials
!     ------------------------
!
      do j = 1, self % noOfKnots(2)
         do i = 1, self % noOfKnots(1)
            grad(:,1) = grad(:,1) + self % points(:,i,j) * dl_i(i) * l_j(j)
            grad(:,2) = grad(:,2) + self % points(:,i,j) *  l_i(i) * dl_j(j)
         end do
      end do

      END SUBROUTINE Compute2DPolyDeriv
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    For debugging purposes, print the information in the surface
!!    interpolant.
!     -----------------------------------------------------------------
!
      SUBROUTINE PrintFacePatch( self )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FacePatch)        :: self
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER :: i,j
!
      PRINT *, "-------------Surface Interpolant--------------"
      DO j = 1, self % noOfKnots(2)
         DO i = 1, self % noOfKnots(1)
            PRINT *, i,j, self % uKnots(i), self % vKnots(j), self % points(:,i,j)
         END DO
      END DO
!
      END SUBROUTINE PrintFacePatch
!
!////////////////////////////////////////////////////////////////////////
!
!
!-------------------------------------------------------------------------
! FaceIs4CorneredQuad returns .TRUE. if the face is represented only by
! the four corners.
!-------------------------------------------------------------------------
!
      LOGICAL FUNCTION FaceIs4CorneredQuad( self )
      IMPLICIT NONE
!
!      ---------
!      Arguments
!      ---------
!
       TYPE(FacePatch)        :: self
!
!      ---------------
!      Local Variables
!      ---------------
!
      IF ( (self % noOfKnots(1) == 2) .AND. (self % noOfKnots(2) == 2) )     THEN
         FaceIs4CorneredQuad = .TRUE.
      ELSE
         FaceIs4CorneredQuad = .FALSE.
      END IF 
      
      END FUNCTION FaceIs4CorneredQuad
!@mark -
!
! /////////////////////////////////////////////////////////////////////
!
!!    Compute the newton divided difference table.
!!    The node values are destroyed while creating the table.
!
! /////////////////////////////////////////////////////////////////////
!
      subroutine divdif( nc, x, table )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                          :: nc
      REAL(KIND=RP), DIMENSION(nc) :: x, table
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                          :: i, j
!
      do 100 i = 2,nc
         do 100 j = nc,i,-1
            table(j) = (table(j) - table(j-1))/(x(j) - x(j-i+1))
 100  continue
!
      return
      END subroutine divdif
!
! /////////////////////////////////////////////////////////////////////
!
!!    Compute the newton form interpolant at the point x
!
! /////////////////////////////////////////////////////////////////////
!
      FUNCTION EvaluateNewtonPolynomial( x, nc, xvals, table )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                          :: nc
      REAL(KIND=RP), DIMENSION(nc) :: xVals, table
      REAL(KIND=RP)                :: x
      REAL(KIND=RP)                :: EvaluateNewtonPolynomial
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP)                :: w, p
      INTEGER                          :: k
      
      p = table(1)
      w = 1.0_RP
      do 100 k = 2,nc
         w = (x - xVals(k-1))*w
         p = p + table(k)*w
 100  continue
      EvaluateNewtonPolynomial = p
!
      END FUNCTION EvaluateNewtonPolynomial
!
! /////////////////////////////////////////////////////////////////////
!
!!    Compute the derivative of the newton form interpolant at 
!!    the point x
!
! /////////////////////////////////////////////////////////////////////
!
      FUNCTION EvaluateNewtonPolynomialDeriv( x, nc, xvals, table )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                      :: nc
      REAL(KIND=RP), DIMENSION(nc) :: xVals, table
      REAL(KIND=RP)                :: x
      REAL(KIND=RP)                :: EvaluateNewtonPolynomialDeriv
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(nc) :: s
      REAL(KIND=RP)                :: w, pp, z
      INTEGER                          :: i,j
      
      pp = table(2)
      EvaluateNewtonPolynomialDeriv = pp
      IF( nc == 2 )     RETURN
      s(1) = x - xVals(1)
      s(2) = x - xVals(2)
      pp = pp + (s(1) + s(2))*table(3)
      EvaluateNewtonPolynomialDeriv = pp
      IF( nc == 3 )     RETURN
      w = s(1)*s(2)
      
      DO i = 4,nc
         z = 0.0_RP
         DO j = 1, i-1
            s(j) = s(j)*(x - xVals(i-1))
            z = z + s(j)
         END DO
         s(i) = w
         z    = z + w
         w    = w*(x - xVals(i-1))
         pp   = pp + z*table(i)
      END DO
      
      EvaluateNewtonPolynomialDeriv = pp
!
      END FUNCTION EvaluateNewtonPolynomialDeriv
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!!    Compute at x the lagrange polynomial L_k of degree n-1    
!!    whose zeros are given by the z(i) 
!---------------------------------------------------------------------
!
      FUNCTION EvaluateLagrangePoly(k,x,n,z)
!
      IMPLICIT NONE
      REAL(KIND=RP) :: EvaluateLagrangePoly
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER, INTENT(IN)                         :: k !! Which poly
      INTEGER, INTENT(IN)                         :: n !! Order+1
      REAL(KIND=RP), DIMENSION(n), INTENT(IN) :: z !! Nodes
      REAL(KIND=RP)              , INTENT(IN) :: x !! Eval Pt
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER           :: j
      REAL(KIND=RP) :: eLP
      
!                                                                       
      IF(k == 1)     THEN 
         eLP = (x - z(2))/(z(k) - z(2)) 
         DO  j = 3,n 
            eLP = eLP*(x - z(j))/(z(k) - z(j))
         END DO
      ELSE 
         eLP = (x - z(1))/(z(k) - z(1)) 
         DO  j = 2,k-1 
            eLP = eLP*(x - z(j))/(z(k) - z(j))
         END DO
         DO j = k+1,n 
            eLP = eLP*(x - z(j))/(z(k) - z(j))
         END DO
      END IF
      EvaluateLagrangePoly= eLp
!                                                                       
      END FUNCTION EvaluateLagrangePoly
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!!    Compute at x the derivative of the lagrange polynomial L_k of
!!    degree n-1 whose zeros are given by the z(i)
!---------------------------------------------------------------------
!
      FUNCTION EvaluateLagrangePolyDeriv(k,x,n,z)
!
      IMPLICIT NONE
      REAL(KIND=RP) :: EvaluateLagrangePolyDeriv
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER, INTENT(IN)                         :: k !! Which poly
      INTEGER, INTENT(IN)                         :: n !! Order+1
      REAL(KIND=RP), DIMENSION(n), INTENT(IN) :: z !! Nodes
      REAL(KIND=RP)              , INTENT(IN) :: x !! Eval Pt
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER           :: l,m
      REAL(KIND=RP) :: hp,poly
!                                                                       
      hp = 0.0_RP
      DO l = 1,n 
         if(l == k)     CYCLE
         poly = 1.0_RP
         DO m = 1,n 
            if(m == l)     CYCLE
            if(m == k)     CYCLE 
            poly = poly*(x - z(m))/(z(k) - z(m))
         END DO
         hp = hp + poly/(z(k) - z(l)) 
      END DO
      EvaluateLagrangePolyDeriv = hp 
!                                                                       
      END FUNCTION EvaluateLagrangePolyDeriv

!
!
! /////////////////////////////////////////////////////////////////////
!
!  **********     
   END MODULE FacePatchClass
!  **********     