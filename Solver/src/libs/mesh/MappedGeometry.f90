#include "Includes.h"
Module MappedGeometryClass
   USE SMConstants
   USE TransfiniteMapClass
   USE NodalStorageClass
   use MeshTypes
   IMPLICIT NONE

      private
      public MappedGeometry, MappedGeometryFace
      public SetMappingsToCrossProduct
      public vDot, vCross, ComputeVorticity
!
!     ---------
!     Constants
!     ---------
!
      LOGICAL  :: useCrossProductMetrics = .false.
!
!     -----
!     Class
!     -----
!
      TYPE MappedGeometry
            INTEGER                                         :: Nx, Ny, Nz                    ! Polynomial order
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: jGradXi, jGradEta, jGradZeta  ! Contravariant vectors times Jacobian!
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: x                             ! Position of points in absolute coordinates
            REAL(KIND=RP), DIMENSION(:,:,:)   , ALLOCATABLE :: jacobian, invJacobian         ! Mapping Jacobian and 1/Jacobian
            real(kind=RP)                                   :: volume
            real(kind=RP), dimension(:,:,:),    allocatable :: dWall                   ! Minimum distance to the nearest wall
            real(kind=RP), dimension(:,:,:,:),  allocatable :: normal                  ! Wall normal, needed for IB
            real(kind=RP), dimension(:,:,:,:) , allocatable :: ncXi, ncEta, ncZeta     ! Normals at the complementary grid nodes
            real(kind=RP), dimension(:,:,:,:) , allocatable :: t1cXi, t1cEta, t1cZeta  ! Tangent vector 1 at the complementary grid nodes
            real(kind=RP), dimension(:,:,:,:) , allocatable :: t2cXi, t2cEta, t2cZeta  ! Tangent vector 2 at the complementary grid nodes
            real(kind=RP), dimension(:,:,:)   , allocatable :: JfcXi, JfcEta, JfcZeta  ! Normalization term of the faces of the complementary grid
            CONTAINS

            PROCEDURE :: construct               => ConstructMappedGeometry
            PROCEDURE :: destruct                => DestructMappedGeometry
      END TYPE MappedGeometry

      type MappedGeometryFace
         real(kind=RP), dimension(:,:,:), allocatable   :: x
         real(kind=RP), dimension(:,:)  , allocatable   :: jacobian   ! |ja^i|: Normalization term of the normal vectors on a face
         real(kind=RP), dimension(:,:,:), allocatable   :: GradXi, GradEta, GradZeta  ! Contravariant vectors
         real(kind=RP), dimension(:,:,:), allocatable   :: normal     ! normal vector on a face
         real(kind=RP), dimension(:,:,:), allocatable   :: t1         ! Tangent vector (along the xi direction)
         real(kind=RP), dimension(:,:,:), allocatable   :: t2         ! Tangent vector 2 (orthonormal to t1 and normal)
         real(kind=RP), dimension(:,:),   allocatable   :: dWall      ! Minimum distance to the nearest wall
         real(kind=RP)                                  :: surface    ! Surface
         real(kind=RP)                                  :: h          ! Element dimension orthogonal to the face
         contains
            procedure :: construct => ConstructMappedGeometryFace
            procedure :: destruct  => DestructMappedGeometryFace
      end type MappedGeometryFace

!
!  ========
   CONTAINS
!  ========
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ConstructMappedGeometry( self, spAxi, spAeta, spAzeta, mapper )
      IMPLICIT NONE
!
!      ---------
!      Arguments
!      ---------
!
      CLASS(MappedGeometry)  , intent(inout) :: self
      TYPE(TransfiniteHexMap), intent(in)    :: mapper
      TYPE(NodalStorage_t)     , intent(in)    :: spAxi
      TYPE(NodalStorage_t)     , intent(in)    :: spAeta
      TYPE(NodalStorage_t)     , intent(in)    :: spAzeta
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER       :: Nx, Ny, Nz, Nmax
      INTEGER       :: i, j, k
      REAL(KIND=RP) :: nrm
      REAL(KIND=RP) :: grad_x(3,3), jGrad(3), x(3)
!
!     -----------
!     Allocations
!     -----------
!
      Nx        = spAxi   % N
      Ny        = spAeta  % N
      Nz        = spAzeta % N
      Nmax      = MAX(Nx,Ny,Nz)
      self % Nx = Nx
      self % Ny = Ny
      self % Nz = Nz

      if (.not.allocated(self % JGradXi)) ALLOCATE( self % JGradXi  (3,0:Nx,0:Ny,0:Nz) )
      if (.not.allocated(self % JGradEta)) ALLOCATE( self % JGradEta (3,0:Nx,0:Ny,0:Nz) )
      if (.not.allocated(self % JGradZeta)) ALLOCATE( self % JGradZeta(3,0:Nx,0:Ny,0:Nz) )
      if (.not.allocated(self % jacobian)) ALLOCATE( self % jacobian   (0:Nx,0:Ny,0:Nz) )
      if (.not.allocated(self % invJacobian)) ALLOCATE( self % invJacobian(0:Nx,0:Ny,0:Nz) )
      if (.not.allocated(self % x)) ALLOCATE( self % x        (3,0:Nx,0:Ny,0:Nz)    )

      if (allocated(self % JGradXi))  self % JGradXi = 0.0_RP 
      if (allocated(self % JGradEta))  self % JGradEta = 0.0_RP 
      if (allocated(self % JGradZeta)) self % JGradZeta= 0.0_RP 
      if (allocated(self % jacobian)) self % jacobian   = 0.0_RP 
      if (allocated(self % invJacobian)) self % invJacobian= 0.0_RP 
      if (allocated(self % x))  self % x = 0.0_RP 
!
!     --------------------------
!     Compute interior locations
!     --------------------------
!
      !write(*,*)'constructing geometry'
     !write(*,*)'spAxi % x', spAxi % x
      !(*,*)'spAeta % x', spAeta % x
      !write(*,*)'spAzeta % x',spAzeta % x
      DO k = 0, Nz
         DO j= 0, Ny
            DO i = 0,Nx
               x = [spAxi % x(i), spAeta % x(j), spAzeta % x(k)]
               !write(*,*)'x befor transfine map', x
               self % x(:,i,j,k) = mapper %  transfiniteMapAt(x)
               !write(*,*)'self % x(:,i,j,k) (after transfine map)', self % x(:,i,j,k)
            END DO
         END DO
      END DO
!
!     ------------
!     Metric terms
!     ------------
!
      IF ( useCrossProductMetrics ) THEN

         CALL computeMetricTermsCrossProductForm(self, spAxi, spAeta, spAzeta, mapper)

      ELSE

         CALL computeMetricTermsConservativeForm(self, spAxi, spAeta, spAzeta, mapper)

      ENDIF
!
!     -----------
!     Cell volume
!     -----------
!
      self % volume = 0.0_RP
      do k = 0, Nz   ; do j = 0, Ny ; do i = 0, Nx
         self % volume = self % volume + spAxi % w(i) * spAeta % w(j) * spAzeta % w(k) * self % jacobian(i,j,k)
      end do         ; end do       ; end do

   END SUBROUTINE ConstructMappedGeometry
!
!////////////////////////////////////////////////////////////////////////
!
      pure SUBROUTINE DestructMappedGeometry(self)
         IMPLICIT NONE
         CLASS(MappedGeometry), intent(inout) :: self

         safedeallocate( self % jGradXi     )
         safedeallocate( self % jGradEta    )
         safedeallocate( self % jGradZeta   )
         safedeallocate( self % jacobian    )
         safedeallocate( self % invJacobian )
         safedeallocate( self % x           )
         safedeallocate( self % dWall       )
         safedeallocate( self % ncXi        )
         safedeallocate( self % ncEta       )
         safedeallocate( self % ncZeta      )
         safedeallocate( self % t1cXi       )
         safedeallocate( self % t1cEta      )
         safedeallocate( self % t1cZeta     )
         safedeallocate( self % t2cXi       )
         safedeallocate( self % t2cEta      )
         safedeallocate( self % t2cZeta     )
         safedeallocate( self % JfcXi       )
         safedeallocate( self % JfcEta      )
         safedeallocate( self % JfcZeta     )
         safedeallocate( self % normal      )

      END SUBROUTINE DestructMappedGeometry
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine computeMetricTermsConservativeForm(self, spAxi, spAeta, spAzeta, mapper)
!
!        *********************************************************************
!              Currently, the invariant form is implemented
!
!              Ja^i_n = -1/2 \hat{x}^i ( Xl \nabla Xm - Xm \nabla Xl )
!                 (i,j,k) and (n,m,l) cyclic
!        *********************************************************************
!
         use PhysicsStorage
         implicit none
         type(MappedGeometry),    intent(inout) :: self
         type(NodalStorage_t),      intent(in)    :: spAxi, spAeta, spAzeta
         type(TransfiniteHexMap), intent(in)    :: mapper
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, m, n, l
         real(kind=RP)  :: grad_x(NDIM,NDIM,0:self % Nx, 0:self % Ny, 0:self % Nz)
         real(kind=RP)  :: xCGL(NDIM,0:self % Nx, 0:self % Ny, 0:self % Nz)
         real(kind=RP)  :: auxgrad(NDIM,NDIM,0:self % Nx, 0:self % Ny, 0:self % Nz)
         real(kind=RP)  :: coordsProduct(NDIM,0:self % Nx,0:self % Ny,0:self % Nz)
         real(kind=RP)  :: Jai(NDIM,0:self % Nx, 0:self % Ny, 0:self % Nz)
         real(kind=RP)  :: Ja1CGL(NDIM,0:self % Nx, 0:self % Ny, 0:self % Nz)
         real(kind=RP)  :: Ja2CGL(NDIM,0:self % Nx, 0:self % Ny, 0:self % Nz)
         real(kind=RP)  :: Ja3CGL(NDIM,0:self % Nx, 0:self % Ny, 0:self % Nz)
         real(kind=RP)  :: JacobianCGL(0:self % Nx, 0:self % Ny, 0:self % Nz)
         real(kind=RP)  :: x(3)
!
!        Compute the mapping gradient in Chebyshev-Gauss-Lobatto points
!        --------------------------------------------------------------
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            x = [spAxi % xCGL(i), spAeta % xCGL(j), spAzeta % xCGL(k)]
            xCGL(:,i,j,k) = mapper % transfiniteMapAt(x)
            grad_x(:,:,i,j,k) = mapper % metricDerivativesAt(x)
         end do         ; end do          ; end do
!
!        *****************************************
!        Compute the x-coordinates of the mappings
!        *****************************************
!
!        Compute coordinates combination
!        -------------------------------
         do k = 0, self % Nz    ; do j = 0, self % Ny  ; do i = 0, self % Nx
            coordsProduct(:,i,j,k) =   xCGL(3,i,j,k) * grad_x(2,:,i,j,k)  &
                                     - xCGL(2,i,j,k) * grad_x(3,:,i,j,k)
         end do            ; end do          ; end do
!
!        Compute its gradient
!        --------------------
         auxgrad = 0.0_RP
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            do l = 0, self % Nx
               auxgrad(:,1,i,j,k) = auxgrad(:,1,i,j,k) + coordsProduct(:,l,j,k) * spAxi % DCGL(i,l)
            end do

            do l = 0, self % Ny
               auxgrad(:,2,i,j,k) = auxgrad(:,2,i,j,k) + coordsProduct(:,i,l,k) * spAeta % DCGL(j,l)
            end do

            do l = 0, self % Nz
               auxgrad(:,3,i,j,k) = auxgrad(:,3,i,j,k) + coordsProduct(:,i,j,l) * spAzeta % DCGL(k,l)
            end do
         end do         ; end do          ; end do
!
!        Compute the curl
!        ----------------
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            Jai(1,i,j,k) = auxgrad(3,2,i,j,k) - auxgrad(2,3,i,j,k)
            Jai(2,i,j,k) = auxgrad(1,3,i,j,k) - auxgrad(3,1,i,j,k)
            Jai(3,i,j,k) = auxgrad(2,1,i,j,k) - auxgrad(1,2,i,j,k)
         end do         ; end do          ; end do
!
!        Assign to the first coordinate of each metrics
!        ----------------------------------------------
         Ja1CGL(1,:,:,:)  = -0.5_RP * Jai(1,:,:,:)
         Ja2CGL(1,:,:,:)  = -0.5_RP * Jai(2,:,:,:)
         Ja3CGL(1,:,:,:)  = -0.5_RP * Jai(3,:,:,:)
!
!        *****************************************
!        Compute the y-coordinates of the mappings
!        *****************************************
!
!        Compute coordinates combination
!        -------------------------------
         do k = 0, self % Nz    ; do j = 0, self % Ny  ; do i = 0, self % Nx
            coordsProduct(:,i,j,k) =   xCGL(1,i,j,k) * grad_x(3,:,i,j,k) &
                                     - xCGL(3,i,j,k) * grad_x(1,:,i,j,k)
         end do            ; end do          ; end do
!
!        Compute its gradient
!        --------------------
         auxgrad = 0.0_RP
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            do l = 0, self % Nx
               auxgrad(:,1,i,j,k) = auxgrad(:,1,i,j,k) + coordsProduct(:,l,j,k) * spAxi % DCGL(i,l)
            end do

            do l = 0, self % Ny
               auxgrad(:,2,i,j,k) = auxgrad(:,2,i,j,k) + coordsProduct(:,i,l,k) * spAeta % DCGL(j,l)
            end do

            do l = 0, self % Nz
               auxgrad(:,3,i,j,k) = auxgrad(:,3,i,j,k) + coordsProduct(:,i,j,l) * spAzeta % DCGL(k,l)
            end do
         end do         ; end do          ; end do
!
!        Compute the curl
!        ----------------
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            Jai(1,i,j,k) = auxgrad(3,2,i,j,k) - auxgrad(2,3,i,j,k)
            Jai(2,i,j,k) = auxgrad(1,3,i,j,k) - auxgrad(3,1,i,j,k)
            Jai(3,i,j,k) = auxgrad(2,1,i,j,k) - auxgrad(1,2,i,j,k)
         end do         ; end do          ; end do
!
!        Assign to the second coordinate of each metrics
!        -----------------------------------------------
         Ja1CGL(2,:,:,:)  = -0.5_RP*Jai(1,:,:,:)
         Ja2CGL(2,:,:,:)  = -0.5_RP*Jai(2,:,:,:)
         Ja3CGL(2,:,:,:)  = -0.5_RP*Jai(3,:,:,:)
!
!        *****************************************
!        Compute the z-coordinates of the mappings
!        *****************************************
!
!        Compute coordinates combination
!        -------------------------------
         do k = 0, self % Nz    ; do j = 0, self % Ny  ; do i = 0, self % Nx
            coordsProduct(:,i,j,k) =   xCGL(2,i,j,k) * grad_x(1,:,i,j,k) &
                                     - xCGL(1,i,j,k) * grad_x(2,:,i,j,k)
         end do            ; end do          ; end do
!
!        Compute its gradient
!        --------------------
         auxgrad = 0.0_RP
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            do l = 0, self % Nx
               auxgrad(:,1,i,j,k) = auxgrad(:,1,i,j,k) + coordsProduct(:,l,j,k) * spAxi % DCGL(i,l)
            end do

            do l = 0, self % Ny
               auxgrad(:,2,i,j,k) = auxgrad(:,2,i,j,k) + coordsProduct(:,i,l,k) * spAeta % DCGL(j,l)
            end do

            do l = 0, self % Nz
               auxgrad(:,3,i,j,k) = auxgrad(:,3,i,j,k) + coordsProduct(:,i,j,l) * spAzeta % DCGL(k,l)
            end do
         end do         ; end do          ; end do
!
!        Compute the curl
!        ----------------
         do k = 0, self % Nz ; do j = 0, self % Ny  ; do i = 0, self % Nx
            Jai(1,i,j,k) = auxgrad(3,2,i,j,k) - auxgrad(2,3,i,j,k)
            Jai(2,i,j,k) = auxgrad(1,3,i,j,k) - auxgrad(3,1,i,j,k)
            Jai(3,i,j,k) = auxgrad(2,1,i,j,k) - auxgrad(1,2,i,j,k)
         end do         ; end do          ; end do
!
!        Assign to the third coordinate of each metrics
!        ----------------------------------------------
         Ja1CGL(3,:,:,:)  = -0.5_RP * Jai(1,:,:,:)
         Ja2CGL(3,:,:,:)  = -0.5_RP * Jai(2,:,:,:)
         Ja3CGL(3,:,:,:)  = -0.5_RP * Jai(3,:,:,:)
!
!        ********************
!        Compute the Jacobian
!        ********************
!
         do k = 0, self % Nz  ; do j = 0, self % Ny   ; do i = 0, self % Nx
            JacobianCGL(i,j,k) = jacobian3D(a1 = grad_x(:,1,i,j,k), &
                                                a2 = grad_x(:,2,i,j,k), &
                                                a3 = grad_x(:,3,i,j,k)   )
         end do               ; end do                ; end do
!
!        **********************
!        Return to Gauss points
!        **********************
!
         self % jGradXi = 0.0_RP
         self % jGradEta = 0.0_RP
         self % jGradZeta = 0.0_RP
         self % jacobian = 0.0_RP

         do k = 0, self % Nz  ; do j = 0, self % Ny  ; do i = 0, self % Nx
            do n = 0, self % Nz ; do m = 0, self % Ny ; do l = 0, self % Nx
               self % jGradXi(:,i,j,k) = self % jGradXi(:,i,j,k) + Ja1CGL(:,l,m,n) &
                                          * spAxi % TCheb2Gauss(i,l) &
                                          * spAeta % TCheb2Gauss(j,m) &
                                          * spAzeta % TCheb2Gauss(k,n)

               self % jGradEta(:,i,j,k) = self % jGradEta(:,i,j,k) + Ja2CGL(:,l,m,n) &
                                          * spAxi % TCheb2Gauss(i,l) &
                                          * spAeta % TCheb2Gauss(j,m) &
                                          * spAzeta % TCheb2Gauss(k,n)

               self % jGradZeta(:,i,j,k) = self % jGradZeta(:,i,j,k) + Ja3CGL(:,l,m,n) &
                                          * spAxi % TCheb2Gauss(i,l) &
                                          * spAeta % TCheb2Gauss(j,m) &
                                          * spAzeta % TCheb2Gauss(k,n)

               self % jacobian(i,j,k) = self % jacobian(i,j,k) + JacobianCGL(l,m,n) &
                                          * spAxi % TCheb2Gauss(i,l) &
                                          * spAeta % TCheb2Gauss(j,m) &
                                          * spAzeta % TCheb2Gauss(k,n)
            end do              ; end do              ; end do
         end do               ; end do               ; end do
         do k = 0, self % Nz  ; do j = 0, self % Ny  ; do i = 0, self % Nx
            self % invJacobian(i,j,k) = 1._RP / self % jacobian(i,j,k)
         end do               ; end do               ; end do
      end subroutine computeMetricTermsConservativeForm
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE computeMetricTermsCrossProductForm(self, spAxi, spAeta, spAzeta, mapper)
!
!     -----------------------------------------------
!     Compute the metric terms in cross product form
!     -----------------------------------------------
!
         use PhysicsStorage
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(MappedGeometry)   , intent(inout) :: self
         TYPE(NodalStorage_t)     , intent(in)    :: spAxi
         TYPE(NodalStorage_t)     , intent(in)    :: spAeta
         TYPE(NodalStorage_t)     , intent(in)    :: spAzeta
         TYPE(TransfiniteHexMap), intent(in)    :: mapper
!
!        ---------------
!        Local Variables
!        ---------------
!
         INTEGER       :: i,j,k,l,m,n
         INTEGER       :: Nx, Ny, Nz
         REAL(KIND=RP) :: grad_x(3,3)
         real(kind=RP)  :: Ja1CGL(NDIM,0:self % Nx, 0:self % Ny, 0:self % Nz)
         real(kind=RP)  :: Ja2CGL(NDIM,0:self % Nx, 0:self % Ny, 0:self % Nz)
         real(kind=RP)  :: Ja3CGL(NDIM,0:self % Nx, 0:self % Ny, 0:self % Nz)
         real(kind=RP)  :: JacobianCGL(0:self % Nx, 0:self % Ny, 0:self % Nz)

         Nx = spAxi % N
         Ny = spAeta % N
         Nz = spAzeta % N

         DO k = 0, Nz
            DO j = 0,Ny
               DO i = 0,Nx
                  grad_x = mapper % metricDerivativesAt([spAxi % xCGL(i), spAeta % xCGL(j), &
                                                              spAzeta % xCGL(k)])

                  CALL vCross( grad_x(:,2), grad_x(:,3), Ja1CGL (:,i,j,k))
                  CALL vCross( grad_x(:,3), grad_x(:,1), Ja2CGL (:,i,j,k))
                  CALL vCross( grad_x(:,1), grad_x(:,2), Ja3CGL(:,i,j,k))

                  JacobianCGL(i,j,k) = jacobian3D(a1 = grad_x(:,1),a2 = grad_x(:,2),a3 = grad_x(:,3))
               END DO
            END DO
         END DO
!
!        **********************
!        Return to Gauss points
!        **********************
!
         self % jGradXi = 0.0_RP
         self % jGradEta = 0.0_RP
         self % jGradZeta = 0.0_RP
         self % jacobian = 0.0_RP

         do k = 0, self % Nz  ; do j = 0, self % Ny  ; do i = 0, self % Nx
            do n = 0, self % Nz ; do m = 0, self % Ny ; do l = 0, self % Nx
               self % jGradXi(:,i,j,k) = self % jGradXi(:,i,j,k) + Ja1CGL(:,l,m,n) &
                                          * spAxi % TCheb2Gauss(i,l) &
                                          * spAeta % TCheb2Gauss(j,m) &
                                          * spAzeta % TCheb2Gauss(k,n)

               self % jGradEta(:,i,j,k) = self % jGradEta(:,i,j,k) + Ja2CGL(:,l,m,n) &
                                          * spAxi % TCheb2Gauss(i,l) &
                                          * spAeta % TCheb2Gauss(j,m) &
                                          * spAzeta % TCheb2Gauss(k,n)

               self % jGradZeta(:,i,j,k) = self % jGradZeta(:,i,j,k) + Ja3CGL(:,l,m,n) &
                                          * spAxi % TCheb2Gauss(i,l) &
                                          * spAeta % TCheb2Gauss(j,m) &
                                          * spAzeta % TCheb2Gauss(k,n)

               self % jacobian(i,j,k) = self % jacobian(i,j,k) + JacobianCGL(l,m,n) &
                                          * spAxi % TCheb2Gauss(i,l) &
                                          * spAeta % TCheb2Gauss(j,m) &
                                          * spAzeta % TCheb2Gauss(k,n)
            end do              ; end do              ; end do
         end do               ; end do               ; end do
         do k = 0, self % Nz  ; do j = 0, self % Ny  ; do i = 0, self % Nx
            self % invJacobian(i,j,k) = 1._RP / self % jacobian(i,j,k)
         end do               ; end do               ; end do


      END SUBROUTINE computeMetricTermsCrossProductForm
!
!////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------
!  Computation of the metric terms on a face
!  -----------------------------------------
   subroutine ConstructMappedGeometryFace(self, Nf, Nelf, Nel, Nel3D, geom, hexMap, side, projType, eSide, rot, sliding, fID,scale)
      use PhysicsStorage
      use InterpolationMatrices
      implicit none
      class(MappedGeometryFace), intent(inout)  :: self
      integer,                   intent(in)     :: Nf(2)    ! Face polynomial order
      integer,                   intent(in)     :: Nelf(2)  ! Element face pOrder (with rotation)
      integer,                   intent(in)     :: Nel(2)   ! Element face pOrder (without rotation)
      integer,                   intent(in)     :: Nel3D(3) ! Element pOrder
      type(MappedGeometry),      intent(in)     :: geom
      type(TransfiniteHexMap),   intent(in)     :: hexMap
      integer,                   intent(in)     :: side
      integer,                   intent(in)     :: projType
      integer,                   intent(in)     :: eSide
      integer,                   intent(in)     :: rot

      logical,   optional,       intent(in)     :: sliding 
      integer,   optional,       intent(in)     :: fID 
      real(kind=RP), optional,   intent(in)     :: scale
!
!     ---------------
!     Local variables
!     ---------------
!
      integer        :: i, j, k, l, m, ii, jj, p, q
      real(kind=RP)  :: xi, eta
      real(kind=RP)  :: dS      (NDIM,0:Nel(1),0:Nel(2))
      real(kind=RP)  :: GradXi  (NDIM,0:Nel(1),0:Nel(2))
      real(kind=RP)  :: GradEta (NDIM,0:Nel(1),0:Nel(2))
      real(kind=RP)  :: GradZeta(NDIM,0:Nel(1),0:Nel(2))
      real(kind=RP)  :: dSrot      (NDIM,0:Nelf(1),0:Nelf(2))
      real(kind=RP)  :: GradXiRot  (NDIM,0:Nelf(1),0:Nelf(2))
      real(kind=RP)  :: GradEtaRot (NDIM,0:Nelf(1),0:Nelf(2))
      real(kind=RP)  :: GradZetaRot(NDIM,0:Nelf(1),0:Nelf(2))
      real(kind=RP)  :: x(3)
      real(kind=RP), allocatable  :: xx(:,:,:)
      real(kind=RP) :: MInt(0:Nelf(1),0:Nf(1),1:2)
      real(kind=RP) ::  xxx(NDIM, 0:Nelf(1), 0:Nf(2),2)
      real(kind=RP) ::  xrot(NDIM, 0:Nelf(1), 0:Nf(2))
      real(kind=RP) ::  xrot2(NDIM, 0:Nelf(1), 0:Nf(2))
      ! (present(sliding)) then 
         !(:,:,1)=transpose(TsetM(Nelf(1), Nf(1),1,1) % T)
        !! MInt(:,:,2)=transpose(TsetM(Nelf(1), Nf(1),3,1) % T)
     ! end if 

      if (.not.allocated(self % jacobian)) allocate( self % jacobian(0:Nf(1), 0:Nf(2)))
      if (.not.allocated(self % x)) allocate( self % x       (NDIM, 0:Nf(1), 0:Nf(2)))
      if (.not.allocated(self % normal))  allocate( self % normal  (NDIM, 0:Nf(1), 0:Nf(2)))
      if (.not.allocated(self % GradXi)) allocate( self % GradXi  (NDIM, 0:Nf(1), 0:Nf(2)))
      if (.not.allocated(self % GradEta)) allocate( self % GradEta (NDIM, 0:Nf(1), 0:Nf(2)))
      if (.not.allocated(self % GradZeta)) allocate( self % GradZeta(NDIM, 0:Nf(1), 0:Nf(2)))
      if (.not.allocated(self % t1))  allocate( self % t1      (NDIM, 0:Nf(1), 0:Nf(2)))
      if (.not.allocated(self % t2))  allocate( self % t2      (NDIM, 0:Nf(1), 0:Nf(2)))
 
      if (allocated(self % jacobian))  self % jacobian=0.0_RP
      if (allocated(self % x))  self % x =0.0_RP
      if (allocated(self % normal))   self % normal=0.0_RP
      if (allocated(self % GradXi))  self % GradXi =0.0_RP
      if (allocated(self % GradEta))  self % GradEta =0.0_RP
      if (allocated(self % GradZeta))  self % GradZeta=0.0_RP
      if (allocated(self % t1))   self % t1     =0.0_RP
      if (allocated(self % t2))   self % t2      =0.0_RP

      dS = 0.0_RP
      GradXi   = 0.0_RP
      GradEta  = 0.0_RP
      GradZeta = 0.0_RP

      select case(side)
         case(ELEFT)
!
!           Get face coordinates
!           --------------------
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               call coordRotation(NodalStorage(Nf(1)) % x(i), NodalStorage(Nf(2)) % x(j), rot, xi, eta)
               x = [-1.0_RP, xi, eta]
               self % x(:,i,j) = hexMap % transfiniteMapAt(x)
            end do ; end do
!  
!           Get surface Jacobian and normal vector
!           --------------------------------------
            do k = 0, Nel3D(3) ; do j = 0, Nel3D(2) ; do i = 0, Nel3D(1)
               dS(:,j,k) = dS(:,j,k) + geom % jGradXi(:,i,j,k) * NodalStorage(Nel3D(1)) % v(i,LEFT)
               GradXi  (:,j,k) = GradXi  (:,j,k) + geom % jGradXi  (:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(1)) % v(i,LEFT)
               GradEta (:,j,k) = GradEta (:,j,k) + geom % jGradEta (:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(1)) % v(i,LEFT)
               GradZeta(:,j,k) = GradZeta(:,j,k) + geom % jGradZeta(:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(1)) % v(i,LEFT)
            end do           ; end do           ; end do
!
!           Swap orientation
!           ----------------
            dS = -dS

         case(ERIGHT)
!
!           Get face coordinates
!           --------------------
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               call coordRotation(NodalStorage(Nf(1)) % x(i), NodalStorage(Nf(2)) % x(j), rot, xi, eta)
               x = [ 1.0_RP, xi, eta ]
               self % x(:,i,j) = hexMap % transfiniteMapAt(x)
            end do ; end do
!
!           Get surface Jacobian and normal vector
!           --------------------------------------
            do k = 0, Nel3D(3) ; do j = 0, Nel3D(2) ; do i = 0, Nel3D(1)
               dS(:,j,k) = dS(:,j,k) + geom % jGradXi(:,i,j,k) * NodalStorage(Nel3D(1)) % v(i,RIGHT)
               GradXi  (:,j,k) = GradXi  (:,j,k) + geom % jGradXi  (:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(1)) % v(i,RIGHT)
               GradEta (:,j,k) = GradEta (:,j,k) + geom % jGradEta (:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(1)) % v(i,RIGHT)
               GradZeta(:,j,k) = GradZeta(:,j,k) + geom % jGradZeta(:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(1)) % v(i,RIGHT)
            end do           ; end do           ; end do

         case(EBOTTOM)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               call coordRotation(NodalStorage(Nf(1)) % x(i), NodalStorage(Nf(2)) % x(j), rot, xi, eta)
               x = [xi, eta,-1.0_RP]
               self % x(:,i,j) = hexMap % transfiniteMapAt(x)
            end do ; end do
!
!           Get surface Jacobian and normal vector
!           --------------------------------------
            do k = 0, Nel3D(3) ; do j = 0, Nel3D(2) ; do i = 0, Nel3D(1)
               dS(:,i,j) = dS(:,i,j) + geom % jGradZeta(:,i,j,k) * NodalStorage(Nel3D(3)) % v(k,BOTTOM)
               GradXi  (:,i,j) = GradXi  (:,i,j) + geom % jGradXi  (:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(3)) % v(k,BOTTOM)
               GradEta (:,i,j) = GradEta (:,i,j) + geom % jGradEta (:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(3)) % v(k,BOTTOM)
               GradZeta(:,i,j) = GradZeta(:,i,j) + geom % jGradZeta(:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(3)) % v(k,BOTTOM)
            end do           ; end do           ; end do
!
!           Swap orientation
!           ----------------
            dS = -dS

         case(ETOP)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               call coordRotation(NodalStorage(Nf(1)) % x(i), NodalStorage(Nf(2)) % x(j), rot, xi, eta)
               x = [xi, eta, 1.0_RP]
               self % x(:,i,j) = hexMap % transfiniteMapAt(x)
            end do ; end do
!
!           Get surface Jacobian and normal vector
!           --------------------------------------
            do k = 0, Nel3D(3) ; do j = 0, Nel3D(2) ; do i = 0, Nel3D(1)
               dS(:,i,j) = dS(:,i,j) + geom % jGradZeta(:,i,j,k) * NodalStorage(Nel3D(3)) % v(k,TOP)
               GradXi  (:,i,j) = GradXi  (:,i,j) + geom % jGradXi  (:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(3)) % v(k,TOP)
               GradEta (:,i,j) = GradEta (:,i,j) + geom % jGradEta (:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(3)) % v(k,TOP)
               GradZeta(:,i,j) = GradZeta(:,i,j) + geom % jGradZeta(:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(3)) % v(k,TOP)
            end do           ; end do           ; end do

         case(EFRONT)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               call coordRotation(NodalStorage(Nf(1)) % x(i), NodalStorage(Nf(2)) % x(j), rot, xi, eta)
               x = [xi, -1.0_RP, eta]
               self % x(:,i,j) = hexMap % transfiniteMapAt(x)
            end do ; end do
!
!           Get surface Jacobian and normal vector
!           --------------------------------------
            do k = 0, Nel3D(3) ; do j = 0, Nel3D(2) ; do i = 0, Nel3D(1)
               dS(:,i,k) = dS(:,i,k) + geom % jGradEta(:,i,j,k) * NodalStorage(Nel3D(2)) % v(j,FRONT)
               GradXi  (:,i,k) = GradXi  (:,i,k) + geom % jGradXi  (:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(2)) % v(j,FRONT)
               GradEta (:,i,k) = GradEta (:,i,k) + geom % jGradEta (:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(2)) % v(j,FRONT)
               GradZeta(:,i,k) = GradZeta(:,i,k) + geom % jGradZeta(:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(2)) % v(j,FRONT)
            end do           ; end do           ; end do
!
!           Swap orientation
!           ----------------
            dS = -dS

         case(EBACK)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               call coordRotation(NodalStorage(Nf(1)) % x(i), NodalStorage(Nf(2)) % x(j), rot, xi, eta)
               x = [xi, 1.0_RP, eta]
               self % x(:,i,j) = hexMap % transfiniteMapAt(x)
            end do ; end do
            !write(*,*)'calculating x,side EBACK',  self % x
!
!           Get surface Jacobian and normal vector
!           --------------------------------------
            do k = 0, Nel3D(3) ; do j = 0, Nel3D(2) ; do i = 0, Nel3D(1)
               dS(:,i,k) = dS(:,i,k) + geom % jGradEta(:,i,j,k) * NodalStorage(Nel3D(2)) % v(j,BACK)
               GradXi  (:,i,k) = GradXi  (:,i,k) + geom % jGradXi  (:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(2)) % v(j,BACK)
               GradEta (:,i,k) = GradEta (:,i,k) + geom % jGradEta (:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(2)) % v(j,BACK)
               GradZeta(:,i,k) = GradZeta(:,i,k) + geom % jGradZeta(:,i,j,k) * geom % invJacobian(i,j,k) * NodalStorage(Nel3D(2)) % v(j,BACK)
            end do           ; end do           ; end do

      end select
!
!     Change the orientation depending on whether left or right elements are used
!     ---------------------------------------------------------------------------
      if ( eSide .eq. 2 ) dS = -dS

      !if ( eSide .eq. 2 ) write(*,*) 'eSide=2 line 700 of mappedgeometry'
!
!     Perform the rotation
!     --------------------
      if ( rot .eq. 0 ) then
         dSRot = dS           ! Considered separated since is very frequent
         GradXiRot   = GradXi
         GradEtaRot  = GradEta
         GradZetaRot = GradZeta

      else
         do j = 0, Nelf(2) ; do i = 0, Nelf(1)
            call leftIndexes2Right(i,j,Nelf(1), Nelf(2), rot, ii, jj)

            dSRot(:,i,j) = dS(:,ii,jj)
            GradXiRot  (:,i,j) = GradXi  (:,ii,jj)
            GradEtaRot (:,i,j) = GradEta (:,ii,jj)
            GradZetaRot(:,i,j) = GradZeta(:,ii,jj)
         end do            ; end do
      end if
     
!
!     Perform p-Adaption
!     ------------------
      select case(projType)
      case (0)
         self % normal = dSRot
         self % GradXi   = GradXiRot
         self % GradEta  = GradEtaRot
         self % GradZeta = GradZetaRot
      case (1)
         self % normal = 0.0_RP
         self % GradXi   = 0.0_RP
         self % GradEta  = 0.0_RP
         self % GradZeta = 0.0_RP
         do j = 0, Nf(2)  ; do l = 0, Nelf(1)   ; do i = 0, Nf(1)
            self % normal(:,i,j) = self % normal(:,i,j) + Tset(Nelf(1), Nf(1)) % T(i,l) * dSRot(:,l,j)
            self % GradXi  (:,i,j) = self % GradXi  (:,i,j) + Tset(Nelf(1), Nf(1)) % T(i,l) * GradXiRot  (:,l,j)
            self % GradEta (:,i,j) = self % GradEta (:,i,j) + Tset(Nelf(1), Nf(1)) % T(i,l) * GradEtaRot (:,l,j)
            self % GradZeta(:,i,j) = self % GradZeta(:,i,j) + Tset(Nelf(1), Nf(1)) % T(i,l) * GradZetaRot(:,l,j)
         end do                  ; end do                   ; end do

      case (2)
         self % normal = 0.0_RP
         self % GradXi   = 0.0_RP
         self % GradEta  = 0.0_RP
         self % GradZeta = 0.0_RP
         do l = 0, Nelf(2)  ; do j = 0, Nf(2)   ; do i = 0, Nf(1)
            self % normal(:,i,j) = self % normal(:,i,j) + Tset(Nelf(2), Nf(2)) % T(j,l) * dSRot(:,i,l)
            self % GradXi  (:,i,j) = self % GradXi  (:,i,j) + Tset(Nelf(2), Nf(2)) % T(j,l) * GradXiRot  (:,i,l)
            self % GradEta (:,i,j) = self % GradEta (:,i,j) + Tset(Nelf(2), Nf(2)) % T(j,l) * GradEtaRot (:,i,l)
            self % GradZeta(:,i,j) = self % GradZeta(:,i,j) + Tset(Nelf(2), Nf(2)) % T(j,l) * GradZetaRot(:,i,l)
         end do                  ; end do                   ; end do

      case (3)
         self % normal = 0.0_RP
         self % GradXi   = 0.0_RP
         self % GradEta  = 0.0_RP
         self % GradZeta = 0.0_RP
         do l = 0, Nelf(2)  ; do j = 0, Nf(2)
            do m = 0, Nelf(1) ; do i = 0, Nf(1)
               self % normal(:,i,j) = self % normal(:,i,j) +   Tset(Nelf(1), Nf(1)) % T(i,m) &
                                         * Tset(Nelf(2), Nf(2)) % T(j,l) &
                                         * dSRot(:,m,l)
               self % GradXi  (:,i,j) = self % GradXi  (:,i,j) +   Tset(Nelf(1), Nf(1)) % T(i,m) &
                                         * Tset(Nelf(2), Nf(2)) % T(j,l) &
                                         * GradXiRot  (:,m,l)
               self % GradEta (:,i,j) = self % GradEta (:,i,j) +   Tset(Nelf(1), Nf(1)) % T(i,m) &
                                         * Tset(Nelf(2), Nf(2)) % T(j,l) &
                                         * GradEtaRot (:,m,l)
               self % GradZeta(:,i,j) = self % GradZeta(:,i,j) +   Tset(Nelf(1), Nf(1)) % T(i,m) &
                                         * Tset(Nelf(2), Nf(2)) % T(j,l) &
                                         * GradZetaRot(:,m,l)
            end do                 ; end do
         end do                  ; end do
      end select

!
!     Perform h/p-Adaption if it's a sliding mesh
!     ------------------------------------------

      if (present(sliding) .and. present(fID)) then 
         if (sliding) then 
         allocate(xx(NDIM, 0:Nf(1), 0:Nf(2)))
        
         xx= self % x 
         xrot=self % x 
         xxx=0.0_RP
         self % x = 0.0_RP 
         self % normal = 0.0_RP
         self % GradXi   = 0.0_RP
         self % GradEta  = 0.0_RP
         self % GradZeta = 0.0_RP

         if(fID==1) then !2
            do l = 0, Nelf(2)  ; do j = 0, Nf(2)   ; do i = 0, Nf(1)    !3;1
              ! self % x(:,i,j)= self % x(:,i,j) + TsetM(Nelf(1), Nf(1), 4, 1) % T(j,l) * xrot(:,i,l)
               self % x(:,i,j)= self % x(:,i,j) + TsetM(Nelf(1), Nf(1), 3, 1) % T(j,l) * xrot(:,i,l)
               self % normal(:,i,j) = self % normal(:,i,j) + TsetM(Nelf(1), Nf(1), 3, 1) % T(j,l) * dSRot(:,i,l)
               self % GradXi  (:,i,j) = self % GradXi  (:,i,j) + TsetM(Nelf(1), Nf(1), 3, 1) % T(j,l) * GradXiRot  (:,i,l)
               self % GradEta (:,i,j) = self % GradEta (:,i,j) + TsetM(Nelf(1), Nf(1),3 , 1) % T(j,l) * GradEtaRot (:,i,l)
               self % GradZeta(:,i,j) = self % GradZeta(:,i,j) + TsetM(Nelf(1), Nf(1), 3, 1) % T(j,l) * GradZetaRot(:,i,l)
            end do                  ; end do                   ; end do
            self % GradXi=self % GradXi*scale
            self % GradEta=self % GradEta *scale
            self % GradZeta=self % GradZeta*scale
         else 
            !do j = 0, Nf(2)  ; do l = 0, Nelf(1)   ; do i = 0, Nf(1)
            do l = 0, Nelf(2)  ; do j = 0, Nf(2)   ; do i = 0, Nf(1)        !1;1
               self % x(:,i,j)= self % x(:,i,j) + TsetM(Nelf(1), Nf(1), 1, 1) % T(j,l) * xx(:,i,l)
               self % normal(:,i,j) = self % normal(:,i,j) + TsetM(Nelf(1), Nf(1), 1, 1) % T(j,l)* dSRot(:,i,l)
               self % GradXi  (:,i,j) = self % GradXi  (:,i,j) + TsetM(Nelf(1), Nf(1), 1, 1) % T(j,l) * GradXiRot  (:,i,l)
               self % GradEta (:,i,j) = self % GradEta (:,i,j) + TsetM(Nelf(1), Nf(1), 1, 1) % T(j,l) * GradEtaRot (:,i,l)
               self % GradZeta(:,i,j) = self % GradZeta(:,i,j) + TsetM(Nelf(1), Nf(1), 1, 1) % T(j,l) * GradZetaRot(:,i,l)

            end do                  ; end do                   ; end do
            self % GradXi=self % GradXi*scale
            self % GradEta=self % GradEta *scale
            self % GradZeta=self % GradZeta*scale
         end if 
         deallocate(xx)
      end if 
      end if 

!
!     Compute
!     -------
      do j = 0, Nf(2)   ; do i = 0, Nf(1)
         self % jacobian(i,j) = norm2(self % normal(:,i,j))
         self % normal(:,i,j) = self % normal(:,i,j) / self % jacobian(i,j)
      end do            ; end do
!
!     Compute tangent vectors
!     -----------------------
      self % t1 = 0.0_RP
      do j = 0, Nf(2)   ; do l = 0, Nf(1)    ; do i = 0, Nf(1)
         self % t1(:,i,j) = self % t1(:,i,j) + NodalStorage(Nf(1)) % D(i,l) * self % x(:,l,j)
      end do            ; end do             ; end do
!
!     Tangent vectors could not be orthogonal to normal vectors (because of truncation errors)
!     ----------------------------------------------------------------------------------------
      do j = 0, Nf(2)   ; do i = 0, Nf(1)
         self % t1(:,i,j)  = self % t1(:,i,j) - dot_product(self % t1(:,i,j), self % normal(:,i,j)) &
                                              * self % normal(:,i,j)

         self % t1(:,i,j)  = self % t1(:,i,j)  / norm2(self % t1(:,i,j))
         call vCross(self % normal(:,i,j), self % t1(:,i,j), self % t2(:,i,j))
      end do            ; end do
!
!     ------------
!     Face surface
!     ------------
!
      self % surface = 0.0_RP
      do j = 0, Nf(2) ; do i = 0, Nf(1)
         self % surface = self % surface + NodalStorage(Nf(1)) % w(i) * NodalStorage(Nf(2)) % w(j) * self % jacobian(i,j)
      end do          ; end do

   end subroutine ConstructMappedGeometryFace
!
!////////////////////////////////////////////////////////////////////////
!
      pure subroutine DestructMappedGeometryFace(self)
         implicit none
         !-------------------------------------------------------------------
         class(MappedGeometryFace), intent(inout) :: self
         !-------------------------------------------------------------------

         safedeallocate(self % x        )
         safedeallocate(self % jacobian )
         safedeallocate(self % GradXi   )
         safedeallocate(self % GradEta  )
         safedeallocate(self % GradZeta )
         safedeallocate(self % normal   )
         safedeallocate(self % t1       )
         safedeallocate(self % t2       )
         safedeallocate(self % dWall    )

      end subroutine DestructMappedGeometryFace

!
!///////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!!     Returns the jacobian of the transformation computed from
!!     the three co-variant coordinate vectors.
!-------------------------------------------------------------------------------
!
      FUNCTION jacobian3D(a1,a2,a3)
!
      USE SMConstants
      IMPLICIT NONE

      REAL(KIND=RP)               :: jacobian3D
      REAL(KIND=RP), DIMENSION(3) :: a1,a2,a3,v
!
      CALL vCross(a2,a3,v)
      jacobian3D = vDot(a1,v)

      END FUNCTION jacobian3D
!
!///////////////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!!    Returns in result the cross product u x v
!-------------------------------------------------------------------------------
!
      SUBROUTINE vCross(u,v,result)
!
      IMPLICIT NONE

      REAL(KIND=RP), DIMENSION(3) :: u,v,result

      result(1) = u(2)*v(3) - v(2)*u(3)
      result(2) = u(3)*v(1) - v(3)*u(1)
      result(3) = u(1)*v(2) - v(1)*u(2)

      END SUBROUTINE vCross
!
!///////////////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!!    Returns the dot product u.v
!-------------------------------------------------------------------------------
!
      FUNCTION vDot(u,v)
!
      IMPLICIT NONE

      REAL(KIND=RP)               :: vDot
      REAL(KIND=RP), DIMENSION(3) :: u,v

      vDot = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)

      END FUNCTION vDot
!
!///////////////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!!    Returns the 2-norm of u
!-------------------------------------------------------------------------------
!
      FUNCTION vNorm(u)
!
      IMPLICIT NONE

      REAL(KIND=RP)               :: vNorm
      REAL(KIND=RP), DIMENSION(3) :: u

      vNorm = SQRT(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))

      END FUNCTION vNorm

      subroutine SetMappingsToCrossProduct()
         implicit none

         useCrossProductMetrics = .true.

      end subroutine SetMappingsToCrossProduct

!
!///////////////////////////////////////////////////////////////////////////////
!
      pure subroutine ComputeVorticity(U_x, U_y, U_z, vorticity)

         real(kind=RP), intent(in)  :: U_x(3)
         real(kind=RP), intent(in)  :: U_y(3)
         real(kind=RP), intent(in)  :: U_z(3)
         real(kind=RP), intent(out) :: vorticity(3)

         vorticity(1) = U_y(3) - U_z(2)
         vorticity(2) = U_z(1) - U_x(3)
         vorticity(3) = U_x(2) - U_y(1)

      end subroutine ComputeVorticity
!
!///////////////////////////////////////////////////////////////////////////////
!
END Module MappedGeometryClass
