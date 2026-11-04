#include "Includes.h"
module RotaryFaceRotations
!
!  ***********************************************************************
!    Rotary periodic faces: the two sides live in frames rotated by the
!    sector angle. Side 2 (or the received MPI side) is rotated into the
!    frame of the other side after prolongation, so the standard face
!    routines apply. Fluxes projected to side 2 are rotated back.
!  ***********************************************************************
!
   use SMConstants
   use MeshTypes,    only: HMESH_MPI
   use FaceClass,    only: Face
   use HexMeshClass, only: HexMesh
   use MPI_Process_Info, only: MPI_Process
   use RotationData, only: rotationParams
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none

   private
   public RotateRotaryFacesSolution, RotateRotaryFacesGradients
   public RotateRotaryMPIFacesSolution, RotateRotaryMPIFacesGradients
   public RotateBackRotaryFacesGradientFlux, RotateBackRotaryFacesFlux
   public CountRotaryFaces

   contains

   subroutine RotateRotaryFacesSolution(self, element_ids)
!
!     Interior faces: rotate the state prolonged to side 2 (the given elements, or all)
!
      implicit none
      class(HexMesh), intent(inout) :: self
      integer, optional, intent(in)    :: element_ids(:)   ! Elements just prolonged (all if absent)
      integer                          :: iEl, eID, k, fID, i, j, nEl

      if (.not. rotationParams % periodicEnabled) return

      nEl = size(self % elements)
      if (present(element_ids)) nEl = size(element_ids)
!$omp do schedule(runtime) private(eID, k, fID, i, j)
      do iEl = 1, nEl
         eID = iEl
         if (present(element_ids)) eID = element_ids(iEl)
         do k = 1, 6
            fID = self % elements(eID) % faceIDs(k)
            if (.not. IsRotarySide2(self % faces(fID), eID)) cycle
            do j = 0, self % faces(fID) % Nf(2) ; do i = 0, self % faces(fID) % Nf(1)
               self % faces(fID) % storage(2) % Q(2:4,i,j) = RotateVectorAroundAxis(self % faces(fID) % storage(2) % Q(2:4,i,j), &
                                                                                   -rotationParams % periodicAngle, rotationParams % axis)
            end do ; end do
         end do
      end do
!$omp end do
   end subroutine RotateRotaryFacesSolution

   subroutine RotateRotaryFacesGradients(self, element_ids)
!
!     Interior faces: rotate the gradients prolonged to side 2 (the given elements, or all)
!
      implicit none
      class(HexMesh), intent(inout) :: self
      integer, optional, intent(in)    :: element_ids(:)   ! Elements just prolonged (all if absent)
      integer                          :: iEl, eID, k, fID, i, j, nEl

      if (.not. rotationParams % periodicEnabled) return

      nEl = size(self % elements)
      if (present(element_ids)) nEl = size(element_ids)
!$omp do schedule(runtime) private(eID, k, fID, i, j)
      do iEl = 1, nEl
         eID = iEl
         if (present(element_ids)) eID = element_ids(iEl)
         do k = 1, 6
            fID = self % elements(eID) % faceIDs(k)
            if (.not. IsRotarySide2(self % faces(fID), eID)) cycle
            do j = 0, self % faces(fID) % Nf(2) ; do i = 0, self % faces(fID) % Nf(1)
               call RotateGradients(self % faces(fID) % storage(2) % U_x(:,i,j), self % faces(fID) % storage(2) % U_y(:,i,j), &
                                    self % faces(fID) % storage(2) % U_z(:,i,j), -rotationParams % periodicAngle)
            end do ; end do
         end do
      end do
!$omp end do
   end subroutine RotateRotaryFacesGradients

   subroutine RotateRotaryMPIFacesSolution(self)
!
!     MPI faces: rotate the received state into the frame of the local side
!
      implicit none
      class(HexMesh), intent(inout) :: self
      integer                          :: iFace, fID, side, i, j
      real(kind=RP)                    :: angle

      if (.not. rotationParams % periodicEnabled) return

!$omp do schedule(runtime) private(fID, side, angle, i, j)
      do iFace = 1, size(self % faces_mpi)
         fID = self % faces_mpi(iFace)
         if (.not. self % faces(fID) % isRotaryPeriodic) cycle
         side  = ReceivedSide(self % faces(fID))
         angle = merge(-1.0_RP, 1.0_RP, side == 2) * rotationParams % periodicAngle
         do j = 0, self % faces(fID) % Nf(2) ; do i = 0, self % faces(fID) % Nf(1)
            self % faces(fID) % storage(side) % Q(2:4,i,j) = RotateVectorAroundAxis(self % faces(fID) % storage(side) % Q(2:4,i,j), &
                                                                                   angle, rotationParams % axis)
         end do ; end do
      end do
!$omp end do
   end subroutine RotateRotaryMPIFacesSolution

   subroutine RotateRotaryMPIFacesGradients(self)
!
!     MPI faces: rotate the received gradients into the frame of the local side
!
      implicit none
      class(HexMesh), intent(inout) :: self
      integer                          :: iFace, fID, side, i, j
      real(kind=RP)                    :: angle

      if (.not. rotationParams % periodicEnabled) return

!$omp do schedule(runtime) private(fID, side, angle, i, j)
      do iFace = 1, size(self % faces_mpi)
         fID = self % faces_mpi(iFace)
         if (.not. self % faces(fID) % isRotaryPeriodic) cycle
         side  = ReceivedSide(self % faces(fID))
         angle = merge(-1.0_RP, 1.0_RP, side == 2) * rotationParams % periodicAngle
         do j = 0, self % faces(fID) % Nf(2) ; do i = 0, self % faces(fID) % Nf(1)
            call RotateGradients(self % faces(fID) % storage(side) % U_x(:,i,j), self % faces(fID) % storage(side) % U_y(:,i,j), &
                                 self % faces(fID) % storage(side) % U_z(:,i,j), angle)
         end do ; end do
      end do
!$omp end do
   end subroutine RotateRotaryMPIFacesGradients

   subroutine RotateBackRotaryFacesGradientFlux(self, no_of_faces, face_ids)
!
!     Interior faces: rotate the gradient flux projected to side 2 back to its frame
!
      implicit none
      class(HexMesh), intent(inout) :: self
      integer,          intent(in)    :: no_of_faces
      integer,          intent(in)    :: face_ids(no_of_faces)   ! Faces just computed
      integer                          :: iFace, fID, i, j

      if (.not. rotationParams % periodicEnabled) return

!$omp do schedule(runtime) private(fID, i, j)
      do iFace = 1, no_of_faces
         fID = face_ids(iFace)
         if (.not. self % faces(fID) % isRotaryPeriodic .or. self % faces(fID) % faceType == HMESH_MPI) cycle
         do j = 0, self % faces(fID) % Nf(2) ; do i = 0, self % faces(fID) % Nf(1)
            call RotateGradients(self % faces(fID) % storage(2) % unStar(:,IX,i,j), self % faces(fID) % storage(2) % unStar(:,IY,i,j), &
                                 self % faces(fID) % storage(2) % unStar(:,IZ,i,j), rotationParams % periodicAngle)
         end do ; end do
      end do
!$omp end do
   end subroutine RotateBackRotaryFacesGradientFlux

   subroutine RotateBackRotaryFacesFlux(self, no_of_faces, face_ids)
!
!     Interior faces: rotate the flux projected to side 2 back to its frame
!
      implicit none
      class(HexMesh), intent(inout) :: self
      integer,          intent(in)    :: no_of_faces
      integer,          intent(in)    :: face_ids(no_of_faces)   ! Faces just computed
      integer                          :: iFace, fID, i, j

      if (.not. rotationParams % periodicEnabled) return
!$omp barrier   ! The flux loop before this one has no barrier
!$omp do schedule(runtime) private(fID, i, j)
      do iFace = 1, no_of_faces
         fID = face_ids(iFace)
         if (.not. self % faces(fID) % isRotaryPeriodic .or. self % faces(fID) % faceType == HMESH_MPI) cycle
         do j = 0, self % faces(fID) % Nf(2) ; do i = 0, self % faces(fID) % Nf(1)
            self % faces(fID) % storage(2) % FStar(2:4,i,j) = RotateVectorAroundAxis(self % faces(fID) % storage(2) % FStar(2:4,i,j), &
                                                                                    rotationParams % periodicAngle, rotationParams % axis)
         end do ; end do
      end do
!$omp end do
   end subroutine RotateBackRotaryFacesFlux

   logical function IsRotarySide2(f, eID)
!     Is element eID side 2 of the interior rotary periodic face f?
      implicit none
      type(Face), intent(in) :: f
      integer,    intent(in) :: eID
      IsRotarySide2 = f % isRotaryPeriodic .and. f % faceType /= HMESH_MPI .and. f % elementIDs(2) == eID
   end function IsRotarySide2

   integer function ReceivedSide(f)
!     Side of an MPI face received from the neighbour rank (the local element has the larger ID)
      implicit none
      type(Face), intent(in) :: f
      ReceivedSide = merge(2, 1, f % elementIDs(1) >= f % elementIDs(2))
   end function ReceivedSide

   integer function CountRotaryFaces(self) result(nRotaryFaces)
!
!     Number of rotary periodic faces in the whole mesh (an MPI face is present on both ranks)
!
      implicit none
      class(HexMesh), intent(in) :: self
      integer                   :: fID, nLocal, nMPILocal, nMPI, ierr

      nLocal    = 0
      nMPILocal = 0
      do fID = 1, size(self % faces)
         if (.not. self % faces(fID) % isRotaryPeriodic) cycle
         nLocal = nLocal + 1
         if (self % faces(fID) % faceType == HMESH_MPI) nMPILocal = nMPILocal + 1
      end do

      nRotaryFaces = nLocal
      nMPI         = nMPILocal
#ifdef _HAS_MPI_
      if (MPI_Process % doMPIAction) then
         call mpi_allreduce(nLocal,    nRotaryFaces, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         call mpi_allreduce(nMPILocal, nMPI,         1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      end if
#endif
      nRotaryFaces = nRotaryFaces - nMPI / 2

   end function CountRotaryFaces
!
!  ------------------------------------------------------------------------
!  Rotation math
!  ------------------------------------------------------------------------
!
   subroutine BuildRotationMatrix(R, angle, axis)
     ! Construct the 3x3 rotation matrix for a given rotation angle and a UNIT axis
     implicit none
     real(kind=RP), intent(out) :: R(3,3)
     real(kind=RP), intent(in)  :: angle
     real(kind=RP), intent(in)  :: axis(3)
     real(kind=RP) :: c, s, t
     real(kind=RP) :: x, y, z

     x = axis(1)
     y = axis(2)
     z = axis(3)

     c = cos(angle)
     s = sin(angle)
     t = 1.0_rp - c

     R(1,1) = t*x*x + c
     R(1,2) = t*x*y - s*z
     R(1,3) = t*x*z + s*y

     R(2,1) = t*x*y + s*z
     R(2,2) = t*y*y + c
     R(2,3) = t*y*z - s*x

     R(3,1) = t*x*z - s*y
     R(3,2) = t*y*z + s*x
     R(3,3) = t*z*z + c

   end subroutine BuildRotationMatrix

   function RotateVectorAroundAxis(v, angle, axis) result(vRot)
     implicit none
     real(kind=RP), intent(in) :: v(3), axis(3)
     real(kind=RP), intent(in) :: angle
     real(kind=RP)             :: vRot(3)
     real(kind=RP)             :: R(3,3)

     call BuildRotationMatrix(R, angle, axis)
     vRot = matmul(R, v)

   end function

   subroutine RotateGradients(U_x, U_y, U_z, rotAngle)
!
!     Rotate the gradient tensor G = [U_x U_y U_z] of one point by rotAngle around the
!     rotation axis: G <- R G R^T
!        1) Rotate the momentum components (rows 2:4) of each column U_x, U_y, U_z
!        2) Rotate the spatial directions, i.e. [U_x(k), U_y(k), U_z(k)] for every variable k
!     Used for the gradients and for the gradient flux (unStar) of rotary periodic faces.
!
      implicit none
      real(kind=RP), intent(inout) :: U_x(:), U_y(:), U_z(:)
      real(kind=RP), intent(in)    :: rotAngle
      real(kind=RP)                :: R(3,3), v(NDIM)
      integer                      :: k

      call BuildRotationMatrix(R, rotAngle, rotationParams % axis)

      U_x(2:4) = matmul(R, U_x(2:4))
      U_y(2:4) = matmul(R, U_y(2:4))
      U_z(2:4) = matmul(R, U_z(2:4))

      do k = 1, size(U_x)
         v      = matmul(R, [U_x(k), U_y(k), U_z(k)])
         U_x(k) = v(1)
         U_y(k) = v(2)
         U_z(k) = v(3)
      end do

   end subroutine RotateGradients

end module RotaryFaceRotations
