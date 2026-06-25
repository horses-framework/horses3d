#include "Includes.h"
module GradientsStabilization
#if defined(CAHNHILLIARD)
   use SMConstants
   use MeshTypes
   use HexMeshClass
   use PhysicsStorage
   use VariableConversion
   use BoundaryConditions, only: BCs
   implicit none

   private
   public   LambdaEstimator_f

   public   StabilizeGradients

   abstract interface
      subroutine LambdaEstimator_f(mesh, fID, i, j, lambda)
         use SMConstants
         use HexMeshClass
         implicit none
         class(HexMesh),   intent(in)  :: mesh
         integer,          intent(in)  :: fID, i, j
         real(kind=RP),    intent(out) :: lambda
      end subroutine LambdaEstimator_f
   end interface
!
!  ========
   contains
!  ========
!
      subroutine StabilizeGradients(mesh, time)
         implicit none
         class(HexMesh)                :: mesh
         real(kind=RP),    intent(in)  :: time
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eID, fID, i, j, k, m
!
!        **************************************************
!        Compute the face stabilization flux in local faces
!        **************************************************
!
!$omp do schedule(runtime) 
         do fID = 1, SIZE(mesh % faces) 
            associate(f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_INTERIOR) 

               if (mesh % faces(fID) % IsMortar==1) then
               associate(unStar=>mesh% faces(fID)%storage(1)%unStar)
                  unStar=0.0_RP
               end associate
               associate(unStar=>mesh% faces(fID)%storage(2)%unStar)
                  unStar=0.0_RP
               end associate
               do m=1,4
                  if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                     call GradientsStabilization_InteriorFace(f=mesh % faces(mesh % faces(fID)%Mortar(m)),masterFace1=f) 
                  end if 
               end do 
               elseif (mesh % faces(fID) % IsMortar==0) then
               call GradientsStabilization_InteriorFace(f) 
               end if 
            case (HMESH_BOUNDARY) 
               call GradientsStabilization_BoundaryFace(f, time) 
 
            end select 
            end associate 
         end do            
!$omp end do 
!
!$omp do schedule(runtime) 
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID))
#ifdef _HAS_MPI_
            if ( e % hasSharedFaces ) cycle
#endif
!
!           Revert the scaling
!           -------------------               
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
               e % storage % c_x(:,i,j,k) = &
                           e % storage % c_x(:,i,j,k) * e % geom % jacobian(i,j,k)
               e % storage % c_y(:,i,j,k) = &
                           e % storage % c_y(:,i,j,k) * e % geom % jacobian(i,j,k)
               e % storage % c_z(:,i,j,k) = &
                           e % storage % c_z(:,i,j,k) * e % geom % jacobian(i,j,k)
            end do         ; end do          ; end do
!
!           Add the surface integrals
!           -------------------------
            call PerformSurfaceIntegrals( e, mesh)
!
!           Perform the scaling
!           -------------------               
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
               e % storage % c_x(:,i,j,k) = &
                           e % storage % c_x(:,i,j,k) / e % geom % jacobian(i,j,k)
               e % storage % c_y(:,i,j,k) = &
                           e % storage % c_y(:,i,j,k) / e % geom % jacobian(i,j,k)
               e % storage % c_z(:,i,j,k) = &
                           e % storage % c_z(:,i,j,k) / e % geom % jacobian(i,j,k)
            end do         ; end do          ; end do

            end associate
         end do
!$omp end do

#ifdef _HAS_MPI_
!$omp do schedule(runtime) 
         do fID = 1, SIZE(mesh % faces) 
            associate(f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_MPI) 
            if (mesh% faces(fID)%IsMortar==1) then 
               associate(unstar=>mesh% faces(fID)%storage(1)%unStar)
                  unstar=0.0_RP
               end associate
               do m=1, 4
                  if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                     call GradientsStabilization_InteriorFace(f=mesh % faces(mesh % faces(fID)%Mortar(m)),masterFace1=f) 
                  end if 
               end do 
            end if 
            call GradientsStabilization_MPIFace(f) 
            end select 
            end associate 
         end do            
!$omp end do 
!$omp single
         if ( mesh % nonconforming ) then
            call mesh % UpdateMPIFacesGradMortarflux(NCOMP)
         end if
   !$omp end single
   
   
   !$omp single
         if ( mesh % nonconforming ) then
            call mesh % GatherMPIFacesGradMortarFlux(NCOMP)
         end if
   !$omp end single
!
!$omp do schedule(runtime) 
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID))
            if ( .not. e % hasSharedFaces ) cycle
!
!           Revert the scaling
!           -------------------               
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
               e % storage % c_x(:,i,j,k) = &
                           e % storage % c_x(:,i,j,k) * e % geom % jacobian(i,j,k)
               e % storage % c_y(:,i,j,k) = &
                           e % storage % c_y(:,i,j,k) * e % geom % jacobian(i,j,k)
               e % storage % c_z(:,i,j,k) = &
                           e % storage % c_z(:,i,j,k) * e % geom % jacobian(i,j,k)
            end do         ; end do          ; end do

!
!           Add the surface integrals
!           -------------------------
            call PerformSurfaceIntegrals( e, mesh)
!
!           Perform the scaling
!           -------------------               
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
               e % storage % c_x(:,i,j,k) = &
                           e % storage % c_x(:,i,j,k) / e % geom % jacobian(i,j,k)
               e % storage % c_y(:,i,j,k) = &
                           e % storage % c_y(:,i,j,k) / e % geom % jacobian(i,j,k)
               e % storage % c_z(:,i,j,k) = &
                           e % storage % c_z(:,i,j,k) / e % geom % jacobian(i,j,k)
            end do         ; end do          ; end do

            end associate
         end do
!$omp end do
#endif
      end subroutine StabilizeGradients

      subroutine GradientsStabilization_InteriorFace(f, masterFace1)
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         type(Face)             :: f
         type(Face),optional    :: masterFace1

         procedure(LambdaEstimator_f)  :: LambdaEstimator
         
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i,j, lm
         real(kind=RP) :: cL(NCOMP), cR(NCOMP), uHat(NCOMP)
         real(kind=RP) :: Hflux(NCOMP,NDIM,0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP) :: vAver(NDIM)

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            cL = f % storage(1) % c(:,i,j)
            cR = f % storage(2) % c(:,i,j)
  
            vAver = AVERAGE( f % storage(1) % v(:,i,j) , f % storage(2) % v(:,i,j))

            Uhat = 0.5_RP * sign2(dot_product(vAver, f % geom % normal(:,i,j)))*(cL - cR) * f % geom % jacobian(i,j) 
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do
         if (f % IsMortar==0) then 
            call f % ProjectGradientFluxToElements(NCOMP, HFlux,(/1,2/),-1)
         end if 
         if (f % IsMortar==2 .and. present(masterFace1)) then 
            call masterFace1 % ProjectMortarGradientFluxToElements(nEqn=NCOMP, slaveFace=f, HFlux=HFlux,whichElements=(/0,2/),factor=-1)
            call f % ProjectGradientFluxToElements(NCOMP, HFlux,(/1,2/),-1)
         end if 
         
      end subroutine GradientsStabilization_InteriorFace   

      subroutine GradientsStabilization_MPIFace(f)
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         type(Face)    :: f
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i,j, thisSide
         real(kind=RP) :: cL(NCOMP), cR(NCOMP), Uhat(NCOMP)
         real(kind=RP) :: Hflux(NCOMP,NDIM,0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP) :: vAver(NDIM)

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            cL = f % storage(1) % c(:,i,j)
            cR = f % storage(2) % c(:,i,j)
   
            vAver = AVERAGE( f % storage(1) % v(:,i,j) , f % storage(2) % v(:,i,j))

            Uhat = 0.5_RP * sign2(dot_product(vAver, f % geom % normal(:,i,j)))*(cL - cR) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectGradientFluxToElements(NCOMP, HFlux, (/thisSide, HMESH_NONE/), -1)

         if (f % IsMortar==2) then 
            !write(*,*) 'this side', thisSide
            call f% Interpolatesmall2biggrad(NCOMP, HFlux)
            
         end if
      end subroutine GradientsStabilization_MPIFace   

      subroutine GradientsStabilization_BoundaryFace(f, time)
         use Physics
         use FaceClass
         implicit none
         type(Face)             :: f
         real(kind=RP)          :: time
         integer                :: i, j
         real(kind=RP)          :: Uhat(NCOMP), cL(NCOMP), cR(NCOMP)
         real(kind=RP)          :: bvExt(NCOMP)

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)

            bvExt =  f % storage(1) % c(:,i,j)
   
            call BCs(f % zone) % bc % PhaseFieldState( f % geom % x(:,i,j), & 
                                time               , &
                                f % geom % normal(:,i,j)      , &
                                bvExt              )
!   
!           -------------------
!           u, v, w, T averages
!           -------------------
!   
            cL = f % storage(1) % Q(:,i,j)
            cR = bvExt
   
            Uhat = 0.5_RP * dot_product(f % storage(1) % v(:,i,j), f % geom % normal(:,i,j))*(cL - cR) * f % geom % jacobian(i,j) / (abs(dot_product(f % storage(1) % v(:,i,j),f % geom % normal(:,i,j))) + epsilon(1.0_RP))
            f % storage(1) % unStar(:,1,i,j) = Uhat * f % geom % normal(1,i,j)
            f % storage(1) % unStar(:,2,i,j) = Uhat * f % geom % normal(2,i,j)
            f % storage(1) % unStar(:,3,i,j) = Uhat * f % geom % normal(3,i,j)

         end do ; end do   
         
      end subroutine GradientsStabilization_BoundaryFace

      subroutine PerformSurfaceIntegrals( e, mesh )
         use ElementClass
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use DGIntegrals
         implicit none
         class(Element)                      :: e
         class(HexMesh)                      :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: faceInt_x(NCOMP, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP)        :: faceInt_y(NCOMP, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP)        :: faceInt_z(NCOMP, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )

         call VectorWeakIntegrals % StdFace(e, NCOMP, &
               mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % unStar, &
               mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % unStar, &
               mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % unStar, &
               mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % unStar, &
               mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % unStar, &
               mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % unStar, &
               faceInt_x, faceInt_y, faceInt_z )

         e % storage % c_x = e % storage % c_x + faceInt_x
         e % storage % c_y = e % storage % c_y + faceInt_y
         e % storage % c_z = e % storage % c_z + faceInt_z

      end subroutine PerformSurfaceIntegrals
#endif
end module GradientsStabilization