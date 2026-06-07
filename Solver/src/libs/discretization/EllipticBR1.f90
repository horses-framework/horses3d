#include "Includes.h"
module EllipticBR2
   use SMConstants
   use Headers
   use MeshTypes
   use ElementClass
   use HexMeshClass
   use PhysicsStorage
   use Physics
   use MPI_Process_Info
   use MPI_Face_Class
   use EllipticDiscretizationClass
   use VariableConversion
   use FluidData
   use BoundaryConditions, only: BCs
   implicit none
!
!
   private
   public   BassiRebay2_t

   type, extends(EllipticDiscretization_t)   :: BassiRebay2_t
      real(kind=RP)        :: eta = 1.0_RP
      contains
         procedure      :: Construct          => BR2_Construct
         procedure      :: ComputeGradient    => BR2_ComputeGradient
         procedure      :: ComputeInnerFluxes => BR2_ComputeInnerFluxes
         procedure      :: RiemannSolver      => BR2_RiemannSolver
         procedure      :: Describe           => BR2_Describe
   end type BassiRebay2_t
!
!  ========
   contains
!  ========
!
      subroutine BR2_Construct(self, controlVariables, eqname)
         use FTValueDictionaryClass
         use mainKeywordsModule
         use MPI_Process_Info
         use PhysicsStorage
         implicit none
         class(BassiRebay2_t)                  :: self
         class(FTValueDictionary), intent(in)  :: controlVariables
         integer,          intent(in)          :: eqname
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)            :: BR2variant
!
!        ----------------------------------------------------------
!        Set the particular procedures to compute the elliptic flux
!        ----------------------------------------------------------
!
         select case (eqName)
         case(ELLIPTIC_NS)
            self % eqName = ELLIPTIC_NS
         
         case(ELLIPTIC_NSSA)
            self % eqName = ELLIPTIC_NSSA

         case(ELLIPTIC_iNS)
            self % eqName = ELLIPTIC_iNS

         case(ELLIPTIC_CH)
            self % eqName = ELLIPTIC_CH

         case(ELLIPTIC_MU)
            self % eqName = ELLIPTIC_MU

         case default
            print*, "Unrecognized equation"
            errorMessage(STD_OUT)
            error stop

         end select

!
!        Request the penalty parameter
!        -----------------------------
         if ( controlVariables % containsKey("penalty parameter") ) then
            self % eta = controlVariables % doublePrecisionValueForKey("penalty parameter")

         else
!            
!           Set 3.0 by default
!           ------------------
            self % eta = 2.0_RP

         end if
            
      end subroutine BR2_Construct
   
      subroutine BR2_Describe(self)
         implicit none
         class(BassiRebay2_t),   intent(in)  :: self
!
!        Display the configuration
!        -------------------------
         if (MPI_Process % isRoot) write(STD_OUT,'(/)')
         call Subsection_Header("Viscous discretization")

         if (.not. MPI_Process % isRoot ) return

         write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","BR2"

         write(STD_OUT,'(30X,A,A30,F10.3)') "->","Penalty parameter: ", self % eta

#ifdef NAVIERSTOKES
         select case (self % eqName)
         case (ELLIPTIC_NS,ELLIPTIC_NSSA)
            select case (grad_vars)
            case(GRADVARS_STATE);   write(STD_OUT,'(30X,A,A30,A)') "->","Gradient variables: ","State"
            case(GRADVARS_ENTROPY); write(STD_OUT,'(30X,A,A30,A)') "->","Gradient variables: ","Entropy"
            case(GRADVARS_ENERGY);  write(STD_OUT,'(30X,A,A30,A)') "->","Gradient variables: ","Energy"
            end select
         end select
#endif
      end subroutine BR2_Describe

      subroutine BR2_ComputeGradient( self , nEqn, nGradEqn, mesh , time , GetGradients, HO_Elements, element_mask, Level)
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use MPI_Process_Info
         implicit none
         class(BassiRebay2_t), intent(in) :: self
         integer,              intent(in) :: nEqn
         integer,              intent(in) :: nGradEqn
         class(HexMesh)                   :: mesh
         real(kind=RP),        intent(in) :: time
         procedure(GetGradientValues_f)   :: GetGradients
         logical, intent(in), optional    :: HO_Elements
         logical, intent(in), optional    :: element_mask(:)
		 integer, intent(in), optional    :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: Nx, Ny, Nz
         integer :: i, j, k, m
         integer :: eID , fID , dimID , eqID, fIDs(6), iFace, iEl
         logical :: HOElements

         if (present(HO_Elements)) then
            HOElements = HO_Elements
         else
            HOElements = .false.
         end if
!
!        ***********************
!        Compute local gradients
!        ***********************
!
         if (HOElements) then
!$omp do schedule(runtime)
            do eID = 1, size(mesh % HO_Elements)
               associate( e => mesh % elements(mesh % HO_Elements(eID)))
               call e % ComputeLocalGradient(nEqn, nGradEqn, GetGradients, .false.)
   !
   !           Prolong to faces
   !           ----------------
               fIDs = e % faceIDs
               if (.not.mesh%nonconforming) then
               call e % ProlongGradientsToFaces(nGradEqn, mesh % faces(fIDs(1)),&
                                                mesh % faces(fIDs(2)),&
                                                mesh % faces(fIDs(3)),&
                                                mesh % faces(fIDs(4)),&
                                                mesh % faces(fIDs(5)),&
                                                mesh % faces(fIDs(6)) )
               else 
               call e % ProlongGradientsToFaces(nGradEqn, mesh % faces(fIDs(1)),&
                                                mesh % faces(fIDs(2)),&
                                                mesh % faces(fIDs(3)),&
                                                mesh % faces(fIDs(4)),&
                                                mesh % faces(fIDs(5)),&
                                                mesh % faces(fIDs(6)),&
                                                faces=mesh % faces)
               end if 

               end associate
            end do
!$omp end do    
         else
!$omp do schedule(runtime)
            do eID = 1, size(mesh % elements)
										  
               associate( e => mesh % elements(eID) )
               call e % ComputeLocalGradient(nEqn, nGradEqn, GetGradients, .false.)
   !
   !           Prolong to faces
   !           ----------------
               fIDs = e % faceIDs
               if (.not.mesh%nonconforming) then
               call e % ProlongGradientsToFaces(nGradEqn, mesh % faces(fIDs(1)),&
                                                mesh % faces(fIDs(2)),&
                                                mesh % faces(fIDs(3)),&
                                                mesh % faces(fIDs(4)),&
                                                mesh % faces(fIDs(5)),&
                                                mesh % faces(fIDs(6)) )
               else 
               call e % ProlongGradientsToFaces(nGradEqn, mesh % faces(fIDs(1)),&
                                                mesh % faces(fIDs(2)),&
                                                mesh % faces(fIDs(3)),&
                                                mesh % faces(fIDs(4)),&
                                                mesh % faces(fIDs(5)),&
                                                mesh % faces(fIDs(6)),&
                                                faces=mesh % faces)  
               end if 
               end associate
            end do
!$omp end do  
         end if
!
!        **********************************************
!        Compute interface solution of non-shared faces
!        **********************************************
!
         if (HOElements) then
!$omp do schedule(runtime) private(fID)
            do iFace = 1, size(mesh % HO_FacesInterior)
               fID = mesh % HO_FacesInterior(iFace)
               if (mesh % faces(fID) % IsMortar==1) then 
                  associate(unStar=>mesh% faces(fID)%storage(1)%unStar)
                     unStar=0.0_RP
                  end associate
                  associate(unStar=>mesh% faces(fID)%storage(2)%unStar)
                     unStar=0.0_RP
                  end associate
                  do m=1,4
                     if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                        call BR2_GradientInterfaceSolution(fma=mesh % faces(fID), nEqn=nEqn, nGradEqn=nGradEqn, GetGradients=GetGradients, &
                        f=mesh % faces(mesh % faces(fID)%Mortar(m)))
                     end if 
                  end do 
               elseif (mesh % faces(fID) % IsMortar==0) then 
               call BR2_GradientInterfaceSolution(mesh % faces(fID), nEqn, nGradEqn, GetGradients)
               end if 
            end do
!$omp end do nowait
         else
!$omp do schedule(runtime) private(fID)
            do iFace = 1, size(mesh % faces_interior)
               fID = mesh % faces_interior(iFace)
               if (mesh % faces(fID) % IsMortar==1) then 
                  associate(unStar=>mesh% faces(fID)%storage(1)%unStar)
                     unStar=0.0_RP
                  end associate
                  associate(unStar=>mesh% faces(fID)%storage(2)%unStar)
                     unStar=0.0_RP
                  end associate
                  do m=1,4
                     if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                        call BR2_GradientInterfaceSolution(fma=mesh % faces(fID), nEqn=nEqn, nGradEqn=nGradEqn, GetGradients=GetGradients, &
                        f=mesh % faces(mesh % faces(fID)%Mortar(m)))
                     end if 
                  end do 
               elseif (mesh % faces(fID) % IsMortar==0) then 
               call BR2_GradientInterfaceSolution(mesh % faces(fID), nEqn, nGradEqn, GetGradients)
               end if 
            end do
!$omp end do nowait
         end if

         if (HOElements) then
!$omp do schedule(runtime) private(fID)
            do iFace = 1, size(mesh % HO_FacesBoundary)
               fID = mesh % HO_FacesBoundary(iFace)
               call BR2_GradientInterfaceSolutionBoundary(mesh % faces(fID), nEqn, nGradEqn, time, GetGradients)
            end do
!$omp end do 
         else
!$omp do schedule(runtime) private(fID)
            do iFace = 1, size(mesh % faces_boundary)
               fID = mesh % faces_boundary(iFace)
               call BR2_GradientInterfaceSolutionBoundary(mesh % faces(fID), nEqn, nGradEqn, time, GetGradients)
            end do
!$omp end do 
         end if
!
!        **********************
!        Compute face integrals
!        **********************
!
         if (HOElements) then
!$omp do schedule(runtime) private(eID) 
            do iEl = 1, size(mesh % HO_ElementsSequential)
               eID = mesh % HO_ElementsSequential(iEl)
               call BR2_ComputeGradientFaceIntegrals(self, nGradEqn, mesh % elements(eID), mesh)
            end do
!$omp end do
         else
!$omp do schedule(runtime) private(eID) 
            do iEl = 1, size(mesh % elements_sequential)
               eID = mesh % elements_sequential(iEl)
               call BR2_ComputeGradientFaceIntegrals(self, nGradEqn, mesh % elements(eID), mesh)
            end do
!$omp end do
         end if
!
!        ******************
!        Wait for MPI faces
!        ******************
!
#ifdef _HAS_MPI_
!$omp single
         if ( MPI_Process % doMPIAction ) then 
            call mesh % GatherMPIFacesSolution(nEqn)
         end if
!$omp end single
!
!        *******************************
!        Compute MPI interface solutions
!        *******************************
!
!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % faces_mpi)
            fID = mesh % faces_mpi(iFace)
            if (mesh% faces(fID)%IsMortar==1) then 
               associate(unstar=>mesh% faces(fID)%storage(1)%unStar)
                  unstar=0.0_RP
               end associate
               do m=1, 4
                  if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                  call BR2_GradientInterfaceSolution(fma=mesh % faces(fID), nEqn=nEqn, nGradEqn=nGradEqn, GetGradients=GetGradients, &
                  f=mesh % faces(mesh % faces(fID)%Mortar(m)))
                  end if 
               end do
            end if 
            call BR2_GradientInterfaceSolutionMPI(mesh % faces(fID), nEqn, nGradEqn, GetGradients)
         end do
!$omp end do 

!$omp single
         if ( mesh % nonconforming ) then
            call mesh % UpdateMPIFacesGradMortarflux(nGradEqn)
         end if
   !$omp end single
   
   
   !$omp single
         if ( mesh % nonconforming ) then
            call mesh % GatherMPIFacesGradMortarFlux(nGradEqn)
         end if
   !$omp end single 
!
!        **************************************************
!        Compute face integrals for elements with MPI faces
!        **************************************************
!
         if (HOElements) then
!$omp do schedule(runtime) private(eID)
            do iEl = 1, size(mesh % HO_ElementsMPI)
               eID = mesh % HO_ElementsMPI(iEl)
               call BR2_ComputeGradientFaceIntegrals(self, nGradEqn, mesh % elements(eID), mesh)
            end do
!$omp end do
         else
!$omp do schedule(runtime) private(eID)
            do iEl = 1, size(mesh % elements_mpi)
               eID = mesh % elements_mpi(iEl)
               call BR2_ComputeGradientFaceIntegrals(self, nGradEqn, mesh % elements(eID), mesh)
            end do
!$omp end do
         end if
#endif

      end subroutine BR2_ComputeGradient
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR2_ComputeGradientFaceIntegrals( self, nGradEqn, e, mesh)
!
!        *******************************************************
!              The surface integrals in the BR2 method considers
!           also the correction of the surface gradients with
!           the stabilizing term:
!              \nabla u^* = \nabla u + \eta r_f([[u]])
!
!           where
!              \int_{e} r_f([[u]])\tau = -0.5\int_{f} [[u]] \tau ds
!
!        *******************************************************
!
         use ElementClass
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use DGIntegrals
         use NodalStorageClass, only: NodalStorage
         implicit none
         class(BassiRebay2_t),   intent(in) :: self
         integer,                intent(in) :: nGradEqn
         class(Element)                     :: e
         class(HexMesh)                     :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i, j, k
         real(kind=RP) :: invjac(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP) :: faceInt_x(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP) :: faceInt_y(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP) :: faceInt_z(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP) :: bv_x(0:e % Nxyz(1),2)
         real(kind=RP) :: bv_y(0:e % Nxyz(2),2)
         real(kind=RP) :: bv_z(0:e % Nxyz(3),2)

         call VectorWeakIntegrals % StdFace(e, nGradEqn, &
               mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % unStar, &
               mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % unStar, &
               mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % unStar, &
               mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % unStar, &
               mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % unStar, &
               mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % unStar, &
               faceInt_x, faceInt_y, faceInt_z )
!
!        Add the integrals weighted with the Jacobian
!        --------------------------------------------
         do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2)    ; do i = 0, e % Nxyz(1)
            invjac(i,j,k) = 1.0_RP / e % geom % jacobian(i,j,k)
            e % storage % U_x(:,i,j,k) = e % storage % U_x(:,i,j,k) - faceInt_x(:,i,j,k) * invjac(i,j,k)
            e % storage % U_y(:,i,j,k) = e % storage % U_y(:,i,j,k) - faceInt_y(:,i,j,k) * invjac(i,j,k)
            e % storage % U_z(:,i,j,k) = e % storage % U_z(:,i,j,k) - faceInt_z(:,i,j,k) * invjac(i,j,k)
         end do                  ; end do                   ; end do
!
!        ******************************************
!        Perform the interface gradients correction
!        ******************************************
!
         associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                   spAeta  => NodalStorage(e % Nxyz(2)), &
                   spAzeta => NodalStorage(e % Nxyz(3)) )
         bv_x = spAxi % b * spAxi % v
         bv_y = spAeta % b * spAeta % v
         bv_z = spAzeta % b * spAzeta % v
         end associate
!
!        ----------------
!>       Xi-contributions
!        ----------------
!

         associate(U_x => mesh % faces(e % faceIDs(ELEFT)) % storage(e % faceSide(ELEFT)) % U_x, &
                   U_y => mesh % faces(e % faceIDs(ELEFT)) % storage(e % faceSide(ELEFT)) % U_y, &
                   U_z => mesh % faces(e % faceIDs(ELEFT)) % storage(e % faceSide(ELEFT)) % U_z, &
                   unStar => mesh % faces(e % faceIDs(ELEFT)) % storage(e % faceSide(ELEFT)) % unStar )

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            U_x(:,j,k) = U_x(:,j,k) - self % eta * unStar(:,1,j,k) * bv_x(i,LEFT) * invjac(i,j,k)
            U_y(:,j,k) = U_y(:,j,k) - self % eta * unStar(:,2,j,k) * bv_x(i,LEFT) * invjac(i,j,k)
            U_z(:,j,k) = U_z(:,j,k) - self % eta * unStar(:,3,j,k) * bv_x(i,LEFT) * invjac(i,j,k)
         end do                 ; end do                ; end do
         end associate

         associate(U_x => mesh % faces(e % faceIDs(ERIGHT)) % storage(e % faceSide(ERIGHT)) % U_x, &
                   U_y => mesh % faces(e % faceIDs(ERIGHT)) % storage(e % faceSide(ERIGHT)) % U_y, &
                   U_z => mesh % faces(e % faceIDs(ERIGHT)) % storage(e % faceSide(ERIGHT)) % U_z, &
                   unStar => mesh % faces(e % faceIDs(ERIGHT)) % storage(e % faceSide(ERIGHT)) % unStar )

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            U_x(:,j,k) = U_x(:,j,k) - self % eta * unStar(:,1,j,k) * bv_x(i,RIGHT) * invjac(i,j,k)
            U_y(:,j,k) = U_y(:,j,k) - self % eta * unStar(:,2,j,k) * bv_x(i,RIGHT) * invjac(i,j,k)
            U_z(:,j,k) = U_z(:,j,k) - self % eta * unStar(:,3,j,k) * bv_x(i,RIGHT) * invjac(i,j,k)
         end do                 ; end do                ; end do
         end associate
!
!        -----------------
!>       Eta-contributions
!        -----------------
!
         associate(U_x => mesh % faces(e % faceIDs(EFRONT)) % storage(e % faceSide(EFRONT)) % U_x, &
                   U_y => mesh % faces(e % faceIDs(EFRONT)) % storage(e % faceSide(EFRONT)) % U_y, &
                   U_z => mesh % faces(e % faceIDs(EFRONT)) % storage(e % faceSide(EFRONT)) % U_z, &
                   unStar => mesh % faces(e % faceIDs(EFRONT)) % storage(e % faceSide(EFRONT)) % unStar )

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            U_x(:,i,k) = U_x(:,i,k) - self % eta * unStar(:,1,i,k) * bv_y(j,LEFT) * invjac(i,j,k)
            U_y(:,i,k) = U_y(:,i,k) - self % eta * unStar(:,2,i,k) * bv_y(j,LEFT) * invjac(i,j,k)
            U_z(:,i,k) = U_z(:,i,k) - self % eta * unStar(:,3,i,k) * bv_y(j,LEFT) * invjac(i,j,k)
         end do                 ; end do                ; end do
         end associate

         associate(U_x => mesh % faces(e % faceIDs(EBACK)) % storage(e % faceSide(EBACK)) % U_x, &
                   U_y => mesh % faces(e % faceIDs(EBACK)) % storage(e % faceSide(EBACK)) % U_y, &
                   U_z => mesh % faces(e % faceIDs(EBACK)) % storage(e % faceSide(EBACK)) % U_z, &
                   unStar => mesh % faces(e % faceIDs(EBACK)) % storage(e % faceSide(EBACK)) % unStar )

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            U_x(:,i,k) = U_x(:,i,k) - self % eta * unStar(:,1,i,k) * bv_y(j,RIGHT) * invjac(i,j,k)
            U_y(:,i,k) = U_y(:,i,k) - self % eta * unStar(:,2,i,k) * bv_y(j,RIGHT) * invjac(i,j,k)
            U_z(:,i,k) = U_z(:,i,k) - self % eta * unStar(:,3,i,k) * bv_y(j,RIGHT) * invjac(i,j,k)
         end do                 ; end do                ; end do
         end associate
!
!        ------------------
!>       Zeta-contributions
!        ------------------
!
         associate(U_x => mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % U_x, &
                   U_y => mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % U_y, &
                   U_z => mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % U_z, &
                   unStar => mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % unStar )

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            U_x(:,i,j) = U_x(:,i,j) - self % eta * unStar(:,1,i,j) * bv_z(k,LEFT) * invjac(i,j,k)
            U_y(:,i,j) = U_y(:,i,j) - self % eta * unStar(:,2,i,j) * bv_z(k,LEFT) * invjac(i,j,k)
            U_z(:,i,j) = U_z(:,i,j) - self % eta * unStar(:,3,i,j) * bv_z(k,LEFT) * invjac(i,j,k)
         end do                 ; end do                ; end do
         end associate

         associate(U_x => mesh % faces(e % faceIDs(ETOP)) % storage(e % faceSide(ETOP)) % U_x, &
                   U_y => mesh % faces(e % faceIDs(ETOP)) % storage(e % faceSide(ETOP)) % U_y, &
                   U_z => mesh % faces(e % faceIDs(ETOP)) % storage(e % faceSide(ETOP)) % U_z, &
                   unStar => mesh % faces(e % faceIDs(ETOP)) % storage(e % faceSide(ETOP)) % unStar )

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            U_x(:,i,j) = U_x(:,i,j) - self % eta * unStar(:,1,i,j) * bv_z(k,RIGHT) * invjac(i,j,k)
            U_y(:,i,j) = U_y(:,i,j) - self % eta * unStar(:,2,i,j) * bv_z(k,RIGHT) * invjac(i,j,k)
            U_z(:,i,j) = U_z(:,i,j) - self % eta * unStar(:,3,i,j) * bv_z(k,RIGHT) * invjac(i,j,k)
         end do                 ; end do                ; end do
         end associate

      end subroutine BR2_ComputeGradientFaceIntegrals
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR2_GradientInterfaceSolution(f, nEqn, nGradEqn, GetGradients, fma)
!
!        ************************************************
!           The BR2 is written in strong form, since it
!           is more efficient for the interpolation to
!           boundaries of the local gradients. Hence,
!           the numerical flux is compensated with the
!           interior solution to yield the interface
!           jumps:
!              U_x += -0.5\int_{S} [[u]]\tau ds
!
!           We compute here the interface fluxes:
!              unStar = -0.5*[[u]]dS
!        ************************************************
!
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         type(Face)                       :: f
         integer, intent(in)              :: nEqn, nGradEqn
         procedure(GetGradientValues_f)   :: GetGradients
         type(Face), optional             :: fma
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: UL(nGradEqn), UR(nGradEqn)
         real(kind=RP) :: Uhat(nGradEqn)
         real(kind=RP) :: Hflux(nGradEqn,NDIM,0:f % Nf(1), 0:f % Nf(2))
         integer       :: i,j
         
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR)
   
            Uhat = 0.5_RP * (UL - UR) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

            if (f % IsMortar==0) then 
         
               call f % ProjectGradientFluxToElements(nGradEqn, HFlux,(/1,2/),1)
            end if 
            if (f % IsMortar==2 .and. present(fma)) then 
                 call fma % ProjectMortarGradientFluxToElements(nEqn=nGradEqn, fma=f, HFlux=HFlux,whichElements=(/1,0/),factor=1) 

                 call f % ProjectGradientFluxToElements(nGradEqn, HFlux,(/0,2/),1)
            end if 
      end subroutine BR2_GradientInterfaceSolution   

      subroutine BR2_GradientInterfaceSolutionMPI(f, nEqn, nGradEqn, GetGradients)
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         type(Face)                       :: f
         integer, intent(in)              :: nEqn, nGradEqn
         procedure(GetGradientValues_f)   :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: UL(nGradEqn), UR(nGradEqn)
         real(kind=RP) :: Uhat(nGradEqn)
         real(kind=RP) :: Hflux(nGradEqn,NDIM,0:f % Nf(1), 0:f % Nf(2))
         integer       :: i,j, thisSide
         
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR)
   
            Uhat = 0.5_RP * (UL - UR) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectGradientFluxToElements(nGradEqn, HFlux,(/thisSide, HMESH_NONE/),1)
         
         if (f % IsMortar==2) then 
            !write(*,*) 'this side', thisSide
            call f% Interpolatesmall2biggrad(nGradEqn, HFlux)
            
         end if 
       
      end subroutine BR2_GradientInterfaceSolutionMPI   

      subroutine BR2_GradientInterfaceSolutionBoundary(f, nEqn, nGradEqn, time, GetGradients)
         use Physics
         use FaceClass
         implicit none
         type(Face)                       :: f
         integer,    intent(in)           :: nEqn
         integer,    intent(in)           :: nGradEqn
         real(kind=RP)                    :: time
         procedure(GetGradientValues_f)   :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i, j
         real(kind=RP) :: Uhat(nGradEqn), UL(nGradEqn), UR(nGradEqn)
         real(kind=RP) :: bvExt(nEqn)

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)

            bvExt =  f % storage(1) % Q(:,i,j)
   
            call BCs(f % zone) % bc % StateForEqn( nEqn, f % geom % x(:,i,j), &
                                time               , &
                                f % geom % normal(:,i,j)      , &
                                bvExt              )
!   
!           -------------------
!           u, v, w, T averages
!           -------------------
!   
            call GetGradients( nEqn, nGradEqn, f % storage(1) % Q(:,i,j), UL )
            call GetGradients( nEqn, nGradEqn, bvExt, UR )
   
            Uhat = 0.5_RP * (UL - UR) * f % geom % jacobian(i,j)
            
            f % storage(1) % unStar(:,1,i,j) = Uhat * f % geom % normal(1,i,j)
            f % storage(1) % unStar(:,2,i,j) = Uhat * f % geom % normal(2,i,j)
            f % storage(1) % unStar(:,3,i,j) = Uhat * f % geom % normal(3,i,j)

         end do ; end do   
         
      end subroutine BR2_GradientInterfaceSolutionBoundary
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR2_ComputeInnerFluxes(self, nEqn, nGradEqn, EllipticFlux, GetViscosity, e, contravariantFlux )
         use ElementClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay2_t) ,     intent (in) :: self
         integer,                   intent(in)  :: nEqn, nGradEqn
         procedure(EllipticFlux_f)              :: EllipticFlux
         procedure(GetViscosity_f)              :: GetViscosity
         type(Element)                          :: e
         real(kind=RP)           , intent (out) :: contravariantFlux(1:nEqn, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)       :: delta
         real(kind=RP)       :: cartesianFlux(1:nEqn, 1:NDIM)
         real(kind=RP)       :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: beta(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: kappa(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         integer             :: i, j, k

#if (!defined(CAHNHILLIARD))

#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
         mu = e % storage % mu_ns(1,:,:,:)
         kappa = e % storage % mu_ns(2,:,:,:)
         beta  = 0.0_RP

#elif defined(NAVIERSTOKES) && (SPALARTALMARAS)
         mu    = e % storage % mu_ns(1,:,:,:)
         kappa = e % storage % mu_ns(2,:,:,:)
         beta  = e % storage % mu_ns(3,:,:,:)
#elif defined(INCNS)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            call GetViscosity(e % storage % Q(INSRHO,i,j,k), mu(i,j,k))      
         end do                ; end do                ; end do

         kappa = 0.0_RP
         beta  = 0.0_RP
#endif

#else /* !(defined(CAHNHILLIARD) */ 

#if defined(NAVIERSTOKES)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            call GetViscosity(e % storage % c(1,i,j,k), mu(i,j,k))      
         end do                ; end do                ; end do
         kappa = 1.0_RP / ( thermodynamics % gammaMinus1 * &
                               POW2( dimensionless % Mach) * dimensionless % Pr ) * mu

         beta  = 0.0_RP
#elif defined(INCNS)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            call GetViscosity(e % storage % c(1,i,j,k), mu(i,j,k))      
         end do                ; end do                ; end do

         kappa = 0.0_RP
         beta  = 0.0_RP
#endif
#endif


         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            call EllipticFlux(nEqn, nGradEqn, e % storage % Q(:,i,j,k), e % storage % U_x(:,i,j,k), &
                              e % storage % U_y(:,i,j,k), e % storage % U_z(:,i,j,k), mu(i,j,k), beta(i,j,k), kappa(i,j,k), cartesianFlux )

            contravariantFlux(:,i,j,k,IX) =     cartesianFlux(:,IX) * e % geom % jGradXi(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradXi(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradXi(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IY) =     cartesianFlux(:,IX) * e % geom % jGradEta(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradEta(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradEta(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IZ) =     cartesianFlux(:,IX) * e % geom % jGradZeta(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradZeta(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)

         end do               ; end do            ; end do

      end subroutine BR2_ComputeInnerFluxes

      subroutine BR2_RiemannSolver ( self , nEqn, nGradEqn, EllipticFlux, f, QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
                                           mu_left, mu_right, nHat , dWall, &
#ifdef MULTIPHASE
sigma, & 
#endif
flux )
         use SMConstants
         use PhysicsStorage
         use Physics
         use FaceClass
         implicit none
         class(BassiRebay2_t)            :: self
         integer,       intent(in)       :: nEqn
         integer,       intent(in)       :: nGradEqn
         procedure(EllipticFlux_f)       :: EllipticFlux
         class(Face),   intent(in)       :: f
         real(kind=RP), intent(in)       :: QLeft(nEqn)
         real(kind=RP), intent(in)       :: QRight(nEqn)
         real(kind=RP), intent(in)       :: U_xLeft(nGradEqn)
         real(kind=RP), intent(in)       :: U_yLeft(nGradEqn)
         real(kind=RP), intent(in)       :: U_zLeft(nGradEqn)
         real(kind=RP), intent(in)       :: U_xRight(nGradEqn)
         real(kind=RP), intent(in)       :: U_yRight(nGradEqn)
         real(kind=RP), intent(in)       :: U_zRight(nGradEqn)
         real(kind=RP), intent(in)       :: mu_left(3)
         real(kind=RP), intent(in)       :: mu_right(3)
         real(kind=RP), intent(in)       :: nHat(NDIM)
         real(kind=RP), intent(in)       :: dWall
#ifdef MULTIPHASE
         real(kind=RP), intent(in)       :: sigma(nEqn)
#endif
         real(kind=RP), intent(out)      :: flux(nEqn)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: Q(nEqn) , U_x(nGradEqn) , U_y(nGradEqn) , U_z(nGradEqn)
         real(kind=RP)     :: flux_vec(nEqn,NDIM), fL(nEqn,NDIM), fR(nEqn,NDIM)
         real(kind=RP)     :: sigma0

         call EllipticFlux(nEqn, nGradEqn, QLeft, U_xLeft, U_yLeft, U_zLeft, mu_left(1), mu_left(2), mu_left(3), fL)
         call EllipticFlux(nEqn, nGradEqn, QRight, U_xRight, U_yRight, U_zRight, mu_right(1), mu_right(2), mu_right(3), fR)

         flux_vec = 0.5_RP * (fL + fR)

         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ) 

#ifdef MULTIPHASE
         sigma0 = 0.5_RP * self % sigma * (maxval(f % Nf))*(maxval(f % Nf)+1) / f % geom % h
         flux = flux - sigma0 * sigma * (QLeft-QRight)
#endif

      end subroutine BR2_RiemannSolver
end module EllipticBR2