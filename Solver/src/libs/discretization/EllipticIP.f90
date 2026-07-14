#include "Includes.h"
module EllipticIP
   use SMConstants
   use Headers
   use MeshTypes
   use ElementClass
   use HexMeshClass
   use Physics
   use PhysicsStorage
   use VariableConversion
   use MPI_Process_Info
   use MPI_Face_Class
   use EllipticDiscretizationClass
   use DGSEMClass
   use FluidData
   use BoundaryConditions           , only: BCs
   use Utilities                    , only: toLower, dot_product
   implicit none
!
!
   private
   public   InteriorPenalty_t, SIPG, IIPG, NIPG
   public   IP_GradientInterfaceSolution, IP_GradientInterfaceSolutionBoundary
   public   IP_GradientInterfaceSolutionMPI

   integer, parameter   :: SIPG = -1
   integer, parameter   :: IIPG = 0
   integer, parameter   :: NIPG = 1

   type, extends(EllipticDiscretization_t)   :: InteriorPenalty_t
      procedure(PenaltyParameter_f), pointer   :: PenaltyParameter
      integer              :: IPmethod = SIPG
      contains
         procedure      :: Construct               => IP_Construct
         procedure      :: ComputeGradient         => IP_ComputeGradient
         procedure      :: ComputeInnerFluxes      => IP_ComputeInnerFluxes
         procedure      :: RiemannSolver           => IP_RiemannSolver
#if defined(NAVIERSTOKES) && !(SPALARTALMARAS)
         procedure      :: RiemannSolver_Jacobians => IP_RiemannSolver_Jacobians
#endif
         procedure      :: Describe                => IP_Describe
   end type InteriorPenalty_t

   abstract interface
      function PenaltyParameter_f(self, f)
         use SMConstants
         use FaceClass
         import InteriorPenalty_t
         implicit none
         class(InteriorPenalty_t)   :: self
         class(Face), intent(in)    :: f
         real(kind=RP)              :: PenaltyParameter_f
      end function PenaltyParameter_f
   end interface
!
!  ========
   contains
!  ========
!
      subroutine IP_Construct(self, controlVariables, eqname)
         use FTValueDictionaryClass
         use mainKeywordsModule
         use MPI_Process_Info
         use PhysicsStorage
         implicit none
         class(InteriorPenalty_t)              :: self
         class(FTValueDictionary), intent(in)  :: controlVariables
         integer,                  intent(in)  :: eqname
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)            :: eqnameaux
         character(len=LINE_LENGTH)            :: IPvariant
!
!        ----------------------------------------------------------
!        Set the particular procedures to compute the elliptic flux
!        ----------------------------------------------------------
!
         select case (eqName)
         case(ELLIPTIC_NS)
            self % eqName = ELLIPTIC_NS
            self % PenaltyParameter => PenaltyParameterNS

         case(ELLIPTIC_NSSA)
            self % eqName = ELLIPTIC_NSSA
            self % PenaltyParameter => PenaltyParameterNS

         case(ELLIPTIC_iNS)
            self % eqName = ELLIPTIC_iNS
            self % PenaltyParameter => PenaltyParameterNS

         case(ELLIPTIC_CH)
            self % eqName = ELLIPTIC_CH
            self % PenaltyParameter => PenaltyParameterCH

         case(ELLIPTIC_MU)
            self % eqName = ELLIPTIC_MU
            self % PenaltyParameter => PenaltyParameterNS

         case default
            print*, "Unrecognized equation"
            errorMessage(STD_OUT)
            error stop

         end select
!
!        Request the penalty parameter
!        -----------------------------
         if ( controlVariables % containsKey("penalty parameter") ) then
            self % sigma = controlVariables % doublePrecisionValueForKey("penalty parameter")

         else
!            
!           Set 1.0 by default
!           ------------------
            self % sigma = 1.0_RP

         end if
!
!        Request the interior penalty variant
!        ------------------------------------
         if ( controlVariables % containsKey("interior penalty variant") ) then
            IPvariant = controlVariables % stringValueForKey("interior penalty variant", LINE_LENGTH)
            call toLower(IPVariant)
      
         else
!
!           Select SIPG by default
!           ----------------------
            IPvariant = "sipg"
   
         end if

         select case (trim(IPvariant))
         case("sipg")
            self % IPmethod = SIPG

         case("iipg")
            self % IPmethod = IIPG

         case("nipg")
            self % IPmethod = NIPG

         case default
            if ( MPI_Process % isRoot ) then
            print*, "Unknown selected IP variant", trim(IPvariant), "."
            print*, "Available options are:"
            print*, "   * SIPG"
            print*, "   * IIPG"
            print*, "   * NIPG"
            errorMessage(STD_OUT)
            error stop
            end if
         end select

      end subroutine IP_Construct

      subroutine IP_Describe(self)
         implicit none
         class(InteriorPenalty_t),  intent(in)  :: self
!
!        Display the configuration
!        -------------------------
         if (MPI_Process % isRoot) write(STD_OUT,'(/)')
         call Subsection_Header("Viscous discretization")

         if (.not. MPI_Process % isRoot ) return

         write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","IP"

         select case(self % IPmethod)
         case(SIPG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Interior penalty variant: ","SIPG"

         case(NIPG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Interior penalty variant: ","NIPG"

         case(IIPG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Interior penalty variant: ","IIPG"
            
         end select

         write(STD_OUT,'(30X,A,A30,F10.3)') "->","Penalty parameter: ", self % sigma
            
      end subroutine IP_Describe

      subroutine IP_ComputeGradient(self, nEqn, nGradEqn, mesh, time, GetGradients, HO_Elements, element_mask, Level)
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use MPI_Process_Info
         implicit none
         class(InteriorPenalty_t), intent(in) :: self
         integer,                  intent(in) :: nEqn, nGradEqn
         class(HexMesh)                       :: mesh
         real(kind=RP),        intent(in)     :: time
         procedure(GetGradientValues_f)       :: GetGradients 
         logical, intent(in), optional        :: HO_Elements
         logical, intent(in), optional        :: element_mask(:)
		 integer, intent(in), optional        :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: Nx, Ny, Nz
         integer :: i, j, k, m, ii, jj
         integer :: eID , fID , dimID , eqID, fIDs(6), iFace, iEl
         logical :: HOElements
         real(kind=RP)   :: hStarAux(nGradEqn, NDIM, 0:mesh%faces(1)% NfRight(1), 0:mesh%faces(1)% NfRight(1))

         if (present(HO_Elements)) then
            HOElements = HO_Elements
         else
            HOElements = .false.
         end if
!
!        *********************************
!        Volume loops and prolong to faces
!        *********************************
!
         if (HOElements) then
!$omp do schedule(runtime)
            do eID = 1, size(mesh % HO_Elements)
               associate( e => mesh % elements(mesh % HO_Elements(eID)) )
               call e % ComputeLocalGradient(nEqn, nGradEqn, GetGradients, .false.)
   !
   !           Prolong to faces
   !           ----------------
               fIDs = e % faceIDs

            if (.not.mesh%SlidingMesh%active) then 
               if (.not.mesh%nonconforming) then 
               call e % ProlongGradientsToFaces(nGradEqn, mesh % faces(fIDs(1)),&
                                                mesh % faces(fIDs(2)),&
                                                mesh % faces(fIDs(3)),&
                                                mesh % faces(fIDs(4)),&
                                                mesh % faces(fIDs(5)),&
                                             mesh % faces(fIDs(6)))
               else 
               call e % ProlongGradientsToFaces(nGradEqn, mesh % faces(fIDs(1)),&
                                             mesh % faces(fIDs(2)),&
                                             mesh % faces(fIDs(3)),&
                                             mesh % faces(fIDs(4)),&
                                             mesh % faces(fIDs(5)),&
                                             mesh % faces(fIDs(6)),&
                                             faces=mesh % faces)
               end if 
            else 
               call e %  ProlongGradientsToFaces(nGradEqn, &
                                             fFR=mesh % faces(fIDs(1)),&
                                             fBK=mesh % faces(fIDs(2)),&
                                             fBOT=mesh % faces(fIDs(3)),&
                                             fR=mesh % faces(fIDs(4)),&
                                             fT=mesh % faces(fIDs(5)),&
                                             fL=mesh % faces(fIDs(6)),&
                                             faces=mesh % mortar_faces )
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
               if (.not.mesh%SlidingMesh%active) then 
                  if (.not.mesh%nonconforming) then 
                  call e % ProlongGradientsToFaces(nGradEqn, mesh % faces(fIDs(1)),&
                                                mesh % faces(fIDs(2)),&
                                                mesh % faces(fIDs(3)),&
                                                mesh % faces(fIDs(4)),&
                                                mesh % faces(fIDs(5)),&
                                                mesh % faces(fIDs(6)))
                  else 
                  call e % ProlongGradientsToFaces(nGradEqn, mesh % faces(fIDs(1)),&
                                                mesh % faces(fIDs(2)),&
                                                mesh % faces(fIDs(3)),&
                                                mesh % faces(fIDs(4)),&
                                                mesh % faces(fIDs(5)),&
                                                mesh % faces(fIDs(6)),&
                                                faces=mesh % faces)
                  end if 
               else 
                  call e %  ProlongGradientsToFaces(nGradEqn, &
                                                fFR=mesh % faces(fIDs(1)),&
                                                fBK=mesh % faces(fIDs(2)),&
                                                fBOT=mesh % faces(fIDs(3)),&
                                                fR=mesh % faces(fIDs(4)),&
                                                fT=mesh % faces(fIDs(5)),&
                                                fL=mesh % faces(fIDs(6)),&
                                                faces=mesh % mortar_faces )
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
               if (mesh % faces(fID) % MortarType == MORTAR_SLIDING) then 
                  associate(unStar=>mesh% faces(fID)%storage(1)%unStar)
                     unStar=0.0_RP
                  end associate
                  associate(unStar=>mesh% faces(fID)%storage(2)%unStar)
                     unStar=0.0_RP
                  end associate 
               end if 
               if (mesh % faces(fID) % MortarType == MORTAR_BIG) then 
                  associate(unStar=>mesh% faces(fID)%storage(1)%unStar)
                     unStar=0.0_RP
                  end associate
                  associate(unStar=>mesh% faces(fID)%storage(2)%unStar)
                     unStar=0.0_RP
                  end associate
                  do m=1,4
                     if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                     call IP_GradientInterfaceSolution(masterFace1=mesh % faces(fID), nEqn=nEqn, nGradEqn=nGradEqn, GetGradients=GetGradients, &
                     f=mesh % faces(mesh % faces(fID)%Mortar(m)))
                     end if 
                  end do 
               elseif(mesh % faces(fID) % MortarType == MORTAR_NONE) then
               call IP_GradientInterfaceSolution(mesh % faces(fID), nEqn, nGradEqn, GetGradients)
               end if 
            end do
!$omp end do 
         else
!$omp do schedule(runtime) private(fID)
            do iFace = 1, size(mesh % faces_interior)
               fID = mesh % faces_interior(iFace)
               if (mesh % faces(fID) % MortarType == MORTAR_SLIDING) then 
                  associate(unStar=>mesh% faces(fID)%storage(1)%unStar)
                     unStar=0.0_RP
                  end associate
                  associate(unStar=>mesh% faces(fID)%storage(2)%unStar)
                     unStar=0.0_RP
                  end associate 
               end if 
               if (mesh % faces(fID) % MortarType == MORTAR_BIG) then 
                  associate(unStar=>mesh% faces(fID)%storage(1)%unStar)
                     unStar=0.0_RP
                  end associate
                  associate(unStar=>mesh% faces(fID)%storage(2)%unStar)
                     unStar=0.0_RP
                  end associate
                  do m=1,4
                     if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                     call IP_GradientInterfaceSolution(masterFace1=mesh % faces(fID), nEqn=nEqn, nGradEqn=nGradEqn, GetGradients=GetGradients, &
                     f=mesh % faces(mesh % faces(fID)%Mortar(m)))
                     end if 
                  end do 
               elseif(mesh % faces(fID) % MortarType == MORTAR_NONE) then
               call IP_GradientInterfaceSolution(mesh % faces(fID), nEqn, nGradEqn, GetGradients)
               end if 
            end do
!$omp end do 
         end if

      if (mesh%SlidingMesh%active) then 
!$omp do schedule(runtime) private(fID)
                  do iFace = 1, size(mesh % mortar_faces)
                     fID = mesh % mortar_faces(iFace)%ID
                     associate(unStar=>mesh % faces(mesh % mortar_faces(fID)%Mortar(1))%storage(1)%unStar)
                        unStar=0.0_RP
                     end associate
                     associate(unStar=>mesh % faces(mesh % mortar_faces(fID)%Mortar(1))%storage(2)%unStar)
                        unStar=0.0_RP
                     end associate

                     associate(unStar=>mesh % faces(mesh % mortar_faces(fID)%Mortar(2))%storage(1)%unStar)
                        unStar=0.0_RP
                     end associate
                     associate(unStar=>mesh % faces(mesh % mortar_faces(fID)%Mortar(2))%storage(2)%unStar)
                        unStar=0.0_RP
                     end associate
                  end do 
!$omp end do        
               end if 

      if (mesh%SlidingMesh%active) then 
!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % mortar_faces)
            fID = mesh % mortar_faces(iFace)%ID
            call IP_GradientInterfaceSolution(f=mesh % mortar_faces(fID), nEqn=nEqn, nGradEqn=nGradEqn, GetGradients=GetGradients,&
            masterFace1=mesh % faces (mesh % mortar_faces(fID)%Mortar(1)), masterFace2=mesh % faces (mesh % mortar_faces(fID)%Mortar(2)), sliding=.true.)
            
         end do 
!$omp end do        
      end if 


         if (HOElements) then
!$omp do schedule(runtime) private(fID)
            do iFace = 1, size(mesh % HO_FacesBoundary)
               fID = mesh % HO_FacesBoundary(iFace)
               call IP_GradientInterfaceSolutionBoundary(mesh % faces(fID), nEqn, nGradEqn, time, GetGradients)
            end do
!$omp end do 
         else
!$omp do schedule(runtime) private(fID)
            do iFace = 1, size(mesh % faces_boundary)
               fID = mesh % faces_boundary(iFace)
               call IP_GradientInterfaceSolutionBoundary(mesh % faces(fID), nEqn, nGradEqn, time, GetGradients)
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
               call IP_ComputeGradientFaceIntegrals(self,nGradEqn, mesh % elements(eID), mesh)
            end do
!$omp end do
         else
!$omp do schedule(runtime) private(eID) 
            do iEl = 1, size(mesh % elements_sequential)
               eID = mesh % elements_sequential(iEl)
               call IP_ComputeGradientFaceIntegrals(self,nGradEqn, mesh % elements(eID), mesh)
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
            if (mesh% faces(fID)%MortarType == MORTAR_BIG) then 
               associate(UnStar=>mesh% faces(fID)%storage(1)%UnStar)
                  UnStar=0.0_RP
               end associate
               do m=1,4
                  if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                     call IP_GradientInterfaceSolution(masterFace1=mesh % faces(fID), nEqn=nEqn, nGradEqn=nGradEqn, GetGradients=GetGradients, &
                     f=mesh % faces(mesh % faces(fID)%Mortar(m)))
                  end if 
               end do
            end if 
            call IP_GradientInterfaceSolutionMPI(mesh % faces(fID), nEqn, nGradEqn, GetGradients)
         end do
!$omp end do 

!$omp single
         if ( mesh % nonconforming ) then
            call mesh % UpdateMPIFacesGradMortarflux(NCONS)
         end if
!$omp end single
   
   
!$omp single
         if ( mesh % nonconforming ) then
            call mesh % GatherMPIFacesGradMortarFlux(NCONS)
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
               call IP_ComputeGradientFaceIntegrals(self,nGradEqn, mesh % elements(eID), mesh)
            end do
!$omp end do
         else
!$omp do schedule(runtime) private(eID)
            do iEl = 1, size(mesh % elements_mpi)
               eID = mesh % elements_mpi(iEl)
               call IP_ComputeGradientFaceIntegrals(self,nGradEqn, mesh % elements(eID), mesh)
            end do
!$omp end do
         end if
#endif

      end subroutine IP_ComputeGradient
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine IP_ComputeGradientFaceIntegrals(self,nGradEqn, e, mesh)
         use ElementClass
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use DGIntegrals
         implicit none
         type(InteriorPenalty_t),         intent(in) :: self
         integer,                    intent(in) :: nGradEqn
         class(Element)                         :: e
         class(HexMesh)                         :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: i, j, k
         real(kind=RP)        :: invjac
         real(kind=RP)        :: faceInt_x(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP)        :: faceInt_y(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP)        :: faceInt_z(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )

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
            invjac = self % IPmethod * e % geom % invJacobian(i,j,k)
            e % storage % U_x(:,i,j,k) = e % storage % U_x(:,i,j,k) + faceInt_x(:,i,j,k) * invjac
            e % storage % U_y(:,i,j,k) = e % storage % U_y(:,i,j,k) + faceInt_y(:,i,j,k) * invjac
            e % storage % U_z(:,i,j,k) = e % storage % U_z(:,i,j,k) + faceInt_z(:,i,j,k) * invjac
         end do                  ; end do                   ; end do
!
      end subroutine IP_ComputeGradientFaceIntegrals
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine IP_GradientInterfaceSolution(f, nEqn, nGradEqn, GetGradients, masterFace1, masterFace2, sliding)
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
         type(Face), optional             :: masterFace1
         type(Face), optional             :: masterFace2
         logical, optional                :: sliding 

!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: UL(nGradEqn), UR(nGradEqn)
         real(kind=RP) :: Uhat(nGradEqn)
         real(kind=RP) :: Hflux(nGradEqn,NDIM,0:f % Nf(1), 0:f % Nf(2))

         integer       :: i,j, lm

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
#ifdef MULTIPHASE
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL, rho_ = f % storage(1) % rho(i,j))
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR, rho_ = f % storage(2) % rho(i,j))
#else
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR)
#endif

#ifdef MULTIPHASE
!           The multiphase solver needs the Chemical potential as first entropy variable
!           ----------------------------------------------------------------------------
            UL(IGMU) = f % storage(1) % mu(1,i,j)
            UR(IGMU) = f % storage(2) % mu(1,i,j)
#endif


            Uhat = 0.5_RP * (UL - UR) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

      if (.not.present(sliding)) then 
         if (f % MortarType == MORTAR_NONE) then 
         call f % ProjectGradientFluxToElements(nGradEqn, HFlux,(/1,2/),1)
         end if 
         if (f % MortarType == MORTAR_SMALL4 .and. present(masterFace1)) then 
            call masterFace1 % ProjectMortarGradientFluxToElements(nEqn=nGradEqn, slaveFace=f, HFlux=HFlux,whichElements=(/0,2/),factor=1)
            call f % ProjectGradientFluxToElements(nGradEqn, HFlux,(/0,2/),1)
         end if 
      else 
         call masterFace1 % ProjectMortarGradientFluxToElements(nEqn=nGradEqn, slaveFace=f, Hflux=HFlux,whichElements=(/1,0/),factor=1, sliding=sliding) 

         call masterFace2 % ProjectMortarGradientFluxToElements(nEqn=nGradEqn, slaveFace=f, Hflux=HFlux,whichElements=(/2,0/),factor=1,sliding=sliding) 
      end if 
         
      end subroutine IP_GradientInterfaceSolution   

      subroutine IP_GradientInterfaceSolutionMPI(f, nEqn, nGradEqn, GetGradients)
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
         integer,    intent(in)           :: nEqn, nGradEqn
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
#ifdef MULTIPHASE
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL, rho_ = f % storage(1) % rho(i,j))
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR, rho_ = f % storage(2) % rho(i,j))
#else
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR)
#endif

#ifdef MULTIPHASE
!           The multiphase solver needs the Chemical potential as first entropy variable
!           ----------------------------------------------------------------------------
            UL(IGMU) = f % storage(1) % mu(1,i,j)
            UR(IGMU) = f % storage(2) % mu(1,i,j)
#endif

   
            Uhat = 0.5_RP * (UL - UR) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectGradientFluxToElements(nGradEqn, HFlux,(/thisSide, HMESH_NONE/),1)
         if (f % MortarType == MORTAR_SMALL4) then 
         
            call f% Interpolatesmall2biggrad(NCONS, HFlux)
            
         end if 
      end subroutine IP_GradientInterfaceSolutionMPI   

      subroutine IP_GradientInterfaceSolutionBoundary(f, nEqn, nGradEqn, time, GetGradients)
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

#if defined(INCNS) || defined(MULTIPHASE)
         if ( trim(BCs(f % zone) % bc % BCType) /= "freeslipwall" ) then
#endif

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
#ifdef MULTIPHASE
            call GetGradients(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), UL, f % storage(1) % rho(i,j))
            call GetGradients(nEqn, nGradEqn, bvExt, UR, f % storage(1) % rho(i,j))
#else
            call GetGradients(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), UL)
            call GetGradients(nEqn, nGradEqn, bvExt, UR)
#endif

#ifdef MULTIPHASE
!           The multiphase solver needs the Chemical potential as first entropy variable
!           ----------------------------------------------------------------------------
            UL(IGMU) = f % storage(1) % mu(1,i,j)
            UR(IGMU) = f % storage(1) % mu(1,i,j)
#endif

   
            Uhat = 0.5_RP * (UL - UR) * f % geom % jacobian(i,j)
            
            f % storage(1) % unStar(:,1,i,j) = Uhat * f % geom % normal(1,i,j)
            f % storage(1) % unStar(:,2,i,j) = Uhat * f % geom % normal(2,i,j)
            f % storage(1) % unStar(:,3,i,j) = Uhat * f % geom % normal(3,i,j)

         end do ; end do   

#if defined(INCNS) || defined(MULTIPHASE)
         else
!
!           *****************************
!           Set W* = W in free slip walls: [[W]] = 0!
!           *****************************

            f % storage(1) % unStar = 0.0_RP 
         end if 
#endif

         
      end subroutine IP_GradientInterfaceSolutionBoundary
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine IP_ComputeInnerFluxes( self , nEqn, nGradEqn, EllipticFlux, GetViscosity, e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         use Physics
         implicit none
         class(InteriorPenalty_t) ,     intent(in)  :: self
         integer,                       intent(in)  :: nEqn, nGradEqn
         procedure(EllipticFlux_f)                  :: EllipticFlux
         procedure(GetViscosity_f)                  :: GetViscosity
         type(Element)                              :: e
         real(kind=RP)           , intent (out)     :: contravariantFlux(1:nEqn, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
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
            call EllipticFlux( nEqn, nGradEqn, e % storage % Q(:,i,j,k), e % storage % U_x(:,i,j,k), &
                               e % storage % U_y(:,i,j,k) , e % storage % U_z(:,i,j,k), mu(i,j,k), beta(i,j,k), kappa(i,j,k), cartesianFlux )
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

      end subroutine IP_ComputeInnerFluxes
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------
!     Function to get the IP penalty parameter
!     ----------------------------------------
      function PenaltyParameterNS(self, f)
         use FaceClass
         implicit none
         class(InteriorPenalty_t)   :: self
         class(Face), intent(in)    :: f
         real(kind=RP)              :: PenaltyParameterNS

         PenaltyParameterNS = 0.5_RP*self % sigma * (maxval(f % Nf)+1)*(maxval(f % Nf)+2) / f % geom % h 

      end function PenaltyParameterNS

      function PenaltyParameterCH(self, f)
         use FaceClass
         implicit none
         class(InteriorPenalty_t)   :: self
         class(Face), intent(in)    :: f
         real(kind=RP)              :: PenaltyParameterCH

         PenaltyParameterCH = 0.5_RP*self % sigma * (maxval(f % Nf)+1)*(maxval(f % Nf)+2) / f % geom % h 
         !PenaltyParameterCH = 0.25_RP * self % sigma * (maxval(f % Nf))  *(maxval(f % Nf)+1) / f % geom % h 

      end function PenaltyParameterCH

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine IP_RiemannSolver ( self , nEqn, nGradEqn, EllipticFlux, f, QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
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
         class(InteriorPenalty_t)        :: self
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
         real(kind=RP), intent(in)       :: mu_left(3), mu_right(3)
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
         real(kind=RP)     :: flux_vec(nEqn,NDIM)
         real(kind=RP)     :: flux_vecL(nEqn,NDIM)
         real(kind=RP)     :: flux_vecR(nEqn,NDIM)
         real(kind=RP)     :: delta, mu
         

         call EllipticFlux(nEqn, nGradEqn, QLeft , U_xLeft , U_yLeft , U_zLeft, mu_left(1), mu_left(2), mu_left(3), flux_vecL )
         call EllipticFlux(nEqn, nGradEqn, QRight , U_xRight , U_yRight , U_zRight, mu_right(1), mu_right(2), mu_right(3), flux_vecR )

         flux_vec = 0.5_RP * (flux_vecL + flux_vecR)

#ifdef NAVIERSTOKES
         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ) - self % PenaltyParameter(f) * dimensionless % mu * (QLeft - QRight)
#else
         mu = 0.5_RP*(mu_left(1) + mu_right(1))
         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ) - self % PenaltyParameter(f) * mu * (QLeft - QRight)
#endif

      end subroutine IP_RiemannSolver
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------------------
!     Subroutine to get the Jacobian (with respect to ∇q⁺ and ∇q⁻) of the numerical
!     contravariant fluxes on a face. Stored in:
!     -> f % storage(side) % dFv_dGradQF(:,:,:,:,i,j)
!                    |                   |_| | | |_|
!                    |                    |  | |  | 
!                    |                    |  | |  |__Coordinate indexes in face 
!                    |                    |  | |_____1 for inner term, 2 for outer term
!                    |                    |  |_______∇q component: 1, 2, 3
!                    |                    |__________Jacobian for this component
!                    |_______________________________1 for ∇q⁺ and 2 for ∇q⁻
!     -----------------------------------------------------------------------------
#if defined(NAVIERSTOKES) && !(SPALARTALMARAS)
      subroutine IP_RiemannSolver_Jacobians( self, f) 
         use FaceClass
         use Physics
         use PhysicsStorage
         implicit none
         !--------------------------------------------
         class(InteriorPenalty_t), intent(in)    :: self
         type(Face)              , intent(inout) :: f
         !--------------------------------------------
         real(kind=RP), DIMENSION(NCONS,NCONS,NDIM,NDIM) :: df_dgradq   ! Cartesian Jacobian tensor
         real(kind=RP), DIMENSION(NCONS,NCONS,NDIM)      :: dfdq_
         real(kind=RP), parameter :: SideSign(2) = (/ 1._RP, -1._RP /)
         real(kind=RP) :: mu, sigma
         integer :: i,j    ! Face coordinate counters
         integer :: n, m ! Index of G_xx
         integer :: side
         !--------------------------------------------
!
!        Initializations
!        ---------------
         mu    = dimensionless % mu             ! TODO: change for Cahn-Hilliard
         sigma = self % PenaltyParameter(f)
         sigma = sigma * mu
         
         do side = 1, 2
            do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
!
!           ************************************************
!           Jacobian with respect to ∇q: dF/d∇q⁺ and dF/d∇q⁻
!           ************************************************
!
!            
!              Definitions
!              -----------
               associate( Q             => f % storage(side) % Q  (:,i,j)                , &
                          U_x           => f % storage(side) % U_x(:,i,j)                , &
                          U_y           => f % storage(side) % U_y(:,i,j)                , &
                          U_z           => f % storage(side) % U_z(:,i,j)                , &
                          nHat          => f % geom % normal(:,i,j)                      , &
                          dFStar_dq     => f % storage(side) % dFStar_dqF(:,:,i,j)       , &
                          dF_dGradQ_out => f % storage(side) % dFv_dGradQF(:,:,:,i,j) )
               
               call ViscousJacobian(Q, U_x, U_y, U_z, df_dgradq, dfdq_)
               
!               
!            For the outer surface integral
!            ******************************
!            
!              Construct face point Jacobians
!              ------------------------------
               dF_dGradQ_out = 0._RP
               do m = 1, NDIM ; do n = 1, NDIM
                  dF_dGradQ_out(:,:,1) = dF_dGradQ_out(:,:,1) + df_dgradq(:,:,n,m) * f % geom % GradXi  (n,i,j) * nHat(m)
                  dF_dGradQ_out(:,:,2) = dF_dGradQ_out(:,:,2) + df_dgradq(:,:,n,m) * f % geom % GradEta (n,i,j) * nHat(m)
                  dF_dGradQ_out(:,:,3) = dF_dGradQ_out(:,:,3) + df_dgradq(:,:,n,m) * f % geom % GradZeta(n,i,j) * nHat(m)
               end do          ; end do
               
!              Multiply by 1/2 (IP scheme) and the jacobian (surface integral) 
!              ---------------------------------------------------------------
               dF_dGradQ_out = dF_dGradQ_out * 0.5_RP * f % geom % jacobian(i,j)
               dFStar_dq   = dFStar_dq - 0.5_RP * f % geom % jacobian(i,j) * dot_product( dfdq_, nHat )
               
!
!           *********************************************
!           Jacobian with respect to q: dF/dq⁺ and dF/dq⁻
!           *********************************************
!
!
!              Penalty contribution (shifts dFStar_dq matrix)
!              ----------------------------------------------
            
               do n = 1, NCONS
                  dFStar_dq(n,n) = dFStar_dq(n,n) + SideSign(side) * sigma * f % geom % jacobian(i,j)
               end do
               end associate
               
            end do              ; end do
            
         end do
         
      end subroutine IP_RiemannSolver_Jacobians
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module EllipticIP