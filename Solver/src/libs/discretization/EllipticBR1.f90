#include "Includes.h"
module EllipticBR1
   use SMConstants
   use Headers
   use MeshTypes
   use Physics
   use VariableConversion
   use PhysicsStorage
   use MPI_Process_Info
   use MPI_Face_Class
   use EllipticDiscretizationClass
   use FluidData
   use BoundaryConditions, only: BCs
   implicit none
!
!
   private
   public   BassiRebay1_t

   type, extends(EllipticDiscretization_t)   :: BassiRebay1_t
      contains
         procedure      :: ComputeGradient           => BR1_ComputeGradient
         procedure      :: LiftGradients             => BR1_LiftGradients
         procedure      :: LiftGradientsHO           => BR1_LiftGradientsHO
         procedure      :: ComputeInnerFluxes        => BR1_ComputeInnerFluxes
         procedure      :: RiemannSolver             => BR1_RiemannSolver
         procedure      :: Describe => BR1_Describe
   end type BassiRebay1_t
!
!  ========
   contains
!  ========
!
      subroutine BR1_Describe(self)
         implicit none
         class(BassiRebay1_t),   intent(in)  :: self  
!
!        Display the configuration
!        -------------------------
         if (MPI_Process % isRoot) write(STD_OUT,'(/)')

         select case (self % eqName)
         case (ELLIPTIC_NS,ELLIPTIC_NSSA,ELLIPTIC_iNS,ELLIPTIC_MU)
            call Subsection_Header("Viscous discretization")
      
         case (ELLIPTIC_CH)
            call Subsection_Header("Cahn--Hilliard discretization")

         end select
   

         if (.not. MPI_Process % isRoot ) return


         write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","BR1"
         write(STD_OUT,'(30X,A,A30,F6.3)') "->","Penalty parameter: ",self % sigma

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

      end subroutine BR1_Describe

      subroutine BR1_ComputeGradient(self, nEqn, nGradEqn, mesh, time, GetGradients, HO_Elements, element_mask, Level)
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t), intent(in) :: self
         integer,              intent(in) :: nEqn, nGradEqn
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
         integer                :: i, j, k
         integer                :: eID , fID , dimID , eqID, fIDs(6), iFace, iEl, locLevel
         logical                :: set_mu
         logical                :: HOElements
         logical                :: compute_element

         if (present(HO_Elements)) then
            HOElements = HO_Elements
         else
            HOElements = .false.
         end if
!
!        ************
!        Volume loops
!        ************
!
#ifdef MULTIPHASE
         select case (self % eqName)
         case(ELLIPTIC_MU)
            set_mu = .true.
         case default
            set_mu = .false.
         end select
#else
         set_mu = .false.
#endif

      if (HOElements) then
!$omp do schedule(runtime) private(eID)
         do i = 1 , size(mesh % HO_Elements)
            eID = mesh % HO_Elements(i)
            compute_element = .true.
            if (present(element_mask)) compute_element = element_mask(eID)
            
            if (compute_element) then
               call mesh % elements(eID) % ComputeLocalGradient(nEqn, nGradEqn, GetGradients, set_mu)
            endif
         end do
!$omp end do nowait
         call self % LiftGradientsHO(nEqn, nGradEqn, mesh, time, GetGradients, element_mask)
      else
	  	 if (present(Level)) then
			 locLevel=Level
!$omp do schedule(runtime) private(eID)
			 do iEl = 1 , mesh % MLRK % MLIter(locLevel,8) 
				eID = mesh % MLRK % MLIter_eIDN(iEl)
				compute_element = .true.
				if (present(element_mask)) compute_element = element_mask(eID)
				
				if (compute_element) then
				   call mesh % elements(eID) % ComputeLocalGradient(nEqn, nGradEqn, GetGradients, set_mu)
				endif
			 end do
!$omp end do nowait
             call self % LiftGradients(nEqn, nGradEqn, mesh, time, GetGradients, element_mask, locLevel)
         else
!$omp do schedule(runtime) private(eID)
			 do iEl = 1 , size(mesh % elements)
				compute_element = .true.
				if (present(element_mask)) compute_element = element_mask(iEl)
				
				if (compute_element) then
				   call mesh % elements(iEl) % ComputeLocalGradient(nEqn, nGradEqn, GetGradients, set_mu)
				endif
			 end do
!$omp end do nowait
             call self % LiftGradients(nEqn, nGradEqn, mesh, time, GetGradients, element_mask)
         end if
      end if
   
      end subroutine BR1_ComputeGradient
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for a strong computation of the gradient
!        ----------------------------------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_LiftGradients(self, nEqn, nGradEqn, mesh, time, GetGradients, element_mask, Level)
!
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t), intent(in) :: self
         integer,              intent(in) :: nEqn, nGradEqn
         class(HexMesh)                   :: mesh
         real(kind=RP),        intent(in) :: time
         procedure(GetGradientValues_f)   :: GetGradients
         logical, intent(in), optional    :: element_mask(:)
		 integer, intent(in), optional        :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i, j, k, m
         integer                :: eID , fID , dimID , eqID, fIDs(6), iFace, iEl, locLevel, locLevelm1
         logical                :: compute_element
         logical, allocatable   :: face_mask(:)

         if (present(element_mask)) then
            allocate(face_mask(size(mesh % faces)))
            !$omp parallel do schedule(runtime)
            do fID = 1, size(mesh % faces)
                associate(f => mesh % faces(fID))
                    face_mask(fID) = .false.
                    if (f % elementIDs(1) > 0) face_mask(fID) = element_mask(f % elementIDs(1))
                    if (f % elementIDs(2) > 0) face_mask(fID) = face_mask(fID) .or. element_mask(f % elementIDs(2))
                end associate
            end do
            !$omp end parallel do
        endif
		
		
		 if (present(Level)) then
            locLevel = Level
			locLevelm1 = max(locLevel-1,1)
!
!        *******************************************
!        Compute Riemann solvers of non-shared faces
!        *******************************************
!
!$omp do schedule(runtime) private(fID)
			 do iFace = 1, mesh % MLRK % MLIter(locLevelm1,3)
				fID = mesh % MLRK % MLIter_fID_Interior(iFace)
				compute_element = .true.
				if (present(element_mask)) compute_element = face_mask(fID)
				
				if (compute_element) then
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
                        call BR1_ComputeElementInterfaceAverage(self=self, masterFace1=mesh % faces(fID), nEqn=nEqn, nGradEqn=nGradEqn, GetGradients=GetGradients, &
                        f=mesh % faces(mesh % faces(fID)%Mortar(m)))
                     end if
                  end do
               elseif (mesh % faces(fID) % MortarType == MORTAR_NONE) then
                  call BR1_ComputeElementInterfaceAverage(self, mesh % faces(fID), nEqn, nGradEqn, GetGradients)
               end if
				endif
			 end do
!$omp end do nowait

!$omp do schedule(runtime) private(fID)
			 do iFace = 1, mesh % MLRK % MLIter(locLevelm1,4)
				fID = mesh % MLRK %  MLIter_fID_Boundary(iFace)
				compute_element = .true.
				if (present(element_mask)) compute_element = face_mask(fID)
				
				if (compute_element) then
				   call BR1_ComputeBoundaryFlux(self, mesh % faces(fID), nEqn, nGradEqn, time, GetGradients)
				endif
			 end do
!$omp end do 
!
!$omp do schedule(runtime) private(eID)
			 do iEl = 1, mesh % MLRK % MLIter(locLevel,9) 
				eID = mesh % MLRK % MLIter_eIDN_Seq(iEl)
				compute_element = .true.
				if (present(element_mask)) compute_element = element_mask(eID)
				
				if (compute_element) then
				   associate(e => mesh % elements(eID))
	!
	!              Add the surface integrals
	!              -------------------------
				   call BR1_GradientFaceLoop( self , nGradEqn, e, mesh)
	!
	!              Prolong gradients
	!              -----------------
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
				endif
			 end do
!$omp end do

#ifdef _HAS_MPI_
!$omp single
			 if ( MPI_Process % doMPIAction ) then 
				call mesh % GatherMPIFacesSolution(nEqn)
			 end if
!$omp end single

!$omp do schedule(runtime) private(fID)
			 do iFace = 1, mesh % MLRK % MLIter(locLevelm1,7)
				fID = mesh % MLRK % MLIter_fID_MPI(iFace)
				if (mesh % faces(fID) % MortarType == MORTAR_SLIDING) then
               associate(unStar=>mesh% faces(fID)%storage(1)%unStar)
                  unStar=0.0_RP
               end associate
               associate(unStar=>mesh% faces(fID)%storage(2)%unStar)
                  unStar=0.0_RP
               end associate
            end if
				if (mesh% faces(fID)%MortarType == MORTAR_BIG) then
               associate(unStar=>mesh% faces(fID)%storage(1)%unStar)
                  unStar=0.0_RP
               end associate
               do m=1,4
                  if (mesh % faces(fID)%Mortar(m) .ne. 0) then
                     call BR1_ComputeElementInterfaceAverage(self=self, masterFace1=mesh % faces(fID), nEqn=nEqn, nGradEqn=nGradEqn, GetGradients=GetGradients, &
                     f=mesh % faces(mesh % faces(fID)%Mortar(m)))
                  end if
               end do
            end if
				call BR1_ComputeMPIFaceAverage(self, mesh % faces(fID), nEqn, nGradEqn, GetGradients)
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
!$omp do schedule(runtime) private(eID) 
			 do iEl = 1, mesh % MLRK % MLIter(locLevel,10)
				eID = mesh % MLRK % MLIter_eIDN_MPI(iEl)
				associate(e => mesh % elements(eID))
	!
	!           Add the surface integrals
	!           -------------------------
				call BR1_GradientFaceLoop(self, nGradEqn, e, mesh)
	!
	!           Prolong gradients
	!           -----------------
				fIDs = e % faceIDs
				if ( .not.mesh % nonconforming ) then
               call e % ProlongGradientsToFaces(nGradEqn, mesh % faces(fIDs(1)),&
                                                mesh % faces(fIDs(2)),&
                                                mesh % faces(fIDs(3)),&
                                                mesh % faces(fIDs(4)),&
                                                mesh % faces(fIDs(5)),&
                                                mesh % faces(fIDs(6)))
            else
               call e % ProlongGradientsToFaces(nGradEqn,fFR=mesh % faces(fIDs(1)),&
                                                fBK=mesh % faces(fIDs(2)),&
                                                fBOT=mesh % faces(fIDs(3)),&
                                                fR=mesh % faces(fIDs(4)),&
                                                fT=mesh % faces(fIDs(5)),&
                                                fL=mesh % faces(fIDs(6)), faces=mesh%faces)
            end if
				end associate
			 end do
!$omp end do
#endif
         else
!
!        *******************************************
!        Compute Riemann solvers of non-shared faces
!        *******************************************
!
!$omp do schedule(runtime) private(fID)
			 do iFace = 1, size(mesh % faces_interior)
				fID = mesh % faces_interior(iFace)
				compute_element = .true.
				if (present(element_mask)) compute_element = face_mask(fID)
				
				if (compute_element) then
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
                     call BR1_ComputeElementInterfaceAverage(self=self, masterFace1=mesh % faces(fID), nEqn=nEqn, nGradEqn=nGradEqn, GetGradients=GetGradients, &
                     f=mesh % faces(mesh % faces(fID)%Mortar(m)))
                     end if 
                  end do 
               elseif (mesh % faces(fID) % MortarType == MORTAR_NONE) then
                  call BR1_ComputeElementInterfaceAverage(self, mesh % faces(fID), nEqn, nGradEqn, GetGradients)
               end if
				endif
			 end do
!$omp end do nowait

          if (mesh%SlidingMesh%active) then 
!$omp do schedule(runtime) private(fID)
                  do iFace = 1, size(mesh % mortar_faces)
                     fID = mesh % mortar_faces(iFace)%ID
                     associate(unStar=>mesh% mortar_faces(fID)%storage(1)%unStar)
                        unStar=0.0_RP
                     end associate
                     associate(unStar=>mesh% mortar_faces(fID)%storage(2)%unStar)
                        unStar=0.0_RP
                     end associate
                  end do 
!$omp end do        
               end if 

      if (mesh%SlidingMesh%active) then 
!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % mortar_faces)
            fID = mesh % mortar_faces(iFace)%ID
            call BR1_ComputeElementInterfaceAverage(self=self, masterFace1=mesh % faces (mesh % mortar_faces(fID)%Mortar(1)),masterFace2=mesh % faces (mesh % mortar_faces(fID)%Mortar(2)),&
             nEqn=nEqn, nGradEqn=nGradEqn, GetGradients=GetGradients, f=mesh % mortar_faces(fID), sliding=.true.)
         end do 
!$omp end do        
      end if


!$omp do schedule(runtime) private(fID)
			 do iFace = 1, size(mesh % faces_boundary)
				fID = mesh % faces_boundary(iFace)
				compute_element = .true.
				if (present(element_mask)) compute_element = face_mask(fID)
				
				if (compute_element) then
				   call BR1_ComputeBoundaryFlux(self, mesh % faces(fID), nEqn, nGradEqn, time, GetGradients)
				endif
			 end do
!$omp end do 
!
!$omp do schedule(runtime) private(eID)
			 do iEl = 1, size(mesh % elements_sequential)
				eID = mesh % elements_sequential(iEl)
				compute_element = .true.
				if (present(element_mask)) compute_element = element_mask(eID)
				
				if (compute_element) then
				   associate(e => mesh % elements(eID))
!
!              Add the surface integrals
!              -------------------------
				   call BR1_GradientFaceLoop( self , nGradEqn, e, mesh)
!
!              Prolong gradients
!              -----------------
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
				endif
			 end do
!$omp end do

#ifdef _HAS_MPI_
!$omp single
			 if ( MPI_Process % doMPIAction ) then 
				call mesh % GatherMPIFacesSolution(nEqn)
			 end if
!$omp end single

!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % faces_mpi)
            fID = mesh % faces_mpi(iFace)
            if (mesh% faces(fID)%MortarType == MORTAR_BIG) then 
               associate(unStar=>mesh% faces(fID)%storage(1)%unStar)
                  unStar=0.0_RP
               end associate
               do m=1,4
                  if (mesh % faces(fID)%Mortar(m) .ne. 0) then 
                     call BR1_ComputeElementInterfaceAverage(self=self, masterFace1=mesh % faces(fID), nEqn=nEqn, nGradEqn=nGradEqn, GetGradients=GetGradients, &
                     f=mesh % faces(mesh % faces(fID)%Mortar(m)))
                  end if 
               end do 
            end if 
            call BR1_ComputeMPIFaceAverage(self, mesh % faces(fID), nEqn, nGradEqn, GetGradients)
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
!$omp do schedule(runtime) private(eID) 
			 do iEl = 1, size(mesh % elements_mpi)
				eID = mesh % elements_mpi(iEl)
				associate(e => mesh % elements(eID))
	!
	!           Add the surface integrals
	!           -------------------------
				call BR1_GradientFaceLoop(self, nGradEqn, e, mesh)
	!
	!           Prolong gradients
	!           -----------------
				fIDs = e % faceIDs
				if ( .not.mesh % nonconforming ) then
               call e % ProlongGradientsToFaces(nGradEqn, mesh % faces(fIDs(1)),&
                                                mesh % faces(fIDs(2)),&
                                                mesh % faces(fIDs(3)),&
                                                mesh % faces(fIDs(4)),&
                                                mesh % faces(fIDs(5)),&
                                                mesh % faces(fIDs(6)))
            else
               call e % ProlongGradientsToFaces(nGradEqn,fFR=mesh % faces(fIDs(1)),&
                                                fBK=mesh % faces(fIDs(2)),&
                                                fBOT=mesh % faces(fIDs(3)),&
                                                fR=mesh % faces(fIDs(4)),&
                                                fT=mesh % faces(fIDs(5)),&
                                                fL=mesh % faces(fIDs(6)), faces=mesh%faces)
            end if
				end associate
			 end do
!$omp end do
#endif
         end if

      end subroutine BR1_LiftGradients

      subroutine BR1_LiftGradientsHO(self, nEqn, nGradEqn, mesh, time, GetGradients, element_mask, Level)
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t), intent(in) :: self
         integer,              intent(in) :: nEqn, nGradEqn
         class(HexMesh)                   :: mesh
         real(kind=RP),        intent(in) :: time
         procedure(GetGradientValues_f)   :: GetGradients
         logical, intent(in), optional    :: element_mask(:)
		 integer, intent(in), optional        :: Level
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i, j, k, m 
         integer                :: eID , fID , dimID , eqID, fIDs(6), iFace, iEl
         logical                :: compute_element
         logical, allocatable   :: face_mask(:)

         if (present(element_mask)) then
            allocate(face_mask(size(mesh % faces)))
            !$omp parallel do schedule(runtime)
            do fID = 1, size(mesh % faces)
                associate(f => mesh % faces(fID))
                    face_mask(fID) = .false.
                    if (f % elementIDs(1) > 0) face_mask(fID) = element_mask(f % elementIDs(1))
                    if (f % elementIDs(2) > 0) face_mask(fID) = face_mask(fID) .or. element_mask(f % elementIDs(2))
                end associate
            end do
            !$omp end parallel do
        endif
!
!        *******************************************
!        Compute Riemann solvers of non-shared faces
!        *******************************************
!
!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % HO_FacesInterior)
            fID = mesh % HO_FacesInterior(iFace)
            compute_element = .true.
            if (present(element_mask)) compute_element = face_mask(fID)
            
            if (compute_element) then
               call BR1_ComputeElementInterfaceAverage(self, mesh % faces(fID), nEqn, nGradEqn, GetGradients)
            endif
         end do
!$omp end do nowait

!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % HO_FacesBoundary)
            fID = mesh % HO_FacesBoundary(iFace)
            compute_element = .true.
            if (present(element_mask)) compute_element = face_mask(fID)
            
            if (compute_element) then
               call BR1_ComputeBoundaryFlux(self, mesh % faces(fID), nEqn, nGradEqn, time, GetGradients)
            endif
         end do
!$omp end do 
!
!$omp do schedule(runtime) private(eID)
         do iEl = 1, size(mesh % HO_ElementsSequential)
            eID = mesh % HO_ElementsSequential(iEl)
            compute_element = .true.
            if (present(element_mask)) compute_element = element_mask(eID)
            
            if (compute_element) then
               associate(e => mesh % elements(eID))
!  
!              Add the surface integrals
!              -------------------------
               call BR1_GradientFaceLoop( self , nGradEqn, e, mesh)
!  
!              Prolong gradients
!              -----------------
               fIDs = e % faceIDs
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

               end associate
            endif
         end do
!$omp end do

#ifdef _HAS_MPI_
!$omp single
         if ( MPI_Process % doMPIAction ) then 
            call mesh % GatherMPIFacesSolution(nEqn)
         end if
!$omp end single

!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % faces_mpi)
            fID = mesh % faces_mpi(iFace)
            call BR1_ComputeMPIFaceAverage(self, mesh % faces(fID), nEqn, nGradEqn, GetGradients)
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
!$omp do schedule(runtime) private(eID) 
         do iEl = 1, size(mesh % HO_ElementsMPI)
            eID = mesh % HO_ElementsMPI(iEl)
            associate(e => mesh % elements(eID))
!
!           Add the surface integrals
!           -------------------------
            call BR1_GradientFaceLoop(self, nGradEqn, e, mesh)
!
!           Prolong gradients
!           -----------------
            fIDs = e % faceIDs
            if ( .not.mesh % nonconforming ) then
            call e % ProlongGradientsToFaces(nGradEqn, mesh % faces(fIDs(1)),&
                                             mesh % faces(fIDs(2)),&
                                             mesh % faces(fIDs(3)),&
                                             mesh % faces(fIDs(4)),&
                                             mesh % faces(fIDs(5)),&
                                             mesh % faces(fIDs(6)))
            else
               call e % ProlongGradientsToFaces(nGradEqn,fFR=mesh % faces(fIDs(1)),&
               fBK=mesh % faces(fIDs(2)),&
               fBOT=mesh % faces(fIDs(3)),&
               fR=mesh % faces(fIDs(4)),&
               fT=mesh % faces(fIDs(5)),&
               fL=mesh % faces(fIDs(6)), faces=mesh%faces)
            end if 
 
            end associate
         end do
!$omp end do
#endif

      end subroutine BR1_LiftGradientsHO
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_GradientFaceLoop( self, nGradEqn, e, mesh )
         use ElementClass
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use DGIntegrals
         implicit none
         class(BassiRebay1_t),   intent(in)  :: self
         integer,                intent(in)  :: nGradEqn
         class(Element)                      :: e
         class(HexMesh)                      :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: faceInt_x(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP)        :: faceInt_y(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP)        :: faceInt_z(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         integer  :: i,j,k

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
            e % storage % U_x(:,i,j,k) = e % storage % U_x(:,i,j,k) + faceInt_x(:,i,j,k) * e % geom % InvJacobian(i,j,k)
            e % storage % U_y(:,i,j,k) = e % storage % U_y(:,i,j,k) + faceInt_y(:,i,j,k) * e % geom % InvJacobian(i,j,k)
            e % storage % U_z(:,i,j,k) = e % storage % U_z(:,i,j,k) + faceInt_z(:,i,j,k) * e % geom % InvJacobian(i,j,k)
         end do                  ; end do                   ; end do

      end subroutine BR1_GradientFaceLoop
!
      subroutine BR1_ComputeElementInterfaceAverage(self, f, nEqn, nGradEqn, GetGradients,masterFace1, masterFace2, sliding)
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         class(BassiRebay1_t),   intent(in)  :: self
         type(Face)                       :: f
         integer,    intent(in)           :: nEqn, nGradEqn
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
         real(kind=RP) :: uStar(nGradEqn)
         real(kind=RP) :: uStar_n(nGradEqn,NDIM,0:f % Nf(1), 0:f % Nf(2))
 

         integer       :: i,j, lm
         integer       :: Sidearray(2)
           do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
#ifdef MULTIPHASE
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL, rho_ = f % storage(1) % rho(i,j))
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR, rho_ = f % storage(2) % rho(i,j))
#else
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR)
#endif

#ifdef MULTIPHASE
            select case (self % eqName)
            case (ELLIPTIC_MU)
!
!              The multiphase solver needs the Chemical potential as first entropy variable
!              ----------------------------------------------------------------------------
               UL(IGMU) = f % storage(1) % mu(1,i,j)
               UR(IGMU) = f % storage(2) % mu(1,i,j)
            end select
#endif

   
            uStar = 0.5_RP * (UR - UL) * f % geom % jacobian(i,j)
            uStar_n(:,IX,i,j) = uStar * f % geom % normal(IX,i,j)
            uStar_n(:,IY,i,j) = uStar * f % geom % normal(IY,i,j)
            uStar_n(:,IZ,i,j) = uStar * f % geom % normal(IZ,i,j)
           end do               ; end do
         
         if (.not.present(sliding)) then 
            if (f % MortarType == MORTAR_NONE) then 
            Sidearray = (/1,2/)
            call f % ProjectGradientFluxToElements(nGradEqn, uStar_n,Sidearray,1)
            end if 
            if (f % MortarType == MORTAR_SMALL4 .and. present(masterFace1)) then 
               Sidearray = (/1,0/)
               call masterFace1 % ProjectMortarGradientFluxToElements(nEqn=nGradEqn, slaveFace=f, Hflux=uStar_n,whichElements=Sidearray,factor=1) 
               Sidearray = (/0,2/)
               call f % ProjectGradientFluxToElements(nGradEqn, uStar_n,Sidearray,1)
            end if 
         else 
            Sidearray = (/1,0/)
            call masterFace1 % ProjectMortarGradientFluxToElements(nEqn=nGradEqn, slaveFace=f, Hflux=uStar_n,whichElements=Sidearray,factor=1, sliding=sliding) 
            Sidearray = (/2,0/)
            call masterFace2 % ProjectMortarGradientFluxToElements(nEqn=nGradEqn, slaveFace=f, Hflux=uStar_n,whichElements=Sidearray,factor=1,sliding=sliding) 
         end if 
         
      end subroutine BR1_ComputeElementInterfaceAverage   

      subroutine BR1_ComputeMPIFaceAverage(self, f, nEqn, nGradEqn, GetGradients)
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         class(BassiRebay1_t),   intent(in)  :: self
         type(Face)                       :: f
         integer, intent(in)              :: nEqn, nGradEqn
         procedure(GetGradientValues_f)   :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: UL(nGradEqn), UR(nGradEqn)
         real(kind=RP) :: uStar(nGradEqn)
         real(kind=RP) :: uStar_n(nGradEqn,NDIM,0:f % Nf(1), 0:f % Nf(2))
         integer       :: i,j
         integer       :: Sidearray(2)
         
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
#ifdef MULTIPHASE
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL, rho_ = f % storage(1) % rho(i,j))
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR, rho_ = f % storage(2) % rho(i,j))
#else
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR)
#endif

#ifdef MULTIPHASE
            select case (self % eqName)
            case (ELLIPTIC_MU)
!
!              The multiphase solver needs the Chemical potential as first entropy variable
!              ----------------------------------------------------------------------------
               UL(IGMU) = f % storage(1) % mu(1,i,j)
               UR(IGMU) = f % storage(2) % mu(1,i,j)
            end select
#endif

   
            uStar = 0.5_RP * (UR - UL) * f % geom % jacobian(i,j)
            uStar_n(:,IX,i,j) = uStar * f % geom % normal(IX,i,j)
            uStar_n(:,IY,i,j) = uStar * f % geom % normal(IY,i,j)
            uStar_n(:,IZ,i,j) = uStar * f % geom % normal(IZ,i,j)
         end do               ; end do

         Sidearray = (/maxloc(f % elementIDs, dim = 1), HMESH_NONE/)
         call f % ProjectGradientFluxToElements(nGradEqn, uStar_n,Sidearray,1)

         if (f % MortarType == MORTAR_SMALL4) then 
            
            call f% Interpolatesmall2biggrad(nGradEqn, uStar_n)
            
         end if 
         
      end subroutine BR1_ComputeMPIFaceAverage   

      subroutine BR1_ComputeBoundaryFlux(self, f, nEqn, nGradEqn, time, GetGradients)
         use Physics
         use FaceClass
         implicit none
         class(BassiRebay1_t),   intent(in)  :: self
         type(Face)                       :: f
         integer, intent(in)              :: nEqn, nGradEqn
         real(kind=RP), intent(in)        :: time
         procedure(GetGradientValues_f)   :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i, j, bcID
         real(kind=RP) :: u_int(nGradEqn), u_star(nGradEqn)

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)

#ifdef MULTIPHASE
            call GetGradients(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), u_int, f % storage(1) % rho(i,j))
#else
            call GetGradients(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), u_int)
#endif

#ifdef MULTIPHASE
            select case (self % eqName)
            case (ELLIPTIC_MU)
!
!              The multiphase solver needs the Chemical potential as first entropy variable
!              ----------------------------------------------------------------------------
               u_int(IGMU) = f % storage(1) % mu(1,i,j)
            end select
#endif
   
            u_star = u_int
            call BCs(f % zone) % bc % GradVarsForEqn( nEqn,&
                                nGradEqn,                  &
                                f % geom % x(:,i,j),       &
                                time               ,       &
                                f % geom % normal(:,i,j),  &
                                f % storage(1) % Q(:,i,j), &
                                u_star, GetGradients           )

            f % storage(1) % unStar(:,1,i,j) = (u_star-u_int) * f % geom % normal(1,i,j) * f % geom % jacobian(i,j)
            f % storage(1) % unStar(:,2,i,j) = (u_star-u_int) * f % geom % normal(2,i,j) * f % geom % jacobian(i,j)    
            f % storage(1) % unStar(:,3,i,j) = (u_star-u_int) * f % geom % normal(3,i,j) * f % geom % jacobian(i,j)

         end do ; end do   

      end subroutine BR1_ComputeBoundaryFlux
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_ComputeInnerFluxes( self , nEqn, nGradEqn, EllipticFlux, GetViscosity, e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t) ,     intent (in) :: self
         integer,                   intent(in)  :: nEqn
         integer,                   intent(in)  :: nGradEqn
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
         real(kind=RP)       :: kappa(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: beta(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         integer             :: i, j, k

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

#elif defined(MULTIPHASE)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            call GetViscosity(e % storage % Q(IMC,i,j,k), mu(i,j,k))      
         end do                ; end do                ; end do
	
         mu = mu + e % storage % mu_NS(1,:,:,:) ! Add Subgrid Viscosity
         kappa = 0.0_RP
         beta  = multiphase % M0_star

#endif


         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            call EllipticFlux( nEqn, nGradEqn, e % storage % Q(:,i,j,k) , e % storage % U_x(:,i,j,k) , & 
                               e % storage % U_y(:,i,j,k) , e % storage % U_z(:,i,j,k), mu(i,j,k), beta(i,j,k), kappa(i,j,k), cartesianFlux)
            
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

      end subroutine BR1_ComputeInnerFluxes

      subroutine BR1_RiemannSolver ( self , nEqn, nGradEqn, EllipticFlux, f, QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
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
         class(BassiRebay1_t)            :: self
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

      end subroutine BR1_RiemannSolver
end module EllipticBR1