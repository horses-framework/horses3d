#include "Includes.h"
module EllipticDiscretizationClass
   use SMConstants
   use MeshTypes
   use Physics
   use PhysicsStorage
   use MPI_Face_Class
   use VariableConversion
   implicit none
!
   private
   public   EllipticDiscretization_t, EllipticFlux_f, GetViscosity_f
   public   ELLIPTIC_NS, ELLIPTIC_NSSA, ELLIPTIC_iNS, ELLIPTIC_CH, ELLIPTIC_MU

   public BaseClass_ComputeGradient

   type EllipticDiscretization_t
      integer                                        :: eqName
      real(kind=RP)                                  :: sigma = 1.0_RP
      contains
         procedure      :: Construct                 => BaseClass_Construct
         procedure      :: ComputeGradient           => BaseClass_ComputeGradient
         procedure      :: ComputeLocalGradients     => BaseClass_ComputeGradient
         procedure      :: LiftGradients             => BaseClass_LiftGradients
         procedure      :: ComputeInnerFluxes        => BaseClass_ComputeInnerFluxes
         procedure      :: RiemannSolver             => BaseClass_RiemannSolver
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
         procedure      :: ComputeInnerFluxJacobian  => BaseClass_ComputeInnerFluxJacobian
         procedure      :: RiemannSolver_Jacobians   => BaseClass_RiemannSolver_Jacobians
#endif
         procedure      :: Describe                  => BaseClass_Describe
   end type EllipticDiscretization_t

   abstract interface
      pure subroutine EllipticFlux_f(nEqn, nGradEqn, Q, U_x, U_y, U_z, mu, beta, kappa, &
#ifdef TRANSPORT
transportD, &
#endif
F)
         use SMConstants
         use PhysicsStorage
         implicit none
         integer,       intent(in)  :: nEqn, nGradEqn
         real(kind=RP), intent(in)  :: Q   (1:nEqn     )
         real(kind=RP), intent(in)  :: U_x (1:nGradEqn)
         real(kind=RP), intent(in)  :: U_y (1:nGradEqn)
         real(kind=RP), intent(in)  :: U_z (1:nGradEqn)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: beta
         real(kind=RP), intent(in)  :: kappa
#ifdef TRANSPORT
         real(kind=RP), intent(in)  :: transportD(1:NDIM, 1:NDIM)
#endif
         real(kind=RP), intent(out) :: F(1:nEqn, 1:NDIM)
      end subroutine EllipticFlux_f
      pure subroutine GetViscosity_f(phi, mu)
         use SMConstants
         implicit none
         real(kind=RP), intent(in)  :: phi
         real(kind=RP), intent(out) :: mu
      end subroutine GetViscosity_f
   end interface

   enum, bind(C)
      enumerator :: ELLIPTIC_NS, ELLIPTIC_NSSA, ELLIPTIC_CH, ELLIPTIC_MU, ELLIPTIC_iNS
   end enum
!
!  ========
   contains
!  ========
!
      subroutine BaseClass_Construct(self, controlVariables, eqName)
         use FTValueDictionaryClass
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         use PhysicsStorage
         implicit none
         class(EllipticDiscretization_t)       :: self
         class(FTValueDictionary), intent(in)  :: controlVariables
         integer, intent(in)                   :: eqName
!
!        Request the penalty parameter
!        -----------------------------
         if ( controlVariables % containsKey("penalty parameter") ) then
            self % sigma = controlVariables % doublePrecisionValueForKey("penalty parameter")

         else
!            
!           Set 0.0 by default
!           ------------------
            self % sigma = 0.0_RP

         end if


         select case (eqName)
         case(ELLIPTIC_NS)
            self % eqName = ELLIPTIC_NS

         case(ELLIPTIC_NSSA)
            self % eqName = ELLIPTIC_NSSA

         case(ELLIPTIC_CH)
            self % eqName = ELLIPTIC_CH

         case(ELLIPTIC_MU)
            self % eqName = ELLIPTIC_MU
         
         case(ELLIPTIC_iNS)
            self % eqName = ELLIPTIC_iNS

         case default
            print*, "Unrecognized equation"
            errorMessage(STD_OUT)
            error stop

         end select
      

      end subroutine BaseClass_Construct

      subroutine BaseClass_Describe(self)
         implicit none
         class(EllipticDiscretization_t), intent(in)  :: self

      end subroutine BaseClass_Describe

      subroutine BaseClass_ComputeGradient(self, nEqn, nGradEqn, mesh, time, GetGradients, HO_Elements, element_mask, Level)
!
!        *****************************************************
!           BaseClass computes Local Gradients by default
!              Do not change.. Used by ComputeTimeDerivativeIsolated
!        *****************************************************
!           
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         class(EllipticDiscretization_t), intent(in) :: self
         integer,             intent(in)  :: nEqn, nGradEqn
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
         integer  :: eID, lID, locLevel
         logical  :: set_mu
         logical  :: compute_element

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

		 if (present(Level)) then
             locLevel = Level
!$omp do schedule(runtime) private(eID)
			 do lID = 1 , mesh % MLRK % MLIter(locLevel,8) 
				eID = mesh % MLRK % MLIter_eIDN(lID)
				compute_element = .true.
				if (present(element_mask)) compute_element = element_mask(eID)
				
				if (compute_element) then
				   call mesh % elements(eID) % ComputeLocalGradient(nEqn, nGradEqn, GetGradients, set_mu)
				endif
			 end do
!$omp end do nowait
         else
!$omp do schedule(runtime)
         do eID = 1 , size(mesh % elements)						  
            compute_element = .true.
            if (present(element_mask)) compute_element = element_mask(eID)
            
            if (compute_element) then
               call mesh % elements(eID) % ComputeLocalGradient(nEqn, nGradEqn, GetGradients, set_mu)
            endif
         end do
!$omp end do nowait
         end if

      end subroutine BaseClass_ComputeGradient

      subroutine BaseClass_LiftGradients(self, nEqn, nGradEqn, mesh, time, GetGradients, element_mask, Level)
!
!        *****************************************************
!        Lift gradients: do nothing here
!        *****************************************************
!           
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         class(EllipticDiscretization_t), intent(in) :: self
         integer,             intent(in)  :: nEqn, nGradEqn
         class(HexMesh)                   :: mesh
         real(kind=RP),        intent(in) :: time
         procedure(GetGradientValues_f)   :: GetGradients
         logical, intent(in), optional    :: element_mask(:)
		 integer, intent(in), optional    :: Level
		 
		 integer  :: locLevel
		 
		 if (present(Level)) then
            locLevel = Level
         else
            locLevel = 1
         end if

      end subroutine BaseClass_LiftGradients

      subroutine BaseClass_ComputeInnerFluxes( self, nEqn, nGradEqn, EllipticFlux, GetViscosity, e, contravariantFlux )
         use ElementClass
         use PhysicsStorage
         implicit none
         class(EllipticDiscretization_t) ,  intent(in) :: self
         integer,                           intent(in) :: nEqn, nGradEqn
         procedure(EllipticFlux_f)                     :: EllipticFlux
         procedure(GetViscosity_f)                     :: GetViscosity
         type(Element)                                 :: e
         real(kind=RP)           ,  intent (out)       :: contravariantFlux(1:nEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         contravariantFlux = 0.0_RP

      end subroutine BaseClass_ComputeInnerFluxes
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------------------------------------------------
!     BaseClass_ComputeInnerFluxJacobian:
!
!     Subroutine to get the Jacobians of the contravariant viscous fluxes on the inner Gauss points of an element
!     ( \tilde{G}_{xx} )
!
!     -> Jacobian with respect to ∇q: a 4th order tensor of the form
!        df_dgradq (:,:,:,:,i,j,k)
!                   |_| | | |_|_| 
!                    |  | |  |
!                    |  | |  |____Coordinate indexes in element 
!                    |  | |_______dim: flux direction: f(1), g(2), h(3) [first  index of \tilde{G} matrix]
!                    |  |_________∇q component: 1, 2, 3                 [second index of \tilde{G} matrix]
!                    |____________Jacobian for this component
!        
!
!     -> Jacobian with respect to q: a 3rd order tensor of the form
!             dFdQ (:,:,i,j,k,dim)
!                   |_| |_|_| |
!                    |   |    |_flux direction: f(1), g(2), h(3)
!                    |   |______Coordinate indexes in element 
!                    |__________Jacobian for this component
!              (added to existing Jacobian)
!     --------------------------------------------------------------------------------------------------
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
      subroutine BaseClass_ComputeInnerFluxJacobian( self, e, df_dgradq, dFdQ) 
         use ElementClass
         use Physics
         use PhysicsStorage
         implicit none
         !--------------------------------------------
         class(EllipticDiscretization_t), intent(in)    :: self
         type(Element)         , intent(inout) :: e                                                                                          !<  This element
         real(kind=RP)         , intent(out)   :: df_dgradq( NCONS, NCONS, NDIM, NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3) )         !>  Contravariant Jacobian with respect to ∇q
         real(kind=RP)         , intent(inout) :: dFdQ     ( NCONS, NCONS, NDIM      , 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3) )   !<> Contravariant Jacobian with respect to q
         !--------------------------------------------
         real(kind=RP), DIMENSION(NCONS,NCONS,NDIM,NDIM) :: df_dgradq_cart ! Cartesian Jacobian tensor with respect to ∇q
         real(kind=RP), DIMENSION(NCONS,NCONS,NDIM)      :: dfdq_cart      ! Cartesian Jacobian tensor with respect to q
#ifdef TRANSPORT
         real(kind=RP), DIMENSION(1:NDIM,1:NDIM)         :: transportD
#endif
         integer :: i,j,k     ! Coordinate counters
         integer :: i1, i2    ! Index of G_xx
         !--------------------------------------------
         
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
#ifdef TRANSPORT
            call e % getTransportDiffusion(i, j, k, transportD)
            call ViscousJacobian ( Q   = e % storage % Q  (:,i,j,k), &
                                   Q_x = e % storage % U_x(:,i,j,k), &
                                   Q_y = e % storage % U_y(:,i,j,k), &
                                   Q_z = e % storage % U_z(:,i,j,k), &
                                   transportD = transportD, &
                                   df_dgradq = df_dgradq_cart, &
                                   df_dq     = dfdq_cart)
#else
            call ViscousJacobian ( Q   = e % storage % Q  (:,i,j,k), &
                                   Q_x = e % storage % U_x(:,i,j,k), &
                                   Q_y = e % storage % U_y(:,i,j,k), &
                                   Q_z = e % storage % U_z(:,i,j,k), &
                                   df_dgradq = df_dgradq_cart, &
                                   df_dq     = dfdq_cart)
#endif
!
!           ***********************************
!           Jacobian with respect to ∇q: dF/d∇q
!           ***********************************
!
!           Fill contravariant Jacobian tensor
!           ----------------------------------
            df_dgradq(:,:,:,:,i,j,k) = 0._RP
!
!           \tilde{G}_{1x}
!           --------------
            do i2 = 1, NDIM ; do i1 = 1, NDIM
               df_dgradq(:,:,i1,1,i,j,k) = df_dgradq(:,:,i1,1,i,j,k) + df_dgradq_cart(:,:,i1,i2) * e % geom % jGradXi  (i2,i,j,k)
            end do          ; end do
!
!           \tilde{G}_{2x}
!           --------------
            do i2 = 1, NDIM ; do i1 = 1, NDIM
               df_dgradq(:,:,i1,2,i,j,k) = df_dgradq(:,:,i1,2,i,j,k) + df_dgradq_cart(:,:,i1,i2) * e % geom % jGradEta (i2,i,j,k)
            end do          ; end do
!
!           \tilde{G}_{3x}
!           --------------
            do i2 = 1, NDIM ; do i1 = 1, NDIM
               df_dgradq(:,:,i1,3,i,j,k) = df_dgradq(:,:,i1,3,i,j,k) + df_dgradq_cart(:,:,i1,i2) * e % geom % jGradZeta(i2,i,j,k)
            end do          ; end do
            
!
!           *********************************
!           Jacobian with respect to q: dF/dq
!           *********************************
!
            
            dFdQ(:,:,IX,i,j,k) = dFdQ(:,:,IX,i,j,k) - (  e % geom % jGradXi  (1,i,j,k) * dfdq_cart(:,:,1) + &
                                                         e % geom % jGradXi  (2,i,j,k) * dfdq_cart(:,:,2) + &
                                                         e % geom % jGradXi  (3,i,j,k) * dfdq_cart(:,:,3) )

            dFdQ(:,:,IY,i,j,k) = dFdQ(:,:,IY,i,j,k) - (  e % geom % jGradEta (1,i,j,k) * dfdq_cart(:,:,1) + &
                                                         e % geom % jGradEta (2,i,j,k) * dfdq_cart(:,:,2) + &
                                                         e % geom % jGradEta (3,i,j,k) * dfdq_cart(:,:,3) )

            dFdQ(:,:,IZ,i,j,k) = dFdQ(:,:,IZ,i,j,k) - (  e % geom % jGradZeta(1,i,j,k) * dfdq_cart(:,:,1) + &
                                                         e % geom % jGradZeta(2,i,j,k) * dfdq_cart(:,:,2) + &
                                                         e % geom % jGradZeta(3,i,j,k) * dfdq_cart(:,:,3) )
            
         end do                ; end do                ; end do
      end subroutine BaseClass_ComputeInnerFluxJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BaseClass_RiemannSolver_Jacobians( self, f) 
         use FaceClass
         use Physics
         use PhysicsStorage
         implicit none
         !--------------------------------------------
         class(EllipticDiscretization_t), intent(in)    :: self
         type(Face)            , intent(inout) :: f
         !--------------------------------------------
!$omp single
         error stop ':: RiemannSolver_Jacobians not implemented for selected Viscous Method'
!$omp end single
      end subroutine BaseClass_RiemannSolver_Jacobians
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#endif
      subroutine BaseClass_RiemannSolver ( self, nEqn, nGradEqn, EllipticFlux, f, QLeft, QRight, U_xLeft, U_yLeft, U_zLeft, U_xRight, U_yRight, U_zRight, &
                                           mu_left, mu_right, nHat , dWall, &
#ifdef MULTIPHASE
sigma, & 
#endif
#ifdef TRANSPORT
transportDLeft, transportDRight, &
#endif
flux )
         use SMConstants
         use SMConstants
         use PhysicsStorage
         use FaceClass
         implicit none
         class(EllipticDiscretization_t) :: self
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
#ifdef TRANSPORT
         real(kind=RP), intent(in)     :: transportDLeft(1:NDIM,1:NDIM)
         real(kind=RP), intent(in)     :: transportDRight(1:NDIM,1:NDIM)
#endif
         real(kind=RP), intent(out)      :: flux(nEqn)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         flux = 0.0_RP

      end subroutine BaseClass_RiemannSolver
end module EllipticDiscretizationClass