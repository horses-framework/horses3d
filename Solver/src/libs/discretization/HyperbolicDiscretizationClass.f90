#include "Includes.h"

#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE) || defined(ACOUSTIC)
module HyperbolicDiscretizationClass
   use SMConstants
#if defined(SPALARTALMARAS)
   use RiemannSolvers_NSSA
#elif defined(NAVIERSTOKES) && (!(TRANSPORT))
   use RiemannSolvers_NS
#elif defined(NAVIERSTOKES) && (TRANSPORT)
   use RiemannSolvers_NSTPT
#elif defined(INCNS)
   use RiemannSolvers_iNS
#elif defined(MULTIPHASE)
   use RiemannSolvers_MU
#elif defined(ACOUSTIC)
   use RiemannSolvers_CAA
#endif
   implicit none

   private
   public HyperbolicDiscretization_t

   type HyperbolicDiscretization_t
      contains
         procedure   :: Initialize               => BaseClass_Initialize
         procedure   :: ComputeInnerFluxes       => BaseClass_ComputeInnerFluxes
#if defined(NAVIERSTOKES) && !defined(SPALARTALMARAS)
         procedure   :: ComputeInnerFluxJacobian => BaseClass_ComputeInnerFluxJacobian
#endif
   end type HyperbolicDiscretization_t

   abstract interface
      pure subroutine HyperbolicFlux0D_f(Q, F, rho_, Qbase &
#ifdef TRANSPORT
      , transportVelocity &
#endif
)
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)           :: Q(1:NCONS     )
         real(kind=RP), intent(out)          :: F(1:NCONS, 1:NDIM)
         real(kind=RP), intent(in), optional :: rho_
#if defined(ACOUSTIC)
         real(kind=RP), intent(in), optional :: Qbase(1:NCONSB)
#else
         real(kind=RP), intent(in), optional :: Qbase(1:NCONS)
#endif
#ifdef TRANSPORT
         real(kind=RP), intent(in), optional :: transportVelocity(1:NDIM)
#endif
      end subroutine HyperbolicFlux0D_f
   end interface
!
!  ========
   contains
!  ========
!
      subroutine BaseClass_Initialize(self, controlVariables)
         use FTValueDictionaryClass
         use mainKeywordsModule
         use Headers
         use PhysicsStorage
         implicit none
         class(HyperbolicDiscretization_t) :: self
         class(FTValueDictionary),  intent(in)   :: controlVariables

!
!        Setup the Riemann solver
!        ------------------------
         call SetRiemannSolver(controlVariables)
!
!        Describe
!        --------
         if (.not. MPI_Process % isRoot ) return

         call Subsection_Header("Hyperbolic discretization")

         write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","Standard"

         call DescribeRiemannSolver

         if ( computeGradients ) then
            write(STD_OUT,'(30X,A,A30,A)') "->","Gradients computation: ", "Enabled."
         else
            write(STD_OUT,'(30X,A,A30,A)') "->","Gradients computation: ", "Disabled."
         end if

      end subroutine BaseClass_Initialize

      subroutine BaseClass_ComputeInnerFluxes( self , e , HyperbolicFlux, contravariantFlux,  useBaseFlow_)
         use ElementClass
         use Physics
         use PhysicsStorage
         implicit none
         class(HyperbolicDiscretization_t), intent(in)  :: self
         type(Element),           intent(in)  :: e
         procedure(HyperbolicFlux0D_f)        :: HyperbolicFlux
         real(kind=RP),           intent(out) :: contravariantFlux(1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
         logical, intent(in), optional        :: useBaseFlow_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: i, j, k
         real(kind=RP)      :: cartesianFlux(1:NCONS, 1:NDIM)
         logical            :: useBaseFlow
#ifdef TRANSPORT
         real(kind=RP)      :: transportVelocity(1:NDIM)
#endif

         if (present(useBaseFlow_)) then
            useBaseFlow = useBaseFlow_
         else
            useBaseFlow = .false.
         end if

         if (useBaseFlow) then

#if defined(ACOUSTIC)
             do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2)    ; do i = 0, e%Nxyz(1)
                call HyperbolicFlux( e % storage % Q(:,i,j,k), cartesianFlux(:,:), Qbase = e % storage % Qbase(:,i,j,k) )

                contravariantFlux(:,i,j,k,IX) =    cartesianFlux(:,IX) * e % geom % jGradXi(IX,i,j,k)  &
                                                 + cartesianFlux(:,IY) * e % geom % jGradXi(IY,i,j,k)  &
                                                 + cartesianFlux(:,IZ) * e % geom % jGradXi(IZ,i,j,k)


                contravariantFlux(:,i,j,k,IY) =   cartesianFlux(:,IX) * e % geom % jGradEta(IX,i,j,k)  &
                                                 + cartesianFlux(:,IY) * e % geom % jGradEta(IY,i,j,k)  &
                                                 + cartesianFlux(:,IZ) * e % geom % jGradEta(IZ,i,j,k)


                contravariantFlux(:,i,j,k,IZ) =   cartesianFlux(:,IX) * e % geom % jGradZeta(IX,i,j,k)  &
                                                 + cartesianFlux(:,IY) * e % geom % jGradZeta(IY,i,j,k)  &
                                                 + cartesianFlux(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)
             end do               ; end do                ; end do

#endif
         else

             do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2)    ; do i = 0, e%Nxyz(1)
#ifdef TRANSPORT
                call e % getTransportVelocity(i, j, k, transportVelocity)
                call HyperbolicFlux( e % storage % Q(:,i,j,k), cartesianFlux(:,:), rho_ = e % storage % rho(i,j,k), transportVelocity = transportVelocity)
#else
                call HyperbolicFlux( e % storage % Q(:,i,j,k), cartesianFlux(:,:), rho_ = e % storage % rho(i,j,k))
#endif

                contravariantFlux(:,i,j,k,IX) =    cartesianFlux(:,IX) * e % geom % jGradXi(IX,i,j,k)  &
                                                 + cartesianFlux(:,IY) * e % geom % jGradXi(IY,i,j,k)  &
                                                 + cartesianFlux(:,IZ) * e % geom % jGradXi(IZ,i,j,k)


                contravariantFlux(:,i,j,k,IY) =   cartesianFlux(:,IX) * e % geom % jGradEta(IX,i,j,k)  &
                                                 + cartesianFlux(:,IY) * e % geom % jGradEta(IY,i,j,k)  &
                                                 + cartesianFlux(:,IZ) * e % geom % jGradEta(IZ,i,j,k)


                contravariantFlux(:,i,j,k,IZ) =   cartesianFlux(:,IX) * e % geom % jGradZeta(IX,i,j,k)  &
                                                 + cartesianFlux(:,IY) * e % geom % jGradZeta(IY,i,j,k)  &
                                                 + cartesianFlux(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)

             end do               ; end do                ; end do
         end if

      end subroutine BaseClass_ComputeInnerFluxes
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------
!     Subroutine to get the Jacobian of the contravariant fluxes
!     -> dFdQ (:,:,i,j,k,dim)
!                 |     |
!              jac|coord|flux in cartesian direction dim
!     ----------------------------------------------------------
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
      subroutine BaseClass_ComputeInnerFluxJacobian( self, e, dFdQ)
         use ElementClass
         use Physics
         use PhysicsStorage
         implicit none
         !--------------------------------------------
         class(HyperbolicDiscretization_t), intent(in)  :: self
         type(Element)          , intent(in)  :: e
         real(kind=RP)          , intent(out) :: dFdQ( NCONS, NCONS, NDIM , 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         !--------------------------------------------
         real(kind=RP), DIMENSION(NCONS,NCONS)  :: dfdq_,dgdq_,dhdq_
         integer                                :: i,j,k
         !--------------------------------------------

         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)

            call InviscidJacobian(e % storage % Q(:,i,j,k),dfdq_,dgdq_,dhdq_)


            dFdQ(:,:,IX,i,j,k) = e % geom % jGradXi  (1,i,j,k) * dfdq_ + &
                                 e % geom % jGradXi  (2,i,j,k) * dgdq_ + &
                                 e % geom % jGradXi  (3,i,j,k) * dhdq_

            dFdQ(:,:,IY,i,j,k) = e % geom % jGradEta (1,i,j,k) * dfdq_ + &
                                 e % geom % jGradEta (2,i,j,k) * dgdq_ + &
                                 e % geom % jGradEta (3,i,j,k) * dhdq_

            dFdQ(:,:,IZ,i,j,k) = e % geom % jGradZeta(1,i,j,k) * dfdq_ + &
                                 e % geom % jGradZeta(2,i,j,k) * dgdq_ + &
                                 e % geom % jGradZeta(3,i,j,k) * dhdq_
         end do                ; end do                ; end do

      end subroutine BaseClass_ComputeInnerFluxJacobian
#endif
end module HyperbolicDiscretizationClass
#endif
