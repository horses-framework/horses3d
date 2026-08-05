module VolumeIntegrals
   use SMConstants
   use PhysicsStorage
   use Physics
   use VariableConversion
   use ElementClass
   use HexMeshClass
   use FluidData
   use NodalStorageClass, only: NodalStorage, NodalStorage_t
#ifdef _HAS_MPI_
   use mpi
   use MPI_Utilities, only: MPI_MinMax
#endif
#include "Includes.h"

   private
   public   VOLUME

#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
   public KINETIC_ENERGY, KINETIC_ENERGY_RATE, KINETIC_ENERGY_BALANCE, ENSTROPHY, VELOCITY
   public ENTROPY, ENTROPY_RATE, INTERNAL_ENERGY, MOMENTUM, SOURCE, PSOURCE, ARTIFICIAL_DISSIPATION
   public RELATIVE_KINETIC_ENERGY, RELATIVE_TOTAL_ENERGY
   public ENTROPY_BALANCE, MATH_ENTROPY
#endif

#if defined(INCNS)
   public MASS, ENTROPY, KINETIC_ENERGY_RATE, ENTROPY_RATE, SOURCE
#endif

#if defined(MULTIPHASE)
   public ENTROPY_RATE, ENTROPY_BALANCE, PHASE2_AREA, PHASE2_XCOG, PHASE2_XVEL, SOURCE
#endif

#if defined(CAHNHILLIARD)
   public FREE_ENERGY
#endif

#if defined(ACOUSTIC)
   public ACOUSTIC_ENERGY, SOURCE
#endif

   public   ScalarVolumeIntegral, VectorVolumeIntegral, GetSensorRange



   enum, bind(C)
      enumerator :: VOLUME
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
      enumerator :: KINETIC_ENERGY, KINETIC_ENERGY_RATE, KINETIC_ENERGY_BALANCE
      enumerator :: ENSTROPHY, VELOCITY, ENTROPY, ENTROPY_RATE, INTERNAL_ENERGY, MOMENTUM, SOURCE, PSOURCE
      enumerator :: RELATIVE_KINETIC_ENERGY, RELATIVE_TOTAL_ENERGY ! Relative to the IB velocity
      enumerator :: ARTIFICIAL_DISSIPATION, ENTROPY_BALANCE, MATH_ENTROPY
#endif
#if defined(INCNS)
      enumerator :: MASS, ENTROPY, KINETIC_ENERGY_RATE, ENTROPY_RATE, SOURCE
#endif
#if defined(MULTIPHASE)
      enumerator :: ENTROPY_RATE, ENTROPY_BALANCE, PHASE2_AREA, PHASE2_XCOG, PHASE2_XVEL, SOURCE
#endif
#if defined(CAHNHILLIARD)
      enumerator :: FREE_ENERGY
#endif
#if defined(ACOUSTIC)
      enumerator :: ACOUSTIC_ENERGY, SOURCE
#endif
   end enum
!
!  ========
   contains
!  ========
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!           SCALAR INTEGRALS PROCEDURES
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      function ScalarVolumeIntegral(mesh, integralType, time) result(val)
!
!        -----------------------------------------------------------
!           This function computes scalar integrals, that is, those
!           in the form:
!                 val = \int v dx
!        -----------------------------------------------------------
!

         implicit none
         class(HexMesh),      intent(in)  :: mesh
         integer,             intent(in)  :: integralType
         real(kind=RP), optional, intent(in) :: time
         real(kind=RP)                    :: val
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: localVal
         integer       :: eID, ierr

         if ((integralType .eq. RELATIVE_KINETIC_ENERGY) .or. (integralType .eq. RELATIVE_TOTAL_ENERGY)) then
            val = ScalarVolumeIntegral_RELATIVE(mesh, integralType, time)
            return
         end if

!
!        Initialization
!        --------------
         val = 0.0_RP
!
!        Loop the mesh
!        -------------
!$omp parallel do reduction(+:val) private(eID) schedule(guided)
         do eID = 1, mesh % no_of_elements
!
!           Compute the integral
!           --------------------
            val = val + ScalarVolumeIntegral_Local(mesh % elements(eID), &
                                                           integralType    )

         end do
!$omp end parallel do

#ifdef _HAS_MPI_
            localVal = val
            call mpi_allreduce(localVal, val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      end function ScalarVolumeIntegral
   
      function ScalarVolumeIntegral_RELATIVE(mesh, integralType, time) result(val)
!
!        -----------------------------------------------------------
!           This function computes scalar integrals, that is, those
!           in the form:
!                 val = \int v dx
!        -----------------------------------------------------------
!

         implicit none
         class(HexMesh),      intent(in)  :: mesh
         integer,             intent(in)  :: integralType
         real(kind=RP),       intent(in)  :: time
         real(kind=RP)                    :: val
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: localVal
         integer       :: eID, ierr
!
!        Initialization
!        --------------
         val = 0.0_RP
!
!        Loop the mesh
!        -------------
!$omp parallel do reduction(+:val) private(eID) schedule(guided)
         do eID = 1, mesh % no_of_elements
!
!           Compute the integral
!           --------------------
            val = val + ScalarVolumeIntegral_Local(mesh % elements(eID), &
                                                           integralType, &
                                                           time = time, &
                                                           IBM = mesh % IBM    )

         end do
!$omp end parallel do

#ifdef _HAS_MPI_
            localVal = val
            call mpi_allreduce(localVal, val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      end function ScalarVolumeIntegral_RELATIVE

      function ScalarVolumeIntegral_Local(e, integralType, time, IBM) result(val)
         use IBMClass, only: IBM_type
         implicit none
         class(Element),      target, intent(in)     :: e
         integer,                     intent(in)     :: integralType
         real(kind=RP), optional, intent(in) :: time
         class(IBM_type), optional, intent(in) :: IBM
         real(kind=RP)                               :: val
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: Nel(3)    ! Element polynomial order
         integer     :: i, j, k
         real(kind=RP)           :: EntropyVars(NCONS)
         real(kind=RP)           :: ViscFlux(NCONS,NDIM)
         real(kind=RP)           :: KinEn(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: U_x(NDIM,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: U_y(NDIM,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: U_z(NDIM,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: uvw(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: inv_rho(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: p_3d(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: grad_Mp(1:NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: M_grad_p(1:NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: correction_term(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: Ma2(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: p, s, ms
         real(kind=RP), pointer  :: Qb(:)
         real(kind=RP)           :: free_en, fchem, entr, area, rho , u , v, w, en, thetaeddy
         real(kind=RP)           :: Strain(NDIM,NDIM)
         real(kind=RP)           :: mu
         real(kind=RP)           :: Q_target(NCONS) ! Used for the IB velocity

         Nel = e % Nxyz

         associate ( wx => NodalStorage(e % Nxyz(1)) % w, &
                     wy => NodalStorage(e % Nxyz(2)) % w, &
                     wz => NodalStorage(e % Nxyz(3)) % w    )
!
!        Initialization
!        --------------
         val = 0.0_RP
!
!        Perform the numerical integration
!        ---------------------------------
         select case ( integralType )

         case ( VOLUME )
!
!           **********************************
!           Computes the volume integral
!              val = \int dV
!           **********************************
!
            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k)
            end do            ; end do           ; end do

#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
         case ( KINETIC_ENERGY )
!
!           ***********************************
!              Computes the kinetic energy
!              integral:
!              K = \int \rho V^2 dV
!           ***********************************
!

            KinEn =         POW2(e % storage % Q(IRHOU,:,:,:))
            KinEn = KinEn + POW2( e % storage % Q(IRHOV,:,:,:) )
            KinEn = KinEn + POW2( e % storage % Q(IRHOW,:,:,:) )
            KinEn = 0.5_RP * KinEn / e % storage % Q(IRHO,:,:,:)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * kinEn(i,j,k)
            end do            ; end do           ; end do

         case ( KINETIC_ENERGY_RATE )
!
!           ***********************************
!              Computes the kinetic energy
!              time derivative:
!              K_t = (d/dt)\int \rho V^2 dV
!           ***********************************
!
            uvw = e % storage % Q(IRHOU,:,:,:) / e % storage % Q(IRHO,:,:,:)
            KinEn = uvw * e % storage % QDot(IRHOU,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(IRHO,:,:,:)

            uvw = e % storage % Q(IRHOV,:,:,:) / e % storage % Q(IRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(IRHOV,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(IRHO,:,:,:)

            uvw = e % storage % Q(IRHOW,:,:,:) / e % storage % Q(IRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(IRHOW,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(IRHO,:,:,:)



            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * kinEn(i,j,k)
            end do            ; end do           ; end do

         case ( KINETIC_ENERGY_BALANCE )
!
!           ***************************************************
!              Computes the kinetic energy
!              balance: will also work the for Pirozzoli scheme
!              with energy gradient variables.
!           ***************************************************
!
            inv_rho = 1.0_RP / e % storage % Q(IRHO,:,:,:)
            uvw = e % storage % Q(IRHOU,:,:,:) * inv_rho
            KinEn = uvw * e % storage % QDot(IRHOU,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(IRHO,:,:,:)

            uvw = e % storage % Q(IRHOV,:,:,:) * inv_rho
            KinEn = KinEn + uvw * e % storage % QDot(IRHOV,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(IRHO,:,:,:)

            uvw = e % storage % Q(IRHOW,:,:,:) * inv_rho
            KinEn = KinEn + uvw * e % storage % QDot(IRHOW,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(IRHO,:,:,:)
!
!           I also need the pressure work
!           -----------------------------
            p_3d = thermodynamics % gammaMinus1*(e % storage % Q(IRHOE,:,:,:) &
                     - 0.5_RP*(POW2(e % storage % Q(IRHOU,:,:,:)) + POW2(e % storage % Q(IRHOV,:,:,:)) + &
                               POW2(e % storage % Q(IRHOW,:,:,:)))*inv_rho)
!
!           This correction term is because the split-form scheme will compute the pressure terms with de-aliased metrics, so I need
!           to replicate that de-aliasing: 0.5<u,M∇p> + 0.5<u,∇(Mp)> = <u,∇(Mp)> + 0.5(<u,M∇p-∇(Mp)>).
!           The first term is computed as <-p,M∇·u>, using the gradients, that already contain the surface term.
!           ------------------------------------------------------------------------------------------------------------------------
            call GetPressureLocalGradient(Nel, p_3d, e % geom % jGradXi, e % geom % jGradEta, e % geom % jGradZeta, grad_Mp, M_grad_p)

            correction_term = 0.5_RP * (  e % storage % Q(IRHOU,:,:,:)*(M_grad_p(IX,:,:,:)-grad_Mp(IX,:,:,:)) &
                                        + e % storage % Q(IRHOV,:,:,:)*(M_grad_p(IY,:,:,:)-grad_Mp(IY,:,:,:)) &
                                        + e % storage % Q(IRHOW,:,:,:)*(M_grad_p(IZ,:,:,:)-grad_Mp(IZ,:,:,:)))*inv_rho

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               call ViscousFlux_ENERGY(NCONS, NGRAD, e % storage % Q(:,i,j,k), e % storage % U_x(:,i,j,k), e % storage % U_y(:,i,j,k), &
                                    e % storage % U_z(:,i,j,k), e % storage % mu_ns(1,i,j,k), 0.0_RP, e % storage % mu_ns(2,i,j,k), ViscFlux)

               val = val + wx(i) * wy(j) * wz(k) * (e % geom % jacobian(i,j,k) * ( kinEn(i,j,k) + &
                     sum(ViscFlux(2:4,1)*e % storage % U_x(2:4,i,j,k)+ViscFlux(2:4,2)*e % storage % U_y(2:4,i,j,k)+ViscFlux(2:4,3)*e % storage % U_z(2:4,i,j,k)) &
                   - p_3d(i,j,k)*(e % storage % U_x(IRHOU,i,j,k) + e % storage % U_y(IRHOV,i,j,k) + e % storage % U_z(IRHOW,i,j,k))) &
                   + correction_term(i,j,k))
            end do            ; end do           ; end do

            val = val + e % storage % artificialDiss


         case ( ENSTROPHY )
!
!           ***************************
!           Computes the flow enstrophy
!           ***************************
!

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               call getVelocityGradients(e % storage % Q(:,i,j,k),&
                                         e % storage % U_x(:,i,j,k), &
                                         e % storage % U_y(:,i,j,k), &
                                         e % storage % U_z(:,i,j,k), &
                                         U_x(:,i,j,k), U_y(:,i,j,k), U_z(:,i,j,k))

               KinEn =   POW2( U_y(IZ,i,j,k) - U_z(IY,i,j,k) ) &
                       + POW2( U_z(IX,i,j,k) - U_x(IZ,i,j,k) ) &
                       + POW2( U_x(IY,i,j,k) - U_y(IX,i,j,k) )

               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * kinEn(i,j,k)
            end do            ; end do           ; end do

         case ( VELOCITY )

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * sqrt(   POW2(e % storage % Q(IRHOU,i,j,k)) &
                                                           + POW2(e % storage % Q(IRHOV,i,j,k)) &
                                                           + POW2(e % storage % Q(IRHOW,i,j,k))  ) &
                                        / e % storage % Q(IRHO,i,j,k) * e % geom % jacobian(i,j,k)
            end do            ; end do           ; end do

         case ( ENTROPY )
!
!           ********************************************
!              Computes the specific entropy integral
!           ********************************************
!
            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               p = Pressure( e % storage % Q(:,i,j,k) )
               s = ( log(p) - thermodynamics % gamma * log(e % storage % Q(IRHO,i,j,k)) )
               val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * s
            end do            ; end do           ; end do

          case ( MATH_ENTROPY )
!
!           ******************************************************************
!              Computes the mathematical entropy as: -\rho s / (\gamma - 1)
!           ******************************************************************
!
            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               p = Pressure( e % storage % Q(:,i,j,k) )
               s = ( log(p) - thermodynamics % gamma * log(e % storage % Q(IRHO,i,j,k)) )
               ms = - e % storage % Q(IRHO,i,j,k) * s / thermodynamics % gammaMinus1
               val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * ms
            end do            ; end do           ; end do

          case ( ENTROPY_RATE )
!
!           ************************************************************
!              Computes the specific entropy integral time derivative
!           ************************************************************
!
            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               call NSGradientVariables_ENTROPY(NCONS, NGRAD, e % storage % Q(:,i,j,k), EntropyVars)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * dot_product(e % storage % QDot(:,i,j,k),EntropyVars)
            end do            ; end do           ; end do


         case (ENTROPY_BALANCE)
!
!           ****************************************************************************
!              Computes the specific entropy integral time derivative minus viscous work
!           ****************************************************************************
!
            select case(grad_vars)
            case(GRADVARS_STATE)
               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  call NSGradientVariables_ENTROPY(NCONS, NGRAD, e % storage % Q(:,i,j,k), EntropyVars)
                  call ViscousFlux_STATE(NCONS, NGRAD, e % storage % Q(:,i,j,k), e % storage % U_x(:,i,j,k), e % storage % U_y(:,i,j,k), &
                                    e % storage % U_z(:,i,j,k), e % storage % mu_ns(1,i,j,k), 0.0_RP, e % storage % mu_ns(2,i,j,k), ViscFlux)
                  val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * (dot_product(e % storage % QDot(:,i,j,k),EntropyVars) + &
                     sum(ViscFlux(:,1)*e % storage % U_x(:,i,j,k)+ViscFlux(:,2)*e % storage % U_y(:,i,j,k)+ViscFlux(:,3)*e % storage % U_z(:,i,j,k)))
               end do            ; end do           ; end do

            case(GRADVARS_ENTROPY)
               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  call NSGradientVariables_ENTROPY(NCONS, NGRAD, e % storage % Q(:,i,j,k), EntropyVars)
                  call ViscousFlux_ENTROPY(NCONS, NGRAD, e % storage % Q(:,i,j,k), e % storage % U_x(:,i,j,k), e % storage % U_y(:,i,j,k), &
                                    e % storage % U_z(:,i,j,k), e % storage % mu_ns(1,i,j,k), 0.0_RP, e % storage % mu_ns(2,i,j,k), ViscFlux)
                  val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * (dot_product(e % storage % QDot(:,i,j,k),EntropyVars) + &
                     sum(ViscFlux(:,1)*e % storage % U_x(:,i,j,k)+ViscFlux(:,2)*e % storage % U_y(:,i,j,k)+ViscFlux(:,3)*e % storage % U_z(:,i,j,k)))
               end do            ; end do           ; end do

               val = val + e % storage % artificialDiss

            case(GRADVARS_ENERGY)
               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  call NSGradientVariables_ENTROPY(NCONS, NGRAD, e % storage % Q(:,i,j,k), EntropyVars)
                  call ViscousFlux_ENERGY(NCONS, NGRAD, e % storage % Q(:,i,j,k), e % storage % U_x(:,i,j,k), e % storage % U_y(:,i,j,k), &
                                    e % storage % U_z(:,i,j,k), e % storage % mu_ns(1,i,j,k), 0.0_RP, e % storage % mu_ns(2,i,j,k), ViscFlux)
                  val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * (dot_product(e % storage % QDot(:,i,j,k),EntropyVars) + &
                     sum(ViscFlux(:,1)*e % storage % U_x(:,i,j,k)+ViscFlux(:,2)*e % storage % U_y(:,i,j,k)+ViscFlux(:,3)*e % storage % U_z(:,i,j,k)))
               end do            ; end do           ; end do

            end select

         case (ARTIFICIAL_DISSIPATION)
            val = val + e % storage % artificialDiss

         case ( INTERNAL_ENERGY )
            !
            !           ***********************************
            !              Computes the internal energy
            !              integral:
            !              \rho e = \int \rho e dV
            !           ***********************************
            !

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * e % storage % Q(IRHOE,i,j,k)
            end do            ; end do           ; end do

         case ( RELATIVE_KINETIC_ENERGY )
!
!           ***********************************
!              Computes the relative kinetic energy
!              integral:
!              K = 1/2 * \int \rho (V-V_s)^2 dV
!            V: flow velocity
!            V_s: immersed body velocity
!            (1/2) * rho * (v-v_s)^2 = (1/2) * rho * v^2 + (1/2) * rho * v_s^2 - rho * v * v_s =
!               = (1/2) * (1/rho) * (rho*v)^2 + (1/2) * rho * v_s^2 - (rho*v) * v_s
!               = (1/2) * (1/rho) * (rho*v)^2 + (1/2) * (1/rho) * (rho*v_s)^2 - (1/rho) * (rho*v) * (rho*v_s)
!           ***********************************
!

            KinEn =         POW2(e % storage % Q(IRHOU,:,:,:))
            KinEn = KinEn + POW2( e % storage % Q(IRHOV,:,:,:) )
            KinEn = KinEn + POW2( e % storage % Q(IRHOW,:,:,:) )
            KinEn = 0.5_RP * KinEn / e % storage % Q(IRHO,:,:,:)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               en = kinEn(i,j,k)
               if( e % isInsideBody(i,j,k) ) then
                  if( IBM% stl(e% STL(i,j,k))% move ) then
                     Q_target = IBM % MaskVelocity( e % storage% Q(:,i,j,k), NCONS, e % STL(i,j,k), e % geom% x(:,i,j,k), time )
                     en = en + (0.5_rp / e % storage % Q(IRHO,i,j,k)) * (POW2(Q_target(IRHOU)) + POW2(Q_target(IRHOV)) + POW2(Q_target(IRHOW))) - &
                        - (1.0_rp / e % storage % Q(IRHO,i,j,k)) * (e % storage % Q(IRHOU,i,j,k) * Q_target(IRHOU) + e % storage % Q(IRHOV,i,j,k) * Q_target(IRHOV) + e % storage % Q(IRHOW,i,j,k) * Q_target(IRHOW))
                  end if
               end if
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * en
            end do            ; end do           ; end do

         case ( RELATIVE_TOTAL_ENERGY )
!
!           ***********************************
!              Computes the relative total energy
!              integral:
!              E_rel =  \int E_rel dV
!              E_rel = E - v_s * (rho*v) + 0.5 * rho * v_s^2
!            V: flow velocity
!            V_s: immersed body velocity
!           ***********************************
!

           do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               en = e % storage % Q(IRHOE,i,j,k)
               if( e % isInsideBody(i,j,k) ) then
                  if( IBM % stl(e % STL(i,j,k)) % move ) then
                     Q_target = IBM % MaskVelocity( e % storage% Q(:,i,j,k), NCONS, e % STL(i,j,k), e % geom% x(:,i,j,k), time )
                     en = en + 0.5_rp / Q_target(IRHO) * (POW2(Q_target(IRHOU)) + POW2(Q_target(IRHOV)) + POW2(Q_target(IRHOW))) - &
                        - (1.0_rp / e % storage % Q(IRHO,i,j,k)) * (e % storage % Q(IRHOU,i,j,k) * Q_target(IRHOU) + e % storage % Q(IRHOV,i,j,k) * Q_target(IRHOV) + e % storage % Q(IRHOW,i,j,k) * Q_target(IRHOW))
                  end if
               end if
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * en 
            end do            ; end do           ; end do
#endif

#if defined(INCNS)
         case (MASS)
            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * e % storage % Q(INSRHO,i,j,k)
            end do            ; end do           ; end do

         case (ENTROPY)
            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               entr =   0.5_RP*(sum(POW2(e % storage % Q(INSRHOU:INSRHOW,i,j,k))))/e % storage % Q(INSRHO,i,j,k) &
                         + 0.5_RP * POW2(e % storage % Q(INSP,i,j,k)) / thermodynamics % rho0c02
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * entr
            end do            ; end do           ; end do

         case (KINETIC_ENERGY_RATE)
!
!           ***********************************
!              Computes the kinetic energy
!              time derivative:
!              K_t = (d/dt)\int \rho V^2 dV
!           ***********************************
!
            uvw = e % storage % Q(INSRHOU,:,:,:) / e % storage % Q(INSRHO,:,:,:)
            KinEn = uvw * e % storage % QDot(INSRHOU,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(INSRHO,:,:,:)

            uvw = e % storage % Q(INSRHOV,:,:,:) / e % storage % Q(INSRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(INSRHOV,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(INSRHO,:,:,:)

            uvw = e % storage % Q(INSRHOW,:,:,:) / e % storage % Q(INSRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(INSRHOW,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(INSRHO,:,:,:)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * kinEn(i,j,k)
            end do            ; end do           ; end do

         case (ENTROPY_RATE)
!
!           ***********************************
!              Computes the kinetic energy
!              time derivative:
!              K_t = (d/dt)\int \rho V^2 dV
!           ***********************************
!
            uvw = e % storage % Q(INSRHOU,:,:,:) / e % storage % Q(INSRHO,:,:,:)
            KinEn = uvw * e % storage % QDot(INSRHOU,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(INSRHO,:,:,:)

            uvw = e % storage % Q(INSRHOV,:,:,:) / e % storage % Q(INSRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(INSRHOV,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(INSRHO,:,:,:)

            uvw = e % storage % Q(INSRHOW,:,:,:) / e % storage % Q(INSRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(INSRHOW,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(INSRHO,:,:,:)

            KinEn = KinEn + (1.0_RP/thermodynamics % rho0c02)*e % storage % Q(INSP,:,:,:)*e % storage % QDot(INSP,:,:,:)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * kinEn(i,j,k)
            end do            ; end do           ; end do


#endif
#ifdef MULTIPHASE
         case (ENTROPY_RATE)
!
!           ***********************************
!              Computes the kinetic energy
!              time derivative:
!              K_t = (d/dt)\int \rho V^2 dV
!           ***********************************
!
            KinEn = e % storage % QDot(IMC,:,:,:) * e % storage % mu(1,:,:,:)
            KinEn = KinEn + e % storage % Q(IMSQRHOU,:,:,:)*e % storage % QDot(IMSQRHOU,:,:,:)
            KinEn = KinEn + e % storage % Q(IMSQRHOV,:,:,:)*e % storage % QDot(IMSQRHOV,:,:,:)
            KinEn = KinEn + e % storage % Q(IMSQRHOW,:,:,:)*e % storage % QDot(IMSQRHOW,:,:,:)
            !KinEn = KinEn + dimensionless % Ma2*e % storage % Q(IMP,:,:,:)*e % storage % QDot(IMP,:,:,:)
            Ma2 = dimensionless % Ma2(1) * e % storage % Q(IMC,:,:,:) + dimensionless % Ma2(2) * (1.0_RP - e % storage % Q(IMC,:,:,:))
            KinEn = KinEn + Ma2*e % storage % Q(IMP,:,:,:)*e % storage % QDot(IMP,:,:,:)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * kinEn(i,j,k)
            end do            ; end do           ; end do


         case (ENTROPY_BALANCE)
!
!           ***********************************
!              Computes the kinetic energy
!              time derivative:
!              K_t = (d/dt)\int \rho V^2 dV
!           ***********************************
!
            KinEn = e % storage % QDot(IMC,:,:,:) * e % storage % mu(1,:,:,:)
            KinEn = KinEn + e % storage % Q(IMSQRHOU,:,:,:)*e % storage % QDot(IMSQRHOU,:,:,:)
            KinEn = KinEn + e % storage % Q(IMSQRHOV,:,:,:)*e % storage % QDot(IMSQRHOV,:,:,:)
            KinEn = KinEn + e % storage % Q(IMSQRHOW,:,:,:)*e % storage % QDot(IMSQRHOW,:,:,:)
            !KinEn = KinEn + dimensionless % Ma2*e % storage % Q(IMP,:,:,:)*e % storage % QDot(IMP,:,:,:)
            Ma2 = dimensionless % Ma2(1) * e % storage % Q(IMC,:,:,:) + dimensionless % Ma2(2) * (1.0_RP - e % storage % Q(IMC,:,:,:))
            KinEn = KinEn + Ma2*e % storage % Q(IMP,:,:,:)*e % storage % QDot(IMP,:,:,:)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)

               mu = max(min(e % storage % Q(IMC,i,j,k),1.0_RP),0.0_RP)
               mu = dimensionless % mu(2) + (dimensionless % mu(1) - dimensionless % mu(2))*mu
               Strain(1,1) = e % storage % U_x(IGU,i,j,k)
               Strain(2,2) = e % storage % U_y(IGV,i,j,k)
               Strain(3,3) = e % storage % U_z(IGW,i,j,k)

               Strain(1,2) = 0.5_RP * (e % storage % U_x(IGV,i,j,k) + e % storage % U_y(IGU,i,j,k))
               Strain(1,3) = 0.5_RP * (e % storage % U_x(IGW,i,j,k) + e % storage % U_z(IGU,i,j,k))
               Strain(2,3) = 0.5_RP * (e % storage % U_y(IGW,i,j,k) + e % storage % U_z(IGV,i,j,k))

               Strain(2,1) = Strain(1,2)
               Strain(3,1) = Strain(1,3)
               Strain(3,2) = Strain(2,3)
!               Strain = 0.0_RP

               val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * (kinEn(i,j,k) + 2.0_RP*mu*sum(Strain*Strain)  &
                        + multiphase % M0 * (POW2(e % storage % U_x(IGMU,i,j,k)) + POW2(e % storage % U_y(IGMU,i,j,k)) + POW2(e % storage % U_z(IGMU,i,j,k)) ) )

            end do            ; end do           ; end do

         case (PHASE2_AREA)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val + wx(i)*wy(j)*wz(k)*e % geom % jacobian(i,j,k)*(1.0_RP-e % storage % Q(IMC,i,j,k))
            end do            ; end do           ; end do

         case (PHASE2_XCOG)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val + wx(i)*wy(j)*wz(k)*e % geom % jacobian(i,j,k)*e % geom % x(IX,i,j,k)*(1.0_RP-e % storage % Q(IMC,i,j,k))
            end do            ; end do           ; end do

            val = val

         case (PHASE2_XVEL)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               rho = dimensionless % rho(1) * e % storage % Q(IMC,i,j,k) + dimensionless % rho(2) * (1.0_RP - e % storage % Q(IMC,i,j,k))
               val = val + wx(i)*wy(j)*wz(k)*e % geom % jacobian(i,j,k)*e % storage % Q(IMSQRHOU,i,j,k)*(1.0_RP-e % storage % Q(IMC,i,j,k)) / sqrt(rho)
            end do            ; end do           ; end do


#endif
#if defined(CAHNHILLIARD)
         case (FREE_ENERGY)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               call QuarticDWP(e % storage % c(1,i,j,k), fchem)
               free_en = fchem + 0.5_RP * POW2(multiphase % eps) *(   POW2(e % storage % c_x(1,i,j,k)) &
                                                                       + POW2(e % storage % c_y(1,i,j,k)) &
                                                                       + POW2(e % storage % c_z(1,i,j,k)) )
               val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * free_en
            end do            ; end do           ; end do

#endif
#if defined(ACOUSTIC)
         case (ACOUSTIC_ENERGY)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                val = val + wx(i)*wy(j)*wz(k)*e % geom % jacobian(i,j,k)* ( &
                      0.5_RP * e % storage % Qbase(ICAARHO,i,j,k) * ( POW2(e %storage % Q(ICAAU,i,j,k)) &
                                                                    + POW2(e %storage % Q(ICAAV,i,j,k)) &
                                                                    + POW2(e %storage % Q(ICAAW,i,j,k)) ) &
                    + 0.5_RP / (thermodynamics % gamma * e % storage % Qbase(ICAAP,i,j,k)) * POW2(e %storage % Q(ICAAP,i,j,k)) )
            end do            ; end do           ; end do
#endif
         end select

         end associate

      end function ScalarVolumeIntegral_Local
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!           VECTOR INTEGRALS PROCEDURES
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      function VectorVolumeIntegral(mesh, integralType, num_of_vars) result(val)
!
!        -----------------------------------------------------------
!           This function computes vector integrals, that is, those
!           in the form:
!                 val = \int \vec{v} dx
!           Implemented integrals are:
!              * VELOCITY
!              * MOMENTUM
!        -----------------------------------------------------------
!
         implicit none
         class(HexMesh),      intent(in)  :: mesh
         integer,             intent(in)  :: integralType
         integer,             intent(in)  :: num_of_vars
         real(kind=RP)                    :: val(num_of_vars)
!
!        ---------------
!        Local variables
!        ---------------
!
         logical       :: fiveVars
         integer       :: eID, ierr
         real(kind=RP) :: localVal(num_of_vars)
         real(kind=RP) :: valAux(num_of_vars)
         real(kind=RP) :: val1, val2, val3, val4, val5

         if (num_of_vars == 5) then    ! Ugly hack.. But only way to make it work with ifort....
            fiveVars = .TRUE.
         else
            fiveVars = .FALSE.
         end if

!
!        Initialization
!        --------------
         val1 = 0.0_RP
         val2 = 0.0_RP
         val3 = 0.0_RP
         val4 = 0.0_RP
         val5 = 0.0_RP
!
!        Loop the mesh
!        -------------
!$omp parallel do reduction(+:val1,val2,val3,val4,val5) private(valAux) schedule(guided)
         do eID = 1, mesh % no_of_elements
!
!           Compute the integral
!           --------------------
            valAux = VectorVolumeIntegral_Local(mesh % elements(eID), integralType  , num_of_vars  )

            val1 = val1 + valAux(1)
            val2 = val2 + valAux(2)
            val3 = val3 + valAux(3)

            if (fiveVars) then
               val4 = val4 + valAux(4)
               val5 = val5 + valAux(5)
            end if
         end do
!$omp end parallel do

         val(1:3) = [val1, val2, val3]

         if (fiveVars) val(4:5) = [val4, val5]

#ifdef _HAS_MPI_
         localVal = val
         call mpi_allreduce(localVal, val, num_of_vars, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      end function VectorVolumeIntegral
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      function VectorVolumeIntegral_Local(e, integralType, num_of_vars) result(val)
         implicit none
         !-arguments---------------------------------------------------
         class(Element),      target, intent(in)     :: e
         integer,                     intent(in)     :: integralType
         integer                    , intent(in)     :: num_of_vars
         real(kind=RP)                               :: val(num_of_vars)
         !-local-variables---------------------------------------------
         integer     :: Nel(3)    ! Element polynomial order
         integer     :: i, j, k
         !-------------------------------------------------------------

         Nel = e % Nxyz

         associate ( wx => NodalStorage(e % Nxyz(1)) % w, &
                     wy => NodalStorage(e % Nxyz(2)) % w, &
                     wz => NodalStorage(e % Nxyz(3)) % w    )
!
!        Initialization
!        --------------
         val = 0.0_RP
!
!        Perform the numerical integration
!        ---------------------------------
         select case ( integralType )

#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
            case ( VELOCITY )

               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  val = val +   wx(i) * wy(j) * wz(k) * e % storage % Q(IRHOU:IRHOW,i,j,k) / e % storage % Q(IRHO,i,j,k) * e % geom % jacobian(i,j,k)
               end do            ; end do           ; end do

            case ( MOMENTUM )

               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  val = val +   wx(i) * wy(j) * wz(k) * e % storage % Q(IRHOU:IRHOW,i,j,k) * e % geom % jacobian(i,j,k)
               end do            ; end do           ; end do

            case ( SOURCE )

               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  val = val +   wx(i) * wy(j) * wz(k) * e % storage % S_NS(1:num_of_vars,i,j,k) * e % geom % jacobian(i,j,k)
               end do            ; end do           ; end do

            case ( PSOURCE )

               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  val = val +   wx(i) * wy(j) * wz(k) * e % storage % S_NSP(1:num_of_vars,i,j,k) * e % geom % jacobian(i,j,k)
               end do            ; end do           ; end do
#endif
#if defined(INCNS)
            case ( SOURCE )

               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  val = val +   wx(i) * wy(j) * wz(k) * e % storage % S_NS(1:num_of_vars,i,j,k) * e % geom % jacobian(i,j,k)
               end do            ; end do           ; end do
#endif
#if defined(MULTIPHASE)
            case ( SOURCE )

               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  val = val +   wx(i) * wy(j) * wz(k) * e % storage % S_NS(1:num_of_vars,i,j,k) * e % geom % jacobian(i,j,k)
               end do            ; end do           ; end do
#endif
#if defined(ACOUSTIC)
            case ( SOURCE )

               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  val = val +   wx(i) * wy(j) * wz(k) * e % storage % S_NS(1:num_of_vars,i,j,k) * e % geom % jacobian(i,j,k)
               end do            ; end do           ; end do
#endif

            case default

               write(STD_OUT,'(A,A)') 'VectorVolumeIntegral :: ERROR: Not defined integral type'
               error stop 99

         end select

         end associate
      end function VectorVolumeIntegral_Local

      subroutine GetPressureLocalGradient(N, p, Ja_xi, Ja_eta, Ja_zeta, grad_Mp, M_grad_p)
         implicit none
         integer, intent(in)        :: N(3)
         real(kind=RP), intent(in)  :: p(0:N(1),0:N(2),0:N(3))
         real(kind=RP), intent(in)  :: Ja_xi(1:NDIM,0:N(1),0:N(2),0:N(3))
         real(kind=RP), intent(in)  :: Ja_eta(1:NDIM,0:N(1),0:N(2),0:N(3))
         real(kind=RP), intent(in)  :: Ja_zeta(1:NDIM,0:N(1),0:N(2),0:N(3))
         real(kind=RP), intent(out) :: grad_Mp(1:NDIM,0:N(1),0:N(2),0:N(3))
         real(kind=RP), intent(out) :: M_grad_p(1:NDIM,0:N(1),0:N(2),0:N(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k, l
         type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta

         spAxi   => NodalStorage(N(1))
         spAeta  => NodalStorage(N(2))
         spAzeta => NodalStorage(N(3))

         grad_Mp = 0.0_RP
         M_grad_p = 0.0_RP
         do k = 0, N(3)   ; do j = 0, N(2) ; do l = 0, N(1) ; do i = 0, N(1)
            grad_Mp(:,i,j,k)  = grad_Mp(:,i,j,k)  + p(l,j,k) * Ja_xi(:,l,j,k) * spAxi % D(i,l)
            M_grad_p(:,i,j,k) = M_grad_p(:,i,j,k) + p(l,j,k) * Ja_xi(:,i,j,k) * spAxi % D(i,l)
         end do           ; end do         ; end do         ; end do

         do k = 0, N(3)   ; do l = 0, N(2) ; do j = 0, N(2) ; do i = 0, N(1)
            grad_Mp(:,i,j,k)  = grad_Mp(:,i,j,k)  + p(i,l,k) * Ja_eta(:,i,l,k) * spAeta % D(j,l)
            M_grad_p(:,i,j,k) = M_grad_p(:,i,j,k) + p(i,l,k) * Ja_eta(:,i,j,k) * spAeta % D(j,l)
         end do           ; end do         ; end do         ; end do

         do l = 0, N(3)   ; do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            grad_Mp(:,i,j,k)  = grad_Mp(:,i,j,k)  + p(i,j,l) * Ja_zeta(:,i,j,l) * spAzeta % D(k,l)
            M_grad_p(:,i,j,k) = M_grad_p(:,i,j,k) + p(i,j,l) * Ja_zeta(:,i,j,k) * spAzeta % D(k,l)
         end do           ; end do         ; end do         ; end do

         nullify (spAxi, spAeta, spAzeta)

      end subroutine GetPressureLocalGradient

      subroutine GetSensorRange(mesh, minSensor, maxSensor)
         implicit none
         class(HexMesh), intent(in)  :: mesh
         real(RP),       intent(out) :: minSensor
         real(RP),       intent(out) :: maxSensor
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: ielem
         integer  :: ierr


         minSensor =  huge(1.0_RP)/10.0_RP
         maxSensor = -huge(1.0_RP)/10.0_RP

!! !$omp parallel do schedule(static) private(ielem)
         ! do ielem = 1, mesh % no_of_elements
         !    minSensor = min(minSensor, mesh % elements(ielem) % storage % sensor)
         !    maxSensor = max(maxSensor, mesh % elements(ielem) % storage % sensor)
         ! end do
!! !$omp end parallel do

!$omp parallel shared(maxSensor, minSensor, mesh) default(private)
!$omp do reduction(max:maxSensor) schedule(runtime)
         DO ielem = 1, mesh % no_of_elements
            maxSensor = max(maxSensor, mesh % elements(ielem) % storage % sensor)
         END DO
!$omp end do
   
!$omp do reduction(min:minSensor) schedule(runtime)
         DO ielem = 1, mesh % no_of_elements
            minSensor = min(minSensor, mesh % elements(ielem) % storage % sensor)
         END DO
!$omp end do
!$omp end parallel

#ifdef _HAS_MPI_
         call MPI_MinMax(minSensor, maxSensor)
#endif

      end subroutine GetSensorRange

end module VolumeIntegrals
