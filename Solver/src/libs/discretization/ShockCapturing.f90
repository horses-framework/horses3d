#include "Includes.h"
module ShockCapturing

   use SMConstants,                only: RP, STD_OUT, LINE_LENGTH, NDIM, IX, IY, IZ
   use PhysicsStorage,             only: NCONS, NGRAD
   use FluidData,                  only: dimensionless
   use Utilities,                  only: toLower
   use DGSEMClass,                 only: ComputeTimeDerivative_f, DGSem
   use HexMeshClass,               only: HexMesh
   use ElementClass,               only: Element
#if !defined(SPALARTALMARAS)
   use LESModels,                  only: Smagorinsky_t
#endif
   use SpectralVanishingViscosity, only: SVV, InitializeSVV
   use DGIntegrals,                only: ScalarWeakIntegrals, ScalarStrongIntegrals

   use ShockCapturingKeywords
   use SCsensorClass, only: SCsensor_t, Set_SCsensor, Destruct_SCsensor

   implicit none

   public :: Initialize_ShockCapturing
   public :: ShockCapturingDriver

   type ArtificialViscosity_t

      integer,  private :: region           !< Sensor region where it acts (1 or 2)
      integer,  private :: updateMethod     !< Method to compute the viscosity coefficient
      real(RP), private :: mu1              !< First viscosity parameter (low)
      real(RP), private :: alpha1           !< Second viscosity parameter (low)
      real(RP), private :: mu2              !< First viscosity parameter (high)
      real(RP), private :: alpha2           !< Second viscosity parameter (high)
      real(RP), private :: mu2alpha         !< Ratio alpha/mu
      logical,  private :: alphaIsPropToMu  !< .true. if alpha/mu is defined
#if !defined (SPALARTALMARAS)
      type(Smagorinsky_t), private :: Smagorinsky  !< For automatic viscosity
#endif

      contains

         procedure :: Initialize => AV_initialize
         procedure :: Viscosity  => AV_viscosity
         procedure :: Describe   => AV_describe

   end type ArtificialViscosity_t

   type SCdriver_t

      logical  :: isActive = .false.  !< On/Off flag

      type(SCsensor_t),             private              :: sensor   !< Sensor
      class(ArtificialViscosity_t), private, allocatable :: method1  !< First SC method
      class(ArtificialViscosity_t), private, allocatable :: method2  !< Second SC method

      contains

         procedure :: Detect           => SC_detect
         procedure :: ComputeViscosity => SC_viscosity
         procedure :: Describe         => SC_describe

         final :: SC_destruct

   end type SCdriver_t

   type, extends(ArtificialViscosity_t) :: SC_NoSVV_t

      integer,                                 private :: fluxType
      procedure(Viscous_Int), nopass, pointer, private :: ViscousFlux => null()

      contains

         procedure :: Initialize => NoSVV_initialize
         procedure :: Viscosity  => NoSVV_viscosity
         procedure :: Describe   => NoSVV_describe

         final :: NoSVV_destruct

   end type SC_NoSVV_t

   type, extends(ArtificialViscosity_t) :: SC_SVV_t

      real(RP), private :: sqrt_mu1
      real(RP), private :: sqrt_alpha1
      real(RP), private :: sqrt_mu2
      real(RP), private :: sqrt_alpha2
      real(RP), private :: sqrt_mu2alpha

      contains

         procedure :: Initialize => SVV_initialize
         procedure :: Viscosity  => SVV_viscosity
         procedure :: Describe   => SVV_describe

   end type SC_SVV_t
!
!  Interfaces
!  ----------
   abstract interface
      pure subroutine Viscous_Int(nEqn, nGradEqn, Q, Q_x, Q_y, Q_z, mu, beta, kappa, F)
         import RP, NDIM
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: nGradEqn
         real(kind=RP), intent(in)  :: Q   (1:nEqn    )
         real(kind=RP), intent(in)  :: Q_x (1:nGradEqn)
         real(kind=RP), intent(in)  :: Q_y (1:nGradEqn)
         real(kind=RP), intent(in)  :: Q_z (1:nGradEqn)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: beta
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(out) :: F(1:nEqn, 1:NDIM)
      end subroutine Viscous_Int
   end interface

   type(SCdriver_t), allocatable :: ShockCapturingDriver
!
!  ========
   contains
!  ========
!
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Initializer
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Initialize_ShockCapturing(self, controlVariables, sem, &
                                        TimeDerivative, TimeDerivativeIsolated)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(SCdriver_t), allocatable, intent(inout) :: self
      class(FTValueDictionary),      intent(in)    :: controlVariables
      class(DGSem),                  intent(inout) :: sem
      procedure(ComputeTimeDerivative_f)           :: TimeDerivative
      procedure(ComputeTimeDerivative_f)           :: TimeDerivativeIsolated
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: method
      integer                       :: minSteps

!
!     Check if shock-capturing is requested
!     -------------------------------------
      allocate(self)

      if (controlVariables % containsKey(SC_KEY)) then
         self % isActive = controlVariables % logicalValueForKey(SC_KEY)
      else
         self % isActive = .false.
      end if

      if (.not. self % isActive) then
         return
      end if
!
!     Shock-capturing methods
!     -----------------------
      if (controlVariables % containsKey(SC_METHOD1_KEY)) then
         method = controlVariables % stringValueForKey(SC_METHOD1_KEY, LINE_LENGTH)
      else
         method = SC_NO_VAL
      end if
      call toLower(method)

      select case (trim(method))
      case (SC_NOSVV_VAL)
         allocate(SC_NoSVV_t :: self % method1)

      case (SC_SVV_VAL)
         allocate(SC_SVV_t :: self % method1)

      case (SC_NO_VAL)
         safedeallocate(self % method1)

      case default
         write(STD_OUT,*) 'ERROR. Unavailable first shock-capturing method. Options are:'
         write(STD_OUT,*) '   * ', SC_NO_VAL
         write(STD_OUT,*) '   * ', SC_NOSVV_VAL
         write(STD_OUT,*) '   * ', SC_SVV_VAL
         error stop

      end select

      if (controlVariables % containsKey(SC_METHOD2_KEY)) then
         method = controlVariables % stringValueForKey(SC_METHOD2_KEY, LINE_LENGTH)
      else
         method = SC_NO_VAL
      end if
      call toLower(method)

      select case (trim(method))
      case (SC_NOSVV_VAL)
         allocate(SC_NoSVV_t :: self % method2)

      case (SC_SVV_VAL)
         allocate(SC_SVV_t :: self % method2)

      case (SC_NO_VAL)
         safedeallocate(self % method2)

      case default
         write(STD_OUT,*) 'ERROR. Unavailable second shock-capturing method. Options are:'
         write(STD_OUT,*) '   * ', SC_NO_VAL
         write(STD_OUT,*) '   * ', SC_NOSVV_VAL
         write(STD_OUT,*) '   * ', SC_SVV_VAL
         error stop

      end select
!
!     Initialize viscous and hyperbolic terms
!     ---------------------------------------
      if (allocated(self % method1)) call self % method1 % Initialize(controlVariables, sem % mesh, 1)
      if (allocated(self % method2)) call self % method2 % Initialize(controlVariables, sem % mesh, 2)
!
!     Sensor 'inertia'
!     ----------------
      if (controlVariables % containsKey(SC_SENSOR_INERTIA_KEY)) then
         minSteps = controlVariables % realValueForKey(SC_SENSOR_INERTIA_KEY)
         if (minSteps < 1) then
            write(STD_OUT,*) 'ERROR. Sensor inertia must be at least 1.'
            error stop
         end if
      else
         minSteps = 1
      end if
!
!     Construct sensor
!     ----------------
      call Set_SCsensor(self % sensor, controlVariables, sem, minSteps, &
                        TimeDerivative, TimeDerivativeIsolated)

   end subroutine Initialize_ShockCapturing
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Base class
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_detect(self, sem, t)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCdriver_t), intent(inout) :: self
      type(DGSem),       intent(inout) :: sem
      real(RP),          intent(in)    :: t


      call self % sensor % Compute(sem, t)

   end subroutine SC_detect
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_viscosity(self, mesh, e, SCflux)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCdriver_t), intent(in)    :: self
      type(HexMesh),     intent(inout) :: mesh
      type(Element),     intent(inout) :: e
      real(RP),          intent(out)   :: SCflux(1:NCONS,       &
                                                 0:e % Nxyz(1), &
                                                 0:e % Nxyz(2), &
                                                 0:e % Nxyz(3), &
                                                 1:NDIM)
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: switch
      logical  :: updated


      switch = e % storage % sensor
      updated = .false.

      if (switch >= 1.0_RP) then
         if (allocated(self % method2)) then
            call self % method2 % Viscosity(mesh, e, switch, SCflux)
            updated = .true.
         end if

      elseif (switch > 0.0_RP) then
         if (allocated(self % method1)) then
            call self % method1 % Viscosity(mesh, e, switch, SCflux)
            updated = .true.
         end if

      end if

      if (.not. updated) then
         SCflux = 0.0_RP
         e % storage % artificialDiss = 0.0_RP
      end if

   end subroutine SC_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_describe(self)
!
!     -------
!     Modules
!     -------
      use MPI_Process_Info, only: MPI_Process
      use Headers,          only: Subsection_Header
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCdriver_t), intent(in) :: self


      if (.not. MPI_Process % isRoot .or. .not. self % isActive) return

      write(STD_OUT, "(/)")
      call Subsection_Header("Shock-Capturing")

      call self % sensor % Describe()

      if (allocated(self % method1)) then
         write(STD_OUT,*) ""
         write(STD_OUT,"(30X,A,A30)") "=>", "First method"
         call self % method1 % Describe()
      end if

      if (allocated(self % method2)) then
         write(STD_OUT,*) ""
         write(STD_OUT,"(30X,A,A30)") "=>", "Second method"
         call self % method2 % Describe()
      end if

   end subroutine SC_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_destruct(self)
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(SCdriver_t), intent(inout) :: self


      self%isActive = .false.

      call Destruct_SCsensor(self % sensor)

      if (allocated(self % method1)) deallocate(self % method1)
      if (allocated(self % method2)) deallocate(self % method2)

   end subroutine SC_destruct
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Shock-capturing methods base class
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine AV_initialize(self, controlVariables, mesh, region)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(ArtificialViscosity_t), intent(inout) :: self
      type(FTValueDictionary),      intent(in)    :: controlVariables
      type(HexMesh),                intent(inout) :: mesh
      integer,                      intent(in)    :: region
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: update

!
!     Viscosity values (mu and alpha)
!     -------------------------------
      if (region == 1) then

         if (controlVariables % containsKey(SC_MU1_KEY)) then
            self % mu1 = controlVariables % doublePrecisionValueForKey(SC_MU1_KEY)
         else
            self % mu1 = 0.0_RP
         end if

         if (controlVariables % containsKey(SC_MU2_KEY)) then
            self % mu2 = controlVariables % doublePrecisionValueForKey(SC_MU2_KEY)
         else
            self % mu2 = self % mu1
         end if

      else

         if (controlVariables % containsKey(SC_MU2_KEY)) then
            self % mu2 = controlVariables % doublePrecisionValueForKey(SC_MU2_KEY)
         else
            self % mu2 = 0.0_RP
         end if
         self % mu1 = self % mu2

      end if

      if (controlVariables % containsKey(SC_ALPHA_MU_KEY)) then

         self % alphaIsPropToMu = .true.
         self % mu2alpha        = controlVariables % doublePrecisionValueForKey(SC_ALPHA_MU_KEY)
         self % alpha1          = self % mu2alpha * self % mu1
         self % alpha2          = self % mu2alpha * self % mu2

      else

         self % alphaIsPropToMu = .false.

         if (region == 1) then

            if (controlVariables % containsKey(SC_ALPHA1_KEY)) then
               self % alpha1 = controlVariables % doublePrecisionValueForKey(SC_ALPHA1_KEY)
            else
               self % alpha1 = 0.0_RP
            end if

            if (controlVariables % containsKey(SC_ALPHA2_KEY)) then
               self % alpha2 = controlVariables % doublePrecisionValueForKey(SC_ALPHA2_KEY)
            else
               self % alpha2 = self % alpha1
            end if

         else

            if (controlVariables % containsKey(SC_ALPHA2_KEY)) then
               self % alpha2 = controlVariables % doublePrecisionValueForKey(SC_ALPHA2_KEY)
            else
               self % alpha2 = 0.0_RP
            end if
            self % alpha1 = self % alpha2

         end if

      end if
!
!     Viscosity update method
!     -----------------------
      if (region == 1) then

         if (controlVariables % containsKey(SC_UPDATE_KEY)) then
            update = controlVariables % StringValueForKey(SC_UPDATE_KEY, LINE_LENGTH)
            call toLower(update)

            select case (trim(update))
            case (SC_CONST_VAL)
               self % updateMethod = SC_CONST_ID

            case (SC_SENSOR_VAL)
               self % updateMethod = SC_SENSOR_ID

#if !defined (SPALARTALMARAS)
            case (SC_SMAG_VAL)

               self % updateMethod = SC_SMAG_ID
               if (.not. self % alphaIsPropToMu) then
                  write(STD_OUT,*) 'ERROR. Alpha must be proportional to mu when using shock-capturing with LES.'
                  error stop
               end if

               ! TODO: Use the default constructor
               self % Smagorinsky % active = .true.
               self % Smagorinsky % requiresWallDistances = .false.
               self % Smagorinsky % WallModel = 0  ! No wall model
               self % Smagorinsky % CS = self % mu1
#endif

            case default
               write(STD_OUT,*) 'ERROR. Unavailable shock-capturing update strategy. Options are:'
               write(STD_OUT,*) '   * ', SC_CONST_VAL
               write(STD_OUT,*) '   * ', SC_SENSOR_VAL
#if !defined (SPALARTALMARAS)
               write(STD_OUT,*) '   * ', SC_SMAG_VAL
#endif
               error stop

            end select

         else
            self % updateMethod = SC_CONST_ID

         end if

      else

         self % updateMethod = SC_CONST_ID

      end if

   end subroutine AV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine AV_viscosity(self, mesh, e, switch, SCflux)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(ArtificialViscosity_t), intent(in)    :: self
      type(HexMesh),                intent(inout) :: mesh
      type(Element),                intent(inout) :: e
      real(RP),                     intent(in)    :: switch
      real(RP),                     intent(out)   :: SCflux(1:NCONS,       &
                                                            0:e % Nxyz(1), &
                                                            0:e % Nxyz(2), &
                                                            0:e % Nxyz(3), &
                                                            1:NDIM)


      SCflux = 0.0_RP

   end subroutine AV_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine AV_describe(self)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(ArtificialViscosity_t), intent(in) :: self

   end subroutine AV_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Simple artificial viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine NoSVV_initialize(self, controlVariables, mesh, region)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
      use PhysicsStorage, only: grad_vars, GRADVARS_STATE, &
                                GRADVARS_ENTROPY, GRADVARS_ENERGY
      use Physics,        only: ViscousFlux_STATE, ViscousFlux_ENTROPY, &
                                ViscousFlux_ENERGY, GuermondPopovFlux_ENTROPY
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_NoSVV_t),       intent(inout) :: self
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(HexMesh),           intent(inout) :: mesh
      integer,                 intent(in)    :: region
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: flux

!
!     Parent initializer
!     ------------------
      call self % ArtificialViscosity_t % Initialize(controlVariables, mesh, region)
      self % region = region
!
!     Set the flux type
!     -----------------
      if (region == 1 .and. controlVariables % containsKey(SC_VISC_FLUX1_KEY)) then
         flux = controlVariables % stringValueForKey(SC_VISC_FLUX1_KEY, LINE_LENGTH)
      elseif (region == 2 .and. controlVariables % containsKey(SC_VISC_FLUX2_KEY)) then
         flux = controlVariables % stringValueForKey(SC_VISC_FLUX2_KEY, LINE_LENGTH)
      else
         flux = SC_PHYS_VAL
      end if

      call toLower(flux)

      select case (trim(flux))
      case (SC_PHYS_VAL); self % fluxType = SC_PHYS_ID
      case (SC_GP_VAL);   self % fluxType = SC_GP_ID
      case default
         write(STD_OUT,'(A,I1,A)') 'ERROR. Artificial viscosity ', region, &
                                   ' not recognized. Options are:'
         write(STD_OUT,*) '   * ', SC_PHYS_VAL
         write(STD_OUT,*) '   * ', SC_GP_VAL
         error stop
      end select

      select case (self % fluxType)
      case (SC_PHYS_ID)
         select case (grad_vars)
         case (GRADVARS_STATE);   self % ViscousFlux => ViscousFlux_STATE
         case (GRADVARS_ENTROPY); self % ViscousFlux => ViscousFlux_ENTROPY
         case (GRADVARS_ENERGY);  self % ViscousFlux => ViscousFlux_ENERGY
         end select

      case (SC_GP_ID)
         select case (grad_vars)
         case (GRADVARS_ENTROPY); self % ViscousFlux => GuermondPopovFlux_ENTROPY
         case default
            write(STD_OUT,*) "ERROR. Guermond-Popov (2014) artificial ",  &
                              "viscosity is only configured for Entropy ", &
                              "gradient variables"
            error stop
         end select

      end select

   end subroutine NoSVV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine NoSVV_viscosity(self, mesh, e, switch, SCflux)
!
!     --------------------------------------------------------------------------
!     TODO: Introduce alpha viscosity, which probably means reimplementing here
!           all the viscous fluxes of `Physics_NS`...
!     --------------------------------------------------------------------------
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_NoSVV_t), intent(in)    :: self
      type(HexMesh),     intent(inout) :: mesh
      type(Element),     intent(inout) :: e
      real(RP),          intent(in)    :: switch
      real(RP),          intent(out)   :: SCflux(1:NCONS,       &
                                                 0:e % Nxyz(1), &
                                                 0:e % Nxyz(2), &
                                                 0:e % Nxyz(3), &
                                                 1:NDIM)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i
      integer  :: j
      integer  :: k
      integer  :: fIDs(6)
      real(RP) :: delta
      real(RP) :: kappa
      real(RP) :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
      real(RP) :: covariantFlux(1:NCONS, 1:NDIM)


      if (switch > 0.0_RP) then
!
!        Compute viscosity
!        -----------------
         select case (self % updateMethod)
         case (SC_CONST_ID)
            mu = merge(self % mu2, self % mu1, switch >= 1.0_RP) * e % hn

         case (SC_SENSOR_ID)
            mu = (self % mu1 * (1.0_RP-switch) + self % mu2 * switch) * e % hn

#if !defined (SPALARTALMARAS)
         case (SC_SMAG_ID)

            delta = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               call self % Smagorinsky % ComputeViscosity(delta, e % geom % dWall(i,j,k), &
                                                          e % storage % Q(:,i,j,k),       &
                                                          e % storage % U_x(:,i,j,k),     &
                                                          e % storage % U_y(:,i,j,k),     &
                                                          e % storage % U_z(:,i,j,k),     &
                                                          mu(i,j,k))
            end do                ; end do                ; end do
#endif

         end select
!
!        Compute the viscous flux
!        ------------------------
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)

            kappa = dimensionless % mu_to_kappa * mu(i,j,k)
            call self % ViscousFlux(NCONS, NGRAD, e % storage % Q(:,i,j,k), &
                                    e % storage % U_x(:,i,j,k),             &
                                    e % storage % U_y(:,i,j,k),             &
                                    e % storage % U_z(:,i,j,k),             &
                                    mu(i,j,k), 0.0_RP, kappa,               &
                                    covariantflux)

            SCflux(:,i,j,k,IX) = covariantFlux(:,IX) * e % geom % jGradXi(IX,i,j,k) &
                               + covariantFlux(:,IY) * e % geom % jGradXi(IY,i,j,k) &
                               + covariantFlux(:,IZ) * e % geom % jGradXi(IZ,i,j,k)


            SCflux(:,i,j,k,IY) = covariantFlux(:,IX) * e % geom % jGradEta(IX,i,j,k) &
                               + covariantFlux(:,IY) * e % geom % jGradEta(IY,i,j,k) &
                               + covariantFlux(:,IZ) * e % geom % jGradEta(IZ,i,j,k)


            SCflux(:,i,j,k,IZ) = covariantFlux(:,IX) * e % geom % jGradZeta(IX,i,j,k) &
                               + covariantFlux(:,IY) * e % geom % jGradZeta(IY,i,j,k) &
                               + covariantFlux(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)

         end do                ; end do                ; end do

      else

         e % storage % artificialDiss = 0.0_RP
         SCflux = 0.0_RP

      end if
!
!     Project to faces
!     ----------------
      fIDs = e % faceIDs
      if (.not. mesh%nonconforming) then 
      call e % ProlongAviscFluxToFaces(NCONS, SCflux, mesh % faces(fIDs(1)), &
                                                      mesh % faces(fIDs(2)), &
                                                      mesh % faces(fIDs(3)), &
                                                      mesh % faces(fIDs(4)), &
                                                      mesh % faces(fIDs(5)), &
                                                      mesh % faces(fIDs(6)) )
      else 
      call e % ProlongAviscFluxToFaces(NCONS, SCflux, mesh % faces(fIDs(1)), &
                                                      mesh % faces(fIDs(2)), &
                                                      mesh % faces(fIDs(3)), &
                                                      mesh % faces(fIDs(4)), &
                                                      mesh % faces(fIDs(5)), &
                                                      mesh % faces(fIDs(6)), faces=mesh % faces)
      end if 

   end subroutine NoSVV_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine NoSVV_describe(self)
!
!     -------
!     Modules
!     -------
      use MPI_Process_Info, only: MPI_Process
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_NoSVV_t), intent(in) :: self


      if (.not. MPI_Process % isRoot) return

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Viscosity update method: "
      select case (self % updateMethod)
         case (SC_CONST_ID);  write(STD_OUT,"(A)") SC_CONST_VAL
         case (SC_SENSOR_ID); write(STD_OUT,"(A)") SC_SENSOR_VAL
         case (SC_SMAG_ID);   write(STD_OUT,"(A)") SC_SMAG_VAL
      end select

      if (self % updateMethod == SC_SMAG_ID) then
#if !defined (SPALARTALMARAS)
         write(STD_OUT,"(30X,A,A30,1pG10.3)") "->", "LES intensity (CS): ", self % Smagorinsky % CS
#endif
      else
         if (self % region == 1) then
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Mu viscosity 1: ", self % mu1
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Mu viscosity 2: ", self % mu2
         else ! self % region == 2
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Mu viscosity: ", self % mu2
         end if
      end if

      if (self % alphaIsPropToMu) then
         write(STD_OUT,"(30X,A,A30,1pG10.3,A)") "->", "Alpha viscosity: ", self % mu2alpha, "x mu"
      else
         if (self % region == 1) then
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Alpha viscosity 1: ", self % alpha1
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Alpha viscosity 2: ", self % alpha2
         else ! self % region == 2
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Alpha viscosity: ", self % alpha2
         end if
      end if

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Dissipation type: "
      select case (self % fluxType)
         case (SC_PHYS_ID); write(STD_OUT,"(A)") SC_PHYS_VAL
         case (SC_GP_ID);   write(STD_OUT,"(A)") SC_GP_VAL
      end select

   end subroutine NoSVV_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine NoSVV_destruct(self)
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(SC_NoSVV_t), intent(inout) :: self


      if (associated(self % ViscousFlux)) nullify(self % ViscousFlux)

   end subroutine NoSVV_destruct
!
!///////////////////////////////////////////////////////////////////////////////
!
!     SVV filtered viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SVV_initialize(self, controlVariables, mesh, region)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_SVV_t),         intent(inout) :: self
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(HexMesh),           intent(inout) :: mesh
      integer,                 intent(in)    :: region


      ! TODO: Implement it also for region 2, but SVV does not seem very useful there...
      if (region == 2) then
         write(STD_OUT,*) "ERROR. SVV viscosity can be used only in the first region of the sensor."
         error stop
      end if
!
!     Parent initializer
!     ------------------
      call self % ArtificialViscosity_t % Initialize(controlVariables, mesh, region)
      self % region = region
!
!     Set the square root of the viscosities
!     --------------------------------------
      self % sqrt_mu1 = sqrt(self % mu1)
      self % sqrt_mu2 = sqrt(self % mu2)

      if (self % alphaIsPropToMu) then
         self % sqrt_mu2alpha = sqrt(self % mu2alpha)
      else
         self % sqrt_alpha1 = sqrt(self % alpha1)
         self % sqrt_alpha2 = sqrt(self % alpha2)
      end if
!
!     Start the SVV module
!     --------------------
      call InitializeSVV(SVV, controlVariables, mesh)

   end subroutine SVV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SVV_viscosity(self, mesh, e, switch, SCflux)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_SVV_t), intent(in)    :: self
      type(HexMesh),   intent(inout) :: mesh
      type(Element),   intent(inout) :: e
      real(RP),        intent(in)    :: switch
      real(RP),        intent(out)   :: SCflux(1:NCONS,       &
                                               0:e % Nxyz(1), &
                                               0:e % Nxyz(2), &
                                               0:e % Nxyz(3), &
                                               1:NDIM)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i
      integer  :: j
      integer  :: k
      integer  :: fIDs(6)
      real(RP) :: delta
      real(RP) :: salpha
      real(RP) :: sqrt_mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
      real(RP) :: sqrt_alpha(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))


      if (switch > 0.0_RP) then
!
!        Compute viscosities
!        -------------------
         select case (self % updateMethod)
         case (SC_CONST_ID)
            sqrt_mu = merge(self % sqrt_mu2,    self % sqrt_mu1,    switch >= 1.0_RP) * e % hn
            salpha  = merge(self % sqrt_alpha2, self % sqrt_alpha1, switch >= 1.0_RP) * e % hn

         case (SC_SENSOR_ID)
            sqrt_mu = (self % sqrt_mu1 * (1.0_RP-switch) + self % sqrt_mu2 * switch) * e % hn
            salpha  = (self % sqrt_alpha1 * (1.0_RP-switch) + self % sqrt_alpha2 * switch) * e % hn

#if !defined (SPALARTALMARAS)
         case (SC_SMAG_ID)

            delta = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               call self % Smagorinsky % ComputeViscosity(delta, e % geom % dWall(i,j,k), &
                                                          e % storage % Q(:,i,j,k),       &
                                                          e % storage % U_x(:,i,j,k),     &
                                                          e % storage % U_y(:,i,j,k),     &
                                                          e % storage % U_z(:,i,j,k),     &
                                                          sqrt_mu(i,j,k))
               sqrt_mu(i,j,k) = sqrt(sqrt_mu(i,j,k))
            end do                ; end do                ; end do
#endif

         end select

         if (self % alphaIsPropToMu) then
            sqrt_alpha = self % sqrt_mu2alpha * sqrt_mu
         else
            sqrt_alpha = salpha
         end if
!
!        Compute the viscous flux
!        ------------------------
         call SVV % ComputeInnerFluxes(e, sqrt_mu, sqrt_alpha, SCflux)

      else

         e % storage % artificialDiss = 0.0_RP
         SCflux = 0.0_RP

      end if
!
!     Project to faces
!     ----------------
      fIDs = e % faceIDs
      if(.not. mesh%nonconforming) then 
      call e % ProlongAviscFluxToFaces(NCONS, SCflux, mesh % faces(fIDs(1)), &
                                                      mesh % faces(fIDs(2)), &
                                                      mesh % faces(fIDs(3)), &
                                                      mesh % faces(fIDs(4)), &
                                                      mesh % faces(fIDs(5)), &
                                                      mesh % faces(fIDs(6)))
      else 
      call e % ProlongAviscFluxToFaces(NCONS, SCflux, mesh % faces(fIDs(1)), &
                                                      mesh % faces(fIDs(2)), &
                                                      mesh % faces(fIDs(3)), &
                                                      mesh % faces(fIDs(4)), &
                                                      mesh % faces(fIDs(5)), &
                                                      mesh % faces(fIDs(6)), faces=mesh % faces  )
      end if 


   end subroutine SVV_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SVV_describe(self)
!
!     -------
!     Modules
!     -------
      use MPI_Process_Info, only: MPI_Process
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_SVV_t), intent(in) :: self


      if (.not. MPI_Process % isRoot) return

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Viscosity update method: "
      select case (self % updateMethod)
         case (SC_CONST_ID);  write(STD_OUT,"(A)") SC_CONST_VAL
         case (SC_SENSOR_ID); write(STD_OUT,"(A)") SC_SENSOR_VAL
         case (SC_SMAG_ID);   write(STD_OUT,"(A)") SC_SMAG_VAL
      end select

      if (self % updateMethod == SC_SMAG_ID) then
#if !defined (SPALARTALMARAS)
         write(STD_OUT,"(30X,A,A30,1pG10.3)") "->", "LES intensity (CS): ", self % Smagorinsky % CS
#endif
      else
         if (self % region == 1) then
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Mu viscosity 1: ", self % mu1
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Mu viscosity 2: ", self % mu2
         else ! self % region == 2
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Mu viscosity: ", self % mu2
         end if
      end if

      if (self % alphaIsPropToMu) then
         write(STD_OUT,"(30X,A,A30,1pG10.3,A)") "->", "Alpha viscosity: ", self % mu2alpha, "x mu"
      else
         if (self % region == 1) then
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Alpha viscosity 1: ", self % alpha1
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Alpha viscosity 2: ", self % alpha2
         else ! self % region == 2
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Alpha viscosity: ", self % alpha2
         end if
      end if

      call SVV % Describe()

   end subroutine SVV_describe

end module ShockCapturing