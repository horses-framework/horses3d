#include "Includes.h"
      Module Physics_NSKeywordsModule
         IMPLICIT NONE
         INTEGER, PARAMETER :: KEYWORD_LENGTH = 132
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: REFERENCE_TEMPERATURE_KEY      = "reference temperature (k)"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: REFERENCE_PRESSURE_KEY         = "reference pressure (pa)"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: MACH_NUMBER_KEY                = "mach number"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: REYNOLDS_NUMBER_KEY            = "reynolds number"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: PRANDTL_NUMBER_KEY             = "prandtl number"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: TURBULENT_PRANDTL_NUMBER_KEY   = "turbulent prandtl number"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: FROUDE_NUMBER_KEY              = "froude number"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: GRAVITY_DIRECTION_KEY          = "gravity direction"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: AOA_THETA_KEY                  = "aoa theta"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: AOA_PHI_KEY                    = "aoa phi"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: FLOW_EQUATIONS_KEY             = "flow equations"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: LESMODEL_KEY                   = "les model"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: COMPUTE_GRADIENTS_KEY          = "compute gradients"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: SUTHERLAND_TEMPERATURE_KEY     = "sutherland temperature"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: SUTHERLAND_REF_TEMPERATURE_KEY = "sutherland reference temperature"

         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(2) :: physics_NSKeywords = [MACH_NUMBER_KEY, FLOW_EQUATIONS_KEY]

         !PARTICLES
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: particlesKey             = "lagrangian particles"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfParticlesKey     = "number of particles"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: particlesPerParcelKey    = "particles per parcel"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: sourceTermKey            = "high order particles source term"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: STOKES_NUMBER_PART_KEY   = "stokes number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: GAMMA_PART_KEY           = "gamma"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PHI_M_PART_KEY           = "phi_m"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: I0_PART_KEY              = "radiation source"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: MIN_BOX_KEY              = "minimum box"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: MAX_BOX_KEY              = "maximum box"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: BC_BOX_KEY               = "bc box"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PART_FILE_KEY            = "particles file"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PART_LOG_FILE_KEY        = "vel and temp from file"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PART_LOG_INJ_KEY         = "injection"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PART_INJ_KEY             = "particles injection"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PART_NUMB_PER_STEP_KEY   = "particles per step"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PART_PERIOD_KEY          = "particles iter period"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: INJ_VEL_KEY              = "particles injection velocity"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: INJ_TEMP_KEY             = "particles injection temperature"

      END MODULE Physics_NSKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!
!    ******
     MODULE PhysicsStorage_NS
!    ******
!
     USE SMConstants
     use FluidData_NS
     use FileReadingUtilities, only: getRealArrayFromString

     IMPLICIT NONE

     private
     public    flowIsNavierStokes, NCONS, NGRAD
     public    IRHO, IRHOU, IRHOV, IRHOW, IRHOE
     public    NPRIM, IPIRHO, IPU, IPV, IPW, IPP, IPT, IPA2
     public    TemperatureReNormalization_Sutherland, S_div_TRef_Sutherland
     public    computeGradients

     public    ConstructPhysicsStorage_NS, DestructPhysicsStorage_NS, DescribePhysicsStorage_NS
     public    CheckPhysicsNSInputIntegrity
     public    GRADVARS_STATE, GRADVARS_ENTROPY, GRADVARS_ENERGY
     public    grad_vars, SetGradientVariables
!
!    ----------------------------
!    Either NavierStokes or Euler
!    ----------------------------
!
     logical, protected :: flowIsNavierStokes = .true.
     logical, protected :: computeGradients   = .true.
!
!    --------------------------
!!   The sizes of the NS system
!    --------------------------
!
     INTEGER, PARAMETER :: NCONS = 5, NGRAD = 5
!
!    -------------------------------------------
!!   The positions of the conservative variables
!    -------------------------------------------
!
     INTEGER, PARAMETER       :: IRHO = 1 , IRHOU = 2 , IRHOV = 3 , IRHOW = 4 , IRHOE = 5
!
!    ----------------------------------------
!!   The positions of the primitive variables
!    ----------------------------------------
!
     INTEGER, PARAMETER       :: NPRIM = 7
     INTEGER, PARAMETER       :: IPIRHO = 1, IPU = 2, IPV = 3, IPW = 4, IPP = 5, IPT = 6, IPA2 = 7
!
     real(kind=RP), protected :: TRef_Sutherland   ! Sutherland's reference temperature
     real(kind=RP), protected :: S_Sutherland      ! Sutherland's temperature
     real(kind=RP), protected :: TemperatureReNormalization_Sutherland ! TRef/SutherlandsTRef
     real(kind=RP), protected :: S_div_Tref_Sutherland
!
!    --------------------------------
!    Choice of the gradient variables
!    --------------------------------
!
     enum, bind(C)
        enumerator :: GRADVARS_STATE, GRADVARS_ENTROPY, GRADVARS_ENERGY
     end enum
     integer, protected :: grad_vars = GRADVARS_STATE

!
!    ------------------------------
!    Thermodynamic specs of the Air
!    ------------------------------
!
     type(Thermodynamics_t), target, private :: ThermodynamicsAir = Thermodynamics_t( &
                                                              "Air", & ! Name
                                                      287.15_RP, & ! R  ! J / (kg * °K)
                                                         1.4_RP, & ! gamma
                                                   sqrt(1.4_RP), & ! sqrtGamma
                                                1.4_RP - 1.0_RP, & ! gammaMinus1
                                     (1.4_RP - 1.0_RP) / 2.0_RP, & ! gammaMinus1Div2
                                     (1.4_RP + 1.0_RP) / 2.0_RP, & ! gammaPlus1Div2
                    (1.4_RP - 1.0_RP) / (2.0_RP * sqrt(1.4_RP)), & ! gammaMinus1Div2sg
                          (1.4_RP - 1.0_RP) / (2.0_RP * 1.4_RP), & ! gammaMinus1Div2g
                                     2.0_RP / (1.4_RP + 1.0_RP), & ! InvGammaPlus1Div2
                                     1.0_RP / (1.4_RP - 1.0_RP), & ! InvGammaMinus1
                                                1.0_RP / 1.4_RP, & ! InvGamma
                                   1.4_RP / ( 1.4_RP - 1.0_RP ), & ! gammaDivGammaMinus1
                       287.15_RP * 1.4_RP / ( 1.4_RP - 1.0_RP ), & ! cp
                                287.15_RP / ( 1.4_RP - 1.0_RP ), & ! cv
                                                        0.0_RP  & ! Bulk viscosity ratio
)
!
!    ========
     CONTAINS
!    ========
!
!     ///////////////////////////////////////////////////////
!
!     --------------------------------------------------
!!    Constructor: Define default values for the physics
!!    variables.
!     --------------------------------------------------
!
      SUBROUTINE ConstructPhysicsStorage_NS( controlVariables, Lref, timeref, success )
      USE FTValueDictionaryClass
      USE Physics_NSKeywordsModule
      use mainKeywordsModule, only: saveLambVectorKey
      use Utilities, only: toLower, almostEqual
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FTValueDictionary)      :: controlVariables
      real(kind=RP), intent(in)    :: Lref
      real(kind=RP), intent(out)   :: timeref
      LOGICAL                      :: success
!
!     ---------------
!     Local variables
!     ---------------
!
      CHARACTER(LEN=KEYWORD_LENGTH) :: keyword
      type(Thermodynamics_t), pointer  :: thermodynamics_
      type(RefValues_t)                :: refValues_
      type(Dimensionless_t)            :: dimensionless_
      real(kind=RP), allocatable       :: array(:)

!
!     --------------------
!     Collect input values
!     --------------------
!
      success = .TRUE.
      CALL CheckPhysicsNSInputIntegrity(controlVariables,success)
      IF(.NOT. success) RETURN
!
!
!     ---------------------
!     Set the gas to be air
!     ---------------------
!
      thermodynamics_ => thermodynamicsAir
!
!     ------------------------
!     Dimensionless quantities
!     ------------------------
!
      dimensionless_ % cp = thermodynamics_ % gamma * thermodynamics_ % InvGammaMinus1
      dimensionless_ % cv = thermodynamics_ % InvGammaMinus1
      dimensionless_ % Mach = controlVariables % doublePrecisionValueForKey(MACH_NUMBER_KEY)

      if ( controlVariables % ContainsKey(PRANDTL_NUMBER_KEY) ) then
         dimensionless_ % Pr   = controlVariables % doublePrecisionValueForKey(PRANDTL_NUMBER_KEY)
      else
         dimensionless_ % Pr = 0.72_RP
      end if
      if ( controlVariables % ContainsKey(TURBULENT_PRANDTL_NUMBER_KEY) ) then
         dimensionless_ % Prt   = controlVariables % doublePrecisionValueForKey(TURBULENT_PRANDTL_NUMBER_KEY)
      else
         dimensionless_ % Prt = dimensionless_ % Pr
      end if
!
!     *********************
!     Select flow equations
!     *********************
!
      keyword = controlVariables % stringValueForKey(FLOW_EQUATIONS_KEY,KEYWORD_LENGTH)
      CALL toLower(keyword)

      IF ( keyword == "euler" )     THEN
         flowIsNavierStokes = .FALSE.
         dimensionless_ % Re = 0.0_RP
         dimensionless_ % mu = 0.0_RP
         dimensionless_ % kappa = 0.0_RP

      ELSE
         flowIsNavierStokes = .TRUE.
!
!        ----------------------------
!        Look for the Reynolds number
!        ----------------------------
!
         IF ( controlVariables % containsKey(REYNOLDS_NUMBER_KEY) )     THEN
            dimensionless_ % Re = controlVariables % doublePrecisionValueForKey(REYNOLDS_NUMBER_KEY)

         ELSE
            PRINT *, "Input file is missing entry for keyword: ", REYNOLDS_NUMBER_KEY
            success = .FALSE.
            RETURN

         END IF
!
!        ------------------------------------------------
!        Set molecular viscosity and thermal conductivity
!        ------------------------------------------------
!
         if ( .not. almostEqual(dimensionless_ % Re, 0.0_RP) ) then
            dimensionless_ % mu   = 1.0_RP / dimensionless_ % Re
            dimensionless_ % kappa = 1.0_RP / ( thermodynamics_ % gammaMinus1 * &
                                                 POW2( dimensionless_ % Mach) * &
                                         dimensionless_ % Re * dimensionless_ % Pr )
         else
            dimensionless_ % mu = 0.0_RP
            dimensionless_ % kappa = 0.0_RP

         end if

         dimensionless_ % mu_to_kappa = 1.0_RP / (thermodynamics_ % gammaMinus1*POW2(dimensionless_ % Mach)* dimensionless_ % Pr)
      END IF
!
!     **************************************
!     Check if state gradients are requested
!     **************************************
!
      if ( controlVariables % containsKey(COMPUTE_GRADIENTS_KEY) ) then
!
!        Do not compute gradients if Euler equations are used and it is specified in the control file
!        --------------------------------------------------------------------------------------------
         if ( .not. flowIsNavierStokes ) then
            if ( .not. controlVariables % logicalValueForKey(COMPUTE_GRADIENTS_KEY) ) then
               computeGradients = .false.

            end if
         end if
      else
!
!        Do not compute gradients if Euler equations and the default option is selected
!        ------------------------------------------------------------------------------
         if ( .not. flowIsNavierStokes ) then
            computeGradients = .false.

         end if
      end if

      ! If Lamb vector is requested, computeGradient should be enabled
      if (controlVariables % logicalValueForKey(saveLambVectorKey)) then
         if (.not. computeGradients) then
            print *, trim(COMPUTE_GRADIENTS_KEY), " should be enabled to compute the Lamb vector."
            error stop
         end if
      end if

      dimensionless_ % gammaM2 = thermodynamics_ % gamma * POW2( dimensionless_ % Mach )
!
!     ********************
!     Set reference values
!     ********************
!
      if ( controlVariables % ContainsKey(REFERENCE_TEMPERATURE_KEY) ) then
         refValues_ % T = controlVariables % DoublePrecisionValueForKey(REFERENCE_TEMPERATURE_KEY)
      else
         refValues_ % T = 520.0_RP*5.0_RP/9.0_RP     ! ≃ 288.88 K
      end if

      if ( controlVariables % ContainsKey(REFERENCE_PRESSURE_KEY) ) then
         refValues_ % rho = controlVariables % DoublePrecisionValueForKey(REFERENCE_PRESSURE_KEY) / (thermodynamics_ % R * refValues_ % T)
      else
         refValues_ % rho = 101325.0_RP / (thermodynamics_ % R * refValues_ % T)
      end if

      refValues_ % V =   dimensionless_ % Mach &
                       * sqrt( thermodynamics_ % gamma * thermodynamics_ % R * refValues_ % T )

      refValues_ % p = refValues_ % rho * POW2( refValues_ % V )

      if ( flowIsNavierStokes ) then
         refValues_ % mu = refValues_ % rho * refValues_ % V * Lref * dimensionless_ % mu
         refValues_ % kappa = refValues_ % mu * thermodynamics_ % cp / dimensionless_ % Pr

      else
         refValues_ % mu = 0.0_RP
         refValues_ % kappa = 0.0_RP

      end if

      timeref = Lref / refValues_ % V
!
!     *******************************************
!     Set the Froude number and gravity direction
!     *******************************************
!
      if ( controlVariables % ContainsKey(FROUDE_NUMBER_KEY) ) then
         dimensionless_ % Fr = controlVariables % DoublePrecisionValueForKey(FROUDE_NUMBER_KEY)

      else
!
!        Default Froude number: earth's gravity
!        --------------------------------------
         dimensionless_ % Fr = refValues_ % V / sqrt(9.81_RP * Lref)

      end if

      if ( controlVariables % ContainsKey(GRAVITY_DIRECTION_KEY) ) then
         allocate(array(1:3))
         array = getRealArrayFromString( controlVariables % StringValueForKey(GRAVITY_DIRECTION_KEY,&
                                                                             KEYWORD_LENGTH))
         dimensionless_ % gravity_dir = array(1:3)

         if ( norm2(dimensionless_ % gravity_dir) < epsilon(1.0_RP)*10.0_RP ) then
!
!           Disable gravity
!           ---------------
            dimensionless_ % gravity_dir = 0.0_RP
            dimensionless_ % invFr2 = 0.0_RP
            dimensionless_ % Fr = huge(1.0_RP)
         else
            dimensionless_ % gravity_dir = dimensionless_ % gravity_dir / norm2(dimensionless_ % gravity_dir)
            dimensionless_ % invFr2 = 1.0_RP / POW2(dimensionless_ % Fr)

         end if
      else
         if ( controlVariables % ContainsKey(FROUDE_NUMBER_KEY) ) then
            print*, "When specifying a Froude number, the gravity direction must be specified"
            print*, "Gravity direction = [x,y,z]"
            errorMessage(STD_OUT)
            error stop

         else
!
!           Gravity is disabled
!           -------------------
            dimensionless_ % gravity_dir = 0.0_RP
            dimensionless_ % Fr = huge(1.0_RP)
            dimensionless_ % invFr2 = 0.0_RP

         end if
      end if
!
!     **********
!     LES Models
!     **********
!
      if ( controlVariables % containsKey(LESMODEL_KEY)) then
         if ( controlVariables % stringValueForKey(LESMODEL_KEY,4) .ne. "none") then
!
!           Enable NS fluxes
!           ----------------
            flowIsNavierStokes = .true.
            computeGradients = .true.

         end if
      end if
!
!     ***************
!     Angle of attack
!     ***************
!
      IF ( controlVariables % containsKey(AOA_PHI_KEY) )     THEN
         refValues_ % AOAPhi = controlVariables % doublePrecisionValueForKey(AOA_PHI_KEY)

      ELSE
!
!        Phi angle of attack is zero by default
!        --------------------------------------
         refValues_ % AOAPhi = 0.0_RP

      END IF

      IF ( controlVariables % containsKey(AOA_THETA_KEY) )     THEN
         refValues_ % AOATheta = controlVariables % doublePrecisionValueForKey(AOA_THETA_KEY)

      ELSE
!
!        Theta angle of attack is zero by default
!        ----------------------------------------
         refValues_ % AOATheta = 0.0_RP

      END IF
!
!     **************************
!     Sutherland's law constants
!     **************************
!
      if ( controlVariables % ContainsKey(SUTHERLAND_TEMPERATURE_KEY) ) then
         S_Sutherland = controlVariables % DoublePrecisionValueForKey(SUTHERLAND_TEMPERATURE_KEY)
      else
         S_Sutherland = 198.6_RP*5.0_RP/9.0_RP ! ≃ 110.4
      end if

      if ( controlVariables % ContainsKey(SUTHERLAND_REF_TEMPERATURE_KEY) ) then
         TRef_Sutherland = controlVariables % DoublePrecisionValueForKey(SUTHERLAND_REF_TEMPERATURE_KEY)
      else
         TRef_Sutherland = refValues_ % T
      end if

      S_div_TRef_Sutherland = S_Sutherland / TRef_Sutherland
      TemperatureReNormalization_Sutherland = refValues_ % T / TRef_Sutherland
!
!     **********************************************************************
!     Set the global (proteted) thermodynamics, dimensionless, and refValues
!     **********************************************************************
!
      call setThermodynamics( thermodynamics_ )
      call setDimensionless( dimensionless_ )
      call setRefValues( refValues_ )

      END SUBROUTINE ConstructPhysicsStorage_NS
!
!     ///////////////////////////////////////////////////////
!
!     -------------------------------------------------
!!    Destructor: Does nothing for this storage
!     -------------------------------------------------
!
      SUBROUTINE DestructPhysicsStorage_NS

      END SUBROUTINE DestructPhysicsStorage_NS
!
!     //////////////////////////////////////////////////////
!
!     -----------------------------------------
!!    Descriptor: Shows the gathered data
!     -----------------------------------------
!
      SUBROUTINE DescribePhysicsStorage_NS()
         USE Headers
         use MPI_Process_Info
         IMPLICIT NONE
         real(kind=RP)  :: pRef

         if ( .not. MPI_Process % isRoot ) return

         pRef = thermodynamics % R * refValues % rho * refValues % T

         write(STD_OUT,'(/,/)')
         if (flowIsNavierStokes) then
            call Section_Header("Loading Navier-Stokes physics")
         else
            call Section_Header("Loading Euler physics")
         end if

         write(STD_OUT,'(/)')
         call SubSection_Header("Fluid data")
         write(STD_OUT,'(30X,A,A22,A10)') "->" , "Gas: " , "Air"
         write(STD_OUT,'(30X,A,A22,F10.3,A)') "->" , "State constant: " , thermodynamics % R, " I.S."
         write(STD_OUT,'(30X,A,A22,F10.3)') "->" , "Specific heat ratio: " , thermodynamics % gamma

         write(STD_OUT,'(/)')
         call SubSection_Header("Reference quantities")
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference Temperature: " , refValues % T, " K."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference pressure: " , pRef, " Pa."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference density: " , refValues % rho , " kg/m^3."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference velocity: " , refValues % V , " m/s."

         if ( flowIsNavierStokes ) then
            write(STD_OUT,'(30X,A,A30,1pG10.3,A)') "->" , "Reference viscosity: ",refValues % mu , " Pa·s."
            write(STD_OUT,'(30X,A,A30,1pG10.3,A)') "->" , "Reference conductivity: ", refValues % kappa, " W/(m·K)."
         end if

         write(STD_OUT,'(/)')
         call SubSection_Header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A27,F10.3)') "->" , "Mach number: " , dimensionless % Mach
         if ( flowIsNavierStokes ) then
            write(STD_OUT,'(30X,A,A27,ES10.3)') "->" , "Reynolds number: " , dimensionless % Re
            write(STD_OUT,'(30X,A,A27,F10.3)') "->" , "Prandtl number: " , dimensionless % Pr
            write(STD_OUT,'(30X,A,A27,F10.3)') "->" , "Turbulent Prandtl number: " , dimensionless % Prt
         end if

         write(STD_OUT,'(30X,A,A27,F10.3)') "->" , "Froude number: " , dimensionless % Fr
         write(STD_OUT,'(30X,A,A27,A,F4.1,A,F4.1,A,F4.1,A)') "->" , "Gravity direction: ","[", &
                                                   dimensionless % gravity_dir(1), ", ", &
                                                   dimensionless % gravity_dir(2), ", ", &
                                                   dimensionless % gravity_dir(3), "]"

      END SUBROUTINE DescribePhysicsStorage_NS
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE CheckPhysicsNSInputIntegrity( controlVariables, success )
         USE FTValueDictionaryClass
         USE Physics_NSKeywordsModule
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(FTValueDictionary) :: controlVariables
         LOGICAL                 :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         CLASS(FTObject), POINTER :: obj
         INTEGER                  :: i
         success = .TRUE.

         DO i = 1, SIZE(physics_NSKeywords)
            obj => controlVariables % objectForKey(physics_NSKeywords(i))
            IF ( .NOT. ASSOCIATED(obj) )     THEN
               PRINT *, "Input file is missing entry for keyword: ",physics_NSKeywords(i)
               success = .FALSE.
            END IF
         END DO

      END SUBROUTINE CheckPhysicsNSInputIntegrity

      subroutine SetGradientVariables(grad_vars_)
         implicit none
         integer, intent(in)  :: grad_vars_

         select case(grad_vars_)
         case(GRADVARS_STATE, GRADVARS_ENTROPY, GRADVARS_ENERGY)
            grad_vars = grad_vars_
         case default
            print*, "Unrecognized option"
            errorMessage(STD_OUT)
            error stop
         end select

      end subroutine SetGradientVariables

!
!    **********
     END MODULE PhysicsStorage_NS
!    **********
