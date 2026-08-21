#include "Includes.h"
module NoSlipWallBCClass
   use SMConstants
   use PhysicsStorage
   use FileReaders,            only: controlFileName
   use FileReadingUtilities,   only: GetKeyword, GetValueAsString ,PreprocessInputLine
   use FTValueDictionaryClass, only: FTValueDictionary
   use GenericBoundaryConditionClass
   use FluidData
   use FileReadingUtilities, only: getRealArrayFromString
   use Utilities, only: toLower, almostEqual
   use VariableConversion, only: GetGradientValues_f
   implicit none
!
!  *****************************
!  Default everything to private
!  *****************************
!
   private
!
!  ****************
!  Public variables
!  ****************
!
!
!  ******************
!  Public definitions
!  ******************
!
   public NoSlipWallBC_t
!
!  ****************************
!  Static variables definitions
!  ****************************
!
!
!  ****************
!  Class definition
!  ****************
!
   type, extends(GenericBC_t) ::  NoSlipWallBC_t
#if defined(NAVIERSTOKES)
      logical           :: isAdiabatic
      real(kind=RP)     :: Twall
      real(kind=RP)     :: invTwall
      real(kind=RP)     :: ewall       ! Wall internal energy
      real(kind=RP)     :: wallType
#endif
#ifdef FLOW
      real(kind=RP)     :: vWall(NDIM)
#endif
#ifdef CAHNHILLIARD
      real(kind=RP)     :: thetaw
#endif
      contains
         procedure         :: Destruct          => NoSlipWallBC_Destruct
         procedure         :: Describe          => NoSlipWallBC_Describe
#ifdef FLOW
         procedure         :: FlowState         => NoSlipWallBC_FlowState
         procedure         :: FlowGradVars      => NoSlipWallBC_FlowGradVars
         procedure         :: FlowNeumann       => NoSlipWallBC_FlowNeumann
#endif
#ifdef CAHNHILLIARD
         procedure         :: PhaseFieldState   => NoSlipWallBC_PhaseFieldState
         procedure         :: PhaseFieldNeumann => NoSlipWallBC_PhaseFieldNeumann
         procedure         :: ChemPotState      => NoSlipWallBC_ChemPotState
         procedure         :: ChemPotNeumann    => NoSlipWallBC_ChemPotNeumann
#endif
   end type NoSlipWallBC_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface NoSlipWallBC_t
      module procedure ConstructNoSlipWallBC
   end interface NoSlipWallBC_t
!
!  *******************
!  Function prototypes
!  *******************
!
   abstract interface
   end interface
!
!  ========
   contains
!  ========
!
!/////////////////////////////////////////////////////////
!
!        Class constructor
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      function ConstructNoSlipWallBC(bname)
!
!        ********************************************************************
!        · Definition of the noSlipWall boundary condition in the control file:
!              #define boundary bname
!                 type             = noSlipWall
!                 velocity         = #value        (only in incompressible NS)
!                 Mach number      = #value        (only in compressible NS)
!                 AoAPhi           = #value
!                 AoATheta         = #value
!                 density          = #value        (only in monophase)
!                 pressure         = #value        (only in compressible NS)
!                 multiphase type  = mixed/layered
!                 phase 1 layer x  > #value
!                 phase 1 layer y  > #value
!                 phase 1 layer z  > #value
!                 phase 1 velocity = #value
!                 phase 2 velocity = #value
!              #end
!        ********************************************************************
!
#ifdef TRANSPORT
         use ScalarBoundaryConditions, only: constructScalarBoundaryConditions
#endif
         implicit none
         type(NoSlipWallBC_t)             :: ConstructNoSlipWallBC
         character(len=*), intent(in) :: bname
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: fid, io
         character(len=LINE_LENGTH) :: boundaryHeader
         character(len=LINE_LENGTH) :: currentLine
         character(len=LINE_LENGTH) :: keyword, keyval
         logical                    :: inside
         type(FTValueDIctionary)    :: bcdict

         open(newunit = fid, file = trim(controlFileName), status = "old", action = "read")

         call bcdict % InitWithSize(16)

         ConstructNoSlipWallBC % BCType = "noslipwall"
         ConstructNoSlipWallBC % bname  = bname
         call toLower(ConstructNoSlipWallBC % bname)

         write(boundaryHeader,'(A,A)') "#define boundary ",trim(bname)
         call toLower(boundaryHeader)
!
!        Navigate until the "#define boundary bname" sentinel is found
!        -------------------------------------------------------------
         inside = .false.
         do 
            read(fid, '(A)', iostat=io) currentLine

            IF(io .ne. 0 ) EXIT

            call PreprocessInputLine(currentLine)
            call toLower(currentLine)

            if ( index(trim(currentLine),"#define boundary") .ne. 0 ) then
               inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(ConstructNoSlipWallBC % bname)) 
            end if
!
!           Get all keywords inside the zone
!           --------------------------------
            if ( inside ) then
               if ( trim(currentLine) .eq. "#end" ) exit

               keyword  = ADJUSTL(GetKeyword(currentLine))
               keyval   = ADJUSTL(GetValueAsString(currentLine))
               call ToLower(keyword)
      
               call bcdict % AddValueForKey(keyval, trim(keyword))

            end if

         end do
!
!        Analyze the gathered data
!        -------------------------
#if defined(NAVIERSTOKES)
         if ( bcdict % ContainsKey("wall type (adiabatic/isothermal)") ) then
            keyval = bcdict % StringValueForKey("wall type (adiabatic/isothermal)", LINE_LENGTH)
            call tolower(keyval)
         else
            keyval = "adiabatic"
         end if

         if ( trim(keyval) .eq. "adiabatic" ) then
            ConstructNoSlipWallBC % isAdiabatic = .true.
            ConstructNoSlipWallBC % ewall = 0.0_RP
            ConstructNoSlipWallBC % Twall = 0.0_RP
            ConstructNoSlipWallBC % invTwall = 0.0_RP
            ConstructNoSlipWallBC % wallType = 0.0_RP
         
         else
            ConstructNoSlipWallBC % isAdiabatic = .false.
            call GetValueWithDefault(bcdict, "wall temperature" , refValues % T, ConstructNoSlipWallBC % Twall)
            ConstructNoSlipWallBC % ewall = ConstructNoSlipWallBC % Twall / (refValues % T*thermodynamics % gammaMinus1*dimensionless % gammaM2)
            ConstructNoSlipWallBC % invTwall = dimensionless % gammaM2*refValues % T / ConstructNoSlipWallBC % Twall
            ConstructNoSlipWallBC % wallType = 1.0_RP
         end if
#ifdef TRANSPORT
         call constructScalarBoundaryConditions(bcdict, ConstructNoSlipWallBC)
#endif
#endif
#ifdef FLOW
         if ( bcdict % ContainsKey("wall velocity") ) then
            ConstructNoSlipWallBC % vWall = getRealArrayFromString( bcdict % StringValueForKey("wall velocity",&
                                                                                           LINE_LENGTH))    
            
         else
            ConstructNoSlipWallBC % vWall = 0.0_RP
         end if
#endif
#ifdef CAHNHILLIARD
         call GetValueWithDefault(bcdict, "contact angle", 90.0_RP, ConstructNoSlipWallBC % thetaw)
#endif

         close(fid)
         call bcdict % Destruct
   
      end function ConstructNoSlipWallBC

      subroutine NoSlipWallBC_Describe(self)
!
!        ***************************************************
!              Describe the inflow boundary condition
         implicit none
         class(NoSlipWallBC_t),  intent(in)  :: self
         write(STD_OUT,'(30X,A,A28,A)') "->", " Boundary condition type: ", "NoSlipWall"
#ifdef FLOW
         write(STD_OUT,'(30X,A,A28,A,F10.2,A,F10.2,A,F10.2,A)') "->", ' Wall velocity: ',"[",self % vWall(1),",",self % vWall(2),",",self % vWall(3),"]"
#endif
#ifdef NAVIERSTOKES
         if ( self % isAdiabatic ) then
            write(STD_OUT,'(30X,A,A28,A)') "->", ' Thermal type: ', "Adiabatic"


         else
            write(STD_OUT,'(30X,A,A28,A)') "->", ' Thermal type: ', "Isothermal"
            write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Wall temperature: ', self % Twall
         end if
#if defined(TRANSPORT)
         call self % ScalarBC % Describe()
#endif
#endif
#ifdef CAHNHILLIARD
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Wall contact angle: ', self % thetaw
#endif
         
      end subroutine NoSlipWallBC_Describe

!
!/////////////////////////////////////////////////////////
!
!        Class destructors
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      subroutine NoSlipWallBC_Destruct(self)
         implicit none
         class(NoSlipWallBC_t)    :: self
#ifdef TRANSPORT
         call self % ScalarBC % Destruct()
#endif
      end subroutine NoSlipWallBC_Destruct
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for compressible Navier--Stokes equations
!        -----------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES)
      subroutine NoSlipWallBC_FlowState(self, x, t, nHat, Q)
!
!        *************************************************************
!           Compute the state variables for a general wall
!
!           · SHOULD BE: Cancel out normal velocity
!           · It cancels the whole velocity because otherwise the IP won't work
!              I need to update the IP (and the analytical Jacobian) to this new
!              whole BC approach.
!        *************************************************************
!
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),          intent(in)    :: x(NDIM)
         real(kind=RP),          intent(in)    :: t
         real(kind=RP),          intent(in)    :: nHat(NDIM)
         real(kind=RP),          intent(inout) :: Q(NCONS)

#if defined (SPALARTALMARAS)
         Q(IRHOTHETA)   = -Q(IRHOTHETA)
#endif
         Q(IRHOU:IRHOW) = 2.0_RP * Q(IRHO)*self % vWall - Q(IRHOU:IRHOW)
!        This boundary condition should be
!        ---------------------------------
         !Q(IRHOU:IRHOW) = Q(IRHOU:IRHOW) - 2.0_RP * sum(Q(IRHOU:IRHOW)*nHat)*nHat


         !Isothermal BC
         Q(IRHOE) = Q(IRHOE) + self % wallType * (Q(IRHO) * self % Twall / (refValues % T * dimensionless % gammaM2 * thermodynamics % gammaMinus1) - Q(IRHOE))
#if defined(TRANSPORT)
         call self % ScalarBC % FlowState(x, t, nHat, Q)
#endif

      end subroutine NoSlipWallBC_FlowState

      subroutine NoSlipWallBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
!
!        **************************************************************
!              Computes the set of gradient variables U* at the wall
!        **************************************************************
!
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),          intent(in)    :: x(NDIM)
         real(kind=RP),          intent(in)    :: t
         real(kind=RP),          intent(in)    :: nHat(NDIM)
         real(kind=RP),          intent(in)    :: Q(NCONS)
         real(kind=RP),          intent(inout) :: U(NGRAD)
         procedure(GetGradientValues_f)        :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: Q_aux(NCONS), U1
         real(kind=RP)  :: e_int
         real(kind=RP)  :: invRho

         invRho = 1.0_RP / Q(IRHO)
         e_int = invRho*(Q(IRHOE) - 0.5_RP*invRho*(POW2(Q(IRHOU))+POW2(Q(IRHOV))+POW2(Q(IRHOW))))
      
         Q_aux(IRHO) = Q(IRHO)
         Q_aux(IRHOU:IRHOW) = Q(IRHO)*self % vWall
         Q_aux(IRHOE) = Q(IRHO)*((1.0_RP-self % wallType)*e_int + self % wallType*self % eWall + 0.5_RP*sum(self % vWall*self % vWall))
#if defined (SPALARTALMARAS)
         Q_aux(IRHOTHETA) = 0.0_RP
#endif

         U1 = U(IRHO)

         call GetGradients(NCONS, NGRAD, Q_aux, U)

         U(IRHO) = U1
#if defined(TRANSPORT)
         call self % ScalarBC % FlowGradVars(x, t, nHat, Q, U, GetGradients)
#endif
      end subroutine NoSlipWallBC_FlowGradVars

      subroutine NoSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
!
!        ***********************************************************
!           Cancel out the temperature flux for adiabatic BCs
!        ***********************************************************
!
         implicit none
         class(NoSlipWallBC_t),   intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCONS)
         real(kind=RP),       intent(in)    :: U_x(NGRAD)
         real(kind=RP),       intent(in)    :: U_y(NGRAD)
         real(kind=RP),       intent(in)    :: U_z(NGRAD)
         real(kind=RP),       intent(inout) :: flux(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: viscWork, heatFlux, invRho, u, v, w

         invRho = 1.0_RP / Q(IRHO)
         u      = invRho * Q(IRHOU)
         v      = invRho * Q(IRHOV)
         w      = invRho * Q(IRHOW)
         viscWork = u*flux(IRHOU) + v*flux(IRHOV) + w*flux(IRHOW)
         heatFlux = flux(IRHOE) - viscWork

         flux(IRHO)  = 0.0_RP
         flux(IRHOE) = sum(self % vWall*flux(IRHOU:IRHOW)) + self % wallType * heatFlux  ! 0 (Adiabatic)/ heatFlux (Isothermal)
#if defined(TRANSPORT)
         call self % ScalarBC % FlowNeumann(x, t, nHat, Q, U_x, U_y, U_z, flux)
#endif
      end subroutine NoSlipWallBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for incompressible Navier--Stokes equations
!        -------------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
!
#if defined(INCNS)
      subroutine NoSlipWallBC_FlowState(self, x, t, nHat, Q)
!
!        *************************************************************
!           Compute the state variables for a general wall
!
!           · Density is computed from the interior state
!           · Wall velocity is set to 2v_wall - v_interior
!           · Pressure is computed from the interior state
!        *************************************************************
!

         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)

         Q(INSRHOU:INSRHOW) = -Q(INSRHOU:INSRHOW)

      end subroutine NoSlipWallBC_FlowState

      subroutine NoSlipWallBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),          intent(in)    :: x(NDIM)
         real(kind=RP),          intent(in)    :: t
         real(kind=RP),          intent(in)    :: nHat(NDIM)
         real(kind=RP),          intent(in)    :: Q(NCONS)
         real(kind=RP),          intent(inout) :: U(NGRAD)
         procedure(GetGradientValues_f)        :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         U(INSRHOU:INSRHOW) = self % vWall
   
      end subroutine NoSlipWallBC_FlowGradVars

      subroutine NoSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
!
!        **********************************************
!           Do nothing: no slip walls are Dirichlet BCs
!        **********************************************
!
         implicit none
         class(NoSlipWallBC_t),  intent(in) :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCONS)
         real(kind=RP),       intent(in)    :: U_x(NCONS)
         real(kind=RP),       intent(in)    :: U_y(NCONS)
         real(kind=RP),       intent(in)    :: U_z(NCONS)
         real(kind=RP),       intent(inout) :: flux(NCONS)
      end subroutine NoSlipWallBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for the multiphase solver
!        -------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#ifdef MULTIPHASE
      subroutine NoSlipWallBC_FlowState(self, x, t, nHat, Q)
!
!        *************************************************************
!           Compute the state variables for a general wall
!
!           · Density is computed from the interior state
!           · Wall velocity is set to 2v_wall - v_interior
!           · Pressure is computed from the interior state
!        *************************************************************
!

         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)

         Q(IMSQRHOU:IMSQRHOW) = -Q(IMSQRHOU:IMSQRHOW)

      end subroutine NoSlipWallBC_FlowState

      subroutine NoSlipWallBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),          intent(in)    :: x(NDIM)
         real(kind=RP),          intent(in)    :: t
         real(kind=RP),          intent(in)    :: nHat(NDIM)
         real(kind=RP),          intent(in)    :: Q(NCONS)
         real(kind=RP),          intent(inout) :: U(NGRAD)
         procedure(GetGradientValues_f)        :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         U(IMSQRHOU:IMSQRHOW) = self % vWall
   
      end subroutine NoSlipWallBC_FlowGradVars

      subroutine NoSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
!
!        ************************************************************************
!           No slip wall is dirichlet on momentum, Neumann on chemical potential
!           -> Cancel out the chemical potential gradient
!        ************************************************************************
         implicit none
         class(NoSlipWallBC_t),  intent(in) :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCONS)
         real(kind=RP),       intent(in)    :: U_x(NCONS)
         real(kind=RP),       intent(in)    :: U_y(NCONS)
         real(kind=RP),       intent(in)    :: U_z(NCONS)
         real(kind=RP),       intent(inout) :: flux(NCONS)

         flux(IMC) = 0.0_RP

      end subroutine NoSlipWallBC_FlowNeumann
#endif

!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for Cahn--Hilliard
!        ------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(CAHNHILLIARD)
      subroutine NoSlipWallBC_PhaseFieldState(self, x, t, nHat, Q)
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine NoSlipWallBC_PhaseFieldState

      subroutine NoSlipWallBC_PhaseFieldNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(NoSlipWallBC_t),  intent(in) :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCOMP)
         real(kind=RP),       intent(in)    :: U_x(NCOMP)
         real(kind=RP),       intent(in)    :: U_y(NCOMP)
         real(kind=RP),       intent(in)    :: U_z(NCOMP)
         real(kind=RP),       intent(inout) :: flux(NCOMP)
!
!        ---------------
!        Local variables
!        ---------------
!  
         real(kind=RP), parameter   :: MIN_ = 1.0e-1_RP
         real(kind=RP)              :: prod
         real(kind=RP)              :: prod12, prod13, prod23, c3

         prod = Q(1) * (1.0_RP - Q(1))

         if ( prod .le. MIN_ ) then
            prod = 0.0_RP
         end if

         flux = -4.0_RP * multiphase % invEps * cos(DEG2RAD*self % thetaw) * prod 

      end subroutine NoSlipWallBC_PhaseFieldNeumann

      subroutine NoSlipWallBC_ChemPotState(self, x, t, nHat, Q)
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine NoSlipWallBC_ChemPotState

      subroutine NoSlipWallBC_ChemPotNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(NoSlipWallBC_t),  intent(in) :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCOMP)
         real(kind=RP),       intent(in)    :: U_x(NCOMP)
         real(kind=RP),       intent(in)    :: U_y(NCOMP)
         real(kind=RP),       intent(in)    :: U_z(NCOMP)
         real(kind=RP),       intent(inout) :: flux(NCOMP)

         flux = 0.0_RP

      end subroutine NoSlipWallBC_ChemPotNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for Acoustic APE equations
!        ---------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(ACOUSTIC)
!        *************************************************************
!           Compute the state variables for a general wall
!
!           · Density is computed from the interior state
!           · Wall velocity is set to 2v_wall - v_interior
!           · Pressure is computed from the interior state
!           Exacat same as for FreeSlipWall, no difference in Acoustic
!        *************************************************************
!
      subroutine NoSlipWallBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(NoSlipWallBC_t), intent(in)  :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: vn
!
!        -----------------------------------------------
!        Generate the external flow along the face, that
!        represents a solid wall.
!        -----------------------------------------------
!
         vn = sum(Q(ICAAU:ICAAW)*nHat)

         Q(ICAARHO) = Q(ICAARHO)
         Q(ICAAU:ICAAW) = Q(ICAAU:ICAAW) - 2.0_RP * vn * nHat
         Q(ICAAP) = Q(ICAAP)

      end subroutine NoSlipWallBC_FlowState
!
!     not use for acoustic
      subroutine NoSlipWallBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
!
         implicit none
         class(NoSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),          intent(in)    :: x(NDIM)
         real(kind=RP),          intent(in)    :: t
         real(kind=RP),          intent(in)    :: nHat(NDIM)
         real(kind=RP),          intent(in)    :: Q(NCONS)
         real(kind=RP),          intent(inout) :: U(NGRAD)
         procedure(GetGradientValues_f)        :: GetGradients
      end subroutine NoSlipWallBC_FlowGradVars
!
      subroutine NoSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(NoSlipWallBC_t), intent(in)    :: self
         real(kind=RP),         intent(in)    :: x(NDIM)
         real(kind=RP),         intent(in)    :: t
         real(kind=RP),         intent(in)    :: nHat(NDIM)
         real(kind=RP),         intent(in)    :: Q(NCONS)
         real(kind=RP),         intent(in)    :: U_x(NCONS)
         real(kind=RP),         intent(in)    :: U_y(NCONS)
         real(kind=RP),         intent(in)    :: U_z(NCONS)
         real(kind=RP),         intent(inout) :: flux(NCONS)

         flux = 0.0_RP

      end subroutine NoSlipWallBC_FlowNeumann
#endif
!
end module NoSlipWallBCClass
