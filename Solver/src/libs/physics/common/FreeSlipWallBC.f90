#include "Includes.h"
module FreeSlipWallBCClass
   use SMConstants
   use PhysicsStorage
   use VariableConversion, only: GetGradientValues_f
   use FileReaders,            only: controlFileName
   use FileReadingUtilities,   only: GetKeyword, GetValueAsString, PreprocessInputLine
   use FTValueDictionaryClass, only: FTValueDictionary
   use GenericBoundaryConditionClass
   use FluidData
   use FileReadingUtilities, only: getRealArrayFromString
   use Utilities, only: toLower, almostEqual
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
   public FreeSlipWallBC_t
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
   type, extends(GenericBC_t) ::  FreeSlipWallBC_t
#ifdef NAVIERSTOKES
      logical           :: isAdiabatic
      real(kind=RP)     :: Twall
      real(kind=RP)     :: ewall       ! Wall internal energy
      real(kind=RP)     :: invTwall
      real(kind=RP)     :: wallType    ! 0/Adia 1/Isothermal
#endif
#ifdef CAHNHILLIARD
      real(kind=RP)     :: thetaw
#endif
      contains
         procedure         :: Destruct          => FreeSlipWallBC_Destruct
         procedure         :: Describe          => FreeSlipWallBC_Describe
#ifdef FLOW
         procedure         :: FlowState         => FreeSlipWallBC_FlowState
         procedure         :: FlowGradVars      => FreeSlipWallBC_FlowGradVars
         procedure         :: FlowNeumann       => FreeSlipWallBC_FlowNeumann
#endif
#ifdef CAHNHILLIARD
         procedure         :: PhaseFieldState   => FreeSlipWallBC_PhaseFieldState
         procedure         :: PhaseFieldNeumann => FreeSlipWallBC_PhaseFieldNeumann
         procedure         :: ChemPotState      => FreeSlipWallBC_ChemPotState
         procedure         :: ChemPotNeumann    => FreeSlipWallBC_ChemPotNeumann
#endif
   end type FreeSlipWallBC_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface FreeSlipWallBC_t
      module procedure ConstructFreeSlipWallBC
   end interface FreeSlipWallBC_t
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
      function ConstructFreeSlipWallBC(bname)
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
         type(FreeSlipWallBC_t)             :: ConstructFreeSlipWallBC
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

         ConstructFreeSlipWallBC % BCType = "freeslipwall"
         ConstructFreeSlipWallBC % bname  = bname
         call toLower(ConstructFreeSlipWallBC % bname)

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
               inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(ConstructFreeSlipWallBC % bname)) 
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
#ifdef NAVIERSTOKES
         if ( bcdict % ContainsKey("wall type (adiabatic/isothermal)") ) then
            keyval = bcdict % StringValueForKey("wall type (adiabatic/isothermal)", LINE_LENGTH)
            call tolower(keyval)
         else
            keyval = "adiabatic"
         end if

         if ( trim(keyval) .eq. "adiabatic" ) then
            ConstructFreeSlipWallBC % isAdiabatic = .true.
            ConstructFreeSlipWallBC % ewall = 0.0_RP
            ConstructFreeSlipWallBC % Twall = 0.0_RP
            ConstructFreeSlipWallBC % invTwall = 0.0_RP
            ConstructFreeSlipWallBC % wallType = 0.0_RP
         else
            ConstructFreeSlipWallBC % isAdiabatic = .false.
            call GetValueWithDefault(bcdict, "wall temperature" , refValues % T, ConstructFreeSlipWallBC % Twall     )
            ConstructFreeSlipWallBC % ewall = ConstructFreeSlipWallBC % Twall / (refValues % T*thermodynamics % gammaMinus1*dimensionless % gammaM2)
            ConstructFreeSlipWallBC % invTwall = dimensionless % gammaM2*refValues % T / ConstructFreeSlipWallBC % Twall
            ConstructFreeSlipWallBC % wallType = 1.0_RP
         end if
#ifdef TRANSPORT
         call constructScalarBoundaryConditions(bcdict, ConstructFreeSlipWallBC)
#endif
#endif
#ifdef CAHNHILLIARD
         call GetValueWithDefault(bcdict, "contact angle", 90.0_RP, ConstructFreeSlipWallBC % thetaw)
#endif

         close(fid)
         call bcdict % Destruct
   
      end function ConstructFreeSlipWallBC

      subroutine FreeSlipWallBC_Describe(self)
!
!        ***************************************************
!              Describe the inflow boundary condition
         implicit none
         class(FreeSlipWallBC_t),  intent(in)  :: self
         write(STD_OUT,'(30X,A,A28,A)') "->", " Boundary condition type: ", "FreeSlipWall"
#ifdef NAVIERSTOKES
         if ( self % isAdiabatic ) then
            write(STD_OUT,'(30X,A,A28,A)') "->", ' Thermal type: ', "Adiabatic"


         else
            write(STD_OUT,'(30X,A,A28,A)') "->", ' Thermal type: ', "Isothermal"
            write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Wall temperature: ', self % Twall * refValues % T
         end if
#if defined(TRANSPORT)
         call self % ScalarBC % Describe()
#endif
#endif
#ifdef CAHNHILLIARD
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Wall contact angle: ', self % thetaw
#endif
         
      end subroutine FreeSlipWallBC_Describe

!
!/////////////////////////////////////////////////////////
!
!        Class destructors
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      subroutine FreeSlipWallBC_Destruct(self)
         implicit none
         class(FreeSlipWallBC_t)    :: self
#ifdef TRANSPORT
         call self % ScalarBC % Destruct()
#endif
      end subroutine FreeSlipWallBC_Destruct
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for compressible Navier--Stokes equations
!        -----------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#ifdef NAVIERSTOKES
      subroutine FreeSlipWallBC_FlowState(self, x, t, nHat, Q)
!
!        *************************************************************
!           Compute the state variables for a general wall
!
!           · Density is computed from the interior state
!           · Wall velocity is set to v_interior - 2v_normal:
!                 -> Maintains tangential speed.
!                 -> Removes normal speed (weakly).
!           · Internal energy is either the interior state for 
!              adiabatic walls, or the imposed for isothermal.
!              eBC = eInt + kWallType (eIso - eInt)
!              where kWallType = 0 for adiabatic and 1 for isothermal.        
!        *************************************************************
!
         implicit none
         class(FreeSlipWallBC_t), intent(in)    :: self
         real(kind=RP),           intent(in)    :: x(NDIM)
         real(kind=RP),           intent(in)    :: t
         real(kind=RP),           intent(in)    :: nHat(NDIM)
         real(kind=RP),           intent(inout) :: Q(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: qNorm, pressure_aux

         qNorm = nHat(IX) * Q(IRHOU) + nHat(IY) * Q(IRHOV) + nHat(IZ) * Q(IRHOW)

         Q(IRHOU:IRHOW) = Q(IRHOU:IRHOW) - 2.0_RP * qNorm * nHat

         !Isothermal BC
         pressure_aux = Q(IRHO) * self % Twall / (refValues % T * dimensionless % gammaM2)
         Q(IRHOE) = Q(IRHOE) + self % wallType*(pressure_aux/thermodynamics % gammaMinus1 + 0.5_RP*(POW2(Q(IRHOU))+POW2(Q(IRHOV))+POW2(Q(IRHOW)))/Q(IRHO) - Q(IRHOE))
#if defined(TRANSPORT)
         call self % ScalarBC % FlowState(x, t, nHat, Q)
#endif

      end subroutine FreeSlipWallBC_FlowState

      subroutine FreeSlipWallBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
!
!        *****************************************************************
!           Only set the temperature, velocity is Neumann, use interior!
!        *****************************************************************
!
         implicit none
         class(FreeSlipWallBC_t),  intent(in)    :: self
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
         real(kind=RP)  :: rhou_n
         real(kind=RP)  :: Q_aux(NCONS)

         Q_aux(IRHO) = Q(IRHO)
         Q_aux(IRHOU:IRHOW) = Q(IRHOU:IRHOW)
         Q_aux(IRHOE) = Q(IRHOE) + self % wallType*(Q(IRHO)*self % eWall+0.5_RP*(POW2(Q(IRHOU))+POW2(Q(IRHOV))+POW2(Q(IRHOW)))/Q(IRHO)-Q(IRHOE))
#if defined(SPALARTALMARAS)
         Q_aux(IRHOTHETA)= Q(IRHOTHETA)
#endif
         call GetGradients(NCONS, NGRAD, Q_aux, U)
#if defined(TRANSPORT)
         call self % ScalarBC % FlowGradVars(x, t, nHat, Q, U, GetGradients)
#endif
      end subroutine FreeSlipWallBC_FlowGradVars

      subroutine FreeSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
!
!        ***********************************************************
!           In momentum, free slip is Neumann. In temperature, 
!           depends on the adiabatic/isothermal choice
!        ***********************************************************
!
         implicit none
         class(FreeSlipWallBC_t),   intent(in) :: self
         real(kind=RP),       intent(in)       :: x(NDIM)
         real(kind=RP),       intent(in)       :: t
         real(kind=RP),       intent(in)       :: nHat(NDIM)
         real(kind=RP),       intent(in)       :: Q(NCONS)
         real(kind=RP),       intent(in)       :: U_x(NGRAD)
         real(kind=RP),       intent(in)       :: U_y(NGRAD)
         real(kind=RP),       intent(in)       :: U_z(NGRAD)
         real(kind=RP),       intent(inout)    :: flux(NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!   
         real(kind=RP)  :: viscWork, heatFlux

         viscWork = (flux(IRHOU)*Q(IRHOU)+flux(IRHOV)*Q(IRHOV)+flux(IRHOW)*Q(IRHOW))/Q(IRHO)
         heatFlux = flux(IRHOE) - viscWork
         flux(IRHO:IRHOW) = 0.0_RP
         flux(IRHOE) = self % wallType * heatFlux  ! 0 (Adiabatic)/ heatFlux (Isothermal)
#if defined(SPALARTALMARAS)
         flux(IRHOTHETA) = 0.0_RP
#endif
#if defined(TRANSPORT)
         call self % ScalarBC % FlowNeumann(x, t, nHat, Q, U_x, U_y, U_z, flux)
#endif
      end subroutine FreeSlipWallBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for incompressible Navier--Stokes equations
!        -------------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#ifdef INCNS
      subroutine FreeSlipWallBC_FlowState(self, x, t, nHat, Q)
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
         class(FreeSlipWallBC_t),  intent(in)    :: self
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
         vn = sum(Q(INSRHOU:INSRHOW)*nHat)

         Q(INSRHO)          = Q(INSRHO)
         Q(INSRHOU:INSRHOW) = Q(INSRHOU:INSRHOW) - 2.0_RP * vn * nHat
         Q(INSP)            = Q(INSP)

      end subroutine FreeSlipWallBC_FlowState

      subroutine FreeSlipWallBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
!
!        **************************************************************
!           Use the interior velocity: Neumann BC!
!        **************************************************************
!
         implicit none
         class(FreeSlipWallBC_t),  intent(in)  :: self
         real(kind=RP),          intent(in)    :: x(NDIM)
         real(kind=RP),          intent(in)    :: t
         real(kind=RP),          intent(in)    :: nHat(NDIM)
         real(kind=RP),          intent(in)    :: Q(NCONS)
         real(kind=RP),          intent(inout) :: U(NGRAD)
         procedure(GetGradientValues_f)        :: GetGradients

      end subroutine FreeSlipWallBC_FlowGradVars

      subroutine FreeSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
!
!        ***************************************************************
!           Set homogeneous Neumann BCs everywhere
!        ***************************************************************
!        
         implicit none
         class(FreeSlipWallBC_t),  intent(in) :: self
         real(kind=RP),       intent(in)      :: x(NDIM)
         real(kind=RP),       intent(in)      :: t
         real(kind=RP),       intent(in)      :: nHat(NDIM)
         real(kind=RP),       intent(in)      :: Q(NCONS)
         real(kind=RP),       intent(in)      :: U_x(NCONS)
         real(kind=RP),       intent(in)      :: U_y(NCONS)
         real(kind=RP),       intent(in)      :: U_z(NCONS)
         real(kind=RP),       intent(inout)   :: flux(NCONS)

         flux = 0.0_RP

      end subroutine FreeSlipWallBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for multiphase solver
!        ---------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#ifdef MULTIPHASE
      subroutine FreeSlipWallBC_FlowState(self, x, t, nHat, Q)
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
         class(FreeSlipWallBC_t),  intent(in)    :: self
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
         vn = sum(Q(IMSQRHOU:IMSQRHOW)*nHat)

         Q(IMC)      = Q(IMC)
         Q(IMSQRHOU) = Q(IMSQRHOU) - 2.0_RP * vn * nHat(IX)
         Q(IMSQRHOV) = Q(IMSQRHOV) - 2.0_RP * vn * nHat(IY)
         Q(IMSQRHOW) = Q(IMSQRHOW) - 2.0_RP * vn * nHat(IZ)
         Q(IMP)      = Q(IMP)

      end subroutine FreeSlipWallBC_FlowState

      subroutine FreeSlipWallBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
!
!        **************************************************************
!           Use the interior velocity: Neumann BC!
!        **************************************************************
!
         implicit none
         class(FreeSlipWallBC_t),  intent(in)  :: self
         real(kind=RP),          intent(in)    :: x(NDIM)
         real(kind=RP),          intent(in)    :: t
         real(kind=RP),          intent(in)    :: nHat(NDIM)
         real(kind=RP),          intent(in)    :: Q(NCONS)
         real(kind=RP),          intent(inout) :: U(NGRAD)
         procedure(GetGradientValues_f)        :: GetGradients
      end subroutine FreeSlipWallBC_FlowGradVars

      subroutine FreeSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(FreeSlipWallBC_t),  intent(in) :: self
         real(kind=RP),       intent(in)      :: x(NDIM)
         real(kind=RP),       intent(in)      :: t
         real(kind=RP),       intent(in)      :: nHat(NDIM)
         real(kind=RP),       intent(in)      :: Q(NCONS)
         real(kind=RP),       intent(in)      :: U_x(NCONS)
         real(kind=RP),       intent(in)      :: U_y(NCONS)
         real(kind=RP),       intent(in)      :: U_z(NCONS)
         real(kind=RP),       intent(inout)   :: flux(NCONS)

         flux = 0.0_RP

      end subroutine FreeSlipWallBC_FlowNeumann
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
      subroutine FreeSlipWallBC_PhaseFieldState(self, x, t, nHat, Q)
         implicit none
         class(FreeSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine FreeSlipWallBC_PhaseFieldState

      subroutine FreeSlipWallBC_PhaseFieldNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(FreeSlipWallBC_t),  intent(in) :: self
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

      end subroutine FreeSlipWallBC_PhaseFieldNeumann

      subroutine FreeSlipWallBC_ChemPotState(self, x, t, nHat, Q)
         implicit none
         class(FreeSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine FreeSlipWallBC_ChemPotState

      subroutine FreeSlipWallBC_ChemPotNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(FreeSlipWallBC_t),  intent(in) :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCOMP)
         real(kind=RP),       intent(in)    :: U_x(NCOMP)
         real(kind=RP),       intent(in)    :: U_y(NCOMP)
         real(kind=RP),       intent(in)    :: U_z(NCOMP)
         real(kind=RP),       intent(inout) :: flux(NCOMP)

         flux = 0.0_RP

      end subroutine FreeSlipWallBC_ChemPotNeumann

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
!        *************************************************************
!
      subroutine FreeSlipWallBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(FreeSlipWallBC_t),  intent(in)  :: self
         real(kind=RP),       intent(in)       :: x(NDIM)
         real(kind=RP),       intent(in)       :: t
         real(kind=RP),       intent(in)       :: nHat(NDIM)
         real(kind=RP),       intent(inout)    :: Q(NCONS)
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

      end subroutine FreeSlipWallBC_FlowState
!
!     not use for acoustic
      subroutine FreeSlipWallBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
!
         implicit none
         class(FreeSlipWallBC_t),  intent(in)  :: self
         real(kind=RP),          intent(in)    :: x(NDIM)
         real(kind=RP),          intent(in)    :: t
         real(kind=RP),          intent(in)    :: nHat(NDIM)
         real(kind=RP),          intent(in)    :: Q(NCONS)
         real(kind=RP),          intent(inout) :: U(NGRAD)
         procedure(GetGradientValues_f)        :: GetGradients
      end subroutine FreeSlipWallBC_FlowGradVars
!
      subroutine FreeSlipWallBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(FreeSlipWallBC_t),  intent(in)    :: self
         real(kind=RP),            intent(in)    :: x(NDIM)
         real(kind=RP),            intent(in)    :: t
         real(kind=RP),            intent(in)    :: nHat(NDIM)
         real(kind=RP),            intent(in)    :: Q(NCONS)
         real(kind=RP),            intent(in)    :: U_x(NCONS)
         real(kind=RP),            intent(in)    :: U_y(NCONS)
         real(kind=RP),            intent(in)    :: U_z(NCONS)
         real(kind=RP),            intent(inout) :: flux(NCONS)

         flux = 0.0_RP

      end subroutine FreeSlipWallBC_FlowNeumann
#endif
!
end module FreeSlipWallBCClass
