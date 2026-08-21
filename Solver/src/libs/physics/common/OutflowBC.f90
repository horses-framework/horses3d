#include "Includes.h"
module OutflowBCClass
   use SMConstants
   use PhysicsStorage
   use FileReaders,            only: controlFileName
   use FileReadingUtilities,   only: GetKeyword, GetValueAsString, PreprocessInputLine
   use FTValueDictionaryClass, only: FTValueDictionary
   use Utilities, only: toLower, almostEqual
   use GenericBoundaryConditionClass
   use FluidData
   use VariableConversion
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
   public OutflowBC_t
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
   type, extends(GenericBC_t) ::  OutflowBC_t
#ifdef FLOW
      real(kind=RP)     :: pExt
#endif
      contains
         procedure         :: Destruct          => OutflowBC_Destruct
         procedure         :: Describe          => OutflowBC_Describe
#ifdef FLOW
         procedure         :: FlowState         => OutflowBC_FlowState
         procedure         :: FlowNeumann       => OutflowBC_FlowNeumann
#endif
#if defined(CAHNHILLIARD)
         procedure         :: PhaseFieldState   => OutflowBC_PhaseFieldState
         procedure         :: PhaseFieldNeumann => OutflowBC_PhaseFieldNeumann
         procedure         :: ChemPotState      => OutflowBC_ChemPotState
         procedure         :: ChemPotNeumann    => OutflowBC_ChemPotNeumann
#endif
   end type OutflowBC_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface OutflowBC_t
      module procedure ConstructOutflowBC
   end interface OutflowBC_t
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
      function ConstructOutflowBC(bname)
!
!        ********************************************************************
!        · Definition of the outflow boundary condition in the control file:
!              #define boundary bname
!                 type             = outflow
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
         type(OutflowBC_t)             :: ConstructOutflowBC
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

         ConstructOutflowBC % BCType = "outflow"
         ConstructOutflowBC % bname  = bname
         call toLower(ConstructOutflowBC % bname)

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
               inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(ConstructOutflowBC % bname)) 
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
#ifdef FLOW
         if ( bcdict % ContainsKey("pressure") ) then
            ConstructOutflowBC % pExt = bcdict % DoublePrecisionValueForKey("pressure")
         else if ( bcdict % ContainsKey("dp") ) then
#if defined(NAVIERSTOKES)
            ConstructOutflowBC % pExt = refValues % p / dimensionless % gammaM2 - bcdict % DoublePrecisionValueForKey("dp")
#elif defined(INCNS) || defined(MULTIPHASE)
            ConstructOutflowBC % pExt = 0.0_RP - bcdict % DoublePrecisionValueForKey("dp")
#endif
         else
#if defined(NAVIERSTOKES)
            ConstructOutflowBC % pExt = refValues % p / dimensionless % gammaM2 
#elif defined(INCNS) || defined(MULTIPHASE)
            ConstructOutflowBC % pExt = 0.0_RP
#endif
         end if

         ConstructOutflowBC % pExt = ConstructOutflowBC % pExt / refValues % p

#ifdef TRANSPORT
         call constructScalarBoundaryConditions(bcdict, ConstructOutflowBC)
#endif
#endif

         close(fid)
         call bcdict % Destruct
   
      end function ConstructOutflowBC

      subroutine OutflowBC_Describe(self)
!
!        ***************************************************
!              Describe the outflow boundary condition
!        ***************************************************
         implicit none
         class(OutflowBC_t),  intent(in)  :: self
#ifdef FLOW
         write(STD_OUT,'(30X,A,A28,A)') "->", " Boundary condition type: ", "Outflow"
         write(STD_OUT,'(30X,A,A28,F10.2)') "->", " Outflow pressure: ", self % pExt * refValues % p
#if defined(TRANSPORT)
         call self % ScalarBC % Describe()
#endif
#endif
         
      end subroutine OutflowBC_Describe

!
!/////////////////////////////////////////////////////////
!
!        Class destructors
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      subroutine OutflowBC_Destruct(self)
         implicit none
         class(OutflowBC_t)    :: self
#ifdef TRANSPORT
         call self % ScalarBC % Destruct()
#endif
      end subroutine OutflowBC_Destruct
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for compressible Navier--Stokes equations
!        -----------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES)
      subroutine OutflowBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(OutflowBC_t),   intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
!   
!        ---------------
!        Local Variables
!        ---------------
!   
         REAL(KIND=RP) :: qDotN, qTanx, qTany, qTanz, p, a, a2, eddy_theta
         REAL(KIND=RP) :: rPlus, entropyConstant, u, v, w, rho, normalMachNo
!         
         associate ( gammaMinus1 => thermodynamics % gammaMinus1, &
                     gamma => thermodynamics % gamma )
         
         qDotN = (nHat(1)*Q(2) + nHat(2)*Q(3) + nHat(3)*Q(4))/Q(1)
         qTanx = Q(2)/Q(1) - qDotN*nHat(1)
         qTany = Q(3)/Q(1) - qDotN*nHat(2)
         qTanz = Q(4)/Q(1) - qDotN*nHat(3)
         
         p            = gammaMinus1*( Q(5) - 0.5_RP*(Q(2)**2 + Q(3)**2 + Q(4)**2)/Q(1) )
         a2           = gamma*p/Q(1)
         a            = SQRT(a2)
         normalMachNo = ABS(qDotN/a)
         
         IF ( normalMachNo <= 1.0_RP )     THEN
!   
!           -------------------------------
!           Quantities coming from upstream
!           -------------------------------
!   
            rPlus           = qDotN + 2.0_RP*a/gammaMinus1
            entropyConstant = p - a2*Q(1)
!   
!           ----------------
!           Resolve solution
!           ----------------
!   
            rho   = -(entropyConstant - self % pExt)/a2
            a     = SQRT(gamma*self % pExt/rho)
            qDotN = rPlus - 2.0_RP*a/gammaMinus1
            u     = qTanx + qDotN*nHat(1)
            v     = qTany + qDotN*nHat(2)
            w     = qTanz + qDotN*nHat(3)
#if defined(SPALARTALMARAS)
            eddy_theta = Q(6)/Q(1)
#endif
            Q(1) = rho
            Q(2) = rho*u
            Q(3) = rho*v
            Q(4) = rho*w
            Q(5) = (self % pExt)/gammaMinus1 + 0.5_RP*rho*(u*u + v*v + w*w)
#if defined(SPALARTALMARAS)
            Q(6) = eddy_theta * rho 
#endif
         END IF
#if defined(TRANSPORT)
         call self % ScalarBC % FlowState(x, t, nHat, Q)
#endif   
      end associate

      end subroutine OutflowBC_FlowState

      subroutine OutflowBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(OutflowBC_t),   intent(in)   :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCONS)
         real(kind=RP),       intent(in)    :: U_x(NGRAD)
         real(kind=RP),       intent(in)    :: U_y(NGRAD)
         real(kind=RP),       intent(in)    :: U_z(NGRAD)
         real(kind=RP),       intent(inout) :: flux(NCONS)

         flux = 0.0_RP
#if defined(TRANSPORT)
         call self % ScalarBC % FlowNeumann(x, t, nHat, Q, U_x, U_y, U_z, flux)
#endif     
      end subroutine OutflowBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for incompressible Navier--Stokes equations
!        -------------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(INCNS)
      subroutine OutflowBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: u, v, w, theta, phi, un

         un = Q(INSRHOU)*nHat(IX) + Q(INSRHOV)*nHat(IY) + Q(INSRHOW)*nHat(IZ)
         
         if ( un .ge. -1.0e-4_RP ) then
         
            Q(INSP) = self % pExt

         else

            theta = refValues % AoATheta * PI / 180.0_RP
            phi   = refValues % AoAPhi   * PI / 180.0_RP
   
            u = cos(theta) * cos(phi)
            v = sin(theta) * cos(phi)
            w = sin(phi)

            Q(INSRHOU) = Q(INSRHO)*u
            Q(INSRHOV) = Q(INSRHO)*v
            Q(INSRHOW) = Q(INSRHO)*w

         end if

      end subroutine OutflowBC_FlowState

      subroutine OutflowBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCONS)
         real(kind=RP),       intent(in)    :: U_x(NCONS)
         real(kind=RP),       intent(in)    :: U_y(NCONS)
         real(kind=RP),       intent(in)    :: U_z(NCONS)
         real(kind=RP),       intent(inout) :: flux(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: drhodn, dudn, dvdn, dwdn

         flux = 0.0_RP

      end subroutine OutflowBC_FlowNeumann
#endif
!
!////////////////////////////////////////////////////////////////////////////
!
!  Subroutines for the Multiphase solver
!  -------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#ifdef MULTIPHASE
      subroutine OutflowBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: u, v, w, theta, phi, un, rho, sqrtRho


         un = Q(IMSQRHOU)*nHat(IX) + Q(IMSQRHOV)*nHat(IY) + Q(IMSQRHOW)*nHat(IZ)
         
         if ( un .ge. -1.0e-4_RP ) then
         
            Q(IMP) = self % pExt

         else

            u = dimensionless % vel_dir(IX)
            v = dimensionless % vel_dir(IY)
            w = dimensionless % vel_dir(IZ)

            rho = dimensionless % rho(1)*Q(IMC) + dimensionless % rho(2)*(1.0_RP - Q(IMC))
            sqrtRho = sqrt(rho)

            Q(IMSQRHOU) = sqrtRho*u
            Q(IMSQRHOV) = sqrtRho*v
            Q(IMSQRHOW) = sqrtRho*w

         end if


      end subroutine OutflowBC_FlowState

      subroutine OutflowBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCONS)
         real(kind=RP),       intent(in)    :: U_x(NCONS)
         real(kind=RP),       intent(in)    :: U_y(NCONS)
         real(kind=RP),       intent(in)    :: U_z(NCONS)
         real(kind=RP),       intent(inout) :: flux(NCONS)

         flux = 0.0_RP

      end subroutine OutflowBC_FlowNeumann
#endif

!
!////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for Cahn--Hilliard: all do--nothing in Outflows
!        ----------------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////////
!
#if defined(CAHNHILLIARD)
      subroutine OutflowBC_PhaseFieldState(self, x, t, nHat, Q)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine OutflowBC_PhaseFieldState

      subroutine OutflowBC_PhaseFieldNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCOMP)
         real(kind=RP),       intent(in)    :: U_x(NCOMP)
         real(kind=RP),       intent(in)    :: U_y(NCOMP)
         real(kind=RP),       intent(in)    :: U_z(NCOMP)
         real(kind=RP),       intent(inout) :: flux(NCOMP)

         flux = 0.0_RP
   
      end subroutine OutflowBC_PhaseFieldNeumann

      subroutine OutflowBC_ChemPotState(self, x, t, nHat, Q)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCOMP)
      end subroutine OutflowBC_ChemPotState

      subroutine OutflowBC_ChemPotNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCOMP)
         real(kind=RP),       intent(in)    :: U_x(NCOMP)
         real(kind=RP),       intent(in)    :: U_y(NCOMP)
         real(kind=RP),       intent(in)    :: U_z(NCOMP)
         real(kind=RP),       intent(inout) :: flux(NCOMP)

         flux = 0.0_RP

      end subroutine OutflowBC_ChemPotNeumann
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
      ! prescribed zero fluctuations of all variables, only makes sence for far field BC, possible after sponge
      subroutine OutflowBC_FlowState(self, x, t, nHat, Q)
         implicit none
         class(OutflowBC_t),  intent(in)    :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(inout) :: Q(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: V, rho, p

         V = 0.0_RP
         rho = 0.0_RP
         p = 0.0_RP

         Q(ICAARHO) = rho
         Q(ICAAU) = V
         Q(ICAAV) = V
         Q(ICAAW) = V
         Q(ICAAP) = p

      end subroutine OutflowBC_FlowState

      subroutine OutflowBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
         implicit none
         class(OutflowBC_t),  intent(in)     :: self
         real(kind=RP),       intent(in)    :: x(NDIM)
         real(kind=RP),       intent(in)    :: t
         real(kind=RP),       intent(in)    :: nHat(NDIM)
         real(kind=RP),       intent(in)    :: Q(NCONS)
         real(kind=RP),       intent(in)    :: U_x(NCONS)
         real(kind=RP),       intent(in)    :: U_y(NCONS)
         real(kind=RP),       intent(in)    :: U_z(NCONS)
         real(kind=RP),       intent(inout) :: flux(NCONS)

         flux = 0.0_RP

      end subroutine OutflowBC_FlowNeumann
#endif
!
end module OutflowBCClass
