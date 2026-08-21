#include "Includes.h"
module ScalarBoundaryConditions
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
!    public DirichletBC_t, NeumannBC_t
   public constructScalarBoundaryConditions
!
!  ****************
!  Class definition
!  ****************
!

    type, extends(GenericBC_t) ::  DirichletBC_t
            real(kind=RP)              :: value ! Impose c = c_D. c_D is the value read from the control file.
        contains
            procedure         :: Destruct          => DirichletBC_Destruct
            procedure         :: Describe          => DirichletBC_Describe
            procedure         :: FlowState         => DirichletBC_FlowState
            procedure         :: FlowGradVars      => DirichletBC_FlowGradVars
            procedure         :: FlowNeumann       => DirichletBC_FlowNeumann
    end type DirichletBC_t

    type, extends(GenericBC_t) ::  NeumannBC_t
            real(kind=RP) :: flux ! Impose J * n = q_N, where J = D * grad(c). q_N is the value read from the control file.
        contains
            procedure         :: Destruct          => NeumannBC_Destruct
            procedure         :: Describe          => NeumannBC_Describe
            procedure         :: FlowState         => NeumannBC_FlowState
            procedure         :: FlowGradVars      => NeumannBC_FlowGradVars
            procedure         :: FlowNeumann       => NeumannBC_FlowNeumann
    end type NeumannBC_t

    !  Constructors are exported with the name of the class
    interface DirichletBC_t
        module procedure ConstructDirichletBC
    end interface DirichletBC_t

    interface NeumannBC_t
        module procedure ConstructNeumannBC
    end interface NeumannBC_t

   contains

    subroutine constructScalarBoundaryConditions(bcdict, FlowBC)
        use FTValueDictionaryClass, only: FTValueDictionary
        use GenericBoundaryConditionClass, only: GenericBC_t
        implicit none
        type(FTValueDIctionary), intent(in)    :: bcdict
        class(GenericBC_t), intent(inout) :: FlowBC

        character(LEN=LINE_LENGTH) :: keyword

        if (.not. bcdict % ContainsKey("scalar")) then
            print *, "You should specify a boundary condition for the scalar equation. Options are:"
            print *, "scalar = dirichlet"
            print *, "scalar = neumann"
            print *, "scalar = periodic"
            errorMessage(STD_OUT)
            error stop
        endif
        
        keyword = bcdict % StringValueForKey("scalar", LINE_LENGTH)
        call toLower(keyword)
        select case (trim(keyword))
        case ("dirichlet")
            allocate(DirichletBC_t :: FlowBC % ScalarBC)
            ! FlowBC % ScalarBC = DirichletBC_t(bcdict)
            select type(bc => FlowBC % ScalarBC)
            type is (DirichletBC_t)
                bc = DirichletBC_t(bcdict)
            end select
        case ("neumann")
            allocate(NeumannBC_t :: FlowBC % ScalarBC)
            ! FlowBC % ScalarBC = NeumannBC_t(bcdict)
            select type(bc => FlowBC % ScalarBC)
            type is (NeumannBC_t)
                bc = NeumannBC_t(bcdict)
            end select
        end select
    end subroutine constructScalarBoundaryConditions

!/////////////////////////////////////////////////////////
!
!        DirichletBC_t
!        -----------------
!
!/////////////////////////////////////////////////////////
!
    function ConstructDirichletBC(bcdict)
        implicit none
        type(DirichletBC_t)             :: ConstructDirichletBC
        type(FTValueDIctionary), intent(in)    :: bcdict      
        ConstructDirichletBC % value = bcdict % getValueOrDefault("scalar value", 0.0_RP)
    end function ConstructDirichletBC

    subroutine DirichletBC_Describe(self)
        implicit none
        class(DirichletBC_t),  intent(in)  :: self
        write(STD_OUT,'(30X,A,A28,A)') "->", " Scalar BC type: ", "Dirichlet"
        write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Dirichlet value: ', self % value         
    end subroutine DirichletBC_Describe

    subroutine DirichletBC_Destruct(self)
        implicit none
        class(DirichletBC_t)    :: self
    end subroutine DirichletBC_Destruct

    subroutine DirichletBC_FlowState(self, x, t, nHat, Q)
        implicit none
        class(DirichletBC_t),  intent(in)    :: self
        real(kind=RP),          intent(in)    :: x(NDIM)
        real(kind=RP),          intent(in)    :: t
        real(kind=RP),          intent(in)    :: nHat(NDIM)
        real(kind=RP),          intent(inout) :: Q(NCONS)
#ifdef TRANSPORT
        ! Q(IC) = self % value
        Q(IC) = 2.0_RP * self % value - Q(IC)
#endif
    end subroutine DirichletBC_FlowState

    subroutine DirichletBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
        implicit none
        class(DirichletBC_t),  intent(in)    :: self
        real(kind=RP),          intent(in)    :: x(NDIM)
        real(kind=RP),          intent(in)    :: t
        real(kind=RP),          intent(in)    :: nHat(NDIM)
        real(kind=RP),          intent(in)    :: Q(NCONS)
        real(kind=RP),          intent(inout) :: U(NGRAD)
        procedure(GetGradientValues_f)        :: GetGradients
        
#ifdef TRANSPORT
        ! U(IC) = self % value
#endif
    end subroutine DirichletBC_FlowGradVars

    subroutine DirichletBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
        implicit none
        class(DirichletBC_t),   intent(in)    :: self
        real(kind=RP),       intent(in)    :: x(NDIM)
        real(kind=RP),       intent(in)    :: t
        real(kind=RP),       intent(in)    :: nHat(NDIM)
        real(kind=RP),       intent(in)    :: Q(NCONS)
        real(kind=RP),       intent(in)    :: U_x(NGRAD)
        real(kind=RP),       intent(in)    :: U_y(NGRAD)
        real(kind=RP),       intent(in)    :: U_z(NGRAD)
        real(kind=RP),       intent(inout) :: flux(NCONS)
    end subroutine DirichletBC_FlowNeumann


!/////////////////////////////////////////////////////////
!
!        NeumannBC_t
!        -----------------
!
!/////////////////////////////////////////////////////////
!
    function ConstructNeumannBC(bcdict)
        implicit none
        type(NeumannBC_t)             :: ConstructNeumannBC
        type(FTValueDIctionary), intent(in)    :: bcdict      
        ConstructNeumannBC % flux = bcdict % getValueOrDefault("scalar value", 0.0_RP)
    end function ConstructNeumannBC

    subroutine NeumannBC_Describe(self)
        implicit none
        class(NeumannBC_t),  intent(in)  :: self
        write(STD_OUT,'(30X,A,A28,A)') "->", " Scalar BC type: ", "Neumann"
        write(STD_OUT,'(30X,A,A28,F10.2)') "->", ' Flux value: ', self % flux
    end subroutine NeumannBC_Describe

    subroutine NeumannBC_Destruct(self)
        implicit none
        class(NeumannBC_t)    :: self
    end subroutine NeumannBC_Destruct

    subroutine NeumannBC_FlowState(self, x, t, nHat, Q)
        implicit none
        class(NeumannBC_t),  intent(in)    :: self
        real(kind=RP),          intent(in)    :: x(NDIM)
        real(kind=RP),          intent(in)    :: t
        real(kind=RP),          intent(in)    :: nHat(NDIM)
        real(kind=RP),          intent(inout) :: Q(NCONS)
    end subroutine NeumannBC_FlowState

    subroutine NeumannBC_FlowGradVars(self, x, t, nHat, Q, U, GetGradients)
        implicit none
        class(NeumannBC_t),  intent(in)    :: self
        real(kind=RP),          intent(in)    :: x(NDIM)
        real(kind=RP),          intent(in)    :: t
        real(kind=RP),          intent(in)    :: nHat(NDIM)
        real(kind=RP),          intent(in)    :: Q(NCONS)
        real(kind=RP),          intent(inout) :: U(NGRAD)
        procedure(GetGradientValues_f)        :: GetGradients
    end subroutine NeumannBC_FlowGradVars

    subroutine NeumannBC_FlowNeumann(self, x, t, nHat, Q, U_x, U_y, U_z, flux)
        implicit none
        class(NeumannBC_t),   intent(in)    :: self
        real(kind=RP),       intent(in)    :: x(NDIM)
        real(kind=RP),       intent(in)    :: t
        real(kind=RP),       intent(in)    :: nHat(NDIM)
        real(kind=RP),       intent(in)    :: Q(NCONS)
        real(kind=RP),       intent(in)    :: U_x(NGRAD)
        real(kind=RP),       intent(in)    :: U_y(NGRAD)
        real(kind=RP),       intent(in)    :: U_z(NGRAD)
        real(kind=RP),       intent(inout) :: flux(NCONS)
#ifdef TRANSPORT
        flux(IC) = self % flux
#endif
    end subroutine NeumannBC_FlowNeumann

end module ScalarBoundaryConditions
