#include "Includes.h"
module PhysicsStorage
   use SMConstants     , only: RP, STD_OUT
   use Headers
   use MPI_Process_Info, only: MPI_Process
#ifdef FLOW
   use FluidData, only: refValues, thermodynamics
#endif
#ifdef CAHNHILLIARD
   use FluidData, only: multiphase
#endif
#if defined(SPALARTALMARAS)
   use PhysicsStorage_NSSA
#elif defined(NAVIERSTOKES) && !defined(TRANSPORT)
   use PhysicsStorage_NS
#elif defined(NAVIERSTOKES) && defined(TRANSPORT)
   use PhysicsStorage_NSTPT
#elif defined(INCNS)
   use PhysicsStorage_iNS
#elif defined(MULTIPHASE)
   use PhysicsStorage_MU
#elif defined(ACOUSTIC)
   use PhysicsStorage_CAA
#endif
#ifdef CAHNHILLIARD
   use PhysicsStorage_CH
#endif
   implicit none

#ifdef FLOW
   private refValues
#endif

#ifdef CAHNHILLIARD
   private multiphase
#endif
!
!  *****************************************************************************
!  These are the different modes supported by the ComputeTimeDerivative function
!        Defining this avoids to create specific procedures for each
!  *****************************************************************************
!
   enum, bind(C)
      enumerator :: CTD_IGNORE_MODE
#ifdef CAHNHILLIARD
      enumerator :: CTD_IMEX_EXPLICIT, CTD_IMEX_IMPLICIT
#endif
      enumerator :: CTD_LAPLACIAN
      enumerator :: CTD_DUMMY
   end enum

   real(kind=RP), protected     :: Lref
   real(kind=RP), protected     :: timeref

   character(len=*), parameter   :: REFERENCE_LENGTH_KEY = "reference length (m)" 
   CHARACTER(LEN=*), parameter   :: LambVectorCubeKey = "lamb vector cube"


#if (!defined(FLOW)) && (defined(CAHNHILLIARD))
   integer, parameter :: NCONS = NCOMP
   integer, parameter :: NGRAD = NCOMP
#endif
   
   
   contains
      subroutine ConstructPhysicsStorage(controlVariables, success)
         USE FTValueDictionaryClass
         TYPE(FTValueDictionary)      :: controlVariables
         LOGICAL                      :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: timeRef_NS, pRef
         

         if ( controlVariables % ContainsKey(REFERENCE_LENGTH_KEY) ) then
            Lref = controlVariables % DoublePrecisionValueForKey(REFERENCE_LENGTH_KEY)
         else
            Lref = 1.0_RP
         end if
!
!        Construct NSE physics
!        ---------------------
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS)) && (!(TRANSPORT))
         call ConstructPhysicsStorage_NS( controlVariables, Lref, timeRef_NS, success )
#elif defined(NAVIERSTOKES) && (TRANSPORT)
         call ConstructPhysicsStorage_NSTPT( controlVariables, Lref, timeRef_NS, success )
#elif defined(SPALARTALMARAS)
         call ConstructPhysicsStorage_NSSA( controlVariables, Lref, timeRef_NS, success )
#elif defined(INCNS)
         call ConstructPhysicsStorage_iNS( controlVariables, Lref, timeRef_NS, success )
#elif defined(MULTIPHASE)
         call ConstructPhysicsStorage_MU( controlVariables, Lref, timeRef_NS, success )
#elif defined(ACOUSTIC)
         call ConstructPhysicsStorage_CAA( controlVariables, Lref, timeRef_NS, success )
#endif

!        Navier--Stokes equations set the reference time
!        -----------------------------------------------
#ifdef FLOW
         timeRef = timeRef_NS
         pRef    = refValues % p
#else
         timeRef = 1.0_RP
         pRef    = 1.0_RP
#endif
!
!        Construct CHE physics
!        ---------------------
#ifdef CAHNHILLIARD
         call ConstructPhysicsStorage_CH( controlVariables, Lref, timeRef, pRef, success )
#endif
!
!        ****************
!        Describe physics      
!        ****************
!
         call DescribePhysicsStorage_Common
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS)) && (!(TRANSPORT))
         call DescribePhysicsStorage_NS()
#elif defined(NAVIERSTOKES) && (TRANSPORT)
         call DescribePhysicsStorage_NSTPT()
#elif defined(SPALARTALMARAS)
         call DescribePhysicsStorage_NSSA()
#elif defined(INCNS)
         call DescribePhysicsStorage_iNS()
#elif defined(MULTIPHASE)
         call DescribePhysicsStorage_MU()
#elif defined(ACOUSTIC)
         call DescribePhysicsStorage_CAA()
         call DescribePhysicsStorage_CAABaseFlow(controlVariables)
#endif

#ifdef CAHNHILLIARD
         call DescribePhysicsStorage_CH(Lref)
#endif
   
      end subroutine ConstructPhysicsStorage

      subroutine DescribePhysicsStorage_Common()
         implicit none
         
         if ( .not. MPI_Process % isRoot ) return 
         
         call Section_Header("Loading common physics")

         write(STD_OUT,'(/,/)')

         write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Reference length (m): " , Lref
         write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Reference time (s): "   , timeRef

      end subroutine DescribePhysicsStorage_Common

      subroutine SetReferenceLength(Lref_)
         implicit none
         real(kind=RP), intent(in)  :: Lref_

         Lref = Lref_

      end subroutine SetReferenceLength

end module PhysicsStorage
