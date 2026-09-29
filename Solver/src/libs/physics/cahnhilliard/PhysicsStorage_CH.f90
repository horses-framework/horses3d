#include "Includes.h"
      Module Physics_CHKeywordsModule
         IMPLICIT NONE 
         INTEGER, PARAMETER :: KEYWORD_LENGTH = 132
!
!        ******************
!        Required arguments
!        ******************
!
         character(len=KEYWORD_LENGTH), parameter    :: TCH_KEY               = "chem characteristic time (s)"
         character(len=KEYWORD_LENGTH), parameter    :: INTERFACE_WIDTH_KEY   = "interface width (m)"
         character(len=KEYWORD_LENGTH), parameter    :: INTERFACE_TENSION_KEY = "interface tension (n/m)"
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(3) :: physics_CHKeywords = [INTERFACE_WIDTH_KEY, &
                                                                              TCH_KEY,   &
                                                                              INTERFACE_TENSION_KEY ]
      END MODULE Physics_CHKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!    
!    ******
     MODULE PhysicsStorage_CH
!    ******
!
     USE SMConstants
     USE Physics_CHKeywordsModule
     use FluidData_CH
     use Utilities,            only: toLower, almostEqual
     
     IMPLICIT NONE

     private
     public    NCOMP

     public    ConstructPhysicsStorage_CH, DestructPhysicsStorage_CH, DescribePhysicsStorage_CH
     public    CheckPhysicsCHInputIntegrity

     integer, parameter    :: NCOMP = 1
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
      SUBROUTINE ConstructPhysicsStorage_CH( controlVariables, Lref, timeRef, pRef, success )
      USE FTValueDictionaryClass
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FTValueDictionary)      :: controlVariables
      real(kind=RP),    intent(in) :: Lref, timeRef, pRef
      LOGICAL                      :: success
!
!     ---------------
!     Local variables
!     ---------------
!
      CHARACTER(LEN=KEYWORD_LENGTH) :: keyword
      type(Multiphase_t)            :: multiphase_ 

      multiphase_ = Multiphase_t()
!
!     --------------------
!     Collect input values
!     --------------------
!
      success = .TRUE.
      CALL CheckPhysicsCHInputIntegrity(controlVariables,success)
      IF(.NOT. success) RETURN 
!
!     *****************************
!     Read multiphase properties
!     *****************************
!
      multiphase_ % eps_wDim   = controlVariables % DoublePrecisionValueForKey(INTERFACE_WIDTH_KEY)
      multiphase_ % tCH_wDim   = controlVariables % DoublePrecisionValueForKey(TCH_KEY)
      multiphase_ % sigma_wDim = controlVariables % DoublePrecisionValueForKey(INTERFACE_TENSION_KEY)


      if ( almostEqual(multiphase_ % sigma_wDim, 0.0_RP) ) then
         multiphase_ % M0_wDim = 0.0_RP
      else if ( .not. almostEqual(multiphase_ % tCH_wDim, 0.0_RP) ) then
         multiphase_ % M0_wDim = POW2(Lref)*(multiphase_ % eps_wDim**2) / (multiphase_ % tCH_wDim * multiphase_ % sigma_wDim)
      else
         multiphase_ % M0_wDim = 0.0_RP
      end if


      multiphase_ % eps   = multiphase_ % eps_wDim / Lref
      multiphase_ % tCH   = multiphase_ % tCH_wDim / timeRef
      multiphase_ % sigma = multiphase_ % sigma_wDim / pRef

      if ( almostEqual(multiphase_ % sigma_wDim, 0.0_RP) ) then
         multiphase_ % M0 = 0.0_RP
      else if ( .not. almostEqual(multiphase_ % tCH, 0.0_RP) ) then
         multiphase_ % M0 = (multiphase_ % eps**2) / (multiphase_ % tCH * multiphase_ % sigma)
      else
         multiphase_ % M0 = 0.0_RP
      end if

!
!     ************************************
!     Set the global (proteted) multiphase
!     ************************************
!
      call setMultiphase( multiphase_ )

      END SUBROUTINE ConstructPhysicsStorage_CH
!
!     ///////////////////////////////////////////////////////
!
!     -------------------------------------------------
!!    Destructor: Does nothing for this storage
!     -------------------------------------------------
!
      SUBROUTINE DestructPhysicsStorage_CH
      
      END SUBROUTINE DestructPhysicsStorage_CH
!
!     //////////////////////////////////////////////////////
!
!     -----------------------------------------
!!    Descriptor: Shows the gathered data
!     -----------------------------------------
!
      SUBROUTINE DescribePhysicsStorage_CH(Lref)
         USE Headers
         use MPI_Process_Info
         IMPLICIT NONE
         real(kind=RP), intent(in)  :: Lref

         if ( .not. MPI_Process % isRoot ) return 

         call Section_Header("Loading Cahn-Hilliard physics")

         write(STD_OUT,'(/,/)')

         call SubSection_Header("Chemical properties")
         write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Chemical characteristic time: " , multiphase % tCH_wDim, " (s)."

         if ( almostEqual(multiphase % tCH_wDim, 0.0_RP) ) then
            write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Mobility diffusion: " , multiphase % M0_wDim, " (m^2/s)."
         else
            write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Mobility diffusion: " , POW2(Lref)/multiphase % tCH_wDim, " (m^2/s)."
         end if
         write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Interface tension: " , multiphase % sigma_wDim, " (N/m)."
         write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Interface width: " , multiphase % eps_wDim, " (m)."

         write(STD_OUT,'(/)')
         call SubSection_Header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Interface tension: " , multiphase % sigma, " [-]."
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Epsilon: " , multiphase % eps, " [-]."
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Mobility: " , multiphase % M0, " [-]."

      END SUBROUTINE DescribePhysicsStorage_CH
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE CheckPhysicsCHInputIntegrity( controlVariables, success )  
         USE FTValueDictionaryClass
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
         
         DO i = 1, SIZE(physics_CHKeywords)
            obj => controlVariables % objectForKey(physics_CHKeywords(i))
            IF ( .NOT. ASSOCIATED(obj) )     THEN
               PRINT *, "Input file is missing entry for keyword: ",physics_CHKeywords(i)
               success = .FALSE. 
               errorMessage(STD_OUT)
               error stop
            END IF  
         END DO  
         
      END SUBROUTINE CheckPhysicsCHInputIntegrity
!
!    **********       
     END MODULE PhysicsStorage_CH
!    **********
