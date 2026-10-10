#include "Includes.h"
module AutosaveClass
   use SMConstants
   use FTValueDictionaryClass
   use Utilities, only: almostEqual, toLower

   private
   public   Autosave_t, AUTOSAVE_BY_TIME, AUTOSAVE_BY_ITERATION, AUTOSAVE_UNDEFINED
   public   AutosaveSlidingMesh_t
   public   Autosave_LambVectorConfigure

   integer, parameter :: AUTOSAVE_UNDEFINED    = 0
   integer, parameter :: AUTOSAVE_BY_TIME      = 1
   integer, parameter :: AUTOSAVE_BY_ITERATION = 2

   type Autosave_t
      logical        :: performAutosave
      logical        :: enable
      integer        :: iter_interval
      real(kind=RP)  :: time_interval
      real(kind=RP)  :: nextAutosaveTime
      integer        :: mode
      contains
         procedure   :: Configure => Autosave_Configure
         procedure   :: Autosave  => Autosave_Autosave
   end type Autosave_t

   type, extends(Autosave_t)  :: AutosaveSlidingMesh_t
         contains
         procedure   :: Configure => AutosaveSlidingMesh_Configure
   end type AutosaveSlidingMesh_t

   contains
!
!////////////////////////////////////////////////////////////////////////////////////
!
!     Autosave procedures
!     -------------------
!
!////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Autosave_Configure(self, controlVariables, t0)
      implicit none
      class(Autosave_t)          :: self
      class(FTValueDictionary)   :: controlVariables
      real(kind=RP), intent(in)  :: t0
!
!     ---------------
!     Local variables
!     ---------------
!
      class(FTObject), pointer   :: obj
      character(len=LINE_LENGTH) :: autosaveMode
      character(len=LINE_LENGTH), parameter :: autosaveModeKey = "autosave mode"
      character(len=LINE_LENGTH), parameter :: autosaveIntervalKey = "autosave interval"
      character(len=LINE_LENGTH), parameter :: autosaveByIteration = "iteration"
      character(len=LINE_LENGTH), parameter :: autosaveByTime = "time"
!
!     Check whether the autosave mode is present
!     ------------------------------------------
      obj => controlVariables % objectForKey(trim(autosaveModeKey))

      if ( associated(obj) ) then
!
!        Present: associate the appropriate mode
!        ---------------------------------------
         autosaveMode = controlVariables % stringValueForKey(trim(autosaveModeKey), requestedLength = LINE_LENGTH)
         call ToLower(autosaveMode)

         if ( trim(autosaveMode) .eq. trim(autosaveByIteration) ) then
            self % mode = AUTOSAVE_BY_ITERATION

         elseif ( trim(autosaveMode) .eq. trim(autosaveByTime) ) then
            self % mode = AUTOSAVE_BY_TIME

         else
            print*, 'Unknown autosave mode "',trim(autosaveMode),'".'
            print*, "Implemented modes are:"
            print*, "   * iteration"
            print*, "   * time"
            errorMessage(STD_OUT)
            error stop
         end if
   
      else 
!
!        Not present: undefined. It can be determined with the interval data kind
!        ------------------------------------------------------------------------
         self % mode = AUTOSAVE_UNDEFINED
   
      end if
!
!     Check whether the autosave interval is present
!     ----------------------------------------------
      obj => controlVariables % objectForKey(trim(autosaveIntervalKey))

      if ( associated(obj) ) then
!
!        Present: assign to time/iter depending on the mode (see below for the undefined case)
!        --------------------------------------------------         
         select case ( self % mode )

         case (AUTOSAVE_BY_ITERATION)
            self % enable = .true.
            self % iter_interval = controlVariables % integerValueForKey(trim(autosaveIntervalKey))
            self % time_interval = huge(1.0_RP)
   
         case (AUTOSAVE_BY_TIME)
            self % enable = .true.
            self % time_interval = controlVariables % doublePrecisionValueForKey(trim(autosaveIntervalKey))
            self % iter_interval = huge(1)
         
         case (AUTOSAVE_UNDEFINED)
!
!           If time_interval is integer, autosave by iteration is selected. Otherwise, autosave by time
!           -------------------------------------------------------------------------------------------
            self % time_interval = controlVariables % doublePrecisionValueForKey(trim(autosaveIntervalKey))

            if ( almostEqual( self % time_interval, fraction(self % time_interval) ) ) then
               self % enable = .true.
               self % mode = AUTOSAVE_BY_ITERATION
               self % iter_interval = controlVariables % integerValueForKey(trim(autosaveIntervalKey))
               self % time_interval = huge(1.0_RP)

            else
               self % enable = .true.
               self % mode = AUTOSAVE_BY_TIME
               self % iter_interval = huge(1)

            end if

         end select

      else  
!
!        Not present: disable autosave
!        -----------------------------
         self % enable = .false.
         self % mode   = AUTOSAVE_UNDEFINED
         self % iter_interval = huge(1)
         self % time_interval = huge(1.0_RP)

      end if
!
!     Reset the last autosave time
!     ----------------------------
      self % nextAutosaveTime = t0 + self % time_interval
      self % performAutosave = .false.
         
   end subroutine Autosave_Configure

   logical function Autosave_Autosave(self,iter)
      class(Autosave_t),   intent(in)     :: self
      integer,             intent(in)     :: iter

      if ( .not. self % enable ) then
         Autosave_Autosave = .false.
         return
      end if

      select case ( self % mode )
      case (AUTOSAVE_BY_ITERATION)
         if ( mod(iter,self % iter_interval) .eq. 0 ) then
            Autosave_Autosave = .true.
         else
            Autosave_Autosave = .false.
         end if
   
      case (AUTOSAVE_BY_TIME)
         Autosave_Autosave = self % performAutosave
      
      end select
   end function Autosave_Autosave

   subroutine Autosave_LambVectorConfigure(self, controlVariables, t0)
      use mainKeywordsModule, only: saveLambVectorKey
      implicit none
      class(Autosave_t)          :: self
      class(FTValueDictionary)   :: controlVariables
      real(kind=RP), intent(in)  :: t0
!
!     ---------------
!     Local variables
!     ---------------
!
      character(len=LINE_LENGTH) :: autosaveMode
      CHARACTER(LEN=*), PARAMETER :: saveLambVectorModeKey = "save lamb vector mode"
      CHARACTER(LEN=*), PARAMETER :: saveLambVectorIntervalKey = "save lamb vector interval"
      character(len=LINE_LENGTH), parameter :: autosaveByIteration = "iteration"
      character(len=LINE_LENGTH), parameter :: autosaveByTime = "time"

      if (.not. controlVariables % logicalValueForKey(saveLambVectorKey)) then
         self % enable = .false.
         self % mode   = AUTOSAVE_UNDEFINED
         self % iter_interval = huge(1)
         self % time_interval = huge(1.0_RP)
         self % nextAutosaveTime = t0 + self % time_interval
         self % performAutosave = .false.
         return
      end if
!
!     Check whether the autosave mode is present
!     ------------------------------------------
      if ( controlVariables % containsKey(trim(saveLambVectorModeKey)) ) then
!
!        Present: associate the appropriate mode
!        ---------------------------------------
         autosaveMode = controlVariables % stringValueForKey(trim(saveLambVectorModeKey), requestedLength = LINE_LENGTH)
         call ToLower(autosaveMode)

         if ( trim(autosaveMode) .eq. trim(autosaveByIteration) ) then
            self % mode = AUTOSAVE_BY_ITERATION

         elseif ( trim(autosaveMode) .eq. trim(autosaveByTime) ) then
            self % mode = AUTOSAVE_BY_TIME

         else
            print*, 'Unknown autosave mode for the Lamb vector"',trim(autosaveMode),'".'
            print*, "Implemented modes are:"
            print*, "   * iteration"
            print*, "   * time"
            errorMessage(STD_OUT)
            error stop
         end if
   
      else 
            print *, trim(saveLambVectorModeKey), " should be specified to save the Lamb vector."
            errorMessage(STD_OUT)
            error stop
      end if
!
!     Check whether the autosave interval is present
!     ----------------------------------------------
      if ( controlVariables % containsKey(trim(saveLambVectorIntervalKey)) ) then
!
!        Present: assign to time/iter depending on the mode (see below for the undefined case)
!        --------------------------------------------------         
         select case ( self % mode )

         case (AUTOSAVE_BY_ITERATION)
            self % enable = .true.
            self % iter_interval = controlVariables % integerValueForKey(trim(saveLambVectorIntervalKey))
            self % time_interval = huge(1.0_RP)
   
         case (AUTOSAVE_BY_TIME)
            self % enable = .true.
            self % time_interval = controlVariables % doublePrecisionValueForKey(trim(saveLambVectorIntervalKey))
            self % iter_interval = huge(1)

         end select

      else  
            print *, trim(saveLambVectorIntervalKey), " should be specified to save the Lamb vector."
            errorMessage(STD_OUT)
            error stop
      end if
!
!     Reset the last autosave time
!     ----------------------------
      self % nextAutosaveTime = t0 + self % time_interval
      self % performAutosave = .false.
         
   end subroutine Autosave_LambVectorConfigure

   subroutine AutosaveSlidingMesh_Configure(self, controlVariables, t0)
      use ParamfileRegions
      implicit none
      class(AutosaveSlidingMesh_t)          :: self
      class(FTValueDictionary)   :: controlVariables
      real(kind=RP), intent(in)  :: t0
!
!     ---------------
!     Local variables
!     ---------------
!
      character(len=LINE_LENGTH) :: in_label, paramFile
      character(len=LINE_LENGTH) :: autosaveMode
      integer, allocatable :: readIteration
      real(kind=rp), allocatable :: readTime
      character(len=LINE_LENGTH), parameter :: autosaveModeKey = "autosave mode"
      character(len=LINE_LENGTH), parameter :: autosaveIntervalKey = "autosave interval"
      character(len=LINE_LENGTH), parameter :: autosaveByIteration = "iteration"
      character(len=LINE_LENGTH), parameter :: autosaveByTime = "time"

      ! Default values
      self % enable = .false.
      self % mode = AUTOSAVE_UNDEFINED
      self % iter_interval = huge(1)
      self % time_interval = huge(1.0_RP)

      ! Read values
      write(in_label , '(A)') "#define slidingmesh"
      call get_command_argument(1, paramFile) 

      call readValueInRegion( trim( paramFile ), autosaveModeKey, autosaveMode, in_label, "#end" )

      if ( trim(autosaveMode) .eq. trim(autosaveByIteration) ) then
         self % enable = .true.
         self % mode = AUTOSAVE_BY_ITERATION
         call readValueInRegion( trim( paramFile ), autosaveIntervalKey, readIteration, in_label, "#end" )
         self % iter_interval = readIteration
      elseif ( trim(autosaveMode) .eq. trim(autosaveByTime) ) then
         self % enable = .true.
         self % mode = AUTOSAVE_BY_TIME
         call readValueInRegion( trim( paramFile ), autosaveIntervalKey, readTime, in_label, "#end" )
         self % time_interval = readTime
      end if

!
!     Reset the last autosave time
!     ----------------------------
      self % nextAutosaveTime = t0 + self % time_interval
      self % performAutosave = .false.
         
   end subroutine AutosaveSlidingMesh_Configure

end module AutosaveClass