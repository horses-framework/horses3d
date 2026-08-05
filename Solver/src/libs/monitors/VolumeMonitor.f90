#include "Includes.h"
module VolumeMonitorClass
   use SMConstants
   use HexMeshClass
   use MonitorDefinitions
   use PhysicsStorage
   use MPI_Process_Info
   implicit none

   private
!~   public VOLUME
#if defined(NAVIERSTOKES)
!~   public KINETIC_ENERGY, KINETIC_ENERGY_RATE, ENSTROPHY
!~   public ENTROPY, ENTROPY_RATE, INTERNAL_ENERGY
#elif defined(CAHNHILLIARD)
!~   public FREE_ENERGY
#endif
   public VolumeMonitor_t

!
!  *******************************
!  Volume monitor class definition
!  *******************************
!
   type VolumeMonitor_t
      logical                         :: active
      integer                         :: ID
      integer                         :: bufferLine
      integer                         :: num_of_vars
      real(kind=RP), allocatable      :: values(:,:)
      character(len=STR_LEN_MONITORS) :: monitorName
      character(len=STR_LEN_MONITORS) :: fileName
      character(len=STR_LEN_MONITORS) :: variable
      contains
         procedure   :: Initialization => VolumeMonitor_Initialization
         procedure   :: Update         => VolumeMonitor_Update
         procedure   :: WriteLabel     => VolumeMonitor_WriteLabel
         procedure   :: WriteValues    => VolumeMonitor_WriteValue
         procedure   :: WriteToFile    => VolumeMonitor_WriteToFile
         procedure   :: getLast        => VolumeMonitor_GetLast
         procedure   :: destruct       => VolumeMonitor_Destruct
         procedure   :: copy           => VolumeMonitor_Assign
         generic     :: assignment(=)  => copy
   end type VolumeMonitor_t
!
!  ========
   contains
!  ========
!
!
!///////////////////////////////////////////////////////////////////////////
!
!           VOLUME MONITOR PROCEDURES
!           -------------------------
!///////////////////////////////////////////////////////////////////////////
!
      subroutine VolumeMonitor_Initialization( self , mesh , ID, solution_file , FirstCall)
!
!        *****************************************************************************
!              This subroutine initializes the volume monitor. The following
!           data is obtained from the case file:
!              -> Name: The monitor name (10 characters maximum)
!              -> Variable: The variable to be monitorized.
!        *****************************************************************************
!
         use ParamfileRegions
         use Utilities, only: toLower
         implicit none
         class(VolumeMonitor_t) :: self
         class(HexMesh)         :: mesh
         integer                :: ID
         character(len=*)       :: solution_file
         logical, intent(in)    :: FirstCall
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=STR_LEN_MONITORS)  :: in_label
         character(len=STR_LEN_MONITORS)  :: fileName
         character(len=STR_LEN_MONITORS)  :: paramFile
         integer                          :: fID
         integer                          :: pos
!
!        Get monitor ID
!        --------------
         self % ID = ID
!
!        Search for the parameters in the case file
!        ------------------------------------------
         write(in_label , '(A,I0)') "#define volume monitor " , self % ID

         call get_command_argument(1, paramFile)
         call readValueInRegion ( trim ( paramFile )  , "name"              , self % monitorName      , in_label , "# end" )
         call readValueInRegion ( trim ( paramFile )  , "variable"          , self % variable         , in_label , "# end" )
!
!        Enable the monitor
!        ------------------
         self % active = .true.
         self % bufferLine = 1
!
!        Select the variable from the available list, and compute auxiliary variables if needed
!        --------------------------------------------------------------------------------------
         self % num_of_vars = 1
         call toLower(self % variable)

#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
         select case ( trim ( self % variable ) )
         case ("kinetic energy")
         case ("kinetic energy rate")
         case ("kinetic energy balance")
         case ("enstrophy")
         case ("entropy")
         case ("entropy rate")
         case ("entropy balance")
         case ("math entropy")
         case ("artificial dissipation")
         case ("internal energy")
         case ("relative kinetic energy")
         case ("relative total energy")
         case ("mean velocity")
         case ("velocity")          ; self % num_of_vars = 3
         case ("momentum")          ; self % num_of_vars = 3
         case ("source")            ; self % num_of_vars = NCONS
         case ("particles source")  ; self % num_of_vars = NCONS
         case ("sensor range")      ; self % num_of_vars = 2
         case ("l2rho")
         case ("l2rhou")
         case ("l2rhoe")

         case default

            if ( len_trim (self % variable) .eq. 0 ) then
               print*, "Variable was not specified for volume monitor " , self % ID , "."
            else
               print*, 'Variable "',trim(self % variable),'" volume monitor ', self % ID, ' not implemented yet.'
               print*, "Options available are:"
               print*, "   * Kinetic energy"
               print*, "   * Kinetic energy rate"
               print*, "   * Kinetic energy balance"
               print*, "   * Enstrophy"
               print*, "   * Entropy"
               print*, "   * Entropy rate"
               print*, "   * Entropy balance"
               print*, "   * Math entropy"
               print*, "   * Artificial dissipation"
               print*, "   * Internal energy"
               print*, "   * Relative kenetic energy"
               print*, "   * Relative total energy"
               print*, "   * Mean velocity"
               print*, "   * Velocity"
               print*, "   * Momentum"
               print*, "   * source"
               print*, "   * particles source"
               print*, "   * sensor range"
               error stop "error stopped."

            end if
         end select
#elif defined(INCNS)
         select case ( trim ( self % variable ) )
         case ("mass")
         case ("entropy")
         case ("kinetic energy rate")
         case ("entropy rate")
         case ("source")            ; self % num_of_vars = NCONS
         case default

            if ( len_trim (self % variable) .eq. 0 ) then
               print*, "Variable was not specified for volume monitor " , self % ID , "."
            else
               print*, 'Variable "',trim(self % variable),'" volume monitor ', self % ID, ' not implemented yet.'
               print*, "Options available are:"
               print*, "   * Mass"
               print*, "   * Entropy"
               print*, "   * Kinetic energy rate"
               print*, "   * Entropy rate"
               print*, "   * source"
               error stop "error stopped."

            end if
         end select

#elif defined(MULTIPHASE)
         select case ( trim ( self % variable ) )
         case ("entropy rate")
         case ("entropy balance")
         case ("phase2-area")
         case ("phase2-xcog")
         case ("phase2-xvel")
         case ("source")            ; self % num_of_vars = NCONS
         case default

            if ( len_trim (self % variable) .eq. 0 ) then
               print*, "Variable was not specified for volume monitor " , self % ID , "."
            else
               print*, 'Variable "',trim(self % variable),'" volume monitor ', self % ID, ' not implemented yet.'
               print*, "Options available are:"
               print*, "   * Entropy rate"
               print*, "   * Entropy balance"
               print*, "   * Phase2-Area"
               print*, "   * Phase2-xCoG"
               print*, "   * Phase2-xVel"
               print*, "   * source"
               error stop "error stopped."

            end if
         end select



#elif defined(CAHNHILLIARD)
         select case ( trim ( self % variable ) )
         case ("free energy")
         case default
            if ( len_trim (self % variable) .eq. 0 ) then
               print*, "Variable was not specified for volume monitor " , self % ID , "."
            else
               print*, 'Variable "',trim(self % variable),'" volume monitor ', self % ID, ' not implemented yet.'
               print*, "Options available are:"
               print*, "   * Free energy"
               error stop "error stopped."

            end if
         end select
#elif defined(ACOUSTIC)
         select case ( trim ( self % variable ) )
         case ("acoustic energy")
         case ("source")            ; self % num_of_vars = NCONS
         case default

            if ( len_trim (self % variable) .eq. 0 ) then
               print*, "Variable was not specified for volume monitor " , self % ID , "."
            else
               print*, 'Variable "',trim(self % variable),'" volume monitor ', self % ID, ' not implemented yet.'
               print*, "Options available are:"
               print*, "   * acoustic energy"
               print*, "   * source"
               error stop "error stopped."

            end if
         end select
#endif


         allocate ( self % values(self % num_of_vars,BUFFER_SIZE) )

!
!        Prepare the file in which the monitor is exported
!        -------------------------------------------------
         write( self % fileName , '(A,A,A,A)') trim(solution_file) , "." , trim(self % monitorName) , ".volume"
!
!        Create file
!        -----------
         if (FirstCall) then
            open ( newunit = fID , file = trim(self % fileName) , status = "unknown" , action = "write" )
!
!        Write the file headers
!        ----------------------
            write( fID , '(A20,A  )') "#Monitor name:      ", trim(self % monitorName)
            write( fID , '(A20,A  )') "#Selected variable: " , trim(self % variable)

            write( fID , * )
            select case ( trim(self % variable) )
               case("velocity")
                  write( fID , '(A10,2X,A24,3(2X,A24))') "#Iteration" , "Time" , 'mean-u', 'mean-v', 'mean-w'
               case("momentum")
                  write( fID , '(A10,2X,A24,3(2X,A24))') "#Iteration" , "Time" , 'mean-rhou', 'mean-rhov', 'mean-rhow'
               case("source")
                  write( fID , '(A10,2X,A24,5(2X,A24))') "#Iteration" , "Time" , 'S1', 'S2', 'S3', 'S4', 'S5'
               case("particles source")
                  write( fID , '(A10,2X,A24,5(2X,A24))') "#Iteration" , "Time" , 'pS1', 'pS2', 'pS3', 'pS4', 'pS5'
               case("sensor range")
                  write( fID , '(A10,2X,A24,2(2X,A24))') "#Iteration" , "Time" , 'min sensor', 'max sensor'
               case default
                  write( fID , '(A10,2X,A24,2X,A24)'   ) "#Iteration" , "Time" , trim(self % variable)
            end select
            close ( fID )
         end if

      end subroutine VolumeMonitor_Initialization

      subroutine VolumeMonitor_Update ( self , mesh , bufferPosition, time )
!
!        *******************************************************************
!           This subroutine updates the monitor value computing it from
!           the mesh. It is stored in the "bufferPosition" position of the
!           buffer.
!        *******************************************************************
!
         use VolumeIntegrals
         implicit none
         class   (  VolumeMonitor_t ) :: self
         class   (  HexMesh       )   :: mesh
         integer                      :: bufferPosition
         real(kind=rp), intent(in)    :: time

         self % bufferLine = bufferPosition
!
!        Compute the volume integral
!        ---------------------------
         select case ( trim(self % variable) )
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
         case ("kinetic energy")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, KINETIC_ENERGY) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("kinetic energy rate")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, KINETIC_ENERGY_RATE) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("kinetic energy balance")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, KINETIC_ENERGY_BALANCE) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("enstrophy")
            self % values(1,bufferPosition) = 0.5_RP * ScalarVolumeIntegral(mesh, ENSTROPHY) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("entropy")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, ENTROPY) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("entropy rate")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, ENTROPY_RATE) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("entropy balance")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, ENTROPY_BALANCE) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("math entropy")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, MATH_ENTROPY) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("artificial dissipation")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, ARTIFICIAL_DISSIPATION) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("internal energy")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, INTERNAL_ENERGY) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("relative kinetic energy")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, RELATIVE_KINETIC_ENERGY, time=time) / ScalarVolumeIntegral(mesh, VOLUME)
         
         case ("relative total energy")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, RELATIVE_TOTAL_ENERGY, time=time) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("mean velocity")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, VELOCITY) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("velocity")
            self % values(:,bufferPosition) = VectorVolumeIntegral(mesh, VELOCITY, self % num_of_vars) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("momentum")
            self % values(:,bufferPosition) = VectorVolumeIntegral(mesh, MOMENTUM, self % num_of_vars) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("source")
            self % values(:,bufferPosition) = VectorVolumeIntegral(mesh, SOURCE, self % num_of_vars) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("particles source")
            self % values(:,bufferPosition) = VectorVolumeIntegral(mesh, PSOURCE, self % num_of_vars) / ScalarVolumeIntegral(mesh, VOLUME)

         case ("sensor range")
            call GetSensorRange(mesh, self % values(1,bufferPosition), self % values(2,bufferPosition))

#elif defined(INCNS)
         case ("mass")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, MASS)

         case ("entropy")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, ENTROPY)

         case ("kinetic energy rate")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, KINETIC_ENERGY_RATE)

         case ("entropy rate")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, ENTROPY_RATE)

         case ("source")
            self % values(:,bufferPosition) = VectorVolumeIntegral(mesh, SOURCE, self % num_of_vars) / ScalarVolumeIntegral(mesh, VOLUME)

#elif defined(MULTIPHASE)
         case ("entropy rate")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, ENTROPY_RATE)

         case ("entropy balance")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, ENTROPY_BALANCE)

         case ("phase2-xcog")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, PHASE2_XCOG)

         case ("phase2-xvel")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, PHASE2_XVEL)

         case ("phase2-area")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, PHASE2_AREA)

         case ("source")
            self % values(:,bufferPosition) = VectorVolumeIntegral(mesh, SOURCE, self % num_of_vars) / ScalarVolumeIntegral(mesh, VOLUME)

#elif defined(CAHNHILLIARD)
         case ("free energy")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, FREE_ENERGY) / ScalarVolumeIntegral(mesh, VOLUME)

#elif defined(ACOUSTIC)
         case ("acoustic energy")
            self % values(1,bufferPosition) = ScalarVolumeIntegral(mesh, ACOUSTIC_ENERGY)

         case ("source")
            self % values(:,bufferPosition) = VectorVolumeIntegral(mesh, SOURCE, self % num_of_vars) / ScalarVolumeIntegral(mesh, VOLUME)

#endif
         end select


      end subroutine VolumeMonitor_Update

      subroutine VolumeMonitor_WriteLabel ( self )
!
!        *************************************************************
!              This subroutine writes the label for the volume
!           monitor, when invoked from the time integrator Display
!           procedure.
!        *************************************************************
!
         implicit none
         class(VolumeMonitor_t)             :: self

         select case ( trim(self % variable) )
            case("velocity")
               write(STD_OUT , '(3(3X,A10))' , advance = "no") 'mean-u', 'mean-v', 'mean-w'
            case("momentum")
               write(STD_OUT , '(3(3X,A10))' , advance = "no") 'mean-rhou', 'mean-rhov', 'mean-rhow'
            case("source")
               write(STD_OUT , '(5(3X,A10))' , advance = "no") 'S1', 'S2', 'S3', 'S4', 'S5'
            case("particles source")
               write(STD_OUT , '(5(3X,A10))' , advance = "no") 'pS1', 'pS2', 'pS3', 'pS4', 'pS5'
            case("sensor range")
               write(STD_OUT , '(2(3X,A10))' , advance = "no") 'min sensor', 'max sensor'
            case default
               write(STD_OUT , '(3X,A10)' , advance = "no") trim(self % monitorName(1 : MONITOR_LENGTH))
         end select
      end subroutine VolumeMonitor_WriteLabel
!
!        *************************************************************
!        WriteValue: This subroutine writes the monitor value for the time
!           integrator Display procedure.
!        *************************************************************
      subroutine VolumeMonitor_WriteValue ( self , bufferLine )
         implicit none
         !-arguments-----------------------------------------------
         class(VolumeMonitor_t)     :: self
         integer                    :: bufferLine
         !-local-variables-----------------------------------------
         integer :: i
         !---------------------------------------------------------

         do i=1, self % num_of_vars
            write(STD_OUT , '(1X,A,1X,ES10.3)' , advance = "no") "|" , self % values (i, bufferLine )
         end do

      end subroutine VolumeMonitor_WriteValue

      subroutine VolumeMonitor_WriteToFile ( self , iter , t , no_of_lines)
!
!        *************************************************************
!              This subroutine writes the buffer to the file.
!        *************************************************************
!
         implicit none
         class(VolumeMonitor_t) :: self
         integer                 :: iter(:)
         real(kind=RP)           :: t(:)
         integer                 :: no_of_lines
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: i
         integer                    :: fID
         character(len=LINE_LENGTH) :: fmt

         if ( MPI_Process % isRoot ) then
            open( newunit = fID , file = trim ( self % fileName ) , action = "write" , access = "append" , status = "old" )

            write(fmt,'(A,I0,A)') '(I10,2X,ES24.16,', size(self % values,1), '(2X,ES24.16))'
            do i = 1 , no_of_lines
               write( fID , fmt ) iter(i) , t(i) , self % values(:,i)

            end do

            close ( fID )
         end if

         if ( no_of_lines .ne. 0 ) self % values(:,1) = self % values(:,no_of_lines)

      end subroutine VolumeMonitor_WriteToFile
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
      function VolumeMonitor_GetLast(self) result(lastValues)
         implicit none
         class(VolumeMonitor_t), intent(in) :: self
         real(kind=RP)                      :: lastValues ( size(self % values,1) )

         lastValues(:) = self % values (:,self % bufferLine)

      end function VolumeMonitor_GetLast
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
      elemental subroutine VolumeMonitor_Destruct (self)
         implicit none
         class(VolumeMonitor_t), intent(inout) :: self

         deallocate (self % values)
      end subroutine VolumeMonitor_Destruct

      elemental subroutine VolumeMonitor_Assign (to, from)
         implicit none
         class(VolumeMonitor_t), intent(inout)  :: to
         type(VolumeMonitor_t) , intent(in)     :: from

         to % active       = from % active
         to % ID           = from % ID

         safedeallocate (to % values)
         allocate ( to % values( size(from % values,1) , size(from % values,2) ) )
         to % values       = from % values

         to % monitorName  = from % monitorName
         to % fileName     = from % fileName
         to % variable     = from % variable
      end subroutine VolumeMonitor_Assign
end module VolumeMonitorClass
