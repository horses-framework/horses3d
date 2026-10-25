#include "Includes.h"
module SurfaceMonitorClass
   use SMConstants
   use HexMeshClass
   use MonitorDefinitions
   use PhysicsStorage
   use MPI_Process_Info
   use FluidData
   use FileReadingUtilities, only: getRealArrayFromString
   implicit none


#if defined(NAVIERSTOKES) || defined(INCNS)  
   private
   public   SurfaceMonitor_t


!
!  ********************************
!  Surface monitor class definition
!  ********************************
!
   type SurfaceMonitor_t
      logical                         :: active
      logical                         :: isDimensionless, IBM = .false.
      integer                         :: ID
      real(kind=RP)                   :: direction(NDIM)
      integer                         :: marker
      real(kind=RP), allocatable      :: referenceSurface
      real(kind=RP), allocatable      :: values(:)
      real(kind=RP)                   :: dynamicPressure
      real(kind=RP)                   :: pGauge = 0.0_RP
      character(len=STR_LEN_MONITORS) :: monitorName
      character(len=STR_LEN_MONITORS) :: fileName
      character(len=STR_LEN_MONITORS) :: variable
      contains
         procedure   :: Initialization => SurfaceMonitor_Initialization
         procedure   :: Update         => SurfaceMonitor_Update
         procedure   :: WriteLabel     => SurfaceMonitor_WriteLabel
         procedure   :: WriteValues    => SurfaceMonitor_WriteValue
         procedure   :: WriteToFile    => SurfaceMonitor_WriteToFile
         procedure   :: destruct       => SurfaceMonitor_Destruct
         procedure   :: copy           => SurfaceMonitor_Assign
         generic     :: assignment(=)  => copy
   end type SurfaceMonitor_t

   contains
!
!/////////////////////////////////////////////////////////////////////////
!
!           SURFACE MONITOR PROCEDURES
!           --------------------------
!/////////////////////////////////////////////////////////////////////////
!
      subroutine SurfaceMonitor_Initialization( self , mesh , ID, solution_file , FirstCall)
!
!        *****************************************************************************
!              This subroutine initializes the surface monitor. The following
!           data is obtained from the case file:
!              -> Name: The monitor name (10 characters maximum)
!              -> Marker: The surface marker in which the monitor will be computed.
!              -> Variable: The variable to be monitorized.
!              -> Reference surface (optional): Reference surface for lift/drag coefficients
!              -> Direction (optional): Direction in which the forces are computed
!              -> Gauge pressure (optional): Integrate p - 1/(gamma M^2) instead of p
!        *****************************************************************************
!  
         use ParamfileRegions
         implicit none
         class(SurfaceMonitor_t) :: self
         class(HexMesh)          :: mesh
         integer                 :: ID
         character(len=*)        :: solution_file
         logical, intent(in)     :: FirstCall
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=STR_LEN_MONITORS)  :: in_label
         character(len=STR_LEN_MONITORS)  :: fileName
         character(len=STR_LEN_MONITORS)  :: paramFile
         character(len=STR_LEN_MONITORS)  :: directionName
         integer, allocatable             :: marker
         character(len=STR_LEN_MONITORS)  :: markerName
         integer                          :: pos, i, STLNum
         integer                          :: fID
         integer                          :: zoneID
         real(kind=RP)                    :: directionValue(NDIM)
         logical, allocatable             :: gaugePressure
!
!        Get monitor ID
!        --------------
         self % ID = ID
!
!        Search for the parameters in the case file
!        ------------------------------------------
         write(in_label , '(A,I0)') "#define surface monitor " , self % ID
         
         call get_command_argument(1, paramFile)
         call readValueInRegion ( trim ( paramFile )  , "name"              , self % monitorName      , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "marker"            , markerName              , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "variable"          , self % variable         , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "reference surface" , self % referenceSurface , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "direction"         , directionName        , in_label , "# end" ) 
         call readValueInRegion ( trim ( paramFile )  , "gauge pressure"    , gaugePressure        , in_label , "# end" ) 
!
!        Gauge pressure
!        --------------
         self % pGauge = 0.0_RP
         if ( allocated(gaugePressure) ) then
            if ( gaugePressure ) then
#if defined(NAVIERSTOKES)
               self % pGauge = 1.0_RP / dimensionless % gammaM2
#else
               print*, '"gauge pressure" is only available for compressible flows (surface monitor ', self % ID, ').'
               error stop "error stopped"
#endif
            end if
         end if
!
!        Enable the monitor
!        ------------------
         self % active = .true.
         allocate ( self % values(BUFFER_SIZE) )
!
!        Get the surface marker
!        ----------------------
         self % marker = -1
         if( mesh% IBM% active ) then                
            do STLNum = 1, mesh% IBM% NumOfSTL
               if( mesh% IBM% Integral(STLNum)% compute ) cycle
               if( trim(mesh% IBM% STLfilename(STLNum)) .eq. trim(markerName) ) then
                 if( .not. mesh% IBM% ComputeBandRegion ) then
                     write(*,'(A)') "Warning: for surface monitors with IBM, 'band region' must be set '.true.'"
                     error stop
                  end if
                  self% marker = STLNum    
                  self% IBM    = .true.
                  mesh% IBM% Integral(STLNum)% compute = .true.
                  mesh% IBM% Integral(STLNum)% ListComputed = .false.
                  if( MPI_Process% isRoot ) then 
                     call mesh% IBM% stlSurfaceIntegrals(STLNum)% ReadTessellation( mesh% IBM% STLfilename(STLNum) )
                     if( mesh% IBM% ClipAxis .ne. 0 ) call mesh% IBM% stlSurfaceIntegrals(STLNum)% Clip( mesh% IBM% minCOORDS, mesh% IBM% maxCOORDS, mesh% IBM% ClipAxis, .false. )
                     call mesh% IBM% stlSurfaceIntegrals(STLNum)% SetIntegrationPoints()
                     call mesh% IBM% stlSurfaceIntegrals(STLNum)% SetIntegration( mesh% IBM% NumOfInterPoints )                        
                  end if
                  call mesh% IBM% SetIntegration( STLNum )                  
                  exit
               end if
            end do
         else
            do zoneID = 1, size(mesh % zones)
               if ( trim(mesh % zones(zoneID) % name) .eq. trim(markerName) ) then            
                  self % marker = zoneID
                  exit
               end if
            end do
         end if

         if ( self % marker .eq. -1 ) then
            self % active = .false.
            write(*,'(A,I0)') "Warning: Marker not specified for surface monitor ", self % ID
            write(*,'(A,I0,A)') "     Surface monitor ", self % ID, " disabled."
         end if         
!
!        Select the variable from the available list, and compute auxiliary variables if needed
!        --------------------------------------------------------------------------------------
!
!        ****************************************
         select case ( trim ( self % variable ) )
!        ****************************************
! 
            case ("mass-flow")
               self % isDimensionless = .false.

            case ("flow")
               self % isDimensionless = .false.

            case ("pressure-force")
               self % isDimensionless = .false.
               if ( len_trim(directionName) .eq. 0 ) then
                  print*, "Direction not specified for pressure-force in surface monitor " , self % ID , "."
                  error stop "error stopped"

               else
                  directionValue = getRealArrayFromString(directionName)
                  if ( size(directionValue) .ne. 3 ) then
                     print*, "Incorrect direction for monitor ", self % ID, "."
   
                  else
                     self % direction = directionValue   

                  end if
               end if

            case ("viscous-force")
               self % isDimensionless = .false.
               if ( len_trim(directionName) .eq. 0 ) then
                  print*, "Direction not specified for pressure-force in surface monitor " , self % ID , "."
                  error stop "error stopped"

               else
                  directionValue = getRealArrayFromString(directionName)
                  if ( size(directionValue) .ne. 3 ) then
                     print*, "Incorrect direction for monitor ", self % ID, "."
   
                  else
                     self % direction = directionValue   

                  end if
               end if

            case ("force")
               self % isDimensionless = .false.

               if ( len_trim(directionName) .eq. 0 ) then
                  print*, "Direction not specified for pressure-force in surface monitor " , self % ID , "."
                  error stop "error stopped"

               else
                  directionValue = getRealArrayFromString(directionName)
                  if ( size(directionValue) .ne. 3 ) then
                     print*, "Incorrect direction for monitor ", self % ID, "."
   
                  else
                     self % direction = directionValue   

                  end if
               end if


            case ("lift")
               self % isDimensionless = .true.

               if ( .not. allocated ( self % referenceSurface ) ) then
                  print*, "Reference surface not specified for lift surface monitor " , self % ID , "."
                  error stop "error stopped"
               end if
               
               if ( len_trim(directionName) .eq. 0 ) then
                  print*, "Direction not specified for lift in surface monitor " , self % ID , "."
                  print*, "    ...  Using [0,1,0] as default."
                  self % direction = [0._RP,1._RP,0._RP]

               else
                  directionValue = getRealArrayFromString(directionName)
                  if ( size(directionValue) .ne. 3 ) then
                     print*, "Incorrect direction for monitor ", self % ID, "."
   
                  else
                     self % direction = directionValue   

                  end if
               end if

               self % dynamicPressure = 0.5_RP * refValues % rho * POW2(refValues % V)* self % referenceSurface

            case ("drag")
               self % isDimensionless = .true.

               if ( .not. allocated ( self % referenceSurface ) ) then
                  print*, "Reference surface not specified for drag surface monitor " , self % ID , "."
                  error stop "error stopped"
               end if
               
               if ( len_trim(directionName) .eq. 0 ) then
                  print*, "Direction not specified for drag in surface monitor " , self % ID , "."
                  print*, "    ...  Using [1,0,0] as default."
                  self % direction = [1._RP,0._RP,0._RP]

               else
                  directionValue = getRealArrayFromString(directionName)
                  if ( size(directionValue) .ne. 3 ) then
                     print*, "Incorrect direction for monitor ", self % ID, "."
   
                  else
                     self % direction = directionValue   

                  end if
               end if

               self % dynamicPressure = 0.5_RP * refValues % rho * refValues % V * refValues % V * self % referenceSurface

            case ("pressure-average")
               self % isDimensionless = .false.

            case default

               if ( len_trim (self % variable) .eq. 0 ) then
                  print*, "Variable was not specified for surface monitor " , self % ID , "."
               else
                  print*, 'Variable "',trim(self % variable),'" surface monitor ', self % ID, ' not implemented yet.'
                  print*, "Options available are:"
                  print*, "   * mass-flow"
                  print*, "   * flow"
                  print*, "   * pressure-force"
                  print*, "   * viscous-force"
                  print*, "   * force"
                  print*, "   * lift"
                  print*, "   * drag"
                  print*, "   * pressure-average"
                  error stop "error stopped."

               end if
!
!        **********
         end select
!        **********
!
!        Prepare the file in which the monitor is exported
!        -------------------------------------------------
         write( self % fileName , '(A,A,A,A)') trim(solution_file) , "." , trim(self % monitorName) , ".surface"  
!
!        Create file
!        -----------
         if (FirstCall) then
            open ( newunit = fID , file = trim(self % fileName) , status = "unknown" , action = "write" ) 
!
!        Write the file headers
!        ----------------------
            write( fID , '(A20,A  )') "Monitor name:      ", trim(self % monitorName)
            write( fID , '(A20,I0 )') "Surface marker:    ", self % marker
            write( fID , '(A20,A  )') "Selected variable: " , trim(self % variable)

            if ( self % isDimensionless ) then
               write(fID , '(A20,ES24.10)') "Dynamic pressure: " , self % dynamicPressure
            end if

            if ( self % pGauge .ne. 0.0_RP ) then
               write(fID , '(A20,ES24.10)') "Gauge pressure:   " , self % pGauge
            end if

            write( fID , * )
            write( fID , '(A10,2X,A24,2X,A24)' ) "Iteration" , "Time" , trim(self % variable)

            close ( fID )
         end if
      end subroutine SurfaceMonitor_Initialization

      subroutine SurfaceMonitor_Update ( self, mesh, bufferPosition, iter, autosave, dt )
!
!        *******************************************************************
!           This subroutine updates the monitor value computing it from
!           the mesh. It is stored in the "bufferPosition" position of the 
!           buffer.
!        *******************************************************************
!
         use SurfaceIntegrals
         use IBMClass
         implicit none
         class   (  SurfaceMonitor_t )   :: self
         class   (  HexMesh       )      :: mesh
         integer                         :: bufferPosition, iter, STLNum
         real(kind=RP)                   :: F(NDIM)
         real(kind=RP)                   :: dt
         logical                         :: autosave
!
!        Pressure subtracted in the pressure integrals of this monitor
!        -------------------------------------------------------------
         pGauge = self % pGauge
         
         select case ( trim ( self % variable ) )

         case ("mass-flow")
            if( self% IBM ) then
               STLNum = self% marker
               call ScalarDataReconstruction( mesh% IBM, mesh% elements, STLNum, MASS_FLOW, iter, autosave, dt )
               self % values(bufferPosition) = mesh% IBM% stlSurfaceIntegrals(STLNum)% ComputeScalarIntegral()
            else
               self % values(bufferPosition) = ScalarSurfaceIntegral(mesh, self % marker, MASS_FLOW, iter)
            end if 
         case ("flow")
            if( self% IBM ) then
               STLNum = self% marker
               call ScalarDataReconstruction( mesh% IBM, mesh% elements, STLNum, FLOW_RATE, iter, autosave, dt ) 
               self % values(bufferPosition) = mesh% IBM% stlSurfaceIntegrals(STLNum)% ComputeScalarIntegral()
            else
               self % values(bufferPosition) = ScalarSurfaceIntegral(mesh, self % marker, FLOW_RATE, iter)
            end if 

         case ("pressure-force")
            if( self% IBM ) then 
               STLNum = self% marker
               call VectorDataReconstruction( mesh% IBM, mesh% elements, STLNum, PRESSURE_FORCE, iter, autosave, dt )
               F = mesh% IBM% stlSurfaceIntegrals(STLNum)% ComputeVectorIntegral()
            else
               F = VectorSurfaceIntegral(mesh, self % marker, PRESSURE_FORCE, iter)
            end if
            F = refValues % rho * POW2(refValues % V) * POW2(Lref) * F
            self % values(bufferPosition) = dot_product(F, self % direction)

         case ("viscous-force")
            if( self% IBM ) then 
               STLNum = self% marker
               call VectorDataReconstruction( mesh% IBM, mesh% elements, STLNum, VISCOUS_FORCE, iter, autosave, dt )
               F = mesh% IBM% stlSurfaceIntegrals(STLNum)% ComputeVectorIntegral()
            else
               F = VectorSurfaceIntegral(mesh, self % marker, VISCOUS_FORCE, iter)
            end if
            F = refValues % rho * POW2(refValues % V) * POW2(Lref) * F
            self % values(bufferPosition) = dot_product(F, self % direction)

         case ("force")
            if( self% IBM ) then 
               STLNum = self% marker
               call VectorDataReconstruction( mesh% IBM, mesh% elements, STLNum, TOTAL_FORCE, iter, autosave, dt )
               F = mesh% IBM% stlSurfaceIntegrals(STLNum) % ComputeVectorIntegral()
            else
               F = VectorSurfaceIntegral(mesh, self % marker, TOTAL_FORCE, iter)
            end if 
            F = refValues % rho * POW2(refValues % V) * POW2(Lref) * F
            self % values(bufferPosition) = dot_product(F, self % direction)

         case ("lift")
            if( self% IBM ) then 
               STLNum = self% marker
               call VectorDataReconstruction( mesh% IBM, mesh% elements, STLNum, TOTAL_FORCE, iter, autosave, dt )
               F = mesh% IBM% stlSurfaceIntegrals(STLNum) % ComputeVectorIntegral()
            else
               F = VectorSurfaceIntegral(mesh, self % marker, TOTAL_FORCE, iter)
            end if 
            F = 2.0_RP * POW2(Lref) * F / self % referenceSurface
            self % values(bufferPosition) = dot_product(F, self % direction)

#if defined (NAVIERSTOKES)
         case ("drag")
            if (flowIsNavierStokes) then
               if( self% IBM ) then 
                  STLNum = self% marker
                  call VectorDataReconstruction( mesh% IBM, mesh% elements, STLNum, TOTAL_FORCE, iter, autosave, dt )
                  F = mesh% IBM% stlSurfaceIntegrals(STLNum) % ComputeVectorIntegral()
               else
                  F = VectorSurfaceIntegral(mesh, self % marker, TOTAL_FORCE, iter)
               end if
            else
               if( self% IBM ) then 
                  STLNum = self% marker
                  call VectorDataReconstruction( mesh% IBM, mesh% elements,STLNum, TOTAL_FORCE, iter, autosave, dt )
                  F = mesh% IBM% stlSurfaceIntegrals(STLNum) % ComputeVectorIntegral()
               else
                  F = VectorSurfaceIntegral(mesh, self % marker, PRESSURE_FORCE, iter)
               end if 
            end if
            F = 2.0_RP * POW2(Lref) * F / self % referenceSurface
            self % values(bufferPosition) = dot_product(F, self % direction)
#endif
!TODO if true
#if defined (INCNS)
         case ("drag")
            if (.true.) then 
               if( self% IBM ) then 
                  STLNum = self% marker
                  call VectorDataReconstruction( mesh% IBM, mesh% elements, STLNum, TOTAL_FORCE, iter, autosave, dt )
                  F = mesh% IBM% stlSurfaceIntegrals(STLNum) % ComputeVectorIntegral()
               else
                  F = VectorSurfaceIntegral(mesh, self % marker, TOTAL_FORCE, iter)
               end if
            else
               if( self% IBM ) then 
                  STLNum = self% marker
                  call VectorDataReconstruction( mesh% IBM, mesh% elements,STLNum, TOTAL_FORCE, iter, autosave, dt )
                  F = mesh% IBM% stlSurfaceIntegrals(STLNum) % ComputeVectorIntegral()
               else
                  F = VectorSurfaceIntegral(mesh, self % marker, PRESSURE_FORCE, iter)
               end if 
            end if
            F = 2.0_RP * POW2(Lref) * F / self % referenceSurface
            self % values(bufferPosition) = dot_product(F, self % direction)
#endif
         case ("pressure-average")
            self % values(bufferPosition) = ScalarSurfaceIntegral(mesh, self % marker, PRESSURE_FORCE, iter) / ScalarSurfaceIntegral(mesh, self % marker, SURFACE, iter)
  
         end select
         

      end subroutine SurfaceMonitor_Update

      subroutine SurfaceMonitor_WriteLabel ( self )
!
!        *************************************************************
!              This subroutine writes the label for the surface
!           monitor, when invoked from the time integrator Display
!           procedure.
!        *************************************************************
!
         implicit none
         class(SurfaceMonitor_t), intent(in)   :: self

         write(STD_OUT , '(3X,A10)' , advance = "no") trim(self % monitorName(1 : MONITOR_LENGTH))

      end subroutine SurfaceMonitor_WriteLabel
   
      subroutine SurfaceMonitor_WriteValue ( self , bufferLine ) 
!
!        *************************************************************
!              This subroutine writes the monitor value for the time
!           integrator Display procedure.
!        *************************************************************
!
         implicit none
         class(SurfaceMonitor_t) :: self
         integer                 :: bufferLine

         write(STD_OUT , '(1X,A,1X,ES10.3)' , advance = "no") "|" , self % values ( bufferLine ) 

      end subroutine SurfaceMonitor_WriteValue 

      subroutine SurfaceMonitor_WriteToFile ( self , iter , t , no_of_lines)
!
!        *************************************************************
!              This subroutine writes the buffer to the file.
!        *************************************************************
!
         implicit none  
         class(SurfaceMonitor_t) :: self
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

         if ( MPI_Process % isRoot ) then

            open( newunit = fID , file = trim ( self % fileName ) , action = "write" , access = "append" , status = "old" )
         
            do i = 1 , no_of_lines
               write( fID , '(I10,2X,ES24.16,2X,ES24.16)' ) iter(i) , t(i) , self % values(i)
   
            end do
        
            close ( fID )
         end if
         
         if ( no_of_lines .ne. 0 ) self % values(1) = self % values(no_of_lines)
      
      end subroutine SurfaceMonitor_WriteToFile
!
      elemental subroutine SurfaceMonitor_Destruct (self)
         implicit none
         class(SurfaceMonitor_t), intent(inout) :: self
         
         safedeallocate (self % values)
         safedeallocate (self % referenceSurface)
      end subroutine SurfaceMonitor_Destruct
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      elemental subroutine SurfaceMonitor_Assign(to, from)
         implicit none
         class(SurfaceMonitor_t), intent(inout) :: to
         type(SurfaceMonitor_t),  intent(in)    :: from
         
         if ( from % active) then
            to % active          = from % active
            to % isDimensionless = from % isDimensionless
            to % ID              = from % ID
            to % direction       = from % direction
            to % marker          = from % marker
            
            safedeallocate(to % referenceSurface)
            if ( allocated(from % referenceSurface) ) then
               allocate (to % referenceSurface)
               to % referenceSurface = from % referenceSurface
            end if
            
            safedeallocate(to % values)
            allocate ( to % values( size(from % values) ) )
            to % values          = from % values
            
            to % dynamicPressure = from % dynamicPressure
            to % pGauge          = from % pGauge
            to % monitorName     = from % monitorName
            to % fileName        = from % fileName
            to % variable        = from % variable
         else
            to % active = .FALSE.
         end if
         
      end subroutine SurfaceMonitor_Assign
#endif      
end module SurfaceMonitorClass
