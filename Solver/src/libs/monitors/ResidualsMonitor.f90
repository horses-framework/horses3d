#include "Includes.h"
module ResidualsMonitorClass
   use SMConstants
   use PhysicsStorage
   use HexMeshClass
   use MonitorDefinitions
   use MPI_Process_Info
#ifdef _HAS_MPI_
      use mpi
#endif
   implicit none
   
   private
   public   Residuals_t
!
!  **************************
!  Residuals class definition
!  **************************
!
   type Residuals_t
      logical                         :: active
      real(kind=RP), allocatable      :: values(:,:)
      real(kind=RP), allocatable      :: CPUtime(:)
      character(len=STR_LEN_MONITORS) :: fileName
	  logical                         :: memory =.FALSE.
	  real(kind=RP), allocatable      :: totMemory(:)
      contains
         procedure   :: Initialization => Residuals_Initialization
         procedure   :: Update         => Residuals_Update
         procedure   :: WriteLabel     => Residuals_WriteLabel
         procedure   :: WriteValues    => Residuals_WriteValue
         procedure   :: WriteToFile    => Residuals_WriteToFile
         procedure   :: destruct       => Residuals_Destruct
         procedure   :: copy           => Residuals_Assign
         generic     :: assignment(=)  => copy
   end type Residuals_t
!
!  ========
   contains
!  ========
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!           RESIDUALS ROUTINES
!           ------------------
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Residuals_Initialization( self, solution_file, FirstCall, monitorMem ) 
!
!        *******************************************************************
!              This subroutine initializes the residuals structure
!        *******************************************************************
!
         implicit none
         class(Residuals_t)  :: self
         character(len=*)    :: solution_file
         logical, intent(in) :: FirstCall
		 logical, intent(in) :: monitorMem
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=STR_LEN_MONITORS)  :: fileName
         integer                          :: fID
         integer                          :: pos
!
!        Enable the monitor
!        ------------------
         self % active = .true.
         allocate ( self % values(NCONS,1:BUFFER_SIZE) , self % CPUtime(BUFFER_SIZE) )
!
!        Get monitor file name
!        ---------------------
         write( self % fileName , '(A,A)') trim(solution_file) , ".residuals"  
!
!        Memory Monitor
!        ---------------------
		 self % memory = monitorMem
		 if (self % memory) then
			allocate( self % totMemory(BUFFER_SIZE))
			self % totMemory = 0.0_RP
		 end if 
!
!        Create file to write the residuals
!        ----------------------------------
         if (FirstCall) then
            open ( newunit = fID , file = trim(self % fileName) , status = "unknown" , action = "write" ) 
            write ( fID , ' ( A                                      ) ' ) "#Residuals file"
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
            write ( fID , ' ( A10,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24 ) ' , advance = "no") "#Iteration" , "Time" , &
                        "Total elapsed Time (s)", "Solver elapsed Time (s)" , "continuity" , "x-momentum" , "y-momentum" , "z-momentum", "energy" , "Max-Residual"
#elif defined(SPALARTALMARAS)
            write ( fID , ' ( A10,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24 ) ' , advance = "no") "#Iteration" , "Time" , &
                        "Total elapsed Time (s)", "Solver elapsed Time (s)" , "continuity" , "x-momentum" , "y-momentum" , "z-momentum", "energy" ," Turbulence Param" , "Max-Residual"
#elif defined(INCNS)
            write ( fID , ' ( A10,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24 ) ' , advance = "no") "#Iteration" , "Time" , &
                        "Elapsed Time (s)" , "dens-transp" , "x-momentum" , "y-momentum" , "z-momentum", "div-v" , "Max-Residual"
#elif defined(MULTIPHASE)
            write ( fID , ' ( A10,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24 ) ' , advance = "no") "#Iteration" , "Time" , &
                        "Elapsed Time (s)", "Solver elapsed Time (s)" , "concentration" , "x-momentum" , "y-momentum" , "z-momentum", "div-v" , "Max-Residual"
#elif defined(CAHNHILLIARD)
            write ( fID , ' ( A10,2X,A24,2X,A24) ' , advance = "no") "#Iteration" , "Time" , "concentration"

#elif defined(ACOUSTIC)
            write ( fID , ' ( A10,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24,2X,A24 ) ', advance = "no" ) "#Iteration" , "Time" , &
                        "Elapsed Time (s)" , "density" , "x-velocity" , "y-velocity" , "z-velocity", "pressure" , "Max-Residual"

#endif
			if (self % memory) write ( fID , ' ( 2X,A24 ) ', advance = "no") "Total Memory(kB)" 
			write ( fID , '(2X)')
!
!        Close file
!        ----------
            close ( fID ) 
         end if
      end subroutine Residuals_Initialization

      subroutine Residuals_Update ( self, mesh, maxResiduals, bufferPosition)
!
!        *********************************************************
!              This subroutine updates the residuals values from
!           those computed in the Monitor procedure
!        *********************************************************
!
         implicit none
         class(Residuals_t)         :: self
         class(HexMesh), intent(in) :: mesh
         real(kind=RP)              :: maxResiduals(NCONS)
         integer                    :: bufferPosition
!
!        ---------------
!        Local variables
!        ---------------
!
		 integer 					:: fID, ios
		 integer                    :: ierr 
		 real(kind=RP) 				:: mem_kb
		 character(len=256) 		:: line
!
!        Update buffer values
!        --------------------      
         self % values( 1:NCONS, bufferPosition ) = maxResiduals
!
!        Get total memory
!        --------------------   
		 if (.not. self % memory) return
		 
		! Default to 0 in case parsing fails
		 mem_kb = 0.0d0

		! Open /proc/self/status and read VmRSS
		open(newunit=fID, file="/proc/self/status", status="old", action="read", iostat=ios)
		if (ios == 0) then
			do
				read(fID, '(A)', iostat=ios) line
				if (ios /= 0) exit
				if (line(1:6) == 'VmRSS:') then
					read(line(7:), *, iostat=ios) mem_kb
					exit
				end if
			end do
			close(fID)
		end if
		self % totMemory(bufferPosition) = mem_kb
#ifdef _HAS_MPI_
		if (MPI_Process % doMPIAction) then
			! Reduce memory usage across all ranks
			call MPI_Allreduce(mem_kb, self % totMemory(bufferPosition), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
		end if
#endif

      end subroutine Residuals_Update

      subroutine Residuals_WriteLabel ( self )
!
!        ************************************************************
!              This subroutine displays the residuals labels for the
!           time integrator Display procedure.
!        ************************************************************
!
         implicit none
         class(Residuals_t)             :: self
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
         write(STD_OUT , '(3X,A10)' , advance = "no") "continuity"
         write(STD_OUT , '(3X,A10)' , advance = "no") "x-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "y-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "z-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "energy"
#ifdef TRANSPORT
         write(STD_OUT , '(3X,A10)' , advance = "no") "concentration"
#endif
#elif defined(SPALARTALMARAS)
         write(STD_OUT , '(3X,A10)' , advance = "no") "continuity"
         write(STD_OUT , '(3X,A10)' , advance = "no") "x-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "y-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "z-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "energy"
         write(STD_OUT , '(3X,A10)' , advance = "no") "turbulence parameter"
#elif defined(INCNS)
         write(STD_OUT , '(3X,A10)' , advance = "no") "dens-transp"
         write(STD_OUT , '(3X,A10)' , advance = "no") "x-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "y-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "z-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "div-v"
#elif defined(MULTIPHASE)
         write(STD_OUT , '(3X,A10)' , advance = "no") "concentration"
         write(STD_OUT , '(3X,A10)' , advance = "no") "x-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "y-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "z-momentum"
         write(STD_OUT , '(3X,A10)' , advance = "no") "div-v"
#elif defined(CAHNHILLIARD)
         write(STD_OUT , '(3X,A10)' , advance = "no") "concentration"
#elif defined(ACOUSTIC)
         write(STD_OUT , '(3X,A10)' , advance = "no") "density"
         write(STD_OUT , '(3X,A10)' , advance = "no") "x-velocity"
         write(STD_OUT , '(3X,A10)' , advance = "no") "y-velocity"
         write(STD_OUT , '(3X,A10)' , advance = "no") "z-velocity"
         write(STD_OUT , '(3X,A10)' , advance = "no") "pressure"
#endif
		if (self % memory) then 
			write(STD_OUT , '(3X,A10)' , advance = "no") "memory(kB)"
		end if 

      end subroutine Residuals_WriteLabel
   
      subroutine Residuals_WriteValue ( self , bufferLine ) 
!
!        ***************************************************************
!              This subroutine displays the residuals values for the 
!           time integrator Display procedure
!        ***************************************************************
!
         implicit none
         class(Residuals_t) :: self
         integer            :: bufferLine
!        ---------------------------------------------------------
         integer            :: eq
      
         do eq = 1 , NCONS
            write(STD_OUT , '(1X,A,1X,ES10.3)' , advance = "no") "|" , self % values(eq , bufferLine)
         end do
		if (self % memory) then 
			write(STD_OUT , '(1X,A,1X,ES10.3)' , advance = "no") "|" , self % totMemory(bufferLine)
		end if 

      end subroutine Residuals_WriteValue 

      subroutine Residuals_WriteToFile ( self , iter , t, TotalSimuTime, SolverSimuTime , no_of_lines)
!
!        *********************************************************************
!              This subroutine exports the results to the monitor file.
!           Just "no_of_lines" buffer lines are written.
!        *********************************************************************
!
         implicit none  
         class(Residuals_t)             :: self
         integer                    :: iter(:)
         real(kind=RP)              :: t(:)
         real(kind=RP)              :: TotalSimuTime(:)
         real(kind=RP)              :: SolverSimuTime(:)
         integer                    :: no_of_lines
!        -------------------------------------------
         integer                    :: i
         integer                    :: fID

#ifdef _HPC_MODE
   return
#endif          
         if ( MPI_Process % isRoot ) then
!
!           Open file
!           ---------
            open( newunit = fID , file = trim ( self % fileName ) , action = "write" , access = "append" , status = "old" )
!
!           Write values
!           ------------
            do i = 1 , no_of_lines
               write(fID, '(I10,3(2X,ES24.16))', advance="no") iter(i), t(i), TotalSimuTime(i), SolverSimuTime(i)
               
			   if (self % memory) then 
					write(fID, 111, advance = "no") self % values(1:NCONS,i), maxval(self % values(1:NCONS,i))
					write(fID , '(2X,ES16.8)') self % totMemory(i)
			   else
					write(fID, 111) self % values(1:NCONS,i), maxval(self % values(1:NCONS,i))
			   end if 
            end do
!
!           Close file
!           ----------
            close ( fID )

         end if

         if ( no_of_lines .ne. 0 ) then
            self % values(1:NCONS,1) = self % values(1:NCONS,no_of_lines)
			if (self % memory) then 
				self % totMemory(1)      = self % totMemory(no_of_lines)
			end if
         end if

#if defined(FLOW)
111 format(6(2X,ES24.16))
#elif (!defined(NAVIERSTOKES) && defined(CAHNHILLIARD))
111 format(2(2X,ES24.16))
#endif
      end subroutine Residuals_WriteToFile
      
      pure subroutine Residuals_Destruct (self)
         implicit none
         class(Residuals_t), intent(inout) :: self
         
         safedeallocate (self % values)
         safedeallocate (self % CPUtime)
		 safedeallocate (self % totMemory)
         self % fileName = ""
         self % active   = .FALSE.
		 self % memory   = .FALSE.
         
      end subroutine Residuals_Destruct
      
      elemental subroutine Residuals_Assign(to, from)
         implicit none
         class(Residuals_t), intent(inout)   :: to
         type(Residuals_t) , intent(in)      :: from
         
         to % active = from % active
         
         safedeallocate (to % values)
         allocate ( to % values( size(from % values,1) , size(from % values,2) ) )
         to % values = from % values
         
         safedeallocate (to % CPUtime)
         allocate ( to % CPUtime( size(from % CPUtime) ) )
         to % CPUtime = from % CPUtime
		 
		 safedeallocate (to % totMemory)
		 if (allocated(from % totMemory)) then
			allocate ( to % totMemory( size(from % totMemory) ) )
			to % totMemory = from % totMemory
	     end if 
         
         to % fileName = from % fileName
		 
		 to % memory = from % memory
         
      end subroutine Residuals_Assign
end module ResidualsMonitorClass
