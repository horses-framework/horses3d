program horses2plt
   use SMConstants
   use getTask
   use Mesh2PltModule
   use Solution2PltModule
   use Solution2VtkHdfModule
   use Solution2FoamModule
   use Stats2PltModule
   use SolutionFile
   use SharedSpectralBasis
   use Headers
   use MPI_Process_Info
   implicit none
   integer                                 :: jobType
   character(len=LINE_LENGTH)              :: meshName
   integer                                 :: no_of_solutions
   character(len=LINE_LENGTH), allocatable :: solutionNames(:)
   integer, allocatable                    :: solutionTypes(:)
   logical                                 :: fixedOrder
   integer                                 :: Nout(3)
   integer                                 :: basis
   integer                                 :: mode
   integer                                 :: iSol
   logical                                 :: useCommandArgs
   logical                                 :: oldStats
   logical 								   :: writeMesh

   call MPI_Process % Init

   call Main_Header("HORSES to TecPlot conversion utility",__DATE__,__TIME__)
!
!  Get the job type
!  ----------------
   call getTaskType(jobType, meshName, no_of_solutions, solutionNames, solutionTypes, fixedOrder, Nout, basis, mode, useCommandArgs, oldStats, writeMesh)
!
!  Construct Spectral basis
!  ------------------------
   call ConstructSpectralBasis
!
!  Perform the conversion to tecplot
!  ---------------------------------
   select case (jobType)
   case (MESH_2_PLT)
      write(STD_OUT,'(/,/)')
      call Section_Header("Mesh conversion")
      call Mesh2Plt(meshName, mode)

   case (SOLUTION_2_PLT)
      do iSol = 1, no_of_solutions

         select case (solutionTypes(iSol))
         case ( SOLUTION_FILE, ZONE_SOLUTION_FILE, ZONE_SOLUTION_AND_DOT_FILE, SOLUTION_AND_SENSOR_FILE )
            write(STD_OUT,'(/,/)')
            call Section_Header("Solution conversion")
            call Solution2Plt(meshName, solutionNames(iSol), fixedOrder, basis, Nout, mode)        

         case ( SOLUTION_AND_GRADIENTS_FILE, SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE)
            write(STD_OUT,'(/,/)')
            call Section_Header("Solution with gradients conversion")
            call Solution2Plt(meshName, solutionNames(iSol), fixedOrder, basis, Nout, mode)        

         case ( STATS_FILE )
            write(STD_OUT,'(/,/)')
            if (oldStats) then
                call Section_Header("Statistics legacy file conversion")
                call Stats2Plt(meshName, solutionNames(iSol), fixedOrder, basis, Nout)
            else
                call Section_Header("Statistics file conversion")
                call Solution2Plt(meshName, solutionNames(iSol), fixedOrder, basis, Nout, mode)        
            end if 

         end select
      end do

   case (SOLUTION_2_VTKHDF)
#ifdef HAS_HDF5
      do iSol = 1, no_of_solutions

         select case (solutionTypes(iSol))
         case ( SOLUTION_FILE, ZONE_SOLUTION_FILE, ZONE_SOLUTION_AND_DOT_FILE, SOLUTION_AND_SENSOR_FILE )
            write(STD_OUT,'(/,/)')
            call Section_Header("Solution conversion")
            call Solution2VtkHdf(meshName, solutionNames(iSol), Nout)

         case ( SOLUTION_AND_GRADIENTS_FILE, SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE)
            write(STD_OUT,'(/,/)')
            call Section_Header("Solution with gradients conversion")
            call Solution2VtkHdf(meshName, solutionNames(iSol), Nout)

         case ( STATS_FILE )
            write(STD_OUT,'(/,/)')
            if (oldStats) then
                write(STD_OUT,'(A)') "Statistics legacy file conversion not supported for VTKHDF"
                error stop
            else
               call Section_Header("Statistics file conversion")
               call Solution2VtkHdf(meshName, solutionNames(iSol), Nout)
            end if

         end select
      end do
#else
      write(STD_OUT, '(A)') "HDF5 support must be enabled to save to VTKHDF format"
      error stop
#endif

   case (SOLUTION_2_FOAM)
         call Solution2Foam (meshName, no_of_solutions, solutionNames, solutionTypes, Nout, writeMesh)  

   case (UNKNOWN_JOB)
      call exit(UNKNOWN_JOB)

   end select
!
!  Destruct Spectral basis
!  -----------------------
   call DestructSpectralBasis
   
   call MPI_Process % Close

end program horses2plt