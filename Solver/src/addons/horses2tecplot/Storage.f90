#include "Includes.h"
module Storage
   use SMConstants
   use PhysicsStorage
   use SolutionFile
   use VariableConversion
   implicit none

   private
   public Mesh_t, Element_t, Boundary_t
   public NVARS, NGRADVARS, hasMPIranks, hasBoundaries, isOldStats
   public partitionFileName, boundaryFileName, flowEq
   public hasExtraGradients, hasMu_NS, hasUt_NS, hasWallY, NSTAT, hasMu_sgs, hasSource

   integer                          :: NVARS, NGRADVARS
   logical                          :: hasMPIranks, hasBoundaries, isOldStats
   logical                          :: hasExtraGradients = .false.
   logical                          :: hasUt_NS = .false.
   logical                          :: hasMu_NS = .false.
   logical                          :: hasWallY     = .false.
   logical                          :: hasMu_sgs = .false.
   logical                          :: hasSource = .false.
   character(len=LINE_LENGTH)       :: boundaryFileName, partitionFileName, flowEq
   integer, parameter               :: NSTAT = 9

   type Element_t
!                                /* Mesh quantities */
      integer                    :: eID
      integer                    :: Nmesh(NDIM)
      integer                    :: mpi_rank = 0
      real(kind=RP), pointer     :: x(:,:,:,:)
!                                /* Solution quantities */
      integer                    :: Nsol(NDIM)
      real(kind=RP), pointer     :: Q(:,:,:,:)
      real(kind=RP), pointer     :: QDot(:,:,:,:)
      real(kind=RP), pointer     :: U_x(:,:,:,:)
      real(kind=RP), pointer     :: U_y(:,:,:,:)
      real(kind=RP), pointer     :: U_z(:,:,:,:)
      real(kind=RP), pointer     :: Q_x(:,:,:,:)
      real(kind=RP), pointer     :: Q_y(:,:,:,:)
      real(kind=RP), pointer     :: Q_z(:,:,:,:)
      real(kind=RP), pointer     :: mu_NS(:,:,:,:)
      real(kind=RP), pointer     :: ut_NS(:,:,:,:)
      real(kind=RP), pointer     :: wallY(:,:,:,:)
      real(kind=RP), pointer     :: mu_sgs(:,:,:,:)
      real(kind=RP), pointer     :: source(:,:,:,:)
      real(kind=RP), pointer     :: stats(:,:,:,:)
      real(kind=RP)              :: sensor
!                                /* Output quantities */
      integer                    :: Nout(NDIM)
      real(kind=RP), pointer     :: xOut(:,:,:,:)
      real(kind=RP), pointer     :: Qout(:,:,:,:)
      real(kind=RP), pointer     :: QDot_out(:,:,:,:)
      real(kind=RP), pointer     :: U_xout(:,:,:,:)
      real(kind=RP), pointer     :: U_yout(:,:,:,:)
      real(kind=RP), pointer     :: U_zout(:,:,:,:)
      real(kind=RP), pointer     :: mu_NSout(:,:,:,:)
      real(kind=RP), pointer     :: ut_NSout(:,:,:,:)
      real(kind=RP), pointer     :: wallYout(:,:,:,:)
      real(kind=RP), pointer     :: mu_sgsout(:,:,:,:)
      real(kind=RP), pointer     :: sourceout(:,:,:,:)
      real(kind=RP), pointer     :: statsout(:,:,:,:)

      real(kind=RP), allocatable :: outputVars(:,:,:,:)
   end type Element_t

   type Boundary_t
      character(len=LINE_LENGTH) :: Name
      integer                    :: no_of_faces
      integer, allocatable       :: elements(:)
      integer, allocatable       :: elementSides(:)
	  logical, allocatable       :: cornerDomain(:)
	  logical, allocatable       :: edgeDomain(:)	
   end type Boundary_t

   type Mesh_t
      integer  :: no_of_elements
      integer  :: nodeType
	  real (kind=RP) 	:: time
      type(Element_t),   allocatable    :: elements(:)
      type(Boundary_t),  allocatable    :: boundaries(:)
      character(len=LINE_LENGTH) :: meshName
      character(len=LINE_LENGTH) :: solutionName
      real(kind=RP)              :: refs(NO_OF_SAVED_REFS)
      logical                    :: hasGradients
      logical                    :: hasSensor
      logical                    :: isSurface
      logical                    :: hasTimeDeriv
      logical                    :: isStatistics
      logical                    :: is2D
      contains
         procedure   :: ReadMesh     => Mesh_ReadMesh
         procedure   :: ReadSolution => Mesh_ReadSolution
   end type Mesh_t

   contains
      subroutine Mesh_ReadMesh(self,meshName)
         use Headers
         implicit none
         class(Mesh_t)         :: self
         character(len=*), intent(in)     :: meshName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer, dimension(:), allocatable     :: arrayDimensions
         integer                                :: fileType, dimensionsSize
         integer                                :: fid, eID, i, pos
         character(len=1024)                    :: msg

         self % meshName = trim(meshName)
!
!        Get mesh node type and type of file
!        -----------------------------------
         self % nodeType = getSolutionFileNodeType(meshName)
         fileType = getSolutionFileType(meshName)
!
!        Get number of elements
!        ----------------------
         self % no_of_elements = getSolutionFileNoOfElements(meshName)
!
!        Allocate elements
!        -----------------
         allocate(self % elements(self % no_of_elements))
!
!        Allocate dimension based on type of file
!        ----------------------------------------
         select case (fileType)
         case (ZONE_MESH_FILE)
             dimensionsSize = 3
         case default
             dimensionsSize = 4
         end select

         allocate(arrayDimensions(dimensionsSize))
!
!        Read coordinates
!        ----------------
         fid = putSolutionFileInReadDataMode(meshName)

         self % is2D = .false.
         do eID = 1, self % no_of_elements
            associate ( e => self % elements(eID) )
            e % eID = eID

            call getSolutionFileArrayDimensions(fid,arrayDimensions)
!
!           Allocate memory for the coordinates
!           -----------------------------------
            ! e % Nmesh(1:3) = arrayDimensions(2:4) - 1
            e % Nmesh(1:dimensionsSize-1) = arrayDimensions(2:dimensionsSize) - 1
!           Use a 0 index for the case of a surface mesh
            if (dimensionsSize .eq. 3) then
               e % Nmesh(3) = 0
            end if
            if (e % Nmesh(3) .eq. 0) then
               self % is2D = .true.
            end if
            allocate( e % x(NDIM,0:e % Nmesh(1),0:e % Nmesh(2),0:e % Nmesh(3)) )
!
!           Read data
!           ---------
            read(fid) e % x

            end associate
         end do
!
!        Close file
!        ----------
         close(fid)

!        Read boundary file (if present)
!        -------------------------------
         if (hasBoundaries) then
               call readBoundaryFile(self % boundaries)
         end if
!
!        Describe the mesh
!        -----------------
		 write(STD_OUT,'(/)')
         write(STD_OUT,'(10X,A,A)') "Loading Mesh File:"
         write(STD_OUT,'(10X,A,A)') "-----------------"
         write(msg,'(A,A,A)') 'Mesh file "',trim(meshName),'":'
         call SubSection_Header(trim(msg))

         write(STD_OUT,'(30X,A,A30,I0)') "->", "Number of elements: ", self % no_of_elements
         select case ( self % nodeType )
         case(1)
            write(STD_OUT,'(30X,A,A30,A)') "->","Discretization nodes: ","Gauss"
         case(2)
            write(STD_OUT,'(30X,A,A30,A)') "->","Discretization nodes: ","Gauss-Lobatto"
         end select

!
!        Describe the boundaries
!        -----------------
         if (hasBoundaries) then
             write(msg,'(A,A,A)') 'Boundary file "',trim(boundaryFileName),'":'
             write(STD_OUT,'(/)')
             call SubSection_Header(trim(msg))
             write(STD_OUT,'(30X,A,A30,I0)') "->","Number of Boundaries: ", size(self % boundaries)

         end if

      end subroutine Mesh_ReadMesh

      subroutine Mesh_ReadSolution(self,solutionName)
         use Headers
         implicit none
         class(Mesh_t)         :: self
         character(len=*), intent(in)     :: solutionName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                        :: no_of_elements
         integer, allocatable           :: arrayDimensions(:)
         integer                        :: fileType, dimensionsSize
         integer                        :: fid, eID, pos
         integer                        :: i,j,k
         integer                        :: iter
         real(kind=RP)                  :: time
         character(len=1024)  :: msg

         self % solutionName = trim(solutionName)
		 write(STD_OUT,'(10X,A,A)') "Loading Solution File:"
         write(STD_OUT,'(10X,A,A)') "---------------------"
		 write(STD_OUT,'(30X,A,A30,A)') "->","Solution File: ", trim(solutionName)													  
!
!        Get the solution file type
!        --------------------------
         self % hasTimeDeriv = .false.
         self % isStatistics = .false.
         select case ( getSolutionFileType(trim(solutionName)) )

         case (SOLUTION_FILE)
            dimensionsSize = 4
            self % hasGradients = .false.
            self % hasSensor    = .false.

         case (SOLUTION_AND_GRADIENTS_FILE)
            dimensionsSize = 4
            self % hasGradients = .true.
            self % hasSensor    = .false.

         case (STATS_FILE)
            dimensionsSize = 4
            self % hasGradients = .false.
            self % hasSensor    = .false.
            self % isStatistics = .true.

         case (ZONE_SOLUTION_FILE)
            dimensionsSize = 3
            self % hasGradients = .false.
            self % hasSensor    = .false.

         case (ZONE_SOLUTION_AND_DOT_FILE)
            dimensionsSize = 3
            self % hasGradients = .false.
            self % hasSensor    = .false.
            self % hasTimeDeriv = .true.

         case (SOLUTION_AND_SENSOR_FILE)
            dimensionsSize = 4
            self % hasGradients = .false.
            self % hasSensor    = .true.

         case (SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE)
            dimensionsSize = 4
            self % hasGradients = .true.
            self % hasSensor    = .true.

         case default
            print*, "File expected to be a solution file"
            errorMessage(STD_OUT)
            error stop
         end select

         self % isSurface = (dimensionsSize .eq. 3)

         self % hasGradients = self % hasGradients .or. hasExtraGradients
!
!        Get node type
!        -------------
         if ( getSolutionFileNodeType(solutionName) .ne. self % nodeType ) then
            print*, "Solution and Mesh node type differs"
            errorMessage(STD_OUT)
            error stop
         end if
!
!        Get number of elements
!        ----------------------
         no_of_elements = getSolutionFileNoOfElements(solutionName)
         if ( self % no_of_elements .ne. no_of_elements ) then
            write(STD_OUT,'(30X,A,I0,A,I0,A)') "The number of elements in the mesh (",self % no_of_elements,&
                                           ") differs to that of the solution (",no_of_elements,")."
            errorMessage(STD_OUT)
            error stop
         end if
         allocate(arrayDimensions(dimensionsSize))
!
!        Get time and iteration
!        ----------------------
         call getSolutionFileTimeAndIteration(trim(solutionName),iter,self % time)
!
!        Read reference values
!        ---------------------
         self % refs = getSolutionFileReferenceValues(trim(solutionName))
!
!        Read coordinates
!        ----------------
         fid = putSolutionFileInReadDataMode(solutionName)

         ! call set_getVelocityGradients(GRADVARS_STATE) ! FIXME: MIGHT BE NEEDED FOR HORSES2PLT
         ! write(STD_OUT,'(15X,A)') " WARNING horses2tecplot.90 :: Velocity Gradients set to default (GRADVARS_STATE)"
      
         if ( .not. isOldStats ) then
         ! if ( .not. self % isStatistics ) then
            do eID = 1, self % no_of_elements
               associate ( e => self % elements(eID) )
               call getSolutionFileArrayDimensions(fid,arrayDimensions)

               call getNVARS(arrayDimensions(1), self % isStatistics)
!   
               ! e % Nsol(1:3) = arrayDimensions(2:4) - 1
               e % Nsol(1:dimensionsSize-1) = arrayDimensions(2:dimensionsSize) - 1
               if (dimensionsSize .eq. 3) e % Nsol(3) = 0
!   
!
!              Allocate memory for the statistics
!              ----------------------------------
               if (self % isStatistics) then
                   allocate( e % stats(1:NSTAT,0:e % Nsol(1), 0:e % Nsol(2), 0:e % Nsol(3)) )
                   read(fid) e % stats
!
               end if
!
!              Allocate memory for the coordinates
!              -----------------------------------            
               allocate( e % Q(1:NVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
!
!              Read data
!              ---------
               read(fid) e % Q

              ! Qdot goes before gradients when present
               if (self % hasTimeDeriv) then
                   allocate( e % Qdot(1:NVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
                   read(fid) e% Qdot
               end if

               if ( self % hasGradients ) then
!
!                 Allocate memory for the gradients
!                 ---------------------------------
                  allocate( e % Q_x(1:NVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) ) ! TODO: Allocate depending on physics
                  allocate( e % Q_y(1:NVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
                  allocate( e % Q_z(1:NVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )

                  allocate( e % U_x(1:NGRADVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
                  allocate( e % U_y(1:NGRADVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
                  allocate( e % U_z(1:NGRADVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
!
!                 Read data
!                 ---------
                  read(fid) e % Q_x
                  read(fid) e % Q_y
                  read(fid) e % Q_z

!                 Call set_getVelocityGradients to make the pointer to the actual subroutine, is needed only for the NS
!                 Set state as is the default option TODO point to the correct one if its possible (oscar note)
!                 ---------------------------
                  call set_getVelocityGradients(GRADVARS_STATE)

                  ! Following block works for NS, CH, NSCH and iNS .... but not iNSCH: change 5 by 6 to use iNSCH (NS won't work)
                  if (NVARS .ge. 5) then
                     do k = 0,e % Nsol(3) ; do j = 0, e % Nsol(2) ; do i = 0, e % Nsol(1)
                        call getVelocityGradients(e % Q(:,i,j,k), e % Q_x(1:5,i,j,k), e % Q_y(1:5,i,j,k), e % Q_z(1:5,i,j,k), &
                                                  e % U_x(1:3,i,j,k), e % U_y(1:3,i,j,k), e % U_z(1:3,i,j,k) )
                     end do               ; end do                ; end do
                     !if (NVARS == 6) then
                     !   e % U_x(6,:,:,:) = e % Q_x(6,:,:,:)
                     !   e % U_y(6,:,:,:) = e % Q_y(6,:,:,:)
                     !   e % U_z(6,:,:,:) = e % Q_z(6,:,:,:)
                     !end if
                  else
                     e % U_x = e % Q_x(1:NGRADVARS,:,:,:)
                     e % U_y = e % Q_y(1:NGRADVARS,:,:,:)
                     e % U_z = e % Q_z(1:NGRADVARS,:,:,:)
                  end if

               end if
               if (self % hasSensor) then
                   read(fid) e % sensor
               end if

               if (hasUt_NS) then
                   allocate( e % ut_NS(1,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
                   read(fid) e % ut_NS
               end if 

               if (hasMu_NS) then
                   allocate( e % mu_NS(1,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
                   read(fid) e % mu_NS
               end if 

               if (hasWallY) then
                   allocate( e % wallY(1,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
                   read(fid) e % wallY
               end if 

               if (hasMu_sgs) then
                   allocate( e % mu_sgs(1,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
                   read(fid) e % mu_sgs
               end if

               if (hasSource) then
                   allocate( e % source(1:NVARS,0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3)) )
                   read(fid) e % source
               end if

               end associate
            end do

         else

            do eID = 1, self % no_of_elements
               associate ( e => self % elements(eID) )
               call getSolutionFileArrayDimensions(fid,arrayDimensions)
!
!              Allocate memory for the statistics
!              ----------------------------------
               e % Nsol(1:3) = arrayDimensions(2:4) - 1
               allocate( e % stats(1:NSTAT,0:e % Nsol(1), 0:e % Nsol(2), 0:e % Nsol(3)) )
!
!              Read data
!              ---------
               read(fid) e % stats
               end associate
            end do
         end if
!
!        Close file
!        ----------
         close(fid)

!        Read mpi ranks (if present)
!        ---------------------------
         if (hasMPIranks) then
               call readPartitionFile(self)
         end if
         
!
!        Describe the solution
!        ---------------------
         write(msg,'(A,A,A)') 'Solution file "',trim(solutionName),'":'
         write(STD_OUT,'(/)')
         call SubSection_Header(trim(msg))

         if ( self % isStatistics ) then
            write(STD_OUT,'(30X,A,A30)') "->","File is statistics file."

         else
            if ( self % hasGradients ) then
               write(STD_OUT,'(30X,A,A40,A)') "->","Solution file contains gradients: ", "yes"
            else
               write(STD_OUT,'(30X,A,A40,A)') "->","Solution file contains gradients: ", "no"
            end if
         end if

         write(STD_OUT,'(30X,A,A30,I0)') "->","Iteration: ", iter
         write(STD_OUT,'(30X,A,A30,ES10.3)') "->","Time: ", self % time
         write(STD_OUT,'(30X,A,A30,F7.3)') "->","Reference velocity: ", self % refs(V_REF)
         write(STD_OUT,'(30X,A,A30,F7.3)') "->","Reference density: ", self % refs(RHO_REF)
         write(STD_OUT,'(30X,A,A30,F7.3)') "->","Reference Temperature: ", self % refs(T_REF)
         write(STD_OUT,'(30X,A,A30,F7.3)') "->","Reference Mach number: ", self % refs(MACH_REF)
      end subroutine Mesh_ReadSolution
      
      
      subroutine readPartitionFile(self)
         implicit none

         !-arguments-----------------------------------------------
         class(Mesh_t)   , intent(inout)  :: self

         !-local-variables-----------------------------------------
         integer                    :: pos, nelem, eID, fID
         !---------------------------------------------------------

         open(newunit = fID, file=trim(partitionFileName),action='read')

            read(fID,*) nelem

            do eID = 1, nelem
               read(fID,*) self % elements(eID) % mpi_rank
            end do

         close(fID)
      end subroutine readPartitionFile
      
      subroutine readBoundaryFile(boundaries)
         implicit none
         !-arguments-----------------------------------------------
         type(Boundary_t), allocatable, intent(inout)  :: boundaries(:)

         !-local-variables-----------------------------------------
         integer                    :: pos, fd, no_of_boundaries,bID
         !---------------------------------------------------------

         open(newunit = fd, file=trim(boundaryFileName),action='read')

            read(fd,*) no_of_boundaries
            allocate ( boundaries(no_of_boundaries) )

            do bID = 1, no_of_boundaries

               read(fd,*) boundaries(bID) % Name
               read(fd,*) boundaries(bID) % no_of_faces

               allocate ( boundaries(bID) % elements    (boundaries(bID) % no_of_faces) )
               allocate ( boundaries(bID) % elementSides(boundaries(bID) % no_of_faces) )
			   allocate ( boundaries(bID) % cornerDomain(boundaries(bID) % no_of_faces) )
			   allocate ( boundaries(bID) % edgeDomain  (boundaries(bID) % no_of_faces) )

               read(fd,*) boundaries(bID) % elements
               read(fd,*) boundaries(bID) % elementSides
            end do

         close(fd)

      end subroutine readBoundaryFile

      Subroutine getNVARS(originalDim, useFlowEq)
!     *****************************************************
!     get NVARS from the hsol array size
!     in case that the array size is not NCONS in horses3d
!     such as the stats, use the definitions from physics.
!     *****************************************************
          implicit none
          integer, intent(in)           :: originalDim
          logical, intent(in)           :: useFlowEq

      NVARS = originalDim

      ! use simple hardcoded values taken from physics.
      ! todo: link from horses3d code
      if (useFlowEq) then
          select case (trim(flowEq))
              case ("ns")
                  NVARS = 5
              case ("nssa")
                  NVARS = 6
               case ("nstpt")
                  NVARS = 6
              case ("ins")
                  NVARS = 5
              case ("ch")
                  NVARS = 1
              case ("mu")
                  NVARS = 5
              case ("caa")
                  NVARS = 5
              case default
                  NVARS = 5
          end select
      end if

      NGRADVARS = NVARS
      select case(NVARS)
      case(1)      ; NGRADVARS = 1
      case(6)      ; NGRADVARS = 3
      case default ; NGRADVARS = 3
      end select
          
      End Subroutine getNVARS
end module Storage
