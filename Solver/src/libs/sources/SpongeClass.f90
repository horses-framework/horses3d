!//////////////////////////////////////////////////////
!
!This class represents the numerical source term as a sponge to complement BC

!Limitations:
! 1) 
!ToDO
! 1) Assert that radious, amplitude, ramp ... are positive

#include "Includes.h"
Module SpongeClass  !
#if defined(FLOW) 
    use SMConstants
    use HexMeshClass
    use SolutionFile
    use PhysicsStorage, only: NCONS
    Implicit None

    private
    public sponge, addSourceSponge, ConstructSponge, DestructSponge, UpdateBaseFlowSponge, WriteBaseFlowSponge, creatRampSponge

    !definition of sponge class 
    type sponge_t
        integer                                                  :: numOfSponges     ! number of sponges 
        integer                                                  :: nElements        ! number of elements in sponge region
        integer                                                  :: nElementsAll     ! number of elements in sponge region in all partitions
        integer, dimension(:), allocatable                       :: elementIndexMap  ! map from eID of mesh to sponge arrays
        real(kind=RP)                                            :: amplitude        ! amplitude of the source term
        real(kind=RP)                                            :: delta            ! temporal filter width
        real(kind=RP)                                            :: rampWidth        ! length of the ramp width
        real(kind=RP), dimension(:), allocatable                 :: fringeStart      ! start of Fringe Region
        real(kind=RP), dimension(:), allocatable                 :: fringeEnd        ! end of Fringe Region
        real(kind=RP), dimension(:), allocatable                 :: deltaRise        ! rampRise of Fringe Region
        real(kind=RP), dimension(:), allocatable                 :: deltaFall        ! rampFall of Fringe Region
        real(kind=RP), dimension(:,:),allocatable                :: x0               ! position of start of the sponge (for cylindrical in the center)
        real(kind=RP), dimension(:), allocatable                 :: radious          ! radious of the ramp zone in cylindrical/cirular sponges
        real(kind=RP), dimension(:,:) ,allocatable               :: axis             ! axis vector of the sponge. In cylindrical axis of the cylinder, in cartesian, the aligned vector
        character(len=STRING_CONSTANT_LENGTH),dimension(:), allocatable       :: shapeType        ! shape of the sponge, either cartesian (aligned with an axis) or cylindrical
        character(len=STRING_CONSTANT_LENGTH)                    :: initialFileName  ! file name to load the initial base flow
        character(len=STRING_CONSTANT_LENGTH)                    :: solutionFileName ! file name to write base flow
        logical                                                  :: readBaseFLowFlag     ! read base flow from file or use instantaneous Q to start
        logical                                                  :: useMovingAverage ! to use Qbase as a moving average
        logical                                                  :: isActive = .false.

    end type sponge_t

  type(sponge_t)                                                :: sponge


    integer,                       parameter :: KEYWORD_LENGTH             = 132
    character(len=KEYWORD_LENGTH), parameter :: SPONGE_CYLINDRICAL_NAME    = "cylindrical"
    character(len=KEYWORD_LENGTH), parameter :: SPONGE_CARTESIAN_NAME      = "cartesian"   
    character(len=KEYWORD_LENGTH), parameter :: SPONGE_FRINGE_NAME         = "fringe"   

    enum, bind(C)
    enumerator :: SPONGE_CYLINDRICAL = 1, SPONGE_CARTESIAN, FRINGE
    end enum

  contains

!/////////////////////////////////////////////////////////////////////////
!           SPONGE PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////

    Subroutine ConstructSponge(self, mesh, controlVariables)
        use FileReadingUtilities, only: getRealArrayFromString
        use FTValueDictionaryClass
        use MPI_Process_Info
        use Headers
        use mainKeywordsModule
        use FileReadingUtilities, only: getFileName
        Implicit None
        type(sponge_t)                                          :: self
        type(HexMesh), intent(inout)                            :: mesh
        type(FTValueDictionary), intent(in)                     :: controlVariables

        !local variables
        character(len=STRING_CONSTANT_LENGTH)                   :: tmp, numOfSponges, coordinates, axis, fileNameControl, solution_file, restart_name, restart_baseflow_name
        logical                                                 :: fileExists
        integer                                                 :: i 
        logical                                                 :: useNumberedKeys !for cases where there is no number key in the control file when using a single sponge

        self % isActive = .false.
        if (.not. controlVariables % logicalValueForKey("use sponge")) return
        
        self % numOfSponges = controlVariables % getValueOrDefault("number of sponges",1)
        self % amplitude = controlVariables % getValueOrDefault("sponge amplitude",1.0_RP)
        self % rampWidth = controlVariables % getValueOrDefault("sponge ramp width",1.0_RP)
        self % delta = controlVariables % getValueOrDefault("sponge temporal width",1.0_RP)
        self % useMovingAverage = controlVariables % logicalValueForKey("sponge use moving average")
        
        
        allocate(self % shapeType(self % numOfSponges))
        allocate(self % radious(self % numOfSponges))
        allocate(self % axis(self % numOfSponges,NDIM))
        allocate(self % x0(self % numOfSponges,NDIM))
        allocate(self % fringeStart(self % numOfSponges),self % fringeEnd(self % numOfSponges),self % deltaRise(self % numOfSponges),self % deltaFall(self % numOfSponges))
 
        do i = 1, self% numOfSponges
            write(tmp, '("sponge shape ",I0)') i
            if (controlVariables % containsKey(trim(tmp))) then
                useNumberedKeys = .true.
            else
                useNumberedKeys = .false.
                write(tmp, '("sponge shape")')
            end if

            self % shapeType(i) = trim(controlVariables % stringValueForKey(tmp, requestedLength = STRING_CONSTANT_LENGTH))
            if(self % shapeType(i) == "cylindrical") then
                if (useNumberedKeys) then
                    write(tmp, '("sponge radious ",I0)') i
                else
                    write(tmp, '("sponge radious")')
                end if
                self % radious(i) = controlVariables % getValueOrDefault(tmp,0.0_RP)
            endif
            if((self % shapeType(i) == "cylindrical").OR.(self % shapeType(i) == "cartesian")) then
                if (useNumberedKeys) then
                    write(tmp, '("sponge axis ",I0)') i
                else
                    write(tmp, '("sponge axis")')
                end if
                axis = controlVariables % stringValueForKey(tmp, requestedLength = STRING_CONSTANT_LENGTH)
                self % axis(i,:) = getRealArrayFromString(trim(axis))
                !normalize axis
                self % axis(i,:) = self % axis(i,:)/sqrt(sum(self % axis(i,:)**2))
                
                if (useNumberedKeys) then
                    write(tmp, '("sponge start position ",I0)') i
                else
                    write(tmp, '("sponge start position")')
                end if
                coordinates = controlVariables % stringValueForKey(tmp, requestedLength = STRING_CONSTANT_LENGTH)
                self % x0(i,:) = getRealArrayFromString(trim(coordinates))
            end if 
            
            if(self % shapeType(i) == "fringe") then
                if (useNumberedKeys) then
                    write(tmp, '("fringe start ",I0)') i
                else
                    write(tmp, '("fringe start")')
                end if
                self % fringeStart(i) = controlVariables % getValueOrDefault(tmp,0.0_RP)
                if (useNumberedKeys) then
                    write(tmp, '("fringe end ",I0)') i
                else
                    write(tmp, '("fringe end")')
                end if
                self % fringeEnd(i) = controlVariables % getValueOrDefault(tmp,0.0_RP)
                if (useNumberedKeys) then
                    write(tmp, '("delta rise ",I0)') i
                else
                    write(tmp, '("delta rise")')
                end if
                self % deltaRise(i) = controlVariables % getValueOrDefault(tmp,0.0_RP)
                if (useNumberedKeys) then
                    write(tmp, '("delta fall ",I0)') i
                else
                    write(tmp, '("delta fall")')
                end if
                self % deltaFall(i) = controlVariables % getValueOrDefault(tmp,0.0_RP)
            endif
        end do 

!       Get the file name of the solution
!       ---------------------------------
        solution_file = controlVariables % stringValueForKey( solutionFileNameKey, requestedLength = STRING_CONSTANT_LENGTH )
        solution_file = trim(getFileName(solution_file))
        write(self % solutionFileName,'(2A)')  TRIM(solution_file),'_baseflow'

!       Get the file name for the base flow
!       -----------------------------------
        self % readBaseFLowFlag = .false.
        restart_name = controlVariables % stringValueForKey( restartFileNameKey, requestedLength = STRING_CONSTANT_LENGTH )
        restart_name = trim(getFileName(restart_name))
        write(restart_baseflow_name,'(2A)')  TRIM(restart_name),'_baseflow.hsol'
        inquire(file=trim(restart_baseflow_name), exist=fileExists)
        if (fileExists) then
            self % initialFileName = trim(restart_baseflow_name)
            self % readBaseFLowFlag = .true.
        else
            if (controlVariables % containsKey("sponge base flow file")) then
                self % initialFileName = trim(controlVariables % stringValueForKey("sponge base flow file", requestedLength = STRING_CONSTANT_LENGTH))
                self % readBaseFLowFlag = .true.
            end if
        end if 

        self % isActive = .true.

        ! create arrays and pre calculate values
        call creatRampSponge(self,mesh)
        call initializeBaseFlow(self,mesh)
        
        if ( .not. MPI_Process % isRoot ) return
        call Subsection_Header("Sponge")
        write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of elements: ", self % nElementsAll
        write(STD_OUT,'(30X,A,A28,F10.2)') "->", "Amplitude: ", self % amplitude
        write(STD_OUT,'(30X,A,A28,F10.2)') "->", "Ramp width: ", self % rampWidth
        write(STD_OUT,'(30X,A,A28,L1)') "->", "Use moving average: ", self % useMovingAverage
        do i = 1, self% numOfSponges
            write(*,*) 
            write(STD_OUT,'(30X,A,A28,I0)') "->", "Sponge: ", i 
            write(STD_OUT,'(30X,A,A28,A)') "->", "Shape: ", self % shapeType(i)
            if(self % shapeType(i) == "cylindrical") then
               write(STD_OUT,'(30X,A,A28,F10.2)') "->", "Ramp radious: ", self % radious(i)
            endif
            if((self % shapeType(i) == "cylindrical").OR.(self % shapeType(i) == "cartesian")) then
                write(STD_OUT,'(30X,A,A28,3(F10.2))') "->", "Axis: ", self % axis(i,:)
                write(STD_OUT,'(30X,A,A28,3(F10.2))') "->", "Ramp start: ", self % x0(i,:)
            end if
            if(self % shapeType(i) == "fringe") then
                write(STD_OUT,'(30X,A,A28,3(F10.2))') "->", "Fringe Start: ", self % fringeStart(i)
                write(STD_OUT,'(30X,A,A28,3(F10.2))') "->", "Fringe End  : ", self % fringeEnd(i)
                write(STD_OUT,'(30X,A,A28,3(F10.2))') "->", "delta Rise  : ", self % deltaRise(i)
                write(STD_OUT,'(30X,A,A28,3(F10.2))') "->", "delta Fall  : ", self % deltaFall(i)
            end if 
        end do 

         if (self % readBaseFLowFlag) write(STD_OUT,'(30X,A,A28,A)') "->", "Initial base flow file: ", self % initialFileName

    End Subroutine ConstructSponge
!
    Subroutine DestructSponge(self)
        Implicit None
        class(sponge_t), intent(inout)                   :: self

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return
        safedeallocate(self % elementIndexMap)
        safedeallocate(self % x0)
        safedeallocate(self % radious)
        safedeallocate(self % axis)
        safedeallocate(self % shapeType)
        safedeallocate(self % fringeStart)
        safedeallocate(self % fringeEnd)
        safedeallocate(self % deltaRise)
        safedeallocate(self % deltaFall)

    End Subroutine DestructSponge
!
    Subroutine creatRampSponge(self, mesh)
        use MPI_Process_Info
#ifdef _HAS_MPI_
        use mpi
#endif
        Implicit None
        type(sponge_t)                                          :: self
        type(HexMesh), intent(inout)                            :: mesh

        !local variables
        real(kind=RP), dimension(:,:,:), allocatable            :: xStar, sigma
        real(kind=RP), dimension(NDIM)                          :: r_vector
        real(kind=RP)                                           :: xData(2), Sf(2)
        logical, dimension(:), allocatable                      :: hasSponge
        integer                                                 :: i, j, k, eID, counter, spongeEID, ierr
        integer, dimension(NDIM)                                :: Nxyz
        integer                                                 :: sponge_number
        integer                                                 :: whichSponge = -1

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

        if (allocated(hasSponge)) deallocate(hasSponge)
        allocate(hasSponge(mesh % no_of_elements))
        hasSponge = .false.
        
        counter = 0
        do sponge_number=1 , self % numOfSponges 

            select case (self % shapeType(sponge_number))
            case (SPONGE_CYLINDRICAL_NAME)
                whichSponge = SPONGE_CYLINDRICAL
            case (SPONGE_CARTESIAN_NAME)
                whichSponge = SPONGE_CARTESIAN
            case (SPONGE_FRINGE_NAME)
                whichSponge = FRINGE
            case default
                print*, "Sponge name not recognized."
                errorMessage(STD_OUT)
                error stop
            end select
                
            do eID = 1, mesh % no_of_elements
                associate(e => mesh % elements(eID))
                    Nxyz = e % Nxyz
                    if (allocated(xStar)) deallocate(xStar)
                    allocate(xStar(0:Nxyz(1),0:Nxyz(2),0:Nxyz(3)))
                    do k = 0, Nxyz(3) ; do j = 0, Nxyz(2) ; do i = 0, Nxyz(1)
                        r_vector(:) = e % geom % x(:,i,j,k) - self % x0(sponge_number,:)
                        select case (whichSponge)
                        case (SPONGE_CYLINDRICAL)
                            ! remove components in the axis direction
                            r_vector(:) = r_vector(:) - sum((e % geom % x(:,i,j,k) - self % x0(sponge_number,:))*self % axis(sponge_number,:))*self % axis(sponge_number,:)
                            ! in this case xStar is actually rdiff ^2
                            xStar(i,j,k) = sum(r_vector*r_vector) - POW2(self % radious(sponge_number))
                        case (SPONGE_CARTESIAN)
                            ! in this case xstar is the distance to the plane
                            xStar(i,j,k) = sum(r_vector(:)*self % axis(sponge_number,:))
                        case (FRINGE)
                            ! in this case xstar is the distance to the plane
                            if ((e % geom % x(1,i,j,k).ge.self % fringeStart(sponge_number)).and.(e % geom % x(1,i,j,k).le.self % fringeEnd(sponge_number))) then
                                xStar(i,j,k) = 1.0_RP
                            end if 
                        end select
                    end do         ; end do          ; end do
                    if (any(xStar .ge. 0.0_RP) .AND. .not.hasSponge(eID)  ) then
                        hasSponge(eID) = .true.
                        counter = counter + 1
                    end if 
                end associate
            end do
        end do

        self % nElements = counter
        if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
        call mpi_allreduce(self % nElements, self % nElementsAll, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        else
            self % nElementsAll = self % nElements
        end if

        ! create mapping array
        if (allocated(self % elementIndexMap)) deallocate(self % elementIndexMap)
        allocate( self  %  elementIndexMap(self % nElements) )

        counter = 0
        do eID = 1, mesh % no_of_elements
            if (hasSponge(eID)) then
                counter = counter + 1
                self % elementIndexMap(counter) = eID
            end if 
        end do

        ! now calculate the ramp amplitude
        do sponge_number=1 , self % numOfSponges 

            select case (self % shapeType(sponge_number))
            case (SPONGE_CYLINDRICAL_NAME)
                whichSponge = SPONGE_CYLINDRICAL
            case (SPONGE_CARTESIAN_NAME)
                whichSponge = SPONGE_CARTESIAN
            case (SPONGE_FRINGE_NAME)
                whichSponge = FRINGE
            end select

            do spongeEID = 1, self % nElements
                eID = self % elementIndexMap(spongeEID)
                associate(e => mesh % elements(eID))
                    Nxyz = e % Nxyz
                    if (allocated(xStar)) deallocate(xStar)
                    if (allocated(sigma)) deallocate(sigma)
                    allocate(xStar(0:Nxyz(1),0:Nxyz(2),0:Nxyz(3)))
                    allocate(sigma(0:Nxyz(1),0:Nxyz(2),0:Nxyz(3)))
                    
                    do k = 0, Nxyz(3) ; do j = 0, Nxyz(2) ; do i = 0, Nxyz(1)
                        r_vector(:) = e % geom % x(:,i,j,k) - self % x0(sponge_number,:)
                        select case (whichSponge)
                        case (SPONGE_CYLINDRICAL)
                            ! remove components in the axis direction
                            r_vector(:) = r_vector(:) - sum((e % geom % x(:,i,j,k) - self % x0(sponge_number,:))*self % axis(sponge_number,:))*self % axis(sponge_number,:)
                            ! xStar is non-dimensionalized by the width of the ramp, since the difference is squared, the width is too
                            xStar(i,j,k) = sqrt(abs(sum(r_vector*r_vector) - POW2(self % radious(sponge_number)))) / self % rampWidth
                        case (SPONGE_CARTESIAN)
                            xStar(i,j,k) = sum(r_vector(:)*self % axis(sponge_number,:))/(self % rampWidth)
                        case (FRINGE)
                            xData(1)=(e % geom % x(1,i,j,k)-self % fringeStart(sponge_number))/self % deltaRise(sponge_number)
                            xData(2)=((e % geom % x(1,i,j,k)-self % fringeEnd(sponge_number))/self % deltaFall(sponge_number))+1
                            
                            if (xData(1).LE.0.0_RP) then
                                Sf(1)=0.0_RP
                            else if (xData(1).GE.1.0_RP) then
                                Sf(1)=1.0_RP
                            else
                                !if (abs(xData(1) - 1.0_RP) > epsilon(xData(1)) .and. abs(xData(1)) > epsilon(xData(1))) then
                                    Sf(1)=1.0_RP/(1.0_RP+exp((1.0_RP/max((xData(1)-1.0_RP),0.01_RP))+(1.0_RP/max(xData(1),0.01_RP))))
                                ! else
                                    ! Sf(1)=1.0_RP
                                ! end if 
                            end if
                            
                            if (xData(2).LE.0.0_RP) then
                                Sf(2)=0.0_RP
                            else if (xData(2).GE.1.0_RP) then
                                Sf(2)=1.0_RP
                            else
                                !if (abs(xData(2) - 1.0_RP) > epsilon(xData(2)) .and. abs(xData(2)) > epsilon(xData(2))) then
                                  Sf(2) = 1.0_RP / (1.0_RP + exp((1.0_RP / max((xData(2) - 1.0_RP),0.01_RP)) + (1.0_RP / max(xData(2),0.01_RP))))
                                ! else
                                  ! Sf(2) = 1.0_RP
                                ! end if
                            end if
                            
                            e % storage % intensitySponge(i,j,k)= self % amplitude*(Sf(1)-Sf(2))+e % storage % intensitySponge(i,j,k)
                        end select
                    end do         ; end do          ; end do
                    if((self % shapeType(sponge_number) == "cylindrical").OR.(self % shapeType(sponge_number) == "cartesian")) then
                         if (any(xStar .ge. 0.0_RP)) then
                            ! limit xStar to [0,1] since after 1 should be constant at the amplitude value
                            xStar = max(0.0_RP,xStar)
                            xStar = min(1.0_RP,xStar)
                            ! Sponge Ramping Function, taken from Beck, A., and Munz, C.-D., Direct Aeroacoustic Simulations Based on High Order Discontinuous Galerkin Schemes, Springer, Cham, 2018, Vol. 579
                            sigma  = 6.0_RP*xStar**5.0_RP - 15.0_RP*xStar**4.0_RP + 10.0_RP*xStar**3.0_RP
                            ! limit sigms <=1 since after 1 should be constant at the amplitude value
                            sigma  = MIN(1.0_RP,sigma)
                            ! pre computed intensity, including amplitude and ramp damping
                            e % storage % intensitySponge(:,:,:) = max(e % storage % intensitySponge(:,:,:), self % amplitude * sigma(:,:,:))
                            end if
                    end if 
                end associate
            end do 
        end do

        ! deallocate(xStar, sigma, hasSponge)
        if (allocated(xStar)) deallocate(xStar)
        if (allocated(sigma)) deallocate(sigma)
        if (allocated(hasSponge)) deallocate(hasSponge)

    End Subroutine creatRampSponge
!
    Subroutine addSourceSponge(self,mesh,Level)
        Implicit None
        type(sponge_t)                                          :: self
        type(HexMesh), intent(inout)                            :: mesh
        integer, intent(in), optional                           :: Level   ! MLRK: only elements with MLevel >= Level
    
        !local variables
        integer                                                 :: i, j, k, eID, spongeEID

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

!$omp do schedule(runtime) private(i,j,k,eID)
        do spongeEID = 1, self % nElements
            eID = self % elementIndexMap(spongeEID)
            if ( present(Level) ) then
                if ( mesh % elements(eID) % MLevel < Level ) cycle
            end if
            associate(e => mesh % elements(eID))
                do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                    e % storage % S_NS(:,i,j,k) = - e % storage % intensitySponge(i,j,k) * (e % storage % Q(:,i,j,k) - e % storage % QbaseSponge(:,i,j,k))
                end do         ; end do          ; end do
            end associate
        end do
!$omp end do

    End Subroutine addSourceSponge
!
    Subroutine initializeBaseFlow(self,mesh)
        use ElementClass
        Implicit None
        type(sponge_t)                                          :: self
        type(HexMesh), intent(inout)                            :: mesh

        !local variables
        integer                                                 :: eID
        
        if (self % readBaseFLowFlag) then
            call readBaseFlow(self,mesh)
        else
!$omp do schedule(runtime) private(eID)
            do eID = 1, mesh % no_of_elements
                associate(e => mesh % elements(eID))
                    e % storage % QbaseSponge(:,:,:,:) = e % storage % Q(:,:,:,:) 
                end associate
            end do
!$omp end do
        end if 
!
    End Subroutine initializeBaseFlow

    ! advance in time base flow as single euler step, taken from Beck 2018
    Subroutine UpdateBaseFlowSponge(self, mesh, dt)
        Implicit None
        type(sponge_t)                                          :: self
        type(HexMesh), intent(inout)                            :: mesh
        real(kind=RP), intent(in)                               :: dt

        !local variables
        integer                                                 :: eID, spongeEID

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

!       Only update for moving average
!       ------------------------
        if (.not. self % useMovingAverage) return

!$omp do schedule(runtime)
        do eID = 1, mesh % no_of_elements
            associate(e => mesh % elements(eID))
                e % storage % QbaseSponge(:,:,:,:) = e % storage % QbaseSponge(:,:,:,:) + (e % storage % Q(:,:,:,:) - e % storage % QbaseSponge(:,:,:,:)) * (dt / self%delta)
            end associate
        end do
!$omp end do

    End Subroutine UpdateBaseFlowSponge
!
    Subroutine  readBaseFlow(self,mesh)
        Implicit None
        type(sponge_t)                                          :: self
        type(HexMesh), intent(inout)                            :: mesh

        !local variables
        INTEGER                        :: fID, eID, fileType, no_of_elements, flag 
        integer(kind=AddrInt)          :: pos, padding
        integer                        :: Nxp1, Nyp1, Nzp1, no_of_eqs, array_rank
        real(kind=RP), allocatable     :: Q(:,:,:,:)
        integer, dimension(NDIM)                                :: Nxyz

!       Get the file type
!       -----------------
        fileType = getSolutionFileType(trim(self % initialFileName))

        select case (fileType)
        case(SOLUTION_FILE)
           padding = 1*NCONS
        case default
           print*, "Only solution file format is accepted"
           errorMessage(STD_OUT)
           stop
        end select
!
!       Read the number of elements
!       ---------------------------
        no_of_elements = getSolutionFileNoOfElements(trim(self % initialFileName))

        if ( no_of_elements .ne. mesh % no_of_allElements ) then
           write(STD_OUT,'(A,A)') "The number of elements stored in the restart file ", &
                                  "do not match that of the mesh file"
           errorMessage(STD_OUT)
           stop
        end if
!
!       Read the terminator indicator
!       -----------------------------
        flag = getSolutionFileDataInitFlag(trim(self % initialFileName))

        if ( flag .ne. BEGINNING_DATA ) then
           print*, "Beginning data flag was not found in the file."
           errorMessage(STD_OUT)
           stop
        end if
!
!       Read elements data
!       ------------------
        fID = putSolutionFileInReadDataMode(trim(self % initialFileName))
        
        do eID = 1, size(mesh % elements)
           associate( e => mesh % elements(eID) )
               Nxyz = e % Nxyz
               if (allocated(Q)) deallocate(Q)
               allocate(Q(NCONS, 0:Nxyz(1), 0:Nxyz(2), 0:Nxyz(3)))
               pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*padding*e % offsetIO * SIZEOF_RP
               read(fID, pos=pos) array_rank
               read(fID) no_of_eqs, Nxp1, Nyp1, Nzp1
               if (      ((Nxp1-1) .ne. e % Nxyz(1)) &
                    .or. ((Nyp1-1) .ne. e % Nxyz(2)) &
                    .or. ((Nzp1-1) .ne. e % Nxyz(3)) &
                    .or. (no_of_eqs .ne. NCONS )       ) then
                  write(STD_OUT,'(A,I0,A)') "Error reading restart file: wrong dimension for element "&
                                              ,eID,"."

                  write(STD_OUT,'(A,I0,A,I0,A,I0,A)') "Element dimensions: ", e % Nxyz(1), &
                                                                        " ,", e % Nxyz(2), &
                                                                        " ,", e % Nxyz(3), &
                                                                        "."

                  write(STD_OUT,'(A,I0,A,I0,A,I0,A)') "Restart dimensions: ", Nxp1-1, &
                                                                        " ,", Nyp1-1, &
                                                                        " ,", Nzp1-1, &
                                                                        "."

                  errorMessage(STD_OUT)
                  stop
               end if

               read(fID) Q
               e % storage % QbaseSponge(:,:,:,:) = Q(:,:,:,:)
          end associate
        end do
        deallocate(Q)
!
!       Close the file
!       --------------
        close(fID)

    End Subroutine  readBaseFlow
!
    Subroutine  WriteBaseFlowSponge(self,mesh,iter,time,last)
        use FluidData, only: thermodynamics, refValues, dimensionless
        Implicit None
        class(sponge_t)                                         :: self
        type(HexMesh), intent(in)                               :: mesh
        integer,             intent(in)        :: iter
        real(kind=RP),       intent(in)        :: time
        logical, intent(in), optional          :: last
!
!       ---------------
!       Local variables
!       ---------------
!
        integer                          :: fid, eID
        integer(kind=AddrInt)            :: pos, padding
        real(kind=RP)                    :: refs(NO_OF_SAVED_REFS) 
        real(kind=RP), allocatable       :: Q(:,:,:,:)
        integer, dimension(NDIM)                                :: Nxyz
        character(len=LINE_LENGTH)                              :: name
        logical                                                 :: notAddIter

!       Check if is activated
!       ------------------------
        if (.not. self % isActive) return

!       Only update for moving average
!       ------------------------
        if (.not. self % useMovingAverage) return

        if (present(last)) then
            notAddIter = last
        else
            notAddIter = .false.
        end if 

        if (notAddIter) then
            write(name,'(2A)')  trim(self % solutionFileName),'.hsol'
        else
            write(name,'(2A,I10.10,A)')  trim(self % solutionFileName),'_',iter,'.hsol'
        end if 
!
!       Gather reference quantities
!       ---------------------------
#if defined(NAVIERSTOKES)
        refs(GAMMA_REF) = thermodynamics % gamma
        refs(RGAS_REF)  = thermodynamics % R
        refs(RHO_REF)   = refValues      % rho
        refs(V_REF)     = refValues      % V
        refs(T_REF)     = refValues      % T
        refs(MACH_REF)  = dimensionless  % Mach
#elif defined(INCNS)
        refs(GAMMA_REF) = 0.0_RP
        refs(RGAS_REF)  = 0.0_RP
        refs(RHO_REF)   = refValues      % rho
        refs(V_REF)     = refValues      % V
        refs(T_REF)     = 0.0_RP
        refs(MACH_REF)  = 0.0_RP
#elif defined(Multiphase)
        refs = 0.0_RP
#elif defined(ACOUSTIC)
        refs(GAMMA_REF) = thermodynamics % gamma
        refs(RGAS_REF)  = thermodynamics % R
        refs(RHO_REF)   = refValues      % rho
        refs(V_REF)     = refValues      % V
        refs(T_REF)     = refValues      % T
        refs(MACH_REF)  = dimensionless  % Mach
#endif
!       Create new file
!       ---------------
        call CreateNewSolutionFile(trim(name),SOLUTION_FILE, mesh % nodeType, mesh % no_of_allElements, iter, time, refs)
        
        padding = NCONS
!
!       Write arrays
!       ------------
        fID = putSolutionFileInWriteDataMode(trim(name))

        do eID = 1, mesh % no_of_elements
            associate( e => mesh % elements(eID) )
                Nxyz = e % Nxyz
                if (allocated(Q)) deallocate(Q)
                allocate(Q(NCONS, 0:Nxyz(1), 0:Nxyz(2), 0:Nxyz(3)))
                pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*padding*e % offsetIO * SIZEOF_RP
                Q(1:NCONS,:,:,:)  = e % storage % QbaseSponge(:,:,:,:)
                call writeArray(fid, Q, position=pos)
            end associate
        end do
        deallocate(Q)

        close(fid)
!
!       Close the file
!       --------------
        call SealSolutionFile(trim(name))

    End Subroutine WriteBaseFlowSponge
!
#endif
End Module SpongeClass
