!//////////////////////////////////////////////////////
!
!This module handle the fictitious surfaces that will be export to save solution of a set of faces, such as slices and FWH acoustic analogy
! The solution is saved before the first time step and at a fixed dt

#include "Includes.h"
Module SurfaceMesh
    use SMConstants
    use HexMeshClass
    use ZoneClass
    use AutosaveClass
    use FTValueDictionaryClass, only: FTValueDictionary
    use MPI_Process_Info
    Implicit None
!   
    private
!
    public SURFACE_TYPE_FWH, SURFACE_TYPE_SLICE, FWH_POSITION, SURFACE_TYPE_ACOUSTICS 
    public surfacesMesh
#if defined(NAVIERSTOKES)
    public getU_tauInSurfaces, getWallDistInSurfaces
#endif
!
    integer, parameter          :: SURFACE_TYPE_FWH          = 1
    integer, parameter          :: SURFACE_TYPE_SLICE        = 2
    integer, parameter          :: SURFACE_TYPE_BC           = 3
    integer, parameter          :: SURFACE_TYPE_ACOUSTICS    = 4
    integer, parameter          :: FWH_POSITION              = 1 ! the fwh surface will be always the first position if exists
!
    ! main class of surfaces mesh
    type SurfaceMesh_t 
        integer                                                 :: numberOfSurfaces
        integer, dimension(:), allocatable                      :: totalFaces       ! number of faces in all partitions for each surface
        integer, dimension(:), allocatable                      :: surfaceTypes     ! type of each surface
        type(Zone_t), dimension(:), allocatable                 :: zones            ! fictitious zones that contains the faces of each surface
        logical, dimension(:), allocatable                      :: surfaceActive    ! flag for each surface
        logical, dimension(:), allocatable                      :: isNoSlip         ! flag use for calculate and save variables of noslip bc
        character(len=LINE_LENGTH), dimension(:), allocatable   :: file_names
        integer, dimension(:,:), allocatable                    :: globalFid        ! for I/O with mpi
        integer, dimension(:,:), allocatable                    :: faceOffset       ! for I/O with mpi
        integer, dimension(:,:), allocatable                    :: elementSide      ! for reading surface from file
        type(Autosave_t)                                        :: autosave         ! save solution parameters
        logical                                                 :: saveGradients    ! flag use for save gradients in bcs and slices
        logical                                                 :: mergeFWHandBC    ! flag to merge 2 surfaces if they have the same BCs
        logical                                                 :: active           ! flag use for save at least one file
        logical                                                 :: saveUt           ! flag use for save friction velocity in bcs
        logical                                                 :: saveTurb         ! flag use for save wall normal distance, and viscosity

        contains

            procedure :: construct        => SurfConstruct
            procedure :: destruct         => SurfDestruct
            procedure :: autosaveConfig   => SurfSaveSolutionConfiguration
            procedure :: saveAllMesh      => SurfSaveAllMesh
            procedure :: saveAllSolution  => SurfSaveAllSolution
            procedure :: loadSolution     => SurfLoadSolution
    end type SurfaceMesh_t 

    type(SurfaceMesh_t)                         :: surfacesMesh
    real(kind=RP)                               :: sliceTolerance
!
   contains 
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!           CLASS PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////////////////////////
!         
    Subroutine SurfConstruct(self, controlVariables, mesh)
!     ************************************************************************
!        The surface is construct as an structure of arrays, all the positions 
!        correspond to the info for a single surface.
!        If there is a FWH surface, it will always be the first one (only one
!        is allowed). Then, the acoustic surface is created (only one is allowed). 
!        Finally, it will be the BC, then each direction the slices.
!     ************************************************************************
!
        use FileReadingUtilities, only: getFileName, getCharArrayFromString, getRealArrayFromString, getPath, removePath
        use mainKeywordsModule
        use Utilities,            only: toLower
        use Headers
        use BoundaryConditions
        use MPI_Process_Info
#ifdef _HAS_MPI_
        use mpi
#endif
        implicit none

        class(SurfaceMesh_t)                                :: self
        type(FTValueDictionary), intent(in)                 :: controlVariables
        type (HexMesh), intent(in)                          :: mesh

        !local variables
        integer                                                 :: numberOfSurfaces, numberOfFaces, numberOfSurfacesBC, surfaceCont
        integer                                                 :: maxNumberFaces, nf, numberOfFacesTotal
        integer, dimension(NDIM)                                :: numberOfSurfacesSlices
        integer                                                 :: i, j, ierr
        integer                                                 :: idSliceX, idSliceY, idSliceZ, idBC
        logical                                                 :: hasFWH, hasSliceX, hasSliceY, hasSliceZ, hasBC, hasAcousticPath
        character(len=LINE_LENGTH)                              :: solution_file, read_str, fileName, zoneName, fullExpression, path
        character(len=LINE_LENGTH), dimension(:), allocatable   :: bcs_names, bcs_namesFWH
        integer, dimension(:), allocatable                      :: FWHelementSide, facesIDs, acousticElementSide
        real(kind=RP), dimension(:), allocatable                :: posSliceX, posSliceY, posSliceZ
        real(kind=RP), dimension(:), allocatable                :: limSliceX, limSliceY, limSliceZ
        logical, dimension(:), allocatable                      :: surfaceActive
        logical                                                 :: isNoSlip, surfaceHasFacesInPartitions

        numberOfSurfaces = 0
        self % active = .false.
!
!       Get the possible surfaces to save
!       ---------------------------------
        hasFWH    = controlVariables % containsKey("acoustic analogy")
        hasSliceX = controlVariables % containsKey("slice x coordinates")
        hasSliceY = controlVariables % containsKey("slice y coordinates")
        hasSliceZ = controlVariables % containsKey("slice z coordinates")
        hasBC     = controlVariables % containsKey("boundaries to save")
        hasAcousticPath = controlVariables % containsKey("acoustic path")
        idSliceX  = 0
        idSliceY  = 0
        idSliceZ  = 0
        idBC      = 0
        numberOfSurfacesBC     = 0
        numberOfSurfacesSlices = 0

				if (hasFWH .and. hasAcousticPath) then
						error stop "acoustic analogy and acoustic path cannot be used at the same time"
				end if
!
!       Get the solution file name
!       --------------------------
        solution_file = controlVariables % stringValueForKey( solutionFileNameKey, requestedLength = LINE_LENGTH )
        solution_file = trim(getFileName(solution_file))
        path = trim(getPath(solution_file))
        ! save in a separated folder
        solution_file = trim(path) // "/surfaces/" // trim(removePath(solution_file))
!
        sliceTolerance = controlVariables % getValueOrDefault("slice tolerance", 1.0e-4_RP)
        self % saveGradients = controlVariables % logicalValueForKey("surface save gradients") .or. controlVariables % logicalValueForKey("save gradients")
        self % saveUt = controlVariables % logicalValueForKey("surface save utau")
        self % saveTurb = controlVariables % logicalValueForKey("surface save turbulent")
!
!       get number of surfaces
!       ----------------------
        if (hasFWH) numberOfSurfaces = numberOfSurfaces + 1

        if (hasAcousticPath) numberOfSurfaces = numberOfSurfaces + 1
!
        if (hasSliceX) then
            read_str = controlVariables % stringValueForKey("slice x coordinates", LINE_LENGTH)
            call toLower(read_str)
            posSliceX = getRealArrayFromString(read_str)
            numberOfSurfacesSlices(1) =  size(posSliceX)
            read_str = controlVariables % stringValueForKey("slice x limits", LINE_LENGTH)
            call toLower(read_str)
            if (trim(read_str) .ne. "") limSliceX = getRealArrayFromString(read_str)
        end if
!
        if (hasSliceY) then
            read_str = controlVariables % stringValueForKey("slice y coordinates", LINE_LENGTH)
            call toLower(read_str)
            posSliceY = getRealArrayFromString(read_str)
            numberOfSurfacesSlices(2) =  size(posSliceY)
            read_str = controlVariables % stringValueForKey("slice y limits", LINE_LENGTH)
            call toLower(read_str)
            if (trim(read_str) .ne. "") limSliceY = getRealArrayFromString(read_str)
        end if
!
        if (hasSliceZ) then
            read_str = controlVariables % stringValueForKey("slice z coordinates", LINE_LENGTH)
            call toLower(read_str)
            posSliceZ = getRealArrayFromString(read_str)
            numberOfSurfacesSlices(3) =  size(posSliceZ)
            read_str = controlVariables % stringValueForKey("slice z limits", LINE_LENGTH)
            call toLower(read_str)
            if (trim(read_str) .ne. "") limSliceZ = getRealArrayFromString(read_str)
        end if 
!
        if (hasBC) then
            read_str = controlVariables % stringValueForKey("boundaries to save", LINE_LENGTH)
            call toLower(read_str)
            call getCharArrayFromString(read_str, LINE_LENGTH, bcs_names)
            numberOfSurfacesBC =  size(bcs_names)
        end if 
        self % mergeFWHandBC = .false.
        if (hasFWH .and. hasBC) then
            if (controlVariables % containsKey("acoustic solid surface")) then
                read_str = controlVariables % stringValueForKey("acoustic solid surface", LINE_LENGTH)
                call toLower(read_str)
                call getCharArrayFromString(read_str, LINE_LENGTH, bcs_namesFWH)
                if (all( bcs_names .eq. bcs_namesFWH )) then
                    self % mergeFWHandBC = .true.
                    numberOfSurfacesBC = 0
                    hasBC = .false.
                end if
            end if
        end if
        numberOfSurfaces = numberOfSurfaces + sum(numberOfSurfacesSlices) + numberOfSurfacesBC
!
!       return if there are no surfaces, the flag active will be false and 
!       all of the arrays wont be allocated. The class procedures will look for the flag before doing anything
        if (numberOfSurfaces .eq. 0) return
!
        self % numberOfSurfaces = numberOfSurfaces
        allocate( self % totalFaces(numberOfSurfaces), self % surfaceTypes(numberOfSurfaces), self % zones(numberOfSurfaces), &
                  self % surfaceActive(numberOfSurfaces), self % file_names(numberOfSurfaces), self%isNoSlip(numberOfSurfaces) )
        self % totalFaces = 0
        self % surfaceActive = .false.
        self % surfaceTypes = 0
!
!       now construct all the surfaces
!       ------------------------------
        surfaceCont = 0
        if (hasFWH) then
            surfaceCont = surfaceCont + 1
            call createFWHSurface(self, controlVariables, mesh, solution_file, FWHelementSide)
            ! for now only fwh is loaded from file, if this changes in future, some changes must be done, allocating the side of all
            ! faces of all the file surfaces
            allocate ( self % elementSide(size(FWHelementSide),1) )
            self % elementSide(:,1) = FWHelementSide
        end if 

        if (hasAcousticPath) then
            write(STD_OUT,*) "Creating acoustic path"
            surfaceCont = surfaceCont + 1
            call createAcousticSurface(self, controlVariables, mesh, solution_file, acousticElementSide, surfaceCont)
            ! for now only one acoustic surface is loaded from file, if this changes in future, some changes must be done, allocating the side of all
            ! faces of all the file surfaces
            allocate ( self % elementSide(size(acousticElementSide),1) )
            self % elementSide(:,1) = acousticElementSide
        end if 
!
        if (hasSliceX) then
            slicex_loop: do i = 1, numberOfSurfacesSlices(1)
                surfaceCont = surfaceCont + 1
                idSliceX = idSliceX + 1
                safedeallocate(facesIDs)
                call getSliceFaces(1, posSliceX(i), mesh, numberOfFaces, facesIDs, limSliceX)
                if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
                    call mpi_allreduce(numberOfFaces, numberOfFacesTotal, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
                else
                    numberOfFacesTotal = numberOfFaces
                end if
                surfaceHasFacesInPartitions = numberOfFacesTotal .gt. 0
                write(fileName,'(A,A,I3.3)')  trim(solution_file), '_sx', idSliceX
                zoneName = "sliceSurface"
                call createSingleSurface(self, surfaceCont, SURFACE_TYPE_SLICE, fileName, zoneName, numberOfFaces, facesIDs, .false., surfaceHasFacesInPartitions)
            end do slicex_loop
        end if 
!
        if (hasSliceY) then
            slicey_loop: do i = 1, numberOfSurfacesSlices(2)
                surfaceCont = surfaceCont + 1
                idSliceY = idSliceY + 1
                safedeallocate(facesIDs)
                call getSliceFaces(2, posSliceY(i), mesh, numberOfFaces, facesIDs, limSliceY)
                if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
                    call mpi_allreduce(numberOfFaces, numberOfFacesTotal, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
                else
                    numberOfFacesTotal = numberOfFaces
                end if
                surfaceHasFacesInPartitions = numberOfFacesTotal .gt. 0
                write(fileName,'(A,A,I3.3)')  trim(solution_file), '_sy', idSliceY
                zoneName = "sliceSurface"
                call createSingleSurface(self, surfaceCont, SURFACE_TYPE_SLICE, fileName, zoneName, numberOfFaces, facesIDs,.false.,surfaceHasFacesInPartitions)
            end do slicey_loop
        end if 
!
        if (hasSliceZ) then
            slicez_loop: do i = 1, numberOfSurfacesSlices(3)
                surfaceCont = surfaceCont + 1
                idSliceZ = idSliceZ + 1
                safedeallocate(facesIDs)
                call getSliceFaces(3, posSliceZ(i), mesh, numberOfFaces, facesIDs,limSliceZ)
                if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
                    call mpi_allreduce(numberOfFaces, numberOfFacesTotal, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
                else
                    numberOfFacesTotal = numberOfFaces
                end if
                surfaceHasFacesInPartitions = numberOfFacesTotal .gt. 0
                write(fileName,'(A,A,I3.3)')  trim(solution_file), '_sz', idSliceZ
                zoneName = "sliceSurface"
                call createSingleSurface(self, surfaceCont, SURFACE_TYPE_SLICE, fileName, zoneName, numberOfFaces, facesIDs,.false.,surfaceHasFacesInPartitions)
            end do slicez_loop
        end if 
!
        if (hasBC .and. .not. self % mergeFWHandBC) then
            bcs_loop: do i = 1, numberOfSurfacesBC
                surfaceCont = surfaceCont + 1
                idBC = idBC + 1
                safedeallocate(facesIDs)
                ! isNoSlip = .false.
                do j = 1, size(mesh % zones)
                    if (trim(bcs_names(i)) .eq. trim(mesh % zones(j) % Name)) exit
                end do
                if (j .gt. size(mesh % zones)) cycle bcs_loop
                ! if ( BCs(j) % bc % BCType .eq. "noslipwall" ) isNoSlip = .true.
                isNoSlip = ( BCs(j) % bc % BCType .eq. "noslipwall" ) 
                numberOfFaces = mesh % zones(j) % no_of_faces
                allocate( facesIDs(numberOfFaces) )
                facesIDs = mesh % zones(j) % faces
                if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
                    call mpi_allreduce(numberOfFaces, numberOfFacesTotal, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
                else
                    numberOfFacesTotal = numberOfFaces
                end if
                surfaceHasFacesInPartitions = numberOfFacesTotal .gt. 0
                write(fileName,'(A,A,I3.3)')  trim(solution_file), '_bc', idBC
                zoneName = "BC_Surface"
                call createSingleSurface(self, surfaceCont, SURFACE_TYPE_BC, fileName, zoneName, numberOfFaces, facesIDs, isNoSlip, surfaceHasFacesInPartitions)
            end do bcs_loop
        end if 
        ! get whether at least one partition has faces of the surface
        if ( (MPI_Process % doMPIAction) ) then
                allocate(surfaceActive(numberOfSurfaces))
#ifdef _HAS_MPI_
            call mpi_allreduce(self % surfaceActive, surfaceActive, numberOfSurfaces, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
            self % surfaceActive = surfaceActive
            deallocate(surfaceActive)
        end if
!
!       Gather the total number of faces and
!       Get the max no of faces between all surfaces en each partition
!       ---------------------------------------------------------------
        maxNumberFaces = 0
        do i = 1, self % numberOfSurfaces
            if (.not. self % surfaceActive(i)) cycle
            if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
                call mpi_allreduce(self % zones(i) % no_of_faces, numberOfFaces, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
            else
                numberOfFaces = self % zones(i) % no_of_faces
            end if
!         
            self % totalFaces(i) = numberOfFaces
!         
            if (maxNumberFaces .lt. self % zones(i) % no_of_faces) maxNumberFaces = self % zones(i) % no_of_faces
        end do
!         
!       Prepare each surface for IO
!       ---------------------------
        allocate( self % globalFid(maxNumberFaces,self % numberOfSurfaces), self % faceOffset(maxNumberFaces,self % numberOfSurfaces) )
        self % globalFid = 0
        self % faceOffset = 0
        do i = 1, self % numberOfSurfaces
            if (.not. self % surfaceActive(i)) cycle
            ! only save for the actual faces of the surface in the partition, the extra positions have a value of 0, by default
            nf = self % zones(i) % no_of_faces
            call SurfacePrepareForIO(self % zones(i), mesh, self % totalFaces(i), self % globalFid(1:nf,i), self % faceOffset(1:nf,i))
        end do
!
        ! change to look if there is at least one true in the array of flags
        self % active = .true.
!
        if ( .not. MPI_Process % isRoot ) return
!
!       create the folder inside the expected RESULTS folder
        write(fullExpression,'(A,A,A)') "mkdir -p ", trim(path), "/surfaces"
        call system(trim(fullExpression))
!
!       Describe the zones
!       ------------------
        if (hasFWH .and. .not. (hasBC .or. hasSliceX .or. hasSliceY .or. hasSliceZ)) return
        call Subsection_Header("Fictitious surfaces zone")
        write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of surfaces: ", self % numberOfSurfaces
        if (.not. all(self % surfaceActive)) write(STD_OUT,'(30X,A,A28,I0)') "->", "Inactive surfaces: ", count(.not. self % surfaceActive, dim=1)
        write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of slices: ", idSliceX + idSliceY + idSliceZ
        if (hasBC) write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of boundaries: ", numberOfSurfacesBC
        if (hasFWH) then 
            write(STD_OUT,'(30X,A,A28,L1)') "->", "Save FWH and BC together: ", self % mergeFWHandBC
            write(STD_OUT ,'(/)')
        end if
!
    End Subroutine SurfConstruct
!
    Subroutine SurfDestruct(self)
!
        implicit none
        class(SurfaceMesh_t), intent(inout)              :: self
!
!       Check if is activated
!       ------------------------
        if (.not. self % active) return
        safedeallocate(self % totalFaces)
        safedeallocate(self % surfaceTypes)
        safedeallocate(self % zones)
        safedeallocate(self % surfaceActive)
        safedeallocate(self % isNoSlip)
        safedeallocate(self % file_names)
        safedeallocate(self % globalFid)
        safedeallocate(self % faceOffset)
        safedeallocate(self % elementSide)
!
    End Subroutine SurfDestruct
!
    Subroutine SurfSaveSolutionConfiguration(self, controlVariables, t0)

        implicit none

        class(SurfaceMesh_t)                                :: self
        type(FTValueDictionary), intent(in)                 :: controlVariables
        real(kind=RP), intent(in)                           :: t0

!       ---------------
!       Local variables
!       ---------------
        real(kind=RP)                                       :: dtSave, dtFW, dtSf
        logical                                             :: saveFWH, saveSurf

!       Check if is activated
!       ------------------------
        if (.not. self % active) then
            self % autosave = Autosave_t(.FALSE.,.FALSE.,huge(1),huge(1.0_RP),nextAutosaveTime=huge(1.0_RP),mode=AUTOSAVE_UNDEFINED)
            return
        end if
!
        saveFWH = controlVariables % containsKey("acoustic save timestep")
        saveSurf = controlVariables % containsKey("surface save timestep")
!
        if (saveSurf .and. saveFWH) then
            write(STD_OUT,'(A)') "Warning: both acoustic and surface save timestep have been specified, the minimum value will be set"
        endif
!
!       configure save type, used for update, write and save file
!       --------------------------
        if (saveSurf .or. saveFWH) then
            dtSf = controlVariables % getValueOrDefault("acoustic save timestep", huge(1.0_RP))
            dtFW = controlVariables % getValueOrDefault("surface save timestep", huge(1.0_RP))
            dtSave = min(dtSf, dtFW)
            self % autosave = Autosave_t(.FALSE.,.TRUE.,huge(1),dtSave,nextAutosaveTime=dtSave+t0,mode=AUTOSAVE_BY_TIME)
        else
            !Save each time step
            write(STD_OUT,'(A)') "Warning: neither acoustic nor surface save timestep have been specified, it will be saved each timestep"
            self % autosave = Autosave_t(.FALSE.,.TRUE.,1,huge(1.0_RP),nextAutosaveTime=huge(1.0_RP),mode=AUTOSAVE_BY_ITERATION)
        end if
    End Subroutine SurfSaveSolutionConfiguration
!
    Subroutine SurfSaveAllMesh(self, mesh, initial_iteration, controlVariables)
        implicit none
!
        class(SurfaceMesh_t)                                :: self
        type(HexMesh), intent(in)                           :: mesh
        integer, intent(in)                                 :: initial_iteration
        type(FTValueDictionary), intent(in)                 :: controlVariables

        !local variables
        integer                                             :: i, nf
        character(len=LINE_LENGTH)                          :: FinalName      !  Final name for particular file
        logical                                             :: saveFWH

        if (.not. self % active) return
        saveFWH = controlVariables % logicalValueForKey("acoustic surface mesh save") .or. self % mergeFWHandBC
        do i = 1, self % numberOfSurfaces
            if (.not. self % surfaceActive(i)) cycle
            !skip fwh if not requested
            if ( (self % surfaceTypes(i) .eq. SURFACE_TYPE_FWH) .and. (.not. saveFWH) ) cycle
            !skip acoustic surface if not requested
            if ( (self % surfaceTypes(i) .eq. SURFACE_TYPE_ACOUSTICS) .and. (.not. saveFWH) ) cycle
            write(FinalName,'(2A,I10.10)') trim(self % file_names(i)),'_',initial_iteration
            nf = self % zones(i) % no_of_faces
            call SurfaceSaveMesh(self % zones(i), mesh, FinalName, self % totalFaces(i), self % globalFid(1:nf,i), self % faceOffset(1:nf,i))
        end do
!
    End Subroutine SurfSaveAllMesh
!
    Subroutine SurfSaveAllSolution(self, mesh, iter, time, controlVariables)
        implicit none
!
        class(SurfaceMesh_t)                                :: self
        type(HexMesh), intent(in)                           :: mesh
        integer, intent(in)                                 :: iter
        real(kind=RP), intent(in)                           :: time              !< Simu time
        type(FTValueDictionary), intent(in)                 :: controlVariables

        !local variables
        integer                                             :: i, nf
        character(len=LINE_LENGTH)                          :: FinalName      !  Final name for particular file
        logical                                             :: saveFWH
        integer, dimension(:), allocatable                  :: elemSide
        logical                                             :: saveUt, saveTurb

        if (.not. self % active) return
        saveFWH = controlVariables % logicalValueForKey("acoustic solution save") .or. self % mergeFWHandBC
        do i = 1, self % numberOfSurfaces
            saveUt = .false.
            saveTurb = .false.
            if (.not. self % surfaceActive(i)) cycle
            !skip fwh if not requested
            if ( (self % surfaceTypes(i) .eq. SURFACE_TYPE_FWH) .and. (.not. saveFWH) ) cycle
            !skip acoustic surface if not requested
            if ( (self % surfaceTypes(i) .eq. SURFACE_TYPE_ACOUSTICS) .and. (.not. saveFWH) ) cycle
            write(FinalName,'(2A,I10.10,A)') trim(self % file_names(i)),'_',iter,'.surf.hsol'
            nf = self % zones(i) % no_of_faces
            safedeallocate(elemSide)
            allocate(elemSide(nf))
            if ((self % surfaceTypes(i) .eq. SURFACE_TYPE_FWH) .or. (self % surfaceTypes(i) .eq. SURFACE_TYPE_ACOUSTICS)) then
                elemSide = self % elementSide(:,1)
            else
                ! for slices and BC use always the first element of the face, this can be changed if needed
                elemSide = 1
            end if
            ! save ut or turb if requested only for BC or FWH merged with BC
            if ( (self % surfaceTypes(i) .eq. SURFACE_TYPE_BC) .or. &
                 (self % surfaceTypes(i) .eq. SURFACE_TYPE_FWH .and. self % mergeFWHandBC) ) then
                if (self%isNoSlip(i)) then
                    saveUt = self % saveUt
                    saveTurb = self % saveTurb
                end if
            end if
            call SurfaceSaveSolution(self % zones(i), mesh, time, iter, FinalName, self % totalFaces(i), &
                                 self % globalFid(1:nf,i), self % faceOffset(1:nf,i), elemSide, &
                                 self % surfaceTypes(i),self % saveGradients, &
                                 self % mergeFWHandBC, saveUt, saveTurb)
        end do
!
    End Subroutine SurfSaveAllSolution
!
    Subroutine SurfLoadSolution(self, fileName, mesh, hasGradients)
!     *******************************************************************
!        Only the FWH surface is loaded, for postprocess calculations
!     *******************************************************************
!
        implicit none
        class(SurfaceMesh_t)                                 :: self
        character(len=*), intent(in)                         :: fileName
        class (HexMesh), intent(inout)                       :: mesh
        logical, intent(in)                                  :: hasGradients

        !local variables
        integer                                             :: nf
!
!       Check if is activated
!       ------------------------
        if (.not. self % active) return
!
        if (self % surfaceTypes(FWH_POSITION) .ne. SURFACE_TYPE_FWH) then
            error stop "only fwh surface can load a solution into the solver"
        end if
!
        nf = self % zones(FWH_POSITION) % no_of_faces
        call SurfaceLoadSolution(self % zones(FWH_POSITION), mesh, fileName, self % globalFid(1:nf,FWH_POSITION), self % faceOffset(1:nf,FWH_POSITION), &
                                 self % elementSide(:,1), hasGradients)
!
    End Subroutine SurfLoadSolution
!
    Subroutine createFWHSurface(self, controlVariables, mesh, solution_file, elementSide)
        use FileReadingUtilities, only: getCharArrayFromString
        use Utilities,            only: toLower
#ifdef _HAS_MPI_
        use mpi
#endif
        implicit none

        class(SurfaceMesh_t)                                :: self
        type(FTValueDictionary), intent(in)                 :: controlVariables
        type (HexMesh), intent(in)                          :: mesh
        character(len=LINE_LENGTH), intent(in)              :: solution_file
        integer, dimension(:), allocatable, intent(out)     :: elementSide

        !local variables
        integer                                             :: i
        integer                                             :: no_of_zones, no_of_face_i, ierr, no_of_faces, numberOfFacesTotal
        integer, dimension(:), allocatable                  :: facesIDs, faces_per_zone, zonesIDs
        character(len=LINE_LENGTH)                          :: zones_str, zones_str2, surface_file, fileName, zoneName
        character(len=LINE_LENGTH), allocatable             :: zones_names(:), zones_temp(:), zones_temp2(:)
        logical                                             :: isNoSlip
        logical                                             :: surfaceHasFacesInPartitions

        if (controlVariables % containsKey("acoustic solid surface")) then
            zones_str = controlVariables % stringValueForKey("acoustic solid surface", LINE_LENGTH)
            zones_str2 = controlVariables % stringValueForKey("acoustic solid surface cont", LINE_LENGTH)
            call toLower(zones_str)
            call toLower(zones_str2)
            call getCharArrayFromString(zones_str, LINE_LENGTH, zones_temp)
            if (zones_str2 .ne. "") then
                no_of_zones = size(zones_temp)
                call getCharArrayFromString(zones_str2, LINE_LENGTH, zones_temp2)
                no_of_zones = no_of_zones + size(zones_temp2)
                allocate(zones_names(no_of_zones))
                zones_names(1:size(zones_temp)) = zones_temp
                zones_names(size(zones_temp)+1:no_of_zones) = zones_temp2
                safedeallocate(zones_temp)
                safedeallocate(zones_temp2)
            else
                no_of_zones = size(zones_temp)
                allocate(zones_names(no_of_zones))
                zones_names = zones_temp
                safedeallocate(zones_temp)
            end if 

    !       Get the zones ids from the mesh and for each, the number of faces
    !       --------------------------
            allocate( faces_per_zone(no_of_zones), zonesIDs(no_of_zones) )
            do i = 1, no_of_zones
                zonesIDs(i) = getZoneID(zones_names(i), mesh)
                if (zonesIDs(i) .eq. -1) then
                    write(*,'(A,A,A)') "Warning: acoustic surface ", trim(zones_names(i)), " not found in the mesh, will be ignored"
                    faces_per_zone(i) = 0
                else
                    faces_per_zone(i) = size(mesh % zones(zonesIDs(i)) % faces)
                end if
            end do 

    !       Get the faces Ids of all zones as a single array
    !       --------------------------
            allocate( facesIDs(SUM(faces_per_zone)) , elementSide(SUM(faces_per_zone)))
            no_of_face_i = 1
            do i = 1, no_of_zones
                if (zonesIDs(i) .eq. -1) cycle
                facesIDs(no_of_face_i:no_of_face_i+faces_per_zone(i)-1) = mesh % zones(zonesIDs(i)) % faces
                no_of_face_i = no_of_face_i + faces_per_zone(i) 
            end do 
            ! the side is always 1 since is at a face of a boundary
            ! eSides = 1
            elementSide = 1
            isNoSlip = .true.

            deallocate(zonesIDs, zones_names) 
        elseif (controlVariables % containsKey("acoustic surface file")) then
            allocate( faces_per_zone(1) )
            surface_file = controlVariables % stringValueForKey("acoustic surface file", LINE_LENGTH)
            call SurfaceLoadSurfaceFromFile(mesh, surface_file, facesIDs, faces_per_zone(1), elementSide)
            isNoSlip = .false.
        else
            error stop "acoustic surface for integration is not defined"
        end if

        no_of_faces = sum(faces_per_zone)

        if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
            call mpi_allreduce(no_of_faces, numberOfFacesTotal, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        else
            numberOfFacesTotal = no_of_faces
        end if
        surfaceHasFacesInPartitions = numberOfFacesTotal .gt. 0
!
!       now create the zone and set type, name and flag
!       --------------------------------------------------
        write(fileName,'(A,A)') trim(solution_file),'_fwh'
        zoneName = "FWH_Surface"
        call createSingleSurface(self, FWH_POSITION, SURFACE_TYPE_FWH, fileName, zoneName, no_of_faces, facesIDs, isNoSlip, surfaceHasFacesInPartitions)
!         
    End Subroutine createFWHSurface

    Subroutine createAcousticSurface(self, controlVariables, mesh, solution_file, elementSide, surfaceCont)
        use FileReadingUtilities, only: getCharArrayFromString
        use Utilities,            only: toLower
#ifdef _HAS_MPI_
        use mpi
#endif
        implicit none

        class(SurfaceMesh_t)                                :: self
        type(FTValueDictionary), intent(in)                 :: controlVariables
        type (HexMesh), intent(in)                          :: mesh
        character(len=LINE_LENGTH), intent(in)              :: solution_file
        integer, dimension(:), allocatable, intent(out)     :: elementSide
        integer, intent(in)                                 :: surfaceCont

        !local variables
        integer                                             :: i
        integer                                             :: no_of_zones, no_of_face_i, ierr, no_of_faces, numberOfFacesTotal
        integer, dimension(:), allocatable                  :: facesIDs, faces_per_zone, zonesIDs
        character(len=LINE_LENGTH)                          :: zones_str, zones_str2, surface_file, fileName, zoneName
        character(len=LINE_LENGTH), allocatable             :: zones_names(:), zones_temp(:), zones_temp2(:)
        logical                                             :: isNoSlip
        logical                                             :: surfaceHasFacesInPartitions

        if (controlVariables % containsKey("acoustic surface file")) then
            allocate( faces_per_zone(1) )
            surface_file = controlVariables % stringValueForKey("acoustic surface file", LINE_LENGTH)
            call SurfaceLoadSurfaceFromFile(mesh, surface_file, facesIDs, faces_per_zone(1), elementSide)
            isNoSlip = .false.
        else
            error stop "acoustic surface for integration is not defined"
        end if

        no_of_faces = sum(faces_per_zone)

        if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
            call mpi_allreduce(no_of_faces, numberOfFacesTotal, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        else
            numberOfFacesTotal = no_of_faces
        end if
        surfaceHasFacesInPartitions = numberOfFacesTotal .gt. 0
!
!       now create the zone and set type, name and flag
!       --------------------------------------------------
        write(fileName,'(A,A)') trim(solution_file),'_acoustic'
        zoneName = "virtual_acoustic_surface"
        call createSingleSurface(self, surfaceCont, SURFACE_TYPE_ACOUSTICS, fileName, zoneName, no_of_faces, facesIDs, isNoSlip, surfaceHasFacesInPartitions)
!         
    End Subroutine createAcousticSurface
!         
    Subroutine createSingleSurface(self, surface_index, surface_type, file_name, zone_name, no_of_faces, facesIDs, isNoSlip, has_faces)
        implicit none

        class(SurfaceMesh_t)                                :: self
        integer, intent(in)                                 :: surface_index, surface_type, no_of_faces
        character(len=LINE_LENGTH), intent(in)              :: file_name, zone_name
        integer, dimension(:), intent(in)                   :: facesIDs
        logical, intent(in)                                 :: isNoSlip
        logical, intent(in)                                 :: has_faces

        self % surfaceTypes(surface_index) = surface_type
        self % file_names(surface_index) = trim(file_name)
        ! not set active nor create zone if there are no faces
        if (.not. has_faces) return
        call self % zones(surface_index) % CreateFictitious(-1, trim(zone_name), no_of_faces, facesIDs)
        self % surfaceActive(surface_index) = .true.
        self % isNoSlip(surface_index) = isNoSlip

    End Subroutine createSingleSurface
!         
!/////////////////////////////////////////////////////////////////////////////////////////////
!           SINGLE SURFACE (ZONE) PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////////////////////////
!         
   Subroutine SurfaceSaveSolution(surface_zone, mesh, time, iter, name, no_of_faces, fGlobID, faceOffset, eSides, surface_type, &
                                  saveGradients, saveBCandFWH, saveUt, saveTurb)

!     *******************************************************************
!        This subroutine saves the solution from the face storage to a binary file
!     *******************************************************************
!
      use FaceClass
      use SolutionFile
      use fluiddata
      use PhysicsStorage

      implicit none

      class (Zone_t), intent(in)                           :: surface_zone
      class (HexMesh), intent(in), target                  :: mesh
      real(kind=RP), intent(in)                            :: time
      integer,intent(in)                                   :: iter, no_of_faces, surface_type
      character(len=*), intent(in)                         :: name
      integer, dimension(:), intent(in)                    :: fGlobID, faceOffset
      integer, dimension(:), intent(in)                    :: eSides
      logical                                              :: saveGradients, saveBCandFWH, saveUt, saveTurb

      ! local variables
      integer                                              :: zoneFaceID, meshFaceID, solution_type
      integer                                              :: ierr
      integer, dimension(6)                                :: meshFaceIDs
      class(Face), pointer                                 :: f
      real(kind=RP), dimension(:,:,:), allocatable         :: Q
      real(kind=RP), dimension(NDIM)                       :: x
      integer                                              :: Nx,Ny
      integer                                              :: fid, padding
      integer(kind=AddrInt)                                :: pos
      real(kind=RP)                                        :: refs(NO_OF_SAVED_REFS) 
      logical                                              :: saveQdot
#if (defined(CAHNHILLIARD)) && (!defined(MULTIPHASE))
      logical                                              :: computeGradients = .true.
#endif

!
!     Gather reference quantities
!     ---------------------------
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
#else
         refs = 0.0_RP
#endif

      saveQdot = .false.
      select case (surface_type)
      case (SURFACE_TYPE_FWH)
          saveQdot = .true.
          solution_type = ZONE_SOLUTION_AND_DOT_FILE
          if (saveGradients .and. saveBCandFWH) then
              padding = NCONS*2 + NGRAD*3
          else
              padding = NCONS*2
          end if 
      case (SURFACE_TYPE_BC)
          ! more variables such as friction velocity, mu, and height of first node will be saved in the future
          solution_type = ZONE_SOLUTION_FILE
          if (saveGradients)  then
              padding = NCONS + NGRAD*3
          else
              padding = NCONS
          end if
      case (SURFACE_TYPE_SLICE)
          ! only solution and gradients will be saved for slices
          solution_type = ZONE_SOLUTION_FILE
          if (saveGradients)  then
              padding = NCONS + NGRAD*3
          else
              padding = NCONS
          end if
      case (SURFACE_TYPE_ACOUSTICS)
        solution_type = ZONE_SOLUTION_FILE
          if (saveGradients)  then
              padding = NCONS + NGRAD*3
          else
              padding = NCONS
          end if
      end select

      if (saveUt) padding = padding + 1
      ! save mu_NS and y, for y+ value calc
      if (saveTurb) padding = padding + 2
!
!     Create new file
!     ---------------
      call CreateNewSolutionFile(trim(name), solution_type, mesh % nodeType, &
                                    no_of_faces, iter, time, refs)
!
!     Write arrays
!     ------------
      fID = putSolutionFileInWriteDataMode(trim(name))
!     Loop the zone to get faces
!     ---------------------------------------
      do zoneFaceID = 1, surface_zone % no_of_faces
          meshFaceID = surface_zone % faces(zoneFaceID)
          f => mesh % faces(meshFaceID)

          Nx = f % Nf(1)
          Ny = f % Nf(2)

          allocate (Q(1:NCONS,0:Nx,0:Ny))
          Q(1:NCONS,:,:)  = f % storage(eSides(zoneFaceID)) % Q(:,:,:)

          ! 4 integers are written: number of dimension, and 3 value of the dimensions
          pos = POS_INIT_DATA + (fGlobID(zoneFaceID)-1)*4_AddrInt*SIZEOF_INT + padding * faceOffset(zoneFaceID) * SIZEOF_RP
          call writeArray(fid, Q, position=pos)

          if (saveQdot) then
              Q(1:NCONS,:,:)  = f % storage(eSides(zoneFaceID)) % Qdot(:,:,:)
              write(fid) Q
          end if
          deallocate(Q)

          ! taken and adapted from saveSolution of HexMesh
          if ( saveGradients .and. computeGradients ) then

               allocate(Q(1:NGRAD,0:Nx,0:Ny))

#ifdef FLOW
               Q(1:NGRAD,:,:) = f % storage(eSides(zoneFaceID)) % U_x
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               Q(1:NGRAD,:,:) = f % storage(eSides(zoneFaceID)) % c_x
#endif
               write(fid) Q

#ifdef FLOW
               Q(1:NGRAD,:,:) = f % storage(eSides(zoneFaceID)) % U_y
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               Q(1:NGRAD,:,:) = f % storage(eSides(zoneFaceID)) % c_y
#endif
               write(fid) Q

#ifdef FLOW
               Q(1:NGRAD,:,:) = f % storage(eSides(zoneFaceID)) % U_z
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               Q(1:NGRAD,:,:) = f % storage(eSides(zoneFaceID)) % c_z
#endif
               write(fid) Q

          if (zoneFaceID .eq. 1) then
          end if 
               deallocate(Q)
          end if
          if (saveUt) then
#if defined(NAVIERSTOKES)
               allocate(Q(1,0:Nx,0:Ny))
               Q(1,:,:)= f % storage(eSides(zoneFaceID)) % u_tau_NS(:,:)
               write(fid) Q
               deallocate(Q)
#endif
          end if 
          if (saveTurb) then
#if defined(NAVIERSTOKES)
               allocate(Q(1,0:Nx,0:Ny))
               Q(1,:,:)= f % storage(1) % mu_NS(1,:,:)
               write(fid) Q
               Q(1,:,:)= f % storage(1) % wallNodeDistance(:,:)
               write(fid) Q
               deallocate(Q)
#endif
          end if 
          safedeallocate(Q)
      end do

     close(fid)
!
!    Close the file
!    --------------
     call SealSolutionFile(trim(name))
	 nullify(f)
!
   End Subroutine SurfaceSaveSolution
!         
   Subroutine SurfaceLoadSolution(surface_zone, mesh, fileName, fGlobID, faceOffset, eSides, hasGradients)

      use SolutionFile
      use PhysicsStorage

      implicit none

      class (Zone_t), intent(in)                           :: surface_zone
      class (HexMesh), intent(inout)                       :: mesh
      character(len=*), intent(in)                         :: fileName
      integer, dimension(:), intent(in)                    :: fGlobID, faceOffset
      integer, dimension(:), intent(in)                    :: eSides
      logical, intent(in)                                  :: hasGradients

      ! local variables
      integer                                              :: zoneFaceID, meshFaceID, fileType
      real(kind=RP), dimension(:,:,:), allocatable         :: QF
      integer                                              :: Nx,Ny
      integer                                              :: fID, pos, padding
      integer                                              :: arrayRank, Neq, Npx, Npy

!     get type and num of vars
!     ------------------------
      fileType = getSolutionFileType(trim(fileName))
      select case (fileType)
      case (ZONE_SOLUTION_AND_DOT_FILE)

      case default
          write(STD_OUT,'(A,I0,A,I0)') "Error, hsol file expected is ", ZONE_SOLUTION_AND_DOT_FILE, ", got: ", fileType
          errorMessage(STD_OUT)
      end select

      padding = NCONS*2
      if (hasGradients) padding = padding + NGRAD*3
!
!     Read elements data
!     ------------------
      fID = putSolutionFileInReadDataMode(trim(fileName))

!
!     Loop the zone to get faces and elements
!     ---------------------------------------
      do zoneFaceID = 1, surface_zone % no_of_faces
          meshFaceID = surface_zone % faces(zoneFaceID)
          ! 4 integers were written: number of dimension, and 3 value of the dimensions
          pos = POS_INIT_DATA + (fGlobID(zoneFaceID)-1)*4*SIZEOF_INT + padding * faceOffset(zoneFaceID) * SIZEOF_RP
          associate(f => mesh % faces(meshFaceID))
              Nx = f % Nf(1)
              Ny = f % Nf(2)
!             verify dimensions of each row
              read(fID, pos=pos) arrayRank
              read(fID) Neq, Npx, Npy
              if (     ((Npx-1) .ne. Nx) &
                  .or. ((Npy-1) .ne. Ny) &
                  .or. (Neq     .ne. NCONS ) ) then
                  write(STD_OUT,'(A,I0,A)') "Error reading fwh file: wrong dimension for face "&
                      ,meshFaceID,"."

                  write(STD_OUT,'(A,I0,A,I0,A)') "Face dimensions: ", Nx, &
                      " ,", Ny, "."

                  write(STD_OUT,'(A,I0,A,I0,A)') "File dimensions: ", Npx -1, &
                      " ,", Npy-1, "."
                  errorMessage(STD_OUT)
                  error stop
              end if

              allocate(QF(1:NCONS,0:Nx,0:Ny))
              read(fID) QF
              f % storage(eSides(zoneFaceID)) % Q = QF
              read(fID) QF
              f % storage(eSides(zoneFaceID)) % Qdot = QF
              safedeallocate(QF)
          end associate
      end do

!     Close the file
!     --------------
      close(fID)
!
   End Subroutine SurfaceLoadSolution
!
   Subroutine SurfacePrepareForIO(surface_zone, mesh, totalNumberOfFaces, globalFaceID, faceOffset)

!     *******************************************************************
!        This subroutine creates the arrays necessary for the face binary file
!     *******************************************************************
!
      use FaceClass
#ifdef _HAS_MPI_
        use mpi
#endif
      implicit none
      class (Zone_t), intent(in)                           :: surface_zone
      class (HexMesh), intent(in), target                  :: mesh
      integer,intent(in)                                   :: totalNumberOfFaces
      integer, dimension(surface_zone % no_of_faces), intent(out)          :: globalFaceID, faceOffset

      ! local variables
      integer                                              :: zoneFaceID, meshFaceID, eID, i
      integer, dimension(:), allocatable                   :: gFid, facesSizes, allfacesSizes, allFacesOffset
      integer, dimension(:,:), allocatable                 :: zoneInfoArray
      ! integer, dimension(:), allocatable                 :: zoneInfoArray
      integer                                              :: ierr, fID
      class(Face), pointer                                 :: faces(:)
      integer, dimension(MPI_Process % nProcs)             :: no_of_faces_p, displs
      integer, dimension(1)                                :: idInGlobal
      integer, dimension(:), allocatable                   :: elementsID_p

      faces => mesh % faces

!     *******************************************************************
!     Get the globalFaceID
!     *******************************************************************

      allocate(gFid(totalNumberOfFaces))
      allocate(zoneInfoArray(totalNumberOfFaces,2))

      if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
          call mpi_gather(surface_zone % no_of_faces,1,MPI_INT,no_of_faces_p,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

      if (MPI_Process % isRoot) then
          displs=0
          do i = 2, MPI_Process % nProcs 
              displs(i) = displs(i-1) + no_of_faces_p(i-1)
          end do
      end if

      allocate(elementsID_p(surface_zone%no_of_faces))

			do i=1, surface_zone % no_of_faces
					eID = faces(surface_zone % faces(i)) % elementIDs(1)
					if (eID == 0) then
						eID = faces(surface_zone % faces(i)) % elementIDs(2)
					end if
					elementsID_p(i) = mesh % elements(eID) % globID
			end do

      ! get the global element ID and the face ID as a single 2D array for the root process, will be used to sort
      call mpi_gatherv(elementsID_p, surface_zone % no_of_faces,MPI_INT, zoneInfoArray(:,1), no_of_faces_p, displs, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      call mpi_gatherv(surface_zone % faces, surface_zone % no_of_faces,MPI_INT, &
                                     zoneInfoArray(:,2), no_of_faces_p, displs, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      ! get the sorted array
      if (MPI_Process % isRoot) then
          gFid = getGlobalFaceIDs(zoneInfoArray, totalNumberOfFaces)
      end if

      ! distribute to all partitions
      call mpi_scatterv(gFid, no_of_faces_p, displs, MPI_INT, globalFaceID, surface_zone % no_of_faces, MPI_INT, 0, MPI_COMM_WORLD, ierr)
#endif
      else
          zoneInfoArray(:,1) = mesh % elements(faces(surface_zone % faces) % elementIDs(1)) % globID
          zoneInfoArray(:,2) = surface_zone % faces
          ! get the sorted array
          globalFaceID = getGlobalFaceIDs(zoneInfoArray, totalNumberOfFaces)
      end if

!     Free memory
!     -----------
      deallocate(gFid, zoneInfoArray)

!     *******************************************************************
!     Get the faceOffset, similar to HexMesh_PrepareForIO, but for faces
!     *******************************************************************

!     Get each face storage size
!     ---------------------
      allocate(facesSizes(totalNumberOfFaces), allfacesSizes(totalNumberOfFaces))

      facesSizes = 0 ! default to use allreduce
      do zoneFaceID = 1, surface_zone % no_of_faces
          ! globalFaceID index
          fID = globalFaceID(zoneFaceID)
          meshFaceID = surface_zone % faces(zoneFaceID)
          facesSizes(fID) = ( faces(meshFaceID) % Nf(1) +1 ) * ( faces(meshFaceID) % Nf(2) +1 )
      end do

      allfacesSizes = 0
      if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
          call mpi_allreduce(facesSizes, allfacesSizes, totalNumberOfFaces, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      else
          allfacesSizes = facesSizes
      end if
    
!     Get all faces offset: the accumulation of allfacesSizes
!     -----------------------
      allocate(allFacesOffset(totalNumberOfFaces))

      allFacesOffset(1) = 0
      do fID = 2, totalNumberOfFaces
          allFacesOffset(fID) = allFacesOffset(fID-1) + allfacesSizes(fID-1)
      end do

!     Assign the results to partitions' array
!     ----------------------------------
      do zoneFaceID = 1, surface_zone % no_of_faces
          fID = globalFaceID(zoneFaceID)
          faceOffset(zoneFaceID) = allFacesOffset(fID)
      end do

!     Free memory
!     -----------
      deallocate(facesSizes, allfacesSizes, allFacesOffset)
	  nullify(faces)
!
   End Subroutine SurfacePrepareForIO
!
   Function getGlobalFaceIDs(zoneInfoArray, N) result(gID)
!     *******************************************************************
!        This function gets a unique identifier of each face of the surface_zone,
!        in a unique order needed for I/O
!     *******************************************************************

      use Utilities, only: QsortWithFriend
      implicit none

      integer,dimension(N,2),intent(in)                    :: zoneInfoArray
      integer,intent(in)                                   :: N
      integer,dimension(N)                                 :: giD

      ! local variables
      integer,dimension(N)                                 :: originalIndex, orderedIndex, toOrderArray
      integer                                              :: bigInt, i, maxF, nDigits

      maxF = maxval(zoneInfoArray(:,2))
      nDigits = 0
      ! get number of digits of maxF
      do while (maxF .ne. 0)
          maxF = maxF / 10
          nDigits = nDigits + 1
      end do
      ! convert the two arrays as a one of integers that can be sort, first by the globID of the element and the by the faceID
      bigInt = 10 ** nDigits
      toOrderArray = bigInt * zoneInfoArray(:,1) + zoneInfoArray(:,2)
      ! create simple arrays of indexes
      originalIndex = [(i, i=1,N)]
      orderedIndex = [(i, i=1,N)]
      call QsortWithFriend(toOrderArray,originalIndex)
      ! get the indexes of orders`
      call QsortWithFriend(originalIndex, orderedIndex)
      gID = orderedIndex
!
   End Function getGlobalFaceIDs
!
   Subroutine SurfaceLoadSurfaceFromFile(mesh, surface_file, facesIDs, numberOfFaces, eSides)

!     *******************************************************************
!        This subroutine reads the faces of the surface from a text file
!     *******************************************************************
!
      use ElementConnectivityDefinitions, only: FACES_PER_ELEMENT, NODES_PER_FACE
#ifdef _HAS_MPI_
        use mpi
#endif
      implicit none

      class(HexMesh), intent(in)                          :: mesh
      character(len=LINE_LENGTH), intent(in)              :: surface_file
      integer, dimension(:), allocatable, intent(out)     :: facesIDs, eSides
      integer, intent(out)                                :: numberOfFaces

      ! local variables
      integer                                             :: fd       ! File unit
      integer                                             :: i, j, k
      integer                                             :: ierr
      integer                                             :: allNumberOfFaces, eID, sumFaces
      integer, dimension(:), allocatable                  :: eIDs, allGeIDs, allfIDs
      integer, dimension(:), allocatable                  :: fIDsTemp, eIDsTemp
      integer, dimension(:,:), allocatable                :: nodeIDs, allNodeIDs, nodeIDsTemp
      integer, dimension(NODES_PER_FACE)                  :: faceNodeID
      logical                                             :: hasMPIExtraInfo
         

!       Only root open and read the file
!       ------------------
        if (MPI_Process % isRoot) then
          open(newunit = fd, file = surface_file )   
          read(fd,*) allNumberOfFaces

          allocate( allfIDs(allNumberOfFaces), allGeIDs(allNumberOfFaces) )

          read(fd,*) hasMPIExtraInfo

          if (hasMPIExtraInfo) then
            allocate(allNodeIDs(4,allNumberOfFaces))
            do i = 1, allNumberOfFaces
              read(fd,*) allGeIDs(i), allfIDs(i), allNodeIDs(:,i)
            end do
          else
            if (MPI_Process % doMPIAction) then
                write(STD_OUT,'(A)') "Error, the surface file does not have the information necessary to create it while running with MPI"
                call exit(99)
            end if 
            do i = 1, allNumberOfFaces
              read(fd,*) allGeIDs(i), allfIDs(i)
            end do
          end if
          close(unit=fd)
        end if

!       For MPI, check the nodeIDs and compare it to each face of the element to find the right one in the partition
!       ------------------
        if ( (MPI_Process % doMPIAction) ) then

!         Broadcast from root, since other process didn't read the file
!         ------------------
#ifdef _HAS_MPI_
          call mpi_barrier(MPI_COMM_WORLD, ierr)
          call mpi_Bcast(allNumberOfFaces, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)

          if (.not. MPI_Process % isRoot) allocate( allfIDs(allNumberOfFaces), allGeIDs(allNumberOfFaces), allNodeIDs(4,allNumberOfFaces) )
          ! call mpi_barrier(MPI_COMM_WORLD, ierr)
          call mpi_Bcast(allGeIDs, allNumberOfFaces, MPI_INT, 0, MPI_COMM_WORLD, ierr)

          call mpi_Bcast(allNodeIDs, 4*allNumberOfFaces, MPI_INT, 0, MPI_COMM_WORLD, ierr)
#endif

!         First compare the globaleID read from the element in the partition and save the normal (not global) eID, and the
!         correspondent read nodeIDs
!         ------------------
          allocate( nodeIDs(4,allNumberOfFaces), eIDs(allNumberOfFaces) )
          j = 0
          do i = 1, allNumberOfFaces
            elems_loop:do eID = 1, mesh%no_of_elements
              if (mesh % elements(eID) % globID .eq. allGeIDs(i)) then
                j = j + 1
                eIDs(j) = eID
                nodeIDs(:,j) = allNodeIDs(:,i)
                exit elems_loop
              end if
            end do elems_loop
          end do
          numberOfFaces = j
          allocate( nodeIDsTemp(4,numberOfFaces), eIDsTemp(numberOfFaces), facesIDs(numberOfFaces), eSides(numberOfFaces) )
          eIDsTemp(:) = eIDs(1:numberOfFaces)
          nodeIDsTemp(:,:) = nodeIDs(:,1:numberOfFaces)
          call move_alloc(eIDsTemp, eIDs)
          call move_alloc(nodeIDsTemp, nodeIDs)

!       Check that all the elements are found in the mesh partitions
!       ------------------
#ifdef _HAS_MPI_
         call mpi_allreduce(numberOfFaces, sumFaces, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
          if (sumFaces .ne. allNumberOfFaces .and. MPI_Process % isRoot) then
              print *, "Error:, not all the elements of the surface file where found in the mesh. File elements: ", allNumberOfFaces, &
              "elements found: ", sumFaces
            call exit(99)
          end if 

!         Now get the actual facesID of the partition for each element based on the nodeIDs
!         ------------------
          element_face_loop:do i = 1, numberOfFaces
            associate( e => mesh % elements(eIDs(i)) )
              do j = 1, FACES_PER_ELEMENT
                do k = 1, NODES_PER_FACE
                  faceNodeID(k) = mesh % nodes(mesh % faces(e % faceIDs(j)) % nodeIDs(k)) % globID
                end do
                if (all( nodeIDs(:,i) .eq. faceNodeID)) then
                  facesIDs(i) = e % faceIDs(j)
                  eSides(i) = e % faceSide(j)
                  cycle element_face_loop
                end if
              end do
                ! if code arrives here, the face was not found
              print *, "Error:, for surface file the face was not found for node ids: ", nodeIDs(:,i), &
                    "and element ID: ", eIDs(i)
              call exit(99)
            end associate
          end do element_face_loop

!       For non MPI, just pass to normal arrays and check if the faces correspond to the element, and get side of face
!       ------------------
        else
          numberOfFaces = allNumberOfFaces
          allocate( facesIDs(numberOfFaces), eIDs(numberOfFaces), eSides(numberOfFaces) )
          facesIDs = allfIDs
          eIDs = allGeIDs
          deallocate(allGeIDs, allfIDs)
          eSides = 0
          do i = 1, numberOfFaces
              ! if findloc is supported by the compiler use this line and comment the if
              ! eSides(i) = findloc(mesh % faces(facesIDs(i)) % elementIDs, eIDs(i), dim=1)
              if ( mesh % faces(facesIDs(i)) % elementIDs(1) .eq. eIDs(i) ) then
                eSides(i) = 1
              elseif ( mesh%faces(facesIDs(i))%elementIDs(2) .eq. eIDs(i) ) then
                eSides(i) = 2
              end if
              if (eSides(i) .eq. 0) then
                print *, "Error: the element ", eIDs(i), " does not correspond to the face ", mesh % faces(facesIDs(i)) % ID, &
                    ". The elements of the face are: " , mesh % faces(facesIDs(i)) % elementIDs, ". The faces of the element are: ", mesh % elements(eIDs(i)) % faceIDs
                call exit(99)
              end if 
          end do
        end if
!
   End Subroutine SurfaceLoadSurfaceFromFile
!
   Subroutine SurfaceSaveMesh( surface_zone, mesh, fileName, no_of_faces, fGlobID, faceOffset )
!
!     *******************************************************************
!        This subroutine saves the mesh from the face geometry to a binary file
!     *******************************************************************
!
      use SolutionFile
      use PhysicsStorage, only: Lref
      use FileReadingUtilities, only: removePath, getFileName
      implicit none
      class (Zone_t), intent(in)                           :: surface_zone
      class (HexMesh), intent(in)                          :: mesh
      character(len=*), intent(in)                         :: fileName
      integer,intent(in)                                   :: no_of_faces
      integer, dimension(:), intent(in)                    :: fGlobID, faceOffset

      ! local variables
      integer                                              :: zoneFaceID, meshFaceID
      integer                                              :: fid
      integer(kind=AddrInt)                                :: pos
      character(len=LINE_LENGTH)                           :: meshName
      real(kind=RP), parameter                             :: refs(NO_OF_SAVED_REFS) = 0.0_RP

!
!     Create file: it will be contained in ./MESH
!     -------------------------------------------
      meshName = "./MESH/" // trim(removePath(getFileName(fileName))) // ".surf.hmesh"
      call CreateNewSolutionFile(trim(meshName), ZONE_MESH_FILE, mesh % nodeType, &
                                    no_of_faces, 0, 0.0_RP, refs)
!
!     Write arrays
!     ------------
       fID = putSolutionFileInWriteDataMode(trim(meshName))

!     Loop the zone to get faces
!     ---------------------------------------
      do zoneFaceID = 1, surface_zone % no_of_faces
          meshFaceID = surface_zone % faces(zoneFaceID)
          associate( f => mesh % faces(meshFaceID) )
              ! 4 integers are written: number of dimension, and 3 value of the dimensions
              pos = POS_INIT_DATA + (fGlobID(zoneFaceID)-1)*4_AddrInt*SIZEOF_INT + 3_AddrInt * faceOffset(zoneFaceID) * SIZEOF_RP
              call writeArray(fid, f % geom % x(:,0:f%Nf(1),0:f%Nf(2))*Lref, position=pos)
          end associate

      end do

     close(fid)
!
!
!    Close the file
!    --------------
      call SealSolutionFile(trim(meshName))

   End Subroutine SurfaceSaveMesh
!         
!/////////////////////////////////////////////////////////////////////////////////////////////
!           HELPER PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////////////////////////
!         
    integer Function getZoneID(zone_name, mesh) result(n)

        character(len=*), intent(in)                        :: zone_name
        type(HexMesh), intent(in)                           :: mesh

        !local variables
        integer                                             :: zoneID

         n = -1
         do zoneID = 1, size(mesh % zones)
            if ( trim(mesh % zones(zoneID) % name) .eq. trim(zone_name) ) then
               n = zoneID
               exit
            end if
         end do

    End Function getZoneID
!
    Subroutine getSliceFaces(direction, coordinate, mesh, numberOfFaces, facesIDs, limitsSlice)

        implicit none
        integer, intent(in)                                 :: direction
        real(kind=RP), intent(in)                           :: coordinate
        type(HexMesh), intent(in)                           :: mesh
        integer, intent(out)                                :: numberOfFaces
        integer, dimension(:), allocatable, intent(out)     :: facesIDs
        real(kind=RP), dimension(:), intent(in)             :: limitsSlice

        ! local variables
        integer                                             :: i, j
        integer, dimension(2)                               :: dirLims
        integer, dimension(:), allocatable                  :: facesIDsTemp
        logical                                             :: useLimits


        useLimits = size(limitsSlice) .eq. 4
        allocate(facesIDsTemp(mesh % no_of_faces))
        j = 0
        do i = 1, mesh % no_of_faces
            associate( x => mesh % faces(i) % geom % x(:,:,:) )
                if ( all(abs(x(direction,:,:)-coordinate) .le. sliceTolerance) ) then
                    ! not use the face if there is no element in mesh associated with it.
                    ! weird check, but necessary for strange bug when the coordinate is 0
                    if (mesh % faces(i) % elementIDs(1) .eq. 0) cycle
                    if (useLimits) then
                        dirLims(1) = direction + 1
                        dirLims(2) = direction + 2
                        if (dirLims(1) .gt. NDIM) dirLims(1) = dirLims(1) - NDIM
                        if (dirLims(2) .gt. NDIM) dirLims(2) = dirLims(2) - NDIM
                        if ( all(x(dirLims(1),:,:) .lt. limitsSlice(1)) ) cycle
                        if ( all(x(dirLims(1),:,:) .gt. limitsSlice(2)) ) cycle
                        if ( all(x(dirLims(2),:,:) .lt. limitsSlice(3)) ) cycle
                        if ( all(x(dirLims(2),:,:) .gt. limitsSlice(4)) ) cycle
                    end if 
                    j = j + 1
                    facesIDsTemp(j) = i
                end if 
            end associate
        end do

        numberOfFaces = j
        allocate(facesIDs(numberOfFaces))
        facesIDs = facesIDsTemp(1:numberOfFaces)
!
    End Subroutine getSliceFaces
!         
!/////////////////////////////////////////////////////////////////////////////////////////////
!           ADDITIONAL VARIABLES PROCEDURES --------------------------
!/////////////////////////////////////////////////////////////////////////////////////////////
!         
#if defined(NAVIERSTOKES)
    Subroutine getU_tauInSurfaces(surfaces, mesh)
        use BoundaryConditions
        use Physics
        use fluiddata
        implicit none
        type(SurfaceMesh_t)                                     :: surfaces
        type(HexMesh), intent(inout)                            :: mesh

        !local variables
        integer                                                 :: surfID, faceID, i, j, meshFaceID
        real(kind=RP)                                           :: u_tau_t1, u_tau_t2
        real(kind=RP)                                           :: qq, u, v, w
        real(kind=RP), dimension(NDIM)                          :: freestream_dir


        if (.not. surfaces % active) return
        if (.not. surfaces % saveUt) return

        u  = cos(refValues % AoAtheta * PI / 180.0_RP)*cos(refValues % AoAphi * PI / 180.0_RP)
        v  = sin(refValues % AoAtheta * PI / 180.0_RP)*cos(refValues % AoAphi * PI / 180.0_RP)
        w  = sin(refValues % AoAphi * PI / 180.0_RP)

        freestream_dir = [u,v,w]

        do surfID = 1, surfaces % numberOfSurfaces
            if (.not. surfaces % surfaceActive(surfID)) cycle
            if ( (surfaces % surfaceTypes(surfID) .ne. SURFACE_TYPE_BC) .or. &
                 (surfaces % surfaceTypes(surfID) .eq. SURFACE_TYPE_FWH .and. .not. surfaces % mergeFWHandBC) ) cycle
             if (.not. surfaces % isNoSlip(surfID)) cycle
             ! save u_tau_NS in each face of the no slip bc zones
!!!$omp parallel do default(shared) private(faceID,u_tau_t1,u_tau_t2,meshFaceID,i,j)
             do faceID = 1, surfaces % zones(surfID) % no_of_faces
                meshFaceID = surfaces % zones(surfID) % faces(faceID)
                associate( Q => mesh % faces(meshFaceID) % storage(1) % Q, &
                           U_x => mesh % faces(meshFaceID) % storage(1) % U_x, &
                           U_y => mesh % faces(meshFaceID) % storage(1) % U_y, &
                           U_z => mesh % faces(meshFaceID) % storage(1) % U_z, &
                           normal => mesh % faces(meshFaceID) % geom % normal, &
                           tangent_1 => mesh % faces(meshFaceID) % geom % t1, &
                           tangent_2 => mesh % faces(meshFaceID) % geom % t2, &
                           u_tau => mesh % faces(meshFaceID) % storage(1) % u_tau_NS )
                    do j=0,mesh % faces(meshFaceID) % Nf(2)   ; do i = 0, mesh % faces(meshFaceID) % Nf(1)
                        call getFrictionVelocity(Q(:,i,j),U_x(:,i,j),U_y(:,i,j),U_z(:,i,j),normal(:,i,j),tangent_1(:,i,j),tangent_2(:,i,j), freestream_dir, u_tau(i,j))
                    end do             ; end do
                end associate

             end do
!!!$omp end parallel do

        end do

    End Subroutine getU_tauInSurfaces
!
    Subroutine getWallDistInSurfaces(surfaces, mesh)
        use BoundaryConditions
        use Physics
        implicit none
        type(SurfaceMesh_t)                                     :: surfaces
        type(HexMesh), intent(inout)                            :: mesh

        !local variables
        integer                                                 :: surfID, faceID, i, j, meshFaceID

        if (.not. surfaces % active) return
        if (.not. surfaces % saveUt) return
        do surfID = 1, surfaces % numberOfSurfaces
            if (.not. surfaces % surfaceActive(surfID)) cycle
            if ( (surfaces % surfaceTypes(surfID) .ne. SURFACE_TYPE_BC) .or. &
                 (surfaces % surfaceTypes(surfID) .eq. SURFACE_TYPE_FWH .and. .not. surfaces % mergeFWHandBC) ) cycle
             if (.not. surfaces % isNoSlip(surfID)) cycle
             ! save dwall in each face of the no slip bc zones
!!!$omp parallel
             do faceID = 1, surfaces % zones(surfID) % no_of_faces
                meshFaceID = surfaces % zones(surfID) % faces(faceID)
                associate( f => mesh % faces(meshFaceID), &
                           dw => mesh % faces(meshFaceID) % storage(1) % wallNodeDistance )
                    call getFaceWallDistance(mesh, f, dw)
                end associate
             end do
!!!$omp end parallel do

        end do

    End Subroutine getWallDistInSurfaces
!
    Subroutine getFaceWallDistance(mesh,f,dWall)
        use FaceClass
        use NodalStorageClass, only: GAUSS, GAUSSLOBATTO
        use ElementConnectivityDefinitions, only: normalAxis, FACES_PER_ELEMENT
        implicit none
        type(HexMesh), intent(in)                                   :: mesh
        type(Face), intent(in)                                      :: f
        real(kind=RP), dimension(0:f % Nf(1),0:f % Nf(2)), intent(out)  :: dWall

        ! local variables
        
        integer                                                     :: eID, efID, i, j
        integer                                                     :: normalDirection, indexArray(2), minIndex, elementIndex
        integer                                                     :: nodeStart, nodeEnd, N
        real(kind=RP), dimension(NDIM)                              :: x0, xN, xf, dWallVector
        real(kind=RP), dimension(2)                                 :: dx
        real(kind=RP), dimension(NDIM,0:f % Nf(1),0:f % Nf(2))      :: x

        eID = f % ElementIDs(1)
        ! get direction of the face respect to the element
        associate ( e => mesh % elements(eID) )
            elem_loop:do efID = 1, FACES_PER_ELEMENT
                if (e % faceIDs(efID) .eq. f % ID) then
                    normalDirection = abs(normalAxis(efID))
                    exit elem_loop
                end if
            end do elem_loop
            select case (normalDirection)
                case (1)
                    N = e % Nxyz(1)
                case (2)
                    N = e % Nxyz(2)
                case (3)
                    N = e % Nxyz(3)
                case default
                   write(STD_OUT,'(A)') "Error: normalDirection not found in axisMap"
                   errorMessage(STD_OUT)
                   error stop 
            end select
        end associate

        ! get nodeType
        select case (mesh % nodeType)
            case (GAUSS)
                nodeStart = 0
                nodeEnd = N
            case (GAUSSLOBATTO)
                nodeStart = 1
                nodeEnd = N - 1
        end select

        ! get x of 1 (or 2) and N (or N-1) node of element
        xf = f % geom % x(:,0,0)
        associate ( e => mesh % elements(eID) )
            select case (normalDirection)
                case (1)
                    x0 = e % geom % x(:,nodeStart,0,0)
                    xN = e % geom % x(:,nodeEnd,0,0)
                case (2)
                    x0 = e % geom % x(:,0,nodeStart,0)
                    xN = e % geom % x(:,0,nodeEnd,0)
                case (3)
                    x0 = e % geom % x(:,0,0,nodeStart)
                    xN = e % geom % x(:,0,0,nodeEnd)
            end select

            !compare to x of wall
            ! get index by min distance
            indexArray = [nodeStart,nodeEnd]
            dx(1) = norm2(xf-x0)
            dx(2) = norm2(xf-xN)
            minIndex = minloc(dx,dim=1)

            elementIndex = indexArray(minIndex)
            ! get all min dist of the face
            select case (normalDirection)
                case (1)
                    x(:,:,:) = e % geom % x(:,elementIndex,:,:)
                case (2)
                    x(:,:,:) = e % geom % x(:,:,elementIndex,:)
                case (3)
                    x(:,:,:) = e % geom % x(:,:,:,elementIndex)
            end select
        end associate
        do j = 0, f % Nf(2)
            do i = 0, f % Nf(1)
                dWallVector(:) = x(:,i,j) - f % geom % x(:,i,j)
                dWall(i,j) = norm2(dWallVector)
            end do
        end do

    End Subroutine getFaceWallDistance
#endif
!
End Module SurfaceMesh
