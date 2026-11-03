!
!////////////////////////////////////////////////////////////////////////////////////////
!
!     MeshInterpolationTool
!
!     Generic interpolation of HORSES3D .hsol files from an input mesh onto a different
!     output mesh. The meshes can have different elements, polynomial orders (also
!     non-conforming) and discretization nodes (Gauss / Gauss-Lobatto).
!
!     It gathers the capabilities of horsesConverter (meshInterpolation task),
!     StatsMeshInterpolation and LambVectorMeshInterpolation:
!        * The layout of the variables stored in each file is detected from the file
!          itself. Solution files (with or without gradients, sensor, LES or source terms),
!          statistics files (NS, iNS and MU) and Lamb vector files can be interpolated.
!        * The position of the output nodes inside the input mesh (mapping) is computed
!          once and reused for all the files that share the same input mesh.
!        * The point search starts at the last element found and uses a background grid
!          of bins. Points that are not found are projected onto the closest input element.
!        * Each input file may use its own input mesh (time-dependent meshes).
!
!     Control file keywords:
!
!        tool type                      = mesh interpolation
!
!        ! Input mesh: .hmesh file with the nodes of the mesh where the data is stored
!        input hmesh file name          = MESH/old.hmesh
!        input hmesh per file           = .false.   ! (optional) .true.: each file "<name>_<iter>.hsol"
!                                                   ! uses the mesh "<input hmesh path>/<input hmesh name>_<iter>.hmesh"
!
!        ! Output mesh: either a .hmesh file (the nodes and the discretization nodes are taken from it)...
!        output hmesh file name         = MESH/new.hmesh
!        ! ...or a mesh file, the polynomial order and the discretization nodes (standard keywords).
!        ! The mesh is constructed with the readers of the solver (curved elements supported), so
!        ! the "#define boundary" blocks of the mesh markers are needed (as in the solver control
!        ! file) and "reference length (m)" must be the one of the simulation (default 1)
!        mesh file name                 = MESH/new.msh
!        polynomial order               = 5           ! or "polynomial order file", as in the solver
!        discretization nodes           = Gauss-Lobatto
!
!        ! Files to interpolate (any combination of the following keywords)
!        input file name                = RESULTS/case_0000100.hsol
!        input file name 2              = RESULTS/case.stats.hsol    ! (optional) "input file name 3", ...
!        input file pattern             = RESULTS/Lamb_*.hsol        ! (optional) shell pattern
!
!        ! Optional keywords
!        interpolation field            = auto      ! auto / solution / stats / lamb vector (checks the file type)
!        output directory               = RESULTS   ! Directory of the interpolated files
!        output file suffix             = .Interpolated
!        outside point tolerance        = 0.1       ! Tolerance (reference coordinates) to accept points slightly outside an element
!        stats solver                   = ns        ! ns / ins / mu. Only needed for ambiguous statistics files
!        interpolation variable blocks  = [9,5,1]   ! Expert: number of variables of each array stored per element
!
!     For compatibility, the keywords of the previous tools are also accepted when the new
!     ones are not present: "stats mesh file name" (input hmesh), "stats qbase file name",
!     "stats sound velocity squared file name", "stats gradient sound velocity squared file name",
!     "stats Lamb vector file name", "Lamb vector dir" and "Lamb vector file name".
!
!     The interpolated files are written as <output directory>/<input name><suffix>.hsol
!     (or <input name><suffix>.stats.hsol for files named *.stats.hsol).
!
!     The input .hmesh file can also be given as an input file: the interpolated node
!     coordinates are written as <output directory>/<input name><suffix>.hmesh, which is
!     the .hmesh of the output mesh (useful to use it later as input mesh).
!
!////////////////////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
module MeshInterpolationTool
    use SMConstants
    use NodalStorageClass
    use SolutionFile, only: NO_OF_SAVED_REFS

    implicit none

    private
    public MeshInterpolation

    ! Description of the data stored in a .hsol file
    type FileLayout_t
        integer               :: fileType
        integer               :: nodeType
        integer               :: no_of_elements
        integer               :: iter
        real(kind=RP)         :: time
        real(kind=RP)         :: refs(NO_OF_SAVED_REFS)
        integer               :: padding                ! Number of variables stored per node
        integer, allocatable  :: blocks(:)              ! Number of variables of each array written per element
        integer               :: sensorBlock = 0        ! The element sensor is stored after this block (0: no sensor)
    end type FileLayout_t

    ! Data of one element. The arrays of the file are stacked in the first dimension of Q
    type ElementData_t
        real(kind=RP), allocatable :: Q(:,:,:,:)
        real(kind=RP)              :: sensor = 0.0_RP
    end type ElementData_t

    ! Position of the nodes of an output element in the input mesh
    type ElementMap_t
        integer      , allocatable :: eID(:,:,:)        ! Input element containing each output node
        real(kind=RP), allocatable :: xi(:,:,:,:)       ! Reference coordinates of each output node in that element
    end type ElementMap_t

    ! Background grid of bins to accelerate the point search
    type BinGrid_t
        real(kind=RP)              :: xmin(NDIM), h(NDIM)
        integer                    :: nb(NDIM)
        integer      , allocatable :: first(:)          ! CSR pointers to elems (size no_of_bins+1)
        integer      , allocatable :: elems(:)          ! Elements overlapping each bin
        real(kind=RP), allocatable :: bbmin(:,:)        ! Inflated bounding box of each element
        real(kind=RP), allocatable :: bbmax(:,:)
    end type BinGrid_t

    enum, bind(C)
        enumerator :: POINT_INSIDE = 0, POINT_TOLERANCE, POINT_PROJECTED
    end enum

    character(len=*), parameter :: TOOL_TMP_FILE = ".horses_mesh_interpolation_files.tmp"

    integer       :: NmaxStorage      = -1
    real(kind=RP) :: outsideTolerance = 0.1_RP

    contains

    subroutine MeshInterpolation(controlVariables)
        use FTValueDictionaryClass
        use HexMeshClass
        implicit none
        TYPE(FTValueDictionary)                 :: controlVariables

        ! Local variables
        type(HexMesh)                           :: M_in, M_out
        type(ElementMap_t), allocatable         :: map(:)
        type(BinGrid_t)                         :: grid
        type(FileLayout_t)                      :: layout
        type(ElementData_t), allocatable        :: data_in(:), data_out(:)
        character(len=LINE_LENGTH), allocatable :: fileNames(:)
        character(len=LINE_LENGTH)              :: field, outDir, suffix
        character(len=LINE_LENGTH)              :: inputMeshName, currentMeshName, outName
        logical                                 :: meshPerFile
        integer                                 :: iFile

        !
        ! Read the options and the files to interpolate
        !
        call readOptions(controlVariables, field, outDir, suffix, meshPerFile)
        call getInputFileNames(controlVariables, field, fileNames)

        write(STD_OUT,'(/,A)') "Mesh interpolation"
        write(STD_OUT,'(A)')   "------------------"
        write(STD_OUT,'(A,A)')  "   Interpolation field:     ", trim(field)
        write(STD_OUT,'(A,I0)') "   Number of input files:   ", size(fileNames)
        write(STD_OUT,'(A,A)')  "   Output directory:        ", trim(outDir)
        write(STD_OUT,'(A,L1)') "   Input hmesh per file:    ", meshPerFile

        !
        ! Construct the output mesh
        !
        call buildOutputMesh(controlVariables, M_out)
        call execute_command_line("mkdir -p " // trim(outDir))

        !
        ! Interpolate all the files
        !
        currentMeshName = ""
        do iFile = 1, size(fileNames)
            write(STD_OUT,'(/,A,I0,A,I0,A,A)') "File ", iFile, " of ", size(fileNames), ": ", trim(fileNames(iFile))

            ! Construct the input mesh and the mapping only if the input mesh changes
            inputMeshName = getInputMeshName(controlVariables, fileNames(iFile), meshPerFile)
            if ( trim(inputMeshName) .ne. trim(currentMeshName) ) then
                call buildInputMesh(inputMeshName, M_in)
                call buildBinGrid(M_in, grid)
                call computeMapping(M_in, M_out, grid, map)
                currentMeshName = inputMeshName
            end if

            call readFileLayout(controlVariables, fileNames(iFile), M_in, field, layout)
            call readFileData(fileNames(iFile), M_in, layout, data_in)
            call interpolateData(M_in, M_out, map, layout, data_in, data_out)

            outName = getOutputFileName(fileNames(iFile), outDir, suffix)
            call writeFileData(outName, M_out, layout, data_out)
            write(STD_OUT,'(A,A)') "   Written: ", trim(outName)

            deallocate(data_in, data_out)
        end do

        write(STD_OUT,'(/,A)') "Mesh interpolation finished"

    end subroutine MeshInterpolation
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!     Options and file names
!     ----------------------
!
!////////////////////////////////////////////////////////////////////////////////////////
!
    subroutine readOptions(controlVariables, field, outDir, suffix, meshPerFile)
        use FTValueDictionaryClass
        use Utilities, only: toLower
        implicit none
        TYPE(FTValueDictionary)                 :: controlVariables
        character(len=LINE_LENGTH), intent(out) :: field, outDir, suffix
        logical                   , intent(out) :: meshPerFile

        field = "auto"
        if ( controlVariables % containsKey("interpolation field") ) then
            field = controlVariables % stringValueForKey("interpolation field", requestedLength = LINE_LENGTH)
            call toLower(field)
        end if
        select case ( trim(field) )
        case ("auto", "solution", "stats", "lamb vector")
        case default
            print *, "Unknown interpolation field: ", trim(field)
            print *, "Options available are: auto, solution, stats, lamb vector"
            errorMessage(STD_OUT)
            error stop
        end select

        outDir = "RESULTS"
        if ( controlVariables % containsKey("output directory") ) then
            outDir = controlVariables % stringValueForKey("output directory", requestedLength = LINE_LENGTH)
        end if

        suffix = ".Interpolated"
        if ( controlVariables % containsKey("output file suffix") ) then
            suffix = controlVariables % stringValueForKey("output file suffix", requestedLength = LINE_LENGTH)
        end if

        meshPerFile = getLogicalValue(controlVariables, "input hmesh per file", .false.)

        if ( controlVariables % containsKey("outside point tolerance") ) then
            outsideTolerance = controlVariables % doublePrecisionValueForKey("outside point tolerance")
        end if

    end subroutine readOptions

    logical function getLogicalValue(controlVariables, keyName, default)
        use FTValueDictionaryClass
        use Utilities, only: toLower
        implicit none
        TYPE(FTValueDictionary)       :: controlVariables
        character(len=*), intent(in)  :: keyName
        logical         , intent(in)  :: default

        character(len=LINE_LENGTH) :: val

        getLogicalValue = default
        if ( .not. controlVariables % containsKey(keyName) ) return

        val = controlVariables % stringValueForKey(keyName, requestedLength = LINE_LENGTH)
        call toLower(val)
        select case ( trim(adjustl(val)) )
        case (".true.", "true", "t", "yes", ".t.")
            getLogicalValue = .true.
        case default
            getLogicalValue = .false.
        end select

    end function getLogicalValue

    subroutine getInputFileNames(controlVariables, field, fileNames)
        use FTValueDictionaryClass
        implicit none
        TYPE(FTValueDictionary)                                :: controlVariables
        character(len=*)                         , intent(in)  :: field
        character(len=LINE_LENGTH), allocatable  , intent(out) :: fileNames(:)

        ! Local variables
        character(len=LINE_LENGTH) :: keyName, val, dir
        integer :: i, n
        character(len=LINE_LENGTH), parameter :: statsKeys(4) = [character(len=LINE_LENGTH) :: &
                                                    "stats qbase file name", &
                                                    "stats sound velocity squared file name", &
                                                    "stats gradient sound velocity squared file name", &
                                                    "stats lamb vector file name"]

        allocate(fileNames(0))

        ! Explicit file names
        if ( controlVariables % containsKey("input file name") ) then
            call appendName(fileNames, controlVariables % stringValueForKey("input file name", requestedLength = LINE_LENGTH))
        end if
        i = 2
        do
            write(keyName,'(A,I0)') "input file name ", i
            if ( .not. controlVariables % containsKey(trim(keyName)) ) exit
            call appendName(fileNames, controlVariables % stringValueForKey(trim(keyName), requestedLength = LINE_LENGTH))
            i = i + 1
        end do

        ! File pattern
        if ( controlVariables % containsKey("input file pattern") ) then
            val = controlVariables % stringValueForKey("input file pattern", requestedLength = LINE_LENGTH)
            call globFiles(val, fileNames)
        end if

        ! Compatibility with the keywords of StatsMeshInterpolation
        if ( (size(fileNames) .eq. 0) .and. ((trim(field) .eq. "stats") .or. (trim(field) .eq. "auto")) ) then
            do i = 1, size(statsKeys)
                if ( controlVariables % containsKey(trim(statsKeys(i))) ) then
                    call appendName(fileNames, controlVariables % stringValueForKey(trim(statsKeys(i)), requestedLength = LINE_LENGTH))
                end if
            end do
        end if

        ! Compatibility with the keywords of LambVectorMeshInterpolation
        if ( (size(fileNames) .eq. 0) .and. controlVariables % containsKey("lamb vector file name") ) then
            dir = "."
            if ( controlVariables % containsKey("lamb vector dir") ) then
                dir = controlVariables % stringValueForKey("lamb vector dir", requestedLength = LINE_LENGTH)
            end if
            val = trim(dir) // "/" // trim(controlVariables % stringValueForKey("lamb vector file name", requestedLength = LINE_LENGTH)) // "*"
            call globFiles(val, fileNames)
        end if

        n = size(fileNames)
        if ( n .eq. 0 ) then
            print *, "No files to interpolate were found. Use:"
            print *, "   input file name    = path/to/file.hsol"
            print *, "   input file name 2  = path/to/file2.hsol   (optional)"
            print *, "   input file pattern = path/to/files_*.hsol (optional)"
            errorMessage(STD_OUT)
            error stop
        end if

    end subroutine getInputFileNames

    subroutine appendName(names, name)
        implicit none
        character(len=LINE_LENGTH), allocatable, intent(inout) :: names(:)
        character(len=*)                       , intent(in)    :: name

        character(len=LINE_LENGTH), allocatable :: aux(:)
        integer :: n

        n = size(names)
        allocate(aux(n+1))
        aux(1:n) = names
        aux(n+1) = adjustl(name)
        call move_alloc(aux, names)

    end subroutine appendName

    subroutine globFiles(pattern, names)
        implicit none
        character(len=*)                       , intent(in)    :: pattern
        character(len=LINE_LENGTH), allocatable, intent(inout) :: names(:)

        integer                    :: fid, io, nFound
        character(len=LINE_LENGTH) :: line

        call execute_command_line("ls -1d " // trim(pattern) // " > " // TOOL_TMP_FILE // " 2> /dev/null")

        nFound = 0
        open(newunit=fid, file=TOOL_TMP_FILE, status="old", action="read", iostat=io)
        if ( io .eq. 0 ) then
            do
                read(fid,'(A)',iostat=io) line
                if ( io .ne. 0 ) exit
                if ( len_trim(line) .eq. 0 ) cycle
                call appendName(names, line)
                nFound = nFound + 1
            end do
            close(fid, status="delete")
        end if

        if ( nFound .eq. 0 ) then
            print *, "WARNING: no files match the pattern ", trim(pattern)
        end if

    end subroutine globFiles

    ! Name of the file without path and extension (.hsol, .stats.hsol or .hmesh)
    subroutine splitFileName(fileName, base, ext)
        use FileReadingUtilities, only: RemovePath
        implicit none
        character(len=*)          , intent(in)  :: fileName
        character(len=LINE_LENGTH), intent(out) :: base
        character(len=LINE_LENGTH), intent(out) :: ext

        integer :: n

        base = RemovePath(trim(fileName))
        ext  = ""
        n = len_trim(base)
        if ( n .gt. 6 ) then
            if ( base(n-5:n) .eq. ".hmesh" ) then
                base = base(1:n-6)
                ext  = ".hmesh"
                return
            end if
        end if
        if ( n .gt. 5 ) then
            if ( base(n-4:n) .eq. ".hsol" ) then
                base = base(1:n-5)
                ext  = ".hsol"
            end if
        end if

        n = len_trim(base)
        if ( n .gt. 6 ) then
            if ( base(n-5:n) .eq. ".stats" ) then
                base = base(1:n-6)
                ext  = ".stats" // trim(ext)
            end if
        end if

    end subroutine splitFileName

    character(len=LINE_LENGTH) function getOutputFileName(fileName, outDir, suffix)
        implicit none
        character(len=*), intent(in) :: fileName, outDir, suffix

        character(len=LINE_LENGTH) :: base, ext

        call splitFileName(fileName, base, ext)
        if ( len_trim(ext) .eq. 0 ) ext = ".hsol"
        getOutputFileName = trim(outDir) // "/" // trim(base) // trim(suffix) // trim(ext)

    end function getOutputFileName

    character(len=LINE_LENGTH) function getInputMeshName(controlVariables, fileName, meshPerFile)
        use FTValueDictionaryClass
        use FileReadingUtilities, only: RemovePath, getFileName
        implicit none
        TYPE(FTValueDictionary)       :: controlVariables
        character(len=*), intent(in)  :: fileName
        logical         , intent(in)  :: meshPerFile

        character(len=LINE_LENGTH) :: meshName, base, ext, dir, tag
        integer                    :: pos

        if ( controlVariables % containsKey("input hmesh file name") ) then
            meshName = controlVariables % stringValueForKey("input hmesh file name", requestedLength = LINE_LENGTH)
        elseif ( controlVariables % containsKey("stats mesh file name") ) then
            meshName = controlVariables % stringValueForKey("stats mesh file name", requestedLength = LINE_LENGTH)
        else
            print *, "The mesh of the input files is not specified. Use:"
            print *, "   input hmesh file name = path/to/mesh.hmesh"
            errorMessage(STD_OUT)
            error stop
        end if

        if ( .not. meshPerFile ) then
            getInputMeshName = meshName
            return
        end if

        ! Iteration tag of the file: text after the last "_" until the next "." (e.g. case_0000000100.Lamb.hsol)
        call splitFileName(fileName, base, ext)
        pos = index(base, '_', BACK=.true.)
        tag = base(pos+1:)
        pos = index(tag, '.')
        if ( pos .gt. 0 ) tag = tag(1:pos-1)

        ! <path of the input hmesh>/<name of the input hmesh>_<tag>.hmesh
        pos = index(meshName, '/', BACK=.true.)
        dir = ""
        if ( pos .gt. 0 ) dir = meshName(1:pos)
        getInputMeshName = trim(dir) // trim(getFileName(RemovePath(trim(meshName)))) // "_" // trim(tag) // ".hmesh"

    end function getInputMeshName
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!     Meshes
!     ------
!
!////////////////////////////////////////////////////////////////////////////////////////
!
    subroutine ensureNodalStorage(Nmax)
        ! Construct the nodal storage of both Gauss and Gauss-Lobatto nodes up to order Nmax
        implicit none
        integer, intent(in) :: Nmax

        integer :: k

        if ( Nmax .le. NmaxStorage ) return

        call InitializeNodalStorage(GAUSS, Nmax)
        call InitializeNodalStorage(GAUSSLOBATTO, Nmax)
        do k = 0, Nmax
            call NodalStorage_Gauss(k) % Construct(GAUSS, k)
            call NodalStorage_GaussLobatto(k) % Construct(GAUSSLOBATTO, k)
        end do
        NmaxStorage = Nmax

    end subroutine ensureNodalStorage

    function nodalStorageFor(nodeType) result(sp)
        implicit none
        integer, intent(in)           :: nodeType
        type(NodalStorage_t), pointer :: sp(:)

        select case (nodeType)
        case (GAUSS)
            sp => NodalStorage_Gauss
        case (GAUSSLOBATTO)
            sp => NodalStorage_GaussLobatto
        case default
            print *, "Unknown node type: ", nodeType
            errorMessage(STD_OUT)
            error stop
        end select

    end function nodalStorageFor

    character(len=LINE_LENGTH) function nodeTypeName(nodeType)
        implicit none
        integer, intent(in) :: nodeType

        select case (nodeType)
        case (GAUSS)
            nodeTypeName = "Gauss"
        case (GAUSSLOBATTO)
            nodeTypeName = "Gauss-Lobatto"
        case default
            nodeTypeName = "Unknown"
        end select

    end function nodeTypeName

    subroutine readHmesh(meshName, mesh)
        ! Read the nodes of the mesh from a .hmesh file and prepare it for IO
        use HexMeshClass
        use StatsMeshInterpolation, only: readMesh_FromHmeshFile
        implicit none
        character(len=*), intent(in) :: meshName
        type(HexMesh)                :: mesh

        integer :: nodeType, Nmax, eID
        logical :: exists

        inquire(file=trim(meshName), exist=exists)
        if ( .not. exists ) then
            print *, "Mesh file ", trim(meshName), " does not exist."
            errorMessage(STD_OUT)
            error stop
        end if

        safedeallocate(mesh % elements)
        call readMesh_FromHmeshFile(trim(meshName), mesh, nodeType, Nmax)
        mesh % nodeType          = nodeType
        mesh % no_of_allElements = mesh % no_of_elements
        do eID = 1, mesh % no_of_elements
            mesh % elements(eID) % eID    = eID
            mesh % elements(eID) % globID = eID
        end do
        call mesh % PrepareForIO()

        Nmax = 0
        do eID = 1, mesh % no_of_elements
            Nmax = max(Nmax, maxval(mesh % elements(eID) % Nxyz))
        end do
        call ensureNodalStorage(Nmax)

    end subroutine readHmesh

    subroutine buildInputMesh(meshName, M_in)
        use HexMeshClass
        implicit none
        character(len=*), intent(in) :: meshName
        type(HexMesh)                :: M_in

        write(STD_OUT,'(A,A)') "   Reading input mesh: ", trim(meshName)
        call readHmesh(meshName, M_in)
        write(STD_OUT,'(A,I0,A,A)') "      ", M_in % no_of_elements, " elements, nodes: ", trim(nodeTypeName(M_in % nodeType))

    end subroutine buildInputMesh

    subroutine buildOutputMesh(controlVariables, M_out)
        use FTValueDictionaryClass
        use HexMeshClass
        use mainKeywordsModule    , only: discretizationNodesKey
        use pAdaptationClass      , only: GetMeshPolynomialOrders
        use StatsMeshInterpolation, only: getNodeType
        use ReadMeshFile          , only: constructMeshFromFile
        use PhysicsStorage        , only: ConstructPhysicsStorage
        use InterpolationMatrices , only: Initialize_InterpolationMatrices, Finalize_InterpolationMatrices
        implicit none
        TYPE(FTValueDictionary)  :: controlVariables
        type(HexMesh)            :: M_out

        character(len=LINE_LENGTH) :: fileName
        integer, allocatable       :: Nx(:), Ny(:), Nz(:)
        integer                    :: Nmax, nodeType, eID
        logical                    :: success

        if ( controlVariables % containsKey("output hmesh file name") ) then
            !
            ! Nodes and discretization nodes read from a .hmesh file
            !
            fileName = controlVariables % stringValueForKey("output hmesh file name", requestedLength = LINE_LENGTH)
            write(STD_OUT,'(A,A)') "   Reading output mesh: ", trim(fileName)
            call readHmesh(fileName, M_out)
            if ( controlVariables % containsKey(discretizationNodesKey) ) then
                call getNodeType(controlVariables, nodeType)
                if ( nodeType .ne. M_out % nodeType ) then
                    print *, "WARNING: the discretization nodes of the control file differ from those of the output hmesh."
                    print *, "         The nodes of the output hmesh are used."
                end if
            end if

        else
            !
            ! Mesh file + polynomial order + discretization nodes
            !
            if ( .not. controlVariables % containsKey("mesh file name") ) then
                print *, "The output mesh is not specified. Use either:"
                print *, "   output hmesh file name = path/to/mesh.hmesh"
                print *, "or"
                print *, "   mesh file name         = path/to/mesh.(msh/h5/mesh)"
                print *, "   polynomial order       = N"
                print *, "   discretization nodes   = Gauss/Gauss-Lobatto"
                errorMessage(STD_OUT)
                error stop
            end if
            fileName = controlVariables % stringValueForKey("mesh file name", requestedLength = LINE_LENGTH)
            write(STD_OUT,'(A,A)') "   Reading output mesh: ", trim(fileName)
            call GetMeshPolynomialOrders(controlVariables, Nx, Ny, Nz, Nmax)
            call getNodeType(controlVariables, nodeType)

            ! The mesh readers of the solver use the reference length of the physics. The other
            ! physical parameters are not used, so default values are given if they are missing
            call addDefaultValue(controlVariables, "flow equations" , "NS")
            call addDefaultValue(controlVariables, "mach number"    , "0.1")
            call addDefaultValue(controlVariables, "reynolds number", "1.0")
            call ConstructPhysicsStorage(controlVariables, success)

            ! Nodal storage and interpolation matrices of the output nodes, used by the mesh readers
            Nmax = max(maxval(Nx), maxval(Ny), maxval(Nz))
            call ensureNodalStorage(Nmax)
            NodalStorage => nodalStorageFor(nodeType)
            CurrentNodes =  nodeType
            call Initialize_InterpolationMatrices(Nmax)

            ! Full mesh construction (curved elements supported in all the mesh formats)
            call constructMeshFromFile(M_out, trim(fileName), nodeType, Nx, Ny, Nz, .true., 0, &
                                       getLogicalValue(controlVariables, "periodic relative tolerance", .false.), success)
            if ( .not. success ) then
                print *, "The output mesh could not be constructed."
                errorMessage(STD_OUT)
                error stop
            end if
            call Finalize_InterpolationMatrices

            M_out % nodeType = nodeType
            if ( M_out % no_of_allElements .eq. 0 ) M_out % no_of_allElements = M_out % no_of_elements
            do eID = 1, M_out % no_of_elements
                if ( M_out % elements(eID) % globID .le. 0 ) M_out % elements(eID) % globID = eID
            end do
            call M_out % PrepareForIO()
        end if

        write(STD_OUT,'(A,I0,A,A)') "      ", M_out % no_of_elements, " elements, nodes: ", trim(nodeTypeName(M_out % nodeType))

    end subroutine buildOutputMesh

    subroutine addDefaultValue(controlVariables, keyName, val)
        use FTValueDictionaryClass
        implicit none
        TYPE(FTValueDictionary)      :: controlVariables
        character(len=*), intent(in) :: keyName, val

        if ( .not. controlVariables % containsKey(keyName) ) call controlVariables % addValueForKey(val, keyName)

    end subroutine addDefaultValue
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!     Point search and mapping
!     ------------------------
!
!////////////////////////////////////////////////////////////////////////////////////////
!
    subroutine buildBinGrid(M_in, grid)
        use HexMeshClass
        implicit none
        type(HexMesh)   , intent(in)    :: M_in
        type(BinGrid_t) , intent(inout) :: grid

        integer       :: eID, d, nElem, nBins, b, i, j, k, lo(NDIM), hi(NDIM)
        integer, allocatable :: counter(:)
        real(kind=RP) :: width(NDIM), xmax(NDIM), L(NDIM), h

        nElem = M_in % no_of_elements
        safedeallocate(grid % bbmin) ; allocate(grid % bbmin(NDIM, nElem))
        safedeallocate(grid % bbmax) ; allocate(grid % bbmax(NDIM, nElem))

        ! Bounding boxes of the elements, inflated to contain the whole (curved) element
        ! even if the nodes are inside the element (Gauss nodes)
        do eID = 1, nElem
            associate( x => M_in % elements(eID) % geom % x )
            do d = 1, NDIM
                grid % bbmin(d, eID) = minval(x(d,:,:,:))
                grid % bbmax(d, eID) = maxval(x(d,:,:,:))
            end do
            end associate
            width = grid % bbmax(:,eID) - grid % bbmin(:,eID)
            grid % bbmin(:,eID) = grid % bbmin(:,eID) - 0.5_RP * maxval(width)
            grid % bbmax(:,eID) = grid % bbmax(:,eID) + 0.5_RP * maxval(width)
        end do

        ! Bins of the size of the average element
        grid % xmin = minval(grid % bbmin, dim=2)
        xmax        = maxval(grid % bbmax, dim=2)
        L = max(xmax - grid % xmin, 1.0e-12_RP * max(maxval(xmax - grid % xmin), 1.0_RP))
        h = (product(L) / real(nElem, RP)) ** (1.0_RP / 3.0_RP)
        do d = 1, NDIM
            grid % nb(d) = max(1, min(nint(L(d) / h), 1024))
        end do
        grid % h = L / real(grid % nb, RP)
        nBins = product(grid % nb)

        ! Count the elements in each bin
        allocate(counter(nBins))
        counter = 0
        do eID = 1, nElem
            lo = binIndex(grid, grid % bbmin(:,eID))
            hi = binIndex(grid, grid % bbmax(:,eID))
            do k = lo(3), hi(3) ; do j = lo(2), hi(2) ; do i = lo(1), hi(1)
                b = binLinearIndex(grid, [i,j,k])
                counter(b) = counter(b) + 1
            end do              ; end do              ; end do
        end do

        ! Fill the CSR structure
        safedeallocate(grid % first) ; allocate(grid % first(nBins+1))
        grid % first(1) = 1
        do b = 1, nBins
            grid % first(b+1) = grid % first(b) + counter(b)
        end do
        safedeallocate(grid % elems) ; allocate(grid % elems(grid % first(nBins+1)-1))
        counter = 0
        do eID = 1, nElem
            lo = binIndex(grid, grid % bbmin(:,eID))
            hi = binIndex(grid, grid % bbmax(:,eID))
            do k = lo(3), hi(3) ; do j = lo(2), hi(2) ; do i = lo(1), hi(1)
                b = binLinearIndex(grid, [i,j,k])
                grid % elems(grid % first(b) + counter(b)) = eID
                counter(b) = counter(b) + 1
            end do              ; end do              ; end do
        end do
        deallocate(counter)

        write(STD_OUT,'(A,I0,A,I0,A,I0,A)') "      Search grid: ", grid % nb(1), " x ", grid % nb(2), " x ", grid % nb(3), " bins"

    end subroutine buildBinGrid

    pure function binIndex(grid, x) result(ijk)
        implicit none
        type(BinGrid_t), intent(in) :: grid
        real(kind=RP)  , intent(in) :: x(NDIM)
        integer                     :: ijk(NDIM)

        ijk = int((x - grid % xmin) / grid % h) + 1
        ijk = max(1, min(ijk, grid % nb))

    end function binIndex

    pure integer function binLinearIndex(grid, ijk)
        implicit none
        type(BinGrid_t), intent(in) :: grid
        integer        , intent(in) :: ijk(NDIM)

        binLinearIndex = ijk(1) + (ijk(2)-1) * grid % nb(1) + (ijk(3)-1) * grid % nb(1) * grid % nb(2)

    end function binLinearIndex

    pure logical function insideBox(grid, eID, x)
        implicit none
        type(BinGrid_t), intent(in) :: grid
        integer        , intent(in) :: eID
        real(kind=RP)  , intent(in) :: x(NDIM)

        insideBox = all(x .ge. grid % bbmin(:,eID)) .and. all(x .le. grid % bbmax(:,eID))

    end function insideBox

    subroutine computeMapping(M_in, M_out, grid, map)
        use HexMeshClass
        implicit none
        type(HexMesh)                  , intent(in)    :: M_in
        type(HexMesh)                  , intent(in)    :: M_out
        type(BinGrid_t)                , intent(in)    :: grid
        type(ElementMap_t), allocatable, intent(inout) :: map(:)

        integer       :: eID, i, j, k, N(NDIM), eGuess, eFound, status
        integer       :: nInside, nTolerance, nProjected
        real(kind=RP) :: xi(NDIM)

        ! The global NodalStorage is used to find the points inside the input elements
        NodalStorage => nodalStorageFor(M_in % nodeType)

        safedeallocate(map)
        allocate(map(M_out % no_of_elements))

        nInside = 0 ; nTolerance = 0 ; nProjected = 0

        !$omp parallel default(shared) private(eID, i, j, k, N, eGuess, eFound, xi, status)
        eGuess = 0
        !$omp do schedule(runtime) reduction(+:nInside, nTolerance, nProjected)
        do eID = 1, M_out % no_of_elements
            N = M_out % elements(eID) % Nxyz
            allocate(map(eID) % eID(0:N(1), 0:N(2), 0:N(3)))
            allocate(map(eID) % xi(NDIM, 0:N(1), 0:N(2), 0:N(3)))
            do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                call locatePoint(M_in, grid, M_out % elements(eID) % geom % x(:,i,j,k), eGuess, eFound, xi, status)
                map(eID) % eID(i,j,k)  = eFound
                map(eID) % xi(:,i,j,k) = xi
                eGuess = eFound
                select case (status)
                case (POINT_INSIDE)
                    nInside = nInside + 1
                case (POINT_TOLERANCE)
                    nTolerance = nTolerance + 1
                case default
                    nProjected = nProjected + 1
                end select
            end do         ; end do         ; end do
        end do
        !$omp end do
        !$omp end parallel

        write(STD_OUT,'(A,I0)') "      Output nodes inside the input mesh:          ", nInside
        write(STD_OUT,'(A,I0)') "      Output nodes found within the tolerance:     ", nTolerance
        write(STD_OUT,'(A,I0)') "      Output nodes projected on the closest element: ", nProjected
        if ( nProjected .gt. 0 ) then
            print *, "WARNING: some output nodes are outside the input mesh. Their values are taken from the closest"
            print *, "         point of the closest input element. Check that both meshes represent the same geometry."
        end if

    end subroutine computeMapping

    subroutine locatePoint(M_in, grid, x, eGuess, eFound, xi, status)
        use HexMeshClass
        implicit none
        type(HexMesh)  , intent(in)  :: M_in
        type(BinGrid_t), intent(in)  :: grid
        real(kind=RP)  , intent(in)  :: x(NDIM)
        integer        , intent(in)  :: eGuess
        integer        , intent(out) :: eFound
        real(kind=RP)  , intent(out) :: xi(NDIM)
        integer        , intent(out) :: status

        integer       :: c, bin, e1, bestE
        real(kind=RP) :: xiTry(NDIM), bestXi(NDIM), bestDist

        status = POINT_INSIDE
        bestE  = 0
        bestDist = huge(1.0_RP)

        ! 1. Last element found (consecutive nodes are usually in the same element)
        if ( eGuess .gt. 0 ) then
            if ( insideBox(grid, eGuess, x) ) then
                if ( M_in % elements(eGuess) % FindPointWithCoords(x, 0, xiTry) ) then
                    eFound = eGuess
                    xi     = xiTry
                    return
                end if
                call keepBest(eGuess, xiTry)
            end if
        end if

        ! 2. Elements of the bin that contains the point
        bin = binLinearIndex(grid, binIndex(grid, x))
        do c = grid % first(bin), grid % first(bin+1) - 1
            e1 = grid % elems(c)
            if ( e1 .eq. eGuess ) cycle
            if ( .not. insideBox(grid, e1, x) ) cycle
            if ( M_in % elements(e1) % FindPointWithCoords(x, 0, xiTry) ) then
                eFound = e1
                xi     = xiTry
                return
            end if
            call keepBest(e1, xiTry)
        end do

        ! 3. Point slightly outside an element
        if ( (bestE .gt. 0) .and. (bestDist .le. 1.0_RP + outsideTolerance) ) then
            eFound = bestE
            xi     = max(-1.0_RP, min(1.0_RP, bestXi))
            status = POINT_TOLERANCE
            return
        end if

        ! 4. Point outside the input mesh: closest point of the closest element
        call projectOnClosestElement(M_in, grid, x, eFound, xi)
        status = POINT_PROJECTED

        contains

            subroutine keepBest(e, xiE)
                integer      , intent(in) :: e
                real(kind=RP), intent(in) :: xiE(NDIM)

                if ( any(xiE .ne. xiE) ) return    ! NaN
                if ( maxval(abs(xiE)) .lt. bestDist ) then
                    bestDist = maxval(abs(xiE))
                    bestE    = e
                    bestXi   = xiE
                end if

            end subroutine keepBest

    end subroutine locatePoint

    subroutine projectOnClosestElement(M_in, grid, x, eFound, xi)
        use HexMeshClass
        implicit none
        type(HexMesh)  , intent(in)  :: M_in
        type(BinGrid_t), intent(in)  :: grid
        real(kind=RP)  , intent(in)  :: x(NDIM)
        integer        , intent(out) :: eFound
        real(kind=RP)  , intent(out) :: xi(NDIM)

        integer       :: ijk(NDIM), r, rFound, i, j, k, c, bin, e1, ii, jj, kk, nodeFound(NDIM), N(NDIM)
        real(kind=RP) :: d2, d2min, xiTry(NDIM)
        logical       :: success

        ! Closest input node, searching in rings of bins around the point
        ijk    = binIndex(grid, x)
        d2min  = huge(1.0_RP)
        eFound = 0
        rFound = -1
        nodeFound = 0
        do r = 0, maxval(grid % nb)
            do k = max(1, ijk(3)-r), min(grid % nb(3), ijk(3)+r)
            do j = max(1, ijk(2)-r), min(grid % nb(2), ijk(2)+r)
            do i = max(1, ijk(1)-r), min(grid % nb(1), ijk(1)+r)
                if ( max(abs(i-ijk(1)), abs(j-ijk(2)), abs(k-ijk(3))) .ne. r ) cycle   ! Only the ring
                bin = binLinearIndex(grid, [i,j,k])
                do c = grid % first(bin), grid % first(bin+1) - 1
                    e1 = grid % elems(c)
                    N  = M_in % elements(e1) % Nxyz
                    do kk = 0, N(3) ; do jj = 0, N(2) ; do ii = 0, N(1)
                        d2 = sum((M_in % elements(e1) % geom % x(:,ii,jj,kk) - x)**2)
                        if ( d2 .lt. d2min ) then
                            d2min     = d2
                            eFound    = e1
                            nodeFound = [ii,jj,kk]
                        end if
                    end do          ; end do          ; end do
                end do
            end do
            end do
            end do
            if ( (eFound .gt. 0) .and. (rFound .lt. 0) ) rFound = r
            if ( (rFound .ge. 0) .and. (r .gt. rFound) ) exit   ! One more ring after the first hit
        end do

        ! Reference coordinates of the closest node
        N  = M_in % elements(eFound) % Nxyz
        xi = [NodalStorage(N(1)) % x(nodeFound(1)), NodalStorage(N(2)) % x(nodeFound(2)), NodalStorage(N(3)) % x(nodeFound(3))]

        ! Improve with the (clamped) solution of the mapping, if it is closer
        success = M_in % elements(eFound) % FindPointWithCoords(x, 0, xiTry)
        if ( all(xiTry .eq. xiTry) ) then
            xiTry = max(-1.0_RP, min(1.0_RP, xiTry))
            if ( sum((mapPoint(M_in % elements(eFound) % geom % x, N, xiTry) - x)**2) .lt. d2min ) xi = xiTry
        end if

    end subroutine projectOnClosestElement

    function mapPoint(xe, N, xi) result(x)
        ! Physical coordinates of the reference point xi of an element with nodes xe
        implicit none
        integer      , intent(in) :: N(NDIM)
        real(kind=RP), intent(in) :: xe(NDIM, 0:N(1), 0:N(2), 0:N(3))
        real(kind=RP), intent(in) :: xi(NDIM)
        real(kind=RP)             :: x(NDIM)

        real(kind=RP) :: lxi(0:N(1)), leta(0:N(2)), lzeta(0:N(3))
        integer       :: i, j, k

        lxi   = NodalStorage(N(1)) % lj(xi(1))
        leta  = NodalStorage(N(2)) % lj(xi(2))
        lzeta = NodalStorage(N(3)) % lj(xi(3))

        x = 0.0_RP
        do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            x = x + xe(:,i,j,k) * lxi(i) * leta(j) * lzeta(k)
        end do         ; end do         ; end do

    end function mapPoint
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!     File reading, interpolation and writing
!     ---------------------------------------
!
!////////////////////////////////////////////////////////////////////////////////////////
!
    subroutine readFileLayout(controlVariables, fileName, M_in, field, layout)
        use FTValueDictionaryClass
        use HexMeshClass
        use SolutionFile
        implicit none
        TYPE(FTValueDictionary)             :: controlVariables
        character(len=*)    , intent(in)    :: fileName
        type(HexMesh)       , intent(in)    :: M_in
        character(len=*)    , intent(in)    :: field
        type(FileLayout_t)  , intent(inout) :: layout

        integer                :: fid, rank, dims(4), eID
        integer(kind=AddrInt)  :: fileSize, headers, dataBytes, totalDOF, nodeBytes
        logical                :: exists, hasSensor

        inquire(file=trim(fileName), exist=exists)
        if ( .not. exists ) then
            print *, "File ", trim(fileName), " does not exist."
            errorMessage(STD_OUT)
            error stop
        end if

        ! Header
        layout % fileType       = getSolutionFileType(trim(fileName))
        layout % nodeType       = getSolutionFileNodeType(trim(fileName))
        layout % no_of_elements = getSolutionFileNoOfElements(trim(fileName))
        layout % refs           = getSolutionFileReferenceValues(trim(fileName))
        call getSolutionFileTimeAndIteration(trim(fileName), layout % iter, layout % time)

        if ( layout % no_of_elements .ne. M_in % no_of_elements ) then
            write(STD_OUT,'(A,I0,A,I0,A)') "The number of elements of the file (", layout % no_of_elements, &
                                           ") does not match that of the input mesh (", M_in % no_of_elements, ")."
            errorMessage(STD_OUT)
            error stop
        end if

        ! Dimensions of the first array of the first element
        fid = putSolutionFileInReadDataMode(trim(fileName))
        read(fid, pos=POS_INIT_DATA) rank
        if ( rank .ne. 4 ) then
            print *, "Unsupported array rank in file ", trim(fileName), ": ", rank
            errorMessage(STD_OUT)
            error stop
        end if
        read(fid) dims
        close(fid)

        ! Total number of variables per node, from the size of the file
        hasSensor = (layout % fileType .eq. SOLUTION_AND_SENSOR_FILE) .or. &
                    (layout % fileType .eq. SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE)
        inquire(file=trim(fileName), size=fileSize)
        totalDOF = 0
        do eID = 1, M_in % no_of_elements
            totalDOF = totalDOF + product(M_in % elements(eID) % Nxyz + 1)
        end do
        headers = int(POS_INIT_DATA - 1, AddrInt) + int(layout % no_of_elements, AddrInt) * 5_AddrInt * SIZEOF_INT &
                + int(SIZEOF_INT, AddrInt)   ! END_OF_FILE indicator
        if ( hasSensor ) headers = headers + int(layout % no_of_elements, AddrInt) * SIZEOF_RP
        dataBytes = fileSize - headers
        nodeBytes = totalDOF * SIZEOF_RP
        if ( mod(dataBytes, nodeBytes) .ne. 0 ) dataBytes = dataBytes + SIZEOF_INT     ! Not sealed file
        if ( (dataBytes .le. 0) .or. (mod(dataBytes, nodeBytes) .ne. 0) ) then
            print *, "The size of the file ", trim(fileName), " is not consistent with the input mesh."
            errorMessage(STD_OUT)
            error stop
        end if
        layout % padding = int(dataBytes / nodeBytes)

        call determineBlocks(controlVariables, layout, dims(1))
        call checkField(field, layout)

        write(STD_OUT,'(A,A,A,I0,A)') "   File type: ", trim(fileTypeName(layout % fileType)), &
                                      ", ", layout % padding, " variables per node"
        write(STD_OUT,'(A,*(I0,:,", "))') "   Arrays per element: ", layout % blocks
        if ( layout % nodeType .ne. M_in % nodeType ) then
            write(STD_OUT,'(A,A,A,A,A)') "   NOTE: the file uses ", trim(nodeTypeName(layout % nodeType)), &
                                         " nodes and the input mesh ", trim(nodeTypeName(M_in % nodeType)), " nodes."
        end if

    end subroutine readFileLayout

    subroutine determineBlocks(controlVariables, layout, nv1)
        ! Number of variables of each array stored per element
        use FTValueDictionaryClass
        use SolutionFile
        use FileReadingUtilities, only: getIntArrayFromString
        use Utilities           , only: toLower
        implicit none
        TYPE(FTValueDictionary)            :: controlVariables
        type(FileLayout_t), intent(inout)  :: layout
        integer           , intent(in)     :: nv1

        integer                    :: r, c, ngrad
        logical                    :: gradients, hasSensor
        character(len=LINE_LENGTH) :: statsSolver

        gradients = (layout % fileType .eq. SOLUTION_AND_GRADIENTS_FILE) .or. &
                    (layout % fileType .eq. SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE) .or. &
                    (layout % fileType .eq. STATS_AND_GRADIENTS_FILE)
        hasSensor = (layout % fileType .eq. SOLUTION_AND_SENSOR_FILE) .or. &
                    (layout % fileType .eq. SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE)

        layout % sensorBlock = 0
        r = layout % padding - nv1

        !
        ! Layout given by the user
        !
        if ( controlVariables % containsKey("interpolation variable blocks") ) then
            layout % blocks = getIntArrayFromString(controlVariables % stringValueForKey("interpolation variable blocks", &
                                                                                         requestedLength = LINE_LENGTH))
            if ( hasSensor ) then
                layout % sensorBlock = 1
                if ( gradients ) layout % sensorBlock = 4
            end if
            if ( (sum(layout % blocks) .ne. layout % padding) .or. (layout % blocks(1) .ne. nv1) ) then
                write(STD_OUT,'(A,I0,A,I0)') "The interpolation variable blocks are not consistent with the file: padding ", &
                                             layout % padding, ", first array ", nv1
                errorMessage(STD_OUT)
                error stop
            end if
            return
        end if

        !
        ! Automatic detection
        !
        select case ( layout % fileType )
        case ( SOLUTION_FILE, SOLUTION_AND_GRADIENTS_FILE, SOLUTION_AND_SENSOR_FILE, SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE )
            ! Q (NCONS) [+ U_x, U_y, U_z (NGRAD)] [+ sensor] [+ mu_NS, mu_sgs] [+ source (NCONS)]
            layout % blocks = [nv1]
            if ( gradients ) then
                if ( isSolutionExtra(r - 3*nv1, nv1) ) then
                    ngrad = nv1
                elseif ( mod(r, 3) .eq. 0 ) then
                    ngrad = r / 3
                else
                    call layoutError()
                end if
                layout % blocks = [layout % blocks, ngrad, ngrad, ngrad]
                r = r - 3*ngrad
            end if
            if ( hasSensor ) layout % sensorBlock = size(layout % blocks)
            if ( .not. isSolutionExtra(r, nv1) ) call layoutError()
            if ( (r .eq. 2) .or. (r .eq. nv1 + 2) ) layout % blocks = [layout % blocks, 1, 1]   ! LES viscosities
            if ( (r .eq. nv1) .or. (r .eq. nv1 + 2) ) layout % blocks = [layout % blocks, nv1]  ! Source term

        case ( STATS_FILE, STATS_AND_GRADIENTS_FILE )
            ! S_ij (9) + Q (NCONS) [+ rho (MU)] [+ U_x, U_y, U_z (NGRAD)], or a single array
            if ( r .eq. 0 ) then
                layout % blocks = [nv1]
                return
            end if

            statsSolver = ""
            if ( controlVariables % containsKey("stats solver") ) then
                statsSolver = controlVariables % stringValueForKey("stats solver", requestedLength = LINE_LENGTH)
                call toLower(statsSolver)
                if ( trim(statsSolver) .eq. "auto" ) statsSolver = ""
            end if

            if ( gradients ) then
                if ( (trim(statsSolver) .eq. "mu") .or. ((len_trim(statsSolver) .eq. 0) .and. (mod(r, 4) .eq. 1)) ) then
                    if ( mod(r-1, 4) .ne. 0 ) call layoutError()
                    c = (r - 1) / 4
                    layout % blocks = [nv1, c, 1, c, c, c]
                else
                    if ( mod(r, 4) .ne. 0 ) call layoutError()
                    c = r / 4
                    layout % blocks = [nv1, c, c, c, c]
                end if
            else
                if ( trim(statsSolver) .eq. "mu" ) then
                    layout % blocks = [nv1, r-1, 1]
                elseif ( (len_trim(statsSolver) .eq. 0) .and. (r .eq. 6) ) then
                    print *, "The layout of the statistics file is ambiguous. Specify the solver that wrote it:"
                    print *, "   stats solver = ns/ins/mu"
                    errorMessage(STD_OUT)
                    error stop
                else
                    layout % blocks = [nv1, r]
                end if
            end if

        case ( MESH_FILE )
            ! Node coordinates (.hmesh): interpolating them gives the .hmesh of the output mesh
            if ( r .ne. 0 ) call layoutError()
            layout % blocks = [nv1]

        case default
            print *, "Unsupported file type: ", trim(fileTypeName(layout % fileType))
            errorMessage(STD_OUT)
            error stop
        end select

        contains

            logical function isSolutionExtra(nExtra, nv)
                ! Possible extra variables of a solution file: none, LES (2), source (NCONS), LES + source
                integer, intent(in) :: nExtra, nv
                isSolutionExtra = (nExtra .eq. 0) .or. (nExtra .eq. 2) .or. (nExtra .eq. nv) .or. (nExtra .eq. nv + 2)
            end function isSolutionExtra

            subroutine layoutError()
                write(STD_OUT,'(A,I0,A,I0,A)') "Unable to detect the layout of the file (", layout % padding, &
                                               " variables per node, first array with ", nv1, " variables)."
                print *, "Specify it with: interpolation variable blocks = [n1, n2, ...]"
                errorMessage(STD_OUT)
                error stop
            end subroutine layoutError

    end subroutine determineBlocks

    subroutine checkField(field, layout)
        use SolutionFile
        implicit none
        character(len=*)  , intent(in) :: field
        type(FileLayout_t), intent(in) :: layout

        logical :: ok

        select case ( trim(field) )
        case ("solution")
            ok = (layout % fileType .eq. SOLUTION_FILE) .or. (layout % fileType .eq. SOLUTION_AND_GRADIENTS_FILE) .or. &
                 (layout % fileType .eq. SOLUTION_AND_SENSOR_FILE) .or. (layout % fileType .eq. SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE)
        case ("stats")
            ok = (layout % fileType .eq. STATS_FILE) .or. (layout % fileType .eq. STATS_AND_GRADIENTS_FILE)
        case ("lamb vector")
            ok = (layout % padding .eq. NDIM)
        case default
            ok = .true.
        end select

        if ( .not. ok ) then
            write(STD_OUT,'(A,A,A,A,A,I0,A)') "The file (", trim(fileTypeName(layout % fileType)), &
                                              ") is not a ", trim(field), " file (", layout % padding, " variables per node)."
            errorMessage(STD_OUT)
            error stop
        end if

    end subroutine checkField

    character(len=LINE_LENGTH) function fileTypeName(fileType)
        use SolutionFile
        implicit none
        integer, intent(in) :: fileType

        select case (fileType)
        case (SOLUTION_FILE)                          ; fileTypeName = "solution"
        case (SOLUTION_AND_GRADIENTS_FILE)            ; fileTypeName = "solution and gradients"
        case (SOLUTION_AND_SENSOR_FILE)               ; fileTypeName = "solution and sensor"
        case (SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE) ; fileTypeName = "solution, gradients and sensor"
        case (STATS_FILE)                             ; fileTypeName = "statistics"
        case (STATS_AND_GRADIENTS_FILE)               ; fileTypeName = "statistics and gradients"
        case (MESH_FILE)                              ; fileTypeName = "mesh"
        case default                                  ; fileTypeName = "unknown"
        end select

    end function fileTypeName

    integer(kind=AddrInt) function elementPosition(e, layout)
        ! Position of the data of an element in a .hsol file
        use ElementClass
        use SolutionFile, only: POS_INIT_DATA
        implicit none
        type(Element)     , intent(in) :: e
        type(FileLayout_t), intent(in) :: layout

        elementPosition = POS_INIT_DATA + (e % globID - 1) * 5_AddrInt * SIZEOF_INT &
                        + int(layout % padding, AddrInt) * e % offsetIO * SIZEOF_RP
        if ( layout % sensorBlock .gt. 0 ) elementPosition = elementPosition + (e % globID - 1) * int(SIZEOF_RP, AddrInt)

    end function elementPosition

    subroutine readFileData(fileName, M_in, layout, data_in)
        use HexMeshClass
        use SolutionFile
        implicit none
        character(len=*)                 , intent(in)  :: fileName
        type(HexMesh)                    , intent(in)  :: M_in
        type(FileLayout_t)               , intent(in)  :: layout
        type(ElementData_t), allocatable , intent(out) :: data_in(:)

        integer :: fid, eID, b, lb, ub, rank, dims(4), N(NDIM)

        allocate(data_in(M_in % no_of_elements))

        fid = putSolutionFileInReadDataMode(trim(fileName))
        do eID = 1, M_in % no_of_elements
            N = M_in % elements(eID) % Nxyz

            read(fid, pos=elementPosition(M_in % elements(eID), layout)) rank
            read(fid) dims
            if ( (rank .ne. 4) .or. (dims(1) .ne. layout % blocks(1)) .or. any(dims(2:4) - 1 .ne. N) ) then
                write(STD_OUT,'(A,I0,A)') "The data of element ", eID, " is not consistent with the input mesh."
                errorMessage(STD_OUT)
                error stop
            end if

            allocate(data_in(eID) % Q(layout % padding, 0:N(1), 0:N(2), 0:N(3)))
            lb = 1
            do b = 1, size(layout % blocks)
                ub = lb + layout % blocks(b) - 1
                read(fid) data_in(eID) % Q(lb:ub,:,:,:)
                if ( b .eq. layout % sensorBlock ) read(fid) data_in(eID) % sensor
                lb = ub + 1
            end do
        end do
        close(fid)

    end subroutine readFileData

    subroutine interpolateData(M_in, M_out, map, layout, data_in, data_out)
        use HexMeshClass
        implicit none
        type(HexMesh)                    , intent(in)  :: M_in, M_out
        type(ElementMap_t)               , intent(in)  :: map(:)
        type(FileLayout_t)               , intent(in)  :: layout
        type(ElementData_t)              , intent(in)  :: data_in(:)
        type(ElementData_t), allocatable , intent(out) :: data_out(:)

        type(NodalStorage_t), pointer :: sp(:)
        real(kind=RP), allocatable    :: lxi(:), leta(:), lzeta(:)
        integer                       :: eID, e1, i, j, k, ii, jj, kk, N(NDIM), Nin(NDIM), NmaxIn
        real(kind=RP)                 :: w

        ! The data is interpolated with the nodes of the file (they may differ from those of the hmesh)
        sp => nodalStorageFor(layout % nodeType)

        NmaxIn = 0
        do eID = 1, M_in % no_of_elements
            NmaxIn = max(NmaxIn, maxval(M_in % elements(eID) % Nxyz))
        end do

        allocate(data_out(M_out % no_of_elements))

        !$omp parallel default(shared) private(eID, e1, i, j, k, ii, jj, kk, N, Nin, w, lxi, leta, lzeta)
        allocate(lxi(0:NmaxIn), leta(0:NmaxIn), lzeta(0:NmaxIn))
        !$omp do schedule(runtime)
        do eID = 1, M_out % no_of_elements
            N = M_out % elements(eID) % Nxyz
            allocate(data_out(eID) % Q(layout % padding, 0:N(1), 0:N(2), 0:N(3)))
            data_out(eID) % Q = 0.0_RP
            data_out(eID) % sensor = -huge(1.0_RP)

            do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                e1  = map(eID) % eID(i,j,k)
                Nin = M_in % elements(e1) % Nxyz
                lxi(0:Nin(1))   = sp(Nin(1)) % lj(map(eID) % xi(1,i,j,k))
                leta(0:Nin(2))  = sp(Nin(2)) % lj(map(eID) % xi(2,i,j,k))
                lzeta(0:Nin(3)) = sp(Nin(3)) % lj(map(eID) % xi(3,i,j,k))

                do kk = 0, Nin(3) ; do jj = 0, Nin(2) ; do ii = 0, Nin(1)
                    w = lxi(ii) * leta(jj) * lzeta(kk)
                    data_out(eID) % Q(:,i,j,k) = data_out(eID) % Q(:,i,j,k) + w * data_in(e1) % Q(:,ii,jj,kk)
                end do            ; end do            ; end do

                ! Element sensor: maximum of the input elements that contain the nodes
                data_out(eID) % sensor = max(data_out(eID) % sensor, data_in(e1) % sensor)
            end do         ; end do         ; end do
        end do
        !$omp end do
        deallocate(lxi, leta, lzeta)
        !$omp end parallel

    end subroutine interpolateData

    subroutine writeFileData(fileName, M_out, layout, data_out)
        use HexMeshClass
        use SolutionFile
        implicit none
        character(len=*)   , intent(in) :: fileName
        type(HexMesh)      , intent(in) :: M_out
        type(FileLayout_t) , intent(in) :: layout
        type(ElementData_t), intent(in) :: data_out(:)

        integer :: fid, eID, b, lb, ub

        call CreateNewSolutionFile(trim(fileName), layout % fileType, M_out % nodeType, M_out % no_of_allElements, &
                                   layout % iter, layout % time, layout % refs)

        fid = putSolutionFileInWriteDataMode(trim(fileName))
        do eID = 1, M_out % no_of_elements
            lb = 1
            do b = 1, size(layout % blocks)
                ub = lb + layout % blocks(b) - 1
                if ( b .eq. 1 ) then
                    call writeArray(fid, data_out(eID) % Q(lb:ub,:,:,:), position=elementPosition(M_out % elements(eID), layout))
                else
                    write(fid) data_out(eID) % Q(lb:ub,:,:,:)
                end if
                if ( b .eq. layout % sensorBlock ) write(fid) data_out(eID) % sensor
                lb = ub + 1
            end do
        end do
        close(fid)

        call SealSolutionFile(trim(fileName))

    end subroutine writeFileData

end module MeshInterpolationTool
