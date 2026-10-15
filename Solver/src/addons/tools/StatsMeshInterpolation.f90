#include "Includes.h"
module StatsMeshInterpolation
    use SMConstants
    use omp_lib

    implicit none

    private 
    public InterpolateStats
    public readMesh_FromHmeshFile, getNodeType ! Also used by MeshInterpolationTool

    type StatsElementStorage_t
        real(rp), allocatable :: Q(:,:,:,:)
        real(rp), allocatable :: a2(:,:,:,:) ! soundVelocitySquared
        real(rp), allocatable :: grada2(:,:,:,:) ! gradSoundVelocitySquared
        real(rp), allocatable :: LambStats(:,:,:,:) ! Lamb vector stats
        ! real(rp), allocatable :: Lamb(:,:,:,:) ! Lamb vector
    end type StatsElementStorage_t

    enum, bind(C)
        enumerator :: statsSolverNS = 1, statsSolveriNS, statsSolverMU
    end enum
    integer :: statsSolverID

    logical :: saveGradients = .false.

    contains

    subroutine InterpolateStats(controlVariables)
        ! Interpolate fields defined in mesh M_in onto the mesh M_out.
        use FTValueDictionaryClass
        use pAdaptationClass          , only: GetMeshPolynomialOrders
        use HexMeshClass
        use StorageClass
        use NodalStorageClass
        use ProbeClass
        use TransfiniteMapClass
        use SolutionFile, only: NO_OF_SAVED_REFS
        use FileReadingUtilities
        implicit None
        TYPE(FTValueDictionary)                :: controlVariables
        !
        ! Local variables
        !
        type(HexMesh)                           :: M_in_aux, M_in, M_out
        CHARACTER(LEN=LINE_LENGTH)              :: fileName
        integer :: nodeType_in, nodeType_out
        integer, allocatable                :: Nx(:), Ny(:), Nz(:)
        integer                             :: Nmax_in, Nmax_out, Nmax
        integer                             :: eID, i, j, k, ii, jj, kk

        type(NodalStorage_t), target, allocatable :: NodalStorage_in(:)
        type(NodalStorage_t), pointer :: NodalStorage_out(:)

        type(StatsElementStorage_t), allocatable :: storage_e_in(:), storage_e_out(:) ! To save the stats
        integer :: no_of_stats
        logical :: allocateSoundVelocity
        type(TransfiniteHexMap), pointer :: hexMap, hex8Map, genHexMap
        type(Probe_t) :: probe
        real(kind=RP) :: refs(NO_OF_SAVED_REFS)
        integer :: iter
        real(rp) :: time

        !
        ! Read mesh M_in
        !
        ! First, read the hmesh file to get the nodes for each element
        fileName = controlVariables % stringValueForKey("stats mesh file name", requestedLength = LINE_LENGTH)
        print *, "Reading mesh M_in from hmesh: ", trim(fileName)
        call readMesh_FromHmeshFile(fileName, M_in_aux, nodeType_in, Nmax_in)
        print *, "DONE"
        ! Second, put the polynomial degree of each element in Nx, Ny, Nz array
        allocate(Nx(M_in_aux % no_of_elements))
        allocate(Ny(M_in_aux % no_of_elements))
        allocate(Nz(M_in_aux % no_of_elements))
        do eID = 1, M_in_aux % no_of_elements
            Nx(eID) = M_in_aux % elements(eID) % Nxyz(1)
            Ny(eID) = M_in_aux % elements(eID) % Nxyz(2)
            Nz(eID) = M_in_aux % elements(eID) % Nxyz(3)
        end do
        ! Second, read the msh file to construct all the geometry related data
        fileName = controlVariables % stringValueForKey("stats geo mesh file name", requestedLength = LINE_LENGTH)
        print *, "Reading mesh M_in from geo: ", trim(fileName)
        call ConstructSimpleMesh(M_in, fileName, AllNx=Nx, AllNy=Ny, AllNz=Nz)
        print *, "DONE"
        ! Assign data from hmesh to M_in
        do eID = 1, M_in % no_of_elements
            M_in % elements(eID) % Nxyz = M_in_aux % elements(eID) % Nxyz
            M_in % elements(eID) % geom % x = M_in_aux % elements(eID) % geom % x
        end do
        M_in % nodeType = nodeType_in
        call M_in % PrepareForIO()
        deallocate(Nx)
        deallocate(Ny)
        deallocate(Nz)
        
        !
        ! Read mesh M_out
        !
        print *, "Reading polynomial order for M_out"
        call GetMeshPolynomialOrders(controlVariables,Nx,Ny,Nz,Nmax_out)
        print *, "DONE"
        fileName = controlVariables % stringValueForKey("mesh file name", requestedLength = LINE_LENGTH)
        print *, "Reading mesh M_out: ", trim(fileName)
        call ConstructSimpleMesh(M_out, fileName, AllNx=Nx, AllNy=Ny, AllNz=Nz)
        print *, "DONE"
        deallocate(Nx)
        deallocate(Ny)
        deallocate(Nz)
        ! Get the node type of the mesh M_out
        call getNodeType(controlVariables, nodetype_out)
        M_out % nodeType = nodeType_out
        
        ! Initialize the nodal storage
        if (M_in % nodeType .ne. M_out % nodeType) then
            allocate(NodalStorage_in(0:Nmax_in))
            allocate(NodalStorage_out(0:Nmax_out))
        else
            Nmax = max(Nmax_in, Nmax_out)
            allocate(NodalStorage_in(0:Nmax))
            NodalStorage_out => NodalStorage_in
        end if
        ! Construct nodal storage data for the input mesh
        do eID = 1, M_in % no_of_elements
            do k = 0, M_in % elements(eID) % Nxyz(3)
                call NodalStorage_in(k) % Construct( M_in % nodeType, k )
            end do
            do j = 0, M_in % elements(eID) % Nxyz(1)
                call NodalStorage_in(j) % Construct( M_in % nodeType, j )
            end do
            do i = 0, M_in % elements(eID) % Nxyz(1)
                call NodalStorage_in(i) % Construct( M_in % nodeType, i )
            end do
        end do
        ! Construct nodal storage data for the output mesh
        do eID = 1, M_out % no_of_elements
            do k = 0, M_out % Nz(eID)
                call NodalStorage_out(k) % Construct( M_out % nodeType, k )
            end do
            do j = 0, M_out % Ny(eID)
                call NodalStorage_out(j) % Construct( M_out % nodeType, j )
            end do
            do i = 0, M_out % Nx(eID)
                call NodalStorage_out(i) % Construct( M_out % nodeType, i )
            end do
        end do


        !
        ! Read stats in M_in
        !
        ! Initialize storage for mesh M_in
        allocate(storage_e_in(M_in % no_of_elements))
        fileName = controlVariables % stringValueForKey("stats qbase file name", requestedLength = LINE_LENGTH)
        call readStats(controlVariables, M_in, storage_e_in, refs, iter, time, fileName)


        ! Initialize storage for mesh M_out
        allocate(storage_e_out(M_out % no_of_elements))
        no_of_stats = size(storage_e_in(1) % Q, dim=1)
        allocateSoundVelocity = allocated(storage_e_in(1) % a2)
        call allocateStorage(M_out, storage_e_out, no_of_stats, allocateSoundVelocity)

        !$omp parallel default(shared) private(hex8Map,genHexMap,hexMap)
        ! Allocate the transfinite mappings
        allocate(hex8Map)
        allocate(genHexMap)
        
        !$omp do schedule(runtime) private(eID,i,j,k,ii,jj,kk,probe)
        do eID = 1, M_out % no_of_elements

            ! Create transfinite mapping
            if (M_out % elements(eID) % SurfInfo % IsHex8) then
               call hex8Map % setCorners(M_out % elements(eID) % SurfInfo % corners)
               hexMap => hex8Map
            else
               CALL genHexMap % destruct()
               CALL genHexMap % constructWithFaces(M_out % elements(eID) % SurfInfo % facePatches)
               hexMap => genHexMap
            end if

            do k = 0, M_out % elements(eID) % Nxyz(3)
                do j = 0, M_out % elements(eID) % Nxyz(2)
                    do i = 0, M_out % elements(eID) % Nxyz(1)
                        !
                        ! Create the node as a probe
                        !
                        call assignProbeCoordinates(probe, NodalStorage_out, M_out % elements(eID) % Nxyz, hexMap, eID, i, j, k)

                        !
                        ! Find element of M_in containing the probe
                        !
                        ! This is needed because the global variable NodalStorage is used inside the function
                        NodalStorage => NodalStorage_in
                        call findProbeInMesh(M_in, probe)

                        !
                        ! Interpolation
                        !
                        ! Allocate the polynomials for interpolation
                        safedeallocate(probe % lxi  ) ; allocate( probe % lxi(0 : M_in % elements(probe % eID) % Nxyz(1)) )
                        safedeallocate(probe % leta ) ; allocate( probe % leta(0 : M_in % elements(probe % eID) % Nxyz(2)) )
                        safedeallocate(probe % lzeta) ; allocate( probe % lzeta(0 : M_in % elements(probe % eID) % Nxyz(3)) )
                        probe % lxi = NodalStorage_in(M_in % elements(probe % eID) % Nxyz(1)) % lj(probe % xi(1))
                        probe % leta = NodalStorage_in(M_in % elements(probe % eID) % Nxyz(2)) % lj(probe % xi(2))
                        probe % lzeta = NodalStorage_in(M_in % elements(probe % eID) % Nxyz(3)) % lj(probe % xi(3))
                        ! Compute high-order interpolation
                        do kk = 0, M_in % elements(probe % eID) % Nxyz(3)
                            do jj = 0, M_in % elements(probe % eID) % Nxyz(2)
                                do ii = 0, M_in % elements(probe % eID) % Nxyz(1)
                                    storage_e_out(eID) % Q(:,i,j,k) = storage_e_out(eID) % Q(:,i,j,k) + storage_e_in(probe % eID) % Q(:,ii,jj,kk) * probe % lxi(ii) * probe % leta(jj) * probe % lzeta(kk)
                                    storage_e_out(eID) % a2(:,i,j,k) = storage_e_out(eID) % a2(:,i,j,k) + storage_e_in(probe % eID) % a2(:,ii,jj,kk) * probe % lxi(ii) * probe % leta(jj) * probe % lzeta(kk)
                                    storage_e_out(eID) % grada2(:,i,j,k) = storage_e_out(eID) % grada2(:,i,j,k) + storage_e_in(probe % eID) % grada2(:,ii,jj,kk) * probe % lxi(ii) * probe % leta(jj) * probe % lzeta(kk)
                                    storage_e_out(eID) % LambStats(:,i,j,k) = storage_e_out(eID) % LambStats(:,i,j,k) + storage_e_in(probe % eID) % LambStats(:,ii,jj,kk) * probe % lxi(ii) * probe % leta(jj) * probe % lzeta(kk)
                                    ! storage_e_out(eID) % Lamb(:,i,j,k) = storage_e_out(eID) % Lamb(:,i,j,k) + storage_e_in(probe % eID) % Lamb(:,ii,jj,kk) * probe % lxi(ii) * probe % leta(jj) * probe % lzeta(kk)
                                end do
                            end do
                        end do   

                    end do
                end do
            end do
        end do
        !$omp end do

        deallocate(hex8Map)
        deallocate(genHexMap)

        !$omp end parallel


        ! Export field values at probes
        call M_out % PrepareForIO()
        fileName = controlVariables % stringValueForKey("solution file name", requestedLength = LINE_LENGTH)
        fileName = "./RESULTS/" // trim(removePath(getFileName(fileName)))
        call saveStats(M_out, storage_e_out, refs, iter, time, fileName)

    end subroutine InterpolateStats

    subroutine ConstructSimpleMesh(mesh, meshFileName, AllNx, AllNy, AllNz)
        use SMConstants
        use HexMeshClass
        use readHDF5
        use readSpecM
        use Read_GMSH, only: CheckGMSHversion
        use readGMSH
        use FileReadingUtilities, only: getFileExtension
        implicit none
        !-arguments----------------------------------------------
        type(HexMesh)               :: mesh
        character(len=*)            :: meshFileName
        integer, intent(in)         :: AllNx(:), AllNy(:), AllNz(:)
        !-local-variables----------------------------------------
        integer                    :: gmsh_version
        character(len=LINE_LENGTH) :: ext

        ext = getFileExtension(trim(meshFileName))
        if (trim(ext)=='h5') then
            call ConstructSimpleMesh_FromHDF5File_(mesh, meshFileName, AllNx=AllNx, AllNy=AllNy, AllNz=AllNz)
        elseif (trim(ext)=='mesh') then
            call ConstructSimpleMesh_FromSpecFile_(mesh, meshFileName, AllNx=AllNx, AllNy=AllNy, AllNz=AllNz)
        elseif (trim(ext)=='msh') then
            call CheckGMSHversion (meshFileName, gmsh_version)
            select case (gmsh_version)
            case (4)
                call ConstructSimpleMesh_FromGMSHFile_v4_(mesh, meshFileName, AllNx=AllNx, AllNy=AllNy, AllNz=AllNz)
            case (2)
                call ConstructSimpleMesh_FromGMSHFile_v2_(mesh, meshFileName, AllNx=AllNx, AllNy=AllNy, AllNz=AllNz)
            case default
                error stop "ReadMeshFile :: Unrecognized GMSH version."
            end select
        else
            error stop 'Mesh file extension not recognized.'
        end if
        mesh % no_of_allElements = mesh % no_of_elements

    end subroutine ConstructSimpleMesh

    subroutine readMesh_FromHmeshFile(meshName, mesh, nodeType, Nmax)
        ! Read mesh coordinates from .hmesh file. Purely sequential code.
        use HexMeshClass
        use SolutionFile
        implicit none
        character(len=*), intent(in)     :: meshName
        class(HexMesh) :: mesh
        integer, intent(out) :: nodeType
        integer, intent(out) :: Nmax

        ! Local variables
        integer, dimension(:), allocatable     :: arrayDimensions
        integer                                :: dimensionsSize
        integer                                :: fid, eID
        integer                                :: no_of_elements, npa, ndof
        integer                                :: lb, ub

        ! Get mesh node type and type of file
        nodeType = getSolutionFileNodeType(meshName)

        ! Get number of elements
        no_of_elements = getSolutionFileNoOfElements(meshName)
        mesh % no_of_elements = no_of_elements

        ! Allocate elements
        allocate(mesh % elements(no_of_elements))
        dimensionsSize = 4
        allocate(arrayDimensions(dimensionsSize))

        ! Read coordinates
        fid = putSolutionFileInReadDataMode(meshName)
        npa = 0
        Nmax = 0
        do eID = 1, no_of_elements

            call getSolutionFileArrayDimensions(fid,arrayDimensions)

            ! Allocate memory for the coordinates
            ! arrayDimensions(2:dimensionsSize) = arrayDimensions(2:dimensionsSize) - 1
            mesh % elements(eID) % Nxyz(1:dimensionsSize-1) = arrayDimensions(2:dimensionsSize) - 1
            allocate( mesh % elements(eID) % geom % x(NDIM,0:mesh % elements(eID) % Nxyz(1),0:mesh % elements(eID) % Nxyz(2),0:mesh % elements(eID) % Nxyz(3)) )

            ! Read data
            read(fid) mesh % elements(eID) % geom % x

            ndof = (mesh % elements(eID) % Nxyz(1) + 1) * (mesh % elements(eID) % Nxyz(2) + 1) * (mesh % elements(eID) % Nxyz(3) + 1)
            npa = npa + ndof

            Nmax = max(Nmax, maxval(mesh % elements(eID) % Nxyz))

        end do

        ! Close file
        close(fid)

    end subroutine readMesh_FromHmeshFile

    subroutine getNodeType(controlVariables, nodetype)
        use mainKeywordsModule, only: discretizationNodesKey
        use NodalStorageClass, only: GAUSS, GAUSSLOBATTO
        use FTValueDictionaryClass
        implicit none
        TYPE(FTValueDictionary)                :: controlVariables
        integer, intent(out) :: nodeType

        if (.not. controlVariables % containsKey(discretizationNodesKey)) then
            call controlVariables % addValueForKey("Gauss",discretizationNodesKey)
        end if

        select case ( trim(controlVariables % stringValueForKey(trim(discretizationNodesKey), requestedLength = LINE_LENGTH)) )
            case("Gauss")
                nodeType = GAUSS
            case("Gauss-Lobatto")
                nodeType = GAUSSLOBATTO
            case default
                print*, "Unknown discretization nodes."
                print*, "Options available are:"
                print*, "   * Gauss"
                print*, "   * Gauss-Lobatto"
                errorMessage(STD_OUT)
                error stop
        end select
    end subroutine getNodeType

    subroutine write_points_vtk(filename, npa, xd_dp_pa, f_pa)
        implicit none
        ! input
        character(len=*), intent(in) :: filename
        integer, intent(in) :: npa
        real(rp), intent(in) :: xd_dp_pa(3,npa)
        real(rp), optional, intent(in) :: f_pa(:,:)

        ! local
        integer :: i
        integer :: nf
        integer :: unit

        if (present(f_pa)) nf = size(f_pa,1)

        open(newunit=unit,file=filename,status='replace',action='write',form='formatted')

        ! header
        write(unit,'(A)') '# vtk DataFile Version 3.0'
        write(unit,'(A)') 'point data'
        write(unit,'(A)') 'ASCII'
        write(unit,'(A)') 'DATASET POLYDATA'

        ! points
        write(unit,'(A,I0,A)') 'POINTS ', npa, ' double'
        do i=1,npa
            write(unit,'(3ES24.16)') xd_dp_pa(1,i),xd_dp_pa(2,i),xd_dp_pa(3,i)
        end do

        ! vertices
        write(unit,'(A,I0,1X,I0)') 'VERTICES ', npa, 2*npa
        do i=0,npa-1
            write(unit,'(2I10)') 1,i
        end do

        if (present(f_pa)) then
            ! point data
            write(unit,'(A,I0)') 'POINT_DATA ', npa

            if (nf == 1) then

                write(unit,'(A)') 'SCALARS field double'
                write(unit,'(A)') 'LOOKUP_TABLE default'

                do i=1,npa
                    write(unit,'(ES24.16)') f_pa(1,i)
                end do

            elseif (nf == 3) then

                write(unit,'(A)') 'VECTORS field double'

                do i=1,npa
                    write(unit,'(3ES24.16)') f_pa(1,i),f_pa(2,i),f_pa(3,i)
                end do

            else

                write(*,*) 'Error: f_pa must have size 1 or 3 in first dimension'
                stop

            endif
        end if

        close(unit)

    end subroutine write_points_vtk

    subroutine allocateStorage(M_out, storage_e_out, no_of_stats, allocateSoundVelocity)
        use HexMeshClass
        implicit none
        type(HexMesh), intent(in)          :: M_out
        type(StatsElementStorage_t)          :: storage_e_out(:)
        integer, intent(in) :: no_of_stats
        logical, intent(in) :: allocateSoundVelocity

        integer :: eID

        do eID = 1, M_out % no_of_elements
            ! Allocate storage for stats
            allocate(storage_e_out(eID) % Q(1:no_of_stats, 0:M_out % elements(eID) % Nxyz(1), 0:M_out % elements(eID) % Nxyz(2), 0:M_out % elements(eID) % Nxyz(3)) )
            storage_e_out(eID) % Q = 0.0_rp
            ! Allocate sound velocity if Navier-Stokes solver, i.e., if a2 or grada2 have been allocated.
            if ( allocateSoundVelocity ) then
                allocate( storage_e_out(eID) % a2(1:1, 0:M_out % elements(eID) % Nxyz(1), 0:M_out % elements(eID) % Nxyz(2), 0:M_out % elements(eID) % Nxyz(3)) )
                storage_e_out(eID) % a2 = 0.0_rp
                allocate( storage_e_out(eID) % grada2(1:NDIM, 0:M_out % elements(eID) % Nxyz(1), 0:M_out % elements(eID) % Nxyz(2), 0:M_out % elements(eID) % Nxyz(3)) )
                storage_e_out(eID) % grada2 = 0.0_rp
            end if
            ! Allocate Lamb vector stats
            allocate(storage_e_out(eID) % LambStats(1:NDIM, 0:M_out % elements(eID) % Nxyz(1), 0:M_out % elements(eID) % Nxyz(2), 0:M_out % elements(eID) % Nxyz(3)) )
            storage_e_out(eID) % LambStats = 0.0_rp
            ! Allocate Lamb vector
            ! allocate(storage_e_out(eID) % Lamb(1:NDIM, 0:M_out % elements(eID) % Nxyz(1), 0:M_out % elements(eID) % Nxyz(2), 0:M_out % elements(eID) % Nxyz(3)) )
            ! storage_e_out(eID) % Lamb = 0.0_rp
        end do
    end subroutine allocateStorage

    subroutine assignProbeCoordinates(probe, NodalStorage_out, Nxyz, hexMap, eID, i, j, k)
        use ProbeClass
        use NodalStorageClass
        use TransfiniteMapClass
        implicit none
        type(Probe_t) :: probe
        type(NodalStorage_t), pointer :: NodalStorage_out(:)
        integer, intent(in) :: Nxyz(NDIM)
        type(TransfiniteHexMap), intent(in) :: hexMap
        integer, intent(in) :: eID, i, j, k

        ! Local variables
        TYPE(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta
        real(rp) :: x(NDIM)

        probe % ID = eID
        spAxi => NodalStorage_out(Nxyz(1))
        spAeta => NodalStorage_out(Nxyz(2))
        spAzeta => NodalStorage_out(Nxyz(3))
        ! Coordinates of the point in the computational space
        x = [spAxi % x(i), spAeta % x(j), spAzeta % x(k)]
        ! Coordinates of the point in the physical space
        probe % x = hexMap % transfiniteMapAt(x)

    end subroutine assignProbeCoordinates

    subroutine findProbeInMesh(M_in, probe)
        ! Find the requested point in the mesh
        use ProbeClass
        use NodalStorageClass
        use HexMeshClass
        implicit none
        type(HexMesh) :: M_in
        type(Probe_t) :: probe

        probe % active = M_in % FindPointWithCoords(probe % x, probe % eID, probe % xi)
        ! Check whether the probe is located in another partition
        call probe % LookInOtherPartitions
        ! Disable the probe if the point is not found
        if ( .not. probe % active ) then
            ! if ( MPI_Process % isRoot ) then
                write(STD_OUT,'(A,I0,A)') "Probe ", probe % ID, " was not successfully initialized."
                print*, "Probe is set to inactive."
            ! end if
            ! return
        end if
    end subroutine findProbeInMesh

    subroutine readStats(controlVariables, mesh, storage_e, refs, iter, time, filename)
        use FTValueDictionaryClass
        use HexMeshClass
        use StorageClass
        use Utilities, only: toLower
        use SolutionFile, only: NO_OF_SAVED_REFS
        implicit none
        TYPE(FTValueDictionary)         :: controlVariables
        CLASS(HexMesh)                  :: mesh
        type(StatsElementStorage_t)          :: storage_e(mesh % no_of_elements)
        real(kind=RP), intent(out) :: refs(NO_OF_SAVED_REFS)
        integer, intent(out) :: iter
        real(rp), intent(out) :: time
        character(len=*)                :: fileName

        !
        ! Local variables
        !
        CHARACTER(LEN=LINE_LENGTH) :: qBaseSolverKey = "stats solver"
        CHARACTER(LEN=LINE_LENGTH) :: qBaseSolverNS  = "ns"
        CHARACTER(LEN=LINE_LENGTH) :: qBaseSolveriNS = "ins"
        CHARACTER(LEN=LINE_LENGTH) :: qBaseSolverMU  = "mu"
        character(len=LINE_LENGTH) :: statsSolver


        ! Check that the user has specified the solver that wrote the stats file
        call toLower(qBaseSolverKey)
        if (.not. controlVariables % containsKey(trim(qBaseSolverKey))) then
            print *, trim(qBaseSolverKey), " not specified. Use:"
            print *, trim(qBaseSolverKey), " = ns/ins/mu"
            errorMessage(STD_OUT)
            error stop
        end if
        ! Load which solver generated the stats file
        statsSolver = controlVariables % stringValueForKey(qBaseSolverKey,requestedLength = LINE_LENGTH)

        if (trim(statsSolver) .eq. trim(qBaseSolverNS)) then
            statsSolverID = statsSolverNS
            call readStats_NS(mesh, storage_e, refs, iter, time, fileName)
            call readStatsSoundVelocity_NS(controlVariables, mesh, storage_e)
            call readStatsGradSoundVelocity_NS(controlVariables, mesh, storage_e)
        elseif (trim(statsSolver) .eq. trim(qBaseSolveriNS)) then
            statsSolverID = statsSolveriNS
            call readStats_iNS(mesh, storage_e, refs, iter, time, fileName)
        elseif (trim(statsSolver) .eq. trim(qBaseSolverMU)) then
            statsSolverID = statsSolverMU
            call readStats_MU(mesh, storage_e, refs, iter, time, fileName)
        else
            print *, "Unknown solver of the stats file ", trim(statsSolver)
        end if

        ! Load Lamb vector stats
        call readStatsLambVector(controlVariables, mesh, storage_e)

        ! Load Lamb vector
        ! call readLambVector(controlVariables, mesh, storage_e)

    end subroutine readStats

    subroutine readStats_NS(mesh, storage_e, refs, iter, time, fileName)
        use HexMeshClass
        use StorageClass
        use SolutionFile
        implicit None
        CLASS(HexMesh)                  :: mesh
        type(StatsElementStorage_t) :: storage_e(mesh % no_of_elements)
        real(kind=RP), intent(out) :: refs(NO_OF_SAVED_REFS)
        integer, intent(out) :: iter
        real(rp), intent(out) :: time
        character(len=*)                :: fileName
        !
        ! Local variables
        !
        INTEGER                        :: fID, eID, fileType, no_of_elements, nodetype
        integer                        :: i, j, k, lb, ub
        integer(kind=AddrInt)          :: pos
        integer                        :: NCONS_NS, NGRAD_NS ! NCONS and NGRAD of the NS solver
        integer                        :: no_stat_s, no_stats_read


        NCONS_NS = 5
        NGRAD_NS = 5

        ! Get the file type
        fileType = getSolutionFileType(trim(fileName))
        if ( (fileType .ne. STATS_FILE) .and. (fileType .ne. STATS_AND_GRADIENTS_FILE) ) then
            print*, "The selected file is not a statistics file"
            errorMessage(STD_OUT)
            error stop
        end if
        if (fileType .eq. STATS_AND_GRADIENTS_FILE) saveGradients = .true.
        
        ! Get the node type
        nodeType = getSolutionFileNodeType(trim(fileName))
        if ( nodeType .ne. mesh % nodeType ) then
            print*, "WARNING: Stats file uses a different discretization nodes than the mesh."
            errorMessage(STD_OUT)
        end if
        
        ! Read the number of elements
        no_of_elements = getSolutionFileNoOfElements(trim(fileName))
        if ( no_of_elements .ne. mesh % no_of_elements ) then
            write(STD_OUT,'(A,A)') "The number of elements stored in the stats file ", &
                "do not match that of the mesh file"
            errorMessage(STD_OUT)
            error stop
        end if

        ! Read reference values
        refs = getSolutionFileReferenceValues(trim(fileName))
        ! Read time and iteration
        call getSolutionFileTimeAndIteration(fileName, iter, time)

        
        ! Read elements data
        fID = putSolutionFileInReadDataMode(trim(fileName))
        no_stat_s = 9
        no_stats_read = no_stat_s + NCONS_NS
        if (fileType .eq. STATS_AND_GRADIENTS_FILE) no_stats_read = no_stats_read + NGRAD_NS*NDIM
        do eID = 1, size(mesh % elements)
            associate( e => mesh % elements(eID) )
                pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*no_stats_read*e % offsetIO*SIZEOF_RP
                pos = pos + 5_AddrInt*SIZEOF_INT ! This is to skip the reading of the dimensions and shape in writeArray

                ! Allocate Q
                allocate(storage_e(eID) % Q(1:no_stats_read, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                
                ! Read and initialize velocity
                read(fID, pos=pos) storage_e(eID) % Q(1:no_stat_s,:,:,:)
                
                ! Read NCONS_NS variables 
                read(fID) storage_e(eID) % Q(no_stat_s+1:no_stat_s+NCONS_NS,:,:,:)

                ! Read gradients
                if (fileType .eq. STATS_AND_GRADIENTS_FILE) then
                    lb = no_stat_s + NCONS_NS + 1
                    ub = lb + NGRAD_NS - 1
                    read(fID) storage_e(eID) % Q(lb:ub,:,:,:)
                    lb = ub + 1
                    ub = lb + NGRAD_NS - 1
                    read(fID) storage_e(eID) % Q(lb:ub,:,:,:)
                    lb = ub + 1
                    ub = lb + NGRAD_NS - 1
                    read(fID) storage_e(eID) % Q(lb:ub,:,:,:)
                end if

            end associate
        end do

        ! Close the file
        close(fID)
    end subroutine readStats_NS

    subroutine readStatsSoundVelocity_NS(controlVariables, mesh, storage_e)
        use HexMeshClass
        use StorageClass
        use SolutionFile
        use FTValueDictionaryClass
        use Utilities, only: toLower
        implicit None
        TYPE(FTValueDictionary)         :: controlVariables
        CLASS(HexMesh)                  :: mesh
        type(StatsElementStorage_t) :: storage_e(mesh % no_of_elements)

        !
        ! Local variables
        !
        INTEGER                        :: fID, eID, fileType, no_of_elements, nodetype
        integer                        :: i, j, k, lb, ub
        integer(kind=AddrInt)          :: pos
        integer                        :: no_stat_s, no_stats_read
        character(len=LINE_LENGTH)     :: fileName
        CHARACTER(LEN=LINE_LENGTH) :: soundVelocityFileNameKey           = "stats sound velocity squared file name"


        ! Check that the user has specified the file to read from
        call toLower(soundVelocityFileNameKey)
        if (.not. controlVariables % containsKey(trim(soundVelocityFileNameKey))) then
            print *, trim(soundVelocityFileNameKey), " not specified. Use:"
            print *, trim(soundVelocityFileNameKey), " = path/to/file.SoundVelocitySquared.stats.hsol"
            errorMessage(STD_OUT)
            error stop
        end if
        fileName = controlVariables % stringValueForKey(soundVelocityFileNameKey,requestedLength = LINE_LENGTH)

        ! Get the file type
        fileType = getSolutionFileType(trim(fileName))
        if ( (fileType .ne. STATS_FILE) .and. (fileType .ne. STATS_AND_GRADIENTS_FILE) ) then
            print *, "The file ", fileName, " is not a stats file."
            errorMessage(STD_OUT)
            error stop
        end if
        
        ! Get the node type
        nodeType = getSolutionFileNodeType(trim(fileName))
        if ( nodeType .ne. mesh % nodeType ) then
            print*, "WARNING: Stats file uses a different discretization nodes than the mesh."
            errorMessage(STD_OUT)
        end if
        
        ! Read the number of elements
        no_of_elements = getSolutionFileNoOfElements(trim(fileName))
        if ( no_of_elements .ne. mesh % no_of_elements ) then
            write(STD_OUT,'(A,A)') "The number of elements stored in the stats file ", &
                "do not match that of the mesh file"
            errorMessage(STD_OUT)
            error stop
        end if

        
        ! Read elements data
        fID = putSolutionFileInReadDataMode(trim(fileName))
        do eID = 1, size(mesh % elements)
            associate( e => mesh % elements(eID) )
                pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*e % offsetIO*SIZEOF_RP
                pos = pos + 5_AddrInt*SIZEOF_INT ! This is to skip the reading of the dimensions and shape in writeArray

                ! Allocate Q
                allocate(storage_e(eID) % a2(1:1, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                
                ! Read and initialize sound velocity
                read(fID, pos=pos) storage_e(eID) % a2(1:1,:,:,:)

            end associate
        end do

        ! Close the file
        close(fID)
    end subroutine readStatsSoundVelocity_NS

    subroutine readStatsGradSoundVelocity_NS(controlVariables, mesh, storage_e)
        use HexMeshClass
        use StorageClass
        use SolutionFile
        use FTValueDictionaryClass
        use Utilities, only: toLower
        implicit None
        TYPE(FTValueDictionary)         :: controlVariables
        CLASS(HexMesh)                  :: mesh
        type(StatsElementStorage_t) :: storage_e(mesh % no_of_elements)

        !
        ! Local variables
        !
        INTEGER                        :: fID, eID, fileType, no_of_elements, nodetype
        integer                        :: i, j, k, lb, ub
        integer(kind=AddrInt)          :: pos
        integer                        :: no_stat_s, no_stats_read
        character(len=LINE_LENGTH)     :: fileName
        CHARACTER(LEN=LINE_LENGTH) :: gradSoundVelocityFileNameKey           = "stats gradient sound velocity squared file name"


        ! Check that the user has specified the file to read from
        call toLower(gradSoundVelocityFileNameKey)
        if (.not. controlVariables % containsKey(trim(gradSoundVelocityFileNameKey))) then
            print *, trim(gradSoundVelocityFileNameKey), " not specified. Use:"
            print *, trim(gradSoundVelocityFileNameKey), " = path/to/file.GradientSoundVelocitySquared.stats.hsol"
            errorMessage(STD_OUT)
            error stop
        end if
        fileName = controlVariables % stringValueForKey(gradSoundVelocityFileNameKey,requestedLength = LINE_LENGTH)

        ! Get the file type
        fileType = getSolutionFileType(trim(fileName))
        if ( (fileType .ne. STATS_FILE) .and. (fileType .ne. STATS_AND_GRADIENTS_FILE) ) then
            print *, "The file ", fileName, " is not a stats file."
            errorMessage(STD_OUT)
            error stop
        end if
        
        ! Get the node type
        nodeType = getSolutionFileNodeType(trim(fileName))
        if ( nodeType .ne. mesh % nodeType ) then
            print*, "WARNING: Stats file uses a different discretization nodes than the mesh."
            errorMessage(STD_OUT)
        end if
        
        ! Read the number of elements
        no_of_elements = getSolutionFileNoOfElements(trim(fileName))
        if ( no_of_elements .ne. mesh % no_of_elements ) then
            write(STD_OUT,'(A,A)') "The number of elements stored in the stats file ", &
                "do not match that of the mesh file"
            errorMessage(STD_OUT)
            error stop
        end if

        
        ! Read elements data
        fID = putSolutionFileInReadDataMode(trim(fileName))
        do eID = 1, size(mesh % elements)
            associate( e => mesh % elements(eID) )
                pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*NDIM*e % offsetIO*SIZEOF_RP
                pos = pos + 5_AddrInt*SIZEOF_INT ! This is to skip the reading of the dimensions and shape in writeArray

                ! Allocate Q
                allocate(storage_e(eID) % grada2(1:NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                
                ! Read and initialize sound velocity
                read(fID, pos=pos) storage_e(eID) % grada2(:,:,:,:)

            end associate
        end do

        ! Close the file
        close(fID)
    end subroutine readStatsGradSoundVelocity_NS

    subroutine readStats_iNS(mesh, storage_e, refs, iter, time, fileName)
        use HexMeshClass
        use StorageClass
        use SolutionFile
        Implicit None
        CLASS(HexMesh)                  :: mesh
        type(StatsElementStorage_t) :: storage_e(mesh % no_of_elements)
        real(kind=RP), intent(out) :: refs(NO_OF_SAVED_REFS)
        integer, intent(out) :: iter
        real(rp), intent(out) :: time
        character(len=*)                :: fileName
        ! real(kind=RP)                   :: soundVelocity

        ! Local variables
        INTEGER                        :: fID, eID, fileType, no_of_elements, nodetype
        integer                        :: i, j, k, lb, ub
        integer(kind=AddrInt)          :: pos
        integer                        :: NCONS_iNS, NGRAD_iNS ! NCONS and NGRAD of the iNS solver
        integer                        :: no_stat_s, no_stats_read

        NCONS_iNS = 5
        NGRAD_iNS = 5

        ! Get the file type
        fileType = getSolutionFileType(trim(fileName))
        print *, "fileType: ", fileType

        if ( (fileType .ne. STATS_FILE) .and. (fileType .ne. STATS_AND_GRADIENTS_FILE) ) then
            print*, "The selected file is not a statistics file"
            errorMessage(STD_OUT)
            error stop
        end if
        if (fileType .eq. STATS_AND_GRADIENTS_FILE) saveGradients = .true.
        
        ! Get the node type
        nodeType = getSolutionFileNodeType(trim(fileName))
        if ( nodeType .ne. mesh % nodeType ) then
            print*, "WARNING: Stats file uses a different discretization nodes than the mesh."
            errorMessage(STD_OUT)
        end if
        
        ! Read the number of elements
        no_of_elements = getSolutionFileNoOfElements(trim(fileName))
        if ( no_of_elements .ne. mesh % no_of_elements ) then
            write(STD_OUT,'(A,A)') "The number of elements stored in the stats file ", &
            "do not match that of the mesh file"
            errorMessage(STD_OUT)
            error stop
        end if

        ! Read reference values
        refs = getSolutionFileReferenceValues(trim(fileName))
        ! Read time and iteration
        call getSolutionFileTimeAndIteration(fileName, iter, time)

        
        ! Read elements data
        fID = putSolutionFileInReadDataMode(trim(fileName))
        no_stat_s = 9
        no_stats_read = no_stat_s + NCONS_iNS
        if (fileType .eq. STATS_AND_GRADIENTS_FILE) no_stats_read = no_stats_read + NGRAD_iNS*NDIM
        do eID = 1, size(mesh % elements)
            associate( e => mesh % elements(eID) )
                pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*no_stats_read*e % offsetIO*SIZEOF_RP
                pos = pos + 5_AddrInt*SIZEOF_INT ! This is to skip the reading of the dimensions and shape in writeArray

                ! Allocate Q
                allocate(storage_e(eID) % Q(1:no_stats_read, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                
                ! Read and initialize velocity
                read(fID, pos=pos) storage_e(eID) % Q(1:no_stat_s,:,:,:)
                
                ! Read NCONS_iNS variables 
                read(fID) storage_e(eID) % Q(no_stat_s+1:no_stat_s+NCONS_iNS,:,:,:)

                ! Read gradients
                if (fileType .eq. STATS_AND_GRADIENTS_FILE) then
                    lb = no_stat_s + NCONS_iNS + 1
                    ub = lb + NGRAD_iNS - 1
                    read(fID) storage_e(eID) % Q(lb:ub,:,:,:)
                    lb = ub + 1
                    ub = lb + NGRAD_iNS - 1
                    read(fID) storage_e(eID) % Q(lb:ub,:,:,:)
                    lb = ub + 1
                    ub = lb + NGRAD_iNS - 1
                    read(fID) storage_e(eID) % Q(lb:ub,:,:,:)
                end if
                
            end associate
        end do

        ! Close the file
        close(fID)
    end subroutine readStats_iNS

    subroutine readStats_MU(mesh, storage_e, refs, iter, time, fileName)
        use HexMeshClass
        use StorageClass
        use SolutionFile
        Implicit None
        CLASS(HexMesh)                  :: mesh
        type(StatsElementStorage_t) :: storage_e(mesh % no_of_elements)
        real(kind=RP), intent(out) :: refs(NO_OF_SAVED_REFS)
        integer, intent(out) :: iter
        real(rp), intent(out) :: time
        character(len=*)                :: fileName
        ! real(kind=RP)                   :: soundVelocity

        ! Local variables
        INTEGER                        :: fID, eID, fileType, no_of_elements, nodetype
        integer                        :: i, j, k, lb, ub
        integer(kind=AddrInt)          :: pos
        integer                        :: NCONS_MU, NGRAD_MU ! NCONS and NGRAD of the MU solver
        integer                        :: no_stat_s, no_stats_read

        NCONS_MU = 5
        NGRAD_MU = 5

        ! Get the file type
        fileType = getSolutionFileType(trim(fileName))
        print *, "fileType: ", fileType

        if ( (fileType .ne. STATS_FILE) .and. (fileType .ne. STATS_AND_GRADIENTS_FILE) ) then
            print*, "The selected file is not a statistics file"
            errorMessage(STD_OUT)
            error stop
        end if
        if (fileType .eq. STATS_AND_GRADIENTS_FILE) saveGradients = .true.
        
        ! Get the node type
        nodeType = getSolutionFileNodeType(trim(fileName))
        if ( nodeType .ne. mesh % nodeType ) then
            print*, "WARNING: Stats file uses a different discretization nodes than the mesh."
            errorMessage(STD_OUT)
        end if
        
        ! Read the number of elements
        no_of_elements = getSolutionFileNoOfElements(trim(fileName))
        if ( no_of_elements .ne. mesh % no_of_elements ) then
            write(STD_OUT,'(A,A)') "The number of elements stored in the stats file ", &
            "do not match that of the mesh file"
            errorMessage(STD_OUT)
            error stop
        end if

        ! Read reference values
        refs = getSolutionFileReferenceValues(trim(fileName))
        ! Read time and iteration
        call getSolutionFileTimeAndIteration(fileName, iter, time)

        
        ! Read elements data
        fID = putSolutionFileInReadDataMode(trim(fileName))
        no_stat_s = 9
        no_stats_read = no_stat_s + NCONS_MU + 1 ! NCONS and density
        if (fileType .eq. STATS_AND_GRADIENTS_FILE) no_stats_read = no_stats_read + NGRAD_MU*NDIM
        do eID = 1, size(mesh % elements)
            associate( e => mesh % elements(eID) )
                pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*no_stats_read*e % offsetIO*SIZEOF_RP
                pos = pos + 5_AddrInt*SIZEOF_INT ! This is to skip the reading of the dimensions and shape in writeArray

                ! Allocate Q
                allocate(storage_e(eID) % Q(1:no_stats_read, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                
                ! Read and initialize velocity
                read(fID, pos=pos) storage_e(eID) % Q(1:no_stat_s,:,:,:)
                
                ! Read NCONS_MU variables 
                read(fID) storage_e(eID) % Q(no_stat_s+1:no_stat_s+NCONS_MU,:,:,:)

                ! Read density
                read(fID) storage_e(eID) % Q(no_stat_s+NCONS_MU+1:no_stat_s+NCONS_MU+1,:,:,:)

                ! Read gradients
                if (fileType .eq. STATS_AND_GRADIENTS_FILE) then
                    lb = no_stat_s + NCONS_MU + 1 + 1
                    ub = lb + NGRAD_MU - 1
                    read(fID) storage_e(eID) % Q(lb:ub,:,:,:)
                    lb = ub + 1
                    ub = lb + NGRAD_MU - 1
                    read(fID) storage_e(eID) % Q(lb:ub,:,:,:)
                    lb = ub + 1
                    ub = lb + NGRAD_MU - 1
                    read(fID) storage_e(eID) % Q(lb:ub,:,:,:)
                end if
                
            end associate
        end do

        ! Close the file
        close(fID)
    end subroutine readStats_MU

    subroutine readStatsLambVector(controlVariables, mesh, storage_e)
        use HexMeshClass
        use StorageClass
        use SolutionFile
        use FTValueDictionaryClass
        use Utilities, only: toLower
        implicit None
        TYPE(FTValueDictionary)         :: controlVariables
        CLASS(HexMesh)                  :: mesh
        type(StatsElementStorage_t) :: storage_e(mesh % no_of_elements)

        !
        ! Local variables
        !
        INTEGER                        :: fID, eID, fileType, no_of_elements, nodetype
        integer                        :: i, j, k, lb, ub
        integer(kind=AddrInt)          :: pos
        integer                        :: no_stat_s, no_stats_read
        character(len=LINE_LENGTH)     :: fileName
        CHARACTER(LEN=LINE_LENGTH) :: LambStatsFileNameKey           = "stats Lamb vector file name"


        ! Check that the user has specified the file to read from
        call toLower(LambStatsFileNameKey)
        if (.not. controlVariables % containsKey(trim(LambStatsFileNameKey))) then
            print *, trim(LambStatsFileNameKey), " not specified. Use:"
            print *, trim(LambStatsFileNameKey), " = path/to/Lamb.stats.hsol"
            errorMessage(STD_OUT)
            error stop
        end if
        fileName = controlVariables % stringValueForKey(LambStatsFileNameKey,requestedLength = LINE_LENGTH)

        ! Get the file type
        fileType = getSolutionFileType(trim(fileName))
        if ( (fileType .ne. STATS_FILE) .and. (fileType .ne. STATS_AND_GRADIENTS_FILE) ) then
            print *, "The file ", fileName, " is not a stats file."
            errorMessage(STD_OUT)
            error stop
        end if
        
        ! Get the node type
        nodeType = getSolutionFileNodeType(trim(fileName))
        if ( nodeType .ne. mesh % nodeType ) then
            print*, "WARNING: Stats file uses a different discretization nodes than the mesh."
            errorMessage(STD_OUT)
        end if
        
        ! Read the number of elements
        no_of_elements = getSolutionFileNoOfElements(trim(fileName))
        if ( no_of_elements .ne. mesh % no_of_elements ) then
            write(STD_OUT,'(A,A)') "The number of elements stored in the stats file ", &
                "do not match that of the mesh file"
            errorMessage(STD_OUT)
            error stop
        end if

        
        ! Read elements data
        fID = putSolutionFileInReadDataMode(trim(fileName))
        do eID = 1, size(mesh % elements)
            associate( e => mesh % elements(eID) )
                pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*NDIM*e % offsetIO*SIZEOF_RP
                pos = pos + 5_AddrInt*SIZEOF_INT ! This is to skip the reading of the dimensions and shape in writeArray

                ! Allocate Q
                allocate(storage_e(eID) % LambStats(1:NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                
                ! Read and initialize sound velocity
                read(fID, pos=pos) storage_e(eID) % LambStats(:,:,:,:)

            end associate
        end do

        ! Close the file
        close(fID)
    end subroutine readStatsLambVector

    ! subroutine readLambVector(controlVariables, mesh, storage_e)
    !     use HexMeshClass
    !     use StorageClass
    !     use SolutionFile
    !     use FTValueDictionaryClass
    !     use Utilities, only: toLower
    !     implicit None
    !     TYPE(FTValueDictionary)         :: controlVariables
    !     CLASS(HexMesh)                  :: mesh
    !     type(StatsElementStorage_t) :: storage_e(mesh % no_of_elements)

    !     !
    !     ! Local variables
    !     !
    !     INTEGER                        :: fID, eID, fileType, no_of_elements, nodetype
    !     integer                        :: i, j, k, lb, ub
    !     integer(kind=AddrInt)          :: pos
    !     integer                        :: no_stat_s, no_stats_read
    !     character(len=LINE_LENGTH)     :: fileName
    !     CHARACTER(LEN=LINE_LENGTH) :: LambFileNameKey           = "Lamb vector file name"


    !     ! Check that the user has specified the file to read from
    !     call toLower(LambFileNameKey)
    !     if (.not. controlVariables % containsKey(trim(LambFileNameKey))) then
    !         print *, trim(LambFileNameKey), " not specified. Use:"
    !         print *, trim(LambFileNameKey), " = path/to/Lamb.hsol"
    !         errorMessage(STD_OUT)
    !         error stop
    !     end if
    !     fileName = controlVariables % stringValueForKey(LambFileNameKey,requestedLength = LINE_LENGTH)

    !     ! Get the file type
    !     fileType = getSolutionFileType(trim(fileName))
    !     if ( (fileType .ne. SOLUTION_FILE) ) then
    !         print *, "The file ", fileName, " is not a solution file."
    !         errorMessage(STD_OUT)
    !         error stop
    !     end if
        
    !     ! Get the node type
    !     nodeType = getSolutionFileNodeType(trim(fileName))
    !     if ( nodeType .ne. mesh % nodeType ) then
    !         print*, "WARNING: Stats file uses a different discretization nodes than the mesh."
    !         errorMessage(STD_OUT)
    !     end if
        
    !     ! Read the number of elements
    !     no_of_elements = getSolutionFileNoOfElements(trim(fileName))
    !     if ( no_of_elements .ne. mesh % no_of_elements ) then
    !         write(STD_OUT,'(A,A)') "The number of elements stored in the stats file ", &
    !             "do not match that of the mesh file"
    !         errorMessage(STD_OUT)
    !         error stop
    !     end if

        
    !     ! Read elements data
    !     fID = putSolutionFileInReadDataMode(trim(fileName))
    !     do eID = 1, size(mesh % elements)
    !         associate( e => mesh % elements(eID) )
    !             pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*NDIM*e % offsetIO*SIZEOF_RP
    !             pos = pos + 5_AddrInt*SIZEOF_INT ! This is to skip the reading of the dimensions and shape in writeArray

    !             ! Allocate Q
    !             allocate(storage_e(eID) % Lamb(1:NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                
    !             ! Read and initialize sound velocity
    !             read(fID, pos=pos) storage_e(eID) % Lamb(:,:,:,:)

    !         end associate
    !     end do

    !     ! Close the file
    !     close(fID)
    ! end subroutine readLambVector

    subroutine saveStats(mesh, storage_e, refs, iter, time, baseName)
        use HexMeshClass
        use SolutionFile
        implicit none
        class(HexMesh),      intent(in)        :: mesh
        class(StatsElementStorage_t), intent(in) :: storage_e(:)
        real(kind=RP),       intent(in)        :: refs(NO_OF_SAVED_REFS)
        integer,             intent(in)        :: iter
        real(kind=RP),       intent(in)        :: time
        character(LEN=*)                       :: baseName

        ! Local variables
        character(LEN=LINE_LENGTH) :: fileName

        write(fileName,'(A,A)') trim(baseName),'.Interpolated.stats.hsol'
        if (statsSolverID .eq. statsSolverNS) then
            call saveStatistics_NS(mesh, storage_e, refs, iter, time, fileName)
        elseif (statsSolverID .eq. statsSolveriNS) then
            call saveStatistics_iNS(mesh, storage_e, refs, iter, time, fileName)
        elseif (statsSolverID .eq. statsSolverMU) then
            call saveStatistics_MU(mesh, storage_e, refs, iter, time, fileName)
        else
            print *, "Unknown stats solver to export"
            error stop
        end if
        
        if (allocated(storage_e(1) % a2)) then
            ! Save sound velocity squared
            write(fileName,'(A,A)') trim(baseName),'.InterpolatedSoundVelocitySquared.stats.hsol'
            call saveStatsSoundVelocitySquared(mesh, storage_e, refs, iter, time, fileName)
        end if
        if (allocated(storage_e(1) % grada2)) then
            ! Save gradient sound velocity squared
            write(fileName,'(A,A)') trim(baseName),'.InterpolatedGradientSoundVelocitySquared.stats.hsol'
            call saveStatsGradientSoundVelocitySquared(mesh, storage_e, refs, iter, time, fileName)
        end if

        ! Save Lamb vector stats
        write(fileName,'(A,A)') trim(baseName),'.InterpolatedLamb.stats.hsol'
        call saveStatsLambVector(mesh, storage_e, refs, iter, time, fileName)

        ! Save Lamb vector
        ! write(fileName,'(A,A)') trim(baseName),'.InterpolatedLamb.hsol'
        ! call saveLambVector(mesh, storage_e, refs, iter, time, fileName)

    end subroutine saveStats

    subroutine saveStatistics_NS(mesh, storage_e, refs, iter, time, name)
        use HexMeshClass
        use SolutionFile
        implicit none
        class(HexMesh),      intent(in)        :: mesh
        class(StatsElementStorage_t), intent(in) :: storage_e(:)
        real(kind=RP),       intent(in)        :: refs(NO_OF_SAVED_REFS)
        integer,             intent(in)        :: iter
        real(kind=RP),       intent(in)        :: time
        character(len=*),    intent(in)        :: name
        
        
        ! Local variables
        integer                          :: fid, eID, fileType
        integer                          :: no_of_written_vars, no_stat_s
        integer                          :: NCONS_NS, NGRAD_NS ! NCONS and NGRAD of the NS solver
        integer(kind=AddrInt)            :: pos
        real(kind=RP), allocatable       :: Q(:,:,:,:)
        
        NCONS_NS = 5
        NGRAD_NS = 5
        
        ! Create new file
        fileType = STATS_FILE
        if ( saveGradients ) fileType = STATS_AND_GRADIENTS_FILE
        call CreateNewSolutionFile(trim(name), fileType, mesh % nodeType, mesh % no_of_allElements, iter, time, refs)
        
        ! Write arrays
        fID = putSolutionFileInWriteDataMode(trim(name))
        ! Set the number of written variables for the correct offset
        no_stat_s = 9 ! S_ij
        no_of_written_vars = no_stat_s + NCONS_NS ! 9 + NCONS
        if ( saveGradients ) no_of_written_vars = no_of_written_vars + NGRAD_NS * NDIM
        do eID = 1, mesh % no_of_elements
            associate( e => mesh % elements(eID) )
                pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*no_of_written_vars*e % offsetIO*SIZEOF_RP
                ! Write S_ij
                call writeArray(fid, storage_e(eID) % Q(1:no_stat_s,:,:,:), position=pos)
                ! Write state variables
                allocate(Q(NCONS_NS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                Q(1:NCONS_NS,:,:,:) = storage_e(eID) % Q(no_stat_s+1:no_stat_s+NCONS_NS,:,:,:)
                write(fid) Q
                deallocate(Q)

                ! Write gradients
                if ( saveGradients ) then
                    allocate(Q(NGRAD_NS,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                    ! UX
                    Q(1:NGRAD_NS,:,:,:) = storage_e(eID) % Q(no_stat_s+NCONS_NS+1:no_stat_s+NCONS_NS+NGRAD_NS,:,:,:)
                    write(fid) Q
                    ! UY
                    Q(1:NGRAD_NS,:,:,:) = storage_e(eID) % Q(no_stat_s+NCONS_NS+1+NGRAD_NS:no_stat_s+NCONS_NS+2*NGRAD_NS,:,:,:)
                    write(fid) Q
                    ! UZ
                    Q(1:NGRAD_NS,:,:,:) = storage_e(eID) % Q(no_stat_s+NCONS_NS+1+2*NGRAD_NS:no_stat_s+NCONS_NS+NDIM*NGRAD_NS,:,:,:)
                    write(fid) Q
                    deallocate(Q)
                end if
            end associate
        end do
        close(fid)
        
        ! Close the file
        call SealSolutionFile(trim(name))

    end subroutine saveStatistics_NS

    subroutine saveStatistics_iNS(mesh, storage_e, refs, iter, time, name)
        use HexMeshClass
        use SolutionFile
        implicit none
        class(HexMesh),      intent(in)        :: mesh
        class(StatsElementStorage_t), intent(in) :: storage_e(:)
        real(kind=RP),       intent(in)        :: refs(NO_OF_SAVED_REFS)
        integer,             intent(in)        :: iter
        real(kind=RP),       intent(in)        :: time
        character(len=*),    intent(in)        :: name
        
        
        ! Local variables
        integer                          :: fid, eID, fileType
        integer                          :: no_of_written_vars, no_stat_s
        integer                          :: NCONS_iNS, NGRAD_iNS ! NCONS and NGRAD of the iNS solver
        integer(kind=AddrInt)            :: pos
        real(kind=RP), allocatable       :: Q(:,:,:,:)
        
        NCONS_iNS = 5
        NGRAD_iNS = 5
        
        ! Create new file
        fileType = STATS_FILE
        if ( saveGradients ) fileType = STATS_AND_GRADIENTS_FILE
        call CreateNewSolutionFile(trim(name), fileType, mesh % nodeType, mesh % no_of_allElements, iter, time, refs)
        
        ! Write arrays
        fID = putSolutionFileInWriteDataMode(trim(name))
        ! Set the number of written variables for the correct offset
        no_stat_s = 9 ! S_ij
        no_of_written_vars = no_stat_s + NCONS_iNS ! 9 + NCONS
        if ( saveGradients ) no_of_written_vars = no_of_written_vars + NGRAD_iNS * NDIM
        do eID = 1, mesh % no_of_elements
            associate( e => mesh % elements(eID) )
                pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*no_of_written_vars*e % offsetIO*SIZEOF_RP
                ! Write S_ij
                call writeArray(fid, storage_e(eID) % Q(1:no_stat_s,:,:,:), position=pos)
                ! Write state variables
                allocate(Q(NCONS_iNS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                Q(1:NCONS_iNS,:,:,:) = storage_e(eID) % Q(no_stat_s+1:no_stat_s+NCONS_iNS,:,:,:)
                write(fid) Q
                deallocate(Q)

                ! Write gradients
                if ( saveGradients ) then
                    allocate(Q(NGRAD_iNS,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                    ! UX
                    Q(1:NGRAD_iNS,:,:,:) = storage_e(eID) % Q(no_stat_s+NCONS_iNS+1:no_stat_s+NCONS_iNS+NGRAD_iNS,:,:,:)
                    write(fid) Q
                    ! UY
                    Q(1:NGRAD_iNS,:,:,:) = storage_e(eID) % Q(no_stat_s+NCONS_iNS+1+NGRAD_iNS:no_stat_s+NCONS_iNS+2*NGRAD_iNS,:,:,:)
                    write(fid) Q
                    ! UZ
                    Q(1:NGRAD_iNS,:,:,:) = storage_e(eID) % Q(no_stat_s+NCONS_iNS+1+2*NGRAD_iNS:no_stat_s+NCONS_iNS+NDIM*NGRAD_iNS,:,:,:)
                    write(fid) Q
                    deallocate(Q)
                end if
            end associate
        end do
        close(fid)
        
        ! Close the file
        call SealSolutionFile(trim(name))

    end subroutine saveStatistics_iNS

    subroutine saveStatistics_MU(mesh, storage_e, refs, iter, time, name)
        use HexMeshClass
        use SolutionFile
        implicit none
        class(HexMesh),      intent(in)        :: mesh
        class(StatsElementStorage_t), intent(in) :: storage_e(:)
        real(kind=RP),       intent(in)        :: refs(NO_OF_SAVED_REFS)
        integer,             intent(in)        :: iter
        real(kind=RP),       intent(in)        :: time
        character(len=*),    intent(in)        :: name
        
        
        ! Local variables
        integer                          :: fid, eID, fileType
        integer                          :: no_of_written_vars, no_stat_s
        integer                          :: NCONS_MU, NGRAD_MU ! NCONS and NGRAD of the MU solver
        integer(kind=AddrInt)            :: pos
        real(kind=RP), allocatable       :: Q(:,:,:,:)
        
        NCONS_MU = 5
        NGRAD_MU = 5
        
        ! Create new file
        fileType = STATS_FILE
        if ( saveGradients ) fileType = STATS_AND_GRADIENTS_FILE
        call CreateNewSolutionFile(trim(name), fileType, mesh % nodeType, mesh % no_of_allElements, iter, time, refs)
        
        ! Write arrays
        fID = putSolutionFileInWriteDataMode(trim(name))
        ! Set the number of written variables for the correct offset
        no_stat_s = 9 ! S_ij
        no_of_written_vars = no_stat_s + NCONS_MU + 1 ! 9 + NCONS + rho
        if ( saveGradients ) no_of_written_vars = no_of_written_vars + NGRAD_MU * NDIM
        do eID = 1, mesh % no_of_elements
            associate( e => mesh % elements(eID) )
                pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*no_of_written_vars*e % offsetIO*SIZEOF_RP
                ! Write S_ij
                call writeArray(fid, storage_e(eID) % Q(1:no_stat_s,:,:,:), position=pos)
                ! Write state variables
                allocate(Q(NCONS_MU, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                Q(1:NCONS_MU,:,:,:) = storage_e(eID) % Q(no_stat_s+1:no_stat_s+NCONS_MU,:,:,:)
                write(fid) Q
                deallocate(Q)
                ! Write the density
                allocate(Q(1, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                Q(1:1,:,:,:) = storage_e(eID) % Q(no_stat_s+NCONS_MU+1:no_stat_s+NCONS_MU+1,:,:,:)
                write(fid) Q
                deallocate(Q)

                ! Write gradients
                if ( saveGradients ) then
                    allocate(Q(NGRAD_MU,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
                    ! UX
                    Q(1:NGRAD_MU,:,:,:) = storage_e(eID) % Q(no_stat_s+NCONS_MU+2:no_stat_s+NCONS_MU+1+NGRAD_MU,:,:,:)
                    write(fid) Q
                    ! UY
                    Q(1:NGRAD_MU,:,:,:) = storage_e(eID) % Q(no_stat_s+NCONS_MU+2+NGRAD_MU:no_stat_s+NCONS_MU+1+2*NGRAD_MU,:,:,:)
                    write(fid) Q
                    ! UZ
                    Q(1:NGRAD_MU,:,:,:) = storage_e(eID) % Q(no_stat_s+NCONS_MU+2+2*NGRAD_MU:no_stat_s+NCONS_MU+1+NDIM*NGRAD_MU,:,:,:)
                    write(fid) Q
                    deallocate(Q)
                end if
            end associate
        end do
        close(fid)
        
        ! Close the file
        call SealSolutionFile(trim(name))

    end subroutine saveStatistics_MU

    subroutine saveStatsSoundVelocitySquared(mesh, storage_e, refs, iter, time, name)
        use HexMeshClass
        use SolutionFile
        implicit none
        class(HexMesh),      intent(in)        :: mesh
        class(StatsElementStorage_t), intent(in) :: storage_e(:)
        real(kind=RP),       intent(in)        :: refs(NO_OF_SAVED_REFS)
        integer,             intent(in)        :: iter
        real(kind=RP),       intent(in)        :: time
        character(len=*),    intent(in)        :: name
        
        ! Local variables
        integer                          :: fid, eID
        integer                          :: no_stat_s
        integer(kind=AddrInt)            :: pos 
        real(kind=RP), allocatable       :: Q(:,:,:,:)

        ! Create new file
        call CreateNewSolutionFile(trim(name),STATS_FILE, mesh % nodeType, mesh % no_of_Allelements, iter, time, refs)
        
        ! Write arrays
        fID = putSolutionFileInWriteDataMode(trim(name))
        do eID = 1, mesh % no_of_elements
            pos = POS_INIT_DATA + (mesh % elements(eID) % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*mesh % elements(eID) % offsetIO*SIZEOF_RP
            call writeArray(fid, storage_e(eID) % a2, position=pos)
        end do
        close(fID)

        ! Close the file
        call SealSolutionFile(trim(name))

    end subroutine saveStatsSoundVelocitySquared

    subroutine saveStatsGradientSoundVelocitySquared(mesh, storage_e, refs, iter, time, name)
        use HexMeshClass
        use SolutionFile
        implicit none
        class(HexMesh),      intent(in)        :: mesh
        class(StatsElementStorage_t), intent(in) :: storage_e(:)
        real(kind=RP),       intent(in)        :: refs(NO_OF_SAVED_REFS)
        integer,             intent(in)        :: iter
        real(kind=RP),       intent(in)        :: time
        character(len=*),    intent(in)        :: name
        
        ! Local variables
        integer                          :: fid, eID
        integer                          :: no_stat_s
        integer(kind=AddrInt)            :: pos 
        real(kind=RP), allocatable       :: Q(:,:,:,:)

        ! Create new file
        call CreateNewSolutionFile(trim(name),STATS_FILE, mesh % nodeType, mesh % no_of_Allelements, iter, time, refs)
        
        ! Write arrays
        fID = putSolutionFileInWriteDataMode(trim(name))
        do eID = 1, mesh % no_of_elements
            pos = POS_INIT_DATA + (mesh % elements(eID) % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*NDIM*mesh % elements(eID) % offsetIO*SIZEOF_RP
            call writeArray(fid, storage_e(eID) % grada2, position=pos)
        end do
        close(fID)

        ! Close the file
        call SealSolutionFile(trim(name))

    end subroutine saveStatsGradientSoundVelocitySquared

    subroutine saveStatsLambVector(mesh, storage_e, refs, iter, time, name)
        use HexMeshClass
        use SolutionFile
        implicit none
        class(HexMesh),      intent(in)        :: mesh
        class(StatsElementStorage_t), intent(in) :: storage_e(:)
        real(kind=RP),       intent(in)        :: refs(NO_OF_SAVED_REFS)
        integer,             intent(in)        :: iter
        real(kind=RP),       intent(in)        :: time
        character(len=*),    intent(in)        :: name
        
        ! Local variables
        integer                          :: fid, eID
        integer                          :: no_stat_s
        integer(kind=AddrInt)            :: pos 
        real(kind=RP), allocatable       :: Q(:,:,:,:)

        ! Create new file
        call CreateNewSolutionFile(trim(name),STATS_FILE, mesh % nodeType, mesh % no_of_Allelements, iter, time, refs)
        
        ! Write arrays
        fID = putSolutionFileInWriteDataMode(trim(name))
        do eID = 1, mesh % no_of_elements
            pos = POS_INIT_DATA + (mesh % elements(eID) % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*NDIM*mesh % elements(eID) % offsetIO*SIZEOF_RP
            call writeArray(fid, storage_e(eID) % LambStats, position=pos)
        end do
        close(fID)

        ! Close the file
        call SealSolutionFile(trim(name))

    end subroutine saveStatsLambVector

    ! subroutine saveLambVector(mesh, storage_e, refs, iter, time, name)
    !     use HexMeshClass
    !     use SolutionFile
    !     implicit none
    !     class(HexMesh),      intent(in)        :: mesh
    !     class(StatsElementStorage_t), intent(in) :: storage_e(:)
    !     real(kind=RP),       intent(in)        :: refs(NO_OF_SAVED_REFS)
    !     integer,             intent(in)        :: iter
    !     real(kind=RP),       intent(in)        :: time
    !     character(len=*),    intent(in)        :: name
        
    !     ! Local variables
    !     integer                          :: fid, eID
    !     integer                          :: no_stat_s
    !     integer(kind=AddrInt)            :: pos 
    !     real(kind=RP), allocatable       :: Q(:,:,:,:)

    !     ! Create new file
    !     call CreateNewSolutionFile(trim(name),SOLUTION_FILE, mesh % nodeType, mesh % no_of_Allelements, iter, time, refs)
        
    !     ! Write arrays
    !     fID = putSolutionFileInWriteDataMode(trim(name))
    !     do eID = 1, mesh % no_of_elements
    !         pos = POS_INIT_DATA + (mesh % elements(eID) % globID-1)*5_AddrInt*SIZEOF_INT + 1_AddrInt*NDIM*mesh % elements(eID) % offsetIO*SIZEOF_RP
    !         call writeArray(fid, storage_e(eID) % Lamb, position=pos)
    !     end do
    !     close(fID)

    !     ! Close the file
    !     call SealSolutionFile(trim(name))

    ! end subroutine saveLambVector

end module StatsMeshInterpolation