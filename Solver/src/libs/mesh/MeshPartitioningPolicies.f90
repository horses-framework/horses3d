#include "Includes.h"
module MeshPartitioningPolicies
    use SMConstants
    implicit none
    
    private
    public GetMETISElementsPartitionByPolicy
    public partitioningNumberOfRegions_KEY
    public setPointerGetRegion
    public setPointerUserDefinedMeshPartitioning

    character(len=LINE_LENGTH), parameter :: partitioningNumberOfRegions_KEY = "partitioning regions"
    character(len=LINE_LENGTH), parameter :: partitioningRegionsPolicy_KEY = "partitioning regions policy"
    character(len=LINE_LENGTH), parameter :: partitioningRegionsPolicyFixedRatio_KEY = "fixed-ratio"
    character(len=LINE_LENGTH), parameter :: partitioningRegionsPolicyProportional_KEY = "proportional"
    character(len=LINE_LENGTH), parameter :: partitioningRegionsPolicyShared_KEY = "shared"
    character(len=LINE_LENGTH), parameter :: partitioningRegionsPolicyCustom_KEY = "custom"
    integer :: partitioningRegionsPolicyType
    integer, parameter :: partitioningRegionsPolicyFixedRatioType = 1
    integer, parameter :: partitioningRegionsPolicyProportionalType = 2
    integer, parameter :: partitioningRegionsPolicySharedType = 3
    integer, parameter :: partitioningRegionsPolicyCustomType = 4
    character(len=LINE_LENGTH), parameter :: partitioningRegionsRatio_KEY = "partitioning regions ratio"

    abstract interface
        pure integer function getRegion_f(e)
            use ElementClass
            implicit none
            type(Element), intent(in) :: e
        end function getRegion_f
    end interface

    abstract interface
        subroutine userDefinedMeshPartitioning_f(mesh, no_of_domains, elementsDomain, nodesDomain, controlVariables)
            use HexMeshClass
            use FTValueDictionaryClass
            implicit none
            type(HexMesh), intent(in)              :: mesh
            integer,       intent(in)              :: no_of_domains
            integer,       intent(out)             :: elementsDomain(mesh % no_of_elements)
            integer,       intent(out)             :: nodesDomain(size(mesh % nodes))
            type(FTValueDictionary), intent(in)    :: controlVariables
        end subroutine userDefinedMeshPartitioning_f
    end interface

    procedure(getRegion_f), pointer :: getRegion => null()
    procedure(userDefinedMeshPartitioning_f), pointer :: userDefinedMeshPartitioning => null()
    
    contains

    subroutine setPointerGetRegion(user_func)
        procedure(getRegion_f) :: user_func
        getRegion => user_func
    end subroutine setPointerGetRegion

    subroutine setPointerUserDefinedMeshPartitioning(user_func)
        procedure(userDefinedMeshPartitioning_f) :: user_func
        userDefinedMeshPartitioning => user_func
    end subroutine setPointerUserDefinedMeshPartitioning

    subroutine GetMETISElementsPartitionByPolicy(mesh, no_of_domains, elementsDomain, nodesDomain, controlVariables)
        use HexMeshClass
        use SMConstants
        use MPI_Process_Info
        use FTValueDictionaryClass
        use Utilities                   , only : toLower
        use FileReadingUtilities        , only : getIntArrayFromString
        implicit none
        type(HexMesh), intent(in)              :: mesh
        integer,       intent(in)              :: no_of_domains
        integer,       intent(out)             :: elementsDomain(mesh % no_of_elements)
        integer,       intent(out)             :: nodesDomain(size(mesh % nodes))
        ! logical,       intent(in)              :: useWeights
        type(FTValueDictionary), intent(in)    :: controlVariables
        ! integer,       intent(in)              :: nLevel 
        ! integer,       intent(in)      	     :: eID_Order(mesh % no_of_elements) !AJRTODO
        ! integer,       intent(in)      	     :: nElementLevel(nLevel)

        ! Local variables
        integer, allocatable :: r_arrayPartMat(:)
        character(len=LINE_LENGTH) :: r_value
        integer :: ra, nra, ea, nea, pe, npe
        integer, allocatable :: ndom_ra(:), ner_ra(:), ra_ea(:), ndr_ra(:)
        real(kind=rp), allocatable :: rat_ra(:)

        nea = mesh % no_of_elements
        npe = 8 ! Number of points per element (=8 for hexahedron)

        !
        ! Read regions information
        !
        ! Check that getRegion pointer has been associated
        if (.not. associated(getRegion)) then
            print *, "getRegion pointer has not been associated. Associate it in the UserDefinedStartup subroutine of the ProblemFile."
            print *, "See an example in the documentation: AJRTODO."
            errorMessage(STD_OUT)
            error stop
        end if
        ! Get the number of regions
        nra = controlVariables % integerValueForKey(partitioningNumberOfRegions_KEY)
        allocate(ner_ra(nra)) ! For each region (ra), the number of elements in such region (ner)
        allocate(ra_ea(nea)) ! Stores, for each global element (ea), its region (ra)
        ner_ra = 0
        ! Count the number of elements in each region
        do ea = 1, nea
            ra = getRegion(mesh % elements(ea))
            ner_ra(ra) = ner_ra(ra) + 1
            ra_ea(ea) = ra
        end do
        if (MPI_Process % isRoot) then
            write(*,*) "From inside METIS partitioning "
            do ra = 1, nra
            write(*,'(A,I0,A,I0)') 'Identified ', ner_ra(ra), ' elements of region ', ra
            end do
        end if

        allocate(ndom_ra(nra)) ! For each region (ra), the number of domains this region should be decomposed into
        if (controlVariables % containsKey(partitioningRegionsPolicy_KEY)) then
            r_value = controlVariables % StringValueForKey(partitioningRegionsPolicy_KEY, requestedLength=LINE_LENGTH)
            call toLower(r_value)
            select case (r_value)
            case (partitioningRegionsPolicyFixedRatio_KEY)
                partitioningRegionsPolicyType = partitioningRegionsPolicyFixedRatioType
            case (partitioningRegionsPolicyProportional_KEY)
                partitioningRegionsPolicyType = partitioningRegionsPolicyProportionalType
            case (partitioningRegionsPolicyShared_KEY)
                partitioningRegionsPolicyType = partitioningRegionsPolicySharedType
            case (partitioningRegionsPolicyCustom_KEY)
                partitioningRegionsPolicyType = partitioningRegionsPolicyCustomType
            case default
                print *, "Unknown partition policy: ", r_value
                print *, "Implemeted policies are:"
                print *, "* ", partitioningRegionsPolicyFixedRatio_KEY
                print *, "* ", partitioningRegionsPolicyProportional_KEY
                print *, "* ", partitioningRegionsPolicyShared_KEY
                print *, "* ", partitioningRegionsPolicyCustom_KEY
                errorMessage(STD_OUT)
                error stop
            end select
        else
            ! Default is proportional
            partitioningRegionsPolicyType = partitioningRegionsPolicyProportionalType
        end if
        if (partitioningRegionsPolicyType == partitioningRegionsPolicyCustomType) then
            ! AJRTODO: Test
            if (.not. associated(userDefinedMeshPartitioning)) then
                print *, "userDefinedMeshPartitioning pointer has not been associated. Associate it in the UserDefinedStartup subroutine of the ProblemFile."
                print *, "See an example in the documentation: AJRTODO."
                errorMessage(STD_OUT)
                error stop
            end if
            call userDefinedMeshPartitioning(mesh, no_of_domains, elementsDomain, nodesDomain, controlVariables)
        else
            allocate(ndr_ra(nra))
            select case (partitioningRegionsPolicyType)
            case (partitioningRegionsPolicyFixedRatioType)
                ! AJRTODO: Test
                rat_ra = getIntArrayFromString(controlVariables % stringValueForKey(partitioningRegionsRatio_KEY, requestedLength=LINE_LENGTH))
                if (size(rat_ra) .ne. nra) then
                    print *, "The size of the vector of ratios does not match the number of regions."
                    print *, "Number of regions: ", nra
                    print *, "Size of vector of ratios: ", size(rat_ra)
                    errorMessage(STD_OUT)
                    error stop
                end if
                call compute_fixed_ratio_partition(nra, no_of_domains, rat_ra, ndr_ra)
            case (partitioningRegionsPolicyProportionalType)
                call compute_proportional_partition(nra, no_of_domains, ner_ra, ndr_ra)
            case (partitioningRegionsPolicySharedType)
                do ra = 1, nra
                    ndr_ra(ra) = no_of_domains
                end do
            end select
            ! Partition the elements in each region
            print *, "ndr_ra: ", ndr_ra
            elementsDomain = -1
            do ra = 1, nra
                call partition_region_METIS(mesh, ra, ra_ea, ner_ra(ra), ndr_ra(ra), elementsDomain)
            end do
            if (any(elementsDomain == -1)) then
                print *, "Some elements have not been assigned to a MPI domain."
                errorMessage(STD_OUT)
                error stop
            end if
            ! Assign the nodes to each MPI domain
            nodesDomain = -1
            do ea = 1, nea
                do pe = 1, npe
                    nodesDomain(mesh % elements(ea) % nodeIDs(pe)) = elementsDomain(ea)
                end do
            end do
            if (any(nodesDomain == -1)) then
                print *, "Some nodes have not been assigned to a MPI domain."
                errorMessage(STD_OUT)
                error stop
            end if
        end if

    end subroutine GetMETISElementsPartitionByPolicy


    !====================================================================
    ! Allocate MPI ranks according to arbitrary positive weights.
    !
    ! Every region with a positive weight receives at least one rank.
    ! The remaining ranks are distributed proportionally using the
    ! largest remainder method.
    !====================================================================
    subroutine apportion_ranks_from_weights(n_regions, n_ranks, weight, region_ranks)
        use SMConstants
        implicit none

        integer, intent(in) :: n_regions
        integer, intent(in) :: n_ranks
        real(kind=rp), intent(in) :: weight(n_regions)

        integer, intent(out) :: region_ranks(n_regions)

        real(kind=rp), parameter :: tol = 1.0e-12_rp

        real(kind=rp) :: ratio(n_regions)
        real(kind=rp) :: exact(n_regions)
        real(kind=rp) :: remainder(n_regions)

        real(kind=rp) :: total_weight

        integer :: n_nonzero
        integer :: n_remaining
        integer :: assigned
        integer :: i, j

        !---------------------------------------------------------------
        ! Compute normalized weights.
        !---------------------------------------------------------------

        total_weight = sum(weight)

        ratio = weight / total_weight

        !---------------------------------------------------------------
        ! Give one rank to every non-empty region.
        !---------------------------------------------------------------

        region_ranks = 0

        do i = 1, n_regions
            if (weight(i) > tol) then
                region_ranks(i) = 1
            end if
        end do

        n_nonzero   = count(weight > tol)
        n_remaining = n_ranks - n_nonzero

        if (n_remaining == 0) return

        !---------------------------------------------------------------
        ! Allocate remaining ranks.
        !---------------------------------------------------------------

        assigned = 0

        do i = 1, n_regions

            if (weight(i) > tol) then

                exact(i) = ratio(i) * dble(n_remaining)

                region_ranks(i) = region_ranks(i) + floor(exact(i))

                remainder(i) = exact(i) - floor(exact(i))

                assigned = assigned + floor(exact(i))

            else

                remainder(i) = -1.0_rp

            end if

        end do

        do while (assigned < n_remaining)

            j = maxloc(remainder, dim=1)

            region_ranks(j) = region_ranks(j) + 1

            remainder(j) = -1.0_rp

            assigned = assigned + 1

        end do

    end subroutine apportion_ranks_from_weights

    !====================================================================
    ! Evenly distribute a set of identical items among a fixed number of
    ! bins.
    !
    ! The allocation is as balanced as possible: the number of items
    ! assigned to any two bins differs by at most one. If the items cannot
    ! be divided exactly, the first bins receive one additional item.
    !
    ! Examples:
    !   10 items, 4 bins -> [3,3,2,2]
    !    3 items, 5 bins -> [1,1,1,0,0]
    !====================================================================
    subroutine apportion_items(n_items, n_bins, items_per_bin)

    implicit none

    integer, intent(in)  :: n_items
    integer, intent(in)  :: n_bins
    integer, intent(out) :: items_per_bin(n_bins)

    integer :: q
    integer :: r

    q = n_items / n_bins
    r = mod(n_items, n_bins)

    items_per_bin = q

    if (r > 0) then
        items_per_bin(1:r) = items_per_bin(1:r) + 1
    end if

    end subroutine

    !====================================================================
    ! FixedRatio policy.
    !====================================================================
    subroutine compute_fixed_ratio_partition(nra, nda, rat_ra, ndr_ra)
        use SMConstants
        implicit none

        integer, intent(in) :: nra ! Number of regions (nra)
        integer, intent(in) :: nda ! Number of total MPI domains (nda)
        real(kind=rp), intent(in) :: rat_ra(nra) ! Ratio of elements for each region (ra)
        integer, intent(out) :: ndr_ra(nra) ! For each region (ra), the number of domains it got assigned (ndr)

        real(kind=rp), parameter :: tol = 1.0e-12_rp

        if (any(rat_ra < 0.0_rp)) then
            error stop "FixedRatio: negative ratios are not allowed."
        end if

        if (abs(sum(rat_ra) - 1.0_rp) > tol) then
            error stop "FixedRatio: ratios must sum to one."
        end if

        if (count(rat_ra > tol) > nda) then
            error stop "FixedRatio: more non-empty regions than MPI ranks."
        end if

        call apportion_ranks_from_weights( &
            nra, nda, rat_ra, ndr_ra)

    end subroutine compute_fixed_ratio_partition

    !====================================================================
    ! Proportional policy.
    !====================================================================
    subroutine compute_proportional_partition(nra, nda, ner_ra, ndr_ra)
        use SMConstants
        implicit none

        integer, intent(in) :: nra ! Number of regions (nra)
        integer, intent(in) :: nda ! Number of total MPI domains (nda)
        integer, intent(in) :: ner_ra(nra) ! For each region (ra), the number of elements in such region (ner)
        integer, intent(out) :: ndr_ra(nra) ! For each region (ra), the number of domains it got assigned (ndr)

        if (any(ner_ra < 0)) then
            error stop "Proportional: negative element counts are not allowed."
        end if

        if (sum(ner_ra) == 0) then
            error stop "Proportional: mesh contains no elements."
        end if

        if (count(ner_ra > 0) > nda) then
            error stop "Proportional: more non-empty regions than MPI ranks."
        end if

        call apportion_ranks_from_weights( &
            nra, nda, real(ner_ra, kind=rp), ndr_ra)

    end subroutine compute_proportional_partition

    !====================================================================
    ! Shared policy. !AJRTODO: Delete?
    !====================================================================
    subroutine compute_shared_partition(nra, nda, ner_ra, ne_ra_da)

    implicit none

    integer, intent(in)  :: nra ! Number of regions (nra)
    integer, intent(in)  :: nda ! Number of total MPI domains (nda)
    integer, intent(in)  :: ner_ra(nra) ! For each region (ra), the number of elements in such region (ner)
    integer, intent(out) :: ne_ra_da(nra,nda) ! The number of elements (ne) assigned from each region (ra) to each MPI domain (da)

    integer :: ra

    if (any(ner_ra < 0)) then
        error stop "Shared: negative element counts."
    end if

    do ra = 1, nra

        call apportion_items( &
            ner_ra(ra), &
            nda, &
            ne_ra_da(ra,:))

    end do

    end subroutine

    !====================================================================
    ! Partition element list with METIS
    !====================================================================
    subroutine partition_region_METIS(mesh, ra, ra_ea, ner, ndr, da_ea)
        use HexMeshClass
        use SMConstants
        implicit none
        type(HexMesh), intent(in) :: mesh
        integer, intent(in) :: ra ! The current region
        integer, intent(in) :: ra_ea(mesh % no_of_elements) ! For each element (ea), the region it belongs to (ra)
        integer, intent(in) :: ner ! The number of elements in the current region
        integer, intent(in) :: ndr ! The number of MPI domains that this region should be divided into
        integer, intent(out) :: da_ea(mesh % no_of_elements) ! For each element (ea), the MPI domain it belongs to (da)

        ! Local variables
        integer, pointer :: vwgt(:) => null(), vsize(:) => null()
        real(kind=RP), pointer :: tpwgt(:) => null()
        real(kind=rp) :: objval
        integer :: nvertex = 8, ncommon = 4
        integer, allocatable :: eptr(:), eind(:), opts(:)
        integer :: ea, nea, npa, elemind, nodeind
        integer, allocatable :: da_er(:), da_pa(:)

        nea = mesh % no_of_elements
        npa = size(mesh % nodes)

        
        ! Build connectivity for the elements in the region
        allocate(eptr(ner+1), eind(nvertex*ner))
        
        elemind = 1
        nodeind = 1
        do ea = 1, nea
            if (ra_ea(ea) == ra) then ! If the element belongs to the current region
                eptr(elemind) = nodeind - 1
                eind(nodeind:nodeind+nvertex-1) = mesh % elements(ea) % nodeIDs - 1
                nodeind = nodeind + nvertex
                elemind = elemind + 1
            end if
        end do
        eptr(elemind) = nodeind - 1

        
        ! METIS options
        allocate(opts(0:39))
        call METIS_SetDefaultOptions(opts)
        opts(1) = -1  ! MIN_EDGE_CUT
        opts(5) = 0   ! No verbosity
        
        ! Partition subset
        allocate(da_er(ner), da_pa(npa)) ! The domain (da) for each element of the region (er) or mesh node (pa)
        call METIS_PartMeshDual(ner, npa, eptr, eind, vwgt, vsize, &
                            ncommon, ndr, tpwgt, opts, objval, da_er, da_pa)
        
        
        ! Assign to domain in the global list of elements
        elemind = 1
        do ea = 1, nea
            if (ra_ea(ea) == ra) then ! If the element belongs to the current region
                da_ea(ea) = da_er(elemind) + 1 ! METIS starts with 0
                elemind = elemind + 1
            end if
        end do
        deallocate(eptr, eind, opts, da_er, da_pa)
    end subroutine partition_region_METIS

end module MeshPartitioningPolicies