subroutine GetMETISElementsPartition(mesh, no_of_domains, elementsDomain, nodesDomain, useWeights, controlVariables, &
										nLevel, eID_Order, nElementLevel)
!
!     *********************************************************************
!        This subroutine performs the mesh partitioning using METIS
!        with two-phase partitioning for MixedRK multiphase simulations
!     *********************************************************************
!        
      use HexMeshClass
      use SMConstants
      use MPI_Process_Info
      use FTValueDictionaryClass
	  use Utilities                   , only : reindexIntegerList, combine_partitions
      use MeshPartitioningPolicies
      implicit none
      type(HexMesh), intent(in)              :: mesh
      integer,       intent(in)              :: no_of_domains
      integer,       intent(out)             :: elementsDomain(mesh % no_of_elements)
      integer,       intent(out)             :: nodesDomain(size(mesh % nodes))
      logical,       intent(in)              :: useWeights
      type(FTValueDictionary), intent(in)    :: controlVariables
	  integer,       intent(in)              :: nLevel 
	  integer,       intent(in)      	     :: eID_Order(mesh % no_of_elements)
	  integer,       intent(in)      	     :: nElementLevel(nLevel)
!
!     ---------------
!     Local Variables
!     ---------------
!
      integer                :: i, j, domain_id, level, starteID
      integer                :: ielem
      integer                :: ne                    ! # elements
      integer                :: nn                    ! # nodes
      integer                :: nvertex               ! # vertices per element
      integer, allocatable   :: eptr(:)               ! index in eind --> eptr(ne+1)
      integer, allocatable   :: eind(:)               ! vertices of each element   --> eind(nvertex*ne)
      integer, pointer       :: vwgt(:) => null()     ! vertices weights
      integer, allocatable   :: vsize(:) 
      integer, parameter     :: ncommon = 4           ! common faces for dual nesh (4 for hexahedrals)
      real(kind=RP), allocatable  :: tpwgt(:)         ! partitions' weights --> tpwgt(no_of_domains)
      integer, allocatable   :: opts(:)               ! metis options
      integer, parameter     :: metis_noptions = 39   ! number of metis options
      integer                :: objval                ! objective calculated value
      integer, parameter     :: MIN_EDGE_CUT = -1     ! option to minimize edge-cut
      integer, parameter     :: MIN_COMM_VOL = 1      ! option to minimize the communication volume
      integer, allocatable, target   :: weights(:)
      logical                :: is_water              ! Flag for water elements
      integer                :: water_count, air_count ! Counters for water and air elements
      logical                :: do_water_air_count    ! Flag to control water/air counting
	  integer, allocatable   :: elementsDomainLevel(:), nodesDomainLevel(:), mapToOld(:)
	  integer, allocatable   :: refEleDomain(:), inputMLRKDomain(:,:)
	  
#ifdef _HAS_METIS_
!
!     **************
!     Initialization    
!     **************
!
      ne = mesh % no_of_elements
      nn = size(mesh % nodes)
      nvertex = 8
	  
	  allocate(inputMLRKDomain(mesh % no_of_elements, nLevel))
	  inputMLRKDomain = 0
	
      ! Check if we should do water/air counting based on explicit method
      do_water_air_count = .false.
      if(trim(controlVariables % stringValueForKey('explicit method', requestedLength = LINE_LENGTH)) == 'MixedRK') then
          do_water_air_count = .true.
      end if

      ! Count water and air elements only if needed
      water_count = 0
      air_count = 0
      if (do_water_air_count) then
          do ielem=1,ne
              is_water = .false.
              if (associated(mesh % elements(ielem) % storage)) then
                    is_water = all(mesh % elements(ielem) % storage % Q(1,:,:,:) < 1.0_RP - 1e-8_RP)
              end if

              if (is_water) then
                  water_count = water_count + 1
              else
                  air_count = air_count + 1
              end if
          end do

          if (MPI_Process % isRoot) then
              write(*,*) "From inside METIS partitioning "
              write(*,'(A,I6,A,I6,A)') 'Identified ', water_count, ' water elements and ', air_count, ' air elements'
          end if
      end if

      !
      ! Partitioning policy
      !
      if (controlVariables % containsKey(partitioningNumberOfRegions_KEY)) then
        call GetMETISElementsPartitionByPolicy(mesh, no_of_domains, elementsDomain, nodesDomain, controlVariables)
        return
      end if
!
!     Construct the reference element partition from single level for MLRK
!     --------------------------------------------------------------------      
	  if (nLevel.gt.1) then
	    allocate(refEleDomain(ne))
		allocate(eptr(ne+1), eind(nvertex*ne))
		
		refEleDomain    =0
!
!     	Gather each element nodes: they are stored using C indexes (starting on 0)
!     	-------------------------
		i = 1
		do ielem=1,ne
!
!        	Save each element nodes
!        	-----------------------
			eind(i:i+nvertex-1) = mesh % elements(ielem) % nodeIDs - 1
!
!        	Save each element ID
!        	--------------------
			eptr(ielem) = i - 1
			i = i + nvertex
		end do
!
!     	Termination: set the last element position in the N+1 entry
!     	-----------
		eptr(ne+1) = i - 1
!
!     	*****************
!     	Set METIS options
!     	*****************
!
		allocate(opts(0:metis_noptions))
		call METIS_SetDefaultOptions(opts)
!
!     	First option chooses the method: edge-cut / communication volume      
!     	-------------------------------
		opts(1) = MIN_EDGE_CUT
!
!     	Disable verbosity
!     	-----------------
		opts(5) = 0
!
!     	*******************************
!     	Calculate weights based on NDOF
!     	*******************************
!
		if (useWeights) then
		  allocate(weights(ne))
		  do ielem=1,ne
			  ! Base weight from polynomial order
			  weights(ielem) = product(mesh % elements(ielem) % Nxyz + 1)
		  end do
		  if (maxval(weights) .ne. minval(weights)) then
			  vwgt => weights
		  endif
		end if 
	    
		call METIS_PartMeshDual(ne, nn,  eptr, eind,  vwgt, vsize, ncommon, no_of_domains, tpwgt,  opts,  objval, refEleDomain, nodesDomain)
!
!     	Recover FORTRAN displacements by adding 1
!     	-----------------------------------------
		refEleDomain = refEleDomain + 1
!
!       Free memory
!       -----------
!
		deallocate (eptr, eind, opts)
		if (associated(vwgt)) nullify(vwgt)       ! vwgt is a pointer to a target. nullify is enough
		if (allocated(weights)) deallocate(weights)
		if (allocated(tpwgt)) deallocate(tpwgt)
		if (allocated(vsize)) deallocate(vsize)
	  end if 
!
!     Perform METIS partitioning based on the element's level - if not MLRK then nLevel=1
!     -----------------------------------------------------------------------------------
	  starteID = 1								   ! counter of the elementID for multi level 
	  do level =1, nLevel
		ne = nElementLevel(level)
		allocate(eptr(ne+1))
		allocate(eind(nvertex*ne))
!
!     	Gather each element nodes: they are stored using C indexes (starting on 0)
!     	-------------------------------------------------------------------------
		i = 1
		j = 0
		do ielem=starteID, ne+starteID-1
		    j=j+1
!
!        	Save each element nodes
!        	-----------------------
			eind(i:i+nvertex-1) = mesh % elements(eID_Order(ielem)) % nodeIDs - 1
!
!        	Save each element ID
!        	--------------------
			eptr(j) = i - 1
			i = i + nvertex
		end do
!
!     	reindexIntegerList as such it is compactly renumbered from 0 to totalNodes-1
!     	----------------------------------------------------------------------------
        allocate(mapToOld(size(mesh % nodes)))
		mapToOld = [(i, i=1,size(mesh % nodes))]
		if (nLevel.gt.1) then
		    deallocate(mapToOld)
			call reindexIntegerList(eind,nn,mapToOld)
		end if
		allocate(elementsDomainLevel(1:ne), nodesDomainLevel(1:nn))
!
!     	Termination: set the last element position in the N+1 entry
!     	-----------------------------------------------------------
		eptr(ne+1) = i - 1
!
!     	*****************
!     	Set METIS options
!     	*****************
!
		allocate(opts(0:metis_noptions))
		call METIS_SetDefaultOptions(opts)
!
!     	First option chooses the method: edge-cut / communication volume      
!     	----------------------------------------------------------------
		opts(1) = MIN_EDGE_CUT
!
!     	Disable verbosity
!     	-----------------
		opts(5) = 0
!
!     	************************************
!     	Two-phase partitioning for MixedRK
!     	************************************
!
		if (do_water_air_count .and. water_count > 0 .and. air_count > 0) then
!
!     		Separate air and water elements and distribute across all domains
!     		This ensures every domain gets both air and water elements
!     		----------------------------------------------------------------
			call partition_air_water_separately(mesh, ne, starteID, eID_Order, no_of_domains, &
												water_count, air_count, elementsDomainLevel)
		else
!
!     		*******************************
!     		Calculate weights based on NDOF
!     		*******************************
!
			if (useWeights) then
			  allocate(weights(ne))
			  do ielem=1,ne
				  ! Base weight from polynomial order
				  weights(ielem) = product(mesh % elements(eID_Order(ielem+starteID-1)) % Nxyz + 1)
			  end do
			  if (maxval(weights) .ne. minval(weights)) then
				  vwgt => weights
			  endif
			end if 
!
!     		Regular METIS partitioning for non-MixedRK cases
!     		-------------------------------------------------
			call METIS_PartMeshDual(ne, nn, eptr, eind, vwgt, vsize, ncommon, no_of_domains, tpwgt, opts, objval, elementsDomainLevel, nodesDomainLevel)
		end if
!
!     	Recover FORTRAN displacements by adding 1
!     	-----------------------------------------
		elementsDomainLevel = elementsDomainLevel + 1
		if (.not. (do_water_air_count .and. water_count > 0 .and. air_count > 0)) then
			nodesDomainLevel = nodesDomainLevel + 1
		end if
!
!     	Send information to inputMLRKDomain - relocate to true elementID
!     	----------------------------------------------------------------
        do i=1, ne
			inputMLRKDomain(eID_Order(starteID+i-1), level) = elementsDomainLevel(i)
		end do 
		starteID = starteID + ne
!
!       Free memory
!       -----------
!
		deallocate (elementsDomainLevel, mapToOld, eptr, eind, opts)
		if (allocated(nodesDomainLevel)) deallocate(nodesDomainLevel)
		if (associated(vwgt)) 	nullify(vwgt)       ! vwgt is a pointer to a target. nullify is enough
		if (allocated(weights)) deallocate(weights)
		if (allocated(tpwgt)) 	deallocate(tpwgt)
		if (allocated(vsize)) 	deallocate(vsize)
	 end do
	 

	 do i=1, mesh % no_of_elements
		elementsDomain(i) = sum(inputMLRKDomain(i,:))
	 end do 
!
!    Combine_partitions between MLRK levels as such it has minimum MPI faces
!    -----------------------------------------------------------------------
	 if (nLevel.gt.1) then
		call combine_partitions(mesh % no_of_elements, nLevel, no_of_domains, refEleDomain, inputMLRKDomain, elementsDomain)
	 end if 
!
!    Update node domains based on final element assignments
!    ------------------------------------------------------
     nodesDomain = 1
     do ielem = 1, mesh % no_of_elements
         do i = 1, nvertex
             nodesDomain(mesh % elements(ielem) % nodeIDs(i)) = elementsDomain(ielem)
         end do
     end do
!
!     **********************
!     Print domain statistics
!     **********************
!
	  if (do_water_air_count .and. MPI_Process % isRoot) then
          do domain_id = 1, no_of_domains
              water_count = 0
              air_count = 0
              do ielem = 1, mesh % no_of_elements
                  if (elementsDomain(ielem) == domain_id) then
                      ! Determine if this is a water element
                      is_water = .false.
                      if (associated(mesh % elements(ielem) % storage)) then
                          is_water = all(mesh % elements(ielem) % storage % Q(1,:,:,:) < 1.0_RP - 1e-8_RP)
                      end if
                      if (is_water) then
                          water_count = water_count + 1
                      else
                          air_count = air_count + 1
                      end if
                  end if
              end do
              write(*,'(A,I3,A,I6,A,I6,A,I6)') 'Domain ', domain_id, ': Water elements = ', &
                   water_count, ', Air elements = ', air_count, ', Total = ', water_count + air_count
          end do
      end if

contains

subroutine partition_air_water_separately(mesh, ne, starteID, eID_Order, no_of_domains, &
                                         water_count, air_count, elementsDomainLevel)
    implicit none
    type(HexMesh), intent(in) :: mesh
    integer, intent(in) :: ne, starteID, no_of_domains, water_count, air_count
    integer, intent(in) :: eID_Order(:)
    integer, intent(out) :: elementsDomainLevel(:)
    
    integer, allocatable :: air_elems(:), water_elems(:)
    integer, allocatable :: air_domains(:), water_domains(:)
    integer :: ielem, w_idx, a_idx, i
    logical :: is_water
    
!   Collect air and water element indices
    allocate(air_elems(air_count), water_elems(water_count))
    w_idx = 0
    a_idx = 0
    
    do ielem = 1, ne
        is_water = .false.
        if (associated(mesh % elements(eID_Order(ielem+starteID-1)) % storage)) then
            is_water = all(mesh % elements(eID_Order(ielem+starteID-1)) % storage % Q(1,:,:,:) < 1.0_RP - 1e-8_RP)
        end if
        
        if (is_water) then
            w_idx = w_idx + 1
            water_elems(w_idx) = ielem
        else
            a_idx = a_idx + 1
            air_elems(a_idx) = ielem
        end if
    end do

!   Partition air elements with METIS
    allocate(air_domains(air_count))
    call partition_subset_with_metis(mesh, air_elems, air_count, starteID, eID_Order, no_of_domains, air_domains)
    
!   Partition water elements with METIS
    allocate(water_domains(water_count))
    call partition_subset_with_metis(mesh, water_elems, water_count, starteID, eID_Order, no_of_domains, water_domains)

!   Simply assign METIS results directly
    elementsDomainLevel = 0
    do i = 1, air_count
        elementsDomainLevel(air_elems(i)) = air_domains(i)
    end do
    do i = 1, water_count
        elementsDomainLevel(water_elems(i)) = water_domains(i)
    end do
    
    deallocate(air_elems, water_elems, air_domains, water_domains)
end subroutine partition_air_water_separately

subroutine partition_subset_with_metis(mesh, elem_subset, subset_count, starteID, eID_Order, no_of_domains, domains)
!   Partition a subset of elements using METIS for optimal spatial locality
    implicit none
    type(HexMesh), intent(in) :: mesh
    integer, intent(in) :: subset_count, starteID, no_of_domains
    integer, intent(in) :: elem_subset(:), eID_Order(:)
    integer, intent(out) :: domains(:)
    
    integer, allocatable :: eptr(:), eind(:), opts(:), subset_domains(:), subset_nodes(:)
    integer, pointer :: vwgt(:) => null(), vsize(:) => null()
    real(kind=RP), pointer :: tpwgt(:) => null()
    integer :: i, j, objval, nvertex = 8, ncommon = 4
    
    if (subset_count == 0) return
    
!   Build connectivity for subset elements only
    allocate(eptr(subset_count+1), eind(nvertex*subset_count))
    
    j = 1
    do i = 1, subset_count
        eptr(i) = j - 1  ! 0-based indexing
        eind(j:j+nvertex-1) = mesh % elements(eID_Order(elem_subset(i)+starteID-1)) % nodeIDs - 1
        j = j + nvertex
    end do
    eptr(subset_count+1) = j - 1
    
!   METIS options
    allocate(opts(0:39))
    call METIS_SetDefaultOptions(opts)
    opts(1) = -1  ! MIN_EDGE_CUT
    opts(5) = 0   ! No verbosity
    
!   Partition subset
    allocate(subset_domains(subset_count), subset_nodes(size(mesh % nodes)))
    call METIS_PartMeshDual(subset_count, size(mesh % nodes), eptr, eind, vwgt, vsize, &
                           ncommon, no_of_domains, tpwgt, opts, objval, subset_domains, subset_nodes)
    
!   Copy results
    domains(1:subset_count) = subset_domains(1:subset_count)
    
    deallocate(eptr, eind, opts, subset_domains, subset_nodes)
end subroutine partition_subset_with_metis


#else
      write(STD_OUT,*) "METIS library was not linked."
      error stop
#endif

   END SUBROUTINE GetMETISElementsPartition
