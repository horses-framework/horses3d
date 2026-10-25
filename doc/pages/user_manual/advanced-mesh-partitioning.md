title: Advanced mesh partitioning
---

It is possible to decide the partitioning method with a keyword specified in the control file, see [Running a simulation](running-a-simulation.html). However, sometimes it is better to have more control over the how the domain is partiotioned into the different MPI domains. For instance, in multi-phase simulations, it might be advantegous to decompose each phase separately. Another example might be simulations featuring sliding meshes, where a MPI domain does not contain elements from the static and the sliding regions.

To give full control to the user of the mesh partitioning, it is possible to define regions in the computational domain. Each element is assigned to a region (a phase in a multi-phase simulation or sliding/static part) and each region is partitioned according to some policy.

The implemented policies, controlled with the keyword `partitioning regions policy`, are:

- Fixed-ratio: the user sets in the control file the percentage of MPI domains assigned to each region.

- Proportional: the number of MPI domains assigned to each region is determined by the number of elements of each region.

- Shared: each region is partitioned into the total number of MPI domains. This is what @AbbBallout implemented for water and air.

- Custom: the user has full control and sets for each element which MPI domain it belongs to.

To implement this, the user has to code one function in the `ProblemFile.f90` of the simulation case. For fixed-ratio, proportional and shared policies, the user should implement a function to assign each element to a region. For the custom policy, the user should implement the full element-to-domain mapping.

## Fixed-ratio partitioning policy

When using this policy, the number of MPI domains assigned to a region is determined with a user-defined parameter setting the percentage of MPI domains. It is guaranteed that all the regions get at least one MPI domain.

This policy is executed when the keyword `partitioning regions policy` is set to `fixed-ratio` in the control file. The user should also specify the number of different regions present in the domain with the keyword `partitioning regions`. Additionally, the user should also specify the percentage of MPI domains assigned to each region with the keyword `partitioning regions ratio`. For instance, in a domain with two regions, and giving 30% of the MPI domains to the first regions and 70% to the second, the keywords would be
```
partitioning regions = 2
partitioning regions policy = fixed-ratio
partitioning regions ratio = [0.3, 0.7]
```

Moreover, the user has to implement a function inside the `ProblemFile.f90` of the simulation case according to the following interface:
```fortran
abstract interface
    pure integer function getRegion_f(e)
        use ElementClass
        implicit none
        type(Element), intent(in) :: e
    end function getRegion_f
end interface
```
This pure function, given an element of the mesh, returns the region to which this element is assigned to. The return value is an integer in the set \( \lbrace 1, ..., N_r \rbrace \), where \( N_r \) is the number of different regions in the domain.

An example of this policy can be found in the test `test/Multiphase/MeshPartitioningPolicies/FixedRatio`. 

## Proportional partitioning policy

When using this policy, the number of MPI domains assigned to a region is computed as \( \frac{\text{# elements in the region}}{\text{total number of elements}} \). It is guaranteed that all the regions get at least one MPI domain.

This policy is executed when the keyword `partitioning regions policy` is set to `proportional` in the control file. The user should also specify the number of different regions present in the domain with the keyword `partitioning regions`. For instance, in a domain with three regions, the keywords would be
```
partitioning regions = 3
partitioning regions policy = proportional
```

Moreover, the user has to implement a function inside the `ProblemFile.f90` of the simulation case according to the following interface:
```fortran
abstract interface
    pure integer function getRegion_f(e)
        use ElementClass
        implicit none
        type(Element), intent(in) :: e
    end function getRegion_f
end interface
```
This pure function, given an element of the mesh, returns the region to which this element is assigned to. The return value is an integer in the set \( \lbrace 1, ..., N_r \rbrace \), where \( N_r \) is the number of different regions in the domain.

An example of this policy can be found in the test `test/Multiphase/MeshPartitioningPolicies/Proportional`. 

## Shared partitioning policy

When using this policy, each region is partitioned across all the available MPI domains. 

This policy is executed when the keyword `partitioning regions policy` is set to `shared` in the control file. The user should also specify the number of different regions present in the domain with the keyword `partitioning regions`. For instance, in a domain with three regions, the keywords would be
```
partitioning regions = 3
partitioning regions policy = shared
```

Moreover, the user has to implement a function inside the `ProblemFile.f90` of the simulation case according to the following interface:
```fortran
abstract interface
    pure integer function getRegion_f(e)
        use ElementClass
        implicit none
        type(Element), intent(in) :: e
    end function getRegion_f
end interface
```
This pure function, given an element of the mesh, returns the region to which this element is assigned to. The return value is an integer in the set \( \lbrace 1, ..., N_r \rbrace \), where \( N_r \) is the number of different regions in the domain.

An example of this policy can be found in the test `test/Multiphase/MeshPartitioningPolicies/Shared`. 

## Custom partitioning policy

This policy is executed when the keyword `partitioning regions policy` is set to `custom` in the control file. 

The user has to implement a subroutine inside the `ProblemFile.f90` of the simulation case according to the following interface:
```fortran
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
```
In this subroutine, the variables `elementsDomain` and `nodesDomain` have to be filled in. The variable `elementsDomain` specifies, for each element of the mesh, the MPI domain it belongs to. Similarly, the variable `nodesDomain` specifies, for each node of the mesh, the MPI domain it belongs to. 

Moreover, in the `UserDefinedStartup` function of the `ProblemFile.f90`, a function pointer has to be set through a call to the subroutine `setPointerUserDefinedMeshPartitioning`.

An example of this policy can be found in the test `test/Multiphase/MeshPartitioningPolicies/Custom`. 