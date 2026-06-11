title: Sliding mesh
---
[TOC]

WRITE HERE THE DOCUMENTATION FOR THE SLIDING MESHES. YOU CAN TAKE INSPIRATION FROM THE STRUCTURE BELOW.

The p-adaptation methods are used when the p-adaptation region is specified in the control file. There are two different types of p-adaptation algorithms: A Truncation Error (TE) algorithm and a Reinforcement Learning (RL) algorithm.

| Keyword               | Description                                                                                                                                                                                              | Default value  |
|-----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------|
| adaptation type       | *CHARACTER*: Can be either "TE" or "RL".                                                                                                                                                                | TE             |
| Nmax                  | *INTEGER(3)*: Maximum polynomial order in each direction for the p-adaptation algorithm (limited to 6 for RL adaptation type).                                                                       | **Mandatory keyword**  |
| Nmin                  | *INTEGER(3)*: Minimum polynomial order in each direction for the p-adaptation algorithm.                                                                                                                | [1,1,1]        |
| conforming boundaries| *CHARACTER*(*): Specifies the boundaries of the geometry that must be forced to be conforming after the p-adaptation process.                                                                           | --             |
| order across faces    | *CHARACTER*: Mathematical expression to specify the maximum polynomial order jump across faces. Currently, only \(N*2/3\) and \(N-1\) are supported.                                                      | \(N-1\)          |
| mode                  | *CHARACTER*: p-Adaptation mode. Can be *static*, *time* or *iteration*. Static p-adaptation is performed once at the beginning of a simulation for steady or unsteady simulations.                    | *static*       |
| interval              | *INTEGER/REAL*: In dynamic p-adaptation cases, this keyword specifies the iteration (integer) or time (real) interval for p-adaptation.                                                                 | *huge number*  |
| restart files         | *LOGICAL*: If .TRUE., the program writes restart files before and after the p-adaptation.                                                                                                               | .FALSE.        |
| adjust nz             | *LOGICAL*: If .TRUE., the order across faces is adjusted in the directions xi, eta, and zeta of the face (being zeta the normal direction). If .FALSE., the order is only adjusted in the xi and eta directions. | .FALSE.        |


## Truncation Error p-Adaptation


This algorithm can perform a p-adaptation to decrease the truncation error below a threshold.

```Markdown
#define p-adaptation
   adaptation type       = TE
   Truncation error type = isolated
   truncation error      = 1.d-2
   Nmax                  = [10,10,10]
   Nmin                  = [2 ,2 ,2 ]
   Conforming boundaries = [InnerCylinder,sphere]
   order across faces    = N*2/3
   increasing            = .FALSE.
   write error files     = .FALSE.
   adjust nz             = .FALSE.
   mode                  = time
   interval              = 1.d0
   restart files         = .TRUE.
   max N decrease        = 1
   padapted mg sweeps pre      = 10
   padapted mg sweeps post     = 12
   padapted mg sweeps coarsest = 20
#end
```

| Keyword                    | Description                                                                                                                                                                                     | Default value       |
|----------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------|
| truncation error type      | *CHARACTER*: Can be either "isolated" or "non-isolated".                                                                                                                                       | isolated            |
| truncation error           | *REAL*: Target truncation error for the p-adaptation algorithm.                                                                                                                                 | **Mandatory keyword** |
| coarse truncation error    | *REAL*: Truncation error used for coarsening.                                                                                                                                                  | same as truncation error |
| increasing                 | *LOGICAL*: If .TRUE. the multi-stage FMG adaptation algorithm is used.                                                                                                                         | .FALSE.             |
| write error files          | *LOGICAL*: If .TRUE., the program writes a file per element containing the directional tau-estimations. The files are stored in the folder ./*TauEstimation*/. When the simulation has several adaptation stages, the new information is just appended. | .FALSE.             |
| max N decrease             | *INTEGER*: Maximum decrease in the polynomial order in every p-adaptation procedure.                                                                                                           | \(N-N_{*min*}\) |
| post smoothing residual    | *REAL*: Specifies the maximum allowable deviation of \(\partial_t q\) after the p-adaptation procedure.                                                                                          | --                  |
| post smoothing method      | *CHARACTER*: Either RK3 or FAS.                                                                                                                                                                | RK3, if the last keyword is activated |
| estimation files           | *CHARACTER*: Name of the folder that contains the error estimations obtained with the multi tau-estimation (see the section [below](#MultiTau)).                                                        | --                  |
| estimation files number    | *INTEGER(2)*: First and last estimation stages to be used for p-adaptation.                                                                                                                    | Mandatory if last keyword is used. |
| padapted \(\ll keyword \gg\) | *MULTIPLE*: Specifies control file keywords that should be replaced after the adaptation procedure. Currently, only 'mg sweeps', 'mg sweeps pre', 'mg sweeps post', and 'mg sweeps coarsest' are supported. | --                  |


## Reinforcement Learning p-Adaptation

This algorithm can perform a p-adaptation based on a trained RL agent.

```Markdown
#define p-adaptation
   adaptation type       = RL
   agent file            = policy_padaptation/p_adaptation_policy
   tolerance             = 1d-2
   Nmax                  = [6, 6, 2]
   Nmin                  = [2, 2, 2]
   Conforming boundaries = [cylinder]
   restart files         = .FALSE.
   mode                  = iteration
   interval              = 10
   threshold             = 2.0
#end
```

| Keyword       | Description                                                                                                                     | Default value       |
|---------------|---------------------------------------------------------------------------------------------------------------------------------|---------------------|
| agent file    | *CHARACTER*: Relative path to the binary file that provides the policy of the RL agent.                                      | **Mandatory keyword** |
| tolerance     | *REAL*: Tolerance for the RL agent. The smaller the tolerance, the wider the adapted region.                                   | 0.01                |
| threshold     | *REAL*: The mesh will be adapted only if the percentage of the elements that require adaptation (in relation to the total number of elements) is above the threshold. | 0.0                 |


## Multiple truncation error estimations 
<a name="MultiTau"></a>

When using Truncation Error based adaptation, a static p-adaptation procedure can be driven by a set of error estimations, which have to be performed beforehand in a simulation with the following block:

```Markdown
#define multi tau-estimation
   truncation error type = isolated
   interval              = 10
   folder                = MultiTau
#end
```
