# Sliding mesh

The sliding mesh feature takes a conforming cylindrical mesh and, if activated, builds a sliding interface at the specified radius. It then performs the prescribed rotation and generates the mortars needed to handle the resulting non-conforming interface between the fixed and rotating regions.

| Keyword               | Description                                                                                                                                                                                              | Default value  |
|-----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------|
| center                | *REAL(2)*: The center of rotation.                                                                                                                                                                        | **Mandatory keyword** |
| radius                | *REAL*: Radius defining the sliding mesh interface.                                                                                                                                                       | **Mandatory keyword** |
| angle                 | *REAL*: Rotation angle (radians) applied to the rotating region at every time step, once the sliding mesh is active.                                                                                     | **Mandatory keyword** |
| rotation axis         | *CHARACTER*: Rotation axis. Possible values are `x`, `y`, or `z`.                                                                                                                                         | **Mandatory keyword** |
| autosave mode         | *CHARACTER*: Can be `iteration` or `time`. Enables autosaving of the sliding mesh state at the specified interval.                                                                                       | Disabled       |
| autosave interval     | *INTEGER/REAL*: Iteration (integer) or time (real) interval for autosaving, depending on `autosave mode`.                                                                                                | --             |

```
#define SlidingMesh
  center = [0.0d0, 0.d0]
  radius = 1.01d0
  angle = 0.07853981d0
  rotation axis = y
#end
```
