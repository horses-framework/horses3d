title: Sliding mesh
---
[TOC]

WRITE HERE THE DOCUMENTATION FOR THE SLIDING MESHES. YOU CAN TAKE INSPIRATION FROM THE STRUCTURE BELOW.

TODO @Hatem

| Keyword               | Description                                                                                                                                                                                              | Default value  |
|-----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------|
| center       | *REAL(2)*: The center of rotation. | Mandatory keyword. |
| radius       | *REAL*: Radius defining sliding mesh interface. | Mandatory keyword. |
| angle       | *REAL*: Angle defining sliding mesh interface. | Mandatory keyword. |
| rotation axis       | *CHARACTER*: Rotation axis. Possible values are `"x"`, `"y"`, or `"z"`. | Mandatory keyword. |


This algorithm can perform a p-adaptation to decrease the truncation error below a threshold.

```Markdown
#define slidingmesh
   centerx       = TE
   centery
   radius
   angle
   rotation axis
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