title: Input and Output Files
---

[TOC]


## Input files

- Control file (\*.control)
- Mesh file (\*.mesh / \*.h5 / \*.msh)
- Polynomial order file (\*.omesh)
- Problem File (ProblemFile.f90)

<small>Notes on the GMSH format (\*.msh) and general workflow using GMSH.</small>

- Curved geometry supported up to polynomial order 5.
- Curved geometry should be generated using the following options: `tools -> options -> mesh -> general -> element order`.
- HORSES3D can read mesh format 4.1 and 2.2 (legacy format).
- The solution to most of the problems mesh reading is to load it in GMSH and export to format 2.2 to have a clean ASCII file.

## Output files

- Solution file (\*.hsol)
- Horses mesh file (\*.hmesh)
- Boundary information (\*.bmesh)
- Partition file (\*.pmesh)
- Polynomial order file (\*.omesh)
- MLRK level file (\*.lmesh), written by 'Multi Level RK3' runs (see [Explicit solvers](explicit-solvers.html))
- Monitor files (\*.volume / \*.surface / \*.residuals)

@Note Refrain from using tabs 
