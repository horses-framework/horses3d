# Changelog

This file documents notable changes in the repository.

## Unreleased

#### Added
- Added integration with codecov to obtain coverage reports for each PR. Coverage
  is collected from the gfortran `DEBUG` CI runs.
- Added the `EXTRA_FFLAGS` option to the solver and problem file Makefiles to pass
  additional compilation flags.

#### Changed
- Removed temporary array creation when matching boundary tags in the GMSH mesh reader.
