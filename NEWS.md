# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- `Polyhedron` is now defined in `Gridap` as `GeneralPolytope{3}<Polytope{3}`. Its internal `data` is renamed to `metadata`.
- `get_data(::Polyhedron)` is renamed to `get_metadata(::Polyhedron)`.
- `simplexify(::Polyhedron)` now returns the same data structure as `simplexify(::Polytope)`. The previous implementation is renamed to `simplexify_interior`.
- `check_graph` is renamed to `check_polytope_graph`.


## [0.2.1] - 2024-05-02

### Added

- Distributed cutter (`cut()`) using `GridapDistributed` and `PartitionedArrays` since [#28](https://github.com/gridap/STLCutters.jl/pull/28).


## [0.2.0] - 2024-02-15

### Added

- `cut_facets()` for the integration of cut cell boundaries with `GridapEmbeded` [#27](https://github.com/gridap/STLCutters.jl/pull/27). `cut_facets()` allows moment fitting integration.

### Changed

- `cut()` interface has been changed to align with `GridapEmbeded`.
- `subtriangulate()` also returns cut facet data.
- The output order of _simplexify_ functions has changed: (1st) coordinates, (2nd) connectivities.
  
### Removed

- `plot(::Polyhedron)` is removed to reduce dependencies [#28](https://github.com/gridap/STLCutters.jl/pull/28).


## [0.1.0] - 2021-09-04

A changelog is not maintained for this version.
