# Distributed

## Introduction

When dealing with large-scale problems, this package can be accelerated through two types of parallelization. The first one is multi-threading, which uses [Julia Threads](https://docs.julialang.org/en/v1/base/multi-threading/) for shared memory parallelization (e.g., `julia -t 4`). This method adds some speed-up. However, it is only efficient for a reduced number of threads.

The second one is a distributed memory computing. For such parallelization, we use MPI (`mpiexec -np 4 julia input.jl`) through [`PartitionedArrays`](https://www.francescverdugo.com/PartitionedArrays.jl/stable) and [`GridapDistributed`](https://gridap.github.io/GridapDistributed.jl/dev/). With MPI we can compute large-scale problems efficiently, up to thousands of cores.

Additionally, the distributed memory implementation is built on top of [GridapEmbedded](https://github.com/gridap/GridapEmbedded.jl). GridapEmbedded provides parallelization tools since [v0.9.2](https://github.com/gridap/GridapEmbedded.jl/releases/tag/v0.9.2).


## Underlying distributed algorithms

The implementation of multi-threading is straightforwardly applied to the embarrassingly parallel loops. However, the distributed memory implementation requires more involved algorithms for the global computations, e.g., the classification of the background cells as inside or outside. These algorithms are described in Chapter 4 of the following PhD thesis.

> Pere A. Martorell. "Unfitted finite element methods for explicit boundary representations". PhD Thesis. Universitat Politècnica de Catalunya. 2024. [hdl.handle.net/10803/690625](https://www.tdx.cat/handle/10803/690625)

## Multi-threading usage

The usage of this package with multi-threading is the same as for the serial case (see [Serial example](@ref)). The user only needs to set the number of threads when initializing Julia. E.g.,

```bash
julia -threads 4 
```

## Distributed memory usage

The distributed usage needs to be set up at the driver level. The user needs to install and load `MPI.jl`, `PartitionedArrays.jl`, and `GridapDistributed.jl`. 

Here, we provide a basic example of solving the Poisson equation with distributed STLCutters with 8 MPI tasks. Run the following command.

```bash
mpiexec -np 8 julia poisson.jl
```

!!! note
    Instead of running `mpiexec -np 8 julia poisson.jl`, it is recommended to use the `mpiexec` function from the `MPI.jl` package.
    ```julia
    using MPI
    mpiexec() do cmd
        run(`$cmd -np 8 $(Base.julia_cmd()) --project=$(Base.active_project()) poisson.jl`)
    end
    ```

Where `poisson.jl` is the following code.

```julia
using STLCutters
using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
parts = (2,2,2)
cells = (10,10,10)
filename = "stl_file_path.stl"
with_mpi() do distribute
  ranks = distribute(LinearIndices((prod(parts),)))
  # Domain and discretization
  geo = STLGeometry(filename)
  pmin,pmax = get_bounding_box(geo)
  model = CartesianDiscreteModel(ranks,parts,pmin,pmax,cells)
  cutgeo = cut(model,geo)
  # Cell aggregation
  model,cutgeo,aggregates = aggregate(AggregateAllCutCells(),cutgeo)
  # Triangulations
  Ω_act = Triangulation(cutgeo,ACTIVE)
  Ω = Triangulation(cutgeo)
  Γ = EmbeddedBoundary(cutgeo)
  nΓ = get_normal_vector(Γ)   
  dΩ = Measure(Ω,2)
  dΓ = Measure(Γ,2)
  # FE spaces
  Vstd = TestFESpace(Ω_act,ReferenceFE(lagrangian,Float64,1))
  V = AgFEMSpace(model,Vstd,aggregates)
  U = TrialFESpace(V)
  # Weak form
  γ = 10.0
  h = (pmax - pmin)[1] / cells[1]
  ud(x) = x[1] - x[2]
  f = 0
  a(u,v) =
    ∫( ∇(v)⋅∇(u) )dΩ +
    ∫( (γ/h)*v*u  - v*(nΓ⋅∇(u)) - (nΓ⋅∇(v))*u )dΓ
  l(v) =
    ∫( v*f )dΩ +
    ∫( (γ/h)*v*ud - (nΓ⋅∇(v))*ud )dΓ
  # Solve
  op = AffineFEOperator(a,l,U,V)
  uh = solve(op)
  writevtk(Ω,"results",cellfields=["uh"=>uh])
end
```
!!! note
    The STL file can be downloaded using [`download_thingi10k`](@ref).

!!! note
    It is recommended to use `with_debug()` instead of `with_mpi()` for debugging in serialized execution, see more details in [`PartitionedArrays`](https://www.francescverdugo.com/PartitionedArrays.jl/stable).


!!! note
    One can consider a different stabilization of the small cut-cell problem instead of AgFEM. Then, the `aggregate` and `AgFEMSpace` need to be removed.
    See more examples in [`GridapEmbedded`](https://github.com/gridap/GridapEmbedded.jl)


!!! warning
    Even though the distributed algorithms are proven to be efficient for large-scale weak scaling tests [Martorell, 2024](https://www.tdx.cat/handle/10803/690625). The performance of this implementation is not tested.


## Usage with `p4est`

The implementation is general for all the `Gridap.jl` triangulations.  Thus, we can use `p4est` through `GridapP4est.jl` in order to exploit its features. With `p4est`, we can take advantage of adaptive mesh refinement (AMR) and redistribution. 

### Redistribute with `p4est`

The following example is using `redistribute` to improve the load balancing of the subdomains.

```julia
using STLCutters
using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using GridapP4est
parts = 8
cells = (2,2,2)
filename = "stl_file_path.stl"
with_mpi() do distribute
  ranks = distribute(LinearIndices((prod(parts),)))
  # Domain and discretization
  geo = STLGeometry(filename)
  pmin,pmax = get_bounding_box(geo)
  coarse_model = CartesianDiscreteModel(pmin,pmax,cells)
  model = OctreeDistributedDiscreteModel(ranks,coarse_model,2)
  cutgeo = cut(model,geo)
  # Redistribute to avoid void subdmoains
  weights = compute_resdistribute_weights(cutgeo)
  model, = GridapDistributed.redistribute(model,weights=weights)
  # Re-compute discretization
  cutgeo = cut(model,geo)
  # Cell aggregation
  model,cutgeo,aggregates = aggregate(AggregateAllCutCells(),cutgeo)
  # Triangulations
  Ω_act = Triangulation(cutgeo,ACTIVE)
  Ω = Triangulation(cutgeo)
  Γ = EmbeddedBoundary(cutgeo)
  nΓ = get_normal_vector(Γ)   
  dΩ = Measure(Ω,2)
  dΓ = Measure(Γ,2)
  # FE spaces
  Vstd = TestFESpace(Ω_act,ReferenceFE(lagrangian,Float64,1))
  V = AgFEMSpace(model,Vstd,aggregates)
  U = TrialFESpace(V)
  # Weak form
  γ = 10.0
  h = (pmax - pmin)[1] / cells[1]
  ud(x) = x[1] - x[2]
  f = 0
  a(u,v) =
    ∫( ∇(v)⋅∇(u) )dΩ +
    ∫( (γ/h)*v*u  - v*(nΓ⋅∇(u)) - (nΓ⋅∇(v))*u )dΓ
  l(v) =
    ∫( v*f )dΩ +
    ∫( (γ/h)*v*ud - (nΓ⋅∇(v))*ud )dΓ
  # Solve
  op = AffineFEOperator(a,l,U,V)
  uh = solve(op)
  writevtk(Ω,"results",cellfields=["uh"=>uh])
end
```

!!! note
    For compatibility reasons, we are not interested on adding `GridapP4est` as a dependency. The interested can add it to its own driver. 

!!! note
    The `compute_resdistribute_weights` is a `GridapEmbedded` function since [#95](https://github.com/gridap/GridapEmbedded.jl/pull/95)

### Adaptive mesh refinement with `p4est`

The `p4est` library is designed to provide AMR. In the following example, we can `adapt` a mesh using the `STLCutters.jl` output and generate a new discretization on the refined mesh. Using `STLCutters.jl` works out-of-the-box with 2:1 balanced adapted meshes. The user can use this discretization to solve partial differential equations and combine them with `redistribute`.

```julia
using STLCutters
using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using GridapP4est
parts = 8
cells = (16,16,16)
filename = "stl_file_path.stl"
with_mpi() do distribute
  ranks = distribute(LinearIndices((prod(parts),)))
  # Domain and discretization
  geo = STLGeometry(filename)
  pmin,pmax = get_bounding_box(geo)
  coarse_model = CartesianDiscreteModel(pmin,pmax,cells)
  model = OctreeDistributedDiscreteModel(ranks,coarse_model,2)
  cutgeo = cut(model,geo)
  # Redistribute to avoid void subdmoains
  weights = compute_adaptive_flags(cutgeo)
  model, = Gridap.Adaptivity.adapt(model,weights=weights)
  model = get_model(model.dmodel)
  # Re-compute discretization
  cutgeo = cut(model,geo)
  Ωin = Triangulation(cutgeo,PHYSICAL_IN)
  Ωout = Triangulation(cutgeo,PHYSICAL_OUT)
  Ωbg = Triangulation(model)
  dΩin = Measure(Ωin,2)
  dΩout = Measure(Ωout,2)
  dΩbg = Measure(Ωbg,2)
  e = ∑( ∫(1)dΩin ) + ∑( ∫(1)dΩout ) - ∑( ∫(1)dΩbg )
  i_am_main(ranks) && println("Volume error: $e")
  writevtk(Ωin,"trian")
  writevtk(Ωbg,"trian_bg")
end
```

!!! note
    The `compute_adaptive_flags` is a `GridapEmbedded` function since [#95](https://github.com/gridap/GridapEmbedded.jl/pull/95)