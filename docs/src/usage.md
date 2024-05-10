# Usage

## Installation

STLCutters is a registered package. You can install it by running:

```julia
# Use ] to enter the Pkg REPL mode
pkg> add STLCutters
```

## Load package

Load the package normally with

```
using STLCutters
```

## Downloading STL files

The STLCutters package works with STL files. Therefore, first, you need to create or download an STL file. The [Thingi10k](https://ten-thousand-models.appspot.com/) dataset is one of the options for testing. You can download any indexed geometry with [`download_thingi10k`](@ref). E.g.,

```
filename = download_thingi10k(293137)
```

## Discretization

Since STLCutters is an extension of [GridapEmbedded](https://github/gridap/GridapEmbedded.jl) it utilizes the same workflow to solve PDEs on embedded domains. In particular, STLCutters extends the [`cut`](@ref STLCutters.cut) function from GridapEmbedded.

We load the STL file with an [`STLGeometry`](@ref) object. E.g.,

```
filename = download_thingi10k(293137)
geo = STLGeometry(filename)
```

Then, we define a `DiscreteModel` around, e.g., a `CartesianDiscreteModel`.

```julia
pmin,pmax = get_bounding_box(geo)
model = CartesianDiscreteModel(pmin,pmax,())
```

Now, we can [`cut`](@ref STLCutters.cut) the `model` with the `STLGeometry`.

```julia
cutgeo = cut(model,geo)
```

!!! note
    Some STL files may not be properly defined. It is recommended to check the requisites before proceeding with [`check_requisites`](@ref). In particular, one checks if the STL is a closed surface, a manifold and has a bounded density of faces.

## Usage with Gridap

Once, the geometry is discretized one can generate the embedded triangulations to solve PDEs with [Gridap](https://github.com/gridap/Gridap.jl), see also [Gridap Tutorials](https://gridap.github.io/Tutorials/stable).

Like in GridapEmbedded, we extract the embedded triangulations as follows.

```julia
Ωact = Triangulation(cutgeo,ACTIVE)
Ω = Triangulation(cutgeo)
Γ = EmbeddedBoundary(cutgeo)
Λ = SkeletonTriangulation(cutgeo)
```

## Serial example

Now, we provide an example of the solution of a Poisson problem on the embedded domain.

```julia
using STLCutters
using GridapEmbedded
cells = (10,10,10)
filename = "stl_file_path.stl"
# Domain and discretization
geo = STLGeometry(filename)
pmin,pmax = get_bounding_box(geo)
model = CartesianDiscreteModel(pmin,pmax,cells)
cutgeo = cut(model,geo)
# Cell aggregation
aggregates = aggregate(AggregateAllCutCells(),cutgeo)
# Triangulations
Ω_act = Triangulation(cutgeo,ACTIVE)
Ω = Triangulation(cutgeo)
Γ = EmbeddedBoundary(cutgeo)
nΓ = get_normal_vector(Γ)   
dΩ = Measure(Ω,2)
dΓ = Measure(Γ,2)
# FE spaces
Vstd = TestFESpace(Ω_act,ReferenceFE(lagrangian,Float64,1))
V = AgFEMSpace(Vstd)
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
```

!!! note
    The STL file can be downloaded using [`download_thingi10k`](@ref).

!!! note
    One can consider a different stabilization of the small cut-cell problem instead of AgFEM. Then, the `aggregate` and `AgFEMSpace` need to be removed.
    See more examples in [`GridapEmbedded`](https://github.com/gridap/GridapEmbedded.jl)
