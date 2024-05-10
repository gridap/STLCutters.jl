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


!!! warning
    Do not mix `@docs` here
