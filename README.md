# STLCutters

STL to cell-wise triangulation to solve FE problems in [Gridap.jl](https://github.com/gridap/Gridap.jl) through [GridapEmbedded.jl](https://github.com/gridap/GridapEmbedded.jl)

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.jcp.2022.111162-blue)](https://authors.elsevier.com/c/1er0r508HsZ58)<!-- Restore on May 21th, 2022(https://doi.org/10.1016/j.jcp.2022.111162)-->
[![CI](https://github.com/gridap/STLCutters.jl/workflows/CI/badge.svg)](https://github.com/gridap/STLCutters.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/gridap/STLCutters.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/gridap/STLCutters.jl)


## Installation

```julia
# Type ] to enter package mode
pkg> add STLCutters
```

## Examples

### Sub-triangulation examples

Use a test geometry, e.g., `47076.stl` (Chichen Itza)
```julia
julia> include("examples/SubTriangulation.jl")
julia> filename = "test/data/47076.stl"
julia> SubTriangulation.main(filename,n=50,output="example1")
```
![Example 1](examples/example1.png)

Download a geometry directly from [Thingi10k](https://ten-thousand-models.appspot.com/), e.g, [37384](https://ten-thousand-models.appspot.com/detail.html?file_id=37384)
```julia
julia> include("examples/SubTriangulation.jl")
julia> filename = SubTriangulation.download(37384)
julia> SubTriangulation.main(filename,n=50,output="example2")
```
![Example 2](examples/example2.png)

### Finite Elements examples

Solve a **Poisson** equation on a test geometry, e.g., `293137.stl` (Low-Poly Bunny)
 ```julia
julia> include("examples/Poisson.jl")
julia> filename = "test/data/293137.stl"
julia> Poisson.main(filename,n=20,output="example3")
```

![Example 3](examples/example3.png)

Solve a **Linear Elasticity** problem on a test geometry, e.g., `550964.stl` (Eiffel Tower in a 5 degree slope)
 ```julia
julia> include("examples/LinearElasticity.jl")
julia> filename = "test/data/550964.stl"
julia> LinearElasticity.main(filename,n=50,force=(tand(5),0,-1),output="example4")
```

![Example 4](examples/example4.png)

Solve an **Incompressible Flow** problem on a test geometry, e.g., `47076.stl` (Chichen Itza)
 ```julia
julia> # ENV["ENABLE_MKL"] = "" ## Uncomment if GridapPardiso.jl requirements are fulfilled
julia> include("examples/Stokes.jl")
julia> filename = "test/data/47076.stl"
julia> Stokes.main(filename,n=10,output="example5")
```

![Example 5](examples/example5.png)
