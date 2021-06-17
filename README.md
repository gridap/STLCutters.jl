# STLCutters

STL to cell-wise triangulation to solve FE problems in [Gridap.jl](https://github.com/gridap/Gridap.jl) through [GridapEmbedded.jl](https://github.com/gridap/GridapEmbedded.jl)

<!--
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pmartorell.github.io/STLCutters.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pmartorell.github.io/STLCutters.jl/dev)
-->
[![CI](https://github.com/pmartorell/STLCutters.jl/workflows/CI/badge.svg)](https://github.com/pmartorell/STLCutters.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/pmartorell/STLCutters.jl/branch/master/graph/badge.svg?token=4mowFw2RKC)](https://codecov.io/gh/pmartorell/STLCutters.jl)


## Installation

```julia
# Type ] to enter package mode
pkg> add git@github.com:pmartorell/STLCutters.jl.git
```

## Examples

### Sub-triangulation examples

Use a test geometry, e.g., 47076.stl (Chichen Itza)
```julia
julia> include("examples/SubTriangulation.jl")
julia> filename = "test/data/47076.stl"
julia> SubTriangulation.main(filename,nmax=50,output="example1")
```
![Example 1](examples/example1.png)

Download a geometry directly from [Thingi10k](https://ten-thousand-models.appspot.com/), e.g, [37384](https://ten-thousand-models.appspot.com/detail.html?file_id=37384)
```julia
julia> include("examples/SubTriangulation.jl")
julia> filename = SubTriangulation.download(37384)
julia> SubTriangulation.main(filename,nmax=50,output="example2")
```
![Example 2](examples/example2.png)

### Poisson examples

Solve a Poisson equation on a test geometry, e.g., 293137.stl (Low-Poly Bunny)
 ```julia
julia> include("examples/Poisson.jl")
julia> filename = "test/data/293137.stl"
julia> Poisson.main(filename,n=20,output="example3")
```

![Example 3](examples/example3.png)
