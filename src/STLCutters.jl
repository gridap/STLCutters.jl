module STLCutters

using LinearAlgebra
using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Helpers


export Point
export Table

export insert_vertices!
export distribute_vertices
export compute_grid
export initial_mesh

include("CellRefinement.jl")
include("FaceInsertion.jl")
include("FaceDistribution.jl")

end # module
