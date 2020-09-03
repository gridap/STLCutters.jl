module STLCutters

using LinearAlgebra
using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Helpers


export Table
export Point
export Segment

export distance
export projection
export have_intersection
export intersection_point
export insert_vertices!
export insert_edges!
export distribute_faces
export compute_grid
export initial_mesh

include("Intersections.jl")
include("CellRefinement.jl")
include("FaceInsertion.jl")
include("FaceDistribution.jl")

end # module
