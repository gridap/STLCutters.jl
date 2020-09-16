module STLCutters

using LinearAlgebra
using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Helpers

import Gridap.ReferenceFEs: get_polytope

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

include("CellMeshes.jl")

include("Intersections.jl")
include("STLs.jl")
include("CellRefinement.jl")
include("FaceInsertion.jl")
include("FaceDistribution.jl")

include("tables/LookupTables.jl")

end # module
