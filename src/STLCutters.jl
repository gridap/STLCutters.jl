module STLCutters

using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Helpers
using GridapEmbedded
using GridapEmbedded.CSG
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters
using LinearAlgebra
using FileIO

import MeshIO

import Gridap.ReferenceFEs: get_polytope
import Gridap.ReferenceFEs: is_simplex 
import Gridap.ReferenceFEs: is_n_cube 
import Gridap.ReferenceFEs: get_vertex_coordinates
import Gridap.ReferenceFEs: num_faces 
import Gridap.ReferenceFEs: num_vertices
import Gridap.ReferenceFEs: num_edges 
import Gridap.ReferenceFEs: num_facets 
import Gridap.ReferenceFEs: num_dims 

import GridapEmbedded.Interfaces: cut
import GridapEmbedded.CSG: get_tree
import GridapEmbedded.CSG: compatible_geometries
import GridapEmbedded.CSG: similar_geometry

export Table
export Point
export Segment
export Triangle
export STLGeometry

export distance
export projection
export have_intersection
export intersection_point
export insert_vertices!
export insert_edges!
export insert_facets!
export distribute_faces
export compute_grid
export initial_mesh

export cut
export surface
export square

include("CellMeshes.jl")

include("Intersections.jl")
include("STLs.jl")
include("CellRefinement.jl")
include("FaceInsertion.jl")
include("FaceDistribution.jl")
include("GridRefinement.jl")
include("SurfaceRefinement.jl")

include("Embedded.jl")

include("tables/LookupTables.jl")

end # module
