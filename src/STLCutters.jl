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
using GraphRecipes
using Plots

import MeshIO

import Gridap: writevtk
import Gridap.ReferenceFEs: get_polytope
import Gridap.ReferenceFEs: is_simplex 
import Gridap.ReferenceFEs: is_n_cube 
import Gridap.ReferenceFEs: get_vertex_coordinates
import Gridap.ReferenceFEs: num_faces 
import Gridap.ReferenceFEs: num_vertices
import Gridap.ReferenceFEs: num_edges 
import Gridap.ReferenceFEs: num_facets 
import Gridap.ReferenceFEs: num_dims 
import Gridap.ReferenceFEs: simplexify 
import Gridap.ReferenceFEs: get_bounding_box 

import GridapEmbedded.Interfaces: cut
import GridapEmbedded.CSG: get_tree
import GridapEmbedded.CSG: compatible_geometries
import GridapEmbedded.CSG: similar_geometry

import Plots: plot

export STLGeometry

export read_stl
export compute_stl_model
export merge_nodes
export merge_and_collapse
export get_bounding_box
export compute_grid
export is_water_tight
export compute_submesh
export volume, volumes
export surface, surfaces
export min_height


include("Intersections.jl")
include("STLs.jl")
include("Polyhedron.jl")

include("Embedded.jl")

include("Tests/Tests.jl")

end # module
