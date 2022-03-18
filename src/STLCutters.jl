module STLCutters

import MeshIO
import Plots
import GraphRecipes

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
using ProgressMeter

import Gridap: writevtk
import Gridap.ReferenceFEs: get_polytope
import Gridap.ReferenceFEs: is_simplex
import Gridap.ReferenceFEs: is_n_cube
import Gridap.ReferenceFEs: get_vertex_coordinates
import Gridap.ReferenceFEs: get_faces
import Gridap.ReferenceFEs: get_face_vertices
import Gridap.ReferenceFEs: num_faces
import Gridap.ReferenceFEs: num_vertices
import Gridap.ReferenceFEs: num_edges
import Gridap.ReferenceFEs: num_facets
import Gridap.ReferenceFEs: num_dims
import Gridap.ReferenceFEs: num_point_dims
import Gridap.ReferenceFEs: get_offset
import Gridap.ReferenceFEs: get_offsets
import Gridap.ReferenceFEs: get_facedims
import Gridap.ReferenceFEs: simplexify
import Gridap.ReferenceFEs: get_bounding_box
import Gridap.Geometry: num_cells
import Gridap.Geometry: get_cell_vertices

import GridapEmbedded.Interfaces: cut
import GridapEmbedded.CSG: get_tree
import GridapEmbedded.CSG: compatible_geometries
import GridapEmbedded.CSG: similar_geometry

import Plots: plot

export STLGeometry
export STLCutter

export subtriangulate
export get_bounding_box

include("SimplexFaces.jl")
include("STLs.jl")
include("Polyhedra.jl")
include("SubTriangulations.jl")

include("Embedded.jl")

end # module
