module STLCutters

import MeshIO
import Plots
import GraphRecipes
import Downloads

using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Helpers
using Gridap.ReferenceFEs
using GridapEmbedded
using GridapEmbedded.CSG
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters
using LinearAlgebra
using FillArrays
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
import Gridap.Geometry: get_polytopes

import GridapEmbedded.Interfaces: cut
import GridapEmbedded.Interfaces: cut_facets
import GridapEmbedded.CSG: get_tree
import GridapEmbedded.CSG: compatible_geometries
import GridapEmbedded.CSG: similar_geometry

import Plots: plot

import Base: split

export STLGeometry
export STLCutter

export subtriangulate
export get_bounding_box
export check_requisites
export download_thingi10k

include("SimplexFaces.jl")
include("STLs.jl")
include("Polyhedra.jl")
include("SubTriangulations.jl")

include("Embedded.jl")

end # module
