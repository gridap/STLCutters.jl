module STLCutters

import MeshIO
import Downloads

using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
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
using GridapDistributed
using PartitionedArrays

using Gridap.Geometry: get_face_to_parent_face
using Gridap.Geometry: get_cell_to_parent_cell
using Gridap.ReferenceFEs: check_polytope_graph
using Gridap.ReferenceFEs: get_graph
using Gridap.ReferenceFEs: get_metadata
using Gridap.ReferenceFEs: isactive
using Gridap.ReferenceFEs: point_eltype
using GridapEmbedded.Distributed: consistent_bgcell_to_inoutcut!
using GridapEmbedded.Distributed: consistent_bgfacet_to_inoutcut!
using GridapEmbedded.Distributed: DistributedEmbeddedDiscretization
using GridapEmbedded.Distributed: AbstractEmbeddedDiscretization
using GridapDistributed: DistributedDiscreteModel
using GridapDistributed: remove_ghost_cells

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
import Gridap.ReferenceFEs: simplexify_interior
import Gridap.ReferenceFEs: simplexify_surface
import Gridap.ReferenceFEs: generate_polytope_data
import Gridap.ReferenceFEs: get_bounding_box
import Gridap.ReferenceFEs: get_facet_orientations
import Gridap.ReferenceFEs: get_facet_normal
import Gridap.ReferenceFEs: get_edge_tangent
import Gridap.ReferenceFEs: get_dimranges
import Gridap.ReferenceFEs: get_dimrange
import Gridap.ReferenceFEs: GraphPolytope
import Gridap.Geometry: num_cells
import Gridap.Geometry: get_cell_vertices
import Gridap.Geometry: get_polytopes
import Gridap.Geometry: Triangulation

import GridapEmbedded.Interfaces: cut
import GridapEmbedded.Interfaces: cut_facets
import GridapEmbedded.Interfaces: get_background_model
import GridapEmbedded.Interfaces: get_geometry
import GridapEmbedded.Interfaces: compute_bgcell_to_inoutcut
import GridapEmbedded.Interfaces: compute_bgfacet_to_inoutcut
import GridapEmbedded.Interfaces: compute_subcell_to_inout
import GridapEmbedded.Interfaces: EmbeddedBoundary
import GridapEmbedded.Interfaces: GhostSkeleton
import GridapEmbedded.Interfaces: cut_facets
import GridapEmbedded.Interfaces: cut_facets
import GridapEmbedded.CSG: get_tree
import GridapEmbedded.CSG: compatible_geometries
import GridapEmbedded.CSG: similar_geometry
import GridapEmbedded.AgFEM: aggregate
import GridapEmbedded.Distributed: change_bgmodel
import GridapEmbedded.Distributed: get_ls_to_bgcell_to_inoutcut
import GridapEmbedded.Distributed: get_ls_to_bgfacet_to_inoutcut
import GridapEmbedded.Distributed: remove_ghost_subfacets


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
include("Distributed.jl")

end # module
