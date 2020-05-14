module STLCutters

using FileIO
using MeshIO
using LinearAlgebra
using WriteVTK
import LinearAlgebra: dot, norm, det, cross


export num_components
export component_type
export mutable
export VectorValue
export get_data

export Point
export Segment
export Triangle
export Tetrahedron
export BoundingBox
export Quadrilater
export Hexahedron
export STL
export SurfaceMesh
export CartesianMesh
export VolumeMesh
export MeshCutter
export CellMesh

export average
export num_vertices
export num_edges
export num_facets
export num_faces
export num_dfaces

export num_dims
export norm
export orthogonal
export distance
export normal
export center
export area
export volume
export measure
export measure_sign
export contains_projection
export have_intersection_point
export have_intersection
export intersection
export projection
export closest_point
export Table
export isactive
export remove!
export compact!
export writevtk
export num_cells
export get_cell
export get_dface_to_vertices
export get_faces
export refine!
export add_vertex!
export is_watter_tight
export compute_cell_to_surface_mesh_faces

include("tables/LookupCutTables.jl");

include("Helpers.jl")

include("MutableVectorValues.jl")

include("VectorValues.jl")

include("Points.jl")

include("Segments.jl")

include("Triangles.jl")

include("Tetrahedrons.jl")

include("BoundingBoxes.jl")

include("Quadrilaters.jl")

include("Hexahedrons.jl")

include("Tables.jl")

include("STLs.jl")

include("SurfaceMeshes.jl")

include("CartesianMeshes.jl")

include("VolumeMeshes.jl")

include("IncrementalSurfaceMeshes.jl")

include("CellMeshes.jl")

include("BulkMeshes.jl")


include("GridapIntegration/GridapIntegration.jl")

end # module
