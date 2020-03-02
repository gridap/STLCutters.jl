module STLCutter

using FileIO
using MeshIO
using LinearAlgebra
using WriteVTK
import LinearAlgebra: dot, norm


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
export Hexahedron
export STL
export SurfaceMesh
export CartesianMesh
export FaceCutter
export CellSubMesh

export average
export num_vertices
export num_edges
export num_facets
export num_faces

export norm
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
export get_dface_to_nfaces
export refine!
export add_vertex!

include("tables/LookupCutTables.jl");

include("Helpers.jl")

include("MutableVectorValues.jl")

include("VectorValues.jl")

include("Points.jl")

include("Segments.jl")

include("Triangles.jl")

include("Tetrahedrons.jl")

include("BoundingBoxes.jl")

include("Hexahedrons.jl")

include("Tables.jl")

include("STLs.jl")

include("SurfaceMeshes.jl")

include("CartesianMeshes.jl")

include("CellSubMeshes.jl")

include("BulkMeshes.jl")

end # module
