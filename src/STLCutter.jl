module STLCutter

using FileIO
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
export RawSTL
export ConformingSTL
export StructuredBulkMesh

export average
export num_vertices
export num_edges
export num_facets
export num_dfaces
export num_dims
export norm
export distance
export normal
export center
export area
export volume
export measure
export measure_sign
export have_intersection
export intersection
export projection
export closest_point
export table_cache
export getlist
export getlist!
export TableOfVectors
export CompressedTable
export AbstractTable
export writevtk
export num_cells
export get_cell

macro check(test)
  quote
    @assert $(esc(test)) $(string(test))
  end
end

macro check(test,msg)
  quote
    @assert $(esc(test)) $msg
  end
end

include("MutableVectorValues.jl")

include("VectorValues.jl")

include("Points.jl")

include("Segments.jl")

include("Triangles.jl")

include("Tetrahedrons.jl")

include("BoundingBoxes.jl")

include("Tables.jl")

include("RawSTLs.jl")

include("ConformingSTLs.jl")

include("BulkMeshes.jl")

end # module
