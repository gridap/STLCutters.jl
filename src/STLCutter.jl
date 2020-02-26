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
export STL
export SurfaceMesh
export CartesianMesh

export average
export num_vertices
export num_edges
export num_facets
export num_faces
export num_dims
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

macro abstractmethod()
  quote
    error("This function belongs to an interface definition and cannot be used.")
  end
end

"""
    tfill(v, ::Val{D}) where D
Returns a tuple of length `D` that contains `D` times the object `v`.
In contrast to `tuple(fill(v,D)...)` which returns the same result, this function is type-stable.
"""
function tfill(v, ::Val{D}) where D
  t = tfill(v, Val{D-1}())
  (v,t...)
end

tfill(v,::Val{0}) = ()
tfill(v,::Val{1}) = (v,)
tfill(v,::Val{2}) = (v,v)
tfill(v,::Val{3}) = (v,v,v)

include("MutableVectorValues.jl")

include("VectorValues.jl")

include("Points.jl")

include("Segments.jl")

include("Triangles.jl")

include("Tetrahedrons.jl")

include("BoundingBoxes.jl")

include("Tables.jl")

include("STLs.jl")

include("SurfaceMeshes.jl")

include("CartesianMeshes.jl")

include("BulkMeshes.jl")

end # module
