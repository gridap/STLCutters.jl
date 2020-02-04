module STLCutter

using FileIO
using LinearAlgebra
import LinearAlgebra: dot, norm

export num_components
export component_type
export mutable
export VectorValue

export Point
export Segment
export Triangle
export Tetrahedron
export HexaCell
export RawSTL
export ConformingSTL

export average
export num_vertices
export num_facets
export num_dims
export norm
export distance
export normal
export center
export volume
export have_intersection
export intersection
export table_cache
export getlist
export getlist!
export TableOfVectors
export AbstractTable

include("MutableVectorValues.jl")

include("VectorValues.jl")

include("Points.jl")

include("Segments.jl")

include("Triangles.jl")

include("Tetrahedrons.jl")

include("HexaCells.jl")

include("Tables.jl")

include("RawSTLs.jl")

include("ConformingSTLs.jl")

end # module
