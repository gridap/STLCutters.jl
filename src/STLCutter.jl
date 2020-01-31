module STLCutter

using FileIO
using LinearAlgebra
import LinearAlgebra: dot, norm

export num_components
export component_type
export mutable
export VectorValue

export Point
export RawSTL
export ConformingSTL
export num_vertices
export num_facets
export num_dims
export norm
export distance
export table_cache
export getlist
export getlist!
export TableOfVectors
export AbstractTable

include("MutableVectorValues.jl")

include("VectorValues.jl")

include("Points.jl")

include("Tables.jl")

include("RawSTLs.jl")

include("ConformingSTLs.jl")

end # module
