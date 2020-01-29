module STLCutter

using FileIO

export Point
export RawSTL
export ConformingSTL
export num_vertices
export num_facets
export num_dims
export distance
export table_cache
export getlist
export getlist!
export TableOfVectors
export AbstractTable

include("Points.jl")

include("Tables.jl")

include("RawSTLs.jl")

include("ConformingSTLs.jl")

end # module
