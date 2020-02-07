module RawSTLsTests

using Test
using STLCutter
import STLCutter: compact_map, map_repeated_vertices

stl = RawSTL(joinpath(@__DIR__,"data/cube.stl"))

@test num_dims(stl) == 3

vertices_map = map_repeated_vertices(stl)
vertices_map = map_repeated_vertices(stl)
vertices_map = compact_map(vertices_map)

@test maximum(vertices_map) == 8

end # module
