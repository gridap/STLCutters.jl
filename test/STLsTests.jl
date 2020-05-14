module STLsTests

using Test
using STLCutters
import STLCutters: delete_repeated_vertices!, delete_repeated_vertices

stl = STL(joinpath(@__DIR__,"data/cube.stl"))
@test num_dims(stl) == 3

stl_1 = delete_repeated_vertices(stl)
@test stl_1 != stl
@test num_vertices(stl_1) == 8
@test num_vertices(stl) == 36
@test maximum(stl.facet_to_vertices) == 36
@test maximum(stl_1.facet_to_vertices) == 8


stl_2 = delete_repeated_vertices!(stl)
@test stl_2 === stl
@test num_vertices(stl) == 8
@test stl_1 == stl_2
@test stl_1 !== stl_2

end # module
