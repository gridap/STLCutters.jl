module STLsTests

using Test
using STLCutter
import STLCutter: delete_repeated_vertices!

stl = STL(joinpath(@__DIR__,"data/cube.stl"))
@show stl
@test num_dims(stl) == 3

stl = delete_repeated_vertices!(stl)
@show stl
@test num_vertices(stl) == 8

end # module
