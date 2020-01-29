module RawSTLsTests

using Test
using STLCutter

stl = RawSTL(joinpath(@__DIR__,"data/cube.stl"))

@test num_dims(stl) == 3

end # module
