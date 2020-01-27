module STLCutterTests

using Test

@testset "Points" begin include("PointsTests.jl") end

@testset "ConformingSTLs" begin include("ConformingSTLsTests.jl") end

end # module
