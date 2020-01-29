module STLCutterTests

using Test

@testset "Points" begin include("PointsTests.jl") end

@testset "TablesTests" begin include("TablesTests.jl") end

@testset "RawSTLs" begin include("RawSTLsTests.jl") end

@testset "ConformingSTLs" begin include("ConformingSTLsTests.jl") end

@testset "GeometriesTests" begin include("GeometriesTests.jl") end

end # module
