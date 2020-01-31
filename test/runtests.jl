module STLCutterTests

using Test

@testset "MutableVectorValues" begin include("MutableVectorValuesTests.jl") end

@testset "VectorValues" begin include("VectorValuesTests.jl") end

@testset "Points" begin include("PointsTests.jl") end

@testset "SegmentsTests" begin include("SegmentsTests.jl") end

@testset "TrianglesTests" begin include("TrianglesTests.jl") end

@testset "HexaCellsTests" begin include("HexaCellsTests.jl") end

@testset "TablesTests" begin include("TablesTests.jl") end

@testset "RawSTLs" begin include("RawSTLsTests.jl") end

@testset "ConformingSTLs" begin include("ConformingSTLsTests.jl") end

#@testset "GeometriesTests" begin include("GeometriesTests.jl") end

end # module
