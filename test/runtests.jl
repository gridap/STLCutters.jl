module STLCutterTests

using Test

@testset "MutableVectorValues" begin include("MutableVectorValuesTests.jl") end

@testset "VectorValues" begin include("VectorValuesTests.jl") end

@testset "Points" begin include("PointsTests.jl") end

@testset "Segments" begin include("SegmentsTests.jl") end

@testset "Triangles" begin include("TrianglesTests.jl") end

@testset "Quadrilaters" begin include("QuadrilatersTests.jl") end

@testset "Hexahedrons" begin include("HexahedronsTests.jl") end

@testset "Tetrahedrons" begin include("TetrahedronsTests.jl") end

@testset "BoundingBoxes" begin include("BoundingBoxesTests.jl") end

@testset "Tables" begin include("TablesTests.jl") end

@testset "STLs" begin include("STLsTests.jl") end

@testset "SurfaceMeshes" begin include("SurfaceMeshesTests.jl") end

@testset "CartesianMeshes" begin include("CartesianMeshesTests.jl") end

@testset "VolumeMeshes" begin include("VolumeMeshesTests.jl") end

@testset "CellMeshes" begin include("CellMeshesTests.jl") end

@testset "IncrementalSurfaceMeshes" begin include("IncrementalSurfaceMeshesTests.jl") end

@testset "BulkMeshes" begin include("BulkMeshesTests.jl") end

end # module
