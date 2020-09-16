module STLCuttersTests

using Test

@testset "CellMeshes" begin include("CellMeshesTests.jl") end

@testset "Intersections" begin include("IntersectionsTests.jl") end
@testset "RefineVertices" begin include("RefineVerticesTests.jl") end
@testset "RefineEdges" begin include("RefineEdgesTests.jl") end

end # module
