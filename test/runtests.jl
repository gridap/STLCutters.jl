module STLCuttersTests

using Test

@testset "CellMeshes" begin include("CellMeshesTests.jl") end

@testset "STLs" begin include("STLsTests.jl") end
@testset "Intersections" begin include("IntersectionsTests.jl") end
@testset "RefineVertices" begin include("RefineVerticesTests.jl") end
@testset "RefineEdges" begin include("RefineEdgesTests.jl") end
@testset "RefineVerticesEdges" begin include("RefineVerticesEdgesTests.jl") end
@testset "RefineFacets" begin include("RefineFacetsTests.jl") end
@testset "GridRefinement" begin include("GridRefinementTests.jl") end
@testset "SurfaceRefinement" begin include("SurfaceRefinementTests.jl") end

end # module
