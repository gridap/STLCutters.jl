module STLCuttersTests

using Test

@testset "STLs" begin include("STLsTests.jl") end
@testset "Intersections" begin include("IntersectionsTests.jl") end
#@testset "Polyhedra" begin include("PolyhedraTests.jl") end

end # module
