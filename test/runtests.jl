module STLCuttersTests

using Test

@testset "STLs" begin include("STLsTests.jl") end
@testset "Intersections" begin include("IntersectionsTests.jl") end
@testset "Polyhedra" begin include("PolyhedraTests.jl") end
@testset "Rotations" begin include("Rotations.jl") end

@testset "Poisson" begin include("Poisson.jl") end

end # module
