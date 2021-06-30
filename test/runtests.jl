module STLCuttersTests

using Test

@testset "STLs" begin include("STLsTests.jl") end
@testset "SimplexFaces" begin include("SimplexFacesTests.jl") end
@testset "Polyhedra" begin include("PolyhedraTests.jl") end
@testset "SubTriangulations" begin include("SubTriangulationsTests.jl") end
@testset "Rotations" begin include("Rotations.jl") end
@testset "SampleGeometries" begin include("SampleGeometriesTests.jl") end

@testset "Integration" begin include("Integration.jl") end
@testset "Poisson" begin include("Poisson.jl") end
@testset "PoissonAgFEM" begin include("PoissonAgFEM.jl") end

include(joinpath(@__DIR__,"..","examples","runexamples.jl"))

end # module
