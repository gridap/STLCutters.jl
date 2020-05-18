module GridapIntegrationTests

using Test

@testset "Integration" begin include("IntegrationTests.jl") end

@testset "Poisson" begin include("PoissonTests.jl") end

@testset "Stokes" begin include("StokesTests.jl") end

end

