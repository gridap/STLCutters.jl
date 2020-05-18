module GridapIntegrationTests

using Test

@testset "Integration" begin include("IntegrationTests.jl") end

@testset "Poisson" begin include("PoissonTests.jl") end

end

