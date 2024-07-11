module SequentialTests

using Test

@time @testset "CutterSeq" begin include("CutterTests.jl") end

@time @testset "PoissonSeq" begin include("PoissonTests.jl") end

end
