module SequentialTests

using Test

@time @testset "CutterSeq" begin include("CutterTests.jl") end

end
