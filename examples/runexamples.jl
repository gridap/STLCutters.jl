module RunExample

using Test

@time @testset "SubTriangulation" begin
  include("SubTriangulation.jl")

  filename = joinpath(@__DIR__,"..","test/data/47076.stl")
  SubTriangulation.main(filename,n=50)

  filename = SubTriangulation.download(37384)
  filename = joinpath(@__DIR__,"..","test/data/37384.stl")
  SubTriangulation.main(filename,n=50)
  # rm(filename)

end

@time @testset "Poisson" begin
  include("Poisson.jl")

  filename = joinpath(@__DIR__,"..","test/data/293137.stl")
  Poisson.main(filename,n=20)
end

@time @testset "LinearElasticity" begin
  include("LinearElasticity.jl")

  filename = joinpath(@__DIR__,"..","test/data/550964.stl")
  LinearElasticity.main(filename,n=50)
end

@time @testset "Stokes" begin
  include("Stokes.jl")

  filename = joinpath(@__DIR__,"..","test/data/47076.stl")
  Stokes.main(filename,n=10)
end

end # module
