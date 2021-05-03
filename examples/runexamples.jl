module RunExample

using Test

@time @testset "SubTriangulation" begin
  include("SubTriangulation.jl")

  filename = joinpath(@__DIR__,"../test/data/47076.stl")
  SubTriangulation.main(filename,nmax=50)
  SubTriangulation.main(filename,nmax=10,nmin=10,kdtree=true)

  filename = joinpath(@__DIR__,"../test/data/47076_sf.obj")
  SubTriangulation.main(filename,nmax=50)

  filename = SubTriangulation.download(293137)
  SubTriangulation.main(filename,nmax=50)
  rm(filename)

  filename = SubTriangulation.download(80084)
  SubTriangulation.main(filename,nmax=20)
  rm(filename)

  filename = SubTriangulation.download(65395)
  SubTriangulation.main(filename,nmax=20)
  rm(filename)

  filename = SubTriangulation.download(77343)
  SubTriangulation.main(filename,nmax=20)
  rm(filename)

  filename = SubTriangulation.download(95436)
  SubTriangulation.main(filename,nmax=20)
  rm(filename)

  filename = SubTriangulation.download(243015)
  SubTriangulation.main(filename,nmax=100)
  rm(filename)

  filename = SubTriangulation.download(57657)
  SubTriangulation.main(filename,nmax=100,tolfactor=10000)
  rm(filename)

  filename = SubTriangulation.download(1452677)
  SubTriangulation.main(filename,nmax=20)
  rm(filename)

end

end # module
