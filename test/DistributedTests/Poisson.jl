module DistributedPoissonTests

using Gridap
using STLCutters
using GridapDistributed
using PartitionedArrays
using GridapEmbedded
using Test

function main(distribute;
  np = (1,1,1),
  nc = (2,2,2),
  geoname = "cube",
  δ = 0.2,
  tolfactor=10^4,
  simplex=false,
  vtk = false,
  verbose = false)

  # Manufactured solution
  u(x) = x[1] + x[2] - x[3]
  f(x) = - Δ(u)(x)
  ud(x) = u(x)

  filename = joinpath(@__DIR__,"..","data","$geoname.stl")
  if !isfile(filename)
    filename = download_thingi10k(id;path="")
  end

  ranks = distribute(LinearIndices((prod(np),)))
  geo = STLGeometry(filename)

  pmin,pmax = get_bounding_box(geo)
  diagonal = pmax-pmin
  pmin = pmin - diagonal*δ
  pmax = pmax + diagonal*δ
  bgmodel = CartesianDiscreteModel(ranks,np,pmin,pmax,nc)
  if simplex
    bgmodel = simplexify(bgmodel,positive=true)
  end

  cutter = STLCutter(;tolfactor)
  cutgeo = cut(cutter,bgmodel,geo)

  model,cutgeo,aggregates = aggregate(AggregateAllCutCells(),cutgeo)

  Ω = Triangulation(cutgeo,PHYSICAL_IN)
  Γ = EmbeddedBoundary(cutgeo)
  n_Γ = get_normal_vector(Γ)
  degree = 2
  dΩ = Measure(Ω,degree)
  dΓ = Measure(Γ,degree)

  # Setup FESpace
  order = 1
  Ω_act = Triangulation(cutgeo,ACTIVE)
  Vstd = TestFESpace(Ω_act,ReferenceFE(lagrangian,Float64,order),conformity=:H1)
  V = AggFESpace(model,Vstd,aggregates)
  U = TrialFESpace(V)

  # Weak form
  γ = 10.0
  h = (pmax - pmin)[1] / nc[1]

  a(u,v) =
    ∫( ∇(v)⋅∇(u) ) * dΩ +
    ∫( (γ/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u ) * dΓ

  l(v) =
    ∫( v*f ) * dΩ +
    ∫( (γ/h)*v*ud - (n_Γ⋅∇(v))*ud ) * dΓ

  # FE problem
  @time op = AffineFEOperator(a,l,U,V)
  @time uh = solve(op)

  e = u - uh

  # Postprocess
  l2(u) = sqrt(sum( ∫( u*u )*dΩ ))
  h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ ))

  el2 = l2(e)
  eh1 = h1(e)
  ul2 = l2(uh)
  uh1 = h1(uh)

  if vtk
    writevtk(Ω,"results",cellfields=["uh"=>uh])
  end

  if verbose
    println("L2 error: ", el2/ul2)
    println("H1 error: ", eh1/uh1)
  end

  @test el2/ul2 < 1.e-8
  @test eh1/uh1 < 1.e-7
end

end # module
