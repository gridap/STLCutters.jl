module LinearElasticity

using IterativeSolvers: cg
using AlgebraicMultigrid

using STLCutters
using Gridap
using GridapEmbedded
using Test

using Gridap.TensorValues

function main(filename;n=20,δ=0.2,force=(0,0,-1)::NTuple{3},output=nothing)

  # Constitutive law
  E = 1e5
  ν = 0.3
  μ = E / (2*(1 + ν))
  λ = (E*ν)/((1+ν)*(1-2*ν))
  σ(ε) = λ*tr(ε)*one(ε) + 2*μ*ε

  f = VectorValue(force)
  ud = VectorValue(0,0,0)

  geo = STLGeometry( filename )

  pmin,pmax = get_bounding_box(geo)
  zmin = pmin[3]+0.001*norm(pmin-pmax)
  diagonal = pmax-pmin
  pmin = pmin - diagonal*δ
  pmax = pmax + diagonal*δ
  pmin = VectorValue(pmin[1],pmin[2],zmin)
  partition = (n,n,n)

  bgmodel = CartesianDiscreteModel(pmin,pmax,partition)

  # Cut the background model

  cutgeo,facet_to_inoutcut = cut(bgmodel,geo)

  strategy = AggregateAllCutCells()
  aggregates = aggregate(strategy,cutgeo,facet_to_inoutcut)

  # Setup integration meshes
  Ω_bg = Triangulation(bgmodel)
  Ω = Triangulation(cutgeo,PHYSICAL)
  Γ = EmbeddedBoundary(cutgeo)

  # Setup normal vectors
  n_Γ = get_normal_vector(Γ)

  # Setup Lebesgue measures
  order = 1
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓ = Measure(Γ,degree)

  # Setup FESpace
  Ω_act = Triangulation(cutgeo,ACTIVE)
  reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
  Vstd = FESpace(Ω_act,reffe,dirichlet_tags=["boundary"],conformity=:H1)

  V = AgFEMSpace(Vstd,aggregates)
  U = TrialFESpace(V,[ud])

  # Weak form

  a(u,v) = ∫( ε(v) ⊙ (σ∘ε(u)) )dΩ

  l(v) = ∫( v⋅f ) * dΩ

  # FE problem
  op = AffineFEOperator(a,l,U,V)

  A = get_matrix(op)
  b = get_vector(op)

  p = aspreconditioner(smoothed_aggregation(A))

  x = cg(A,b,verbose=true,Pl=p,reltol=1e-12)

  uh = FEFunction(U,x)

  if output !== nothing
    colors = color_aggregates(aggregates,bgmodel)
    writevtk(Ω_bg,output*"_trian",
      celldata=["aggregate"=>aggregates,"color"=>colors],
      cellfields=["uh"=>uh])
    writevtk(Ω,output*"_results",cellfields=["uh"=>uh,"σ"=>σ∘ε(uh),"ε"=>ε(uh)])
    writevtk(geo,output*"_geo")
  end

end

end # module
