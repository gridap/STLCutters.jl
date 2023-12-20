module LinearElasticity

using IterativeSolvers: cg
using AlgebraicMultigrid: smoothed_aggregation
using AlgebraicMultigrid: aspreconditioner

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

  @test check_requisites(geo,bgmodel)

  # Cut the background model

  cutgeo = cut(bgmodel,geo)

  strategy = AggregateAllCutCells()
  aggregates = aggregate(strategy,cutgeo)

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

  x = cg(A,b,verbose=true,Pl=p,reltol=1e-10)

  uh = FEFunction(U,x)

  σh = project(norm∘σ∘ε(uh),Ω_act,dΩ,Float64,order,aggregates)

  if output !== nothing
    colors = color_aggregates(aggregates,bgmodel)
    writevtk(Ω_bg,output*"_trian",
      celldata=["aggregate"=>aggregates,"color"=>colors],
      cellfields=["uh"=>uh])
    writevtk(Ω,output*"_results",cellfields=["uh"=>uh,"σ"=>σh,"ε"=>ε(uh)])
    writevtk(Γ,output*"_boundary",cellfields=["uh"=>uh,"σ"=>σh,"ε"=>ε(uh)])
    writevtk(geo,output*"_geo")
  end

end

function project(q,model,dΩ,T,order,aggregates)
  reffe = ReferenceFE(lagrangian,T,order)
  Vstd = FESpace(model,reffe,conformity=:H1)
  V = AgFEMSpace(Vstd,aggregates)
  a(u,v) = ∫( u⊙v )*dΩ
  l(v) = ∫( v⊙q )*dΩ
  op = AffineFEOperator(a,l,V,V)
  A = get_matrix(op)
  b = get_vector(op)
  p = aspreconditioner(smoothed_aggregation(A))
  x = cg(A,b,verbose=false,Pl=p,reltol=1e-8)
  qh = FEFunction(V,x)
  qh
end

end # module
