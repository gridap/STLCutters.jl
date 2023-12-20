module Stokes

using Test
using STLCutters
using Gridap
using GridapEmbedded
using Gridap.ReferenceFEs

const ENABLE_MKL = haskey(ENV,"ENABLE_MKL")

if ENABLE_MKL
  import Pkg
  Pkg.add("GridapPardiso")
  Pkg.add("SparseMatricesCSR")
  using GridapPardiso
  using SparseMatricesCSR
end

function main(filename;n,output=nothing)

  geo = STLGeometry( filename )

  δ = 0.2
  pmin,pmax = get_bounding_box(geo)
  L = pmax-pmin
  zmin = pmin[3]+0.001*norm(L)
  pmin = pmin - L*δ - (1-δ)*L.*Point(0,1,0)
  pmax = pmax + L*δ + (2-δ)*L.*Point(0,1,0)
  pmin = VectorValue(pmin[1],pmin[2],zmin)

  # Background model
  cells = (n,3*n,n)
  bgmodel = CartesianDiscreteModel(pmin,pmax,cells)

  @test check_requisites(geo,bgmodel)

  labels = get_face_labeling(bgmodel)
  wall_tags = union(get_faces(HEX)[ [21,22,25,26] ]...)
  add_tag_from_tags!(labels,"wall",wall_tags)
  add_tag_from_tags!(labels,"inlet",23)

  # Cut the background model
  cutgeo = cut(bgmodel,geo)

  Ωb = Triangulation(bgmodel)
  Ω = Triangulation(cutgeo,PHYSICAL_OUT)
  Ωa = Triangulation(cutgeo,ACTIVE_OUT)
  Γ = EmbeddedBoundary(cutgeo)
  n_Γ = -get_normal_vector(Γ)

  k = 2
  reffe_u = ReferenceFE(lagrangian,VectorValue{3,Float64},k)
  reffe_p = ReferenceFE(lagrangian,Float64,k-1,space=:P)
  Vstd = TestFESpace(Ωa,reffe_u,dirichlet_tags=["inlet","wall"])
  Qstd = FESpace(Ωa,reffe_p,conformity=:L2)

  threshold = 0.5
  strategy = AggregateCutCellsByThreshold(threshold)
  aggregates = aggregate(strategy,cutgeo,geo,OUT)
  V = AgFEMSpace(Vstd,aggregates)
  Q = AgFEMSpace(Qstd,aggregates)

  function u_1d(x,xmin,xmax)
    L = xmax - xmin
    y = (x-xmin)/L
    4*y-4*y^2
  end
  vmax = 0.2*VectorValue(0,1,0)
  u_in(x) = vmax*u_1d(x[1],pmin[1],pmax[1])*u_1d(x[3],pmin[3],pmax[3])
  u_wall = VectorValue(0.,0.,0.)

  U = TrialFESpace(V,[u_in,u_wall])
  P = Q

  X = MultiFieldFESpace([U,P])
  Y = MultiFieldFESpace([V,Q])

  # perhaps for just a picture 2*k is enough also in cut cells
  dΩ = Measure(Ω,2*k) #2*(k-1),3*2*(k-1)
  dΓ = Measure(Γ,2*k) #2*3*k)

  β = 10
  h =  minimum(Tuple(pmax-pmin)./cells)
  τ = β*k^2/h

  a((u,p),(v,q)) =
    ∫( ∇(u)⊙∇(v) - p*(∇⋅v) - (∇⋅u)*q )dΩ +
    ∫(
      τ*v⋅u  -
      v⋅(n_Γ⋅∇(u)) -
      (n_Γ⋅∇(v))⋅u +
      (p*n_Γ)⋅v +
      (q*n_Γ)⋅u )dΓ

  l((v,q)) = 0

  if ENABLE_MKL
    assem = SparseMatrixAssembler(SparseMatrixCSR{1,Float64,Int},Vector{Float64},X,Y)
    op = AffineFEOperator(a,l,X,Y,assem)
    ls = PardisoSolver()
    solver = LinearFESolver(ls)
    uh,ph = solve(solver,op)
  else
    op = AffineFEOperator(a,l,X,Y)
    uh,ph = solve(op)
  end

  ph = project(ph,Ωa,dΩ,Float64,k,aggregates)

  if !isnothing(output)
    writevtk(Ωb,output*"_Ωb")
    writevtk(Ω,output*"_Ω",order=2,cellfields=["uh"=>uh,"ph"=>ph])
    writevtk(Ωa,output*"_Ωa",cellfields=["uh"=>uh,"ph"=>ph])
    writevtk(Γ,output*"_Γ",cellfields=["n"=>n_Γ,"uh"=>uh,"ph"=>ph])
  end

end

function project(q,model,dΩ,T,order,aggregates)
  reffe = ReferenceFE(lagrangian,T,order)
  Vstd = FESpace(model,reffe,conformity=:H1)
  V = AgFEMSpace(Vstd,aggregates)
  a(u,v) = ∫( u⊙v )*dΩ
  l(v) = ∫( v⊙q )*dΩ
  op = AffineFEOperator(a,l,V,V)
  qh = solve(op)
  qh
end

end # module
