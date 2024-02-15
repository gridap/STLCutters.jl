module PoisonAgFEMTests

using IterativeSolvers: cg
using AlgebraicMultigrid

using STLCutters
using Gridap
using GridapEmbedded
using Test

using STLCutters: compute_stl_model
using STLCutters: read_stl, merge_nodes, get_bounding_box

# Manufactured solution
u(x) = x[1] + x[2] - x[3]
f(x) = - Δ(u)(x)
ud(x) = u(x)

stlpath = joinpath(@__DIR__,"data/Bunny-LowPoly.stl")

@time geo = STLGeometry( stlpath )

n = 28
δ = 0.2
pmin,pmax = get_bounding_box(geo)
diagonal = pmax-pmin
pmin = pmin - diagonal*δ
pmax = pmax + diagonal*δ
partition = (n,n,n)

bgmodel = CartesianDiscreteModel(pmin,pmax,partition)

# Cut the background model
cutter = STLCutter()
@time cutgeo = cut(cutter,bgmodel,geo)

strategy = AggregateAllCutCells()
aggregates = aggregate(strategy,cutgeo)

# Setup integration meshes
Ω_bg = Triangulation(bgmodel)
Ω = Triangulation(cutgeo,PHYSICAL)
Γd = EmbeddedBoundary(cutgeo)

# Setup normal vectors
n_Γd = get_normal_vector(Γd)

#writevtk(Ω,"trian_O")
#writevtk(Γd,"trian_Gd",cellfields=["normal"=>n_Γd])
#writevtk(Triangulation(bgmodel),"bgtrian")

# Setup Lebesgue measures
order = 1
degree = 2*order
dΩ = Measure(Ω,degree)
dΓd = Measure(Γd,degree)
#dΓg = Measure(Γg,degree)

vol = sum( ∫(1)*dΩ  )
surf = sum( ∫(1)*dΓd )
@show vol,surf

# Setup FESpace
Ω_act = Triangulation(cutgeo,ACTIVE)
Vstd = FESpace(Ω_act,ReferenceFE(lagrangian,Float64,order),conformity=:H1)

@time V = AgFEMSpace(Vstd,aggregates) #,Vser)
U = TrialFESpace(V)
# Weak form
γd = 10.0
h = (pmax - pmin)[1] / partition[1]

a(u,v) =
  ∫( ∇(v)⋅∇(u) ) * dΩ  +
  ∫( (γd/h)*v*u  - v*(n_Γd⋅∇(u)) - (n_Γd⋅∇(v))*u ) * dΓd

l(v) =
  ∫( v*f ) * dΩ +
  ∫( (γd/h)*v*ud - (n_Γd⋅∇(v))*ud ) * dΓd

# FE problem
@time op = AffineFEOperator(a,l,U,V)

A = get_matrix(op)
b = get_vector(op)

@time p = aspreconditioner(smoothed_aggregation(A))

@time  x = cg(A,b,verbose=true,Pl=p,reltol=1e-10)

@time uh = FEFunction(U,x)

@time e = u - uh

# Postprocess
l2(u) = sqrt(sum( ∫( u*u )*dΩ ))
h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ ))

el2 = l2(e)
eh1 = h1(e)
ul2 = l2(uh)
uh1 = h1(uh)

#colors = color_aggregates(aggregates,bgmodel)
#writevtk(Ω_bg,"trian",celldata=["aggregate"=>aggregates,"color"=>colors],cellfields=["uh"=>uh])
#writevtk(Ω,"results",cellfields=["uh"=>uh])

@test el2/ul2 < 1.e-9
@test eh1/uh1 < 1.e-9

end # module
