module PoisonTests

using STLCutters
using Gridap
import Gridap: ∇
using GridapEmbedded
using Test

using STLCutters: compute_stl_model
using STLCutters: read_stl, merge_nodes, get_bounding_box

# Manufactured solution
u(x) = x[1] + x[2] - x[3]
∇u(x) = VectorValue( 1, 1, -1)
Δu(x) = 0
f(x) = - Δu(x)
ud(x) = u(x)
∇(::typeof(u)) = ∇u


@time X,T,N = read_stl(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))
#X,T,N = read_stl(joinpath(@__DIR__,"data/cube.stl"))
#@time X,T,N = read_stl(joinpath(@__DIR__,"data/441708_sf.obj"))
@time stl = compute_stl_model(T,X)
@time stl = merge_nodes(stl)
# writevtk(stl,"geo")
n = 20
δ = 0.2
pmin,pmax = get_bounding_box(stl)
diagonal = pmax-pmin
pmin = pmin - diagonal*δ
pmax = pmax + diagonal*δ
partition = (n,n,n)

geo = STLGeometry(stl)
bgmodel = CartesianDiscreteModel(pmin,pmax,partition)

# Cut the background model

@time cutgeo, = cut(bgmodel,geo)

# Setup integration meshes
Ω = Triangulation(cutgeo)
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

vol = sum( ∫(1)*dΩ  )
surf = sum( ∫(1)*dΓd )
@show vol,surf

# Setup FESpace
model = DiscreteModel(cutgeo)
V = TestFESpace(model,ReferenceFE(lagrangian,Float64,order),conformity=:H1)
U = TrialFESpace(V) 

# Weak form
γd = 10.0
γg = 0.1
h = (pmax - pmin)[1] / partition[1]

a(u,v) =
  ∫( ∇(v)⋅∇(u) ) * dΩ +
  ∫( (γd/h)*v*u  - v*(n_Γd⋅∇(u)) - (n_Γd⋅∇(v))*u ) * dΓd 

l(v) =
  ∫( v*f ) * dΩ +
  ∫( (γd/h)*v*ud - (n_Γd⋅∇(v))*ud ) * dΓd

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

#writevtk(Ω,"results",cellfields=["uh"=>uh])

@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

end # module
