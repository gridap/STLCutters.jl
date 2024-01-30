#module IntegrationTests

using STLCutters
using Gridap
import Gridap: ∇
using GridapEmbedded
using Gridap.ReferenceFEs
using Test

using STLCutters: compute_stl_model
using STLCutters: read_stl, merge_nodes, get_bounding_box

# Manufactured solution
u0(x) = x[1] + x[2] - x[3]

X,T,N = read_stl(joinpath(@__DIR__,"data/cube.stl"))
stl = compute_stl_model(X,T)
stl = merge_nodes(stl)
# writevtk(stl,"geo")
n = 10
δ = 0.2
pmin,pmax = get_bounding_box(stl)
diagonal = pmax-pmin
pmin = pmin - diagonal*δ
pmax = pmax + diagonal*δ
partition = (n,n,n)

geo = STLGeometry(stl)
bgmodel = CartesianDiscreteModel(pmin,pmax,partition)

# Cut the background model

cutgeo = cut(bgmodel,geo)

# Setup integration meshes
Ω = Triangulation(cutgeo,PHYSICAL)
Γd = EmbeddedBoundary(cutgeo)
Γg = GhostSkeleton(cutgeo)

# Setup normal vectors
n_Γd = get_normal_vector(Γd)
n_Γg = get_normal_vector(Γg)

#writevtk(Ω,"trian_O")
#writevtk(Γd,"trian_Gd",cellfields=["normal"=>n_Γd])
#writevtk(Γg,"trian_Gg",cellfields=["normal"=>n_Γg])
#writevtk(Triangulation(bgmodel),"bgtrian")

# Setup Lebesgue measures
order = 1
degree = 2*order
dΩ = Measure(Ω,degree)
dΓd = Measure(Γd,degree)
dΓg = Measure(Γg,degree)

vol = sum( ∫(1)*dΩ  )
surf = sum( ∫(1)*dΓd )

@test vol ≈ 1
@test surf ≈ 6

# Setup FESpace
Ω_act = Triangulation(cutgeo,ACTIVE,geo)
V = TestFESpace(Ω_act,ReferenceFE(lagrangian,Float64,order),conformity=:H1)
v = FEFunction(V,rand(num_free_dofs(V)))
u = interpolate(u0,V)

# Check divergence theorem
a = sum( ∫( ∇(v)⋅∇(u) ) * dΩ )
b = sum( ∫( v⋅n_Γd⋅∇(u) ) * dΓd )
@test abs( a-b ) < 1e-9

# Moment fitted
Ω_act_in = Triangulation(cutgeo,ACTIVE_IN,geo)
Ω_act_out = Triangulation(cutgeo,ACTIVE_OUT,geo)
dΩᵐ_in = Measure(Ω_act_in,Quadrature(momentfitted,cutgeo,degree,in_or_out=IN))
dΩᵐ_out = Measure(Ω_act_out,Quadrature(momentfitted,cutgeo,degree,in_or_out=OUT))

f = x -> x[1] + 1
f = 1
# @test ∑(∫(f)dΩᵐ_in) ≈ ∑(∫(f)dΩ)

# Simplex background
#

bgmodel = CartesianDiscreteModel(pmin,pmax,partition)
bgmodel = simplexify(bgmodel,positive=true)

# Cut the background model

cutgeo = cut(bgmodel,geo)

# Setup integration meshes
Ω = Triangulation(cutgeo,PHYSICAL)
Γd = EmbeddedBoundary(cutgeo)
Γg = GhostSkeleton(cutgeo)

# Setup normal vectors
n_Γd = get_normal_vector(Γd)
n_Γg = get_normal_vector(Γg)

#writevtk(Ω,"trian_O")
#writevtk(Γd,"trian_Gd",cellfields=["normal"=>n_Γd])
#writevtk(Γg,"trian_Gg",cellfields=["normal"=>n_Γg])
#writevtk(Triangulation(bgmodel),"bgtrian")

# Setup Lebesgue measures
order = 1
degree = 2*order
dΩ = Measure(Ω,degree)
dΓd = Measure(Γd,degree)
dΓg = Measure(Γg,degree)

vol = sum( ∫(1)*dΩ  )
surf = sum( ∫(1)*dΓd )

@test vol ≈ 1
@test surf ≈ 6

# Setup FESpace
Ω_act = Triangulation(cutgeo,ACTIVE,geo)
V = TestFESpace(Ω_act,ReferenceFE(lagrangian,Float64,order),conformity=:H1)
v = FEFunction(V,rand(num_free_dofs(V)))
u = interpolate(u0,V)

# Check divergence theorem
a = sum( ∫( ∇(v)⋅∇(u) ) * dΩ )
b = sum( ∫( v⋅n_Γd⋅∇(u) ) * dΓd )
@test abs( a-b ) < 1e-9

# Moment fitted

#end
