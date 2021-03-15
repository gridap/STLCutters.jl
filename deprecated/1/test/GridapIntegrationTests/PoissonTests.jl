module PoissonTests

using STLCutters.GridapIntegration

using Gridap
import Gridap: ∇
using GridapEmbedded
using Test

# Manufactured solution
u(x) = x[1] + x[2] # - x[3]
∇u(x) = VectorValue( 1, 1) #, -1)
Δu(x) = 0
f(x) = - Δu(x)
ud(x) = u(x)
∇(::typeof(u)) = ∇u

n = 10
partition = (n,n)
geo = square(name="wall")

#partition = (n,n,n)
#geo = STLGeometry("test/data/Bunny-LowPoly.stl",name="wall")
#test/data/Bunny-LowPoly.stl

box = 1.5*get_metadata(geo)

bgmodel = CartesianDiscreteModel(box.pmin,box.pmax,partition)

# Cut the background model
cutgeo = cut(bgmodel,geo)

# Setup integration meshes
trian_Ω = Triangulation(cutgeo)
trian_Γd = EmbeddedBoundary(cutgeo)
trian_Γg = GhostSkeleton(cutgeo,"wall")

# Setup normal vectors
n_Γd = get_normal_vector(trian_Γd)
n_Γg = get_normal_vector(trian_Γg)

# Setup cuadratures
order = 1
quad_Ω = CellQuadrature(trian_Ω,2*order)
quad_Γd = CellQuadrature(trian_Γd,2*order)
quad_Γg = CellQuadrature(trian_Γg,2*order)

# Setup FESpace
model = DiscreteModel(cutgeo)
V = TestFESpace(model=model,valuetype=Float64,reffe=:Lagrangian,order=order,conformity=:H1)
U = TrialFESpace(V)

# Weak form
γd = 10.0
γg = 0.1
h = (box.pmax - box.pmin)[1] / partition[1]
a_Ω(u,v) = ∇(v)⋅∇(u)
l_Ω(v) = v⋅f
a_Γd(u,v) = (γd/h)⋅v⋅u  - v⋅(n_Γd⋅∇(u)) - (n_Γd⋅∇(v))⋅u
l_Γd(v) = (γd/h)⋅v⋅ud - (n_Γd⋅∇(v))⋅ud
a_Γg(v,u) = (γg⋅h)⋅jump(n_Γg⋅∇(v))⋅jump(n_Γg⋅∇(u))

# FE problem
t_Ω = AffineFETerm(a_Ω,l_Ω,trian_Ω,quad_Ω)
t_Γd = AffineFETerm(a_Γd,l_Γd,trian_Γd,quad_Γd)
t_Γg = LinearFETerm(a_Γg,trian_Γg,quad_Γg)
op = AffineFEOperator(U,V,t_Ω,t_Γd,t_Γg)
uh = solve(op)

# Postprocess
uh_Ω = restrict(uh,trian_Ω)
writevtk(trian_Ω,"results",cellfields=["uh"=>uh_Ω])

tol = 1e-9
e = u - uh_Ω
el2 = sqrt(sum(integrate(e⋅e,trian_Ω,quad_Ω)))
@test el2 < tol
eh1 = sqrt(sum(integrate(e⋅e+a_Ω(e,e),trian_Ω,quad_Ω)))
@test eh1 < tol

end # module
