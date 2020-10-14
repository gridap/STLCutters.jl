module PoisonTests

using STLCutters
using Gridap
import Gridap: ∇
using GridapEmbedded
using Test

using STLCutters: compute_stl_model

# Manufactured solution
u(x) = x[1] + x[2] # - x[3]
∇u(x) = VectorValue( 1, 1) #, -1)
Δu(x) = 0
f(x) = - Δu(x)
ud(x) = u(x)
∇(::typeof(u)) = ∇u

n = 10
partition = (n,n)

pmin,pmax = Point(0.0,0.0),Point(1.0,1.0)

vertices = [ Point(0.15,0.15),
             Point(0.15,0.85),
             Point(0.85,0.85), 
             Point(0.85,0.15) ]
faces = [[1,2],[2,3],[3,4],[4,1]]

stl = compute_stl_model(Table(faces),vertices)

geo = STLGeometry(stl,name="wall")

#partition = (n,n,n)
#geo = STLGeometry("test/data/Bunny-LowPoly.stl",name="wall")
#test/data/Bunny-LowPoly.stl


bgmodel = CartesianDiscreteModel(pmin,pmax,partition)

# Cut the background model

cutgeo = cut(bgmodel,geo)

# Setup integration meshes
trian_Ω = Triangulation(cutgeo)
trian_Γ = EmbeddedBoundary(cutgeo)
trian_Γg = GhostSkeleton(cutgeo,"wall")

# Setup normal vectors
n_Γ = get_normal_vector(trian_Γ)
n_Γg = get_normal_vector(trian_Γg)

# Setup cuadratures
order = 1
quad_Ω = CellQuadrature(trian_Ω,2*order)
quad_Γ = CellQuadrature(trian_Γ,2*order)
quad_Γg = CellQuadrature(trian_Γg,2*order)

# Setup FESpace
model = DiscreteModel(cutgeo)
V = TestFESpace(model=model,valuetype=Float64,reffe=:Lagrangian,order=order,conformity=:H1)
U = TrialFESpace(V)

# Weak form
const γd = 10.0
const γg = 0.1
h = (pmax - pmin)[1] / partition[1]
a_Ω(u,v) = ∇(v)⋅∇(u)
l_Ω(v) = v*f
a_Γ(u,v) = (γd/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u
l_Γ(v) = (γd/h)*v*ud - (n_Γ⋅∇(v))*ud
a_Γg(v,u) = (γg*h)*jump(n_Γg⋅∇(v))*jump(n_Γg⋅∇(u))

# FE problem
t_Ω = AffineFETerm(a_Ω,l_Ω,trian_Ω,quad_Ω)
t_Γ = AffineFETerm(a_Γ,l_Γ,trian_Γ,quad_Γ)
t_Γg = LinearFETerm(a_Γg,trian_Γg,quad_Γg)
op = AffineFEOperator(U,V,t_Ω,t_Γ,t_Γg)
uh = solve(op)

# Postprocess
uh_Ω = restrict(uh,trian_Ω)
writevtk(trian_Ω,"results",cellfields=["uh"=>uh_Ω])

tol = 1e-9
e = u - uh_Ω
el2 = sqrt(sum(integrate(e*e,trian_Ω,quad_Ω)))
@test el2 < tol
eh1 = sqrt(sum(integrate(e*e+a_Ω(e,e),trian_Ω,quad_Ω)))
@test eh1 < tol

end # module
