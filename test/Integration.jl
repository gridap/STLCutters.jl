module IntegrationTests

using STLCutters
using Gridap
import Gridap: ∇
using Gridap.Geometry
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

bgmodel = CartesianDiscreteModel(pmin,pmax,partition)

# Cut the background model

cutgeom = cut(bgmodel,geo)
model_in = DiscreteModel(cutgeom)

degree = 2

trian_in = Triangulation(cutgeom)
test_triangulation(trian_in)
quad_in = CellQuadrature(trian_in,degree)
vol = sum(integrate(1,trian_in,quad_in))
@test vol ≈ (0.85-0.15)^2

trian_Γ = EmbeddedBoundary(cutgeom)
test_triangulation(trian_Γ)
quad_Γ = CellQuadrature(trian_Γ,degree)
surf = sum(integrate(1,trian_Γ,quad_Γ))
@test surf ≈ (0.85-0.15)*2^2

trian_Γg_in = GhostSkeleton(cutgeom)

n_Γ = get_normal_vector(trian_Γ)


V_in = TestFESpace(model=model_in,valuetype=Float64,reffe=:Lagrangian,order=1,conformity=:H1)

v_in = FEFunction(V_in,rand(num_free_dofs(V_in)))

v_in_Γ = restrict(v_in,trian_Γ)

v_in_in = restrict(v_in,trian_in)

v_in_Γg_in = restrict(v_in,trian_Γg_in)

trian = Triangulation(bgmodel)

# Check divergence theorem
u(x) = x[1] + x[2]
u_in = interpolate(V_in,u)
u_in_Γ = restrict(u_in,trian_Γ)
u_in_in = restrict(u_in,trian_in)
a = sum( integrate(∇(v_in_in)⋅∇(u_in_in),trian_in,quad_in) )
b = sum( integrate(v_in_Γ⋅n_Γ⋅∇(u_in_Γ),trian_Γ,quad_Γ) )
@test abs(a-b) < 1e-9

end # module
