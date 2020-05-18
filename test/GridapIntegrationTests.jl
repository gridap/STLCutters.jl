module testing

using Test

using STLCutters.GridapIntegration

using Gridap
using GridapEmbedded
using Gridap.Geometry

import Gridap: ∇

# Manufactured solution
u(x) = x[1] + x[2]# - x[3]
∇u(x) = VectorValue( 1, 1) #, -1)
Δu(x) = 0
f(x) = - Δu(x)
ud(x) = u(x)
∇(::typeof(u)) = ∇u



n = 2
outputfile = "out"

partition = (n,n)
geo = square(name="square")
#partition = (n,n,n)
#geo = STLGeometry("test/data/cube.stl",name="cube")
box = 1.5*get_metadata(geo)
bgmodel = CartesianDiscreteModel(box.pmin,box.pmax,partition)




# Forcing data
#ud = 1
#f = 10

# Cut the background model
cutgeo = cut(bgmodel,geo)

# Setup integration meshes
trian_Ω = Triangulation(cutgeo)
trian_Γd = EmbeddedBoundary(cutgeo)
trian_Γg = GhostSkeleton(cutgeo,"square")

# Setup normal vectors
n_Γd = get_normal_vector(trian_Γd)
n_Γg = get_normal_vector(trian_Γg)

writevtk(trian_Ω,"trian_O")
writevtk(trian_Γd,"trian_Gd",cellfields=["normal"=>n_Γd])
writevtk(trian_Γg,"trian_Gg",cellfields=["normal"=>n_Γg])
#writevtk(Triangulation(bgmodel),"bgtrian")

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
a_Ω(u,v) = ∇(v)*∇(u)
l_Ω(v) = v*f
a_Γd(u,v) = (γd/h)*v*u  - v*(n_Γd*∇(u)) - (n_Γd*∇(v))*u
l_Γd(v) = (γd/h)*v*ud - (n_Γd*∇(v))*ud
a_Γg(v,u) = (γg*h)*jump(n_Γg*∇(v))*jump(n_Γg*∇(u))

# FE problem
t_Ω = AffineFETerm(a_Ω,l_Ω,trian_Ω,quad_Ω)
t_Γd = AffineFETerm(a_Γd,l_Γd,trian_Γd,quad_Γd)
t_Γg = LinearFETerm(a_Γg,trian_Γg,quad_Γg)
op = AffineFEOperator(U,V,t_Ω,t_Γd,t_Γg)
uh = solve(op)

#writevtk(trian_Ω,"trian_Omega",celldata=["bgcell"=>collect(Int64,get_cell_id(trian_Ω))])
#writevtk(trian_Γd,"trian_Gamma_d",celldata=["bgcell"=>collect(Int64,get_cell_id(trian_Γd))])
#writevtk(bgmodel,"bgmodel")

# Postprocess
uh_Ω = restrict(uh,trian_Ω)
if outputfile !== nothing
  writevtk(trian_Ω,outputfile,cellfields=["uh"=>uh_Ω])
end

uh_Γd = restrict(uh,trian_Γd)
writevtk(trian_Γd,"u_trian_Gd",cellfields=["uh"=>uh_Γd])


tol = 1.0e-9
e = u - uh_Ω
@show el2 = sqrt(sum(integrate(e*e,trian_Ω,quad_Ω)))

uh_Γ = restrict(uh,trian_Γd)
e = u - uh_Γ
@show el2_Γ = sqrt(sum(integrate(e*e,trian_Γd,quad_Γd)))

uh = interpolate(U,u)

uh_Ω = restrict(uh,trian_Ω)
e = u - uh_Ω
@show el2 = sqrt(sum(integrate(e*e,trian_Ω,quad_Ω)))


uh_Γ = restrict(uh,trian_Γd)
e = u - uh_Γ
@show el2_Γ = sqrt(sum(integrate(e*e,trian_Γd,quad_Γd)))


#using Gridap.Integration
#q_0 = get_coordinates(quad_Ω)
#q_to_x = get_cell_map(trian_Ω)
#x_0 = evaluate(q_to_x,q_0)
#writevtk(x_0,"x_0")
#writevtk(q_0,"q_0")


trian_in = trian_Ω
quad_in = quad_Ω
quad_Γ = quad_Γd
trian_Γ = trian_Γd
trian_Γg_in = trian_Γg
model_in = model
n_Γ = n_Γd

V_in = TestFESpace(model=model_in,valuetype=Float64,reffe=:Lagrangian,order=1,conformity=:H1)

v_in = FEFunction(V_in,rand(num_free_dofs(V_in)))

v_in_Γ = restrict(v_in,trian_Γ)

v_in_in = restrict(v_in,trian_in)

v_in_Γg_in = restrict(v_in,trian_Γg_in)

# Check divergence theorem
u(x) = x[1] + x[2] #+ x[3]
u_in = interpolate(V_in,u)
u_in_Γ = restrict(u_in,trian_Γ)
u_in_in = restrict(u_in,trian_in)
a = sum( integrate(∇(v_in_in)*∇(u_in_in),trian_in,quad_in) )
b = sum( integrate(v_in_Γ*n_Γ*∇(u_in_Γ),trian_Γ,quad_Γ) )

@show abs(a-b)
#@test abs(a-b) < 1.0e-9

val = collect(integrate(∇(v_in_in)*∇(u_in_in),trian_in,quad_in))
val2 = collect(integrate(∇(u_in_in),trian_in,quad_in))
val2 = collect(integrate(u_in_in,trian_in,quad_in))
val3 = collect(integrate(u_in_Γ,trian_Γ,quad_Γ))

writevtk(trian_in,"grad_v_in_in",order=2,cellfields=["∇v_in_in"=>∇(v_in_in),"∇u_in_in"=>∇(u_in_in),"u_in_in"=>u_in_in])
writevtk(trian_in,"vals",order=2,cellfields=["v"=>val,"v2"=>val2])
writevtk(trian_Γ,"vals_g",order=2,cellfields=["v"=>val3])
writevtk(trian_Γ,joinpath("trian_G"),order=2,cellfields=["u_in"=>u_in_Γ,"normal"=>n_Γ])
#writevtk(trian_Γg_in,joinpath("trian_Gg_in"),cellfields=["v_in"=>mean(v_in_Γg_in)])

v_in_Γg_in = restrict(v_in,trian_Γg)


#writevtk(Triangulation(model),joinpath("trian"),order=2,cellfields=["v_in"=>v_in])
#writevtk(trian_Γ,joinpath("trian_G"),order=2,cellfields=["v_in"=>v_in_Γ,"normal"=>n_Γ])
writevtk(trian_Γg,joinpath("trian_Gg_in"),cellfields=["v_in"=>mean(v_in_Γg_in)])


_u = 1
U = TrialFESpace(V_in)
u_in = interpolate(U,_u)
u_in_Γ = restrict(u_in,trian_Γ)


u_in_in = restrict(u_in,trian_in)

@show sum( integrate(u_in_in,trian_in,quad_in) ) 
@show sum( integrate(u_in_Γ,trian_Γ,quad_Γ) ) , surface(geo) 


end # module
