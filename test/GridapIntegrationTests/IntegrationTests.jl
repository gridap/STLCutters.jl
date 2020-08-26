module IntegrationTests

using STLCutters.GridapIntegration

using Gridap
import Gridap: ∇
using GridapEmbedded
using Gridap.Geometry
using Test

n = 10

geos = [ square(name="wall"), STLGeometry(joinpath(@__DIR__,"../data/cube.stl"),name="wall") ]
partitions = [ (n,n), (n,n,n) ]
volumes = [ 1.0, 1.0 ]
surfaces = [ 4.0, 6.0 ]


for i in 1:length(geos)
  partition = partitions[i]
  geo = geos[i]

  box = 1.5*get_metadata(geo)

  bgmodel = CartesianDiscreteModel(box.pmin,box.pmax,partition)

  # Cut the background model
  cutgeom = cut(bgmodel,geo)

  model_in = DiscreteModel(cutgeom)

  degree = 2

  trian_in = Triangulation(cutgeom)
  test_triangulation(trian_in)
  quad_in = CellQuadrature(trian_in,degree)
  vol = sum(integrate(1,trian_in,quad_in))
  @test abs(volumes[i] - vol) < 1.0e-3

  trian_Γ = EmbeddedBoundary(cutgeom)
  test_triangulation(trian_Γ)
  quad_Γ = CellQuadrature(trian_Γ,degree)
  surf = sum(integrate(1,trian_Γ,quad_Γ))
  @test abs(surf - surfaces[i]) < 1.0e-3

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

end

end # module
