module STLsTests

using Test
using STLCutters
using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs

using STLCutters: read_stl
using STLCutters: compute_stl_model
using STLCutters: split_disconnected_parts
using STLCutters: save_as_stl

using STLCutters: is_water_tight
using STLCutters: merge_nodes

X,T,N = read_stl(joinpath(@__DIR__,"data/cube.stl"))

stl = compute_stl_model(T,X)
@test !is_water_tight(stl)
#writevtk(stl.grid,"cube",cellfields=["normals"=>N])


stl = merge_nodes(stl)
@test is_water_tight(stl)
#writevtk(stl.grid,"cube",cellfields=["normals"=>N])

f = save_as_stl(stl,"_cube")
X1,T1,N1 = read_stl(f)
@test X1 == X
@test T1 == T
@test N1 == N
rm(f)

X,T,N = read_stl(joinpath(@__DIR__,"data/65904.stl"))
stl = compute_stl_model(T,X)
stl = merge_nodes(stl)

stls = split_disconnected_parts(stl)

@test length(stls) == 3
@test sum(num_cells.(stls)) == num_cells(stl)
@test sum(num_vertices.(stls)) == num_vertices(stl)

f = save_as_stl(stl,"_65904")

rm.(f)

end # module
