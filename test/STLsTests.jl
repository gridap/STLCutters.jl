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

using STLCutters: STL

X,T,N = read_stl(joinpath(@__DIR__,"data/cube.stl"))

stlmodel = compute_stl_model(T,X)
@test !is_water_tight(stlmodel)
#writevtk(stl.grid,"cube",cellfields=["normals"=>N])

stl = STL(stlmodel)

@test num_faces(stl) == num_faces(stlmodel)
@test num_vertices(stl) == num_vertices(stlmodel)
@test num_edges(stl) == num_edges(stlmodel)
@test num_facets(stl) == num_facets(stlmodel)

stlmodel = merge_nodes(stlmodel)
@test is_water_tight(stlmodel)
#writevtk(stl.grid,"cube",cellfields=["normals"=>N])

stl = STL(stlmodel)

@test num_faces(stl) == num_faces(stlmodel)

f = save_as_stl(stlmodel,"_cube")
X1,T1,N1 = read_stl(f)
@test X1 == X
@test T1 == T
@test N1 == N
rm(f)

X,T,N = read_stl(joinpath(@__DIR__,"data/65904.stl"))
stlmodel = compute_stl_model(T,X)
stlmodel = merge_nodes(stlmodel)

stls = split_disconnected_parts(stlmodel)

@test length(stls) == 3
@test sum(num_cells,stls) == num_cells(stlmodel)
@test sum(num_vertices,stls) == num_vertices(stlmodel)

f = save_as_stl(stls,"_65904")

rm.(f)

end # module
