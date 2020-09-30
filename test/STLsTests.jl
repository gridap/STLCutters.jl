module STLsTests

using Test
using STLCutters
using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs

using STLCutters: read_stl
using STLCutters: compute_stl_model

using STLCutters: is_water_tight
using STLCutters: merge_nodes

X,T,N = read_stl(joinpath(@__DIR__,"data/cube.stl"))

stl = compute_stl_model(T,X)
@test !is_water_tight(stl)
writevtk(stl.grid,"cube",cellfields=["normals"=>N])


stl = merge_nodes(stl)
@test is_water_tight(stl)
writevtk(stl.grid,"cube",cellfields=["normals"=>N])

end # module
