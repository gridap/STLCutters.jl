module STLsTests

using Test
using STLCutters
using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs

using STLCutters: read_stl
using STLCutters: compute_stl_grid

using STLCutters: is_water_tight
using STLCutters: merge_nodes

X,T,N = read_stl("test/data/cube.stl")

stl = compute_stl_grid(T,X)
topo = GridTopology(stl)
@test !is_water_tight(topo)
#writevtk(stl,"cube",cellfields=["normals"=>N])


stl = merge_nodes(stl)
topo = GridTopology(stl)
@test is_water_tight(topo)
#writevtk(stl,"cube",cellfields=["normals"=>N])

end # module
