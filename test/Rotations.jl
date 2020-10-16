module GridRefinementTests
  
using Test
using Gridap
using STLCutters
using DelimitedFiles

using Gridap.Geometry

using STLCutters: compute_stl_model
using STLCutters: refine_grid 
using STLCutters: volume,volumes
using STLCutters: FACE_CUT, FACE_IN, FACE_OUT
using STLCutters: get_bounding_box 
using STLCutters: read_stl 
using STLCutters: merge_nodes 



R(θ) = TensorValue(cos(θ),sin(θ),-sin(θ),cos(θ))
origin = Point(0.5,0.5)

bg_mesh = compute_linear_grid(QUAD4)
bg_mesh = CartesianGrid( Point(0.0,0.0),(0.2,0.2),(5,5))

vertices = [
  Point(0.2,0.2),
  Point(0.2,0.8),
  Point(0.8,0.8),
  Point(0.8,0.2) ]
faces = [[1,2],[2,3],[3,4],[4,1]]

θs = π.*[ exp10.(-10:-1); 1 ./ (5:-1:2) ]

for (i,θ) in enumerate(θs)

_vertices = origin .+ R(θ).⋅(vertices .- origin)
stl = compute_stl_model(Table(faces),_vertices)
out = refine_grid(bg_mesh,stl)
T,X,reffes,cell_types,cell_to_io,cell_to_bgcell,bgcell_to_ioc = out 
submesh = UnstructuredGrid(X,Table(T),reffes,cell_types)
vol = volume(submesh) + sum(volumes(bg_mesh).*(bgcell_to_ioc.≠FACE_CUT))
@test volume(bg_mesh) ≈ vol
vol1 = sum(volumes(submesh).*(cell_to_io.==FACE_IN))
vol2 = sum(volumes(bg_mesh).*(bgcell_to_ioc.==FACE_IN))
@test vol1+vol2 ≈ 0.6^2
writevtk(submesh,"submesh_$i",cellfields=["io"=>cell_to_io,"bgcell"=>cell_to_bgcell])
writevtk(bg_mesh,"bgmesh_$i",cellfields=["io"=>bgcell_to_ioc])
writevtk(get_grid(stl),"_stl_$i")
end

end # module
