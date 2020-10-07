module GridRefinementTests
  
using Test
using Gridap
using STLCutters

using Gridap.Geometry

using STLCutters: compute_stl_model
using STLCutters: refine_grid 
using STLCutters: volume,volumes
using STLCutters: FACE_CUT, FACE_IN, FACE_OUT


vertices = [
  Point(-0.1,0.5),
  Point(0.1,0.8),
  Point(0.4,0.4),
  Point(0.7,0.7),
  Point(0.8,0.3),
  Point(1.1,0.3) ]
faces = [[1,2],[2,3],[3,4],[4,5],[5,6]]

stl = compute_stl_model(Table(faces),vertices)
bg_mesh = compute_linear_grid(QUAD4)
bg_mesh = CartesianGrid( Point(0.,0.), (0.3,0.3), (3,3) )

out = refine_grid(bg_mesh,stl)
T,X,reffes,cell_types,cell_to_io,cell_to_bgcell,bgcell_to_ioc = out 

submesh = UnstructuredGrid(X,Table(T),reffes,cell_types)
vol = volume(submesh) + sum(volumes(bg_mesh).*(bgcell_to_ioc.≠FACE_CUT))
@test volume(bg_mesh) ≈ vol

writevtk(submesh,"submesh",cellfields=["io"=>cell_to_io,"bgcell"=>cell_to_bgcell])
writevtk(bg_mesh,"bgmesh",cellfields=["io"=>bgcell_to_ioc])
writevtk(stl,"stl")

vertices = [
  Point(-2.0,-2.0,3.0),
  Point(0.4,-2.0,0.3),
  Point(0.4,2.0,0.3),
  Point(3.0,3.0,3.0) ]
facets = [[1,2,3],[2,4,3]]

stl = compute_stl_model(Table(facets),vertices)
bg_mesh = compute_linear_grid(HEX8)
bg_mesh = CartesianGrid( Point(0.,0.,0.), (0.3,0.3,0.3), (3,3,3) )

out = refine_grid(bg_mesh,stl)
T,X,reffes,cell_types,cell_to_io,cell_to_bgcell,bgcell_to_ioc = out 

submesh = UnstructuredGrid(X,Table(T),reffes,cell_types)
vol = volume(submesh) + sum(volumes(bg_mesh).*(bgcell_to_ioc.≠FACE_CUT))
@test volume(bg_mesh) ≈ vol

writevtk(submesh,"submesh",cellfields=["io"=>cell_to_io,"bgcell"=>cell_to_bgcell])
writevtk(bg_mesh,"bgmesh",cellfields=["io"=>bgcell_to_ioc])
writevtk(stl,"stl")


vertices = [
  Point(-0.1,0.5,0.5),
  Point(0.5,1.1,0.5),
  Point(1.1,0.5,0.5),
  Point(0.5,0.1,0.6),
  Point(1.1,1.1,0.6),
  Point(-0.1,1.1,0.5),
  Point(-0.1,-0.1,0.3),
  Point(0.5,-0.1,0.2),
  Point(1.1,-0.1,0.3) ]
facets = [[1,2,3],[1,3,4],[1,6,2],[3,2,5],[1,4,7],[4,8,7],[4,9,8],[3,9,4]]

stl = compute_stl_model(Table(facets),vertices)
bg_mesh = compute_linear_grid(HEX8)
bg_mesh = CartesianGrid( Point(0.,0.,0.), (0.3,0.3,0.3), (3,3,3) )

out = refine_grid(bg_mesh,stl)
T,X,reffes,cell_types,cell_to_io,cell_to_bgcell,bgcell_to_ioc = out 

submesh = UnstructuredGrid(X,Table(T),reffes,cell_types)
vol = volume(submesh) + sum(volumes(bg_mesh).*(bgcell_to_ioc.≠FACE_CUT))
@test volume(bg_mesh) ≈ vol

writevtk(submesh,"submesh",cellfields=["io"=>cell_to_io,"bgcell"=>cell_to_bgcell])
writevtk(bg_mesh,"bgmesh",cellfields=["io"=>bgcell_to_ioc])
writevtk(stl,"stl")


end # module
