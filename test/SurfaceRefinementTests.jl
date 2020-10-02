module SurfaceRefinementTests

using Test
using Gridap
using STLCutters


using STLCutters: compute_stl_model
using STLCutters: refine_stl_face
using STLCutters: compute_face_to_cells
using STLCutters: refine_surface 


stl_vertices = [
  Point(-0.5,0.5),
  Point(0.5,0.5) ]
stl_faces = Table([[1,2]])
stl = compute_stl_model(stl_faces,stl_vertices)
origin = Point(0,0)
sizes = (1,1)
partition = (1,1)
grid = CartesianGrid(origin,sizes,partition)
T,X = refine_stl_face(grid,1,stl,3)

stl_vertices = [
  Point(-0.5,0.5,0.5),
  Point(0.5,0.5,0.5),
  Point(0.5,-0.5,0.5) ]
stl_faces = Table([[1,2,3]])
stl = compute_stl_model(stl_faces,stl_vertices)

origin = Point(0,0,0)
sizes = (1,1,1)
partition = (1,1,1)
grid = CartesianGrid(origin,sizes,partition)
T,X = refine_stl_face(grid,1,stl,7)

origin = Point(0,0,0)
sizes = (0.2,0.2,0.2)
partition = (3,3,3)
grid = CartesianGrid(origin,sizes,partition)
stlf_to_c,c_to_stlf = compute_face_to_cells(grid,stl)
stl,f_to_bgcell = refine_surface(grid,stl,stlf_to_c)
writevtk(grid,"bg_mesh")
writevtk(get_grid(stl),"surface",cellfields=["bgcell"=>f_to_bgcell])

# TODO: Testing

end # module
