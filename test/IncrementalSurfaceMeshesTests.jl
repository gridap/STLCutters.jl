module IncrementalSurfaceMeshesTests 

using STLCutter

using Test

using STLCutter: cut_surface_mesh, surface, expand

using STLCutter: get_facet_coordinates, global_dface, face_dimension, local_dface, @check

points = [
  Point( 0.0, 0.0, 0.5 ), 
  Point( 1.0, 0.0, 0.5 ), 
  Point( 0.0, 1.0, 0.5 ), 
  Point( 1.0, 1.0, 0.5 ) ] 

c2v = Table( [ [ 1, 2, 3 ], [ 2, 3, 4 ] ] ) 

sm = SurfaceMesh( points, c2v )

box = BoundingBox(
  Point( 0.2, 0.2, 0.2 ),
  Point( 0.7, 0.7, 0.7 ) )

box = BoundingBox(
  - Point( 0.2, 0.2, 0.2 ),
    Point( 1.7, 1.7, 1.7 ) )

stl = STL(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))
#stl = STL(joinpath(@__DIR__,"data/cube.stl"))

sm = SurfaceMesh(stl)

box = BoundingBox(sm)

box = expand(box,0.1)

bg_mesh = CartesianMesh( box, 10 )

new_sm, = cut_surface_mesh(sm,bg_mesh)

writevtk(sm,"sm")
writevtk(bg_mesh,"bg_mesh")

writevtk(new_sm,"new_sm")

@test surface(sm) ≈ surface(new_sm)


end # module


#TODO:
# 
# [x] Use caches to avoid recursive memory memory allocation
# [x] Move to /test and /src folders
# [ ] Try to code it all dimension agnostic
# [x] Force same orientation in facets than their "parent"
# [x] Get (compute) the outputs: 
#   * [x] SurfaceMesh
#   * [x] sm_faces ↦ bg_faces (or inverse, if needed)
# [ ] Use the information in cell_mesh cutter
# [-] Test and debug with different surface mesh and bg mesh combinations
# [-] Add tests to the /test file
#
# DEBUG:
#
#  [x] Test not passing: overlapped faces
#
