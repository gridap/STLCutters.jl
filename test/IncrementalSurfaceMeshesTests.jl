module IncrementalSurfaceMeshesTests 

using STLCutter

using Test

using STLCutter: cut_surface_mesh

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

bg_mesh = CartesianMesh( box, 2 )

new_sm, d_to_df_to_smf = cut_surface_mesh(sm,bg_mesh)

writevtk(sm,"sm")
writevtk(bg_mesh,"bg_mesh")

display(d_to_df_to_smf)
writevtk(new_sm,"new_sm")

end # module


#TODO:
# 
# [ ] Use caches to avoid recursive memory memory allocation
# [x] Move to /test and /src folders
# [ ] Try to code it all dimension agnostic
# [x] Force same orientation in facets than their "parent"
# [x] Get (compute) the outputs: 
#   * [x] SurfaceMesh
#   * [x] sm_faces ↦ bg_faces (or inverse, if needed)
# [ ] Use the information in cell_mesh cutter
# [ ] Test and debug with different surface mesh and bg mesh combinations
# [ ] Add tests to the /test file
#
