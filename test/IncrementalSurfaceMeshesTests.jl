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

#stl = STL(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))
stl = STL(joinpath(@__DIR__,"data/cube.stl"))
stl = STL(joinpath(@__DIR__,"data/wine_glass.stl"))

sm = SurfaceMesh(stl)

box = BoundingBox(sm)

box = expand(box,0.1)

bg_mesh = CartesianMesh( box, 20 )

new_sm, d_to_bg_df_to_smf, new_sm_face_to_sm_face = cut_surface_mesh(sm,bg_mesh)

writevtk(sm,"sm")
writevtk(bg_mesh,"bg_mesh")

writevtk(new_sm,"new_sm")

@test surface(sm) ≈ surface(new_sm)

@test is_watter_tight(sm)

end # module

