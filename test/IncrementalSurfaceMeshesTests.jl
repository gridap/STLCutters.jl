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
#stl = STL(joinpath(@__DIR__,"data/wine_glass.stl"))

sm = SurfaceMesh(stl)

box = BoundingBox(sm)

box = expand(box,0.1)

bg_mesh = CartesianMesh( box, 20 )

new_sm, d_to_bg_df_to_smf, new_sm_face_to_sm_face = cut_surface_mesh(sm,bg_mesh)

#summed_surface = zeros( num_facets(sm) )
#D = num_dims(sm)
#for facet in 1:num_facets(new_sm)
#  facet_coordinates = get_facet_coordinates(new_sm,facet)
#  face = global_dface(new_sm,D-1,facet)
#  old_face = new_sm_face_to_sm_face[face]
#  @check face_dimension(sm,old_face) == D-1
#  old_facet = local_dface(sm,old_face,D-1)
#
#  summed_surface[old_facet] += measure(facet_coordinates)
#end
#
#
#for facet in 1:num_facets(sm)
#  facet_coordinates = get_facet_coordinates(sm,facet)
#  if measure(facet_coordinates) ≉ summed_surface[facet]
#    @show facet
#    @show facet_coordinates
#    @show measure(facet_coordinates), summed_surface[facet]
#  end 
#end
# 
#i = 0
#for facet in 1:num_facets(new_sm)
#  global i
#  facet_coordinates = get_facet_coordinates(new_sm,facet)
#  face = global_dface(new_sm,D-1,facet)
#  old_face = new_sm_face_to_sm_face[face]
#  @check face_dimension(sm,old_face) == D-1
#  old_facet = local_dface(sm,old_face,D-1)
#  if old_facet == 5
#    @show facet, facet_coordinates
#    i += 1
#    @show i, facet
#    writevtk( facet_coordinates, "t$i" )
#  end
#end
#
#bg_cell_to_sm_face = d_to_bg_df_to_smf[D+1]
#@show global_dface(new_sm,D-1,78)
#@show global_dface(new_sm,D-1,79)
#@show global_dface(new_sm,D-1,81)
#for bg_cell in 1:num_cells(bg_mesh)
#  for lface in 1:length(bg_cell_to_sm_face,bg_cell)
#    sm_face = bg_cell_to_sm_face[bg_cell,lface]
#    if sm_face ∈ 
#      [ global_dface(new_sm,D-1,78),
#        global_dface(new_sm,D-1,79),
#        global_dface(new_sm,D-1,81) ]
#      @show bg_cell
#    end
#  end
#end

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
