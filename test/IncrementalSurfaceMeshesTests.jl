module IncrementalSurfaceMeshesTests 

using STLCutter

using Test

using STLCutter: IncrementalSurfaceMesh

using STLCutter: num_faces, num_vertices, get_vertex_coordinates, get_dface_to_vertices, UNSET 

using STLCutter: find_closest_background_face, cells_containg_nfaces!, add_vertex!,  fetch_boundary_nfaces!, complete_boundary!, connect_boundary_faces!, add_faces!

using STLCutter: is_dface_in_background_facet, is_nface_around_dface

using STLCutter: compute_surface_mesh_face_to_cells 

using STLCutter: background_face_to_surface_mesh_faces

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
vm = VolumeMesh( bg_mesh )
sm_face_to_bg_cells = compute_surface_mesh_face_to_cells(bg_mesh,sm)

display(sm_face_to_bg_cells)

point_type = eltype( get_vertex_coordinates(sm) )
D = num_dims(sm)
T = Float64

ism = IncrementalSurfaceMesh{D,T}()
bg_cells = Int[]
nfaces = Int[]
vertices = Int[]

cell_to_facets = get_faces(vm,D,D-1)
facet_to_cells = get_faces(vm,D-1,D)

for vertex in 1:num_vertices(sm)
  point = get_vertex_coordinates(sm,vertex)
  bg_d, bg_face = find_closest_background_face( vm, point, sm_face_to_bg_cells, vertex )
  add_vertex!(ism,bg_d,bg_face,0,vertex,point)
end

for sm_d in 1:D-1
  for sm_dface in 1:num_faces(sm,sm_d)
    n = sm_d-1
    cells_containg_nfaces!(bg_cells,vm,sm,sm_d,sm_dface,ism,n)
    while length(bg_cells) > 0
      bg_cell = pop!(bg_cells)
      fetch_boundary_nfaces!(ism,n,nfaces,vm,bg_cell,sm,sm_d,sm_dface)
      num_nfaces = length(nfaces)
      complete_boundary!(ism,n,nfaces,vertices,vm,bg_cell,sm,sm_d,sm_dface)
      connect_boundary_faces!(vm,bg_cell,ism,n,nfaces,sm_d,sm_dface)
      add_faces!(ism,n,nfaces,sm,sm_d,sm_dface,vm,bg_cell)
      for i in num_nfaces+1:length(nfaces)
        nface = nfaces[i]
        for lfacet in 1:length(cell_to_facets,bg_cell)
          bg_facet = cell_to_facets[bg_cell,lfacet]
          if is_dface_in_background_facet(ism,n,nface,vm,bg_facet)
            for lcell in 1:length(facet_to_cells,bg_facet)
              _bg_cell = facet_to_cells[bg_facet,lcell]
              if _bg_cell != bg_cell && _bg_cell ∉ bg_cells
                push!(bg_cells,_bg_cell)
              end
            end
          end
        end
      end
    end
  end
end

@show get_dface_to_vertices(sm,0)
@show get_dface_to_vertices(sm,1)
@show get_dface_to_vertices(sm,2)

@show ism.vertex_to_bg_d
@show ism.vertex_to_bg_dface

@show get_dface_to_vertices(ism,0)
@show get_dface_to_vertices(ism,1)
@show get_dface_to_vertices(ism,2)



writevtk(sm,"sm")
writevtk(ism,"ism")
writevtk(bg_mesh,"bg_mesh")

sm = SurfaceMesh(ism.vertex_coordinates,ism.d_to_dface_to_vertices)
d_to_df_to_smf = background_face_to_surface_mesh_faces(ism,vm)

display(d_to_df_to_smf)
writevtk(sm,"new_sm")

end # module




#TODO:
# 
# [ ] Use caches to avoid recursive memory memory allocation
# [x] Move to /test and /src folders
# [ ] Try to code it all dimension agnostic
# [ ] Force same orientation in facets than their "parent"
# [x] Get (compute) the outputs: 
#   * [x] SurfaceMesh
#   * [x] sm_faces ↦ bg_faces (or inverse, if needed)
# [ ] Use the information in cell_mesh cutter
# [ ] Test and debug with different surface mesh and bg mesh combinations
# [ ] Add tests to the /test file
#
