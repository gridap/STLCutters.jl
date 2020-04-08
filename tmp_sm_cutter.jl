module sm_cutter

using STLCutter
using STLCutter: compute_surface_mesh_face_to_cells, get_vertex_coordinates, get_cell_coordinates, UNSET


function find_closest_background_face(bg_mesh::VolumeMesh,point::Point,sm_face_to_bg_cells::Table,vertex::Integer)
  min_dist = tol()
  closest_cell = UNSET
  for i in 1:length(sm_face_to_bg_cells,vertex)
    cell = sm_face_to_bg_cells[vertex,i]
    cell_coordinates = get_cell_coordinates(bg_mesh, cell )
    dist = distance(cell_coordinates,point)
    if dist < min_dist
      min_dist = dist
      closest_cell = cell
    end
  end
  closest_cell != UNSET || return (UNSET,UNSET)
  min_dist = tol()
  closest_face = UNSET
  for d in 0:D
    c_to_df = get_faces(bg_mesh,D,d)
    for ldface in 1:length(c_to_df,closest_cell)
      dface = c_to_df[closest_cell,ldface]
      dist = distance(bg_mesh,d,dface, point)
      if dist < min_dist
        min_dist = dist
        closest_face = dface
      end
    end
    if closest_face != UNSET
      return d, closest_face
    end
  end
  @assert false
end

tol() = 1e-10

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

bg_mesh = CartesianMesh( box, 1 ) 
vm = VolumeMesh( bg_mesh )
sm_face_to_bg_cells = compute_surface_mesh_face_to_cells(bg_mesh,sm)

display(sm_face_to_bg_cells)

point_type = eltype( get_vertex_coordinates(sm) )
D = num_dims(sm)

vertex_coordinates = point_type[]
d_to_dface_to_vertices = [ Table{Int}() for d in 0:D-1 ]
d_to_dface_to_sm_n = [ Int[] for d in 0:D-1 ]
d_to_dface_to_sm_nface = [ Int[] for d in 0:D-1 ]

vertex_to_bg_d = Int[]
vertex_to_bg_dface = Int[]

for vertex in 1:num_vertices(sm)
  point = get_vertex_coordinates(sm,vertex)
  _d, bg_face = find_closest_background_face( vm, point, sm_face_to_bg_cells, vertex )

  push!( vertex_coordinates, point )
  push!( vertex_to_bg_dface, bg_face )
  push!( vertex_to_bg_d, _d )
end

@show vertex_to_bg_d



#for d in 1:D-1
#  for SMFace in faces(d,surface_mesh)
#    
#    n = d-1
#    bg_cells = cells_containg_nfaces( bg_mesh, _sm, n )
#    
#    while length(bg_cells) > 0
#      bg_cell = pop!(bg_cells)
#        
#      ## 1) Fetch sm_face-bg_cell intersection already computed
#      ##    ( ∂Face_SM ∩ Cell_BgM ) ∩ ( 2) of neighbor cells )
#
#      nfaces = fetch_boundary_nfaces( bg_cell, _sm, n )
#
#      ## 2) Complete intersections in cell boundary, to define a close domain
#      ##    Face_SM ∩ ∂Cell_BgM 
#      
#      append!(nfaces, complete_boundary(bg_cell,_sm,SMFace) )
#
#      ## 3) Now we have a convex domain on n-dim space bounded by nfaces
#      ##    Create symplex mesh locally (dface = nface ∪ vertex):
#
#      add_faces!( sm, nfaces, n )
#      
#      ## *) Feed depth-first traversal
#
#      for _bg_cell in cells_touching( bg_mesh, nfaces ) #new nfaces actually, what is computed in 2)
#        push!( bg_cells, _bg_cell )
#      end
#
#    end
#  end
#end



end # module
