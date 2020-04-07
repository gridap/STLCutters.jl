module sm_cutter

using STLCutter
using STLCutter: compute_surface_mesh_face_to_cells

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

bg_mesh = CartesianMesh( box, 1 ) 

sm_face_to_bg_cells = compute_surface_mesh_face_to_cells(bg_mesh,sm)

display(sm_face_to_bg_cells)

point_type = eltype( get_vertices(sm) )

vertex_coordinates = point_type[]
d_to_dface_to_vertices = [ Table{Int}() for d in 0:D-1 ]
d_to_dface_to_sm_n = [ Int[] for d in 0:D-1 ]
d_to_dface_to_sm_nface = [ Int[] for d in 0:D-1 ]

vertex_to_bg_d = Int[]
vertex_to_bg_dface = Int[]

for vertex in 1:num_vertices(sm)
  point = get_vertex_coordinates(sm,vertex)
  bg_face, _d = find_closest_background_face( bg_mesh, point )

  push!( vertex_coordinates, point )
  push!( vertex_to_bg_dface, bg_face )
  push!( vertex_to_bg_d, _d )
end


function find_closest_background_face(bg_mesh::CartesianMesh,point::Point)
  min_dist = tol()
  closest_cell = UNSET
  for i in length(sm_face_to_bg_cells,vertex)
    cell = sm_face_to_bg_cells[vertex,i]
    cell_coordinates = get_cell(bg_mesh, cell )
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
    for ldface in 1:num_dfaces(bg_mesh,d)
      dface = get_face(bg_mesh,ldface,d)
      dist = distance(bg_mesh,dface,d, point)
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


function distance(box::BoundingBox,p::Point)
  if have_intersection(box,p) 
    min_dist = 0.0
  else
    min_dist = typemax( 0.0 )
    for iface in 1:num_faces(box)
      face = get_face(box,iface)
      dist = distance(face,p)
      if dist < min_dist
        min_dist = dist
      end
    end
  end
  min_dist
end


end # module
