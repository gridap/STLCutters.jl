module sm_cutter

using STLCutter
using STLCutter: compute_surface_mesh_face_to_cells, get_vertex_coordinates, get_cell_coordinates, get_dface_to_vertices, UNSET
import STLCutter: closest_point, distance, get_dface_to_vertices, num_faces

struct IncrementalSurfaceMesh{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  d_to_dface_to_vertices::Vector{Table{Int}}
  d_to_dface_to_sm_n::Vector{Vector{Int8}}
  d_to_dface_to_sm_nface::Vector{Vector{Int}}
  vertex_to_bg_d::Vector{Int8}
  vertex_to_bg_dface::Vector{Int}
end

function IncrementalSurfaceMesh{D,T}() where {D,T}
  vertex_coordinates = Point{D,T}[]
  d_to_dface_to_vertices = [ Table{Int}() for d in 0:D-1 ]
  d_to_dface_to_sm_n = [ Int8[] for d in 0:D-1 ]
  d_to_dface_to_sm_nface = [ Int[] for d in 0:D-1 ]
  vertex_to_bg_d = Int8[]
  vertex_to_bg_dface = Int[]
  IncrementalSurfaceMesh(
    vertex_coordinates,
    d_to_dface_to_vertices,
    d_to_dface_to_sm_n,
    d_to_dface_to_sm_nface,
    vertex_to_bg_d,
    vertex_to_bg_dface )
end

function add_vertex!(sm::IncrementalSurfaceMesh,bg_d::Integer,bg_face::Integer,point::Point,sm_d::Integer,sm_face::Integer)
  push!( sm.vertex_coordinates, point )
  push!( sm.vertex_to_bg_dface, bg_face )
  push!( sm.vertex_to_bg_d, bg_d )

  vertex = length( sm.vertex_coordinates )
  add_face!(sm,0,[vertex],sm_d,sm_face)
end

function add_face!(sm::IncrementalSurfaceMesh,d::Integer,vertices::Vector,sm_d::Integer,sm_dface::Integer)
  push!( sm.d_to_dface_to_vertices[d+1], vertices )
  push!( sm.d_to_dface_to_sm_n[d+1], sm_d )
  push!( sm.d_to_dface_to_sm_nface[d+1], sm_dface )
end

function get_background_face(sm::IncrementalSurfaceMesh,vertex::Integer)
  sm.vertex_to_bg_d[vertex], sm.vertex_to_bg_dface[vertex]
end

function get_dface_to_vertices(sm::IncrementalSurfaceMesh,d::Integer)
  sm.d_to_dface_to_vertices[d+1]
end

function num_faces(sm::IncrementalSurfaceMesh,d::Integer)
  length( get_dface_to_vertices(sm,d) )
end

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

function cells_containg_nfaces!(
  cells::Vector,
  m::VolumeMesh,
  sm::SurfaceMesh,
  sm_d::Integer,
  sm_face::Integer,
  n::Integer,
  _sm::IncrementalSurfaceMesh)

  resize!(cells,0)
  D = num_dims(m)
  sm_df_to_sm_nf = get_faces(sm,sm_d,n)
  sm_nf_to_sm_v = get_dface_to_vertices(sm,n)
  for lnface in 1:length( sm_df_to_sm_nf, sm_d )
    sm_nface = sm_df_to_sm_nf[ sm_face, lnface ]
    for lvertex in 1:length( sm_nf_to_sm_v, sm_nface )
      sm_vertex = sm_nf_to_sm_v[ sm_nface, lvertex ]
      d,dface = get_background_face(_sm,sm_vertex)
      if dface != UNSET
        df_to_c = get_faces(m,d,D)
        for i in 1:length(df_to_c,dface)
          cell = df_to_c[dface,i]
          if cell ∉ cells
            push!(cells,cell)
          end
        end
      end
    end
  end
  cells
end

function fetch_boundary_nfaces!(
  nfaces::Vector,n::Integer,_sm::IncrementalSurfaceMesh,
  m::VolumeMesh,bg_cell::Integer)

  resize!(nfaces,0)
  for nface in 1:num_faces(_sm,n)
    if is_nface_in_bg_cell(nface,n,_sm,m,bg_cell)
      push!(nfaces,nface)
    end
  end
  nfaces
end

function is_nface_in_bg_cell(
  nface::Integer,n::Integer,_sm::IncrementalSurfaceMesh,
  m::VolumeMesh,bg_cell::Integer)

  D = num_dims(m)
  nface_to_vertices = get_dface_to_vertices(_sm,n)
  for lvertex in 1:length(nface_to_vertices,nface)
    vertex = nface_to_vertices[nface,lvertex]
    bg_d, bg_dface = get_background_face(_sm,vertex)
    bg_c_to_bg_df = get_faces(m,D,bg_d)
    for ldface in 1:length(bg_c_to_bg_df,bg_cell)
      _bg_dface = bg_c_to_bg_df[bg_cell,ldface]
      if _bg_dface == bg_dface
        return true
      end
    end
  end
  false
end

#function complete_boundary( )
#  sm_d
#  n
#
#  vertices = get_vertices(nfaces)
#
#  while length(vertices) > 0
#    vertex = pop!(vertices)
#    bg_dface = vertex_to_bg_d
#    bg_d = vertex_to_bg_dface
#    _find_next_intersection()
#  end
#  
#
#end
#
#
#function _find_next_intersection(bg_cell...)
#  min_dist = tol()
#  closest_bg_face = UNSET
#
#  c_to_f = get_faces(m,D,D-1)
#  for d in 0:D-sm_d
#    f_to_df = get_faces(m,D,D-1)
#    for lfacet in 1:length(c_to_f,bg_cell)
#      bg_facet = c_to_f[bg_cell,lfacet]
#
#      if sm_d == 1 || is_around(m,facet,bg_face,bg_d) 
#        for ldface in 1:lenth(f_to_df,bg_facet)
#          _bg_dface = f_to_df[bg_facet,ldface]
#
#          if contains_any_nface( vertices, d, _bg_dface) 
#            
#            dist = distance(m,d,_bg_dface,sm,sm_d,sm_dface)
#          if dist < min_dist
#            min_dist = dist
#            closest_bg_face = _bg_dface
#          end
#        end
#      end
#    end
#    if closest_bg_face != UNSET
#      _d = d
#      break
#    end
#  end
#
#  if closest_bg_face != UNSET
#    point = closest_point(bg_m,_d,closest_bg_face,sm,sm_d,sm_dface)
#    return closest_bg_face,_d,point
#  end
#  
#  return UNSET,UNSET, zero( eltype(get_vertex_coordinates(m)) )
#end
#
#function get_vertices(nfaces)
#  vertices = Int[]
#  for nface in nfaces
#    for lvertex in 1:length(nface_to_vertices,nface)
#      vertex = nface_to_vertices[nface,lvertex]
#      if vertex ∉ vertices
#        push!(vertices,vertex)
#      end
#    end
#  end
#  vertices
#end

@generated function distance(
  m::VolumeMesh{D},d::Integer,dface::Integer,
  sm::SurfaceMesh{D},n::Integer,sm_nface::Integer) where D

  msg1 = "\$d-face does not exist at \$(typeof(m))"
  msg2 = "\$n-face does not exist at \$(typeof(sm))"
  msg3 = "distance() from \$d-face in \$(typeof(m)) against \$n-face in \$(typeof(sm)) not implemented"
  
  body = ""
  for d in 0:D
    sm_body = ""
    body *= "if d == $d \n"
    body *= "  f$d = get_face_coordinates(m,Val{$d}(),dface) \n"
    for n in 0:min(D-d,D-1)
      sm_body *= "if n == $n \n"
      sm_body *= "    sm_f$n = get_face_coordinates(sm,Val{$n}(),sm_nface) \n"
      sm_body *= "    distance(f$d,sm_f$n)\n  "
      sm_body *= "else"
    end
    error = "if n > $(D-1) \n"
    error *= "    throw(ArgumentError(\"$msg2\"))\n  "
    error *= "else \n"
    error *= "    throw(ArgumentError(\"$msg3\")) \n  "
    error *= "end"
    condition = sm_body * error
    body *= "  $condition \n"
    body *= "else"
  end
  error = "  throw(ArgumentError(\"$msg1\"))\n"
  str = body * '\n' * error * "end"

  Meta.parse(str)
end

@generated function closest_point(
  m::VolumeMesh{D},d::Integer,dface::Integer,
  sm::SurfaceMesh{D},n::Integer,sm_nface::Integer) where D

  msg1 = "\$d-face does not exist at \$(typeof(m))"
  msg2 = "\$n-face does not exist at \$(typeof(sm))"
  msg3 = "closest_point() from \$d-face in \$(typeof(m)) against \$n-face in \$(typeof(sm)) not implemented"
  
  body = ""
  for d in 0:D
    sm_body = ""
    body *= "if d == $d \n"
    body *= "  f$d = get_face_coordinates(m,Val{$d}(),dface) \n"
    for n in 0:min(D-d,D-1)
      sm_body *= "if n == $n \n"
      sm_body *= "    sm_f$n = get_face_coordinates(sm,Val{$n}(),sm_nface) \n"
      sm_body *= "    closest_point(f$d,sm_f$n)\n  "
      sm_body *= "else"
    end
    error = "if n > $(D-1) \n"
    error *= "    throw(ArgumentError(\"$msg2\"))\n  "
    error *= "else \n"
    error *= "    throw(ArgumentError(\"$msg3\")) \n  "
    error *= "end"
    condition = sm_body * error
    body *= "  $condition \n"
    body *= "else"
  end
  error = "  throw(ArgumentError(\"$msg1\"))\n"
  str = body * '\n' * error * "end"

  Meta.parse(str)
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
T = Float64

_sm = IncrementalSurfaceMesh{D,T}()

vertex_coordinates = point_type[]
d_to_dface_to_vertices = [ Table{Int}() for d in 0:D-1 ]
d_to_dface_to_sm_n = [ Int[] for d in 0:D-1 ]
d_to_dface_to_sm_nface = [ Int[] for d in 0:D-1 ]

vertex_to_bg_d = Int[]
vertex_to_bg_dface = Int[]

bg_cells = Int[]
nfaces = Int[]


for vertex in 1:num_vertices(sm)
  point = get_vertex_coordinates(sm,vertex)
  _d, bg_face = find_closest_background_face( vm, point, sm_face_to_bg_cells, vertex )
  
  add_vertex!(_sm,_d,bg_face,point,0,vertex)

  push!( vertex_coordinates, point )
  push!( vertex_to_bg_dface, bg_face )
  push!( vertex_to_bg_d, _d )
  push!( d_to_dface_to_vertices[0+1], [vertex] )
  push!( d_to_dface_to_sm_n[0+1], 0 )
  push!( d_to_dface_to_sm_nface[0+1], vertex )
end

@show vertex_to_bg_d


for d in 1:D-1
  for sm_dface in 1:num_faces(sm,d)
    n = d-1
    cells_containg_nfaces!(bg_cells,vm,sm,d,sm_dface,n,_sm)
    @show bg_cells
    while length(bg_cells) > 0
      bg_cell = pop!(bg_cells)
      fetch_boundary_nfaces!(nfaces,n,_sm,vm,bg_cell)
      @show nfaces
    end
  end
end

end # module
