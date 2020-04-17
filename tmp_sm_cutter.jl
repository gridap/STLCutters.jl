module sm_cutter

using STLCutter
using WriteVTK
using STLCutter: compute_surface_mesh_face_to_cells, get_vertex_coordinates, get_cell_coordinates, get_dface_to_vertices, get_face_coordinates, UNSET, @check
import STLCutter: closest_point, distance, get_dface_to_vertices, num_faces, num_vertices, writevtk, get_vertex_coordinates 

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

function add_vertex!(
  sm::IncrementalSurfaceMesh,
  bg_d::Integer,bg_face::Integer,point::Point,
  sm_d::Integer,sm_face::Integer)

  push!( sm.vertex_coordinates, point )
  push!( sm.vertex_to_bg_dface, bg_face )
  push!( sm.vertex_to_bg_d, bg_d )

  vertex = length( sm.vertex_coordinates )
  add_face!(sm,0,[vertex],sm_d,sm_face)
  vertex
end

function add_face!(sm::IncrementalSurfaceMesh,d::Integer,vertices::Vector,sm_d::Integer,sm_dface::Integer)
  push!( sm.d_to_dface_to_vertices[d+1], vertices )
  push!( sm.d_to_dface_to_sm_n[d+1], sm_d )
  push!( sm.d_to_dface_to_sm_nface[d+1], sm_dface )
  length( sm.d_to_dface_to_vertices[d+1] )
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

num_vertices(sm::IncrementalSurfaceMesh) = num_faces(sm,0)

get_vertex_coordinates(sm::IncrementalSurfaceMesh) = sm.vertex_coordinates

function get_surface_mesh_face(ism::IncrementalSurfaceMesh,d::Integer,dface::Integer)
  sm_d = ism.d_to_dface_to_sm_n[d+1][dface]
  sm_dface = ism.d_to_dface_to_sm_nface[d+1][dface]
  sm_d, sm_dface
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
  ism::IncrementalSurfaceMesh,
  n::Integer)

  resize!(cells,0)
  D = num_dims(m)
  sm_df_toism_nf = get_faces(sm,sm_d,n)
  sm_nf_toism_v = get_dface_to_vertices(sm,n)
  for lnface in 1:length( sm_df_toism_nf, sm_d )
    sm_nface = sm_df_toism_nf[ sm_face, lnface ]
    for lvertex in 1:length( sm_nf_toism_v, sm_nface )
      sm_vertex = sm_nf_toism_v[ sm_nface, lvertex ]
      d,dface = get_background_face(ism,sm_vertex)
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
  ism::IncrementalSurfaceMesh,n::Integer,nfaces::Vector,
  m::VolumeMesh,bg_cell::Integer,
  sm::SurfaceMesh,sm_d::Integer,sm_dface::Integer)

  resize!(nfaces,0)
  for nface in 1:num_faces(ism,n)
    if is_nface_in_surface_mesh_dface(ism,n,nface,sm,sm_d,sm_dface)
      if is_nface_in_bg_cell(nface,n,ism,m,bg_cell)
        push!(nfaces,nface)
      end
    end
  end
  nfaces
end

function is_nface_in_bg_cell(
  nface::Integer,n::Integer,ism::IncrementalSurfaceMesh,
  m::VolumeMesh,bg_cell::Integer)

  D = num_dims(m)
  nface_to_vertices = get_dface_to_vertices(ism,n)
  for lvertex in 1:length(nface_to_vertices,nface)
    vertex = nface_to_vertices[nface,lvertex]
    bg_d, bg_dface = get_background_face(ism,vertex)
    if !is_dface_in_cell(m,bg_d,bg_dface,bg_cell)
      return false
    end
  end
  true
end

function is_dface_in_cell(vm::VolumeMesh,d,dface,cell)
  D = num_dims(vm)
  is_nface_in_dface(vm,d,dface,D,cell)
end

function is_nface_in_surface_mesh_dface(
  ism::IncrementalSurfaceMesh,n::Integer,nface::Integer,
  sm::SurfaceMesh,sm_d::Integer,sm_dface::Integer)

  _d, _dface = get_surface_mesh_face(ism,n,nface)
  return is_nface_in_dface(sm,_d,_dface,sm_d,sm_dface)
end

function is_nface_in_dface(sm::SurfaceMesh,n::Integer,nface::Integer,d::Integer,dface)
  @check n ≤ d

  dface_to_nfaces = get_faces(sm,d,n)
  for lnface in 1:length(dface_to_nfaces,dface)
    _nface = dface_to_nfaces[dface,lnface]
    if _nface == nface
      return true
    end
  end
  false
end

##
function complete_boundary!(
  ism::IncrementalSurfaceMesh,n::Integer,nfaces::Vector,vertices::Vector,
  m::VolumeMesh,bg_cell::Integer,
  sm::SurfaceMesh,sm_d::Integer,sm_dface::Integer)

  vertices = get_vertices!(ism,vertices,n,nfaces)

  while length(vertices) > 0
    vertex = pop!(vertices)
    d,dface = get_background_face(ism,vertex)
    _d, _dface, point = _find_next_intersection(m,bg_cell,d,dface,sm,sm_d,sm_dface,ism,nfaces,n)
    if _dface != UNSET
      _vertex = add_vertex!(ism,_d,_dface,point,sm_d,sm_dface)
      push!( vertices, _vertex )

      ## TODO: Too hardcoded, think in any-dimension
      if n == 0
        nface = _vertex
      else
        nface = add_face!(ism, 1, [vertex,_vertex], sm_d, sm_dface ) 
      end
      push!( nfaces, nface )
    end
  end
end


function _find_next_intersection(
  m::VolumeMesh,bg_cell::Integer,bg_d::Integer,bg_face::Integer,
  sm::SurfaceMesh,sm_d::Integer,sm_dface::Integer,
  ism::IncrementalSurfaceMesh,nfaces::Vector,n::Integer)

  min_dist = tol()
  
  closest_bg_face = UNSET
  _d = UNSET
  D = num_dims(m)
  cell_to_facets = get_faces(m,D,D-1)
  for d in 0:D-sm_d
    facet_to_dfaces = get_faces(m,D-1,d)
    for lfacet in 1:length(cell_to_facets,bg_cell)
      bg_facet = cell_to_facets[bg_cell,lfacet]
      if sm_d == 1 || is_facet_around_dface(m,bg_facet,bg_d,bg_face) 
        for ldface in 1:length(facet_to_dfaces,bg_facet)
          _bg_dface = facet_to_dfaces[bg_facet,ldface]
          if !contains_any_nface(m,d,_bg_dface,ism,n,nfaces)
            dist = distance(m,d,_bg_dface,sm,sm_d,sm_dface)
            if dist < min_dist
              min_dist = dist
              closest_bg_face = _bg_dface
            end
          end
        end
      end
    end
    if closest_bg_face != UNSET
      _d = d
      break
    end
  end

  if closest_bg_face != UNSET
    point = closest_point(m,_d,closest_bg_face,sm,sm_d,sm_dface)
    return _d, closest_bg_face, point
  end
  
  return UNSET,UNSET, zero( eltype(get_vertex_coordinates(m)) )
end

function get_dfaces!(ism::IncrementalSurfaceMesh,n::Integer,nfaces::Vector,d::Integer,dfaces::Vector)
  ## TODO: Avoid global search, storing information
  resize!(dfaces,0)
  for dface in 1:num_faces(ism,d)
    if is_dface_in_any_nface(ism,d,dface,n,nfaces)
      push!(dfaces,dface)
    end
  end
  dfaces
end

function is_dface_in_any_nface(ism::IncrementalSurfaceMesh,d::Integer,dface::Integer,n::Integer,nfaces::Vector)
  @check d ≤ n
  if d == n
    if dface in nfaces
      return true
    else
      return false
    end
  end
  for nface in nfaces
    if is_nface_in_dface(ism,d,dface,n,nface)
      return true
    end
  end
  false
end

function is_nface_in_dface(ism::IncrementalSurfaceMesh,n::Integer,nface::Integer,d::Integer,dface::Integer)
  @check n ≤ d
  dface_to_vertices = get_dface_to_vertices(ism,n)
  for lvertex in 1:length(dface_to_vertices,nface)
    vertex = dface_to_vertices[nface,lvertex]
    if !is_vertex_in_dface(ism,vertex,d,dface)
      return false
    end
  end
  true
end

function get_vertices!(ism::IncrementalSurfaceMesh,vertices::Vector,n::Integer,nfaces::Vector)
  resize!(vertices,0)
  nface_to_vertices = get_dface_to_vertices(ism,n)
  for nface in nfaces
    for lvertex in 1:length(nface_to_vertices,nface)
      vertex = nface_to_vertices[nface,lvertex]
      if vertex ∉ vertices
        push!(vertices,vertex)
      end
    end
  end
  vertices
end

function is_facet_around_dface(m::VolumeMesh,facet::Integer,d::Integer,dface::Integer)
  D = num_dims(m)
  facet_to_vertices = get_faces(m,D-1,0)
  dface_to_vertices = get_faces(m,d,0)
  if d == D-1 && dface == facet
    return false
  end
  for lvertex in 1:length(facet_to_vertices,facet)
    vertex = facet_to_vertices[facet,lvertex]
    for _lvertex in 1:length(dface_to_vertices,dface)
      _vertex = dface_to_vertices[dface,_lvertex]
      if vertex == _vertex
        return true
      end
    end
  end
  false
end

function contains_any_nface(
  m::VolumeMesh,d::Integer,dface::Integer,
  ism::IncrementalSurfaceMesh,n::Integer,nfaces::Vector)
  
  nface_to_vertices = get_dface_to_vertices(ism,n)
  for nface in nfaces
    for lvertex in 1:length(nface_to_vertices,nface)
      vertex = nface_to_vertices[nface,lvertex]
      _d, _dface = get_background_face(ism,vertex)
      if is_nface_in_dface(m,_d,_dface,d,dface)
        return true
      end
    end
  end
  false
end

function is_nface_in_dface(m::VolumeMesh,n::Integer,nface::Integer,d::Integer,dface::Integer)
  if n ≤ d
    dface_to_nfaces = get_faces(m,d,n)
    for lnface in 1:length(dface_to_nfaces,dface)
      _nface = dface_to_nfaces[dface,lnface]
      if nface == _nface
        return true
      end
    end
  end
  false
end


function add_faces!(
  ism::IncrementalSurfaceMesh,n::Integer,nfaces::Vector,
  sm::SurfaceMesh,sm_d::Integer,sm_dface::Integer,
  vm::VolumeMesh,bg_cell::Integer)
  
  _nfaces = Int[]
  fetch_boundary_nfaces!(ism,n,_nfaces,vm,bg_cell,sm,sm_d,sm_dface)
  if length(nfaces) == 0
    return
  end
  nface_to_vertices = get_dface_to_vertices(ism,n)
  vertex = nface_to_vertices[ _nfaces[1], 1 ]
  lvertices = Int[]
  dfaces = Int[]
@show nfaces
  for d in reverse(0:sm_d-1)
   # dfaces = get_dfaces!(ism,n,nfaces,d,dfaces)
    dfaces = fetch_boundary_nfaces!(ism,d,dfaces,vm,bg_cell,sm,sm_d,sm_dface)
    @show sm_d,d,dfaces
    dface_to_vertices = get_dface_to_vertices(ism,d)
    resize!( lvertices, d+2 )
    for dface in dfaces
      if !is_dface_connected_to_vertex(ism,n,nfaces,d,dface,vertex)
        lvertices[1] = vertex
        for i in 1:length(dface_to_vertices,dface)
          lvertices[i+1] = dface_to_vertices[dface,i]
        end
        @show d+1,lvertices
        add_face!(ism,d+1,lvertices,sm_d,sm_dface)
      end
    end
  end
end

function is_dface_connected_to_vertex(
  ism::IncrementalSurfaceMesh,
  n::Integer,nfaces::Vector,
  d::Integer,dface::Integer,
  vertex::Integer)

  dface_to_vertices = get_dface_to_vertices(ism,d)
  for nface in nfaces
    if is_vertex_in_dface(ism,vertex,n,nface)
      if is_nface_in_dface(ism,d,dface,n,nface)
        return true
      end
    end
  end
  false
end

function is_vertex_in_dface(ism::IncrementalSurfaceMesh,vertex::Integer,n::Integer,nface::Integer)
  nface_to_vertices = get_dface_to_vertices(ism,n)
  for lvertex in 1:length(nface_to_vertices,nface)
    _vertex = nface_to_vertices[nface,lvertex]
    if vertex == _vertex
      return true
    end
  end
  false
end


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


function is_dface_in_background_facet(
  ism::IncrementalSurfaceMesh,d::Integer,dface::Integer,
  vm::VolumeMesh,facet::Integer)

  D = num_dims(vm)
  dface_to_vertices = get_dface_to_vertices(ism,d)
  for lvertex in 1:length(dface_to_vertices,dface)
    vertex = dface_to_vertices[dface,lvertex]
    bg_d, bg_dface = get_background_face(ism,vertex)
    if bg_d > D-1 || !is_dface_in_facet(vm,bg_d,bg_dface,facet)
      return false
    end
  end
  true
end

function is_dface_in_facet(vm::VolumeMesh,d::Integer,dface::Integer,facet::Integer) 
  D = num_dims(vm)
  @check d ≤ D-1 
  facet_to_dfaces = get_faces(vm,D-1,d)
  for ldface in 1:length(facet_to_dfaces,facet)
    _dface = facet_to_dfaces[facet,ldface]
    if _dface == dface
      return true
    end
  end
  false
end

function connect_boundary_faces!(
  vm::VolumeMesh,bg_cell::Integer,
  ism::IncrementalSurfaceMesh,n::Integer,nfaces::Vector,
  sm_d::Integer,sm_dface::Integer)

  D = num_dims(vm)
  cell_to_facets = get_faces(vm,D,D-1)
  # in more dimensions D > 3 this would go from d = 2 to D-1 bg_faces in the current bg_cell
  # need to write something more generic
  # -> whiteboard algorithm 

  
  n == D-2 || return
  
  vertices = Int[] ## <- Save load caches
  _vertices = Int[] ## <- Save load caches
  
  vertices = get_vertices!(ism,vertices,n,nfaces)
  for lfacet in 1:length(cell_to_facets,bg_cell) 
    bg_facet = cell_to_facets[bg_cell,lfacet]
    _vertices = get_vertices_in_background_facet(ism,vertices,vm,bg_facet,_vertices)
  
    @check length(_vertices) ≤ n+1
    if length(_vertices) == n+1
      if !are_vertices_in_any_dface(ism,_vertices,n,nfaces)
        nface = add_face!(ism,n,_vertices,sm_d,sm_dface)
        push!(nfaces,nface)
      end
    end
  end
end


function get_vertices_in_background_facet(
  ism::IncrementalSurfaceMesh,vertices::Vector,
  vm::VolumeMesh,bg_facet::Integer,
  _vertices::Vector)

  resize!(_vertices,0)
  for vertex ∈ vertices
    if is_vertex_in_background_facet(ism,vertex,vm,bg_facet)
      push!(_vertices,vertex)
    end
  end
  _vertices
end

function is_vertex_in_background_facet(ism::IncrementalSurfaceMesh,vertex::Integer,vm::VolumeMesh,facet::Integer)
  is_dface_in_background_facet(ism,0,vertex,vm,facet)
end


function are_vertices_in_any_dface(ism::IncrementalSurfaceMesh,vertices::Vector,d::Integer,dfaces::Vector)
  for dface in dfaces
    if are_vertices_in_dface(ism,vertices,d,dface)
      return true
    end
  end
  false
end

function are_vertices_in_dface(ism::IncrementalSurfaceMesh,vertices::Vector,d::Integer,dface::Integer)
  dface_to_vertices = get_dface_to_vertices(ism,d)
  for lvertex in 1:length(dface_to_vertices,dface)
    vertex = dface_to_vertices[dface,lvertex]
    if vertex ∉ vertices
      return false
    end
  end
  true
end

function writevtk(sm::IncrementalSurfaceMesh{D,T},file_base_name) where {D,T}
  d_to_vtk_type_id = Dict(0=>1,1=>3,2=>5,3=>10)
  num_points = num_vertices(sm)
  points = zeros(T,D,num_points)
  for (i ,p ) in enumerate( get_vertex_coordinates(sm) ), d in 1:D
    points[d,i] = p[d]
  end
  cells = MeshCell{Vector{Int64}}[]
  for d in 2 #0:D-1
    dface_to_vertices = get_dface_to_vertices(sm,d)
    vtk_type = VTKCellType(d_to_vtk_type_id[d])
    for i in 1:num_faces(sm,d)
      vertices = [ dface_to_vertices[i,j] for j in 1:vtk_type.nodes ]
      push!( cells, MeshCell(vtk_type,vertices) )
    end
  end
  vtkfile = vtk_grid(file_base_name,points,cells)
  vtk_save(vtkfile)
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
  add_vertex!(ism,bg_d,bg_face,point,0,vertex)
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

end # module




#TODO:
# 
# Use caches to avoid recursive memory memory allocation
# Move to /test and /src folders
# Try to code it all dimension agnostic
# Force same orientation in facets than their "parent"
# Get (compute) the outputs: 
#   * SurfaceMesh
#   * sm_faces ↦ bg_faces (or inverse, if needed)
# Use the information in cell_mesh cutter
# Test and debug with different surface mesh and bg mesh combinations
#
