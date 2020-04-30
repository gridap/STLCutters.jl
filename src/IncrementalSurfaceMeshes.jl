
struct IncrementalSurfaceMesh{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  d_to_dface_to_vertices::Vector{Table{Int}}
  d_to_dface_to_sm_n::Vector{Vector{Int8}}
  d_to_dface_to_sm_nface::Vector{Vector{Int}}
  sm_face_to_ds::Vector{Vector{Int8}}
  sm_face_to_dfaces::Vector{Vector{Int}}
  vertex_to_bg_d::Vector{Int8}
  vertex_to_bg_dface::Vector{Int}
end

function IncrementalSurfaceMesh{D,T}(num_faces) where {D,T}
  vertex_coordinates = Point{D,T}[]
  d_to_dface_to_vertices = [ Table{Int}() for d in 0:D-1 ]
  d_to_dface_to_sm_n = [ Int8[] for d in 0:D-1 ]
  d_to_dface_to_sm_nface = [ Int[] for d in 0:D-1 ]
  sm_face_to_ds = [ Int8[] for face in 1:num_faces ]
  sm_face_to_dfaces = [ Int[] for face in 1:num_faces ]
  vertex_to_bg_d = Int8[]
  vertex_to_bg_dface = Int[]
  IncrementalSurfaceMesh(
    vertex_coordinates,
    d_to_dface_to_vertices,
    d_to_dface_to_sm_n,
    d_to_dface_to_sm_nface,
    sm_face_to_ds,
    sm_face_to_dfaces,
    vertex_to_bg_d,
    vertex_to_bg_dface )
end

function IncrementalSurfaceMesh(sm::SurfaceMesh{D,T}) where {D,T}
  IncrementalSurfaceMesh{D,T}(num_faces(sm))
end

## Getters

function num_faces(sm::IncrementalSurfaceMesh,d::Integer)
  length( get_dface_to_vertices(sm,d) )
end

num_facets(sm::IncrementalSurfaceMesh{D}) where D = num_faces(sm,D-1)

num_vertices(sm::IncrementalSurfaceMesh) = num_faces(sm,0)

function get_background_face(sm::IncrementalSurfaceMesh,vertex::Integer)
  sm.vertex_to_bg_d[vertex], sm.vertex_to_bg_dface[vertex]
end

function get_dface_to_vertices(sm::IncrementalSurfaceMesh,d::Integer)
  sm.d_to_dface_to_vertices[d+1]
end

get_vertex_coordinates(sm::IncrementalSurfaceMesh) = sm.vertex_coordinates

function get_surface_mesh_face(ism::IncrementalSurfaceMesh,d::Integer,dface::Integer)
  sm_d = ism.d_to_dface_to_sm_n[d+1][dface]
  sm_dface = ism.d_to_dface_to_sm_nface[d+1][dface]
  sm_d, sm_dface
end

function get_faces_in_surface_mesh_face(
  ism::IncrementalSurfaceMesh,
  sm::SurfaceMesh,sm_d::Integer,sm_dface::Integer)

  sm_face = global_dface(sm,sm_d,sm_dface)
  ds = ism.sm_face_to_ds[sm_face]
  dfaces = ism.sm_face_to_dfaces[sm_face]
  ds, dfaces
end

function get_vertex_coordinates(sm::IncrementalSurfaceMesh,i::Integer)
  get_face_coordinates(sm,Val{0}(),i)
end

function get_edge_coordinates(sm::IncrementalSurfaceMesh,i::Integer)
  get_face_coordinates(sm,Val{1}(),i)
end

function get_facet_coordinates(sm::IncrementalSurfaceMesh{D},i::Integer) where D
  get_face_coordinates(sm,Val{D-1}(),i)
end

function get_face_coordinates(sm::IncrementalSurfaceMesh,::Val{d}, i::Integer) where d
  error("Not implemented for dimension = $d")
end

function get_face_coordinates(sm::IncrementalSurfaceMesh,::Val{0}, i::Integer)
  sm.vertex_coordinates[i]
end

function get_face_coordinates(sm::IncrementalSurfaceMesh,::Val{1}, i::Integer)
  df_to_v = get_dface_to_vertices(s,1)
  v = get_vertex_coordinates(s)
  Segment(v[df_to_v[i,1]],v[df_to_v[i,2]])
end

function get_face_coordinates(sm::IncrementalSurfaceMesh,::Val{2}, i::Integer)
  df_to_v = get_dface_to_vertices(sm,2)
  v = get_vertex_coordinates(sm)
  Triangle( v[df_to_v[i,1]], v[df_to_v[i,2]], v[df_to_v[i,3]] )
end

function facet_normal(ism::IncrementalSurfaceMesh,i::Integer)
  facet_coordinates = get_facet_coordinates(ism,i)
  facet_normal = normal(facet_coordinates)
  facet_normal / norm(facet_normal)
end

## Geometric queries

tolerance() = 1e-6 ## TODO: Add argument

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

## Printers

function writevtk(sm::IncrementalSurfaceMesh{D,T},file_base_name) where {D,T}
  d_to_vtk_type_id = Dict(0=>1,1=>3,2=>5,3=>10)
  num_points = num_vertices(sm)
  points = zeros(T,D,num_points)
  for (i ,p ) in enumerate( get_vertex_coordinates(sm) ), d in 1:D
    points[d,i] = p[d]
  end
  cells = MeshCell{Vector{Int64}}[]
  for d in 0:D-1
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

##

function cut_surface_mesh(sm::SurfaceMesh{D,T},m::CartesianMesh) where {D,T}
  vm = VolumeMesh( m )
  sm_face_to_bg_cells = compute_surface_mesh_face_to_cells(m,sm)
  ism = IncrementalSurfaceMesh(sm)
  
  bg_cells = Int[]
  nfaces = Int[]
  vertices = Int[]

  nfaces_cache = Int[]
  dfaces_cache = Int[]
  vertices_cache = Int[]
  lvertices_cache = Int[]

  for vertex in 1:num_vertices(sm)
    point = get_vertex_coordinates(sm,vertex)
    bg_d, bg_face = find_closest_background_face( vm, point, sm_face_to_bg_cells, vertex )
    point = closest_point( vm,bg_d,bg_face, sm,0,vertex )
    add_vertex!(ism,bg_d,bg_face,sm,0,vertex,point)
  end

  for sm_d in 1:D-1
    for sm_dface in 1:num_faces(sm,sm_d)
      n = sm_d-1
      cells_containg_nfaces!(bg_cells,vm,sm,sm_d,sm_dface,ism,n)
      while length(bg_cells) > 0
        bg_cell = popfirst!(bg_cells)
        fetch_boundary_nfaces!(ism,n,nfaces,vm,bg_cell,sm,sm_d,sm_dface)
        num_nfaces = length(nfaces)
        complete_boundary!(ism,n,nfaces,vertices,vm,bg_cell,sm,sm_d,sm_dface,vertices_cache,lvertices_cache)
        connect_boundary_faces!(vm,bg_cell,ism,n,nfaces,sm,sm_d,sm_dface, vertices_cache, lvertices_cache)
        add_faces!(ism,sm,sm_d,sm_dface,vm,bg_cell,vertices,dfaces_cache,nfaces_cache,vertices_cache )
        for i in num_nfaces+1:length(nfaces)
          nface = nfaces[i]
          bg_d, bg_dface = get_background_face(ism,n,nface,vm)
          dface_to_cells = get_faces(vm,bg_d,D)
          for lcell in 1:length(dface_to_cells,bg_dface)
            _bg_cell = dface_to_cells[bg_dface,lcell]
            if _bg_cell != bg_cell && _bg_cell ∉ bg_cells
              push!(bg_cells,_bg_cell)
            end
          end
        end
      end
    end
  end
  correct_face_orientation!(ism,sm)
  face_to_sm_face = get_face_to_surface_mesh_face(ism,sm)
  d_to_bg_df_to_face = background_face_to_surface_mesh_faces(ism,vm)
  sm = SurfaceMesh(ism.vertex_coordinates,ism.d_to_dface_to_vertices)
  sm, d_to_bg_df_to_face, face_to_sm_face
end

function add_vertex!(
  ism::IncrementalSurfaceMesh,
  bg_d::Integer,bg_face::Integer,
  sm::SurfaceMesh,sm_d::Integer,sm_face::Integer,
  point::Point)

  push!( ism.vertex_coordinates, point )
  push!( ism.vertex_to_bg_dface, bg_face )
  push!( ism.vertex_to_bg_d, bg_d )

  vertex = length( ism.vertex_coordinates )
  add_face!(ism,0,[vertex],sm,sm_d,sm_face)
  vertex
end

function add_face!(
  ism::IncrementalSurfaceMesh,d::Integer,vertices::Vector,
  sm::SurfaceMesh,sm_d::Integer,sm_dface::Integer)

  sm_face = global_dface(sm,sm_d,sm_dface)
  push!( ism.d_to_dface_to_vertices[d+1], vertices )
  push!( ism.d_to_dface_to_sm_n[d+1], sm_d )
  push!( ism.d_to_dface_to_sm_nface[d+1], sm_dface )
  push!( ism.sm_face_to_ds[sm_face],d)
  push!( ism.sm_face_to_dfaces[sm_face], num_faces(ism,d) )
  num_faces(ism,d) 
end

function find_closest_background_face(
  bg_mesh::VolumeMesh,
  point::Point,
  sm_face_to_bg_cells::Table,
  vertex::Integer)

  D = num_dims(bg_mesh)
  min_dist = tolerance()
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
  min_dist = tolerance()
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
  sm_df_to_sm_nf = get_faces(sm,sm_d,n)
  sm_nf_to_sm_v = get_dface_to_vertices(sm,n)
  for lnface in 1:length( sm_df_to_sm_nf, sm_d )
    sm_nface = sm_df_to_sm_nf[ sm_face, lnface ]
    for lvertex in 1:length( sm_nf_to_sm_v, sm_nface )
      sm_vertex = sm_nf_to_sm_v[ sm_nface, lvertex ]
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
  for sm_n in 0:sm_d
    dface_to_nfaces = get_faces(sm,sm_d,sm_n)
    for lnface in 1:length(dface_to_nfaces,sm_dface)
      sm_nface = dface_to_nfaces[sm_dface,lnface]
      ds, dfaces = get_faces_in_surface_mesh_face(ism,sm,sm_n,sm_nface)
      for i in 1:length(ds)
        if ds[i] == n &&
          is_nface_in_bg_cell(ism,n,dfaces[i],m,bg_cell) &&
          dfaces[i] ∉ nfaces
          
          if sm_n < sm_d || is_nface_in_bg_cell_boundary(ism,n,dfaces[i],m,bg_cell)
            push!(nfaces, dfaces[i])
          end
        end
      end
    end
  end
  nfaces
end

function fetch_boundary_vertices!(
  ism::IncrementalSurfaceMesh,vertices::Vector,
  m::VolumeMesh,bg_cell::Integer,
  sm::SurfaceMesh,sm_d::Integer,sm_dface::Integer)

  fetch_boundary_nfaces!(ism,0,vertices,m,bg_cell,sm,sm_d,sm_dface)
end

function complete_boundary!(
  ism::IncrementalSurfaceMesh,n::Integer,nfaces::Vector,vertices::Vector,
  m::VolumeMesh,bg_cell::Integer,
  sm::SurfaceMesh,sm_d::Integer,sm_dface::Integer,
  _vertices::Vector,lvertices::Vector)

  vertices = fetch_boundary_vertices!(ism,vertices,m,bg_cell,sm,sm_d,sm_dface)
  _vertices = fetch_boundary_vertices!(ism,_vertices,m,bg_cell,sm,sm_d,sm_dface)
  
  while length(vertices) > 0
    vertex = pop!(vertices)
    d,dface = get_background_face(ism,vertex)
    _d, _dface, point = _find_next_intersection(m,bg_cell,d,dface,sm,sm_d,sm_dface,ism,nfaces,n,_vertices)
    if _dface != UNSET
      _vertex = add_vertex!(ism,_d,_dface,sm,sm_d,sm_dface,point)
      push!( vertices, _vertex )
      push!( _vertices, _vertex )
      if n == 0
        nface = _vertex
      else
        @check n == 1
        resize!(lvertices,n+1)
        lvertices[1] = vertex
        lvertices[2] = _vertex
        nface = add_face!( ism, n, lvertices, sm, sm_d, sm_dface ) 
      end
      push!( nfaces, nface )
    end
  end
end

function connect_boundary_faces!(
  vm::VolumeMesh,bg_cell::Integer,
  ism::IncrementalSurfaceMesh,n::Integer,nfaces::Vector,
  sm::SurfaceMesh,sm_d::Integer,sm_dface::Integer,
  vertices::Vector,lvertices::Vector)

  D = num_dims(vm)
  n == D-2 || return
 
  cell_to_facets = get_faces(vm,D,D-1)
  resize!(vertices,0)
  resize!(lvertices,0)
  
  vertices = fetch_boundary_vertices!(ism,vertices,vm,bg_cell,sm,sm_d,sm_dface)
  for lfacet in 1:length(cell_to_facets,bg_cell) 
    bg_facet = cell_to_facets[bg_cell,lfacet]
    lvertices = get_vertices_in_background_facet(ism,vertices,vm,bg_facet,lvertices)
    if length(lvertices) == n+1
      if !are_vertices_in_any_dface(ism,lvertices,n,nfaces)
        nface = add_face!(ism,n,lvertices,sm,sm_d,sm_dface)
        push!(nfaces,nface)
      end
    end
  end
end

function add_faces!(
  ism::IncrementalSurfaceMesh,
  sm::SurfaceMesh,sm_d::Integer,sm_dface::Integer,
  vm::VolumeMesh,bg_cell::Integer,
  vertices::Vector,dfaces::Vector,nfaces::Vector,lvertices::Vector)
  
  vertices = fetch_boundary_vertices!(ism,vertices,vm,bg_cell,sm,sm_d,sm_dface)
  length(vertices) > 0 || return
  vertex = vertices[1]

  for d in 0:sm_d-1
    n = d+1
    dfaces = fetch_boundary_nfaces!(ism,d,dfaces,vm,bg_cell,sm,sm_d,sm_dface)
    nfaces = fetch_boundary_nfaces!(ism,n,nfaces,vm,bg_cell,sm,sm_d,sm_dface)
    dface_to_vertices = get_dface_to_vertices(ism,d)
    resize!( lvertices, n+1 )
    for dface in dfaces
      if !is_vertex_in_dface(ism,vertex,d,dface)
        lvertices[1] = vertex
        for i in 1:length(dface_to_vertices,dface)
          lvertices[i+1] = dface_to_vertices[dface,i]
        end
        if !are_vertices_in_any_dface(ism,lvertices,n,nfaces)
          add_face!(ism,n,lvertices,sm,sm_d,sm_dface)
        end
      end
    end
  end
end

function correct_face_orientation!(ism::IncrementalSurfaceMesh,sm::SurfaceMesh) 
  D = num_dims(sm)
  facet_to_vertices = get_dface_to_vertices(ism,D-1)
  for facet in 1:num_facets(ism)
    _facet_normal = facet_normal(ism,facet)
    sm_d, sm_dface = get_surface_mesh_face(ism::IncrementalSurfaceMesh,D-1,facet)
    @check sm_d == D-1
    sm_facet = sm_dface
    sm_facet_normal = get_facet_normal(sm,sm_facet)
    if sm_facet_normal ⋅ _facet_normal < 0
      n = length( facet_to_vertices, facet )
      vertex = facet_to_vertices[facet,n]
      facet_to_vertices[facet,n] = facet_to_vertices[facet,n-1]
      facet_to_vertices[facet,n-1] = vertex
    end
  end
  ism
end

function surface_mesh_face_to_background_face(ism::IncrementalSurfaceMesh,vm::VolumeMesh)
  sm_face_to_bg_d = Int[]
  sm_face_to_bg_dface = Int[]
  D = num_dims(vm)
  for d in 0:D-1
    for dface in 1:num_faces(ism,d)
      bg_d, bg_dface = get_background_face(ism,d,dface,vm)
      push!(sm_face_to_bg_d,bg_d)
      push!(sm_face_to_bg_dface,bg_face)
    end
  end
  sm_face_to_bg_d, sm_face_to_bg_dface
end

function background_face_to_surface_mesh_faces(ism::IncrementalSurfaceMesh,vm::VolumeMesh) 
  D = num_dims(vm)
  d_to_bg_dfaces = [ Int[] for d in 0:D ] 
  d_to_sm_faces = [ Int[] for d in 0:D ]
  d_to_bg_dface_to_sm_faces = Array{Table{Int}}(undef,D+1)
  sm_face = 0
  for sm_d in 0:D-1
    for sm_dface in 1:num_faces(ism,sm_d)
      sm_face += 1
      bg_d, bg_dface = get_background_face(ism,sm_d,sm_dface,vm)
      push!(d_to_bg_dfaces[bg_d+1],bg_dface)
      push!(d_to_sm_faces[bg_d+1],sm_face)
    end
  end
  for d in 0:D
    d_to_bg_dface_to_sm_faces[d+1] = Table(d_to_sm_faces[d+1],d_to_bg_dfaces[d+1],num_faces(vm,d))
  end
  d_to_bg_dface_to_sm_faces
end

function get_face_to_surface_mesh_face(ism::IncrementalSurfaceMesh,sm::SurfaceMesh)
  D = num_dims(sm)
  face_to_sm_face= zeros(Int, sum( [ num_faces(ism,d) for d in 0:D-1 ] ) )
  face = 0
  for d in 0:D-1, dface in 1:num_faces(ism,d)
    face += 1
    sm_d, sm_dface = get_surface_mesh_face(ism,d,dface)
    sm_face = global_dface(sm,sm_d,sm_dface)
    face_to_sm_face[face] = sm_face
  end
  face_to_sm_face 
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

function _find_next_intersection(
  m::VolumeMesh,bg_cell::Integer,bg_d::Integer,bg_face::Integer,
  sm::SurfaceMesh,sm_d::Integer,sm_dface::Integer,
  ism::IncrementalSurfaceMesh,nfaces::Vector,n::Integer,vertices::Vector)

  min_dist = tolerance()

  closest_bg_face = UNSET
  _d = UNSET
  D = num_dims(m)
  cell_to_facets = get_faces(m,D,D-1)
  for d in 0:D-sm_d
    facet_to_dfaces = get_faces(m,D-1,d)
    for lfacet in 1:length(cell_to_facets,bg_cell)
      bg_facet = cell_to_facets[bg_cell,lfacet]
      if sm_d == 1 || is_dface_in_facet(m,bg_d,bg_face,bg_facet) 
        for ldface in 1:length(facet_to_dfaces,bg_facet)
          _bg_dface = facet_to_dfaces[bg_facet,ldface]
          if !contains_any_nface(m,d,_bg_dface,ism,n,nfaces)
            if !is_any_vertex_on_background_face(ism,vertices,m ,d,_bg_dface)
              dist = distance(m,d,_bg_dface,sm,sm_d,sm_dface)
              if dist < min_dist
                min_dist = dist
                closest_bg_face = _bg_dface
              end
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


function is_any_vertex_on_background_face(
  ism::IncrementalSurfaceMesh,
  vertices::Vector,
  vm::VolumeMesh,bg_d::Integer,bg_dface::Integer)

  for vertex in vertices
    d,dface = get_background_face(ism,vertex)
    if is_nface_in_dface(vm,d,dface,bg_d,bg_dface)
      return true
    end
  end
  false
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

function get_background_face(ism::IncrementalSurfaceMesh,d::Integer,dface::Integer,vm::VolumeMesh)
  D = num_dims(vm)
  dface_to_vertices = get_dface_to_vertices(ism,d)
  max_bg_d = 0
  for lvertex in 1:length(dface_to_vertices,dface)
    vertex = dface_to_vertices[dface,lvertex]
    bg_d, bg_dface = get_background_face(ism,vertex)
    max_bg_d = max( max_bg_d, bg_d )
  end
  for bg_n in max_bg_d:D
    _bg_nface = UNSET
    vertex = dface_to_vertices[dface,1]
    bg_d, bg_dface = get_background_face(ism,vertex)
    bg_dface_to_nfaces = get_faces(vm,bg_d,bg_n)
    for lnface in 1:length(bg_dface_to_nfaces,bg_dface)
      bg_nface = bg_dface_to_nfaces[bg_dface,lnface]
      found = true
      for lvertex in 1:length(dface_to_vertices,dface)
        vertex = dface_to_vertices[dface,lvertex]
        _bg_d, _bg_dface = get_background_face(ism,vertex)
        if !is_nface_in_dface(vm,_bg_d,_bg_dface,bg_n,bg_nface)
          found = false
          break
        end
      end
      if found
        if _bg_nface == UNSET
          _bg_nface = bg_nface
        else
          _bg_nface = UNSET
          break
        end
      end
    end
    if _bg_nface != UNSET
      return bg_n, _bg_nface
    end
  end
  @assert false
end

## Queries

function is_nface_in_bg_cell(
  ism::IncrementalSurfaceMesh,n::Integer,nface::Integer,
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

function is_nface_in_bg_cell_boundary(
  ism::IncrementalSurfaceMesh,n::Integer,nface::Integer,
  m::VolumeMesh,bg_cell::Integer)

  D = num_dims(m)
  nface_to_vertices = get_dface_to_vertices(ism,n)
  for lvertex in 1:length(nface_to_vertices,nface)
    vertex = nface_to_vertices[nface,lvertex]
    bg_d, bg_dface = get_background_face(ism,vertex)
    if bg_d == D && !is_dface_in_cell(m,bg_d,bg_dface,bg_cell)
      return false
    end
  end
  true
end

function is_nface_in_surface_mesh_dface(
  ism::IncrementalSurfaceMesh,n::Integer,nface::Integer,
  sm::SurfaceMesh,sm_d::Integer,sm_dface::Integer)

  _d, _dface = get_surface_mesh_face(ism,n,nface)
  return is_nface_in_dface(sm,_d,_dface,sm_d,sm_dface)
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

function is_vertex_in_background_facet(
  ism::IncrementalSurfaceMesh,vertex::Integer,
  vm::VolumeMesh,bg_facet::Integer)

  is_dface_in_background_facet(ism,0,vertex,vm,bg_facet)
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
