
const FACE_UNDEF = 0
const FACE_IN = 1
const FACE_OUT = 2
const FACE_BOUNDARY = 3
const FACE_CUT = -1

struct CellMeshCache
  d_to_num_dfaces::Vector{Int}
  d_to_dface_to_new_dfaces::Vector{Vector{Vector{Int}}}
  d_to_dface_to_cell_dface::Vector{Vector{Int}}
  cell_to_lfacet_to_orientation::Table{Int}
  vertex_to_surface_mesh_faces::Vector{Vector{Int}}
  facet_to_surface_mesh_facet::Vector{Int}
end

struct CutterCache
  d_to_ldface_to_dface::Table{Int}
  m_n_to_mface_to_nfaces::Matrix{Table{Int}}
  cell_to_lfacet_to_orientation::Table{Int}
end

struct MeshCache
  cell::CellMeshCache
  cutter::CutterCache
  vector_cache_1::Vector{Int}
  vector_cache_2::Vector{Int}
  vector_cache_3::Vector{Int}
end

struct CellMesh{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  m_n_to_mface_to_nfaces::Matrix{Table{Int}}
  d_to_dface_to_in_out_boundary::Table{Int}
  cache::MeshCache
end

function CellMesh(box::BoundingBox{D,T}) where {D,T}
  mesh = CellMesh{D,T}()
  reset!(mesh,box)
end

function CellMesh{D,T}() where {D,T}
  coordinates = Point{D,T}[]
  m_n_to_mf_to_nf = Matrix{Table{Int}}(undef,D+1,D+1)
  d_df_to_io = Table{Int}()
  for d in 0:D, n in 0:D
    m_n_to_mf_to_nf[d+1,n+1] = Table{Int}()
  end
  cache = MeshCache(D)
  CellMesh{D,T}( coordinates, m_n_to_mf_to_nf, d_df_to_io, cache )
end

function reset!(m::CellMesh{D,T},box::BoundingBox{D,T}) where {D,T} 
  table_m_n_to_mface_to_nfaces = D_to_m_n_to_mface_to_nfaces_for_hexD_to_tetD[D]
  resize!(m.vertex_coordinates,num_vertices(box))
  for (i,v) in enumerate( get_vertices(box) )
    m.vertex_coordinates[i] = v
  end
  for d in 0:D, n in 0:d
    df_to_nf = table_m_n_to_mface_to_nfaces[d+1,n+1]
    resize!( get_faces(m,d,n), 0 )
    append!( get_faces(m,d,n), df_to_nf )
  end
  resize!(m.d_to_dface_to_in_out_boundary, 0 )

  reset_d_to_dface_to_in_out_boundary!(m,m.cache)
  reset_cache!(m.cache,m)
  m
end

function add_vertex!(mesh::CellMesh,d::Integer,face::Integer,point::Point,sm_face::Integer)
  if face != UNSET
    add_surface_mesh_face_to_vertex!(mesh.cache,d,face,sm_face)
    _add_vertex!(mesh,d,face,point)
  else
    UNSET
  end
end

function compact!(mesh::CellMesh)
  compact_cache!(mesh.cache,mesh)
  compact_mesh!(mesh)
  reset_d_to_dface_to_in_out_boundary!(mesh,mesh.cache)
  compute_facet_to_cells!(mesh)
  mesh
end

function compute_in_out!(mesh::CellMesh,sm::SurfaceMesh)
  cache = mesh.cache
  compute_facet_to_surface_mesh_facet!(mesh,cache,sm)
  define_cells_touching_boundary!(mesh,cache,sm)
  define_boundary_faces!(mesh,cache)

  queue = get_vector!(cache)
  resize!(queue,num_cells(mesh))

  D = num_dims(mesh)
  c_to_f = get_faces(mesh,D,D-1)
  f_to_c = get_faces(mesh,D-1,D)
  for cell in 1:num_cells(mesh)
    if is_cell_defined(mesh,cell)
      cell_io = get_cell_in_out(mesh,cell)
      head = 1
      tail = 1
      queue[head] = cell
      while head ≤ tail
        head_cell = queue[head]
        head += 1
        for lfacet in 1:length(c_to_f,head_cell)
          facet = c_to_f[head_cell,lfacet]
          if !is_face_boundary(mesh,D-1,facet)
            for i in 1:length(f_to_c,facet)
              cell_around = f_to_c[facet,i]
              if !is_cell_defined(mesh,cell_around)
                set_cell_in_out!(mesh,cell_around,cell_io)
                tail += 1
                queue[tail] = cell_around
              end
            end
          end
        end
      end
    end
  end

  define_faces!(mesh)
  fix_boundary_facets!(mesh)
  mesh
end

function compute_cell_mesh!(
  mesh::CellMesh,
  sm::SurfaceMesh,
  vm::VolumeMesh,cell::Integer,
  d_to_bg_dface_to_sm_faces::Vector{<:Table})

  D = num_dims(vm)

  for sm_d in 0:D-1
    for d in sm_d:D
      cell_to_dfaces = get_faces(vm,D,d)
      dface_to_sm_faces = d_to_bg_dface_to_sm_faces[d+1]
      for ldface in 1:length(cell_to_dfaces,cell)
        dface = cell_to_dfaces[cell,ldface]
        for i in 1:length(dface_to_sm_faces,dface)
          sm_face = dface_to_sm_faces[dface,i]
          if face_dimension(sm,sm_face) == sm_d 
            cut_cell_mesh!(mesh,sm,sm_face,d,ldface)
          end
        end
      end
    end
  end
  compact!(mesh)
  compute_in_out!(mesh,sm)
end

# Intersections

tolerance(m::CellMesh) = 1e-7

function cut_cell_mesh!(mesh::CellMesh,sm::SurfaceMesh,sm_face::Integer,d::Integer,ldface::Integer)
  sm_d = face_dimension(sm,sm_face)
  if sm_d == 0
    point = get_vertex_coordinates(sm,sm_face)
    _d, iface = find_closest_face(mesh,point,d,ldface)
    vertex = add_vertex!(mesh,_d,iface,point,sm_face)
  else
    n = sm_d-1
    df_to_nf = get_faces(mesh,sm_d,sm_d-1)
    nf_to_v = get_dface_to_vertices(mesh,sm_d-1)

    nfaces = find_dfaces_capturing_surface_mesh_face!(mesh,sm_d-1,sm,sm_face)

#    if length(nfaces) == 0
#      _d, iface, point = find_intersection_point_on_boundary_face(mesh,sm,sm_face)
#      vertex = add_vertex!(mesh,_d,iface,point,sm_face)
#      kface = vertex
#      for k in 0:n-1
#        kface != UNSET || break
#        _d, iface, point = find_next_boundary_point(mesh,k,kface,sm,sm_face)
#        vertex = add_vertex!(mesh,_d,iface,point,sm_face)
#        kface = find_container_dface(mesh,k+1,k,kface,vertex)
#      end
#      if kface != UNSET
#        push!(nfaces,kface)
#      end
#    end

    while length(nfaces) > 0
      nface = popfirst!(nfaces)
      _d, iface, point = find_next_point(mesh,nface,sm,sm_face,d,ldface) 
      vertex = add_vertex!(mesh,_d,iface,point,sm_face)
      if vertex != UNSET
        for _nface in reverse(1:num_dfaces(mesh,n))
          if isactive(mesh,n,_nface)
            if _nface != nface && 
              _is_vertex_in_face(mesh,vertex,n,_nface)  && 
              _is_dface_capturing_surface_mesh_face(mesh,n,_nface,sm,sm_face)

              push!(nfaces,_nface)
            end
          end
        end
      end
    end
  end
end

function _is_vertex_in_face(mesh,vertex,d,dface)
  df_to_v = get_dface_to_vertices(mesh,d)
  for lvertex in 1:length(df_to_v,dface)
    if vertex == df_to_v[dface,lvertex]
      return true
    end
  end
  false
end

function find_closest_face(mesh::CellMesh,point::Point,d::Integer,ldface::Integer)
  cell_cache = mesh.cache.cell
  min_distance = tolerance(mesh)
  closest_face = UNSET
  D = num_dims(mesh)
  for n in 0:d
    dface_to_nfaces = get_faces(mesh,d,n)
    for dface in 1:num_dfaces(mesh,d)
      if isactive(mesh,d,dface)
        if d == D || get_cell_dface(cell_cache,d,dface) == ldface
          for lnface in 1:length(dface_to_nfaces,dface)
            nface = dface_to_nfaces[dface,lnface]
            dist = distance(mesh,n,nface,point)
            if dist ≤ min_distance
              min_distance = dist
              closest_face = nface
            end
          end
        end
      end
    end
    if closest_face != UNSET
      return (n,closest_face)
    end
  end
  @check false
  (UNSET,UNSET)
end

function find_closest_face(mesh::CellMesh,point::Point)
  min_distance = tolerance(mesh)
  iface = UNSET
  for d in 0:num_dims(mesh)
    for i in 1:num_dfaces(mesh,d)
      if isactive(mesh,d,i)
        dist = distance(mesh,d,i,point)
        if dist ≤ min_distance
          min_distance = dist
          iface = i
        end
      end
    end
    if iface != UNSET
      return (d,iface)
    end
  end
  (UNSET,UNSET)
end

function find_dfaces_capturing_surface_mesh_face!(mesh::CellMesh,d::Integer,sm::SurfaceMesh,sm_face::Integer)
  dfaces = get_vector_cache!(mesh)
  df_to_v = get_dface_to_vertices(mesh,d)
  for dface in 1:num_dfaces(mesh,d)
    if isactive(mesh,d,dface) && _is_dface_capturing_surface_mesh_face(mesh,d,dface,sm,sm_face)
      push!(dfaces,dface)
    end
  end
  dfaces
end

function find_intersection_point_on_boundary_face(mesh::CellMesh,sm::SurfaceMesh,sm_face::Integer)
  D = num_dims(mesh)
  sm_dim = face_dimension(sm,sm_face)
  for d in 0:D-sm_dim
    for n in 0:d
      df_to_nf = get_faces(mesh,d,n)
      iface = UNSET
      min_distance = tolerance(mesh)
      for dface in 1:num_dfaces(mesh,d)
        if isactive(mesh,d,dface)
          if is_dface_in_cell_mesh_boundary(mesh,d,dface)
            for lnface in 1:length(df_to_nf,dface)
              nface = df_to_nf[dface,lnface]
              dist = distance(mesh,n,nface,sm,sm_face) 
              if dist < min_distance
                min_distance = dist
                iface = nface
              end
            end
          end
        end
      end
      if iface != UNSET
        point = closest_point(mesh,n,iface,sm,sm_face)
        return (n,iface,point)
      end
    end
  end
  point = zero( eltype(get_vertex_coordinates(mesh)) )
  (UNSET,UNSET,point)
end


function find_next_point(
  mesh::CellMesh,dface::Integer,
  sm::SurfaceMesh,sm_face::Integer,
  cell_d::Integer,cell_dface::Integer)

  cell_cache = mesh.cache.cell
  D = num_dims(mesh)
  sm_d = face_dimension(sm,sm_face)
  d = sm_d-1
  dim = cell_d-d-1
  iface = UNSET
  min_distance = tolerance(mesh)
  
  for face in 1:num_dfaces(mesh,cell_d)
    if isactive(mesh,cell_d,face)
      if cell_d == D || get_cell_dface(cell_cache,cell_d,face) == cell_dface
        if _is_dface_in_nface_boundary(mesh,cell_d,face,d,dface)
          opposite_face = _get_opposite_face(mesh,cell_d,face,d,dface)
          if !_is_touching_surface_mesh_face(mesh,dim,opposite_face,sm,sm_face)
            dist = distance(mesh,dim,opposite_face,sm,sm_face) 
            if dist < min_distance
              min_distance = dist
              iface = opposite_face
            end
          end
        end
      end
    end
  end
  if iface != UNSET
    opposite_face = iface
    iface = UNSET
    min_distance = tolerance(mesh)
    for n in 0:dim
      df_to_nf = get_faces(mesh,dim,n)
      for lnface in 1:length(df_to_nf,opposite_face)
        nface = df_to_nf[opposite_face,lnface]
        dist = distance(mesh,n,nface,sm,sm_face)
        if dist < min_distance
          min_distance = dist
          iface = nface
        end
      end
      if iface != UNSET
        point = closest_point(mesh,n,iface,sm,sm_face)
        return (n,iface,point)
      end
    end
  end
  point = zero( eltype(get_vertex_coordinates(mesh)) )
  (UNSET,UNSET,point)
end

function find_next_boundary_point(mesh::CellMesh,d::Integer,dface::Integer,sm::SurfaceMesh,sm_face::Integer)
  D = num_dims(mesh)
  sm_d = face_dimension(sm,sm_face)
  for k in d+1:D-1
    for n in d+1:k
      dim = n-d-1
      kf_to_nf = get_faces(mesh,k,n)
      iface = UNSET
      min_distance = tolerance(mesh)
      for kface in 1:num_dfaces(mesh,k)
        if isactive(mesh,k,kface)
          if is_dface_in_cell_mesh_boundary(mesh,k,kface)
            for lnface in 1:length(kf_to_nf,kface)
              nface = kf_to_nf[kface,lnface]
              if _is_dface_in_nface_boundary(mesh,n,nface,d,dface)
                opposite_face = _get_opposite_face(mesh,n,nface,d,dface)
                if !_is_touching_surface_mesh_face(mesh,dim,opposite_face,sm,sm_face)
                  dist = distance(mesh,dim,opposite_face,sm,sm_face)
                  if dist < min_distance
                    min_distance = dist
                    iface = opposite_face
                  end
                end
              end
            end
          end
        end
      end
      if iface != UNSET
        point = closest_point(mesh,dim,iface,sm,sm_face)
        return (dim,iface,point)
      end
    end
  end
  point = zero(get_vertex_coordinates(mesh,1))
  (UNSET,UNSET,point)
end

function find_container_dface(mesh::CellMesh,d::Integer,n::Integer,nface::Integer,vertex::Integer)
  if vertex == UNSET
    return UNSET
  end
  @check d == n+1
  df_to_nf = get_faces(mesh,d,n)
  df_to_v = get_dface_to_vertices(mesh,d)
  for dface in reverse(1:num_dfaces(mesh,d))
    if isactive(mesh,d,dface)
      found = false
      for lnface in 1:length(df_to_nf,dface)
        if nface == df_to_nf[dface,lnface]
          found = true
          break
        end
      end
      if found
        for lvertex in 1:length(df_to_v,dface)
          if vertex == df_to_v[dface,lvertex]
            return dface
          end
        end
      end
    end
  end
  @check false
  UNSET
end

function _is_dface_in_nface_boundary(mesh::CellMesh,n::Integer,nface::Integer,d::Integer,dface::Integer)
  @check n > d
  nf_to_df = get_faces(mesh,n,d)
  for ldface in 1:length(nf_to_df,nface)
    if nf_to_df[nface,ldface] == dface
      return true
    end
  end
  false
end

function _get_opposite_face(mesh::CellMesh,d::Integer,dface::Integer,n::Integer,nface::Integer)
  opposite_face_dim = d-n-1
  df_to_of = get_faces(mesh,d,opposite_face_dim)
  of_to_v = get_dface_to_vertices(mesh,opposite_face_dim)
  nf_to_v = get_dface_to_vertices(mesh,n)

  for loface in 1:length(df_to_of,dface)
    opposite = true
    oface = df_to_of[dface,loface]
    for oface_lvertex in 1:length(of_to_v,oface), nface_lvertex in 1:length(nf_to_v,nface)
      if of_to_v[oface,oface_lvertex] == nf_to_v[nface,nface_lvertex]
        opposite = false
        break
      end
    end
    if opposite
      return oface
    end
  end
  @check false
  return UNSET
end

function _is_dface_capturing_surface_mesh_face(
  mesh::CellMesh,d::Integer,dface::Integer,
  sm::SurfaceMesh,sm_face::Integer)
  df_to_v = get_dface_to_vertices(mesh,d)
  for lvertex in 1:length(df_to_v,dface)
    vertex = df_to_v[dface,lvertex]
    if !_is_surface_mesh_face_on_vertex(mesh,vertex,sm,sm_face)
      return false
    end
  end
  true
end

function _is_touching_surface_mesh_face(
  mesh::CellMesh,d::Integer,dface::Integer,
  sm::SurfaceMesh,sm_face::Integer)
  df_to_v = get_dface_to_vertices(mesh,d)
  for lvertex in 1:length(df_to_v,dface)
    vertex = df_to_v[dface,lvertex]
    if _is_surface_mesh_face_on_vertex(mesh,vertex,sm,sm_face)
      return true
    end
  end
  false
end

function _is_surface_mesh_face_on_vertex(mesh::CellMesh,vertex::Integer,sm::SurfaceMesh,sm_face::Integer)
  d = face_dimension(sm,sm_face)
  sm_dface = local_dface(sm,sm_face,d)
  for n in 0:d
    df_to_nf = get_faces(sm,d,n)
    for sm_f in get_surface_mesh_faces_on_vertex(mesh,vertex)
      for sm_lnface in 1:length(df_to_nf,sm_dface)
        if sm_f == global_dface(sm,n,df_to_nf[sm_dface,sm_lnface])
          return true  
        end
      end
    end
  end
  false
end

function get_surface_mesh_faces_on_vertex(mesh::CellMesh,vertex::Integer)
  cache = get_cache(mesh)
  cache.cell.vertex_to_surface_mesh_faces[vertex]
end

# Getters

num_dims(mesh::CellMesh{D}) where D = D

function num_vertices(m::CellMesh) 
  length(m.vertex_coordinates)
end

function num_dfaces(m::CellMesh,d::Int)
  if d == 0
    num_vertices(m)
  else
    length( get_faces(m,d,0) ) 
  end
end

num_edges(m::CellMesh) = num_dfaces(m,1)

num_facets(m::CellMesh{D}) where D = num_dfaces(m,D-1)

num_cells(m::CellMesh{D}) where D = num_dfaces(m,D)

function get_faces(m::CellMesh,d::Int,n::Int)
  m.m_n_to_mface_to_nfaces[d+1,n+1]
end

function get_dface_to_vertices(m::CellMesh,d::Int)
  get_faces(m,d,0)
end

function get_cell_to_vertices(m::CellMesh{D}) where D
  get_faces(m,D,0)
end

function get_cell_to_facets(m::CellMesh{D}) where D
  get_faces(m,D,D-1)
end

function get_facet_to_vertices(m::CellMesh{D}) where D
  get_faces(m,D-1,0)
end

function get_vertex_coordinates(m::CellMesh) 
  m.vertex_coordinates
end

function get_face_coordinates(m::CellMesh,::Val{0},i::Integer)
  v = get_vertex_coordinates(m)
  v[i]
end

function get_face_coordinates(m::CellMesh,::Val{1},i::Integer)
  df_to_v = get_dface_to_vertices(m,1)
  v = get_vertex_coordinates(m)
  Segment( v[ df_to_v[i,1] ], v[ df_to_v[i,2] ] )
end

function get_face_coordinates(m::CellMesh,::Val{2},i::Integer)
  df_to_v = get_dface_to_vertices(m,2)
  v = get_vertex_coordinates(m)
  Triangle( v[ df_to_v[i,1] ], v[ df_to_v[i,2] ], v[ df_to_v[i,3] ] )
end

function get_face_coordinates(m::CellMesh,::Val{3},i::Integer)
  df_to_v = get_dface_to_vertices(m,3)
  v = get_vertex_coordinates(m)
  Tetrahedron( v[ df_to_v[i,1]], v[ df_to_v[i,2] ], v[ df_to_v[i,3] ], v[ df_to_v[i,4] ] )
end

function get_face_coordinates(m::CellMesh,::Val{d},i::Integer) where d
  throw(ArgumentError("get_face(::CellMesh,::Val($d),::Integer) not implemented"))
end

function isactive(m::CellMesh,d::Integer,i::Integer)
  d != 0 || return true
  isactive(get_dface_to_vertices(m,d),i)
end

function get_vertex_coordinates(m::CellMesh,i::Integer)
  get_face_coordinates(m,Val{0}(),i)
end

function get_edge_coordinates(m::CellMesh,i::Integer)
  get_face_coordinates(m,Val{1}(),i)
end

function get_facet_coordinates(m::CellMesh{D},i::Integer) where D
  get_face_coordinates(m,Val{D-1}(),i)
end

function get_cell_coordinates(m::CellMesh{D},i::Integer) where D
  get_face_coordinates(m,Val{D}(),i)
end

function facet_normal(m::CellMesh{D},i::Integer) where D
  f = get_facet_coordinates(m,i)
  normal(f)
end

function face_center(m::CellMesh,::Val{d},i::Integer) where d
  f = get_face_coordinates(m,Val{d}(),i)
  center(f)
end

function facet_center(m::CellMesh{D},i::Integer) where D
  face_center(m,Val{D-1}(),i)
end

function get_vertex_in_out_boundary(mesh::CellMesh,cell::Integer)
  get_face_in_out_boundary(mesh,0,cell)
end

function get_cell_in_out(mesh::CellMesh{D},cell::Integer) where D
  get_face_in_out_boundary(mesh,D,cell)
end

function is_face_interior(mesh::CellMesh,d::Integer,face::Integer)
  get_face_in_out_boundary(mesh,d,face) == FACE_IN
end

function is_face_exterior(mesh::CellMesh,d::Integer,face::Integer)
  get_face_in_out_boundary(mesh,d,face) == FACE_OUT
end

function is_face_boundary(mesh::CellMesh,d::Integer,face::Integer)
  get_face_in_out_boundary(mesh,d,face) == FACE_BOUNDARY
end

function is_face_defined(mesh::CellMesh,d::Integer,face::Integer)
  get_face_in_out_boundary(mesh,d,face) != FACE_UNDEF
end

function is_vertex_interior(mesh::CellMesh,vertex::Integer)
  is_face_interior(mesh,0,vertex)
end

function is_vertex_exterior(mesh::CellMesh,vertex::Integer)
  is_face_exterior(mesh,0,vertex)
end

function is_vertex_boundary(mesh::CellMesh,vertex::Integer)
  is_face_boundary(mesh,0,vertex)
end

function is_facet_boundary(mesh::CellMesh{D},facet::Integer) where D
  is_face_boundary(mesh,D-1,facet)
end

function is_cell_interior(mesh::CellMesh{D},cell::Integer) where D
  is_face_interior(mesh,D,cell)
end

function is_cell_exterior(mesh::CellMesh{D},cell::Integer) where D
  is_face_exterior(mesh,D,cell)
end

function is_cell_defined(mesh::CellMesh{D},cell::Integer) where D
  is_face_defined(mesh,D,cell)
end

function are_all_faces_defined(mesh::CellMesh)
  for d in 0:num_dims(mesh), dface in 1:num_dfaces(mesh,d)
    if !is_face_defined(mesh,d,dface)
      return false
    end
  end
  true
end

function is_surface_mesh_captured(mesh::CellMesh)
  for i in 1:num_facets(mesh)
    if is_facet_boundary(mesh,i)
      return true
    end
  end
  false
end

# Setters

function set_face_in_out_boundary!(mesh::CellMesh,d::Integer,face::Integer,val)
  mesh.d_to_dface_to_in_out_boundary[d+1,face] = val
end

function get_face_in_out_boundary(mesh::CellMesh,d::Integer,face::Integer)
  mesh.d_to_dface_to_in_out_boundary[d+1,face]
end

function set_face_as_interior!(mesh::CellMesh,d::Integer,face::Integer)
  set_face_in_out_boundary!(mesh,d,face,FACE_IN)
end

function set_face_as_exterior!(mesh::CellMesh,d::Integer,face::Integer)
  set_face_in_out_boundary!(mesh,d,face,FACE_OUT)
end

function set_face_as_boundary!(mesh::CellMesh,d::Integer,face::Integer)
  set_face_in_out_boundary!(mesh,d,face,FACE_BOUNDARY)
end

function set_cell_in_out!(mesh::CellMesh{D},cell::Integer,val) where D
  @check val != FACE_BOUNDARY
  set_face_in_out_boundary!(mesh,D,cell,val)
end

function set_cell_as_interior!(mesh::CellMesh{D},cell::Integer) where D
  set_face_as_interior!(mesh,D,cell)
end

function set_cell_as_exterior!(mesh::CellMesh{D},cell::Integer) where D
  set_face_as_exterior!(mesh,D,cell)
end

function set_cell_as_boundary!(mesh::CellMesh{D},cell::Integer) where D
  set_face_as_boundary!(mesh,D,cell)
end

function set_facet_as_interior!(mesh::CellMesh{D},facet::Integer) where D
  set_face_as_interior!(mesh,D-1,facet)
end

function set_facet_as_exterior!(mesh::CellMesh{D},facet::Integer) where D
  set_face_as_exterior!(mesh,D-1,facet)
end

function set_facet_as_boundary!(mesh::CellMesh{D},facet::Integer) where D
  set_face_as_boundary!(mesh,D-1,facet)
end

# Geometrical queries

@generated function distance(m::CellMesh{D},d::Int,i::Int,p::Point{D}) where D
  body = ""
  for d in 0:D
    body *= "if d == $d \n"
    body *= "  f$d = get_face_coordinates(m,Val{$d}(),i) \n"
    body *= "  distance(p,f$d) \n"
    body *= "else"
  end
  error = "  throw(ArgumentError(\"\$d-face does not exist\"))\n"
  str = body * '\n' * error * "end"
  Meta.parse(str)
end

function distance(m::CellMesh,d::Integer,face::Integer,sm::SurfaceMesh,sm_face::Integer)
  sm_d = face_dimension(sm,sm_face)
  sm_dface = local_dface(sm,sm_face,sm_d)
  distance(m,d,face,sm,sm_d,sm_dface)
end

@generated function distance(
  m::CellMesh{D},d::Integer,dface::Integer,
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

# closest_point_in_mesh_to_surface_mesh
function closest_point(m::CellMesh,d::Integer,face::Integer,sm::SurfaceMesh,sm_face::Integer)
  sm_d = face_dimension(sm,sm_face)
  sm_dface = local_dface(sm,sm_face,sm_d)
  closest_point(m,d,face,sm,sm_d,sm_dface)
end

@generated function closest_point(
  m::CellMesh{D},d::Integer,dface::Integer,
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

# Printers

function writevtk(m::CellMesh{D,T},file_base_name) where {D,T}
  d_to_vtk_type_id = Dict(0=>1,1=>3,2=>5,3=>10)
  num_points = num_vertices(m)
  points = zeros(T,D,num_points+num_dfaces(m,D-1))
  for (i ,p ) in enumerate(get_vertex_coordinates(m)), d in 1:D
    points[d,i] = p[d]
  end

  # TODO: normal to cell, not to point
  for i in 1:num_facets(m)
    center = facet_center(m,i)
    for d in 1:D
      points[d,i+num_points] = center[d]
    end
  end

  cells = MeshCell{Vector{Int64}}[]
  for d in 0:num_dims(m)
    dface_to_vertices = get_dface_to_vertices(m,d)
    vtk_type = VTKCellType(d_to_vtk_type_id[d])
    for i in 1:num_dfaces(m,d)
      vertices = [ dface_to_vertices[i,j] for j in 1:vtk_type.nodes ]
      push!( cells, MeshCell(vtk_type,vertices) )
    end
  end

  normals = zeros(3,num_points+num_facets(m))
  for i in 1:num_facets(m)
    normal = facet_normal(m,i)
    normal = normal / norm(normal)
    for d in 1:D
      normals[d,i+num_points] = normal[d]
    end
  end

  vtkfile = vtk_grid(file_base_name,points,cells)
  vtkfile["facet_normals",VTKPointData()] = normals
  vtkfile["IO",VTKCellData()] = _get_in_out_face_data(m,0:num_dims(m))

  vtk_save(vtkfile)
end

# Helpers

function CellMeshCache(D::Integer)
  d_num_df = [ 0 for d in 0:D ]
  d_to_df_to_new_df = [ Vector{Int}[] for d in 1:D ]
  d_to_df_to_cdf = [ Int[] for d in 0:D-1 ]
  c_lf_to_o = Table{Int}()
  v_to_smf = Vector{Int}[]
  f_to_sm_f = Int[]
  CellMeshCache( d_num_df, d_to_df_to_new_df, d_to_df_to_cdf, c_lf_to_o, v_to_smf, f_to_sm_f )
end

function CutterCache(D::Integer)
  d_to_ldf_to_df = Table{Int}()
  m_n_to_mf_to_nf = Matrix(undef,D+1,D+1)
  c_ldf_to_o = Table{Int}()
  for i in 1:D+1, j in 1:D+1
    m_n_to_mf_to_nf[i,j] = Table{Int}()
  end
  CutterCache( d_to_ldf_to_df, m_n_to_mf_to_nf, c_ldf_to_o )
end

function MeshCache(D::Integer)
  cell_cache = CellMeshCache(D) 
  cutter_cache = CutterCache(D)
  vector1 = []
  vector2 = []
  vector3 = []
  MeshCache( cell_cache, cutter_cache, vector1, vector2, vector3 )
end

function reset_cell_cache!(cache::CellMeshCache,mesh::CellMesh)
  D = num_dims(mesh)
  for d in 0:D
    set_num_dfaces!(cache,d, num_dfaces(mesh,d) )
  end
  table_c_lf_to_o = D_to_cell_lfacet_to_orientation_for_hexD_to_tetD[D] 
  table_d_to_df_to_cdf = D_to_d_to_subdface_to_dface_for_hexD_to_tetD[D]
  for d in 1:D
    df_to_new_df = get_dface_to_new_dfaces(cache,d)
    if length(df_to_new_df) < num_dfaces(cache,d)
      resize!( df_to_new_df, num_dfaces(cache,d) )
    end
    for i in 1:num_dfaces(mesh,d)
      if !isassigned(df_to_new_df,i)
        df_to_new_df[i] = []
      end
      resize!( df_to_new_df[i],0)
    end
  end
  resize!( cache.cell_to_lfacet_to_orientation, 0 )
  append!( cache.cell_to_lfacet_to_orientation, table_c_lf_to_o )
  v_to_smf = cache.vertex_to_surface_mesh_faces
  for i in 1:num_vertices(mesh)
    if length(v_to_smf) < i
      push!(v_to_smf,[])
    else
      resize!(v_to_smf[i],0) 
    end
  end
  for d in 0:D-1
    df_to_cdf = table_d_to_df_to_cdf[d+1]
    resize!( cache.d_to_dface_to_cell_dface[d+1], 0 )
    append!( cache.d_to_dface_to_cell_dface[d+1], df_to_cdf )
  end
  cache
end

function reset_cutter_cache!(cutter::CutterCache)
  resize!(cutter.d_to_ldface_to_dface,0)
  resize!(cutter.cell_to_lfacet_to_orientation,0)
  for t in cutter.m_n_to_mface_to_nfaces
    resize!(t,0)
  end
  cutter
end

function reset_cache!(cache::MeshCache,mesh::CellMesh)
  reset_cell_cache!(cache.cell,mesh)
  reset_cutter_cache!(cache.cutter)
  cache
end

get_cache(mesh::CellMesh) = mesh.cache

function get_vector_cache!(mesh::CellMesh) 
  cache = get_cache(mesh)
  resize!(cache.vector_cache_1,0)
end

get_vector!(a::MeshCache) = resize!(a.vector_cache_2,0)

get_vector_bis!(a::MeshCache) = resize!(a.vector_cache_3,0)

num_dims(cache::CellMeshCache) = length(cache.d_to_dface_to_new_dfaces)

num_dfaces(cache::CellMeshCache,d::Integer) = cache.d_to_num_dfaces[d+1]

set_num_dfaces!(cache::CellMeshCache,d::Integer,val) = cache.d_to_num_dfaces[d+1] = val

get_orientation(a::CellMeshCache,cell::Integer,lfacet::Integer) = a.cell_to_lfacet_to_orientation[cell,lfacet]

function  set_orientation!(a::CellMeshCache,cell::Integer,lfacet::Integer,val)
  a.cell_to_lfacet_to_orientation[cell,lfacet] = val
end

get_dface_to_new_dfaces(a::CellMeshCache,d::Integer) = a.d_to_dface_to_new_dfaces[d]

get_new_faces(a::CellMeshCache,d::Integer,i::Integer) = a.d_to_dface_to_new_dfaces[d][i]

get_cell_dface(a::CellMeshCache,d::Integer,i::Integer) = a.d_to_dface_to_cell_dface[d+1][i]

function is_dface_in_cell_mesh_boundary(mesh::CellMesh,d::Integer,dface::Integer)
  cache = get_cache(mesh)
  cell_cache = cache.cell
  get_cell_dface(cache.cell,d,dface) != 0
end

function _add_vertex!(m::CellMesh{D,T},p::Point{D,T}) where {D,T}
  v = get_vertex_coordinates(m)
  push!( v, p )
  num_vertices(m)
end

function remove_dface!(m::CellMesh,d::Integer,i::Int)
  for n in 0:d
    df_to_nf = get_faces(m,d,n)
    remove!(df_to_nf,i)
  end
end

num_dfaces(cutter::CutterCache,d::Integer) = length(cutter.d_to_ldface_to_dface,d+1)

num_dims(cutter::CutterCache) = length(cutter.d_to_ldface_to_dface)-1

num_mesh_dims(a::CutterCache) = size(a.m_n_to_mface_to_nfaces,1)-1

function get_dface(cutter::CutterCache,d::Integer,i::Integer) 
  cutter.d_to_ldface_to_dface[d+1,i]
end

function get_faces(cutter::CutterCache,d::Integer,n::Integer) 
  cutter.m_n_to_mface_to_nfaces[d+1,n+1]
end

function set_dface!(cutter::CutterCache,d::Integer,i::Integer,val)
  cutter.d_to_ldface_to_dface[d+1,i] = val
end

function set_orientation!(cutter::CutterCache,cell::Integer,lfacet::Integer, val )
  cutter.cell_to_lfacet_to_orientation[cell,lfacet] = val
end

function _add_vertex!(mesh::CellMesh,dim::Integer,face::Integer,point::Point)
  if dim == 0
    return face
  end
  cache = mesh.cache
  @check isactive(mesh,dim,face)
  for d in dim:num_dims(mesh)
    df_to_nf = get_faces(mesh,d,dim)
    for iface in 1:length(df_to_nf)
      if isactive(df_to_nf,iface)
        for lface in 1:length(df_to_nf,iface)
          if df_to_nf[iface,lface] == face
            refine_dface!( cache, mesh, d,iface, dim,lface )
            append_mesh!( mesh, cache )
            remove_dface!( mesh, d,iface )
            break
          end
        end
      end
    end
  end
 _add_vertex!(mesh,point)
end

function append_mesh!(mesh::CellMesh,caches::MeshCache)
  dim = num_dims(caches.cutter)
  add_faces!(mesh,caches.cutter)
  mesh
end

function add_faces!(mesh::CellMesh,cutter::CutterCache)
  dim = num_dims(cutter)
  for d in 0:dim, n in 0:d
    append!( get_faces(mesh,d,n), get_faces(cutter,d,n) )
  end
  mesh
end

function refine_dface!(caches::MeshCache,mesh::CellMesh,dim::Integer,face::Integer,ldim::Integer,lface::Integer)
  cell_cache = caches.cell
  cutter_cache = caches.cutter

  ## Load cutter cache
  load_tables!(cutter_cache,dim,ldim,lface)
  setup_cell_to_lface_to_orientation!(cutter_cache,cell_cache,face)
  setup_d_to_tdface_to_dface!(cutter_cache,cell_cache,mesh,face,ldim,lface)
  setup_connectivities!(cutter_cache)
  remove_repeated_faces!(cutter_cache,mesh)

  ## Update cell cache
  update_dface_to_new_dfaces!(cell_cache,cutter_cache,face)
  update_cell_to_lface_to_orientation!(cell_cache,cutter_cache)
  update_dface_to_cell_dface!(cell_cache,cutter_cache,face)

  caches
end

function load_tables!(cutter::CutterCache,dim::Integer,ldim::Integer,lid::Integer)
  case = D_to_d_to_dface_to_case_for_cut_tetD[dim][ldim][lid]
  tables_nf_to_nF = D_to_case_to_d_to_subdface_to_dface_for_cut_tetD[dim][case]
  tables_nf_to_mf = D_to_case_to_m_n_to_mface_to_nfaces_for_cut_tetD[dim][case]
  tables_c_to_lf_to_o = D_to_case_to_cell_lfacet_to_orientation_for_cut_tetD[dim][case]

  resize!( cutter.d_to_ldface_to_dface, 0 )
  append!( cutter.d_to_ldface_to_dface, tables_nf_to_nF )
  for d in 0:dim, n in 0:d
    resize!( get_faces(cutter,d,n), 0 )
    append!( get_faces(cutter,d,n), tables_nf_to_mf[d+1,n+1] )
  end
  resize!( cutter.cell_to_lfacet_to_orientation, 0 )
  if num_mesh_dims(cutter) == dim
    append!( cutter.cell_to_lfacet_to_orientation, tables_c_to_lf_to_o )
  end
  cutter
end

function setup_cell_to_lface_to_orientation!(cutter::CutterCache,cache::CellMeshCache,face::Integer)
  D = num_dims(cache)
  d = num_dims(cutter)
  if d == D
    cell = face
    tcell_to_tfacets = get_faces(cutter,D,D-1)
    tfacet_to_facet = 
    for tcell in 1:length(tcell_to_tfacets)
      for ltfacet in 1:length(tcell_to_tfacets,tcell)
        tfacet = tcell_to_tfacets[tcell,ltfacet]
        if get_dface(cutter,D-1,tfacet) != 0
          lfacet = get_dface(cutter,D-1,tfacet)
          set_orientation!(cutter,tcell,ltfacet, get_orientation(cache,cell,lfacet) )
        end
      end
    end
  end
  cutter
end

function setup_d_to_tdface_to_dface!(
  cutter::CutterCache,
  cache::CellMeshCache,
  mesh::CellMesh,
  face::Integer,
  cut_dim::Integer,
  cut_dim_lface::Integer)
  

  dim = num_dims(cutter) 
  dim_face = face
  dim_f_to_cut_dim_f = get_faces(mesh,dim,cut_dim)
  cut_dim_face = dim_f_to_cut_dim_f[dim_face,cut_dim_lface]

  for d in 0:dim
    dim_f_to_df = get_faces(mesh,dim,d)
    cut_dim_f_to_df = get_faces(mesh,cut_dim,d)
    n_dfaces = num_dfaces(mesh,d)
    for tdface in 1:num_dfaces(cutter,d)
      if get_dface(cutter,d,tdface) == 0
        if d == 0 || d > cut_dim
          n_dfaces += 1
          set_dface!(cutter,d,tdface, n_dfaces)
        else
          _d, _lface = _get_lface_containing_ntface(cutter,d,tdface)
          if _lface != UNSET
 
            dim_f_to__df = get_faces(mesh,dim,_d)
            _df_to_df = get_faces(mesh,_d,d)

            _face = dim_f_to__df[dim_face,_lface]
            @check !isactive(mesh,_d,_face)
            _dface = UNSET
            for new_face in get_new_faces(cache,_d,_face)
              for ldface in 1:length(_df_to_df,new_face)
                dface = _df_to_df[new_face,ldface]
                if _are_dface_and_tdface_coincidents(mesh,cutter,d,dface,tdface)
                  _dface = dface
                  break
                end
              end
              if _dface != UNSET
                break
              end
            end
            @check _dface != UNSET
            set_dface!(cutter,d,tdface, _dface)
          else
            n_dfaces += 1
            set_dface!(cutter,d,tdface, n_dfaces)
          end
        end
      else
        dface = dim_f_to_df[dim_face,get_dface(cutter,d,tdface)]
        if !isactive(mesh,d,dface)
          @check d ≥ cut_dim
          for new_dface in get_new_faces(cache,d,dface)
            if _are_dface_and_tdface_coincidents(mesh,cutter,d,new_dface,tdface)
              dface = new_dface
              break
            end
          end
        end
        set_dface!(cutter,d,tdface, dface)
      end
    end
  end
end



function setup_connectivities!(cutter::CutterCache)
  dim = num_dims(cutter)
  for d in 0:dim, n in 0:d
    df_to_nf = get_faces(cutter,d,n)
    for i in 1:length(df_to_nf), j in 1:length(df_to_nf,i)
      df_to_nf[i,j] = get_dface(cutter,n,df_to_nf[i,j])
    end
  end
  cutter
end

function remove_repeated_faces!(cutter::CutterCache,mesh::CellMesh)
  dim = num_dims(cutter)
  for d in 0:dim-1, n in 0:d
    df_to_nf = get_faces(cutter,d,n)
    for i in 1:num_dfaces(cutter,d)
      if get_dface(cutter,d,i) <= length( get_faces(mesh,d,n) )
        remove!(df_to_nf,i)
      end
    end
  end
  cutter
end

function update_dface_to_new_dfaces!(cache::CellMeshCache,cutter::CutterCache,id::Integer)
  dim = num_dims(cutter)
  for d in 1:dim
    df_to_new_df = get_dface_to_new_dfaces(cache,d)
    n_dfaces = num_dfaces(cache,d)
    for i in 1:num_dfaces(cutter,d)
      dface = get_dface(cutter,d,i)
      if dface > num_dfaces(cache,d)
        set_num_dfaces!(cache,d,dface)
        if length(df_to_new_df) < dface
          resize!(df_to_new_df,dface)
        end
      end
      if !isassigned(df_to_new_df,dface)
        df_to_new_df[dface] = []
      elseif dface > n_dfaces
        resize!(df_to_new_df[dface],0)
      end
    end
  end
  for i in 1:num_dfaces(cutter,dim)
    push!( get_new_faces(cache,dim,id), get_dface(cutter,dim,i) )
  end
  cache 
end

function update_cell_to_lface_to_orientation!(cache::CellMeshCache,cutter::CutterCache)
  append!( cache.cell_to_lfacet_to_orientation, cutter.cell_to_lfacet_to_orientation )
  cache
end

function update_dface_to_cell_dface!(cache::CellMeshCache,cutter::CutterCache,face::Integer)
  dim = num_dims(cutter)
  D = num_dims(cache)
  for d in 0:dim-1
    df_to_cdf = cache.d_to_dface_to_cell_dface[d+1]
    for i in 1:num_dfaces(cutter,d)
      if get_dface(cutter,d,i) > length(df_to_cdf)
        @check  get_dface(cutter,d,i) == length(df_to_cdf) + 1
        push!(df_to_cdf,0)
      end
    end
  end
  if dim < D
    df_to_cdf = cache.d_to_dface_to_cell_dface[dim+1]
    for i in 1:num_dfaces(cutter,dim)
      push!(df_to_cdf, df_to_cdf[face] ) 
    end
  end
  cache 
end

function _get_lface_containing_ntface(cutter::CutterCache,n::Integer,nface::Integer)
  for d in n+1:num_dims(cutter)-1
    df_to_nf = get_faces(cutter,d,n)
    df_to_v = get_faces(cutter,d,0)
    nf_to_v = get_faces(cutter,n,0)
    for dface in 1:num_dfaces(cutter,d)
      if get_dface(cutter,d,dface) != 0
        for lnface in 1:length(df_to_nf,dface)
          _nface = df_to_nf[dface,lnface]
          if _nface == nface
            return d, get_dface(cutter,d,dface)
          end
        end
      end
    end
  end
  (UNSET,UNSET)
end


function _find_coincident_dface(
  mesh::CellMesh,
  cache::CellMeshCache,
  cutter::CutterCache,
  d::Integer,
  dface::Integer,
  tdface::Integer)

  for new_dface in get_new_faces(cache,d,dface)
    if _are_dface_and_tdface_coincidents(mesh,cutter,d,new_dface,tdface)
      return new_dface
    end
  end
  @check false
  return UNSET
end

function _are_dface_and_tdface_coincidents(
  mesh::CellMesh,
  cutter::CutterCache,
  d::Integer,
  dface::Integer,
  tdface::Integer)
  
  df_to_v = get_faces(mesh,d,0)
  tdf_to_tv = get_faces(cutter,d,0)
  for lvertex in 1:length(df_to_v,dface)
    coincident = false
    for tlvertex in 1:length(tdf_to_tv,tdface)
      if df_to_v[dface,lvertex] == get_dface( cutter, 0, tdf_to_tv[tdface,tlvertex] )
        coincident = true
        break
      end
    end
    if !coincident
      return false
    end
  end
  true
end

function push_surface_mesh_face!(cache::CellMeshCache,sm_face::Integer)
  n = num_dfaces(cache,0)
  set_num_dfaces!(cache,0,n+1)
  if length(cache.vertex_to_surface_mesh_faces) == n
    push!( cache.vertex_to_surface_mesh_faces, [] )
  else
    resize!( cache.vertex_to_surface_mesh_faces[n+1], 0 )
  end
  push!( cache.vertex_to_surface_mesh_faces[n+1], sm_face )
end

function add_surface_mesh_face_to_vertex!(cache::CellMeshCache,vertex::Integer,sm_face::Integer)
  push!( cache.vertex_to_surface_mesh_faces[vertex], sm_face )
end

function add_surface_mesh_face_to_vertex!(cache::MeshCache,d::Integer,face::Integer,sm_face::Integer)
  if d == 0
    add_surface_mesh_face_to_vertex!(cache.cell,face,sm_face)
  else
    push_surface_mesh_face!(cache.cell,sm_face)
  end
end

function compact_cache!(cache::MeshCache,mesh::CellMesh)
  cell_cache = cache.cell
  D = num_dims(mesh)
  for d in 0:D
    set_num_dfaces!( cell_cache, d, num_dfaces(mesh,d) )
  end
  for d in 1:D
    df_to_new_df = get_dface_to_new_dfaces(cache.cell,d)
    for dface in 1:num_dfaces(mesh,d)
      resize!(df_to_new_df[dface],0)
    end
  end
  for d in 0:D-1
    iface = 0
    for i in 1:num_dfaces(mesh,d)
      iface += 1
      if !isactive(mesh,d,i)
        deleteat!(cell_cache.d_to_dface_to_cell_dface[d+1],iface)
        iface -= 1
      end
    end
  end
  for i in 1:num_cells(mesh)
    if !isactive(mesh,D,i)
      remove!(cell_cache.cell_to_lfacet_to_orientation,i)
    end
  end
  compact!(cell_cache.cell_to_lfacet_to_orientation)
  cache
end

function compact_mesh!(mesh::CellMesh)
  D = num_dims(mesh)
  for d in 0:D
    n_dfaces = 0
    df_to_df = get_faces(mesh,d,d)
    for i in 1:length(df_to_df)
      if isactive(df_to_df,i)
        @check length(df_to_df,i) == 1
        n_dfaces += 1
        df_to_df[i,1] = n_dfaces
      end
    end
  end
  for d in 1:D, n in 0:d-1
    nf_to_nf = get_faces(mesh,n,n)
    df_to_nf = get_faces(mesh,d,n)
    for i in 1:length(df_to_nf)
      if isactive(df_to_nf,i)
        for j in 1:length(df_to_nf,i)
          df_to_nf[i,j] = nf_to_nf[ df_to_nf[i,j], 1 ]
        end
      end
    end
  end
  for d in 0:D, n in 0:d
    compact!( get_faces(mesh,d,n) )
  end
  mesh
end

function _find_surface_mesh_facet(mesh::CellMesh,cache::MeshCache,facet::Integer,sm::SurfaceMesh)
  intersection_cache = get_vector!(cache)
  intersection_cache_bis = get_vector_bis!(cache) 
  v_to_smf = cache.cell.vertex_to_surface_mesh_faces
  D = num_dims(mesh)
  f_to_v = get_dface_to_vertices(mesh,D-1)

  for lvertex in 1:length(f_to_v,facet)
    resize!(intersection_cache_bis,0)
    vertex = f_to_v[facet,lvertex]
    for sm_face in v_to_smf[vertex]
      d = face_dimension(sm,sm_face)
      sm_dface = local_dface(sm,sm_face,d)
      smdf_to_smf = get_faces(sm,d,D-1)
      for i in 1:length(smdf_to_smf,sm_dface)
        sm_facet = smdf_to_smf[sm_dface,i]
        push!(intersection_cache_bis,sm_facet)
      end
    end
    if lvertex == 1
      resize!(intersection_cache,0)
      append!(intersection_cache, intersection_cache_bis )
    else
      for i in 1:length(intersection_cache)
        if intersection_cache != UNSET
          if intersection_cache[i] ∉ intersection_cache_bis
             intersection_cache[i] = UNSET
          end
        end
      end
    end
  end
  for sm_facet in intersection_cache
    if sm_facet != UNSET
      return sm_facet
    end
  end
  UNSET
end

function compute_facet_to_surface_mesh_facet!(mesh::CellMesh,cache::MeshCache,sm::SurfaceMesh)
  D = num_dims(mesh)
  facet_to_surface_mesh_facet = cache.cell.facet_to_surface_mesh_facet
  resize!(facet_to_surface_mesh_facet,num_facets(mesh))
  for facet in 1:num_facets(mesh)
    if isactive(mesh,D-1,facet)
      facet_to_surface_mesh_facet[facet] = _find_surface_mesh_facet(mesh,cache,facet,sm)
    end
  end
end

function define_cells_touching_boundary!(mesh::CellMesh,cache::MeshCache,sm::SurfaceMesh)
  D = num_dims(mesh)
  c_to_f = get_faces(mesh,D,D-1)
  f_to_smf = cache.cell.facet_to_surface_mesh_facet

  for cell in 1:num_cells(mesh)
    if isactive(mesh,D,cell)
      for lfacet in 1:length(c_to_f,cell)
        facet = c_to_f[cell,lfacet]
        if f_to_smf[facet] != UNSET
          set_cell_in_out!(mesh,cell, _define_cell(mesh,cache,sm,cell,lfacet) )
          break
        end
      end
    end
  end
end

function _define_cell(mesh::CellMesh,cache::MeshCache,sm::SurfaceMesh,cell::Integer,lfacet::Integer)
  D = num_dims(mesh)
  cell_cache = cache.cell
  f_to_smf = cell_cache.facet_to_surface_mesh_facet
  c_to_lf_to_o = cell_cache.cell_to_lfacet_to_orientation
  c_to_f = get_faces(mesh,D,D-1)

  icell = cell
  ifacet = c_to_f[icell,lfacet]
  @check f_to_smf[ifacet] != UNSET
  cell = get_cell_coordinates(mesh,icell)
  facet = get_facet_coordinates(mesh,ifacet)
  
  facet_normal = normal(facet)
  facet_normal = facet_normal / norm(facet_normal)

  sm_facet_normal = get_facet_normal(sm,f_to_smf[ifacet])

  @check relative_orientation(facet,cell) == c_to_lf_to_o[icell,lfacet]
  orientation = ( facet_normal ⋅ sm_facet_normal ) * c_to_lf_to_o[icell,lfacet]

  if orientation > 0
    FACE_IN
  elseif orientation < 0
    FACE_OUT
  else
    @check false "orientation == 0, face not defined"
    FACE_UNDEF
  end
end

function reset_d_to_dface_to_in_out_boundary!(mesh::CellMesh,cache::MeshCache)
  D = num_dims(mesh)
  lens = get_vector!(cache)
  resize!(lens,D+1)
  for d in 0:D
    lens[d+1] = num_dfaces(mesh,d)
  end
  d_to_df_to_io = mesh.d_to_dface_to_in_out_boundary
  resize!(d_to_df_to_io,lens)
  fill!(d_to_df_to_io,FACE_UNDEF)
  mesh
end

function define_boundary_faces!(mesh::CellMesh,cache::MeshCache)
  D = num_dims(mesh)
  d = D-1
  f_to_smf = cache.cell.facet_to_surface_mesh_facet
  for facet in 1:num_facets(mesh)
    if f_to_smf[facet] != UNSET
      for d in 0:D-1
        f_to_df = get_faces(mesh,D-1,d)
        for ldf in 1:length(f_to_df,facet)
          dface = f_to_df[facet,ldf]
          set_face_as_boundary!(mesh,d,dface)
        end
      end
    end
  end
  mesh
end

function define_faces!(mesh::CellMesh)
  D = num_dims(mesh)
  for cell in 1:num_cells(mesh)
    cell_io = get_cell_in_out(mesh,cell)
    for d in 0:D-1
      c_to_df = get_faces(mesh,D,d)
      for ldface in 1:length(c_to_df,cell)
        dface = c_to_df[cell,ldface]
        if !is_face_defined(mesh,d,dface)
          set_face_in_out_boundary!(mesh,d,dface,cell_io)
        end
      end
    end
  end
  mesh
end

function fix_boundary_facets!(mesh::CellMesh)
  cache = get_cache(mesh)
  f_to_smf = cache.cell.facet_to_surface_mesh_facet
  c_to_f = get_cell_to_facets(mesh)
  for facet in 1:num_facets(mesh)
    set_facet_as_exterior!(mesh,facet)
  end
  for cell in 1:num_cells(mesh)
    if is_cell_interior(mesh,cell)
      for lfacet in 1:length(c_to_f,cell)
        facet = c_to_f[cell,lfacet]
        if f_to_smf[facet] != UNSET
          set_facet_as_boundary!(mesh,facet)
          _fix_cell_to_facet_orientation!(mesh,cell,lfacet)
        end
      end
    end
  end
end

function _fix_cell_to_facet_orientation!(mesh::CellMesh,cell::Integer,lfacet::Integer)
  cell_cache = get_cache(mesh).cell
  c_to_f = get_cell_to_facets(mesh)
  facet = c_to_f[cell,lfacet]
  if get_orientation(cell_cache,cell,lfacet) < 0
    set_orientation!(cell_cache,cell,lfacet, 1 )
    _flip_facet!(mesh,facet)
  end
end

function _flip_facet!(mesh::CellMesh,facet::Integer)
  f_to_v = get_facet_to_vertices(mesh)
  n_vxf = length(f_to_v,facet)
  if  n_vxf > 1
    _v = f_to_v[facet,n_vxf]
    f_to_v[facet,n_vxf] = f_to_v[facet,n_vxf-1] 
    f_to_v[facet,n_vxf-1] = _v
  end
end


function compute_facet_to_cells!(mesh::CellMesh)
  D = num_dims(mesh)
  compute_nface_to_dfaces_dual!(mesh,D,D-1)
end

function compute_nface_to_dfaces_dual!(mesh::CellMesh,n::Integer,d::Integer)
  n_dfaces = num_dfaces(mesh,d)
  n_nfaces = num_dfaces(mesh,n)
  nf_to_df = get_faces(mesh,n,d)
  lens = get_vector!(mesh.cache)
  resize!(lens,n_dfaces)
  fill!(lens,0)
  for i in 1:n_nfaces, j in 1:length(nf_to_df,i)
    dface = nf_to_df[i,j]
    lens[dface] += 1
  end
  df_to_nf = get_faces(mesh,d,n)
  resize!(df_to_nf,lens)
  fill!(df_to_nf,UNSET)
  fill!(lens,0)
  for i in 1:n_nfaces, j in 1:length(nf_to_df,i)
    dface = nf_to_df[i,j]
    lens[dface] += 1
    df_to_nf[dface,lens[dface]] = i
  end
  mesh
end

function _get_in_out_face_data(m::CellMesh,range)
  data = Int[]
  for d in range, i in 1:num_dfaces(m,d)
    push!(data,get_face_in_out_boundary(m,d,i))
  end
  data
end

function is_any_face_repeated(mesh::CellMesh)
  D = num_dims(mesh)
  for d in 0:D
    df_to_v = get_dface_to_vertices(mesh,d)
    for dface in 1:length(df_to_v)
      if isactive(df_to_v,dface)
        face_found = true
        for _dface in 1:length(df_to_v)
          if _dface != dface
            for lvertex in 1:length(df_to_v,dface)
              vertex = df_to_v[dface,lvertex]
              vertex_found = false
              for _lvertex in 1:length(df_to_v,_dface)
                _vertex = df_to_v[_dface,_lvertex]
                if _vertex == vertex
                  vertex_found = true
                  break
                end
              end
              if !vertex_found
                face_found = false
              end
            end
            if face_found
              return true
            end
          end
        end
      end
    end
  end
  false
end
