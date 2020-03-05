
const FACE_UNDEF = 0
const FACE_IN = 1
const FACE_OUT = 2
const FACE_BOUNDARY = 3

struct CellMesh{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  m_n_to_mface_to_nfaces::Matrix{Table{Int}}
  d_to_dface_to_in_out_boundary::Table{Int}
end

struct CellMeshCache
  d_to_dface_to_new_dfaces::Vector{Vector{Vector{Int}}}
  d_to_dface_to_cell_dface::Vector{Vector{Int}}
  cell_to_lfacet_to_orientation::Table{Int}
  d_to_dface_to_surface_mesh_faces::Vector{Vector{Vector{Int}}}
end

struct CutterCache
  d_to_ldface_to_dface::Table{Int}
  m_n_to_mface_to_nfaces::Matrix{Table{Int}}
  cell_to_lfacet_to_orientation::Table{Int}
end

struct MeshCaches
  cell::CellMeshCache
  cutter::CutterCache
end

function CellMesh(box::BoundingBox{D,T}) where {D,T}
  table_m_n_to_mf_to_nf = D_to_m_n_to_mface_to_nfaces_for_hexD_to_tetD[D]
  
  coordinates = [ v for v in get_vertices(box) ]
  m_n_to_mf_to_nf = Matrix{Table{Int}}(undef,D+1,D+1)
  d_df_to_io = Table{Int}()

  for d in 0:D, n in 0:d
    df_to_nf = table_m_n_to_mf_to_nf[d+1,n+1]
    m_n_to_mf_to_nf[d+1,n+1] = Table( df_to_nf )
  end


  CellMesh{D,T}( coordinates, m_n_to_mf_to_nf, d_df_to_io )
end

function CellMeshCache(mesh::CellMesh{D}) where D
  table_c_lf_to_o = D_to_cell_lfacet_to_orientation_for_hexD_to_tetD[D] 
  table_d_to_df_to_cdf = D_to_d_to_subdface_to_dface_for_hexD_to_tetD[D]
  
  d_to_df_to_new_df = [ Vector{Int}[] for d in 1:D ]
  d_to_df_to_cdf = Vector{Vector{Int}}(undef,D)
  c_lf_to_o = Table( table_c_lf_to_o )
  d_to_df_to_smf = [ Vector{Int}[] for d in 0:D-1 ]

  for d in 1:D
    d_to_df_to_new_df[d] = [ Int[] for i in 1:num_dfaces(mesh,d) ]
  end

  for d in 0:D-1
    d_to_df_to_smf[d+1] = [ Int[] for i in 1:num_dfaces(mesh,d) ]
  end

  for d in 0:D-1
    df_to_cdf = table_d_to_df_to_cdf[d+1]
    d_to_df_to_cdf[d+1] = []
    append!( d_to_df_to_cdf[d+1], df_to_cdf )
  end

  CellMeshCache( d_to_df_to_new_df, d_to_df_to_cdf, c_lf_to_o, d_to_df_to_smf )
end

CutterCache(::T) where T<:CellMesh = CutterCache(T)

function CutterCache(::Type{<:CellMesh{D}}) where D
  d_to_ldf_to_df = Table{Int}()
  m_n_to_mf_to_nf = Matrix(undef,D+1,D+1)
  c_ldf_to_o = Table{Int}()
  for i in 1:D+1, j in 1:D+1
    m_n_to_mf_to_nf[i,j] = Table{Int}()
  end
  CutterCache( d_to_ldf_to_df, m_n_to_mf_to_nf, c_ldf_to_o )
end

function MeshCaches(mesh::CellMesh)
  cell_cache = CellMeshCache(mesh) 
  cutter_cache = CutterCache(mesh)
  MeshCaches( cell_cache, cutter_cache )
end

get_cache(mesh::CellMesh) = MeshCaches(mesh)

function num_vertices(m::CellMesh) 
  length(m.vertex_coordinates)
end

function num_dfaces(m::CellMesh,d::Int)
  if d == 0
    length(m.vertex_coordinates)
  else
    length( get_faces(m,d,0) ) 
  end
end

function get_faces(m::CellMesh,d::Int,n::Int)
  m.m_n_to_mface_to_nfaces[d+1,n+1]
end

function get_dface_to_vertices(m::CellMesh,d::Int)
  get_faces(m,d,0)
end

function get_vertex_coordinates(m::CellMesh) 
  m.vertex_coordinates
end

function get_face(m::CellMesh,::Val{0},i::Integer)
  v = get_vertex_coordinates(m)
  v[i]
end

function get_face(m::CellMesh,::Val{1},i::Integer)
  df_to_v = get_dface_to_vertices(m,1)
  v = get_vertex_coordinates(m)
  Segment( v[ df_to_v[i,1] ], v[ df_to_v[i,2] ] )
end

function get_face(m::CellMesh,::Val{2},i::Integer)
  df_to_v = get_dface_to_vertices(m,2)
  v = get_vertex_coordinates(m)
  Triangle( v[ df_to_v[i,1] ], v[ df_to_v[i,2] ], v[ df_to_v[i,3] ] )
end

function get_face(m::CellMesh,::Val{3},i::Integer)
  df_to_v = get_dface_to_vertices(m,3)
  v = get_vertex_coordinates(m)
  Tetrahedron( v[ df_to_v[i,1]], v[ df_to_v[i,2] ], v[ df_to_v[i,3] ], v[ df_to_v[i,4] ] )
end

function get_face(m::CellMesh,::Val{d},i::Integer) where d
  throw(ArgumentError("get_face(::CellMesh,::Val($d),::Integer) not implemented"))
end

function isactive(m::CellMesh,d::Integer,i::Integer)
  d != 0 || return true
  isactive(get_dface_to_vertices(m,d),i)
end

function get_vertex(m::CellMesh,i::Integer)
  get_face(m,Val{0}(),i)
end

function get_edge(m::CellMesh,i::Integer)
  get_face(m,Val{1}(),i)
end

function get_facet(m::CellMesh{D},i::Integer) where D
  get_face(m,Val{D-1}(),i)
end

function get_cell(m::CellMesh{D},i::Integer) where D
  get_face(m,Val{D}(),i)
end

function facet_normal(m::CellMesh{D},i::Integer) where D
  f = get_facet(m,i)
  normal(f)
end


function face_center(m::CellMesh,::Val{d},i::Integer) where d
  f = get_face(m,Val{d}(),i)
  center(f)
end

function facet_center(m::CellMesh{D},i::Integer) where D
  face_center(m,Val{D-1}(),i)
end

num_dims(cache::CellMeshCache) = length(cache.d_to_dface_to_new_dfaces)

num_dfaces(cutter::CutterCache,d::Integer) = length(cutter.d_to_ldface_to_dface,d+1)

num_dims(cutter::CutterCache) = length(cutter.d_to_ldface_to_dface)-1

function get_dface(cutter::CutterCache,d::Integer,i::Integer) 
  cutter.d_to_ldface_to_dface[d+1,i]
end

function set_dface!(cutter::CutterCache,d::Integer,i::Integer,val)
  cutter.d_to_ldface_to_dface[d+1,i] = val
end

function get_faces(cutter::CutterCache,d::Integer,n::Integer) 
  cutter.m_n_to_mface_to_nfaces[d+1,n+1]
end

function set_orientation!(cutter::CutterCache,cell::Integer,lfacet::Integer, val )
  cutter.cell_to_lfacet_to_orientation[cell,lfacet] = val
end

@generated function distance(m::CellMesh{D},d::Int,i::Int,p::Point{D}) where D
  body = ""
  for d in 0:D
    body *= "if d == $d \n"
    body *= "  f$d = get_face(m,Val{$d}(),i) \n"
    body *= "  distance(p,f$d) \n"
    body *= "else"
  end
  error = "  throw(ArgumentError(\"\$d-face does not exist\"))\n"
  str = body * '\n' * error * "end"
  Meta.parse(str)
end

function initialize!(m::CellMesh{D,T},box::BoundingBox{D,T}) where {D,T} 
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

  m
end


function initialize!(cache::CellMeshCache,mesh::CellMesh{D}) where D
  table_c_lf_to_o = D_to_cell_lfacet_to_orientation_for_hexD_to_tetD[D] 
  table_d_to_df_to_cdf = D_to_d_to_subdface_to_dface_for_hexD_to_tetD[D]
  
  for d in 1:D
    df_to_new_df = get_dface_to_new_dfaces(cache,d)
    resize!( df_to_new_df, num_dfaces(mesh,d) )
    for i in 1:num_dfaces(mesh,d)
      resize!( df_to_new_df[i],0)
    end
  end

  resize!( cache.cell_to_lfacet_to_orientation, 0 )
  append!( cache.cell_to_lfacet_to_orientation, table_c_lf_to_o )
 
  for d in 0:D-1
    df_to_smf = cache.d_to_dface_to_surface_mesh_faces[d+1]
    resize!( df_to_smf, num_dfaces(mesh,d) )
    for i in 1:num_dfaces(mesh,d)
      resize!(df_to_smf[i],0) 
    end
  end
  
  for d in 0:D-1
    df_to_cdf = table_d_to_df_to_cdf[d+1]
    resize!( cache.d_to_dface_to_cell_dface[d+1], 0 )
    append!( cache.d_to_dface_to_cell_dface[d+1], df_to_cdf )
  end

  cache
end

function initialize!(cutter::CutterCache)
  resize!(cutter.d_to_ldface_to_dface,0)
  resize!(cutter.cell_to_lfacet_to_orientation,0)
  for t in cutter.m_n_to_mface_to_nfaces
    resize!(t,0)
  end
  cutter
end

function initialize_cache!(caches::MeshCaches,mesh::CellMesh)
  initialize!(caches.cell,mesh)
  initialize!(caches.cutter)
  caches
end

get_orientation(a::CellMeshCache,cell::Integer,lfacet::Integer) = a.cell_to_lfacet_to_orientation[cell,lfacet]

get_dface_to_new_dfaces(a::CellMeshCache,d) = a.d_to_dface_to_new_dfaces[d]

get_new_faces(a::CellMeshCache,d,i) = a.d_to_dface_to_new_dfaces[d][i]

function add_vertex!(m::CellMesh{D,T},p::Point{D,T}) where {D,T}
  v = get_vertex_coordinates(m)
  push!( v, p )
end

function remove_dface!(m::CellMesh,d::Integer,i::Int)
  for n in 0:d
    df_to_nf = get_faces(m,d,n)
    remove!(df_to_nf,i)
  end
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
  append!( cutter.cell_to_lfacet_to_orientation, tables_c_to_lf_to_o )

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

function setup_d_to_ldface_to_dface!(
  cutter::CutterCache,
  cache::CellMeshCache,
  mesh::CellMesh,
  face::Integer,
  ldim::Integer,
  lface::Integer)

  dim = num_dims(cutter) 
  for n in 0:dim
    df_to_nf = get_faces(mesh,dim,n)
    n_dfaces = num_dfaces(mesh,n)
    for iface in 1:num_dfaces(cutter,n)
      if get_dface(cutter,n,iface) == 0
        n_dfaces += 1
        set_dface!( cutter, n,iface, n_dfaces )
      else
        if n == ldim && get_dface(cutter,n,iface) == lface
          set_dface!( cutter, n,iface, _find_coincident_dface(mesh,cache,cutter, dim,face, ldim,iface) )
        else
          set_dface!( cutter, n,iface, df_to_nf[face,get_dface(cutter,n,iface)] )
        end
      end
    end
  end
  cutter
end

function _find_coincident_dface(
  mesh::CellMesh,
  cache::CellMeshCache,
  cutter::CutterCache,
  dim::Integer,
  face::Integer,
  ldim::Integer,
  tdface::Integer)

  ldface = get_dface(cutter,ldim,tdface)
  dface = get_faces(mesh,dim,ldim)[face,ldface]
  df_to_v = get_faces(mesh,ldim,0)
  tdf_to_tv = get_faces(cutter,ldim,0)
  for new_dface in get_new_faces(cache,ldim,dface)
    found = true
    for lvertex in 1:length(df_to_v,dface)
      coincident = false
      for tlvertex in 1:length(tdf_to_tv,tdface)
        if df_to_v[new_dface,lvertex] == get_dface( cutter, 0, tdf_to_tv[tdface,tlvertex] )
          coincident = true
          break
        end
      end
      if !coincident
        found = false
        break
      end
    end
    if found
      return new_dface
    end
  end
  @check false
  return UNSET
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
    n_dfaces = length(df_to_new_df)
    for i in 1:num_dfaces(cutter,d)
      dface = get_dface(cutter,d,i)
      if dface > length( df_to_new_df )
        resize!(df_to_new_df,dface)
      end
      if !isassigned(df_to_new_df,dface)
        df_to_new_df[dface] = Int[]
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


function add_faces!(mesh::CellMesh,cutter::CutterCache)
  dim = num_dims(cutter)
  for d in 0:dim, n in 0:d
    append!( get_faces(mesh,d,n), get_faces(cutter,d,n) )
  end
  mesh
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


function refine_dface!(caches::MeshCaches,mesh::CellMesh,dim::Integer,face::Integer,ldim::Integer,lface::Integer)
  cell_cache = caches.cell
  cutter_cache = caches.cutter

  ## Load cutter cache
  load_tables!(cutter_cache,dim,ldim,lface)
  setup_cell_to_lface_to_orientation!(cutter_cache,cell_cache,face)
  setup_d_to_ldface_to_dface!(cutter_cache,cell_cache,mesh,face,ldim,lface)
  setup_connectivities!(cutter_cache)
  remove_repeated_faces!(cutter_cache,mesh)

  ## Update cell cache
  update_dface_to_new_dfaces!(cell_cache,cutter_cache,face)
  update_cell_to_lface_to_orientation!(cell_cache,cutter_cache)
  update_dface_to_cell_dface!(cell_cache,cutter_cache,face)

  caches
end

function append_mesh!(mesh::CellMesh,caches::MeshCaches)
  dim = num_dims(caches.cutter)
  add_faces!(mesh,caches.cutter)
  mesh
end

function add_vertex!(mesh::CellMesh{D},caches::MeshCaches,dim::Integer,face::Integer,point::Point) where D
  if dim == 0
    return mesh
  end
  for d in dim:D
    df_to_nf = get_faces(mesh,d,dim)
    for iface in 1:length(df_to_nf)
      if isactive(df_to_nf,iface)
        for lface in 1:length(df_to_nf,iface)
          if df_to_nf[iface,lface] == face
            refine_dface!( caches, mesh, d,iface, dim,lface )
            append_mesh!( mesh, caches )
            remove_dface!( mesh, d,iface )
            break
          end
        end
      end
    end
  end
  add_vertex!(mesh,point)
  mesh
end

function compact!(mesh::CellMesh{D}) where D
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

function writevtk(m::CellMesh{D,T},file_base_name) where {D,T}
  d_to_vtk_type_id = Dict(0=>1,1=>3,2=>5,3=>10)
  num_points = num_vertices(m)
  points = zeros(T,D,num_points+num_dfaces(m,D-1))
  for (i ,p ) in enumerate(get_vertex_coordinates(m)), d in 1:D
    points[d,i] = p[d]
  end
  for i in 1:num_dfaces(m,D-1)
    center = facet_center(m,i)
    for d in 1:D
      points[d,i+num_points] = center[d]
    end
  end

  cells = MeshCell{Vector{Int64}}[]
  for d in D
    dface_to_vertices = get_dface_to_vertices(m,d)
    vtk_type = VTKCellType(d_to_vtk_type_id[d])
    for i in 1:num_dfaces(m,d)
      vertices = [ dface_to_vertices[i,j] for j in 1:vtk_type.nodes ]
      push!( cells, MeshCell(vtk_type,vertices) )
    end
  end

  normals = zeros(3,num_points+num_dfaces(m,D-1))
  for i in 1:num_dfaces(m,D-1)
    normal = facet_normal(m,i)
    normal = normal / norm(normal)
    for d in 1:D
      normals[d,i+num_points] = normal[d]
    end
  end

  vtkfile = vtk_grid(file_base_name,points,cells)
  vtkfile["facet_normals",VTKPointData()] = normals

  vtk_save(vtkfile)
end

