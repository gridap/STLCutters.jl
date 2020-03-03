
struct CellMesh{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  m_n_to_mface_to_nfaces::Matrix{Table{Int}}
  d_to_dface_to_new_dfaces::Vector{Vector{Vector{Int}}}
end

struct CellMeshCache

end

struct CutterCache{D}
  d_to_ldface_to_dface::Table{Int}
  m_n_to_mface_to_nfaces::Matrix{Table{Int}}
end

struct MeshCaches
  cell::CellMeshCache
  cutter::CutterCache
end

CellMesh(b::BoundingBox) = CellMesh(Hexahedron(b))

function CellMesh(c::Hexahedron{N,D,T}) where {N,D,T}
  coordinates = [ v for v in get_vertices(c) ]
  m_n_to_mf_to_nf = Matrix{Table{Int}}(undef,D+1,D+1)

  table_m_n_to_mf_to_nf = D_to_m_n_to_mface_to_nfaces_for_hexD_to_tetD[D]

  for d in 0:D, n in 0:d
    df_to_nf = table_m_n_to_mf_to_nf[d+1,n+1]
    m_n_to_mf_to_nf[d+1,n+1] = Table( df_to_nf )
  end

  d_to_df_to_new_df = [ Vector{Int}[] for d in 1:D ]

  for d in 1:D
    d_to_df_to_new_df[d] = [ Int[] for i in 1:length( m_n_to_mf_to_nf[d+1,0+1] ) ]
  end

  CellMesh{D,T}( coordinates, m_n_to_mf_to_nf, d_to_df_to_new_df )
end

CutterCache(::T) where T<:CellMesh = CutterCache(T)

function CutterCache(::Type{<:CellMesh{D}}) where D
  d_to_ldf_to_df = Table{Int}()
  m_n_to_mf_to_nf = Matrix(undef,D+1,D+1)
  for i in 1:D+1, j in 1:D+1
    m_n_to_mf_to_nf[i,j] = Table{Int}()
  end
  CutterCache{D}( d_to_ldf_to_df, m_n_to_mf_to_nf )
end

function num_vertices(m::CellMesh) 
  length(m.vertex_coordinates)
end

function num_dfaces(m::CellMesh,d::Int)
  if d == 0
    return length(m.vertex_coordinates)
  end
  length( get_faces(m,d,0) ) 
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

function initialize!(m::CellMesh{D,T},c::Hexahedron{N,D,T}) where {N,D,T} 
  table_m_n_to_mface_to_nfaces = D_to_m_n_to_mface_to_nfaces_for_hexD_to_tetD[D]
  resize!(m.vertex_coordinates,N)
  for (i,v) in enumerate( get_vertices(c) )
    m.vertex_coordinates[i] = v
  end
  for d in 0:D, n in 0:d
    df_to_nf = table_m_n_to_mface_to_nfaces[d+1,n+1]
    resize!( get_faces(m,d,n), 0 )
    append!( get_faces(m,d,n), df_to_nf )
  end
  for d in 1:D
    resize!( m.d_to_dface_to_new_dfaces[d], num_dfaces(m,d) )
    for i in 1:num_dfaces(m,d)
      resize!( m.d_to_dface_to_new_dfaces[d][i],0)
    end
  end
  m
end

initialize!(m::CellMesh,b::BoundingBox) = initialize!(m,Hexahedron(b))

function add_vertex!(m::CellMesh{D,T},p::Point{D,T}) where {D,T}
  v = get_vertex_coordinates(m)
  push!( v, p )
end

function delete_dface!(m::CellMesh,d::Integer,i::Int)
  for n in 0:d
    df_to_nf = get_faces(m,d,n)
    remove!(df_to_nf,i)
  end
end

function load_tables!(cutter::CutterCache,dim::Integer,ldim::Integer,lid::Integer)
  case = D_to_d_to_dface_to_case_for_cut_tetD[dim][ldim][lid]
  tables_nf_to_nF = D_to_case_to_d_to_subdface_to_dface_for_cut_tetD[dim][case]
  tables_nf_to_mf = D_to_case_to_m_n_to_mface_to_nfaces_for_cut_tetD[dim][case]
  resize!( cutter.d_to_ldface_to_dface, 0 )
  append!( cutter.d_to_ldface_to_dface, tables_nf_to_nF )
  for d in 0:dim, n in 0:d
    resize!( get_faces(cutter,d,n), 0 )
    append!( get_faces(cutter,d,n), tables_nf_to_mf[d+1,n+1] )
  end
end

function setup_dfaces!(
  cutter::CutterCache,
  mesh::CellMesh,
  id::Integer,
  ldim::Integer,
  lid::Integer)

  dim = num_dims(cutter) 
  for n in 0:dim
    df_to_nf = get_faces(mesh,dim,n)
    n_dfaces = num_dfaces(mesh,n)
    for i in 1:num_dfaces(cutter,n)
      if get_dface(cutter,n,i) == 0
        n_dfaces += 1
        set_dface!( cutter, n, i, n_dfaces )
      else
        if n == ldim && get_dface(cutter,n,i) == lid
          set_dface!( cutter, n, i, _find_coincident_dface(mesh,cutter,dim,id,ldim,i) )
        else
          set_dface!( cutter, n, i, df_to_nf[ id, get_dface(cutter,n,i)  ] )
        end
      end
    end
  end
  cutter
end

function _find_coincident_dface(m::CellMesh,c::CutterCache,dim::Integer,id::Integer,ldim::Integer,ldf::Integer)
  lid = get_dface(c,ldim,ldf)
  df_id = get_faces(m,dim,ldim)[id,lid]
  df_to_v = get_faces(m,ldim,0)
  ldf_to_lv = get_faces(c,ldim,0)
  for new_df in m.d_to_dface_to_new_dfaces[ldim][df_id]
    found = true
    for i in 1:length(df_to_v,df_id)
      coincident = false
      for j in 1:length(ldf_to_lv,ldf)
        if df_to_v[new_df,i] == get_dface( c, 0, ldf_to_lv[ldf,j] )
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
      return new_df
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
      if get_dface(cutter,d,i) <= num_dfaces(mesh,d)
        remove!(df_to_nf,i)
      end
    end
  end
  cutter
end

function update_dface_to_new_dfaces!(mesh::CellMesh,cutter::CutterCache,id::Integer)
  dim = num_dims(cutter)
  for d in 1:dim
    df_to_new_df = mesh.d_to_dface_to_new_dfaces[d]
    for i in 1:num_dfaces(cutter,d)
      dface = get_dface(cutter,d,i)
      if dface > length( df_to_new_df )
        resize!(df_to_new_df,dface)
      end
      if !isassigned(df_to_new_df,dface)
        df_to_new_df[dface] = Int[]
      elseif dface > num_dfaces(mesh,d)
        resize!(df_to_new_df[dface],0)
      end
    end
  end
  for i in 1:num_dfaces(cutter,dim)
    push!( mesh.d_to_dface_to_new_dfaces[dim][id], get_dface(cutter,dim,i) )
  end
  mesh
end

function add_faces!(mesh::CellMesh,cutter::CutterCache)
  dim = num_dims(cutter)
  for d in 0:dim, n in 0:d
    append!( get_faces(mesh,d,n), get_faces(cutter,d,n) )
  end
  mesh
end

function refine!(cutter::CutterCache,mesh::CellMesh,dim::Integer,id::Integer,ldim::Integer,lid::Integer)
  load_tables!(cutter,dim,ldim,lid)
  setup_dfaces!(cutter,mesh,id,ldim,lid)
  setup_connectivities!(cutter)
  remove_repeated_faces!(cutter,mesh)
  cutter
end


function update!(mesh::CellMesh,cutter::CutterCache,id::Integer)
  dim = num_dims(cutter)
  update_dface_to_new_dfaces!(mesh,cutter,id)
  add_faces!(mesh,cutter)
  delete_dface!(mesh,dim,id)
  mesh
end

function refine!(mesh::CellMesh{D},cutter::CutterCache,dim::Integer,id::Integer) where D
  if dim == 0
    return mesh
  end
  for d in dim:D
    df_to_nf = get_faces(mesh,d,dim)
    for i in 1:length(df_to_nf)
      if isactive(df_to_nf,i)
        for j in 1:length(df_to_nf,i)
          if df_to_nf[i,j] == id 
            refine!(cutter,mesh,d,i,dim,j)
            update!(mesh,cutter,i)
            break
          end
        end
      end
    end
  end
  mesh
end

function writevtk(m::CellMesh{D,T},file_base_name) where {D,T}
  d_to_vtk_type_id = Dict(0=>1,1=>3,2=>5,3=>10)
  num_points = num_vertices(m)
  points = zeros(T,D,num_points)
  for (i ,p ) in enumerate(get_vertex_coordinates(m)), d in 1:D
    points[d,i] = p[d]
  end
  cells = MeshCell{Vector{Int64}}[]
  for d in 0:D
    dface_to_vertices = get_dface_to_vertices(m,d)
    num_dfaces = length(dface_to_vertices)
    vtk_type = VTKCellType(d_to_vtk_type_id[d])
    for i in 1:num_dfaces
      vertices = [ dface_to_vertices[i,j] for j in 1:vtk_type.nodes ]
      push!( cells, MeshCell(vtk_type,vertices) )
    end
  end
  vtkfile = vtk_grid(file_base_name,points,cells)
  vtk_save(vtkfile)
end

