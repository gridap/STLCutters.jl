
struct CellSubMesh{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  dface_to_nfaces::Matrix{Table{Int}}
  dface_to_new_dfaces::Vector{Vector{Vector{Int}}}
end

struct FaceCutter{D}
  dfaces::Table{Int}
  dface_to_nfaces::Matrix{Table{Int}}
end

CellSubMesh(b::BoundingBox) = CellSubMesh(Hexahedron(b))

function CellSubMesh(c::Hexahedron{N,D,T}) where {N,D,T}
  coordinates = [ v for v in get_vertices(c) ]
  dface_to_nfaces = Matrix{Table{Int}}(undef,D+1,D+1)

  for d in 0:D, n in 0:d
    df_to_nf = hex_to_tet_nface_to_mface[D][d+1][n+1]
    dface_to_nfaces[d+1,n+1] = Table( df_to_nf )
  end

  nf_to_new_df = [ Vector{Int}[] for d in 1:D ]

  for d in 1:D
    nf_to_new_df[d] = [ Int[] for i in 1:length( dface_to_nfaces[d+1,0+1] ) ]
  end

  CellSubMesh{D,T}(coordinates,dface_to_nfaces,nf_to_new_df)
end

FaceCutter(::T) where T<:CellSubMesh = FaceCutter(T)

function FaceCutter(::Type{<:CellSubMesh{D}}) where D
  dfaces = zero(Table{Int})
  dface_to_nfaces = Matrix(undef,D+1,D+1)
  for i in 1:D+1, j in 1:D+1
    dface_to_nfaces[i,j] = zero(Table{Int})
  end
  FaceCutter{D}(dfaces,dface_to_nfaces)
end

function num_vertices(m::CellSubMesh) 
  length(m.vertex_coordinates)
end

function num_dfaces(m::CellSubMesh,d::Int)
  if d == 0
    num_vertices(m)
  else
    length(m.dface_to_nfaces[d+1,0+1])
  end
end

function dface_to_nfaces(m::CellSubMesh,d::Int,n::Int)
  m.dface_to_nfaces[d+1,n+1]
end

function dface_vertices(m::CellSubMesh,d::Int) 
  m.dface_to_nfaces[d+1,0+1]
end

function get_vertex_coordinates(m::CellSubMesh) 
  m.vertex_coordinates
end

function get_dface(m::CellSubMesh,::Val{0},i::Integer)
  v = get_vertex_coordinates(m)
  v[i]
end

function get_dface(m::CellSubMesh,::Val{1},i::Integer)
  df_to_v = dface_vertices(m,1)
  v = get_vertex_coordinates(m)
  Segment( v[ df_to_v[i,1] ], v[ df_to_v[i,2] ] )
end

function get_dface(m::CellSubMesh,::Val{2},i::Integer)
  df_to_v = dface_vertices(m,2)
  v = get_vertex_coordinates(m)
  Triangle( v[ df_to_v[i,1] ], v[ df_to_v[i,2] ], v[ df_to_v[i,3] ] )
end

function get_dface(m::CellSubMesh,::Val{3},i::Integer)
  df_to_v = dface_vertices(m,3)
  v = get_vertex_coordinates(m)
  Tetrahedron( v[ df_to_v[i,1]], v[ df_to_v[i,2] ], v[ df_to_v[i,3] ], v[ df_to_v[i,4] ] )
end

function get_dface(m::CellSubMesh,::Val{d},i::Integer) where d
  throw(ArgumentError("get_dface(::CellSubMesh,::Val($d),::Integer) not implemented"))
end

function isactive(m::CellSubMesh,d::Integer,i::Integer)
  d != 0 || return true
  isactive(dface_vertices(m,d),i)
end

function get_vertex(m::CellSubMesh,i::Integer)
  get_dface(m,Val{0}(),i)
end

function get_edge(m::CellSubMesh,i::Integer)
  get_dface(m,Val{1}(),i)
end

function get_facet(m::CellSubMesh{D},i::Integer) where D
  get_dface(m,Val{D-1}(),i)
end

function get_cell(m::CellSubMesh{D},i::Integer) where D
  get_dface(m,Val{D}(),i)
end

num_dfaces(cutter::FaceCutter,d::Integer) = length(cutter.dfaces,d+1)

num_dims(cutter::FaceCutter) = length(cutter.dfaces)-1

get_dface(cutter::FaceCutter,d::Integer,i::Integer) = cutter.dfaces[d+1,i]

set_dface!(cutter::FaceCutter,d::Integer,i::Integer,val) = cutter.dfaces[d+1,i] = val

dface_to_nfaces(cutter::FaceCutter,d::Integer,n::Integer) = cutter.dface_to_nfaces[d+1,n+1]

@generated function distance(m::CellSubMesh{D},d::Int,i::Int,p::Point{D}) where D
  body = ""
  for d in 0:D
    body *= "if d == $d \n"
    body *= "  f$d = get_dface(m,Val{$d}(),i) \n"
    body *= "  distance(p,f$d) \n"
    body *= "else"
  end
  error = "  throw(ArgumentError(\"\$d-face does not exist\"))\n"
  str = body * '\n' * error * "end"
  Meta.parse(str)
end

function initialize!(m::CellSubMesh{D,T},c::Hexahedron{N,D,T}) where {N,D,T}
  resize!(m.vertex_coordinates,N)
  for (i,v) in enumerate( get_vertices(c) )
    m.vertex_coordinates[i] = v
  end
  for d in 0:D, n in 0:d
    df_to_nf = hex_to_tet_nface_to_mface[D][d+1][n+1]
    resize!( m.dface_to_nfaces[d+1,n+1], 0 )
    push!( m.dface_to_nfaces[d+1,n+1], df_to_nf )
  end
  for d in 1:D
    resize!( m.dface_to_new_dfaces[d], num_dfaces(m,d) )
    for i in 1:num_dfaces(m,d)
      resize!( m.dface_to_new_dfaces[d][i],0)
    end
  end
  m
end

initialize!(m::CellSubMesh,b::BoundingBox) = initialize!(m,Hexahedron(b))

function add_vertex!(m::CellSubMesh{D,T},p::Point{D,T}) where {D,T}
  v = get_vertex_coordinates(m)
  push!( v, p )
end

function delete_dface!(m::CellSubMesh,d::Integer,i::Int)
  for n in 0:d
    df_to_nf = dface_to_nfaces(m,d,n)
    remove!(df_to_nf,i)
  end
end

function load_tables!(cutter::FaceCutter,dim::Integer,ldim::Integer,lid::Integer)
  case = cut_tet_dface_to_case[dim][ldim][lid]
  tables_nf_to_nF = cut_tet_nsubface_to_nface[dim][case]
  tables_nf_to_mf = cut_tet_nface_to_mface[dim][case]
  resize!(cutter.dfaces,0)
  push!(cutter.dfaces,tables_nf_to_nF)
  for d in 0:dim, n in 0:d
    resize!( cutter.dface_to_nfaces[d+1,n+1], 0 )
    push!( cutter.dface_to_nfaces[d+1,n+1], tables_nf_to_mf[d+1][n+1] )
  end
end

function setup_dfaces!(
  cutter::FaceCutter,
  mesh::CellSubMesh,
  id::Integer,
  ldim::Integer,
  lid::Integer)

  dim = num_dims(cutter) 
  for n in 0:dim
    df_to_nf = dface_to_nfaces(mesh,dim,n)
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

function _find_coincident_dface(m::CellSubMesh,c::FaceCutter,dim::Integer,id::Integer,ldim::Integer,ldf::Integer)
  lid = get_dface(c,ldim,ldf)
  df_id = dface_to_nfaces(m,dim,ldim)[id,lid]
  df_to_v = dface_to_nfaces(m,ldim,0)
  ldf_to_lv = dface_to_nfaces(c,ldim,0)
  for new_df in m.dface_to_new_dfaces[ldim][df_id]
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
        
function setup_connectivities!(cutter::FaceCutter)
  dim = num_dims(cutter)
  for d in 0:dim, n in 0:d
    df_to_nf = dface_to_nfaces(cutter,d,n)
    for i in 1:length(df_to_nf), j in 1:length(df_to_nf,i)
      df_to_nf[i,j] = get_dface(cutter,n,df_to_nf[i,j])
    end
  end
  cutter
end

function remove_repeated_faces!(cutter::FaceCutter,mesh::CellSubMesh)
  dim = num_dims(cutter)
  for d in 0:dim-1, n in 0:d
    df_to_nf = dface_to_nfaces(cutter,d,n)
    for i in 1:num_dfaces(cutter,d)
      if get_dface(cutter,d,i) <= num_dfaces(mesh,d)
        remove!(df_to_nf,i)
      end
    end
  end
  cutter
end

function update_dface_to_new_dfaces!(mesh::CellSubMesh,cutter::FaceCutter,id::Integer)
  dim = num_dims(cutter)
  for d in 1:dim
    df_to_new_df = mesh.dface_to_new_dfaces[d]
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
    push!( mesh.dface_to_new_dfaces[dim][id], get_dface(cutter,dim,i) )
  end
  mesh
end

function add_faces!(mesh::CellSubMesh,cutter::FaceCutter)
  dim = num_dims(cutter)
  for d in 0:dim, n in 0:d
    push!( dface_to_nfaces(mesh,d,n), dface_to_nfaces(cutter,d,n) )
  end
  mesh
end

function refine!(cutter::FaceCutter,mesh::CellSubMesh,dim::Integer,id::Integer,ldim::Integer,lid::Integer)
  load_tables!(cutter,dim,ldim,lid)
  setup_dfaces!(cutter,mesh,id,ldim,lid)
  setup_connectivities!(cutter)
  remove_repeated_faces!(cutter,mesh)
  cutter
end


function update!(mesh::CellSubMesh,cutter::FaceCutter,id::Integer)
  dim = num_dims(cutter)
  update_dface_to_new_dfaces!(mesh,cutter,id)
  add_faces!(mesh,cutter)
  delete_dface!(mesh,dim,id)
  mesh
end

function refine!(mesh::CellSubMesh{D},cutter::FaceCutter,dim::Integer,id::Integer) where D
  if dim == 0
    return mesh
  end
  for d in dim:D
    df_to_nf = dface_to_nfaces(mesh,d,dim)
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

function writevtk(m::CellSubMesh{D,T},file_base_name) where {D,T}
  d_to_vtk_type_id = Dict(0=>1,1=>3,2=>5,3=>10)
  num_points = num_vertices(m)
  points = zeros(T,D,num_points)
  for (i ,p ) in enumerate(get_vertex_coordinates(m)), d in 1:D
    points[d,i] = p[d]
  end
  cells = MeshCell{Vector{Int64}}[]
  for d in 0:D
    dface_to_vertices = dface_vertices(m,d)
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

