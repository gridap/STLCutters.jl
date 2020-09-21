
module CellMeshes

using Gridap
using Gridap.Arrays
using Gridap.Geometry


include("tables/LookupCutTables.jl")

const FACE_UNDEF = 0
const FACE_IN = 1
const FACE_OUT = 2
const FACE_BOUNDARY = 3
const FACE_CUT = -1


struct Plane{D,T}
  origin::Point{D,T}
  normal::VectorValue{D,T}
end

struct BoundingBox{D,T}
  pmin::Point{D,T}
  pmax::Point{D,T}
end

@generated function get_vertices(b::BoundingBox{D,T}) where {D,T}
  N = 2^D
  d = Dict( 0 => "pmin", 1 => "pmax" )
  v_str = [ "" for i in 1:N ]
  for i in 1:N
    bin = digits(i-1,base=2,pad=D)
    data = join([ "b.$(d[b])[$i]," for (i,b) in enumerate(bin) ])
    v_str[i] = "Point{$D,$T}($data),"
  end
  vertices = join(v_str)
  str = "($vertices)"
  Meta.parse(str)
end

num_vertices(::T) where T<:BoundingBox = num_vertices(T)

num_vertices(::Type{<:BoundingBox{D}}) where D = 2^D

num_dims(::Type{<:BoundingBox{D}}) where D = D

num_dims(::T) where T<:BoundingBox = num_dims(T)

Base.length(a::Table,i::Integer) = a.ptrs[i+1] - a.ptrs[i]

function Base.setindex!(a::Table,val,i::Integer,j::Integer)
  @boundscheck 0 < j ≤ length(a,i) || throw(BoundsError(a,(i,j)))
  p = a.ptrs[i]-1
  a.data[p+j] = val
end

function Base.getindex(a::Table,i::Integer,j::Integer)
  @boundscheck 0 < j ≤ length(a,i) || throw(BoundsError(a,(i,j)))
  p = a.ptrs[i]-1
  a.data[p+j]
end

function Base.push!(a::Table,b::Vector)
  append!(a.data,b)
  push!(a.ptrs,length(a.data)+1)
  a
end

function Base.append!(a::Table,b::Matrix)
  n = size(b,2)
  m = size(b,1)
  for i in 1:n
    for j in 1:m
      push!(a.data,b[j,i])
    end
    push!(a.ptrs,length(a.data)+1)
  end
  a
end

function Base.append!(a::Table,b::Table)
  for i in 1:length(b)
    for j in 1:length(b,i)
      push!(a.data,b[i,j])
    end
    push!(a.ptrs,length(a.data)+1)
  end
  a
end

function Base.append!(a::Table,b::Vector{<:Vector})
  for i in 1:length(b)
    push!(a,b[i])
  end
  a
end

function Base.empty!(a::Table)
  resize!(a.data,0)
  resize!(a.ptrs,1)
  a.ptrs[1] = 1
end

function Base.fill!(a::Table,val)
  fill!(a.data,val)
  a
end

function Base.resize!(a::Table,len::Vector)
  n = length(a)
  resize!(a.ptrs,1)
  append!(a.ptrs,len)
  length_to_ptrs!(a.ptrs)
  resize!(a.data,a.ptrs[end]-1)
  a
end

function Base.deleteat!(a::Table,masks::Vector)
  @assert length(a) == length(masks)
  k = 0
  for i in 1:length(a)
    if !Bool(masks[i])
      for j in 1:length(a,i)
        k +=1
        a.data[k] = a[i,j]
      end
    end
  end
  resize!(a.data,k)
  k = 0
  for i in 1:(length(a.ptrs)-1)
    if !Bool(masks[i])
      k +=1
      n = a.ptrs[i+1] - a.ptrs[i]
      a.ptrs[k+1] = a.ptrs[k] + n
    end
  end
  resize!(a.ptrs,k+1)
end

## CellMeshes

struct CellMeshCache
  d_to_num_dfaces::Vector{Int}
  d_to_inactive_dfaces::Vector{Vector{Bool}}
  d_to_dface_to_new_dfaces::Vector{Vector{Vector{Int}}}
  d_to_dface_to_cell_dface::Vector{Vector{Int}}
  cell_to_lfacet_to_orientation::Table{Int,Vector{Int},Vector{Int32}}
  vertex_to_surface_mesh_faces::Vector{Vector{Int}}
  facet_to_surface_mesh_facet::Vector{Int}
end

struct CutterCache
  d_to_ldface_to_dface::Table{Int,Vector{Int},Vector{Int32}}
  m_n_to_mface_to_nfaces::Matrix{Table{Int,Vector{Int},Vector{Int32}}}
  cell_to_lfacet_to_orientation::Table{Int,Vector{Int},Vector{Int32}}
end

struct LevelSetCache{D,T}
  levelsets::Vector{Plane{D,T}}
  vertex_to_levelset_to_inoutcut::Table{Int8,Vector{Int8},Vector{Int32}}
  levelset_region_to_inout::Dict{Int,Int8}
  vector_cache_int::Vector{Int}
  vector_cache_float::Vector{Float64}
end

struct MeshCache
  cell::CellMeshCache
  cutter::CutterCache
  levelset::LevelSetCache
  vector_cache_1::Vector{Int}
  vector_cache_2::Vector{Int}
  vector_cache_3::Vector{Int}
end

struct CellMesh{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  m_n_to_mface_to_nfaces::Matrix{Table{Int,Vector{Int},Vector{Int32}}}
  d_to_dface_to_in_out_boundary::Table{Int8,Vector{Int8},Vector{Int32}}
  cache::MeshCache
end

function CellMesh(box::BoundingBox{D,T}) where {D,T}
  mesh = CellMesh{D,T}()
  reset!(mesh,box)
end

function CellMesh(X::Vector{Point{D,T}},K,p::Polytope{D}) where {D,T}
  mesh = CellMesh{D,T}()
  reset!(mesh,X,K)
end

function CellMesh{D,T}() where {D,T}
  coordinates = Point{D,T}[]
  m_n_to_mf_to_nf = Matrix{Table{Int,Vector{Int},Vector{Int32}}}(undef,D+1,D+1)
  d_df_to_io = empty_table(Int,Int32,0)
  for d in 0:D, n in 0:D
    m_n_to_mf_to_nf[d+1,n+1] = empty_table(Int,Int32,0)
  end
  cache = MeshCache(D)
  CellMesh{D,T}( coordinates, m_n_to_mf_to_nf, d_df_to_io, cache )
end

function reset!(m::CellMesh{D,T},box::BoundingBox{D,T}) where {D,T} 
  resize!(m.vertex_coordinates,num_vertices(box))
  for (i,v) in enumerate( get_vertices(box) )
    m.vertex_coordinates[i] = v
  end
  _reset!(m)
end

function reset!(m::CellMesh,X::Vector{<:Point},K)
  resize!(m.vertex_coordinates,length(K))
  for (i,v) in enumerate(K)
    m.vertex_coordinates[i] = X[v]
  end
  _reset!(m)
end

function _reset!(m::CellMesh{D,T}) where {D,T} 
  table_m_n_to_mface_to_nfaces = D_to_m_n_to_mface_to_nfaces_for_hexD_to_tetD[D]
  for d in 0:D, n in 0:d
    df_to_nf = table_m_n_to_mface_to_nfaces[d+1,n+1]
    empty!( get_faces(m,d,n) )
    append!( get_faces(m,d,n), df_to_nf )
  end
  empty!(m.d_to_dface_to_in_out_boundary)

  reset_d_to_dface_to_in_out_boundary!(m,m.cache)
  reset_cache!(m.cache,m)
  m
end

function CellMeshCache(D::Integer)
  d_num_df = [ 0 for d in 0:D ]
  d_to_active_df = [ Bool[] for d in 0:D ]
  d_to_df_to_new_df = [ Vector{Int}[] for d in 1:D ]
  d_to_df_to_cdf = [ Int[] for d in 0:D-1 ]
  c_lf_to_o = empty_table(Int,Int32,0)
  v_to_smf = Vector{Int}[]
  f_to_sm_f = Int[]
  CellMeshCache( d_num_df, d_to_active_df, d_to_df_to_new_df, d_to_df_to_cdf, c_lf_to_o, v_to_smf, f_to_sm_f )
end

function CutterCache(D::Integer)
  d_to_ldf_to_df = empty_table(Int,Int32,0)
  m_n_to_mf_to_nf = Matrix(undef,D+1,D+1)
  c_ldf_to_o = empty_table(Int,Int32,0)
  for i in 1:D+1, j in 1:D+1
    m_n_to_mf_to_nf[i,j] = empty_table(Int,Int32,0)
  end
  CutterCache( d_to_ldf_to_df, m_n_to_mf_to_nf, c_ldf_to_o )
end

function MeshCache(D::Integer)
  cell_cache = CellMeshCache(D) 
  cutter_cache = CutterCache(D)
  levelset_cache = LevelSetCache{D,Float64}()
  vector1 = []
  vector2 = []
  vector3 = []
  MeshCache( cell_cache, cutter_cache, levelset_cache, vector1, vector2, vector3 )
end

function reset_cell_cache!(cache::CellMeshCache,mesh::CellMesh)
  D = num_dims(mesh)
  for d in 0:D
    set_num_dfaces!(cache,d, num_dfaces(mesh,d) )
    resize!(cache.d_to_inactive_dfaces[d+1], num_dfaces(mesh,d) )
    fill!(cache.d_to_inactive_dfaces[d+1],false)
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
  empty!( cache.cell_to_lfacet_to_orientation )
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
    empty!( cache.d_to_dface_to_cell_dface[d+1] )
    append!( cache.d_to_dface_to_cell_dface[d+1], df_to_cdf )
  end
  cache
end

function reset_cutter_cache!(cutter::CutterCache)
 empty!(cutter.d_to_ldface_to_dface)
 empty!(cutter.cell_to_lfacet_to_orientation)
  for t in cutter.m_n_to_mface_to_nfaces
    empty!(t)
  end
  cutter
end

function reset_cache!(cache::MeshCache,mesh::CellMesh)
  reset_cell_cache!(cache.cell,mesh)
  reset_cutter_cache!(cache.cutter)
  cache
end

get_cache(mesh::CellMesh) = mesh.cache

get_levelset_cache(mesh::CellMesh) = mesh.cache.levelset

function LevelSetCache{D,T}() where {D,T}
  ls = Plane{D,T}[]
  v_to_ls_to_ioc = empty_table(Int8,Int32,0)
  ls_region_to_io = Dict{Int,Int8}()
  int_cache = Int[]
  float_cache = Float64[]
  LevelSetCache(
    ls,
    v_to_ls_to_ioc,
    ls_region_to_io,
    int_cache,
    float_cache)
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

get_vector_cache!(mesh::CellMesh) = empty!(get_cache(mesh).vector_cache_1)

get_vector!(a::MeshCache) = resize!(a.vector_cache_2,0)

get_vector_bis!(a::MeshCache) = resize!(a.vector_cache_3,0)

num_dims(cache::CellMeshCache) = length(cache.d_to_dface_to_new_dfaces)

num_dfaces(cache::CellMeshCache,d::Integer) = cache.d_to_num_dfaces[d+1]

set_num_dfaces!(cache::CellMeshCache,d::Integer,val) = cache.d_to_num_dfaces[d+1] = val

get_dface_to_new_dfaces(a::CellMeshCache,d::Integer) = a.d_to_dface_to_new_dfaces[d]

get_new_faces(a::CellMeshCache,d::Integer,i::Integer) = a.d_to_dface_to_new_dfaces[d][i]

get_cell_dface(a::CellMeshCache,d::Integer,i::Integer) = a.d_to_dface_to_cell_dface[d+1][i]

@inline function isactive(m::CellMesh,d::Integer,i::Integer)
  !m.cache.cell.d_to_inactive_dfaces[d+1][i]
end

num_dims(mesh::CellMesh{D}) where D = D

@inline function num_vertices(m::CellMesh) 
  length(m.vertex_coordinates)
end

@inline function num_dfaces(m::CellMesh,d::Integer)
  if d == 0
    num_vertices(m)
  else
    length( get_faces(m,d,0) ) 
  end
end

@inline num_faces(m::CellMesh,d::Integer) = num_dfaces(m,d)

function num_faces(m::CellMesh{D}) where D
  n = 0
  for d in 0:D
    n += num_faces(m,d)
  end
  n
end

num_edges(m::CellMesh) = num_dfaces(m,1)

num_facets(m::CellMesh{D}) where D = num_dfaces(m,D-1)

num_cells(m::CellMesh{D}) where D = num_dfaces(m,D)

@inline function get_faces(m::CellMesh,d::Integer,n::Integer)
  m.m_n_to_mface_to_nfaces[d+1,n+1]
end

@inline function get_dface_to_vertices(m::CellMesh,d::Integer)
  get_faces(m,d,0)
end

@inline function get_cell_to_vertices(m::CellMesh{D}) where D
  get_faces(m,D,0)
end

@inline function get_cell_to_facets(m::CellMesh{D}) where D
  get_faces(m,D,D-1)
end

@inline function get_facet_to_vertices(m::CellMesh{D}) where D
  get_faces(m,D-1,0)
end

@inline function get_edge_to_vertices(m::CellMesh)
  get_faces(m,1,0)
end

@inline function get_vertex_coordinates(m::CellMesh) 
  m.vertex_coordinates
end

@inline function get_vertex_coordinates(m::CellMesh,i::Integer) 
  m.vertex_coordinates[i]
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

get_orientation(a::CellMeshCache,cell::Integer,lfacet::Integer) = a.cell_to_lfacet_to_orientation[cell,lfacet]

function  set_orientation!(a::CellMeshCache,cell::Integer,lfacet::Integer,val)
  a.cell_to_lfacet_to_orientation[cell][lfacet] = val
end

function set_orientation!(cutter::CutterCache,cell::Integer,lfacet::Integer, val )
  cutter.cell_to_lfacet_to_orientation[cell,lfacet] = val
end

function add_vertex!(mesh::CellMesh,d::Integer,face::Integer,point::Point,sm_face::Integer)
  if face != UNSET
    add_surface_mesh_face_to_vertex!(mesh.cache,d,face,sm_face)
    _add_vertex!(mesh,d,face,point)
  else
    UNSET
  end
end

function _add_vertex!(m::CellMesh{D,T},p::Point{D,T}) where {D,T}
  v = get_vertex_coordinates(m)
  push!( v, p )
  num_vertices(m)
end

function remove_dface!(m::CellMesh,d::Integer,i::Int)
  m.cache.cell.d_to_inactive_dfaces[d+1][i] = true
end

function _add_vertex!(mesh::CellMesh,dim::Integer,face::Integer,point::Point)
  if dim == 0
    return face
  end
  cache = mesh.cache
  @assert isactive(mesh,dim,face)
  for d in dim:num_dims(mesh)
    df_to_nf = get_faces(mesh,d,dim)
    for iface in 1:length(df_to_nf)
      if isactive(mesh,d,iface)
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
  for d in 0:dim
    masks = get_masks(mesh,d)
    cache = get_vector_cache!(mesh)
    resize!(cache,length(get_faces(cutter,d,d)))
    fill!(cache,false)
    append!(masks,cache)
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
  update_dface_to_new_dfaces!(cell_cache,cutter_cache,face) # allocations
  update_cell_to_lface_to_orientation!(cell_cache,cutter_cache)
  update_dface_to_cell_dface!(cell_cache,cutter_cache,face) # casual allocations
  caches
end

function load_tables!(cutter::CutterCache,dim::Integer,ldim::Integer,lid::Integer)
  case = D_to_d_to_dface_to_case_for_cut_tetD[dim][ldim][lid]
  tables_nf_to_nF = D_to_case_to_d_to_subdface_to_dface_for_cut_tetD[dim][case]
  tables_nf_to_mf = D_to_case_to_m_n_to_mface_to_nfaces_for_cut_tetD[dim][case]
  tables_c_to_lf_to_o = D_to_case_to_cell_lfacet_to_orientation_for_cut_tetD[dim][case]

  empty!( cutter.d_to_ldface_to_dface )
  append!( cutter.d_to_ldface_to_dface, tables_nf_to_nF )
  for d in 0:dim, n in 0:d
    empty!( get_faces(cutter,d,n) )
    append!( get_faces(cutter,d,n), tables_nf_to_mf[d+1,n+1] )
  end
  empty!( cutter.cell_to_lfacet_to_orientation )
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
            @assert !isactive(mesh,_d,_face)
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
            @assert _dface != UNSET
            set_dface!(cutter,d,tdface, _dface)
          else
            n_dfaces += 1
            set_dface!(cutter,d,tdface, n_dfaces)
          end
        end
      else
        dface = dim_f_to_df[dim_face,get_dface(cutter,d,tdface)]
        if !isactive(mesh,d,dface)
          @assert d ≥ cut_dim
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
  masks = get_vector_cache!(mesh)
  dim = num_dims(cutter)
  for d in 0:dim-1, n in 0:d
    df_to_nf = get_faces(cutter,d,n)
    resize!(masks,num_dfaces(cutter,d))
    fill!(masks,false)
    for i in 1:num_dfaces(cutter,d)
      if get_dface(cutter,d,i) <= length( get_faces(mesh,d,n) )
        masks[i] = true
      end
    end
    deleteat!(df_to_nf,masks)
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

function update_dface_to_cell_dface!(cache::CellMeshCache,cutter::CutterCache,face::Integer)
  dim = num_dims(cutter)
  D = num_dims(cache)
  for d in 0:dim-1
    df_to_cdf = cache.d_to_dface_to_cell_dface[d+1]
    for i in 1:num_dfaces(cutter,d)
      if get_dface(cutter,d,i) > length(df_to_cdf)
        @assert  get_dface(cutter,d,i) == length(df_to_cdf) + 1
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

function get_masks(cache::CellMeshCache,d::Integer)
  cache.d_to_inactive_dfaces[d+1]
end

function get_masks(mesh::CellMesh,d::Integer)
  get_masks(mesh.cache.cell,d)
end

function compact_mesh!(mesh::CellMesh)
  D = num_dims(mesh)
  for d in 0:D
    n_dfaces = 0
    df_to_df = get_faces(mesh,d,d)
    for i in 1:length(df_to_df)
      if isactive(mesh,d,i)
        @assert length(df_to_df,i) == 1
        n_dfaces += 1
        df_to_df[i,1] = n_dfaces
      end
    end
  end
  for d in 1:D, n in 0:d-1
    nf_to_nf = get_faces(mesh,n,n)
    df_to_nf = get_faces(mesh,d,n)
    for i in 1:length(df_to_nf)
      if isactive(mesh,d,i)
        for j in 1:length(df_to_nf,i)
          df_to_nf[i,j] = nf_to_nf[ df_to_nf[i,j], 1 ]
        end
      end
    end
  end
  for d in 0:D
    for n in 0:d
      deleteat!( get_faces(mesh,d,n), get_masks(mesh,d) )
    end
    deleteat!( get_masks(mesh,d), get_masks(mesh,d) )
  end
  d_to_dfaces_to_ioc = mesh.d_to_dface_to_in_out_boundary
  empty!(d_to_dfaces_to_ioc)
  dfaces_to_ioc = get_vector_cache!(mesh)
  for d in 0:D
    resize!(dfaces_to_ioc,num_faces(mesh,d))
    fill!(dfaces_to_ioc,UNSET)
    push!(d_to_dfaces_to_ioc,dfaces_to_ioc)
  end
  mesh
end

function cut_levelsets!(cell_mesh::CellMesh)
  atol = 1e-13
  ls_cache = get_levelset_cache(cell_mesh)
  distances = get_float_vector(ls_cache)
  levelsets = get_levelsets(ls_cache)
  for (ls,plane) in enumerate(levelsets)
    resize!(distances,num_vertices(cell_mesh))
    for vertex in 1:num_vertices(cell_mesh)
      point = get_vertex_coordinates(cell_mesh,vertex)
      dist = signed_distance(point,plane)
      if abs(dist) < atol
        dist = 0.0
      end
      distances[vertex] = dist
      set_vertex_ioc_to_ls!(ls_cache,vertex,ls,dist)
    end
    for edge in 1:num_edges(cell_mesh)
      if isactive(cell_mesh,1,edge)
        if is_edge_cut(cell_mesh,edge,distances)
          point = edge_intersection(cell_mesh,edge,distances)
          add_vertex_to_edge!(cell_mesh,edge,point,ls)
        end
      end
    end
  end
  cell_mesh
end

function add_vertex_to_edge!(mesh::CellMesh,edge::Integer,point::Point,ls::Integer)
  _add_vertex!(mesh,1,edge,point)
  add_vertex_ioc_to_ls!(mesh,edge,ls)
end

function set_vertex_ioc_to_ls!(ls_cache::LevelSetCache,vertex::Integer,ls::Integer,dist::Float64)
  v_to_ls_to_ioc = ls_cache.vertex_to_levelset_to_inoutcut
  if dist == 0
    v_to_ls_to_ioc[vertex,ls] = FACE_CUT
  elseif dist > 0
    v_to_ls_to_ioc[vertex,ls] = FACE_IN
  else # dist < 0
    v_to_ls_to_ioc[vertex,ls] = FACE_OUT
  end
end

function add_vertex_ioc_to_ls!(mesh::CellMesh,edge::Integer,ls::Integer)
  ls_cache = get_levelset_cache(mesh)
  e_to_v = get_edge_to_vertices(mesh)
  v_to_ls_to_ioc = ls_cache.vertex_to_levelset_to_inoutcut
  ls_to_ioc = get_int_vector(ls_cache)
  resize!(ls_to_ioc,num_levelsets(ls_cache))
  fill!(ls_to_ioc,UNSET)
  v1 = e_to_v[edge,1]
  v2 = e_to_v[edge,2]
  for i in 1:ls-1
    _ioc = UNSET
    ioc1 = v_to_ls_to_ioc[v1,i]
    ioc2 = v_to_ls_to_ioc[v2,i]
    if ioc1 == FACE_CUT
      _ioc = ioc2
    elseif ioc2 == FACE_CUT
      _ioc = ioc1
    elseif ioc1 ≠ ioc2
      _ioc = FACE_CUT
    else # ioc1 == ioc2
      _ioc = ioc1
    end
    ls_to_ioc[i] = _ioc
  end
  ls_to_ioc[ls] = FACE_CUT
  push!(v_to_ls_to_ioc,ls_to_ioc)
end


function is_edge_cut(mesh::CellMesh,edge::Integer,distances)
  e_to_v = get_edge_to_vertices(mesh)
  v1 = e_to_v[edge,1]
  v2 = e_to_v[edge,2]
  dist1 = distances[v1]
  dist2 = distances[v2]
  if dist1 == 0 || dist2 == 0
    false
  elseif sign(dist1) ≠ sign(dist2)
    true
  else
    false
  end
end

function edge_intersection(mesh::CellMesh,edge::Integer,distances)
  @assert is_edge_cut(mesh,edge,distances)
  e_to_v = get_edge_to_vertices(mesh)
  v1 = e_to_v[edge,1]
  v2 = e_to_v[edge,2]
  dist1 = distances[v1]
  dist2 = distances[v2]
  len = abs(dist1) + abs(dist2)
  p1 = get_vertex_coordinates(mesh,v1)
  p2 = get_vertex_coordinates(mesh,v2)
  (abs(dist2)/len)*p1 + (abs(dist1)/len)*p2
end

function signed_distance(p::Point,plane::Plane)
  o = plane.origin
  n = plane.normal
  (p-o)⋅n
end

function get_in_region(cache::LevelSetCache,vertex::Integer)
  r = 0
  v_to_ls_to_ioc = cache.vertex_to_levelset_to_inoutcut
  for ls in 1:num_levelsets(cache)
    ioc = v_to_ls_to_ioc[vertex,ls]
    if ioc ∈ (FACE_IN,FACE_CUT)
      r |= 1<<(ls-1)
    end
  end
  r + 1
end

function get_out_region(cache::LevelSetCache,vertex::Integer)
  r = 0
  v_to_ls_to_ioc = cache.vertex_to_levelset_to_inoutcut
  for ls in 1:num_levelsets(cache)
    ioc = v_to_ls_to_ioc[vertex,ls]
    if ioc == FACE_IN
      r |= 1<<(ls-1)
    end
  end
  r + 1
end

function is_vertex_cut(cache::LevelSetCache,vertex::Integer)
  v_to_ls_to_ioc = cache.vertex_to_levelset_to_inoutcut
  for ls in 1:num_levelsets(cache)
    ioc = v_to_ls_to_ioc[vertex,ls]
    if ioc == FACE_CUT
      return true
    end
  end
  false
end

function define_regions!(mesh::CellMesh)
  ls_cache = get_levelset_cache(mesh) 
  for vertex in 1:num_vertices(mesh)
    if is_vertex_cut(ls_cache,vertex)
      r_in = get_in_region(ls_cache,vertex)
      r_out = get_out_region(ls_cache,vertex)
      set_region_in!(ls_cache,r_in)
      set_region_out!(ls_cache,r_out)
    end
  end
  ls_cache.levelset_region_to_inout
end

function set_region_in!(cache::LevelSetCache,r::Integer)
  cache.levelset_region_to_inout[r] = FACE_IN
end

function set_region_out!(cache::LevelSetCache,r::Integer)
  cache.levelset_region_to_inout[r] = FACE_OUT
end


function is_face_on_boudary(mesh::CellMesh,d::Integer,dface::Integer)
  df_to_v = get_dface_to_vertices(mesh,d)
  ls_cache = get_levelset_cache(mesh)
  v_to_ls_to_ioc = ls_cache.vertex_to_levelset_to_inoutcut
  for ls in 1:num_levelsets(ls_cache)
    ioc = FACE_CUT
    for lv in 1:length(df_to_v,dface)
      vertex = df_to_v[dface,lv]
      if v_to_ls_to_ioc[vertex,ls] ≠ FACE_CUT
        @assert ioc == FACE_CUT || ioc == v_to_ls_to_ioc[vertex,ls]
        ioc = v_to_ls_to_ioc[vertex,ls]
      end
    end
    if ioc == FACE_CUT
      return true
    end
  end
  false
end

function get_face_region(mesh::CellMesh,d::Integer,dface::Integer)
  r = 0
  df_to_v = get_dface_to_vertices(mesh,d)
  ls_cache = get_levelset_cache(mesh)
  v_to_ls_to_ioc = ls_cache.vertex_to_levelset_to_inoutcut
  for ls in 1:num_levelsets(ls_cache)
    ioc = FACE_CUT
    for lv in 1:length(df_to_v,dface)
      vertex = df_to_v[dface,lv]
      if v_to_ls_to_ioc[vertex,ls] ≠ FACE_CUT
        @assert ioc == FACE_CUT || ioc == v_to_ls_to_ioc[vertex,ls]
        ioc = v_to_ls_to_ioc[vertex,ls]
      end
    end
    if ioc == FACE_IN
      r |= 1<<(ls-1)
    end
  end
  r + 1
end

get_cell_region(mesh::CellMesh{D},cell) where D = get_face_region(mesh,D,cell)

function define_faces!(mesh::CellMesh)
  D = num_dims(mesh)
  ls_cache = get_levelset_cache(mesh)
  d_to_dfaces_to_ioc = mesh.d_to_dface_to_in_out_boundary
  r_to_io = ls_cache.levelset_region_to_inout
  for d in 0:D
    for dface in 1:num_faces(mesh,d)
      if is_face_on_boudary(mesh,d,dface)
        d_to_dfaces_to_ioc[d+1,dface] = FACE_CUT
      else
        r = get_face_region(mesh,d,dface)
        d_to_dfaces_to_ioc[d+1,dface] = r_to_io[r]
      end
    end
  end
end

function reset!(cache::LevelSetCache,mesh::CellMesh,levelsets)
  empty!(cache.levelsets)
  empty!(cache.vertex_to_levelset_to_inoutcut)
  empty!(cache.levelset_region_to_inout)
  append!(cache.levelsets,levelsets)
  ls_to_ioc = get_int_vector(cache)
  resize!(ls_to_ioc,length(levelsets))
  fill!(ls_to_ioc,UNSET)
  for vertex in 1:num_vertices(mesh)
    push!(cache.vertex_to_levelset_to_inoutcut,ls_to_ioc)
  end
  cache
end

get_int_vector(cache::LevelSetCache) = cache.vector_cache_int

get_float_vector(cache::LevelSetCache) = cache.vector_cache_float

get_levelsets(cache::LevelSetCache) = cache.levelsets

num_levelsets(cache::LevelSetCache) = length(cache.levelsets)

function get_cell_to_inout(mesh::CellMesh{D}) where D
  mesh.d_to_dface_to_in_out_boundary[D+1]
end

function get_cell_to_inout(mesh::CellMesh{D},cell) where D
  mesh.d_to_dface_to_in_out_boundary[D+1,cell]
end

function get_facet_to_inout(mesh::CellMesh{D}) where D
  mesh.d_to_dface_to_in_out_boundary[D]
end

function get_facet_to_inout(mesh::CellMesh{D},facet) where D
  mesh.d_to_dface_to_in_out_boundary[D,facet]
end

function compute_cell_mesh!(mesh::CellMesh)
  cut_levelsets!(mesh)
  compact_mesh!(mesh)
  define_regions!(mesh)
  define_faces!(mesh)
  mesh
end

function compute_cell_mesh!(mesh::CellMesh,levelsets::Vector{<:Plane})
  ls_cache = mesh.cache.levelset
  reset!(ls_cache,mesh,levelsets)
  compute_cell_mesh!(mesh)
end

end # module
