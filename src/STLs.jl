
struct STLTopology{Dc,Dp,T} <: GridTopology{Dc,Dp}
  vertex_to_coordinates::Vector{Point{Dp,T}}
  n_m_to_nface_to_mfaces::Matrix{Table{Int32,Vector{Int32},Vector{Int32}}}
  polytope::Polytope{Dc}
  facedims::Vector{Int8}
  offsets::Vector{Int32}
end

const STL = STLTopology

function STL(model::DiscreteModel)
  topo = get_grid_topology(model)
  STL(topo)
end

function STL(topo::GridTopology{Dc,Dp}) where {Dc,Dp}
  Dc == Dp-1 || error("STL must be a surface")
  v_coords = get_vertex_coordinates(topo)
  polytope = only(get_polytopes(topo))
  facedims = get_facedims(topo)
  offsets = Int32.(get_offsets(topo))
  T = Table{Int32,Vector{Int32},Vector{Int32}}
  n_m_to_nf_to_mf = Matrix{T}(undef,Dc+1,Dc+1)
  for n in 0:Dc, m in 0:Dc
    n_m_to_nf_to_mf[n+1,m+1] = get_faces(topo,n,m)
  end
  STL(
    v_coords,
    n_m_to_nf_to_mf,
    polytope,
    facedims,
    offsets)
end

num_dims(stl::STL{Dc}) where Dc = Dc

num_point_dims(stl::STL{Dc,Dp}) where {Dc,Dp} = Dp

get_faces(stl::STL,n::Integer,m::Integer) = stl.n_m_to_nface_to_mfaces[n+1,m+1]

get_vertex_coordinates(stl::STL) = stl.vertex_to_coordinates

get_polytope(stl::STL) = stl.polytope

get_polytopes(stl::STL) = Fill(get_polytope(stl),1)

get_offsets(stl::STL) = stl.offsets

get_offset(stl::STL,d::Integer) = get_offsets(stl)[d+1]

get_facedims(stl::STL) = stl.facedims

get_face_vertices(stl::STL,d::Integer) = get_faces(stl,d,0)

get_cell_vertices(stl::STL{Dc}) where Dc = get_face_vertices(stl,Dc)

num_faces(stl::STL,d::Integer) = length(get_faces(stl,d,d))

num_faces(stl::STL) = length(get_facedims(stl))

num_vertices(stl::STL) = num_faces(stl,0)

num_edges(stl::STL) = num_faces(stl,0)

num_facets(stl::STL{D}) where D = num_faces(stl,D-1)

num_cells(stl::STL{Dc}) where Dc = num_faces(stl,Dc)

function get_face_vertices(stl::STL{D}) where D
  @notimplementedif D ≤ 1
  f_to_v = lazy_append(get_face_vertices(stl,0),get_face_vertices(stl,1))
  for d in 2:D
    f_to_v = lazy_append(f_to_v,get_face_vertices(stl,d))
  end
  f_to_v
end

function get_cell!(c,stl::STL{Dc},i::Integer) where Dc
  get_dface!(c,stl,i,Val{Dc}())
end

function get_dface!(c,stl::STL,i::Integer,::Val{d}) where d
  T = get_face_vertices(stl,d)
  X = get_vertex_coordinates(stl)
  p = get_polytope(stl)
  get_dface!(c,T,X,p,i,Val{d}())
end

function get_cell(stl::STL,i::Integer)
  c = get_cell_cache(stl)
  get_cell!(c,stl,i)
end

function get_cell_cache(stl::STL{Dc}) where Dc
  get_dface_cache(stl,Dc)
end

function get_dface_cache(stl::STL,d::Integer)
  T = get_face_vertices(stl,d)
  array_cache(T)
end

function get_cell!(cache,T,X,p::Polytope{D},i::Integer) where D
  get_dface!(cache,T,X,p,i,Val{D}())
end

function get_dface!(cache,T,X,p::Polytope,i::Integer,::Val{d}) where d
  @notimplementedif !is_simplex(p)
  get_simplex_dface!(cache,T,X,i,Val{d}())
end

function get_simplex_dface!(cache,T,X,i::Integer,::Val{d}) where d
  nodes = getindex!(cache,T,i)
  vertices = ntuple( i -> X[nodes[i]], Val{d+1}() )
  simplex_face( vertices )
end

function is_open_surface(stl::GridTopology)
  Dc = num_dims(stl)
  e_to_f = get_faces(stl,Dc-1,Dc)
  1 ≤ maximum(length,e_to_f) ≤ 2 || return false
  1 ≤ minimum(length,e_to_f) ≤ 2 || return false
  true
end

# Read STL from file

function read_stl(filename::String)
  mesh = _load(filename)
  vertex_to_coordinates = MeshIO.decompose(MeshIO.Point3,mesh)
  facet_to_vertices = MeshIO.decompose(MeshIO.TriangleFace,mesh)
  vertex_to_normals = MeshIO.decompose_normals(mesh)
  vertex_to_coordinates = Vector{Point{3,Float64}}( vertex_to_coordinates )
  facet_to_vertices = Table( Vector{Vector{Int}}( facet_to_vertices ) )
  vertex_to_normals = Vector{VectorValue{3,Float64}}( vertex_to_normals )
  facet_to_normals = vertex_to_normals[1:3:end]
  vertex_to_coordinates,facet_to_vertices,facet_to_normals
end

function _load(filename)
  if Sys.islinux()
    F = _file_format(filename)
    load(File{F}(filename))
  else
    load(filename)
  end
end

function _file_format(filename)
  ext =  lowercase(last(splitext(filename)))
  stl,obj = ".stl",".obj"
  ext == stl && !_is_binary(filename) && return format"STL_ASCII"
  ext == stl && _is_binary(filename) && return format"STL_BINARY"
  ext == obj && return format"OBJ"
  msg = "Invalid file format: $filename"
  @unreachable msg
end

function _get_encoding(filename)
  @notimplementedif !Sys.islinux()
  cmd = `file --mime-encoding $filename`
  read(cmd,String)
end

_is_binary(filename) = contains(_get_encoding(filename),"binary")

function Base.convert(::Type{Point{D,T}},x::MeshIO.Point{D}) where {D,T}
  Point{D,T}(x.data)
end

function Base.convert(::Type{Vector{T}},x::MeshIO.TriangleFace) where T
  convert.(T,collect(x))
end

# Save STL to file

function save_as_stl(stl::DiscreteModel{Dc,Dp},filename) where {Dc,Dp}
  (Dc == 2 && Dp == 3 ) || error("Geometry incompatible with STL format")
  filename *= ".stl"
  file = File{format"STL_BINARY"}(filename)
  facetype = MeshIO.GLTriangleFace
  pointtype = MeshIO.Point3f
  normaltype = MeshIO.Vec3f
  faces = Array{facetype}(undef, num_cells(stl))
  vertices = Array{pointtype}(undef, num_cells(stl) * 3)
  normals = Array{normaltype}(undef, num_cells(stl) * 3)
  for i in 1:num_cells(stl)
    faces[i] = facetype( (i-1)*3+1,(i-1)*3+2, (i-1)*3+3 )
    n = normal(get_cell(stl,i))
    normals[ (i-1)*3+1 ] = normaltype(Tuple(n)...)
    normals[ (i-1)*3+2 ] = normals[ (i-1)*3+1 ]
    normals[ (i-1)*3+3 ] = normals[ (i-1)*3+1 ]
    coords = get_cell_coordinates(stl)[i]
    vertices[ (i-1)*3+1 ] = pointtype(Tuple(coords[1])...)
    vertices[ (i-1)*3+2 ] = pointtype(Tuple(coords[2])...)
    vertices[ (i-1)*3+3 ] = pointtype(Tuple(coords[3])...)
  end
  mesh = MeshIO.Mesh(MeshIO.meta(vertices; normals=normals), faces)
  save(file,mesh)
  filename
end

function save_as_stl(stls::Vector{<:DiscreteModel},filename)
  files = []
  for (i,stl) in enumerate(stls)
    f = save_as_stl(stl,filename*"_$i")
    push!(files,f)
  end
  files
end

# STL as Grid and DiscreteModel

function get_cell_cache(grid::Grid)
  T = get_cell_node_ids(grid)
  array_cache(T)
end

function get_cell!(c,grid::Grid,i)
  T = get_cell_node_ids(grid)
  X = get_node_coordinates(grid)
  p = get_polytope( only(get_reffes(grid)) )
  get_cell!(c,T,X,p,i)
end
function get_facet_cache(model::DiscreteModel)
  Dc = num_dims(model)
  get_dface_cache(model,Dc-1)
end

function get_dface_cache(model::DiscreteModel,d)
  topo = get_grid_topology(model)
  T = get_face_vertices(topo,d)
  array_cache(T)
end

function get_facet!(c,model::DiscreteModel{Dc},i) where Dc
  get_dface!(c,model,i,Val{Dc-1}())
end

function get_dface!(c,model::DiscreteModel,i,::Val{d}) where d
  topo = get_grid_topology(model)
  T = get_face_vertices(topo,d)
  X = get_vertex_coordinates(topo)
  p = get_polytope( only(get_reffes(model)) )
  get_dface!(c,T,X,p,i,Val{d}())
end

function get_dface(model::DiscreteModel,i,::Val{d}) where d
  c = get_dface_cache(model,d)
  get_dface!(c,model,i,Val{d}())
end

function get_vertex(stl::DiscreteModel,vertex::Integer)
  get_vertex_coordinates(get_grid_topology(stl))[vertex]
end

function get_edge(stl::DiscreteModel,edge::Integer)
  get_dface(stl,edge,Val{1}())
end

function get_facet(stl::DiscreteModel{Dc},facet::Integer) where Dc
  get_dface(stl,facet,Val{Dc-1}())
end

function get_cell(stl::DiscreteModel{Dc},cell::Integer) where Dc
  get_dface(stl,cell,Val{Dc}())
end

function Base.eps(T::Type{<:AbstractFloat},grid::Grid)
  pmin,pmax = get_bounding_box(grid)
  vmax = max(abs.(Tuple(pmin))...,abs.(Tuple(pmax))...)
  eps(T(vmax))
end

Base.eps(grid::Grid) = eps(Float64,grid)

Base.eps(model::DiscreteModel) = eps(get_grid(model))

function compute_stl_model(
  cell_to_vertices::Table,
  vertex_to_coordinates::Vector{<:Point{D}}) where D

  p = Polytope( tfill(TET_AXIS,Val{D-1}()) )
  compute_model(cell_to_vertices,vertex_to_coordinates,p)
end

function compute_model(
  cell_to_vertices::Table,
  vertex_to_coordinates::Vector{<:Point},
  p::Polytope)

  grid = compute_grid(cell_to_vertices,vertex_to_coordinates,p)
  UnstructuredDiscreteModel(grid)
end

function compute_grid(
  cell_to_nodes::AbstractArray,
  node_to_coordinates::Vector{<:Point},
  p::Polytope)

  T = Table(cell_to_nodes)
  X = node_to_coordinates
  reffes = [LagrangianRefFE(p)]
  cell_types = fill(1,length(T))
  UnstructuredGrid(X,T,reffes,cell_types)
end

function merge_and_collapse(stl::DiscreteModel;atol=eps(Float32,stl))
  stl = merge_nodes(stl;atol)
  collapse_small_facets!(stl;atol)
  stl = merge_nodes(stl;atol)
  stl = delete_duplicated_faces(stl)
  preprocess_small_facets(stl;atol)
end

function collapse_small_facets!(stl::DiscreteModel;atol)
  D = num_dims(stl)
  stl_topo = get_grid_topology(stl)
  coords = get_vertex_coordinates(stl_topo)
  f_to_v = get_faces(stl_topo,D,0)
  f_to_e = get_faces(stl_topo,D,D-1)
  e_to_v = get_faces(stl_topo,D-1,0)
  fv_c = array_cache(f_to_v)
  fe_c = array_cache(f_to_e)
  ev_c = array_cache(e_to_v)
  c = get_facet_cache(stl)
  for f in 1:num_cells(stl)
    for v in getindex!(fv_c,f_to_v,f)
      vertex = coords[v]
      for e in getindex!(fe_c,f_to_e,f)
        vertices = getindex!(ev_c,e_to_v,e)
        v ∉ vertices || continue
        edge = get_facet!(c,stl,e)
        if distance(vertex,edge) < atol
          _vertex = projection(vertex,edge)
          coords[v] = _vertex
        end
      end
    end
  end
  stl
end

function merge_nodes(stl::DiscreteModel;atol=eps(Float32,stl))
  X,T = delete_repeated_vertices(stl;atol)
  compute_stl_model(T,X)
end

function delete_repeated_vertices(stl::DiscreteModel;atol)
  vertex_to_equal = _map_equal_vertices(stl;atol)
  topo = get_grid_topology(stl)
  X = get_vertex_coordinates(topo)
  T = get_cell_vertices(topo)
  X,T = _delete_vertices(X,T,vertex_to_equal)
  _delete_empty_cells!(T)
  T = Table(T)
  X,T
end

function _map_equal_vertices(stl::DiscreteModel;atol)
  _map_equal_vertices_from_cloud(stl;atol)
end

function _map_equal_vertices_from_cloud(stlmodel::DiscreteModel;atol)
  stl = get_grid(stlmodel)
  group_to_vertices =  _group_vertices(stl;atol=atol*10)
  vertices_map = collect(1:num_nodes(stl))
  for _vertices in group_to_vertices
    for i in 1:length(_vertices), j in i+1:length(_vertices)
      if _vertices[i] != _vertices[j]
        vi = get_node_coordinates(stl)[_vertices[i]]
        vj = get_node_coordinates(stl)[_vertices[j]]
        if vertices_map[_vertices[i]] == _vertices[i] &&
           vertices_map[_vertices[j]] == _vertices[j]

          if distance(vi,vj) < atol
            vertices_map[_vertices[j]] = vertices_map[_vertices[i]]
          end
        end
      end
    end
  end
  for i in 1:length(vertices_map)
    nmax = 10
    n = 0
    v = i
    while vertices_map[v] != v
      v = vertices_map[v]
      n += 1
      @assert n < nmax
    end
    if vertices_map[i] != v
      vertices_map[i] = v
    end
  end
  vertices_map
end

function _delete_vertices(X,T,vertex_to_equal_vertex)
  u = unique(vertex_to_equal_vertex)
  m = vertex_to_equal_vertex .== 1:length(vertex_to_equal_vertex)
  m = cumsum(m) .* m
  vertices_map = m[vertex_to_equal_vertex]
  @assert count(iszero,vertices_map) == 0
  X = X[u]
  T = map( i -> vertices_map[i], T )
  X,T
end

function _delete_empty_cells!(T)
  filter!( f -> length(f) == _num_uniques(f) , T )
end

function _num_uniques(a::AbstractVector)
  c = 0
  for i in 1:length(a)
    if a[i] ∉ view(a,1:i-1)
      c += 1
    end
  end
  c
end

function _group_vertices(stl::Grid{Dc,D};atol) where {Dc,D}
  min_length = _compute_min_length(stl;atol)
  digits = -Int(floor(log10(min_length)))
  pmin, pmax = get_bounding_box(stl)
  cells = NTuple{D,Int}[]
  vertices = Int[]
  m = CartesianIndices( tfill(2,Val{D}()) )
  for (iv,v) in enumerate(get_node_coordinates(stl))
    p = v-pmin
    for i in m
      f = d -> i.I[d] == 1 ? floor(p[d];digits) : ceil(p[d];digits)
      c = ntuple( f, Val{D}() )
      c = Int.(round.(exp10(digits).*c)).+1
      push!(cells,c)
      push!(vertices,iv)
    end
  end
  u = unique(cells)
  dict = Dict( zip(u,1:length(u)) )
  cells = map( i -> dict[i], cells )
  cell_to_vertices = [ Int[] for _ in 1:maximum(cells) ]
  for (cell,vertex) in zip(cells,vertices)
    push!(cell_to_vertices[cell],vertex)
  end
  cell_to_vertices
end

function _compute_min_length(mesh::Grid;atol)
  min_length = Inf
  for nodes in get_cell_node_ids(mesh)
    for i in 1:length(nodes), j in i+1:length(nodes)
      vi = get_node_coordinates(mesh)[nodes[i]]
      vj = get_node_coordinates(mesh)[nodes[j]]
      dist = distance(vi,vj)
      if dist < min_length && dist > atol
        min_length = dist
      end
    end
  end
  min_length
end

function preprocess_small_facets(stl::DiscreteModel;atol)
  @notimplementedif !is_edge_manifold(stl) "The geometry is not edge manifold"
  @notimplementedif !is_water_tight(stl) "The geometry is not watter tight"
  max_iters = 100
  for i in 1:max_iters
    stl,incomplete = _preprocess_small_facets(stl;atol)
    incomplete || return stl
  end
  @warn "Unable to fix small facets in $max_iters iterations"
  stl
end

function _preprocess_small_facets(stl::DiscreteModel{Dc};atol) where Dc
  @notimplementedif Dc ≠ 2
  hang_v,cut_e = get_hanging_vertices_and_edges(stl;atol)
  _stl = stl
  incomplete = false
  if !isempty(hang_v)
    incomplete = true
    touched = falses(num_cells(stl))
    topo = get_grid_topology(stl)
    f_to_v = get_faces(topo,Dc,0)
    e_to_v = get_faces(topo,Dc-1,0)
    e_to_f = get_faces(topo,Dc-1,Dc)
    fv = array_cache(f_to_v)
    ev = array_cache(e_to_v)
    ef = array_cache(e_to_f)
    c = get_facet_cache(stl)
    X = get_vertex_coordinates(topo)
    T = Vector( f_to_v )
    Tnew = Vector{Int32}[]
    for (v,e) in zip(hang_v,cut_e)
      edge = get_facet!(c,stl,e)
      p = projection(X[v],edge)
      X[v] = p
      faces_around = getindex!(ef,e_to_f,e)
      if any(f->touched[f],faces_around)
        continue
      end
      for f in faces_around
        touched[f] = true
        ks = split_face!(fv,ev,f_to_v,e_to_v,f,e,v)
        push!(Tnew,ks...)
      end
    end
    deleteat!(T,touched)
    append!(T,Tnew)
    _delete_empty_cells!(T)
    _stl = compute_stl_model(Table(T),X)
    _stl = merge_nodes(_stl;atol)
    _stl = delete_duplicated_faces(_stl)
    if !is_water_tight(_stl)
      _stl = stl
      incomplete = false
    end
  end
  _stl,incomplete
end

function delete_duplicated_faces(stl::DiscreteModel)
  topo = get_grid_topology(stl)
  D = num_dims(stl)
  v_to_f = get_faces(topo,0,D)
  f_to_v = get_faces(topo,D,0)
  vf = array_cache(v_to_f)
  fvi = array_cache(f_to_v)
  fvj = array_cache(f_to_v)
  is_duplicated = falses(num_cells(stl))
  for v in 1:num_vertices(topo)
    facets = getindex!(vf,v_to_f,v)
    for i in 1:length(facets), j in i+1:length(facets)
      fi = facets[i]
      fj = facets[j]
      fi ≠ fj || continue
      vi = getindex!(fvi,f_to_v,fi)
      vj = getindex!(fvj,f_to_v,fj)
      if vi ⊆ vj
        is_duplicated[fi] = is_duplicated[fj] = true
      end
    end
  end
  facets = findall(d->!d,is_duplicated)
  T = f_to_v[facets]
  X = get_vertex_coordinates(topo)
  compute_stl_model(T,X)
end

function split_face!(fv,ev,f_to_v,e_to_v,f,e,v)
  facet = getindex!(fv,f_to_v,f)
  edge = getindex!(ev,e_to_v,e)
  replace( facet, edge[1] => v ),
  replace( facet, edge[2] => v )
end

function get_hanging_vertices_and_edges(stl::DiscreteModel{Dc};atol) where Dc
  @notimplementedif Dc ≠ 2
  topo = get_grid_topology(stl)
  vertex_coordinates = get_vertex_coordinates(topo)
  f_to_e = get_faces(topo,Dc,Dc-1)
  f_to_v = get_faces(topo,Dc,0)
  e_to_v = get_faces(topo,Dc-1,0)
  ce = array_cache(f_to_e)
  cv = array_cache(f_to_v)
  _cv = array_cache(e_to_v)
  c = get_facet_cache(stl)
  hanging_vertices = Int32[]
  cut_edges = Int32[]
  for f in 1:num_cells(stl)
    min_dist = Inf
    hang_v = UNSET
    cut_e = UNSET
    for e in getindex!(ce,f_to_e,f)
      for v in getindex!(cv,f_to_v,f)
        _vertices = getindex!(_cv,e_to_v,e)
        v ∉ _vertices || continue
        edge = get_facet!(c,stl,e)
        point = vertex_coordinates[v]
        dist = distance(edge,point)
        if dist < min_dist
          min_dist = dist
          hang_v = v
          cut_e = e
        end
      end
    end
    if min_dist < atol
      push!(hanging_vertices,hang_v)
      push!(cut_edges,cut_e)
    end
  end
  hanging_vertices,cut_edges
end

function split_disconnected_parts(stl::DiscreteModel)
  Dc = num_dims(stl)
  f_to_v = get_faces(get_grid_topology(stl),Dc,0)
  v_to_f = get_faces(get_grid_topology(stl),0,Dc)
  v_to_part = fill(UNSET,num_vertices(stl))
  f_to_part = fill(UNSET,num_cells(stl))
  num_parts = 0
  stack = Int[]
  for f in 1:num_cells(stl)
    f_to_part[f] == UNSET || continue
    num_parts += 1
    f_to_part[f] = num_parts
    empty!(stack)
    push!(stack,f)
    while !isempty(stack)
      f_curr = pop!(stack)
      for v in f_to_v[f_curr]
        v_to_part[v] = num_parts
        for f_neig in v_to_f[v]
          f_curr ≠ f_neig || continue
          f_to_part[f_neig] == UNSET || continue
          f_to_part[f_neig] = num_parts
          push!(stack,f_neig)
        end
      end
    end
  end
  @assert all(!isequal(UNSET),v_to_part)
  @assert all(!isequal(UNSET),f_to_part)
  coords = get_vertex_coordinates(get_grid_topology(stl))
  stls = typeof(stl)[]
  for part in 1:num_parts
    faces = findall(isequal(part),f_to_part)
    vertices = findall(isequal(part),v_to_part)
    v_to_part_v = fill(UNSET,length(v_to_part))
    v_to_part_v[vertices] = 1:length(vertices)
    _f_to_v = Table(map(i->v_to_part_v[i],f_to_v[faces]))
    _coords = coords[vertices]
    _stl = compute_stl_model(_f_to_v,_coords)
    push!(stls,_stl)
  end
  stls
end

function check_requisites(stl::DiscreteModel,bgmodel::DiscreteModel;
  verbose=false,max_num_facets=10000)

  !verbose || println(join(fill('-',40)))
  check_requisites(stl;verbose)
  check_facet_density(stl,bgmodel;verbose,max_num_facets)
  true
end

function check_requisites(stl::DiscreteModel;verbose=false)
  if !is_surface(stl)
    error("Is not a surface")
  end
  if !is_vertex_manifold(stl)
    error("Is not vertex manifold")
  end
  if !is_edge_manifold(stl)
    error("Is not edge manifold")
  end
  if !is_water_tight(stl)
    error("Is not water tight")
  end
  if has_sharp_edges(stl)
    error("Has sharp edges")
  end
  true
end

function check_facet_density(stl::DiscreteModel,bgmodel::DiscreteModel;
  verbose=false,max_num_facets)

  max_fv = max_num_facets_per_vertex(stl)
  max_fbc = max_num_facets_per_bgcell(stl,bgmodel)
  !verbose || println("Maximum num facets per vertex: $max_fv")
  !verbose || println("Maximum num facets per bgcell: $max_fbc")
  if max_fbc > max_num_facets
    fulfill = false
    if max_fv > max_num_facets
      error("Unable to run geometry $max_num_facets")
    else
      error("Please refine your mesh")
    end
  end
  true
end

function is_vertex_manifold(stlmodel::DiscreteModel{2,3})
  stl = get_grid_topology(stlmodel)
  D = num_dims(stl)
  v_to_e = get_faces(stl,0,1)
  v_to_f = get_faces(stl,0,D)
  f_to_v = get_faces(stl,D,0)
  fc = array_cache(v_to_f)
  vc = array_cache(f_to_v)
  for v in 1:num_vertices(stl)
    vfacets = getindex!(fc,v_to_f,v)
    !isempty(vfacets) || continue
    f0 = vfacets[1]
    fnext = f0
    nf = 0
    while true
      nf += 1
      i = findfirst(isequal(v), getindex!(vc,f_to_v,fnext) )
      inext = i == D+1 ? 1 : i+1
      vnext = getindex!(vc,f_to_v,fnext)[inext]
      faces = getindex!(fc,v_to_f,vnext)
      i = findfirst( f -> f ≠ fnext && v ∈ getindex!(vc,f_to_v,f), faces )
      fnext = isnothing(i) ? UNSET : faces[i]
      fnext ≠ UNSET || break
      fnext ≠ f0 || break
    end
    if nf ≠ length(vfacets)
      return false
    end
  end
  true
end

function is_edge_manifold(stlmodel::DiscreteModel{2,3})
  stl = get_grid_topology(stlmodel)
  e_to_f = get_faces(stl,1,2)
  c = array_cache(e_to_f)
  maximum(e->length(getindex!(c,e_to_f,e)),1:length(e_to_f)) ≤ 2
end


function max_num_facets_per_bgcell(stlmodel,bgmodel)
  stl = STL(stlmodel)
  grid = get_grid(bgmodel)
  c_to_stlf = compute_cell_to_facets(grid,stl)
  maximum(length,c_to_stlf)
end

function max_num_facets_per_vertex(stlmodel)
  stl = get_grid_topology(stlmodel)
  v_to_f = get_faces(stl,0,2)
  c = array_cache(v_to_f)
  maximum(v->length(getindex!(c,v_to_f,v)),1:length(v_to_f))
end

function has_sharp_edges(stlmodel)
  stl = STL(stlmodel)
  Π = get_facet_planes(stl)
  e_to_f = get_faces(stl,1,2)
  c = array_cache(e_to_f)
  for e in 1:length(e_to_f)
    facets = getindex!(c,e_to_f,e)
    Π1 = Π[facets[1]]
    Π2 = Π[facets[2]]
    n1 = normal(Π1)
    n2 = normal(Π2)
    n1 ⋅ n2 ≉ -1 || return true
  end
  false
end

function get_bounding_box(mesh::Grid)
  pmin = get_node_coordinates(mesh)[1]
  pmax = get_node_coordinates(mesh)[1]
  for vertex in get_node_coordinates(mesh)
    pmin = Point(min.(Tuple(pmin),Tuple(vertex)))
    pmax = Point(max.(Tuple(pmax),Tuple(vertex)))
  end
  pmin,pmax
end

function get_bounding_box(model::DiscreteModel)
  get_bounding_box(get_grid_topology(model))
end

function get_bounding_box(msh::GridTopology)
  vertices = get_vertex_coordinates(msh)
  get_bounding_box(vertices)
end

function get_bounding_box(vertices::AbstractVector{<:VectorValue})
  pmin = pmax = vertices[1]
  for vertex in vertices
    pmin = Point( min.( Tuple(pmin), Tuple(vertex) ) )
    pmax = Point( max.( Tuple(pmax), Tuple(vertex) ) )
  end
  pmin,pmax
end

is_water_tight(model::DiscreteModel) = is_water_tight(get_grid_topology(model))

function is_water_tight(top::GridTopology{Dc,Dp}) where {Dc,Dp}
  Dc == Dp - 1 || error("The Geometry is not a surface")
  c_to_f = get_faces(top,Dc-1,Dc)
  maximum(length,c_to_f) == minimum(length,c_to_f) == 2
end

is_surface(model::DiscreteModel) = is_surface(get_grid_topology(model))

function is_surface(::GridTopology{Dc,Dp}) where {Dc,Dp}
  Dc == Dp-1
end

function measure(a::Grid,mask)
  m = BigFloat(0.0,precision=128)
  p = get_polytope(only(get_reffes(a)))
  T = get_cell_node_ids(a)
  X = get_node_coordinates(a)
  c = array_cache(T)
  for i in 1:num_cells(a)
    if mask[i]
      f = get_cell!(c,T,X,p,i)
      m += measure(f)
    end
  end
  Float64(m)
end

function measure(a::Grid,cell_to_inoutcut,inoutcut::Symbol)
  inout_dict = Dict{Symbol,Int8}(
    :in => FACE_IN,
    :out => FACE_OUT,
    :cut => FACE_CUT )
  measure(a,cell_to_inoutcut,inout_dict[inoutcut])
end

function measure(a::Grid,cell_to_val,val)
  measure(a,lazy_map(isequal(val),cell_to_val))
end

function measure(a::Grid)
  measure(a,lazy_map(i->true,1:num_cells(a)))
end

volume(a::Grid{D,D},args...) where D = measure(a,args...)

function surface(a::Grid{Df,Dp},args...) where {Df,Dp}
  @notimplementedif Df ≠ Dp-1
  measure(a,args...)
end

surface(a::DiscreteModel) = surface(get_grid(a))

function measures(a::Grid,cell_to_inoutcut)
  measure(a,cell_to_inoutcut,FACE_IN),
  measure(a,cell_to_inoutcut,FACE_OUT),
  measure(a,cell_to_inoutcut,FACE_CUT)
end

volumes(a::Grid{D,D},args...) where D = measures(a,args...)

function surfaces(a::Grid{Df,Dp},args...) where {Df,Dp}
  @notimplementedif Df ≠ Dp-1
  measures(a,args...)
end

function min_height(grid::Grid)
  min_h = Inf
  c = get_cell_cache(grid)
  for i in 1:num_cells(grid)
    cell = get_cell!(c,grid,i)
    h = min_height(cell)
    if h < min_h
      min_h = h
    end
  end
  min_h
end

min_height(model::DiscreteModel) = min_height(get_grid(model))

function max_length(grid::Grid)
  max_len = 0.0
  c = get_cell_cache(grid)
  for i in 1:num_cells(grid)
    cell = get_cell!(c,grid,i)
    max_len = max(max_len,max_length(cell))
  end
  max_len
end

max_length(model::DiscreteModel) = max_length(get_grid(model))

# Cartesian Grid

function get_bounding_box(grid::CartesianGrid)
  desc = get_cartesian_descriptor(grid)
  pmin = desc.origin
  pmax = desc.origin + VectorValue(desc.partition .* desc.sizes)
  pmin,pmax
end

function measure(a::CartesianGrid)
  pmin,pmax = get_bounding_box(a)
  measure(pmin,pmax)
end

function measure(a::CartesianGrid{D,T,typeof(identity)},mask) where {D,T}
  cell_vol = prod( get_cartesian_descriptor(a).sizes )
  cell_vol * count(mask)
end

function compute_cell_to_facets(model_a::DiscreteModel,grid_b)
  grid_a = get_grid(model_a)
  compute_cell_to_facets(grid_a,grid_b)
end


function compute_cell_to_facets(grid::CartesianGrid,stl)
  CELL_EXPANSION_FACTOR = 1e-3
  desc = get_cartesian_descriptor(grid)
  @assert length(get_reffes(grid)) == 1
  p = get_polytope(get_cell_reffe(grid)[1])
  @notimplementedif desc.map !== identity
  cell_to_stl_facets = [ Int32[] for _ in 1:num_cells(grid) ]
  n = Threads.nthreads()
  thread_to_cells = [ Int32[] for _ in 1:n ]
  thread_to_stl_facets = [ Int32[] for _ in 1:n ]
  coords = get_node_coordinates(grid)
  cell_to_nodes = get_cell_node_ids(grid)
  c = [ ( get_cell_cache(stl), array_cache(cell_to_nodes) ) for _ in 1:n ]
  δ = CELL_EXPANSION_FACTOR
  Threads.@threads for stl_facet in 1:num_cells(stl)
    i = Threads.threadid()
    cc,nc = c[i]
    f = get_cell!(cc,stl,stl_facet)
    pmin,pmax = get_bounding_box(f)
    for cid in get_cells_around(desc,pmin,pmax)
      cell = LinearIndices(desc.partition)[cid.I...]
      nodes = getindex!(nc,cell_to_nodes,cell)
      _pmin = coords[nodes[1]]
      _pmax = coords[nodes[end]]
      Δ = (_pmax - _pmin) * δ
      _pmin = _pmin - Δ
      _pmax = _pmax + Δ
      if voxel_intersection(f,_pmin,_pmax,p)
        push!(thread_to_cells[i],cell)
        push!(thread_to_stl_facets[i],stl_facet)
      end
    end
  end
  cell_to_stl_facets = [ Int32[] for _ in 1:num_cells(grid) ]
  for (cells,stl_facets) in zip(thread_to_cells,thread_to_stl_facets)
    for (cell,stl_facet) in zip(cells,stl_facets)
      push!(cell_to_stl_facets[cell],stl_facet)
    end
  end
  cell_to_stl_facets
end

function compute_cell_to_facets(a::UnstructuredGrid,b)
  tmp = cartesian_bounding_model(a)
  tmp_to_a = compute_cell_to_facets(tmp,a)
  tmp_to_b = compute_cell_to_facets(tmp,b)
  a_to_tmp = inverse_index_map(tmp_to_a,num_cells(a))
  a_to_b = compose_index_map(a_to_tmp,tmp_to_b)
  a_to_b # TODO: filter by true intersections (should work without filter
end

function cartesian_bounding_model(grid::Grid)
  h = max_length(grid)
  pmin,pmax = get_bounding_box(grid)
  s = Tuple(pmax-pmin)
  partition = Int.( ceil.( s ./ h  ))
  CartesianDiscreteModel(pmin,pmax,partition)
end

function inverse_index_map(
  a_to_b::AbstractVector{<:AbstractVector},
  nb=maximum(maximum(a_to_b)))

  na = length(a_to_b)
  as = Int32[]
  bs = Int32[]
  cache = array_cache(a_to_b)
  for a in 1:na
    for b in getindex!(cache,a_to_b,a)
      push!(as,a)
      push!(bs,b)
    end
  end
  assemble_sparse_map(bs,as,nb)
end

function compose_index_map(
  a_to_b::AbstractVector{<:AbstractVector},
  b_to_c::AbstractVector{<:AbstractVector})

  na = length(a_to_b)
  as = Int32[]
  cs = Int32[]
  cache_ab = array_cache(a_to_b)
  cache_bc = array_cache(b_to_c)
  for a in 1:na
    for b in getindex!(cache_ab,a_to_b,a)
      for c in getindex!(cache_bc,b_to_c,b)
        push!(as,a)
        push!(cs,c)
      end
    end
  end
  assemble_sparse_map(as,cs,na)
end

function assemble_sparse_map(as,bs,na=maximum(as))
  a_to_b = map(_->Int32[],1:na)
  for (a,b) in zip(as,bs)
    if b ∉ a_to_b[a]
      push!(a_to_b[a],b)
    end
  end
  a_to_b
end

function get_cells_around(desc::CartesianDescriptor{D},pmin::Point,pmax::Point) where D
  cmin,_ = get_cell_bounds(desc,pmin)
  _,cmax = get_cell_bounds(desc,pmax)
  cmin = CartesianIndices(desc.partition)[cmin]
  cmax = CartesianIndices(desc.partition)[cmax]
  ranges = ntuple( i -> cmin.I[i] : cmax.I[i], Val{D}() )
  CartesianIndices( ranges )
end

function get_cell_bounds(desc::CartesianDescriptor,p::Point)
  function _get_cell(cell)
    cell = Int.(cell)
    cell = max.(cell,1)
    cell = min.(cell,desc.partition)
    LinearIndices(desc.partition)[cell...]
  end
  tol = 0.1
  coords = Tuple(p-desc.origin)./desc.sizes
  cell = floor.(coords).+1
  cell_min = cell .- ( (coords.-(floor.(coords))) .< tol )
  cell_max = cell .+ ( (coords.-(floor.(coords))) .> (1-tol) )
  _get_cell(cell_min),_get_cell(cell_max)
end

# Bisectors

function bisector_plane!(
   cache,
   stl::STL{Dc,Dp},
   d::Integer,
   dface::Integer,
   Πf::AbstractArray) where {Dc,Dp}

  @notimplementedif Dc ≠ Dp-1
  @notimplementedif d ≠ Dc-1
  c,fc = cache
  e_to_f = get_faces(stl,d,Dc)
  facets = getindex!(c,e_to_f,dface)
  length(facets) == 2 || return Πf[ only(facets) ]
  edge = get_dface!(fc,stl,dface,Val{Dc-1}())
  Π1 = Πf[ facets[1] ]
  Π2 = Πf[ facets[2] ]
  bisector_plane(edge,Π1,Π2)
end

function bisector_plane(
   stl::STL,
   d::Integer,
   dface::Integer,
   Πf::AbstractArray)

  c = bisector_plane_cache(stl,d)
  bisector_plane(c,stl,d,dface,Πf)
end

function bisector_plane_cache(stl::STL,d::Integer)
  Dc = num_dims(stl)
  e_to_f = get_faces(stl,d,Dc)
  c = array_cache(e_to_f)
  fc = get_dface_cache(stl,d)
  c,fc
end
