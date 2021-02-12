
function read_stl(filename::String)
  if !unknown(query(filename))
    mesh = load(filename)
  else
    F = _file_format(filename)
    mesh = load(File(F,filename))
  end
  vertex_to_coordinates = MeshIO.decompose(MeshIO.Point{3,Float64}, mesh )
  facet_to_vertices = MeshIO.decompose(MeshIO.Face{3,Int},mesh)
  vertex_to_normals = MeshIO.decompose(MeshIO.Normal{3,Float64},mesh)
  vertex_to_coordinates = Vector{Point{3,Float64}}( vertex_to_coordinates )
  facet_to_vertices = Table( Vector{Vector{Int}}( facet_to_vertices ) )
  vertex_to_normals = Vector{VectorValue{3,Float64}}( vertex_to_normals )
  facet_to_normals = vertex_to_normals[1:3:end]
  vertex_to_coordinates,facet_to_vertices,facet_to_normals
end

function _file_format(filename)
  if endswith(filename,r".stl|.STL")
    cmd = `file --mime-encoding $filename`
    out = read(cmd,String)
    if contains(out,"binary")
      F = format"STL_BINARY"
    elseif contains(out,"ascii")
      F = format"STL_ASCII"
    else
      msg = "Invalid file encoding: $out"
      @unreachable msg 
    end
  elseif  endswith(filename,r".obj")
    F = format"OBJ"
  else
    msg = "Invalid file extension: $filename"
    @unreachable msg
  end
end

function Base.convert(::Type{Point{D,T}},x::MeshIO.Point{D}) where {D,T} 
  Point{D,T}(x.data)
end

function Base.convert(::Type{VectorValue{D,T}},x::MeshIO.Normal{D})  where {D,T}
  VectorValue{D,T}(x.data)
end

Base.convert(::Type{Vector},x::MeshIO.Face) = collect(x.data)

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

function get_dface(stl::DiscreteModel,iface::Integer,::Val{d}) where d
  @notimplementedif length( get_polytopes(stl) ) ≠ 1 
  @notimplementedif !is_simplex( get_polytopes(stl)[1] )
  0 < iface ≤ num_faces(stl,d) || 
    throw(AssertionError(
      "get_dface() :: access at $d-face $iface of $(num_faces(stl,d))"))
  stl_topology = get_grid_topology(stl)
  face_vertices = get_faces(stl_topology,d,0)[iface]
  X = get_vertex_coordinates(stl_topology)
  v = ntuple( i -> X[face_vertices[i]], Val{d+1}() )
  simplex_face(v)
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

function merge_nodes(stl::DiscreteModel)
  X,T = delete_repeated_vertices(get_grid(stl))
  compute_stl_model(T,X)
end

function delete_repeated_vertices(stl::Grid)
  group_to_vertices =  _group_vertices(stl)
  vertices_map = collect(1:num_nodes(stl))
  for _vertices in group_to_vertices
    for i in 1:length(_vertices), j in i+1:length(_vertices)
      vi = get_node_coordinates(stl)[_vertices[i]]
      vj = get_node_coordinates(stl)[_vertices[j]]
      if vertices_map[_vertices[i]] == _vertices[i] && 
         vertices_map[_vertices[j]] == _vertices[j]

        if distance(vi,vj) < TOL
          vertices_map[_vertices[j]] = _vertices[i]
        end
      end
    end
  end
  u = unique(vertices_map)
  m = vertices_map .== 1:length(vertices_map)
  m = cumsum(m) .* m
  vertices_map = m[vertices_map]
  X = get_node_coordinates(stl)[u]
  T = get_cell_nodes(stl)
  T = Table( map( i -> vertices_map[i], T ) )
  X,T
end

function _group_vertices(stl::Grid{Dc,D}) where {Dc,D}
  min_length = _compute_min_length(stl)
  num_digits = -Int(floor(log10(min_length)))
  pmin, pmax = get_bouding_box(stl)
  cells = Int[]
  vertices = Int[]
  ranks = ceil.(Tuple(pmax-pmin),digits=num_digits)
  ranks = Int.(round.(exp10(num_digits).*ranks)).+1
  lids = LinearIndices( ranks )
  m = CartesianIndices( tfill(2,Val{D}()) )
  for (iv,v) in enumerate(get_node_coordinates(stl))
    p = v-pmin
    for i in m
      f = d -> i.I[d] == 1 ? floor(p[d],digits=num_digits) : ceil(p[d],digits=num_digits)
      r = ntuple( f, Val{D}() )
      r = Int.(round.(exp10(num_digits).*r)).+1
      cell = lids[r...]
      push!(cells,cell)
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

function get_bouding_box(mesh::Grid)
  pmin = get_node_coordinates(mesh)[1]
  pmax = get_node_coordinates(mesh)[1]
  for vertex in get_node_coordinates(mesh)
    pmin = Point(min.(Tuple(pmin),Tuple(vertex)))
    pmax = Point(max.(Tuple(pmax),Tuple(vertex)))
  end
  pmin,pmax
end

function _compute_min_length(mesh::Grid)
  min_length = Inf
  for nodes in get_cell_nodes(mesh)
    for i in 1:length(nodes), j in i+1:length(nodes)
      vi = get_node_coordinates(mesh)[nodes[i]]
      vj = get_node_coordinates(mesh)[nodes[j]]
      dist = distance(vi,vj)
      if dist < min_length
        min_length = dist
      end
    end
  end
  min_length
end

is_water_tight(model::DiscreteModel) = is_water_tight(get_grid_topology(model))

is_open_surface(model::DiscreteModel) = is_open_surface(get_grid_topology(model))

function is_water_tight(top::GridTopology{Dc}) where Dc
  l = length.( get_faces(top,Dc-1,Dc) )
  maximum(l) == minimum(l) == 2
end

function is_surface(::GridTopology{Dc,Dp}) where {Dc,Dp}
  Dc == Dp-1
end

function is_open_surface(top::GridTopology{Dc}) where Dc
  is_surface(top) || return false
  l = length.( get_faces(top,Dc-1,Dc) )
  (1 ≤ maximum(l) ≤ 2) && (1 ≤ minimum(l) ≤ 2)
end

function have_intersection(
  grid::Grid,
  cell::Integer,
  m::DiscreteModel,
  face::Integer
  ;atol=nothing)

  dispatch_face(have_intersection,grid,cell,m,face,atol=atol)
end

function is_on_boundary(
  grid::Grid,
  cell::Integer,
  m::DiscreteModel,
  face::Integer
  ;atol::Real)

  dispatch_face(is_on_boundary,grid,cell,m,face,atol=atol)
end

function get_bounding_box(m::DiscreteModel,face::Integer)
  dispatch_face(get_bounding_box,m,face)
end


function dispatch_face(fun,grid::Grid,cell::Integer,args...;kargs...)
  c = get_cell(grid,cell)
  dispatch_face(fun,c,args...;kargs...)
end

function dispatch_face(fun,a,b::DiscreteModel,face::Integer;kargs...)
  topo = get_grid_topology(b)
  d = get_facedims(topo)[face]
  dface = face - get_dimrange(topo,d)[1] + 1
  dispatch_face(fun,a,b,d,dface;kargs...)
end

function dispatch_face(fun,a::DiscreteModel,face::Integer;kargs...)
  topo = get_grid_topology(a)
  d = get_facedims(topo)[face]
  dface = face - get_dimrange(topo,d)[1] + 1
  dispatch_face(fun,a,d,dface;kargs...)
end

@generated(
function dispatch_face(
  fun,a,b::DiscreteModel{D},d::Integer,df::Integer;kargs...) where D

  str = ""
  for d in 0:D
    str *= "if d == $d \n"
    str *= "  f$d =  get_dface(b,df,Val{$d}()) \n"
    str *= "  fun(a,f$d;kargs...) \n"
    str *= "else"
  end
  str *= "\n  @notimplemented \nend"
  Meta.parse(str)
end
)

@generated(
function dispatch_face(
  fun,a::DiscreteModel{D},d::Integer,df::Integer;kargs...) where D

  str = ""
  for d in 0:D
    str *= "if d == $d \n"
    str *= "  f$d =  get_dface(a,df,Val{$d}()) \n"
    str *= "  fun(f$d;kargs...) \n"
    str *= "else"
  end
  str *= "\n  @notimplemented \nend"
  Meta.parse(str)
end
)

function measure(a::Grid)
  m = 0.0
  for i in 1:num_cells(a)
    c = get_cell(a,i)
    m += measure(c)
  end
  m
end

function measures(a::Grid,num::Integer,map)
  m = zeros(num)
  for i in 1:num_cells(a)
    c = get_cell(a,i)
    m[map[i]] += measure(c)
  end
  m
end

function measures(a::Grid)
  measures(a,num_cells(a),1:num_cells(a))
end

volume(a::Grid{D,D}) where D = measure(a)

volumes(a::Grid{D,D},args...) where D = measures(a,args...)

function surface(a::Grid{Df,Dp}) where {Df,Dp}
  @notimplementedif Df ≠ Dp-1
  measure(a)
end

function surfaces(a::Grid{Df,Dp},num::Integer,map) where {Df,Dp}
  @notimplementedif Df ≠ Dp-1
  measures(a,num,map)
end

surface(a::DiscreteModel) = surface(get_grid(a))

surfaces(a::DiscreteModel,args...) = surfaces(get_grid(a),args...)

function get_bounding_box(stl::DiscreteModel)
  vertices = get_vertex_coordinates(get_grid_topology(stl))
  pmin = pmax = vertices[1]
  for vertex in vertices
    pmin = Point( min.( Tuple(pmin), Tuple(vertex) ) )
    pmax = Point( max.( Tuple(pmax), Tuple(vertex) ) )
  end
  pmin,pmax
end
