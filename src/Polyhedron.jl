
## Constants

const OPEN = -1

struct Polyhedron{Dp,Tp,Td}
  vertices::Vector{Point{Dp,Tp}}
  edge_vertex_graph::Vector{Vector{Int}}
  isopen::Bool
  data::Td
end

struct PolyhedronData
  vertex_to_planes::Vector{Vector{Int}}
  vertex_to_original_faces::Vector{Vector{Int}}
  vertex_to_parent_vertex::Vector{Int}
  plane_to_vertex_to_distances::Vector{Vector{Float64}}
  plane_to_ids::Vector{Int}
end

struct CartesianPlane{D,T}
  d::Int8
  value::T
  positive::Bool
end

## Constructors

function Polyhedron(stl::DiscreteModel{Dc,Dp}) where {Dc,Dp}
  𝓖 = compute_graph(stl)
  X = get_node_coordinates(stl)
  isopen = is_open_surface(stl)
  p = Polyhedron(X,𝓖;isopen)
  set_original_faces!(p,stl)
  p
end

function Polyhedron(
  vertices::AbstractVector{<:Point},
  graph::Vector{<:Vector}
  ;isopen=false::Bool)

  Polyhedron(vertices,graph,isopen,polyhedron_data(length(vertices)))
end

function Polyhedron(
  vertices::AbstractVector{<:Point},
  graph::Vector{<:Vector},
  data)

  isopen = false
  Polyhedron(vertices,graph,isopen,data)
end


function Polyhedron(p::Polytope{2},vertices::AbstractVector{<:Point})
  if p == TRI
    e_v_graph = [[2,3],[3,1],[1,2]]
  elseif p == QUAD
    e_v_graph = [[2, 3],[4, 1],[1, 4],[3, 2]]
  else
    @unreachable
  end
  data = polyhedron_data(p)
  Polyhedron(vertices,e_v_graph,data)
end

function Polyhedron(p::Polytope{3},vertices::AbstractVector{<:Point})
  if p == TET
    e_v_graph = [[2,4,3],[3,4,1],[1,4,2]]
  elseif p == HEX
    e_v_graph = [
      [5, 2, 3],
      [6, 4, 1],
      [7, 1, 4],
      [8, 3, 2],
      [1, 7, 6],
      [2, 5, 8],
      [3, 8, 5],
      [4, 6, 7] ]
  else
    @unreachable
  end
  data = polyhedron_data(p)
  Polyhedron(vertices,e_v_graph,data)
end

function Polyhedron(p::Polytope)
  Polyhedron(p,get_vertex_coordinates(p))
end

function CartesianPlane(p::Point{D,T},d::Integer,orientation::Integer) where {D,T}
  @notimplementedif abs( orientation ) ≠ 1
  CartesianPlane{D,T}(Int8(d),p[d],orientation>0)
end

function compute_distances!(p::Polyhedron,Π,faces)
  data = get_data(p)
  v_to_f = get_data(p).vertex_to_original_faces
  Π_to_v_to_d = get_plane_distances(data)
  Π_ids = get_plane_ids(data)
  num_Π = length(Π_ids)
  append!(Π_ids,faces)
  append!(Π_to_v_to_d, [ zeros(num_vertices(p)) for _ in faces ] )
  for (i,(Πi,fi)) in enumerate(zip(Π,faces))
    v_to_d = Π_to_v_to_d[i+num_Π]
    for v in 1:num_vertices(p)
      if fi ∈ v_to_f[v]
        dist = 0.0
      else
        dist = signed_distance(p[v],Πi)
      end
      v_to_d[v] = dist
    end
  end
end

function get_plane_distances(data::PolyhedronData)
  data.plane_to_vertex_to_distances
end

function get_plane_distances(data::PolyhedronData,Π)
  i = findfirst(isequal(Π),get_plane_ids(data))
  get_plane_distances(data)[i]
end

function get_plane_ids(data::PolyhedronData)
  data.plane_to_ids
end

## Polyhedra operations

function clip(poly::Polyhedron,Π;inside=true)
  p = poly
  for Πi in Π
    p⁻,p⁺ = split(p,Πi)
    p = inside ? p⁻ : p⁺
    !isnothing(p) || break
  end
  p
end

function split(p::Polyhedron,Π)
  distances = get_plane_distances(get_data(p),Π)
  smin = Inf
  smax = -Inf
  for v in 1:num_vertices(p)
    isactive(p,v) || continue
    smin = min(smin,distances[v])
    smax = max(smax,distances[v])
  end
  if smin ≥ 0
    return nothing,p
  end
  if smax ≤ 0
    return p,nothing
  end
  new_vertices = empty( get_vertex_coordinates(p) )
  in_graph = deepcopy( get_graph(p) )
  out_graph = deepcopy( get_graph(p) )
  data = deepcopy( get_data(p) )
  D = num_dims(p)
  for v in 1:num_vertices(p)
    isactive(p,v) || continue
    distances[v] > 0 && continue
    for (i,vneig) in enumerate( get_graph(p)[v] )
      vneig ∉ (UNSET,OPEN) || continue
      distances[vneig] > 0 || continue
      vertex = compute_intersection(p[v],distances[v],p[vneig],distances[vneig])
      push!( new_vertices, vertex )
      add_vertex!(data,v,vneig,Π)

      push!( in_graph, fill(UNSET,D) )
      in_graph[v][i] = num_vertices(p) + length(new_vertices)
      in_graph[end][1] = v
      
      push!( out_graph, fill(UNSET,D) )
      ineig = findfirst( isequal(v), get_graph(p)[vneig] )
      out_graph[vneig][ineig] = num_vertices(p) + length(new_vertices)
      out_graph[end][1] = vneig
    end
  end
  complete_graph!(in_graph,num_vertices(p))
  complete_graph!(out_graph,num_vertices(p))
  disconnect_graph!(in_graph,num_vertices(p),distances,true)
  disconnect_graph!(out_graph,num_vertices(p),distances,false)
  add_open_vertices!(in_graph,p)
  add_open_vertices!(out_graph,p)
  vertices = [p.vertices;new_vertices]
  in_data = data
  out_data = deepcopy(data)
  Polyhedron(vertices,in_graph,isopen(p),in_data), 
  Polyhedron(vertices,out_graph,isopen(p),out_data)
end

function flip(p::Polyhedron)
  graph = flip_graph(get_graph(p))
  Polyhedron(get_vertex_coordinates(p),graph,get_data(p))
end

function flip_graph(graph)
  map( i -> reverse(i), graph )
end

function decompose(surf::Polyhedron,cell::Polyhedron,rfaces)
  i = findfirst(i->has_original_reflex_face(surf,i), rfaces )
  if i === nothing
    R = [surf],[cell]
  else
    rf = rfaces[i]
    s⁻,s⁺ = split(surf,rf)
    k⁻,k⁺ = split(cell,rf)
    S = [s⁻,s⁺]
    K = [k⁻,k⁺]
    if length(rfaces) == i
      R = S,K
    else
      Sr,Kr = empty(S),empty(K)
      rfaces = view(rfaces,i+1:length(rfaces))
      for (s,k) in zip(S,K)
        Si,Ki = decompose(s,k,rfaces)
        append!(Sr,Si)
        append!(Kr,Ki)
      end
      R = Sr,Kr
    end
  end
  R
end

function simplexify(poly::Polyhedron{3})
  vstart = fill(UNSET,num_vertices(poly))
  stack = Int[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    vstart[v] == UNSET || continue
    vstart[v] = v
    empty!(stack)
    push!(stack,v)
    while !isempty(stack)
      vcurrent = pop!(stack)
      for vneig in get_graph(poly)[vcurrent]
        if vstart[vneig] == UNSET
          vstart[vneig] = v
          push!(stack,vneig)
        end
      end
    end
  end
  v_to_pv = get_data(poly).vertex_to_parent_vertex
  istouch = map( i -> zeros(Bool,length(i)), get_graph(poly) )
  T = Vector{Int}[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    for i in 1:length(get_graph(poly)[v])
      !istouch[v][i] || continue
      istouch[v][i] = true
      vcurrent = v
      vnext = get_graph(poly)[v][i]
      while vnext != v
        inext = findfirst( isequal(vcurrent), get_graph(poly)[vnext] )
        inext = ( inext % length( get_graph(poly)[vnext] ) ) + 1
        istouch[vnext][inext] = true
        vcurrent = vnext
        vnext = get_graph(poly)[vnext][inext]
        if v ∉ (vstart[v],vcurrent,vnext)
          parents = map(i->v_to_pv[i],(v,vstart[v],vcurrent,vnext))
          zero_vol = false
          for i in 1:length(parents), j in i+1:length(parents)
            if parents[i] == parents[j]
              zero_vol = true
              break
            end
          end
          !zero_vol || continue
          k = [vstart[v],v,vcurrent,vnext]
          push!(T,k)
        end
      end
    end
  end
  T,get_vertex_coordinates(poly)
end

function simplexify(polys::AbstractVector{<:Polyhedron})
  T = Vector{Int}[]
  X = empty(get_vertex_coordinates(first(polys)))
  for poly in polys
    Ti,Xi = simplexify(poly)
    append!(T, map(i->i.+length(X),Ti) )
    append!(X,Xi)
  end
  T,X
end

## Getters

num_dims(::Polyhedron{D}) where D = D

num_vertices(a::Polyhedron) = length(a.vertices)

get_vertex_coordinates(a::Polyhedron) = a.vertices

Base.getindex(a::Polyhedron,i::Integer) = a.vertices[i]

get_graph(a::Polyhedron) = a.edge_vertex_graph

get_data(a::Polyhedron) = a.data

Base.isopen(a::Polyhedron) = a.isopen

get_vertex_to_planes(a::PolyhedronData) = a.vertex_to_planes

function isactive(p::Polyhedron,vertex::Integer)
  !isempty( get_graph(p)[vertex] )
end

function edge_mesh(p::Polyhedron)
  edge_to_vertices = Vector{Int}[]
  for v in 1:num_vertices(p)
    for vneig in get_graph(p)[v]
      if vneig > v
        push!(edge_to_vertices,[v,vneig])
      end
    end
  end
  X = get_vertex_coordinates(p)
  T = Table(edge_to_vertices)
  UnstructuredGrid(X,T,[SEG2],fill(1,length(T)))
end

function get_original_reflex_faces(p::Polyhedron{D},stl::DiscreteModel) where D
  get_original_faces(p,stl,Val{D-2}())
end

function get_original_facets(p::Polyhedron{D},stl::DiscreteModel) where D
  get_original_faces(p,stl,Val{D-1}())
end

function get_original_faces(p::Polyhedron,stl::DiscreteModel,::Val{d}) where d
  faces = collect_original_faces(p,stl,d)
  filter!(f->has_original_face(p,f,Val{d}()),faces)
  faces
end

function collect_original_faces(p::Polyhedron,stl::DiscreteModel,d::Integer)
  v_to_f = get_data(p).vertex_to_original_faces
  facedims = get_facedims(get_grid_topology(stl))
  faces = Int[]
  for v in 1:num_vertices(p)
    if isactive(p,v)
      for f ∈ v_to_f[v]
        if facedims[f] == d && f ∉ faces
          push!(faces,f)
        end
      end
    end
  end
  faces
end

function has_original_reflex_face(p::Polyhedron{D},face::Integer) where D
  has_original_face(p,face,Val{D-2}())
end

function has_original_facet(p::Polyhedron{D},face::Integer) where D
  has_original_face(p,face,Val{D-1}())
end

function has_original_face(p::Polyhedron,face::Integer,::Val{0})
  v_to_f = get_data(p).vertex_to_original_faces
  for v in 1:num_vertices(p)
    if isactive(p,v) && face ∈ v_to_f[v]
      return true
    end
  end
  false
end

function has_original_face(p::Polyhedron,face::Integer,::Val{1})
  v_to_f = get_data(p).vertex_to_original_faces
  v_to_v = get_data(p).vertex_to_parent_vertex
  for v in 1:num_vertices(p)
    if isactive(p,v) && face ∈ v_to_f[v]
      for vneig in get_graph(p)[v]
        vneig ∉ (OPEN,UNSET) || continue
        if face ∈ v_to_f[vneig] && v_to_v[v] ≠ v_to_v[vneig]
          return true
        end
      end
    end
  end
  false
end

function has_original_face(p::Polyhedron,face::Integer,::Val{2})
  d = 2
  v_to_f = get_data(p).vertex_to_original_faces
  v_to_v = get_data(p).vertex_to_parent_vertex
  for v in 1:num_vertices(p)
    if isactive(p,v) && face ∈ v_to_f[v]
      for vneig in get_graph(p)[v]
        vneig ∉ (OPEN,UNSET) || continue
        if face ∈ v_to_f[vneig]
          num_v = 1
          vcurrent = v
          vnext = vneig
          while vnext ≠ v 
            if v_to_v[vnext] ∉ (v_to_v[v],v_to_v[vcurrent])
              num_v += 1
            end
            vcurrent,vnext = vnext,next_vertex(p,vcurrent,vnext)
            vnext ∉ (UNSET,OPEN) || break
            face ∈ v_to_f[vnext] || break
          end
          if vnext == v && num_v ≥ d+1
            return true
          end
        end
      end
    end
  end
  false
end

## Setup


function compute_cell_to_facets(grid::CartesianGrid,stl::DiscreteModel)
  desc = get_cartesian_descriptor(grid)
  @notimplementedif desc.map !== identity
  cell_to_stl_facets = [ Int[] for _ in 1:num_cells(grid) ]
  δ = 0.1
  for stl_facet in 1:num_cells(stl)
    f = get_cell(stl,stl_facet)
    pmin,pmax = get_bounding_box(f)
    for cid in get_cells_around(desc,pmin,pmax)
      cell = LinearIndices(desc.partition)[cid.I...]
      _pmin = get_cell_coordinates(grid)[cell][1]
      _pmax = get_cell_coordinates(grid)[cell][end]
      Δ = (_pmax - _pmin) * δ
      _pmin = _pmin - Δ
      _pmax = _pmax + Δ
      if fast_intersection(f,_pmin,_pmax)
        push!(cell_to_stl_facets[cell],stl_facet)
      end
    end
  end
  cell_to_stl_facets
end

function get_cells_around(desc::CartesianDescriptor{D},pmin::Point,pmax::Point) where D
  cmin,_ = get_cell_bounds(desc,pmin)
  _,cmax = get_cell_bounds(desc,pmax)
  cmin = CartesianIndices(desc.partition)[cmin]
  cmax = CartesianIndices(desc.partition)[cmax]
  ranges = ntuple( i -> cmin.I[i]:cmax.I[i], Val{D}() )
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

## Helpers

function polyhedron_data(num_vertices::Integer)
  v_to_Π = [ Int[] for _ in 1:num_vertices ]
  v_to_of = [ Int[] for _ in 1:num_vertices ]
  v_to_v = collect(1:num_vertices)
  Π_to_v_to_d = Vector{Int}[]
  Π_to_id = Int[]
  PolyhedronData( v_to_Π, v_to_of, v_to_v, Π_to_v_to_d, Π_to_id  )
end

function polyhedron_data(p::Polytope)
  v_to_Π = deepcopy( get_faces(p,0,num_dims(p)-1) )
  v_to_of = [ Int[] for _ in 1:num_vertices(p) ]
  v_to_v = collect(1:num_vertices(p))
  Π_to_v_to_d = Vector{Int}[]
  Π_to_id = Int[]
  PolyhedronData( v_to_Π, v_to_of, v_to_v, Π_to_v_to_d, Π_to_id  )
end

function signed_distance(point::Point{D},Π::CartesianPlane{D}) where D
  Π.positive ? point[Π.d] - Π.value : Π.value - point[Π.d]
end

function complete_graph!(edge_graph,num_vertices::Integer)
  for v in num_vertices+1:length(edge_graph)
    vnext = next_vertex(edge_graph,num_vertices,v)
    vnext ∉ (UNSET,OPEN) || continue
    edge_graph[v][end] = vnext
    edge_graph[vnext][2] = v
  end
end

function add_open_vertices!(graph,poly::Polyhedron)
  if isopen(poly)
    for v in num_vertices(poly)+1:length(graph)
      if UNSET ∉ graph[v]
        push!(graph[v],graph[v][end])
        graph[v][end-1] = OPEN
      end
    end
  end
end

function disconnect_graph!(edge_graph,num_vertices,distances,mask::Bool)
  for i in 1:num_vertices
    if mask == ( distances[i] > 0 )
      edge_graph[i] = empty!( edge_graph[i] )
    end
  end
end

function next_vertex(edge_graph::Vector,num_vertices::Integer,vstart::Integer)
  vcurrent = vstart
  vnext = first( edge_graph[vcurrent] )
  while vnext ≤ num_vertices && vnext ∉ (UNSET,OPEN)
    i = findfirst( isequal(vcurrent), edge_graph[vnext] )
    vcurrent = vnext
    inext = ( i % length( edge_graph[vcurrent] ) ) + 1
    vnext = edge_graph[vcurrent][ inext ]
  end
  vnext
end

function compute_intersection(p1::Point,d1::Real,p2::Point,d2::Real)
  ( p1*(abs(d2))+p2*(abs(d1)) ) / ( abs(d1) + abs(d2) )
end

add_vertex!(data::Nothing,a...) = data 

function add_vertex!(data::PolyhedronData,v1::Integer,v2::Integer,Πid::Integer)
  v_to_Π = data.vertex_to_planes
  v_to_f = data.vertex_to_original_faces
  v_to_pv = data.vertex_to_parent_vertex
  Πs = intersect( v_to_Π[v1], v_to_Π[v2] )
  fs = intersect( v_to_f[v1], v_to_f[v2] )
  push!( Πs, Πid )
  push!( v_to_Π, Πs )
  push!( v_to_f, fs )
  Π_to_v_to_d = get_plane_distances(data)
  Π_to_id = get_plane_ids(data)
  i = findfirst(isequal(Πid),Π_to_id)
  d1 = Π_to_v_to_d[i][v1] # Distances may be in the input
  d2 = Π_to_v_to_d[i][v2]
  pv = iszero(d1) ? v_to_pv[v1] : length(v_to_pv)+1 
  pv = iszero(d2) ? v_to_pv[v2] : pv 
  push!(v_to_pv,pv)
  @assert d1 ≠ d2
  for i in 1:length(Π_to_v_to_d)
    v_to_d = Π_to_v_to_d[i] 
    dist = Π_to_id[i] == Πid ? 0.0 : (d1*v_to_d[v2]-d2*v_to_d[v1])/(d1-d2)
    push!(v_to_d,dist)
  end
  data
end

function compute_graph(stl::DiscreteModel{2})
#  @notimplementedif !is_surface(stl)
  @notimplementedif length(get_reffes(stl)) > 1
  @notimplementedif get_polytope(first(get_reffes(stl))) ≠ TRI
  D = 2
  topo = get_grid_topology(stl)
  v_to_e = get_faces(topo,0,1)
  v_to_f = get_faces(topo,0,D)
  f_to_v = get_faces(topo,D,0)
  graph = [ Int[] for _ in 1:num_vertices(stl) ]
  for v in 1:num_vertices(topo)
    f0 = first(v_to_f[v])
    fnext = f0
    while true
      i = findfirst(isequal(v),f_to_v[fnext])
      inext = i == D+1 ? 1 : i+1
      vnext = f_to_v[fnext][inext]
      push!(graph[v],vnext)
      _f = UNSET
      for f in v_to_f[vnext]
        if f ≠ fnext && v ∈ f_to_v[f]
          _f = f
          break
        end
      end
      fnext = _f
      fnext ≠ UNSET || break
      fnext ≠ f0 || break
    end
    if fnext == UNSET
      if length(graph[v]) < length(v_to_e[v])
        v0 = first(graph[v])
        fnext = f0
        while fnext ≠ UNSET
          i = findfirst(isequal(v),f_to_v[fnext])
          inext = i == 1 ? D+1 : i-1
          vnext = f_to_v[fnext][inext]
          pushfirst!(graph[v],vnext)
          _f = UNSET
          for f in v_to_f[vnext]
            if f ≠ fnext && v ∈ f_to_v[f]
              _f = f
              break
            end
          end
          fnext = _f
        end
        @assert length(graph[v]) == length(v_to_e[v])
      end
      push!(graph[v],OPEN)
    end
  end
  graph
end

function set_original_faces!(data::PolyhedronData,stl::DiscreteModel)
  Dc = num_dims(stl)
  topo = get_grid_topology(stl)
  v_to_f = data.vertex_to_original_faces
  for v in 1:length(v_to_f)
    for d in 0:Dc
      append!( v_to_f[v], map(i->i.+get_offset(topo,d),get_faces(topo,0,d)[v]) )
    end
  end
  data
end

function set_original_faces!(p::Polyhedron,a...)
  set_original_faces!(get_data(p),a...)
  p
end

function restrict(poly::Polyhedron,stl::DiscreteModel,stl_facets)
  nodes = Int[]
  for f in stl_facets
    for n in get_cell_nodes(stl)[f]
      if n ∉ nodes
        push!(nodes,n)
      end
    end
  end
  restrict(poly,nodes)
end

function restrict(p::Polyhedron,nodes)
  graph = get_graph(p)[nodes]
  f = i -> i ∈ nodes ? findfirst(isequal(i),nodes) : OPEN
  graph = map(i->map(f,i),graph)
  data = restrict(get_data(p),nodes)
  vertices = get_vertex_coordinates(p)[nodes]
  isopen = true
  Polyhedron(vertices,graph,isopen,data)
end

function restrict(data::PolyhedronData,nodes)
  v_to_Π = data.vertex_to_planes[nodes]
  v_to_of = data.vertex_to_original_faces[nodes]
  v_to_v = collect(1:length(nodes))
  Π_to_v_to_d = Vector{Int}[]
  Π_to_id = Int[]
  PolyhedronData( v_to_Π, v_to_of, v_to_v, Π_to_v_to_d, Π_to_id  )
end

function next_vertex(p::Polyhedron,vprevious::Integer,vcurrent::Integer)
  i = findfirst( isequal(vprevious), get_graph(p)[vcurrent] )
  i = ( i % length( get_graph(p)[vcurrent] ) ) + 1
  get_graph(p)[vcurrent][ i ]
end

function get_cell_planes(p::Polytope,pmin::Point,pmax::Point)
  @notimplementedif !is_n_cube(p)
  D = num_dims(p)
  [ CartesianPlane(isodd(i)*pmin+iseven(i)*pmax,D-((i-1)÷2),(-1)^i) for i in 1:2*D ]
end

function is_convex(stl::DiscreteModel{Dc,Dp},d::Integer,dface::Integer) where {Dc,Dp}
  @notimplementedif Dc ≠ Dp-1
  @notimplementedif d ≠ Dc-1
  facets = get_faces(get_grid_topology(stl),d,Dc)[dface]
  length(facets) == 2 || return false
  f1 = get_cell(stl,facets[1])
  f2 = get_cell(stl,facets[2])
  n1 = normal(f1)
  c1 = center(f1)
  Π1 = Plane(c1,n1)
  c2 = center(f2)
  signed_distance(c2,Π1) ≤ 0
end

function bisector_plane(stl::DiscreteModel{Dc,Dp},d::Integer,dface::Integer) where {Dc,Dp}
  @notimplementedif Dc ≠ Dp-1
  @notimplementedif d ≠ Dc-1
  facets = get_faces(get_grid_topology(stl),d,Dc)[dface]
  if length(facets) != 2 
    f1 = get_cell(stl,facets[1])
    return Plane(center(f1),normal(f1))
  end
  f1 = get_cell(stl,facets[1])
  f2 = get_cell(stl,facets[2])
  e = get_dface(stl,dface,Val{Dc-1}())
  n1 = normal(f1)
  n2 = normal(f2)
  n = (n1-n2)/2
  c = center(e)
  Plane(c,n)
end

function vertical_plane(
  stl::DiscreteModel{Dc,Dp},
  d::Integer,
  dface::Integer,
  dx::Integer) where {Dc,Dp}

  e = get_dface(stl,dface,Val{Dc-1}())
  c = center(e)
  o = tfill( 0, Val{Dp}() )
  vx = VectorValue( Base.setindex(o,1,dx) )

  # vectors = get_base(e)
  # orthogonal( vectors...,vx )
  @notimplementedif Dp ≠ 3
  v = e[2]-e[1]
  v /= norm(v)
  n = v × vx
  n /= norm(v)
  Plane(c,n)
end
  
function get_reflex_planes(stl::DiscreteModel;bisector=false) 
  Dc = num_dims(stl)
  Dp = num_point_dims(stl)
  T_Π = typeof( bisector_plane(stl,Dc-1,1) )
  Π_r = Vector{T_Π}(undef,num_faces(stl,Dc-1))
  for f in 1:num_faces(stl,Dc-1)
    if !bisector
      Π_r[f] = vertical_plane(stl,Dc-1,f,Dp)
    else
      Π_r[f] = bisector_plane(stl,Dc-1,f)
    end
  end
  Π_r
end

function get_facet_plane(stl::DiscreteModel,i)
  f = get_cell(stl,i)
  n = normal(f)
  c = center(f)
  Plane(c,n)
end

function get_facet_planes(stl::DiscreteModel)
  [ get_facet_plane(stl,i) for i in 1:num_cells(stl) ]
end

function get_convex_faces(stl::DiscreteModel)
  Dc = num_dims(stl)
  f_to_isconvex = zeros(Bool,num_faces(stl,Dc-1))
  for f in 1:num_faces(stl,Dc-1)
    f_to_isconvex[f] = is_convex(stl,Dc-1,f)
  end
  f_to_isconvex
end

function get_disconnected_parts(poly::Polyhedron)
  v_to_part = fill(UNSET,num_vertices(poly))
  stack = Int[]
  num_parts = 0
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    v_to_part[v] == UNSET || continue
    num_parts += 1
    empty!(stack)
    push!(stack,v)
    while !isempty(stack)
      vcurrent = pop!(stack)
      for vneig in get_graph(poly)[vcurrent]
        if vneig ∉ (UNSET,OPEN) && v_to_part[vneig] == UNSET
          v_to_part[vneig] = num_parts
          push!(stack,vneig)
        end
      end
    end
  end
  v_to_part,num_parts
end


function get_disconnected_faces(poly::Polyhedron,stl::DiscreteModel,d::Integer)
  v_to_f = get_data(poly).vertex_to_original_faces
  facedims = get_facedims(get_grid_topology(stl))
  v_to_part,num_parts = get_disconnected_parts(poly)
  part_to_faces = [ Int[] for _ in 1:num_parts ]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    part = v_to_part[v]
    for f in v_to_f[v]
      if facedims[f] == d && f ∉ part_to_faces[part]
        push!(part_to_faces[part],f)
      end
    end
  end
  part_to_faces
end

