
## Constants

const OPEN = -1

struct Polyhedron{Dp,Tp,Td}
  vertices::Vector{Point{Dp,Tp}}
  edge_vertex_graph::Vector{Vector{Int}}
  data::Td
end

struct PolyhedronData
  vertex_to_planes::Vector{Vector{Int}}
  vertex_to_original_faces::Vector{Vector{Int}}
  vertex_to_isopen::Vector{Bool}
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
  p = Polyhedron(X,𝓖,isopen=true)
  set_original_faces!(p,stl)
  p
end

function Polyhedron(
  vertices::AbstractVector{<:Point},
  graph::Vector{<:Vector}
  ;isopen=false::Bool)

  Polyhedron(vertices,graph,polyhedron_data(length(vertices);isopen))
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

## TMP

#struct DistancePlaneList{T,V<:AbstractVector{Vector{T}}}
#  plane_to_vertex_to_distance::V
#end
#
#struct DistancePlane{T}
#  i::Int
#  list::DistancePlaneList{T}
#end
#
#Base.length(a::DistancePlaneList) = length(a.plane_to_vertex_to_distance)
#
#num_planes(a::DistancePlaneList) = length(a)
#
#num_vertices(a::DistancePlaneList) = length(a.plane_to_vertex_to_distance[end])
#
#Base.getindex(a::DistancePlaneList,i::Integer) = DistancePlane(i,a)
#
#function Base.getindex(a::DistancePlaneList,indices::AbstractVector)
#  dists = [ copy(get_distances(a)[i]) for i in indices ]
#  DistancePlaneList(dists)
#end
#
#
#function Base.view(a::DistancePlaneList,indices)
#  dists = view( get_distances(a), indices )
#  DistancePlaneList(dists)
#end
#
#function Base.iterate(a::DistancePlaneList,state=1)
#  if state > length(a)
#    nothing
#  else
#    (a[state],state+1)
#  end
#end
#
#get_distances(a::DistancePlaneList) = a.plane_to_vertex_to_distance
#
#get_distances(a::DistancePlane) = get_distances(a.list)[a.i]
#
#num_vertices(a::DistancePlane) = length(get_distances(a))
#
#function precompute_planes(p::Polyhedron,Πr,stl_e)
#  dists = compute_distances(p,Πr,stl_e)
#  DistancePlaneList(dists)
#end
#
#function add_vertex!(a::DistancePlane,v1,v2)
#  add_vertex!(get_distances(a.list),a.i,v1,v2)
#end
#
#function add_vertex!(a::DistancePlaneList,v1,v2,d1,d2)
#  add_vertex!(get_distances(a),UNSET,v1,v2,d1,d2)
#end
#
#function compute_distances(p::Polyhedron,Π,stl_edges)
#  v_to_e = get_data(p).vertex_to_original_faces
#  Π_to_v_to_d = [ zeros(num_vertices(p)) for _ in 1:length(Π) ]
#  for (i,Πi) in enumerate(Π)
#    for v in 1:num_vertices(p)
#      if isactive(p,v)
#        dist = signed_distance(p[v],Πi)
#        if stl_edges[i] ∈ v_to_e[v]
#          @assert abs(dist) < 1e-9
#          dist = 0.0
#        end
#        Π_to_v_to_d[i][v] = dist
#      end
#    end
#  end
#  Π_to_v_to_d
#end
#
#function add_vertex!(Π_to_v_to_d,Πid,v1::Integer,v2::Integer,d1=nothing,d2=nothing)
#  d1 = !isnothing(d1) ? d1 : Π_to_v_to_d[Πid][v1]
#  d2 = !isnothing(d2) ? d2 : Π_to_v_to_d[Πid][v2]
#  for i in 1:length(Π_to_v_to_d)
#    distances = Π_to_v_to_d[i] 
#    d = Πid == i ? 0.0 : ( d1*distances[v2] - d2*distances[v1] ) / ( d1 - d2 )
#    push!(distances,d)
#  end
#  Π_to_v_to_d
#end
#
#function Base.cat(a::DistancePlaneList,b::DistancePlaneList)
#  @assert length(a) == length(b)
#  dists = [ [get_distances(a)[i];get_distances(b)[i]] for i in 1:length(a) ]
#  DistancePlaneList( dists )
#end
#
#is_precomputed(::Plane) = false
#
#is_precomputed(::CartesianPlane) = false
#
#is_precomputed(::DistancePlane) = true

## Polyhedra operations

function clip(poly::Polyhedron,Π)
  p = poly
  for Πi in Π
    p,_ = split(p,Πi)
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
  add_open_vertices!(in_graph,data,num_vertices(p))
  add_open_vertices!(out_graph,data,num_vertices(p))
  vertices = [p.vertices;new_vertices]
  in_data = data
  out_data = deepcopy(data)
  Polyhedron(vertices,in_graph,in_data), Polyhedron(vertices,out_graph,out_data)
end

function Base.cat(p1::Polyhedron,p2::Polyhedron)
  vertices = vcat( get_vertex_coordinates(p1), get_vertex_coordinates(p2) )
  graph = cat_graph(get_graph(p1),get_graph(p2))
  data = cat_data( get_data(p1), get_data(p2) )
  Polyhedron(vertices,graph,data)
end

function flip(p::Polyhedron)
  graph = flip_graph(get_graph(p))
  Polyhedron(get_vertex_coordinates(p),graph,get_data(p))
end

function flip_graph(graph)
  map( i -> reverse(i), graph )
end

function _connect(p::Polyhedron{D}) where D
  data = get_data(p)
  v_to_Π = get_vertex_to_planes(data)
  e_to_Π = NTuple{D-1,Int}[]
  e_to_v = Vector{Int}[]
  for v in 1:num_vertices(p)
    isactive(p,v) || continue
    if length(v_to_Π[v]) ≥ D-1
      for i in 1:length(v_to_Π[v]), j in i+1:length(v_to_Π[v])
        Πi = v_to_Π[v][i]
        Πj = v_to_Π[v][j]
        Πe = Πi < Πj ? (Πi,Πj) : (Πj,Πi)
        if Πe ∉ e_to_Π
          push!(e_to_Π,Πe)
          push!(e_to_v,[])
        end
        e = findfirst(isequal(Πe),e_to_Π)
        push!(e_to_v[e],v)
      end
    end
  end
  v_to_isopen = data.vertex_to_isopen
  cut_edge = UNSET
  min_num_cuts = typemax(Int)
  for e in 1:length(e_to_v)
    length(e_to_v[e]) ≥ 2 || continue
    sort!(e_to_v[e], by = i -> p[i] ) 
    @show v_to_isopen[ e_to_v[e] ]
    if v_to_isopen[ first(e_to_v[e]) ]
      reverse!(e_to_v[e],1,2)
      @assert !v_to_isopen[ first(e_to_v[e]) ]
    end
    if v_to_isopen[ last(e_to_v[e]) ]
      reverse!(e_to_v[e],length(e_to_v[e])-1,length(e_to_v[e]))
      @assert !v_to_isopen[ last(e_to_v[e]) ]
    end
    if length(e_to_v[e]) > 2
      if length(e_to_v[e]) < min_num_cuts 
        min_num_cuts = length(e_to_v[e])
        cut_edge = e
      end
    end
  end
  @assert cut_edge ≠ UNSET # -> return only surface
  isinside = falses(num_vertices(p))
  istouch = falses(num_vertices(p))
  e1 = e_to_v[cut_edge][1]
  e2 = e_to_v[cut_edge][end]
  v = e_to_v[cut_edge][2]
  vstart = e1
  isinside[vstart] = _is_edge_oriented(p,v,e1,e2)
  istouch[vstart] = true
  stack = Int[]
  push!(stack,vstart)
  while !isempty(stack)
    v = pop!(stack)
    for vneig in get_graph(p)[v]
      e = findfirst( i -> (v,vneig) ⊆ i, e_to_v)
      if !istouch[vneig] 
        @show length(e_to_v[e])
        isinside[vneig] = iseven(length(e_to_v[e])) == isinside[v]
        istouch[vneig] = true
        push!(stack,vneig)
      end
    #  @assert  isinside[vneig] == (iseven(length(e_to_v[e])) == isinside[v])
    end
  end
  @show istouch, isinside
  graph = get_graph(p)
  for e in 1:length(e_to_v)
    if length(e_to_v[e]) > 2
      e1 = e_to_v[e][1]
      istart = isinside[e1] ? 1 : 2
      for i in istart:2:length(e_to_v[e])-1
        v1 = e_to_v[e][i]
        v2 = e_to_v[e][i+1]
        v1_neig = i == 1 ? e2 : UNSET
        v2_neig = i == length(e_to_v[e])-1 ? e1 : UNSET
        link!(graph,v1,v1_neig,v2,v2_neig)
      end
    end
  end
  for v in 1:num_vertices(p)
    if !isinside[v] && istouch[v]
      empty!(graph[v])
    end
  end
  @show e_to_v
  @show e_to_Π
  p
end

function _is_edge_oriented(poly::Polyhedron,v::Integer,e1::Integer,e2::Integer)
  data = get_data(poly)
  v_to_Π = get_vertex_to_planes(data)
  vnext = next_vertex(poly,UNSET,v)
  fs = v_to_Π[v][ findfirst( in(v_to_Π[vnext]), v_to_Π[v] ) ] 
  vnext = next_vertex(poly,e1,e2)
  fc = v_to_Π[e2][ findfirst( in(v_to_Π[vnext]), v_to_Π[e2] ) ]
  fs == fc
end


function merge(p::Polytope,cell::Polyhedron,poly::Polyhedron)
  D = num_dims(p)
  data = get_data(poly)
  v_to_Π = get_vertex_to_planes(data)
  v_to_df = fill(UNSET,num_vertices(poly))
  e_to_v = [ Int[] for _ in 1:num_edges(p) ]
  for v in 1:num_vertices(poly)
    !isempty(get_graph(poly)[v]) || continue
    d = D - length( v_to_Π[v] )
    for f in v_to_Π[v], df in get_faces(p,D-1,d)[f]
      if get_faces(p,d,D-1)[df] ⊆ v_to_Π[v]
        v_to_df[v] = get_dimrange(p,d)[df]
        if d == 1
          push!(e_to_v[df],v)
        end
        break
      end
    end
  end
  graph_in = cat_graph(get_graph(cell),get_graph(poly))
  graph_out = cat_graph(get_graph(cell),get_graph(poly),revert=true)
  vertices = vcat( get_vertex_coordinates(cell), get_vertex_coordinates(poly) )
  data_in = cat_data( get_data(cell), get_data(poly) )
  data_out = deepcopy(data_in)
  min_num_cuts = typemax(0)
  cut_edge = UNSET
  for e in 1:num_edges(p)
    num_cuts = length(e_to_v[e])
    if num_cuts > 0 && num_cuts < min_num_cuts
      min_num_cuts = num_cuts
      cut_edge = e
    end
  end
  if cut_edge ≠ UNSET
    for e in 1:num_edges(p)
      if length(e_to_v[e]) > 1
        d = ((e-1)>>(D-1)) + 1
        sort!(e_to_v[e], by = i -> poly[i][d] ) 
      end
    end
    isinside = falses(num_vertices(cell)) # These can be stored on Int
    istouch = falses(num_vertices(cell))
    vstart = get_face_vertices(p,1)[cut_edge][1]
    isinside[vstart] = is_edge_oriented(p,cell,cut_edge,poly,first(e_to_v[cut_edge]))
    istouch[vstart] = true
    stack = Int[]
    push!(stack,vstart)
    while length(stack) > 0
      vcurrent = pop!(stack)
      for eneig in get_faces(p,0,1)[vcurrent]
        for vneig in get_faces(p,1,0)[eneig]
          if vneig ≠ vcurrent && !istouch[vneig]
            istouch[vneig] = true
            isinside[vneig] = isinside[vcurrent] == iseven(length(e_to_v[eneig]))
            push!(stack,vneig)
          end
        end
      end
    end
    for e in 1:num_edges(p)
      if length(e_to_v[e]) > 0
        edge = get_face_vertices(p,1)[e]
        istarts = (0,1)
        istarts = isinside[edge[1]] ? istarts : reverse(istarts)
        graphs = (graph_in,graph_out)
        for (istart,graph) in zip(istarts,graphs)
          for i in istart:2:length(e_to_v[e])
            v1 = i == 0 ? edge[1] : e_to_v[e][i] + num_vertices(cell)
            v2 = i == length(e_to_v[e]) ? edge[2] : e_to_v[e][i+1] + num_vertices(cell)
            v1_neig = i == 0 ? edge[2] : UNSET
            v2_neig = i == length(e_to_v[e]) ? edge[1] : UNSET
            link!(graph,v1,v1_neig,v2,v2_neig)
          end
        end
        data_in.vertex_to_isopen[first(e_to_v[e])+num_vertices(p)] = false
        data_in.vertex_to_isopen[last(e_to_v[e])+num_vertices(p)] = false
        data_out.vertex_to_isopen[first(e_to_v[e])+num_vertices(p)] = false
        data_out.vertex_to_isopen[last(e_to_v[e])+num_vertices(p)] = false
      end
    end
  else
    isinside = fill( is_cell_inside(p,cell,poly), num_vertices(cell) )
  end
  ## Diconnect unused vertices
  for v in 1:num_vertices(p)
    if !isinside[v]
      empty!(graph_in[v])
    end
  end
  for v in 1:num_vertices(p)
    if isinside[v]
      empty!(graph_out[v])
    end
  end
  delete_open_vertices!(graph_in,data_in)
  delete_open_vertices!(graph_out,data_out)
  #update_open_data!(graph_in,data_in)
  #update_open_data!(graph_out,data_out)
  Polyhedron(vertices,graph_in,data_in), Polyhedron(vertices,graph_out,data_out)
end

function decompose(surf::Polyhedron,cell::Polyhedron,rfaces)
  i = findfirst(i->has_stl_edge(surf,i), rfaces )
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


#
#    if length(stl_edges) == i
#      R = T
#    else
#      R = empty(T)
#      stl_e = view(stl_edges,i+1:length(stl_edges))
#      for k in T
#        Π_n = Π[i+1:length(stl_edges)]
#        Tn = _decompose(k,stl_e,Π_n,e_to_isconvex;isinside)
#        append!(R,Tn)
#      end
#    end
#  end



function decompose(poly::Polyhedron,Πr,e_to_isconvex;isinside)
  stl_edges = get_stl_edges(poly)
  decompose(poly,stl_edges,Πr,e_to_isconvex;isinside)
end

function decompose(poly::Polyhedron,stl_edges,Πr,e_to_isconvex;isinside)
  i = findfirst(i -> (e_to_isconvex[i] != isinside) && has_stl_edge(poly,i), stl_edges )
  if i === nothing 
    R = [poly]
  else
    ei = stl_edges[i]
    Πi = Πr[ei]
    P = split(poly,Πi,-ei)
    T = typeof(poly)[]
    for k in P
      !isnothing(k) || continue
      connect!(k,-ei)
      push!(T,k)
    end
    if length(stl_edges) == i
      R = T
    else
      R = empty(T)
      stl_e = view(stl_edges,i+1:length(stl_edges))
      for k in T
        Tn = decompose(k,stl_e,Πr,e_to_isconvex;isinside)
        append!(R,Tn)
      end
    end
  end
  R
end

function _decompose(poly::Polyhedron,stl_edges,Π,e_to_isconvex;isinside)
  i = findfirst(i -> (e_to_isconvex[i] != isinside) && has_stl_edge(poly,i), stl_edges )
  if i === nothing 
    R = [poly]
  else
    ei = stl_edges[i]
    Πi = Π[i]
    p⁻,p⁺ = split(poly,Πi,-ei)
    T = [p⁻,p⁺]
    @assert nothing ∉ T
    if length(stl_edges) == i
      R = T
    else
      R = empty(T)
      stl_e = view(stl_edges,i+1:length(stl_edges))
      for k in T
        Π_n = Π[i+1:length(stl_edges)]
        Tn = _decompose(k,stl_e,Π_n,e_to_isconvex;isinside)
        append!(R,Tn)
      end
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

function connect!(p::Polyhedron,Πid)
  !is_connected(p) || return p
  data = get_data(p)
  v_to_Π = get_vertex_to_planes(data)
  e_to_Π = Int[]
  e_to_v = Vector{Int}[]
  for v in 1:num_vertices(p)
    if isactive(p,v) && Πid ∈ v_to_Π[v]
      for Πe ∈ v_to_Π[v]
        if Πe ≠ Πid && Πe ∈ e_to_Π
          i = findfirst(isequal(Πe),e_to_Π)
          push!(e_to_v[i],v)
        else
          push!(e_to_Π,Πe)
          push!(e_to_v,[v])
        end
      end
    end
  end
  graph = get_graph(p)
  for e in 1:length(e_to_v)
    if length(e_to_v[e]) > 2
      @notimplementedif isodd( length(e_to_v[e]) )
      vs = sort!( e_to_v[e], by = i -> p[i] )
     # link!(graph,vs[1],vs[end],vs[2],UNSET)
     # link!(graph,vs[end],vs[1],vs[end-1],UNSET)
      for i in 1:2:length(vs)-1
        v0 =  i == 1 ? vs[end] : UNSET
        v1 = i+1 == length(vs) ? vs[1] : UNSET
        link!(graph,vs[i],v0,vs[i+1],v1)
      end
      delete_open_vertices!(graph,vs)
    end
  end
  update_open_data!(graph,data)
  p
end


## Getters

num_dims(::Polyhedron{D}) where D = D

num_vertices(a::Polyhedron) = length(a.vertices)

get_vertex_coordinates(a::Polyhedron) = a.vertices

Base.getindex(a::Polyhedron,i::Integer) = a.vertices[i]

get_graph(a::Polyhedron) = a.edge_vertex_graph

get_data(a::Polyhedron) = a.data

get_vertex_to_planes(a::PolyhedronData) = a.vertex_to_planes

function isactive(p::Polyhedron,vertex::Integer)
  !isempty( get_graph(p)[vertex] )
end

function is_connected(p::Polyhedron)
  for v in 1:num_vertices(p)
    if isactive(p,v) && is_open_part(p,v)
      return false
    end
  end
  true
end

function get_reflect_faces(p::Polyhedron,stl::DiscreteModel)
  @notimplementedif num_dims(p) ≠ 3
  Dr = num_dims(stl)-1
  stack = Int[]
  istouch = zeros(Bool,num_vertices(p))
  v_to_f = get_data(p).vertex_to_original_faces
  facedims = get_facedims(get_grid_topology(stl))
  rf = Int[]
  for v in 1:num_vertices(p)
    if !istouch[v] && isactive(p,v)
      empty!(stack)
      push!(stack,v)
      while !isempty(stack)
        vcurrent = pop!(stack)
        !istouch[vcurrent] || continue
        istouch[vcurrent] = true
        for vneig in get_graph(p)[vcurrent]
          vneig ∉ (OPEN,UNSET) || continue
          if !istouch[vneig]
            push!(stack,vneig)
            for f in v_to_f[vcurrent]
              if f ∈ v_to_f[vneig] && f ∉ rf && facedims[f] == Dr
                push!(rf,f)
              end
            end
          end
        end
      end
    end
  end
  rf
end

function has_stl_edge(p::Polyhedron,edge::Integer)
  v_to_f = get_data(p).vertex_to_original_faces
  for v in 1:num_vertices(p)
    if isactive(p,v) && edge ∈ v_to_f[v]
      for vneig in get_graph(p)[v]
        vneig ≠ OPEN || continue
        if edge ∈ v_to_f[vneig] && p[v] ≠ p[vneig]
          return true
        end
      end
    end
  end
  false
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
  i = 0
  for f in faces
    if has_original_face(p,f,Val{d}())
      i += 1
      faces[i] = f
    end
  end
  resize!(faces,i)
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
  for v in 1:num_vertices(p)
    if isactive(p,v) && face ∈ v_to_f[v]
      for vneig in get_graph(p)[v]
        vneig ∉ (OPEN,UNSET) || continue
        if face ∈ v_to_f[vneig] && p[v] ≠ p[vneig]
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
  for v in 1:num_vertices(p)
    if isactive(p,v) && face ∈ v_to_f[v]
      for vneig in get_graph(p)[v]
        vneig ∉ (OPEN,UNSET) || continue
        if face ∈ v_to_f[vneig]
          num_v = 1
          vcurrent = v
          vnext = vneig
          while vnext ≠ v 
            (p[vcurrent] == p[vnext] || p[v] == p[vnext]) || (num_v += 1)
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

function polyhedron_data(num_vertices::Integer;isopen=false)
  v_to_Π = [ Int[] for _ in 1:num_vertices ]
  v_to_of = [ Int[] for _ in 1:num_vertices ]
  v_to_o = fill(isopen,num_vertices)
  Π_to_v_to_d = Vector{Int}[]
  Π_to_id = Int[]
  PolyhedronData( v_to_Π, v_to_of, v_to_o, Π_to_v_to_d, Π_to_id  )
end

function polyhedron_data(p::Polytope)
  v_to_Π = deepcopy( get_faces(p,0,num_dims(p)-1) )
  v_to_of = [ Int[] for _ in 1:num_vertices(p) ]
  v_to_o = fill(false,num_vertices(p))
  Π_to_v_to_d = Vector{Int}[]
  Π_to_id = Int[]
  PolyhedronData( v_to_Π, v_to_of, v_to_o, Π_to_v_to_d, Π_to_id  )
end

function cat_data(d1::PolyhedronData,d2::PolyhedronData)
  v_to_Π = [ d1.vertex_to_planes; d2.vertex_to_planes ]
  v_to_of = [ d1.vertex_to_original_faces; d2.vertex_to_original_faces ]
  v_to_o = [ d1.vertex_to_isopen; d2.vertex_to_isopen ]
  PolyhedronData( v_to_Π, v_to_of, v_to_o )
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

function add_open_vertices!(graph,data::PolyhedronData,num_vertices::Integer)
  v_to_isopen = data.vertex_to_isopen
  v_to_Π = data.vertex_to_planes
  for v in num_vertices+1:length(graph)
    if v_to_isopen[v] && UNSET ∉ graph[v]
      push!(graph[v],graph[v][end])
      graph[v][end-1] = OPEN
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
  v_to_o = data.vertex_to_isopen
  Πs = intersect( v_to_Π[v1], v_to_Π[v2] )
  fs = intersect( v_to_f[v1], v_to_f[v2] )
  @assert v_to_o[v1] == v_to_o[v1]
  oc = v_to_o[v1]
  push!( Πs, Πid )
  push!( v_to_Π, Πs )
  push!( v_to_f, fs )
  push!( v_to_o, oc )
  Π_to_v_to_d = get_plane_distances(data)
  Π_to_id = get_plane_ids(data)
  i = findfirst(isequal(Πid),Π_to_id)
  d1 = Π_to_v_to_d[i][v1]
  d2 = Π_to_v_to_d[i][v2]
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

function compute_graph(stl::DiscreteModel,face_list)
  nodes = Int[]
  for face in face_list
    for node in get_cell_nodes(stl)[face]
      if node ∉ nodes
        push!(nodes,node)
      end
    end
  end
  num_vertices = length(nodes)
  graph = [ Int[] for _ in 1:num_vertices ]
  for face in face_list
    face_nodes = get_cell_nodes(stl)[face]
    for (i,n) in enumerate( face_nodes )
      n⁺ = face_nodes[ i < length(face_nodes) ? i+1 : 1 ]
      n⁻ = face_nodes[ i > 1 ? i-1 : length(face_nodes) ]
      v = findfirst( isequal(n), nodes )
      v⁺ = findfirst( isequal(n⁺), nodes )
      v⁻ = findfirst( isequal(n⁻), nodes )
      if v⁻ ∉ graph[v]
        push!(graph[v],v⁻)
      end
      if v⁺ ∉ graph[v]
        push!(graph[v],v⁺)
      end
    end
  end
  orient_graph!(graph)
  set_open_nodes!(graph)
  graph, nodes
end

function orient_graph!(graph)
  for v in 1:length(graph)
    if length(graph[v]) > 3
      i = 2
      while i < length(graph[v])
        v_i = graph[v][i]
        found = false
        for v_j in graph[v_i]
          if v_j ≠ graph[v][i-1] && v_j ∈ graph[v]
            k = findfirst( isequal(v_j) , graph[v] )
            graph[v][k] = graph[v][i+1]
            graph[v][i+1] = v_j
            found = true
            break
          end
        end
        found || break
        i += 1
      end
      if i < length(graph[v])
        i_last = i
        i = 1
        inext = length(graph[v])
        iprev = 2
        while inext > i_last
          v_i = graph[v][i]
          found = false
          for v_j in graph[v_i]
            ip = i == length(graph[v]) ? 1 : i+1
            if v_j ≠ graph[v][iprev] && v_j ∈ graph[v]
              k = findfirst( isequal(v_j) , graph[v] )
              graph[v][k] = graph[v][inext]
              graph[v][inext] = v_j
              found = true
              break
            end
          end
          found || break
          i = i == 1 ? length(graph[v]) : i-1
          iprev = i == length(graph[v]) ? 1 : i-1
          inext -= inext
        end
      end
    end
  end
  graph
end

function set_open_nodes!(graph)
  for v in 1:length(graph)
    @notimplementedif length( graph[v] ) < 2
    if length( graph[v] ) == 2
      push!( graph[v], OPEN )
    else
      ibreak = UNSET
      for (i,v_i) in enumerate( graph[v] )
        inext = ( i % length( graph[v] ) ) + 1
        j = findfirst( isequal(v), graph[v_i] )
        jnext = j == 1 ? length( graph[v_i] ) : j-1
        if graph[v][inext] != graph[v_i][jnext]
          ibreak = i
          break
        end
      end
      if ibreak ≠ UNSET
        push!(graph[v],UNSET)
        for i in reverse(ibreak+1:length(graph[v])-1)
          graph[v][i+1] = graph[v][i]
        end
        graph[v][ibreak+1] = OPEN
      end
    end
  end
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
  Polyhedron(vertices,graph,data)
end

function restrict(data::PolyhedronData,nodes)
  v_to_Π = data.vertex_to_planes[nodes]
  v_to_of = data.vertex_to_original_faces[nodes]
  v_to_o = data.vertex_to_isopen[nodes]
  Π_to_v_to_d = Vector{Int}[]
  Π_to_id = Int[]
  PolyhedronData( v_to_Π, v_to_of, v_to_o, Π_to_v_to_d, Π_to_id  )
end

function next_vertex(p::Polyhedron,vprevious::Integer,vcurrent::Integer)
  i = findfirst( isequal(vprevious), get_graph(p)[vcurrent] )
  i = ( i % length( get_graph(p)[vcurrent] ) ) + 1
  get_graph(p)[vcurrent][ i ]
end

function orient_edge_endpoints(
  p::Polytope,
  cell::Polyhedron,
  edge::Integer,
  poly::Polyhedron,
  vertex::Integer)

  vnext = next_vertex(poly,UNSET,vertex)
  data = get_data(poly)
  v_to_Π = get_vertex_to_planes(data)
  f = v_to_Π[vertex][ findfirst( in(v_to_Π[vnext]), v_to_Π[vertex] ) ] # f[next] ∩ f[v]
  e = get_faces(p,1,0)[edge]
  vnext = next_vertex(cell,e[1],e[2])

  vnext ∈ get_face_vertices(p,num_dims(p)-1)[f] ? (e[1],e[2]) : (e[2],e[1])
end

function is_edge_oriented(
  p::Polytope,
  cell::Polyhedron,
  edge::Integer,
  poly::Polyhedron,
  vertex::Integer)

  vnext = next_vertex(poly,UNSET,vertex)
  data = get_data(poly)
  v_to_Π = get_vertex_to_planes(data)
  f = v_to_Π[vertex][ findfirst( in(v_to_Π[vnext]), v_to_Π[vertex] ) ] # f[next] ∩ f[v]
  e = get_faces(p,1,0)[edge]
  vnext = next_vertex(cell,e[1],e[2])
  get_face_vertices(p,num_dims(p)-1)[f] ∋ vnext
end

function link!(graph,v1,vnext1,v2,vnext2)
  i1 = findfirst( isequal(vnext1), graph[v1] )
  i2 = findfirst( isequal(vnext2), graph[v2] )
  graph[v1][i1] = v2
  graph[v2][i2] = v1
end

function is_cell_inside(p::Polytope{D},cell::Polyhedron{D},poly::Polyhedron{D}) where D
  data = get_data(poly)
  v_to_Π = get_vertex_to_planes(data)
  vstart = UNSET
  Πid = UNSET
  # TODO: select vstart s.t. min_length is maximum (for more than one cut face)
  for v in 1:num_vertices(poly)
    if !isempty(get_graph(poly)[v]) && !isempty(v_to_Π[v])
      @assert length(v_to_Π[v]) == 1
      Πid = first(v_to_Π[v])
      vstart = v
      break
    end
  end
  d_Π = D-((Πid-1)÷2)
  i = 1
  pmin,pmax = poly[vstart],poly[vstart]
  mins,maxs = tfill(i,Val{D}()), tfill(i,Val{D}())
  vnext = UNSET
  for vneig in reverse( get_graph(poly)[vstart] )
    if vneig ≠ OPEN && !isempty(v_to_Π[vneig]) && first(v_to_Π[vneig]) == Πid
      vnext = vneig
      break
    end
  end
  vcurrent = vstart
  while vnext ≠ vstart
    i += 1
    _vnext = next_vertex(poly,vcurrent,vnext)
    vcurrent = vnext
    vnext = _vnext
    @assert first(v_to_Π[vnext]) == Πid
    for d in 1:D
      if d ≠ d_Π 
        if poly[vcurrent][d] < pmin[d]
          pmin = Base.setindex(pmin,poly[vcurrent][d],d)
          mins = Base.setindex(mins,i,d)
        end
        if poly[vcurrent][d] > pmax[d]
          pmax = Base.setindex(pmax,poly[vcurrent][d],d)
          maxs = Base.setindex(maxs,i,d)
        end
      end
    end
  end
  d1 = d_Π ≠ 1 ? 1 : 2
  d2 = d_Π ≠ d1+1 ? d1+1 : d1+2
  max_i = max( max(mins...), max(maxs...) ) 
  min_i = min( mins[d1], mins[d2], maxs[d1], maxs[d2] )

  # TODO: rewrite this clearer: xmin-ymax-xmax-ymin -> clockwise
  if mins[d1] ≠ mins[d2]
    if mins[d1] == max_i
      cut_orientation = mins[d2] != min_i
    elseif mins[d1] == min_i
      cut_orientation = mins[d2] == max_i
    else
      cut_orientation = mins[d1] > mins[d2]
    end
  elseif maxs[d1] ≠ maxs[d2]
    if maxs[d1] == max_i
      cut_orientation = maxs[d2] == min_i
    elseif maxs[d1] == min_i
      cut_orientation = maxs[d2] != max_i
    else
      cut_orientation = maxs[d1] < maxs[d2]
    end
  end

  e0 = first( get_faces(p,D-1,1)[Πid] )
  d_e = ((e0-1)>>(D-1)) + 1
  @assert d1 == d_e

  e = get_faces(p,1,0)[e0]
  vnext = next_vertex(cell,e[1],e[2])
  cell_orientation = vnext ∉ get_face_vertices(p,D-1)[Πid]

  cut_orientation != cell_orientation 
end

function get_cell_planes(p::Polytope,pmin::Point,pmax::Point)
  @notimplementedif !is_n_cube(p)
  D = num_dims(p)
  [ CartesianPlane(isodd(i)*pmin+iseven(i)*pmax,D-((i-1)÷2),(-1)^i) for i in 1:2*D ]
end

function delete_open_vertices!(graph,p::Polytope)
  delete_open_vertices!(graph,1:num_vertices(p))
end

function delete_open_vertices!(graph,vrange::AbstractVector)
  stack = Int[]
  istouch = zeros(Bool,length(graph))
  for v in vrange
    if !istouch[v] && !isempty(graph[v])
      empty!(stack)
      push!(stack,v)
      while !isempty(stack)
        vcurrent = pop!(stack)
        for vneig in graph[vcurrent]
          if vneig ≠ OPEN && !istouch[vneig]
            istouch[vneig] = true
            push!(stack,vneig)
          end
        end
      end
    end
  end
  if any( istouch )
    for v in 1:length(graph)
      if istouch[v] && OPEN ∈ graph[v]
        i = findfirst(isequal(OPEN),graph[v])
        deleteat!(graph[v],i)
      end
    end
  else
    for v in 1:length(graph)
      if OPEN ∈ graph[v]
        i = findfirst(isequal(OPEN),graph[v])
        deleteat!(graph[v],i)
      end
    end
  end
  graph
end

function delete_open_vertices!(graph,data::PolyhedronData)
  stack = Int[]
  v_to_isopen = data.vertex_to_isopen
  openpoly = true
  for v in 1:length(graph)
    if !v_to_isopen[v] && !isempty(graph[v])
      openpoly = false 
      empty!(stack)
      push!(stack,v)
      while !isempty(stack)
        vcurrent = pop!(stack)
        for vneig in graph[vcurrent]
          if vneig ≠ OPEN && v_to_isopen[vneig]
            v_to_isopen[vneig] = false
            push!(stack,vneig)
          end
        end
      end
    end
  end
  if !openpoly
    for v in 1:length(graph)
      if !v_to_isopen[v] && OPEN ∈ graph[v]
        i = findfirst(isequal(OPEN),graph[v])
        deleteat!(graph[v],i)
      end
    end
  else
    for v in 1:length(graph)
      if OPEN ∈ graph[v]
        i = findfirst(isequal(OPEN),graph[v])
        deleteat!(graph[v],i)
      end
    end
    fill!(v_to_isopen,false)
  end
  graph
end

function update_open_data!(graph,data)
  stack = Int[]
  openpoly = false
  v_to_isopen = data.vertex_to_isopen
  for v in 1:length(graph)
    if !v_to_isopen[v]
      empty!(stack)
      push!(stack,v)
      while !isempty(stack)
        vcurrent = pop!(stack)
        for vneig in graph[vcurrent]
          if v_to_isopen[vneig]
            v_to_isopen[vneig] = false
            push!(stack,vneig)
          end
        end
      end
    elseif OPEN ∈ graph[v]
      openpoly = true
    end
  end
  if !openpoly
    fill!(v_to_isopen,false)
  end
  graph
end

function is_open_part(p::Polyhedron,vertex::Integer)
  data = get_data(p)
  data.vertex_to_isopen[vertex]
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

function cat_graph(G1,G2;revert=false)
  f = i -> map( j -> j ∈ (UNSET,OPEN) ? j : j+length(G1), i )
  if revert
    r = i -> reverse(i)
    f = r∘f
  end
  [ deepcopy(G1); map(f,G2) ]
end

