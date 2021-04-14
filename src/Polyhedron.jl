
## Constants

const OPEN = -1

const FACE_CUT = -1
const FACE_IN = 1
const FACE_OUT = 2

struct Polyhedron{Dp,Tp,Td}
  vertices::Vector{Point{Dp,Tp}}
  edge_vertex_graph::Vector{Vector{Int32}}
  isopen::Bool # TODO: revise, it may mean is_solid/is_surface
  data::Td
end

struct PolyhedronData
  vertex_to_planes::Vector{Vector{Int32}}
  vertex_to_original_faces::Vector{Vector{Int32}}
  vertex_to_parent_vertex::Vector{Int32}
  vertex_to_parent_edge::Vector{Tuple{Int32,Int32}}
  plane_to_vertex_to_distances::Vector{Vector{Float64}}
  plane_to_ref_plane::Vector{Int32}
  plane_to_ids::Vector{Int32}
end

struct CartesianPlane{D,T}
  d::Int8
  value::T
  positive::Bool
end

## Constructors

function Polyhedron(stl::STL)
  𝓖 = compute_graph(stl)
  X = get_vertex_coordinates(stl)
  isopen = is_open_surface(stl)
  p = Polyhedron(X,𝓖;isopen)
  set_original_faces!(p,stl)
  p
end

function Polyhedron(
  vertices::AbstractVector{<:Point},
  graph::Vector{<:Vector},
  isopen,
  data)

  Polyhedron(collect(vertices),graph,isopen,data)
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
  e_v_graph = map(i->Int32.(i),e_v_graph)
  data = polyhedron_data(p)
  Polyhedron(vertices,e_v_graph,data)
end

function Polyhedron(p::Polytope{3},vertices::AbstractVector{<:Point})
  if p == TET
    e_v_graph = [[2,4,3],[3,4,1],[1,4,2],[1,2,3]]
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
  e_v_graph = map(i->Int32.(i),e_v_graph)
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

function check_graph(p::Polyhedron)
  check_graph(get_graph(p))
end

function check_graph(graph::AbstractVector{<:AbstractVector})
  for v in 1:length(graph)
    !isempty(graph[v]) || continue
    for vneig in graph[v]
      vneig ∉ (OPEN,UNSET) || continue
      v ∈ graph[vneig] || return false
    end
  end
  true
end

function compact!(p::Polyhedron)
  ids = findall(i->!isactive(p,i),1:num_vertices(p))
  old_to_new = fill(UNSET,num_vertices(p))
  new = 0
  for i in 1:num_vertices(p)
    isactive(p,i) || continue
    new += 1
    old = i
    old_to_new[old] = new
  end
  vertices = get_vertex_coordinates(p)
  graph = get_graph(p)
  deleteat!(vertices,ids)
  deleteat!(graph,ids)
  f = i -> i ∈ (OPEN,UNSET) ? i : old_to_new[i]
  map!(i->map!(f,i,i),graph,graph)
  compact!(p.data,ids,old_to_new)
  p
end

function compact!(data::PolyhedronData,ids,old_to_new)
  v_to_Π = data.vertex_to_planes
  v_to_of = data.vertex_to_original_faces
  v_to_pv = data.vertex_to_parent_vertex
  v_to_pe = data.vertex_to_parent_edge
  Π_to_v_to_dist = data.plane_to_vertex_to_distances
  Π_to_ref_Π = data.plane_to_ref_plane
  Π_to_id = data.plane_to_ids

  deleteat!(v_to_Π,ids)
  deleteat!(v_to_of,ids)
  deleteat!(v_to_pv,ids)
  for (i,pv) in enumerate(v_to_pv)
    if old_to_new[pv] == UNSET
      old_to_new[pv] = i
    end
    v_to_pv[i] = old_to_new[pv]
  end
  deleteat!(v_to_pe,ids)
  f = i -> i == UNSET ? i : old_to_new[i]
  for (i,pe) in enumerate(v_to_pe)
    v_to_pe[i] = (f(pe[1]),f(pe[2])) 
  end
  for v_to_dist in Π_to_v_to_dist
    deleteat!(v_to_dist,ids)
  end
  data
end

function Base.copy(data::PolyhedronData)
  v_to_Π = copy(data.vertex_to_planes)
  v_to_of = copy(data.vertex_to_original_faces)
  v_to_pv = copy(data.vertex_to_parent_vertex)
  v_to_pe = copy(data.vertex_to_parent_edge)
  Π_to_v_to_dist = _copy(data.plane_to_vertex_to_distances)
  Π_to_ref_Π = data.plane_to_ref_plane
  Π_to_id = data.plane_to_ids
  PolyhedronData( v_to_Π, v_to_of, v_to_pv, v_to_pe, Π_to_v_to_dist, Π_to_ref_Π, Π_to_id )
end

function compute_distances!(p::Polyhedron,Π,faces)
  data = get_data(p)
  v_to_f = get_data(p).vertex_to_original_faces
  Π_to_v_to_d = get_plane_distances(data)
  Π_ids = get_plane_ids(data)
  Π_to_rΠ = data.plane_to_ref_plane
  num_Π = length(Π_ids)
  append!(Π_ids,faces)
  append!(Π_to_rΠ,fill(UNSET,length(faces)))
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

function has_plane(data::PolyhedronData,Π)
  Π ∈ get_plane_ids(data)
end

## Polyhedra operations

function clip(poly::Polyhedron,Π;inside=true,inout=trues(length(Π)))
  p = poly
  for (i,Πi) in enumerate(Π)
    side = inside == inout[i] ? :left : :right
    p = clip(p,Πi,side)
    !isnothing(p) || break
  end
  p
end

function clip(p::Polyhedron,Π,side)
  @assert side ∈ (:right,:left)
  p⁻,p⁺ = split(p,Π,side)
  side == :left ? p⁻ : p⁺ 
end

function split(p::Polyhedron,Π,side=:both)
  @assert side ∈ (:both,:left,:right)
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
  if smax < 0
    if has_coplanars(get_data(p),Π)
      add_plane!(get_data(p),Π)
    end
    return p,nothing
  end
  new_vertices = empty( get_vertex_coordinates(p) )
  in_graph = _copy(get_graph(p))
  if side == :both
    out_graph = _copy(get_graph(p))
  end
  ≶ = side ∈ (:both,:left) ? (≥) : (<)
  data = copy( get_data(p) )
  D = num_dims(p)
  for v in 1:num_vertices(p)
    isactive(p,v) || continue
    distances[v] ≶ 0 && continue
    for (i,vneig) in enumerate( get_graph(p)[v] )
      vneig ∉ (UNSET,OPEN) || continue
      distances[vneig] ≶ 0 || continue
      vertex = compute_intersection(p[v],distances[v],p[vneig],distances[vneig])
      push!( new_vertices, vertex )
      add_vertex!(data,v,vneig,Π)
      push!( in_graph, fill(UNSET,D) )
      in_graph[v][i] = num_vertices(p) + length(new_vertices)
      in_graph[end][1] = v
      if side == :both
        push!( out_graph, fill(UNSET,D) )
        ineig = findfirst( isequal(v), get_graph(p)[vneig] )
        out_graph[vneig][ineig] = num_vertices(p) + length(new_vertices)
        out_graph[end][1] = vneig
      end
    end
  end
  if side == :both
    complete_graph!(out_graph,num_vertices(p))
    disconnect_graph!(out_graph,num_vertices(p),distances,(<))
    add_open_vertices!(out_graph,p)
    out_vertices = [p.vertices;new_vertices]
    out_data = copy(data)
    update_data!(out_data,out_graph,num_vertices(p))
    p_out = compact!(Polyhedron(out_vertices,out_graph,isopen(p),out_data))
  else
    p_out = nothing
  end
  complete_graph!(in_graph,num_vertices(p))
  disconnect_graph!(in_graph,num_vertices(p),distances,(≶))
  add_open_vertices!(in_graph,p)
  in_vertices = [p.vertices;new_vertices]
  in_data = data
  update_data!(in_data,in_graph,num_vertices(p))
  p_in = compact!(Polyhedron(in_vertices,in_graph,isopen(p),in_data))
  p⁻ = side ≠ :right ? p_in : p_out
  p⁺ = side ≠ :right ? p_out : p_in
  p⁻,p⁺
end

function decompose(surf::Polyhedron,cell::Polyhedron,rfaces,empty_facets,stl::STL)
  i = findfirst(i->has_original_reflex_face(surf,i), rfaces )
  if i === nothing
    R = [surf],[cell]
  else
    rf = rfaces[i]
    if !has_coplanars(surf.data,rf) || !contains_coplanars(surf.data,rf) 
      s⁻,s⁺ = split(surf,rf)
      k⁻,k⁺ = split(cell,rf)
      S = [s⁻,s⁺]
      K = [k⁻,k⁺]
      if any(isnothing,S) || any(!has_facets,S)
        S,K = split_reflex_face(S,K,surf,cell,stl,rf,empty_facets)
      end
      if any(isnothing,K)
        j = findfirst(!isnothing,K)
        S,K = [S[j]],[K[j]]
      end
    else
      S,K = [surf],[cell]
    end
    if length(rfaces) == i
      R = S,K
    else
      Sr,Kr = empty(S),empty(K)
      rfaces = view(rfaces,i+1:length(rfaces))
      for (s,k) in zip(S,K)
        Si,Ki = decompose(s,k,rfaces,empty_facets,stl)
        append!(Sr,Si)
        append!(Kr,Ki)
      end
      R = Sr,Kr
    end
  end
  R
end

function simplexify(poly::Polyhedron{3})
  !isopen(poly) || return simplexify_surface(poly)
  vstart = fill(UNSET,num_vertices(poly))
  stack = Int32[]
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
  istouch = map( i -> falses(length(i)), get_graph(poly) )
  vertex_coordinates = get_vertex_coordinates(poly)
  T = Vector{Int32}[]
  X = vertex_coordinates
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    for i in 1:length(get_graph(poly)[v])
      !istouch[v][i] || continue
      istouch[v][i] = true
      vcurrent = v
      vnext = get_graph(poly)[v][i]
      while vnext != v
        inext = findfirst( isequal(vcurrent), get_graph(poly)[vnext] )
        !isnothing(inext) || break
        inext = ( inext % length( get_graph(poly)[vnext] ) ) + 1
        istouch[vnext][inext] = true
        vcurrent = vnext
        vnext = get_graph(poly)[vnext][inext]
        @assert vcurrent ≠ vnext
        vcurrent ≠ vnext || break
        if v ∉ (vstart[v],vcurrent,vnext)
          k = [vstart[v],v,vcurrent,vnext]
          push!(T,k)
        end
      end
    end
  end
  T,X
end

function simplexify_surface(poly::Polyhedron{3})
  stack = Int32[]
  v_to_pv = get_data(poly).vertex_to_parent_vertex
  istouch = map( i -> falses(length(i)), get_graph(poly) )
  T = Vector{Int32}[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    for i in 1:length(get_graph(poly)[v])
      !istouch[v][i] || continue
      istouch[v][i] = true
      vcurrent = v
      vnext = get_graph(poly)[v][i]
      vnext ∉ (OPEN,UNSET) || continue
      while vnext != v
        inext = findfirst( isequal(vcurrent), get_graph(poly)[vnext] )
        inext = ( inext % length( get_graph(poly)[vnext] ) ) + 1
        istouch[vnext][inext] = true
        vcurrent = vnext
        vnext = get_graph(poly)[vnext][inext]
        vnext ∉ (OPEN,UNSET) || break
        if v ∉ (vcurrent,vnext)
          k = [v,vcurrent,vnext]
          push!(T,k)
        end
      end
    end
  end
  T,get_vertex_coordinates(poly)
end

function simplexify_boundary(poly::Polyhedron{3},stl::STL)
  D = 3
  stack = Int32[]
  v_to_pv = get_data(poly).vertex_to_parent_vertex
  v_to_Π = get_data(poly).vertex_to_planes
  facedims = get_facedims(stl)
  facet_list = Int32[]
  istouch = map( i -> falses(length(i)), get_graph(poly) )
  vertex_coordinates = get_vertex_coordinates(poly)
  T = Vector{Int32}[]
  X = empty(vertex_coordinates)
  f_to_stlf = Int32[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    for i in 1:length(get_graph(poly)[v])
      !istouch[v][i] || continue
      istouch[v][i] = true
      vcurrent = v
      vnext = get_graph(poly)[v][i]
      vnext ∉ (OPEN,UNSET) || continue
      empty!(facet_list)
      append!( facet_list, v_to_Π[v] )
      filter!( i -> i > 0, facet_list)
      filter!( i -> facedims[i] == num_dims(stl), facet_list)
      filter!( i -> i ∈ v_to_Π[vnext], facet_list )
      while vnext != v
        inext = findfirst( isequal(vcurrent), get_graph(poly)[vnext] )
        inext = ( inext % length( get_graph(poly)[vnext] ) ) + 1
        istouch[vnext][inext] = true
        vcurrent = vnext
        vnext = get_graph(poly)[vnext][inext]
        vnext ∉ (OPEN,UNSET) || break
        filter!( i -> i ∈ v_to_Π[vnext], facet_list )
        !isempty(facet_list) || break
        if v ∉ (vcurrent,vnext)
          kx = (v,vcurrent,vnext)
          k = (1:D) .+ length(X)
          for x in kx
            push!(X,vertex_coordinates[x])
          end
          push!(T,k)
          push!(f_to_stlf,first(facet_list))
        end
      end
    end
  end
  T,X,f_to_stlf
end

function simplexify(polys::AbstractVector{<:Polyhedron{Dp,Tp}}) where {Dp,Tp}
  T = Vector{Int32}[]
  X = Point{Dp,Tp}[]
  for poly in polys
    Ti,Xi = simplexify(poly)
    append!(T, map(i->i.+length(X),Ti) )
    append!(X,Xi)
  end
  T,X
end

function simplexify_boundary(
  polys::AbstractVector{<:Polyhedron{Dp,Tp}},
  stl::STL) where {Dp,Tp}

  T = Vector{Int32}[]
  X = Point{Dp,Tp}[]
  f_to_stlf = Int32[]
  for poly in polys
    Ti,Xi,f_to_f_i = simplexify_boundary(poly,stl)
    append!(T, map(i->i.+length(X),Ti) )
    append!(X,Xi)
    append!(f_to_stlf,f_to_f_i)
  end
  T,X,f_to_stlf
end

## Getters

num_dims(::Polyhedron{D}) where D = D

num_vertices(a::Polyhedron) = length(a.vertices)

get_vertex_coordinates(a::Polyhedron) = a.vertices

Base.getindex(a::Polyhedron,i::Integer) = a.vertices[i]

@inline get_graph(a::Polyhedron) = a.edge_vertex_graph

get_data(a::Polyhedron) = a.data

Base.isopen(a::Polyhedron) = a.isopen

get_vertex_to_planes(a::PolyhedronData) = a.vertex_to_planes

function isactive(p::Polyhedron,vertex::Integer)
  !isempty( get_graph(p)[vertex] )
end

function plot(p::Polyhedron,filename=nothing)
  vertices = Int32[]
  for i in 1:num_vertices(p)
    if isactive(p,i)
      push!(vertices,i)
    end
  end
  vertex_to_node = Dict( zip( vertices, 1:length(vertices) ) )
  g = empty( get_graph(p) )
  for (i,vneigs) in enumerate( get_graph(p) )
    isactive(p,i) || continue
    nodes = Int32[]
    for (j,vneig) in enumerate(vneigs)
      if vneig ∈ (UNSET,OPEN)
        push!(vertices,UNSET)
        n = length(vertices)
      else
        n = vertex_to_node[vneig]
      end
      push!(nodes,n)
    end
    push!(g,nodes)
  end
  v_to_pv = p.data.vertex_to_parent_vertex
  names = String[]
  for v in vertices
    if v == UNSET
      name = ""
    else
      name =  v_to_pv[v] == v ? "$v" : "$v($(v_to_pv[v]))"
    end
    push!(names,name)
  end
  kwargs = (fontsize=10,node_size=0,line=(:dot,0.5,1),thickness_scaling=0.5)
  GraphRecipes.graphplot(g,names=names,curves=false;kwargs...)
  isnothing(filename) || Plots.savefig(filename)
  Plots.plot!()
end

function writevtk(p::Polyhedron,filename;kwargs...)
  writevtk(edge_mesh(p),filename;kwargs...)
end

function writevtk(ps::Array{<:Polyhedron},filename)
  for (i,p) in enumerate(ps)
    writevtk(p,filename*"$i")
  end
end

function edge_mesh(p::Polyhedron)
  edge_to_vertices = Vector{Int32}[]
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

function get_original_reflex_faces(p::Polyhedron{D},stl::STL;empty=true) where D
  get_original_faces(p,stl,Val{D-2}();empty)
end

function get_original_facets(p::Polyhedron{D},stl::STL;empty=false) where D
  get_original_faces(p,stl,Val{D-1}();empty)
end

function get_original_faces(p::Polyhedron,stl::STL,::Val{d};empty) where d
  faces = collect_original_faces(p,stl,d)
  filter!(f->has_original_face(p,f,Val{d}();empty),faces)
  faces
end

function collect_original_faces(p::Polyhedron,stl::STL,d::Integer)
  v_to_f = get_data(p).vertex_to_original_faces
  facedims = get_facedims(stl)
  faces = Int32[]
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

function has_original_reflex_face(p::Polyhedron{D},face::Integer;empty=true) where D
  has_original_face(p,face,Val{D-2}();empty)
end

function has_original_facet(p::Polyhedron{D},face::Integer;empty=false) where D
  has_original_face(p,face,Val{D-1}();empty)
end

function has_original_face(p::Polyhedron,face::Integer,::Val{0};empty=false)
  v_to_f = get_data(p).vertex_to_original_faces
  for v in 1:num_vertices(p)
    if isactive(p,v) && face ∈ v_to_f[v]
      return true
    end
  end
  false
end

function has_original_face(p::Polyhedron,face::Integer,::Val{1};empty=true)
  v_to_f = get_data(p).vertex_to_original_faces
  v_to_v = get_data(p).vertex_to_parent_vertex
  for v in 1:num_vertices(p)
    if isactive(p,v) && face ∈ v_to_f[v]
      for vneig in get_graph(p)[v]
        vneig ∉ (OPEN,UNSET) || continue
        if face ∈ v_to_f[vneig] && ( empty || v_to_v[v] ≠ v_to_v[vneig] )
          return true
        end
      end
    end
  end
  false
end

function has_original_face(p::Polyhedron,face::Integer,::Val{2};empty=false)
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
            if v_to_v[vnext] ∉ (v_to_v[v],v_to_v[vcurrent]) || empty
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

function has_faces(p::Polyhedron,::Val{0})
  for v in 1:num_vertices(p)
    isactive(p,v) || continue
    return true
  end
  false
end

function has_faces(p::Polyhedron,::Val{1})
  for v in 1:num_vertices(p)
    isactive(p,v) || continue
    for vneig in get_graph(p)[v]
      @assert isactive(p,vneig)
      vneig ∉ (OPEN,UNSET) || continue
      v_to_v[v] ≠ v_to_v[vneig] || continue
      return true
    end
    @unreachable
  end
  false
end

function has_faces(p::Polyhedron,::Val{2})
  d = 2
  v_to_v = get_data(p).vertex_to_parent_vertex
  for v in 1:num_vertices(p)
    isactive(p,v) || continue
    for vneig in get_graph(p)[v]
      vneig ∉ (OPEN,UNSET) || continue
      num_v = 1
      vcurrent = v
      vnext = vneig
      while vnext ≠ v 
        if v_to_v[vnext] ∉ (v_to_v[v],v_to_v[vcurrent])
          num_v += 1
        end
        vcurrent,vnext = vnext,next_vertex(p,vcurrent,vnext)
        vnext ∉ (UNSET,OPEN) || break
        if vnext == v && num_v ≥ d+1
          return true
        end
      end
    end
  end
  false
end

has_vertices(p::Polyhedron) = has_faces(p,Val{0}())

has_edges(p::Polyhedron) = has_faces(p,Val{1}())

has_facets(p::Polyhedron{D}) where D = has_faces(p,Val{D-1}())


## Setup


function compute_cell_to_facets(grid::CartesianGrid,stl::STL)
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
  δ = 0.1
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
      if fast_intersection(f,_pmin,_pmax,p)
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
  v_to_Π = [ Int32[] for _ in 1:num_vertices ]
  v_to_of = [ Int32[] for _ in 1:num_vertices ]
  v_to_v = collect(1:num_vertices)
  v_to_e = fill((UNSET,UNSET),num_vertices)
  Π_to_v_to_d = Vector{Int32}[]
  Π_to_rΠ = Int32[]
  Π_to_id = Int32[]
  PolyhedronData( v_to_Π, v_to_of, v_to_v, v_to_e, Π_to_v_to_d, Π_to_rΠ, Π_to_id  )
end

function polyhedron_data(p::Polytope)
  v_to_Π = map( i-> Int32.(i), -get_faces(p,0,num_dims(p)-1) )
  v_to_of = [ Int32[] for _ in 1:num_vertices(p) ]
  v_to_v = Int32.(collect(1:num_vertices(p)))
  v_to_e = fill((Int32(UNSET),Int32(UNSET)),num_vertices(p))
  Π_to_v_to_d = Vector{Int32}[]
  Π_to_rΠ = Int32[]
  Π_to_id = Int32[]
  PolyhedronData( v_to_Π, v_to_of, v_to_v, v_to_e, Π_to_v_to_d, Π_to_rΠ, Π_to_id  )
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

function correct_graph!(graph,num_vertices::Integer)
  for v in num_vertices+1:length(graph)
    i = 1
    while i ≤ length(graph[v])
      if graph[v][i] == v
        deleteat!(graph[v],i)
      else 
        i += 1
      end
    end
    i = 1
    while i ≤ length(graph[v]) && length(graph[v]) > 1
      _i = i == length(graph[v]) ? 1 : i+1
      if graph[v][i] == graph[v][_i]
        deleteat!(graph[v],i)
      else 
        i += 1
      end
    end
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

function disconnect_graph!(edge_graph,num_vertices,distances,(≶)::Function)
  for i in 1:num_vertices
    if distances[i] ≶ 0
      edge_graph[i] = empty!( edge_graph[i] )
    end
  end
end

function update_data!(data,graph,num_vertices)
  v_to_pv = data.vertex_to_parent_vertex
  v_to_pe = data.vertex_to_parent_edge
  Π_to_v_to_d = get_plane_distances(data)
  for v in num_vertices+1:length(graph)
    v_to_pv[v] == v || continue
    for vnext in num_vertices+1:length(graph)
      v_to_pv[vnext] == vnext || continue
      if v_to_pe[v] == v_to_pe[vnext] && v_to_pe[v] ≠ (UNSET,UNSET)
        v_to_pv[vnext] = v
        for v_to_d in Π_to_v_to_d
          v_to_d[v] = v_to_d[v_to_pv[v]]
        end
      end
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

function _intersect(a::AbstractVector{T},b::AbstractVector{T}) where T
  n = count(in(a),b)
  r = Vector{T}(undef,n)
  n = 0
  for i in a
    if i ∈ b
      n += 1
      r[n] = i
    end
  end
  r
end

function _copy(a::AbstractVector{<:AbstractVector}) 
  [ copy(i) for i in a ]
end

function add_vertex!(data::PolyhedronData,v1::Integer,v2::Integer,Πid::Integer)
  v_to_Π = data.vertex_to_planes
  v_to_f = data.vertex_to_original_faces
  v_to_pv = data.vertex_to_parent_vertex
  v_to_pe = data.vertex_to_parent_edge
  Πs = _intersect( v_to_Π[v1], v_to_Π[v2] )
  fs = _intersect( v_to_f[v1], v_to_f[v2] )
  pe = (v_to_pv[v1],v_to_pv[v2])
  push!( Πs, Πid )
  push!( v_to_Π, Πs )
  push!( v_to_f, fs )
  push!( v_to_pe, pe)
  Π_to_v_to_d = get_plane_distances(data)
  Π_to_id = get_plane_ids(data)
  i = findfirst(isequal(Πid),Π_to_id)
  d1 = Π_to_v_to_d[i][v1] # Distances may be in the input
  d2 = Π_to_v_to_d[i][v2]
  pv = iszero(d1) ? v_to_pv[v1] : length(v_to_pv)+1 
  pv = iszero(d2) ? v_to_pv[v2] : pv 
  push!(v_to_pv,pv)
  @assert d1 ≠ d2
  ref_i = data.plane_to_ref_plane[i]
  for i in 1:length(Π_to_v_to_d)
    v_to_d = Π_to_v_to_d[i] 
    if pv ≠ length(v_to_pv)
      dist = v_to_d[pv]
    elseif ref_i == UNSET || data.plane_to_ref_plane[i] ∈ (UNSET,i)
      if Π_to_id[i] == Πid 
        dist = 0.0
      else
        dist = (d1*v_to_d[v2]-d2*v_to_d[v1])/(d1-d2)
      end
    else
      j = data.plane_to_ref_plane[i]
      @assert abs(j) < i
      dist = sign(j) * Π_to_v_to_d[abs(j)][end]
    end
    push!(v_to_d,dist)
  end
  data
end

function compute_graph(stl::STL{2})
  @notimplementedif get_polytope(stl) ≠ TRI
  D = 2
  v_to_e = get_faces(stl,0,1)
  v_to_f = get_faces(stl,0,D)
  f_to_v = get_faces(stl,D,0)
  ec = array_cache(v_to_e)
  fc = array_cache(v_to_f)
  vc = array_cache(f_to_v)
  graph = [ Int32[] for _ in 1:num_vertices(stl) ]
  for v in 1:num_vertices(stl)
    vfacets = getindex!(fc,v_to_f,v)
    !isempty(vfacets) || continue
    f0 = first(vfacets)
    fnext = f0
    while true
      i = findfirst(isequal(v), getindex!(vc,f_to_v,fnext) )
      inext = i == D+1 ? 1 : i+1
      vnext = getindex!(vc,f_to_v,fnext)[inext]
      push!(graph[v],vnext)
      faces = getindex!(fc,v_to_f,vnext)
      i = findfirst( f -> f ≠ fnext && v ∈ getindex!(vc,f_to_v,f), faces )
      fnext = isnothing(i) ? UNSET : faces[i]
      fnext ≠ UNSET || break
      fnext ≠ f0 || break
    end
    if fnext == UNSET
      n_edges = length( getindex!(ec,v_to_e,v) )
      if length(graph[v]) < n_edges
        v0 = first(graph[v])
        fnext = f0
        while fnext ≠ UNSET
          i = findfirst(isequal(v), getindex!(vc,f_to_v,fnext) )
          inext = i == 1 ? D+1 : i-1
          vnext = getindex!(vc,f_to_v,fnext)[inext]
          pushfirst!(graph[v],vnext)
          faces = getindex!(fc,v_to_f,vnext)
          i = findfirst( f -> f ≠ fnext && v ∈ getindex!(vc,f_to_v,f), faces )
          fnext = isnothing(i) ? UNSET : faces[i]
        end
        @assert length(graph[v]) == n_edges
      end
      push!(graph[v],OPEN)
    end
  end
  check_graph(graph) || error("Unable to build vertex-edge graph")
  graph
end

function set_original_faces!(data::PolyhedronData,stl::STL)
  Dc = num_dims(stl)
  v_to_f = data.vertex_to_original_faces
  for d in 0:Dc
    offset = get_offset(stl,d)
    _v_to_f = get_faces(stl,0,d)
    c = array_cache(_v_to_f)
    for v in 1:length(v_to_f)
      append!( v_to_f[v], map( i -> i.+offset, getindex!(c,_v_to_f,v) ) )
    end
  end
  data
end

function set_original_faces!(p::Polyhedron,a...)
  set_original_faces!(get_data(p),a...)
  p
end

function restrict(poly::Polyhedron,stl::STL,stl_facets)
  cell_to_nodes = get_cell_vertices(stl)
  c = array_cache(cell_to_nodes)
  nodes = Int32[]
  for f in stl_facets
    for n in getindex!(c,cell_to_nodes,f)
      if n ∉ nodes
        push!(nodes,n)
      end
    end
  end
  sort!(nodes)
  restrict(poly,nodes)
end

function restrict(p::Polyhedron,nodes)
  graph = get_graph(p)[nodes]
  f = i -> Int32( i ∈ nodes ? findfirst(isequal(i),nodes) : OPEN )
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
  v_to_e = fill((UNSET,UNSET),length(nodes))
  Π_to_rΠ = copy(data.plane_to_ref_plane)
  Π_to_id = copy(data.plane_to_ids)
  Π_to_v_to_d = [ data.plane_to_vertex_to_distances[i][nodes] for i in 1:length(Π_to_id) ]
  PolyhedronData( v_to_Π, v_to_of, v_to_v, v_to_e, Π_to_v_to_d, Π_to_rΠ, Π_to_id  )
end

function next_vertex(p::Polyhedron,vprevious::Integer,vcurrent::Integer)
  i = findfirst( isequal(vprevious), get_graph(p)[vcurrent] )
  i = ( i % length( get_graph(p)[vcurrent] ) ) + 1
  get_graph(p)[vcurrent][ i ]
end

function get_cell_planes(p::Polytope,pmin::Point,pmax::Point)
  @notimplementedif !is_n_cube(p)
  D = num_dims(p)
  [ CartesianPlane(isodd(i)*pmin+iseven(i)*pmax,D-((i-1)÷2),1) for i in 1:2*D ],
  -(1:2*D),iseven.(1:2*D)
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
  bisector_plane(e,Plane(f1),Plane(f2))
end

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

function bisector_plane(edge::Face{1,3},Π1::Plane,Π2::Plane)
  n1 = normal(Π1)
  n2 = normal(Π2)
  n1 ⋅ n2 ≉ -1 || error("Edge too sharp")
  n = n1-n2
  if norm(n) < 1
    v = edge[2]-edge[1]
    v /= norm(v)
    _n = n1+n2
    @assert norm(_n) > 1
    _n /= norm(_n)
    n = _n × v
    @assert norm(n) ≈ 1
  end
  n /= norm(n)
  @assert norm(n) ≈ 1
  c = center(edge)
  Plane(c,n)
end

#function get_reflex_planes(stl::DiscreteModel) 
#  Dc = num_dims(stl)
#  Dp = num_point_dims(stl)
#  T_Π = typeof( bisector_plane(stl,Dc-1,1) )
#  Π_r = Vector{T_Π}(undef,num_faces(stl,Dc-1))
#  for f in 1:num_faces(stl,Dc-1)
#    Π_r[f] = bisector_plane(stl,Dc-1,f)
#  end
#  Π_r
#end

function get_reflex_planes(stl::STL,Πf) 
  Dc = num_dims(stl)
  Dp = num_point_dims(stl)
  c = bisector_plane_cache(stl,Dc-1)
  T_Π = typeof( bisector_plane!(c,stl,Dc-1,1,Πf) )
  Π_r = Vector{T_Π}(undef,num_faces(stl,Dc-1))
  for f in 1:num_faces(stl,Dc-1)
    Π_r[f] = bisector_plane!(c,stl,Dc-1,f,Πf)
  end
  Π_r
end

function get_facet_plane!(c,stl::STL,i)
  f = get_cell!(c,stl,i)
  Plane(f)
end

function get_facet_planes(stl::STL)
  c = get_cell_cache(stl)
  [ get_facet_plane!(c,stl,i) for i in 1:num_cells(stl) ]
end

function group_facing_facets(poly::Polyhedron,facets,part_to_facets;inside)
  length(part_to_facets) > 1 || return [facets]
  f_to_part = fill(UNSET,length(facets))
  for (i,f) in enumerate(facets)
    f_to_part[i] = findfirst( p -> f ∈ p, part_to_facets )
  end
  parts = unique(f_to_part)
  length(parts) > 1 || return [facets]
  for i in 1:length(facets)
    p = findfirst(isequal(f_to_part[i]),parts)
    f_to_part[i] = p
  end
  p_to_p_to_facing = [ trues(length(parts)) for _ in 1:length(parts) ]
  for (i,fi) in enumerate(facets), (j,fj) in enumerate(facets)
    fi ≠ fj || continue
    f_to_part[i] ≠ f_to_part[j] || continue 
    if !is_facet_in_facet(poly,fj,fi;inside)
      p_to_p_to_facing[f_to_part[i]][f_to_part[j]] = false
    end
  end
  p_to_group = fill(UNSET,length(parts))
  group = Int32[]
  ids = Int32[]
  num_groups = 0
  for p in 1:length(parts)
    p_to_group[p] == UNSET || continue
    num_groups += 1
    empty!(group)
    push!(group,p)
    for i in 1:length(parts)
      if p_to_p_to_facing[p][i] && p_to_p_to_facing[i][p] && i ∉ group
        push!(group,i)
      end
    end
    empty!(ids)
    for (i,p_i) in enumerate(group)
      for p_j in group
        if !p_to_p_to_facing[p_j][p_i]
          push!(ids,i)
          break
        end
      end
    end
    deleteat!(group,ids)
    for p_i in group
      @assert p_to_group[p_i] == UNSET
      p_to_group[p_i] = num_groups
    end
  end
  group_to_facets = [ Int32[] for _ in 1:num_groups ]
  for (i,f) in enumerate(facets)
    g = p_to_group[f_to_part[i]]
    push!(group_to_facets[g],f)
  end
  group_to_facets
end


function get_disconnected_facets(poly::Polyhedron,stl::STL)
  v_to_pv = poly.data.vertex_to_parent_vertex
  v_to_f = poly.data.vertex_to_original_faces
  facets = get_original_facets(poly,stl,empty=true)
  facet_to_part = fill(UNSET,length(facets))
  stack = Int32[]
  num_parts = 0
  for (i,facet) in enumerate(facets)
    facet_to_part[i] == UNSET || continue
    num_parts += 1
    facet_to_part[i] = num_parts
    empty!(stack)
    push!(stack,facet)
    while !isempty(stack)
      current_facet = pop!(stack)
      v,vnext = get_facet_onset(poly,current_facet) # reduce calls
      vcurrent = v
      while true
        if v_to_pv[vcurrent] ≠ v_to_pv[vnext]
          for fi in v_to_f[vcurrent]
            fi ≠ current_facet || continue
            for fj in v_to_f[vnext]
              fj ≠ current_facet || continue
              if fi == fj
                fi ∈ facets || continue
                _i = findfirst(isequal(fi),facets)
                facet_to_part[_i] == UNSET || continue
                facet_to_part[_i] = num_parts
                push!(stack,fi)
              end
            end
          end
        end
        vcurrent,vnext = vnext,next_vertex(poly,vcurrent,vnext)
        vcurrent ≠ v || break
      end
    end
  end
  part_to_facets = [ Int32[] for _ in 1:num_parts ]
  for (i,facet) in enumerate(facets)
    push!(part_to_facets[facet_to_part[i]],facet)
  end
  part_to_facets
end

function get_facet_onset(poly::Polyhedron,facet::Integer)
  v_to_pv = poly.data.vertex_to_parent_vertex
  v_to_f = poly.data.vertex_to_original_faces
  for v in 1:num_vertices(poly)
    if isactive(poly,v) && facet ∈ v_to_f[v]
      for vneig in get_graph(poly)[v]
        vneig ∉ (OPEN,UNSET) || continue
        if facet ∈ v_to_f[vneig]
          vcurrent = v
          vnext = vneig
          while vnext ≠ v
            vcurrent,vnext = vnext,next_vertex(poly,vcurrent,vnext)
            vnext ∉ (UNSET,OPEN) || break
            facet ∈ v_to_f[vnext] || break
          end
          if vnext == v
            return v,vneig
          end
        end
      end
    end
  end
end

function is_facet_in_facet(poly::Polyhedron,facet,plane;inside,atol=0)
  if has_coplanars(poly.data,plane)
    if are_coplanar(poly.data,facet,plane)
      return true
    end
  end
  v_to_f = get_data(poly).vertex_to_original_faces
  distances = get_plane_distances(get_data(poly),plane)
  smax = -Inf
  smin = Inf
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    if facet ∈ v_to_f[v] && plane ∉ v_to_f[v]
      smin = min(smin,distances[v])
      smax = max(smax,distances[v])
    end
  end
 # @assert !(smin == smax == 0)
  ( smin ≥ 0 && smax ≥ atol ) && return !inside
  ( smin ≤ 0 && smax ≤ -atol ) && return inside
  false
end

function is_reflex(poly::Polyhedron,stl::STL,reflex_face;inside)
  Dc = num_dims(stl)
  rf = reflex_face - get_offset(stl,Dc-1)
  length(get_faces(stl,Dc-1,Dc)[rf]) == 2 || return false
  f1 = get_faces(stl,Dc-1,Dc)[rf][1] + get_offset(stl,Dc)
  f2 = get_faces(stl,Dc-1,Dc)[rf][2] + get_offset(stl,Dc)
  has_plane(poly.data,f1) || return false
  has_plane(poly.data,f2) || return false
  is_facet_in_facet(poly,f1,f2;inside) || return true
  is_facet_in_facet(poly,f2,f1;inside) || return true
  false
end

function refine(
  K::Polyhedron,
  Γ::Polyhedron,
  stl::STL,
  reflex_faces::AbstractVector,
  empty_facets::AbstractVector=[],
  ;inside::Bool)

  reflex_faces = filter( f -> is_reflex(Γ,stl,f;inside), reflex_faces )
  Γn,Kn = decompose(Γ,K,reflex_faces,empty_facets,stl)
  Kn_clip = empty(Kn)
  for (i,(Γi,Ki)) in enumerate(zip(Γn,Kn))
    facets = get_original_facets(Γi,stl)
    if length(facets) > 1
      ids = findall(in(empty_facets),facets)
      if length(ids) < length(facets)
        deleteat!(facets,ids)
      end
    end
    facets = add_missing_facets(Γi,stl,facets,reflex_faces,empty_facets)
    @assert !isempty(facets)
    part_to_facets = get_disconnected_facets(Γi,stl)
    group_to_facets = group_facing_facets(Γi,facets,part_to_facets;inside)
    for facets in group_to_facets
      Ki_clip = clip(Ki,facets;inside)
      !isnothing(Ki_clip) || continue
      push!(Kn_clip,Ki_clip)
    end
  end
  Kn_clip
end

function filter_face_planes(
  stl::STL,
  reflex_planes::AbstractVector,
  reflex_faces::AbstractVector{<:Integer},
  facet_planes::AbstractVector,
  facets::AbstractVector{<:Integer})
  
  Dc = num_dims(stl)
  Πr = view(reflex_planes,lazy_map(i->i-get_offset(stl,Dc-1),reflex_faces))
  Πf = view(facet_planes,lazy_map(i->i-get_offset(stl,Dc),facets))
  lazy_append(Πr,Πf),lazy_append(reflex_faces,facets)
end


function get_cell_nodes_to_inout(polys_in,polys_out,p::Polytope)
  node_to_inout = fill(UNSET,num_vertices(p))
  
  complete_nodes_to_inout!(node_to_inout,polys_in,FACE_IN,p)
  complete_nodes_to_inout!(node_to_inout,polys_out,FACE_OUT,p)
  #@assert UNSET ∉ node_to_inout
  node_to_inout
end

function get_cell_facets_to_inoutcut(polys_in,polys_out,p::Polytope)
  facet_to_inoutcut = fill(UNSET,num_facets(p))
  
  complete_facets_to_inoutcut!(facet_to_inoutcut,polys_in,FACE_IN,p)
  complete_facets_to_inoutcut!(facet_to_inoutcut,polys_out,FACE_OUT,p)
  @assert UNSET ∉ facet_to_inoutcut
  facet_to_inoutcut
end

function complete_nodes_to_inout!(node_to_inout,polys,inout,p::Polytope)
  D = num_dims(p)
  for poly in polys
    v_to_Π = poly.data.vertex_to_planes
    v_to_v = poly.data.vertex_to_parent_vertex
    for v in 1:num_vertices(poly)
      isactive(poly,v) || continue
    #  v = inout == FACE_CUT ? v : v_to_v[v]
      if count(Π -> -2*D ≤ Π < 0 ,v_to_Π[v]) == D
        i = 0
        node = 0
        for _ in 1:D
          i = findnext(Π -> Π < 0,v_to_Π[v],i+1)
          f = -v_to_Π[v][i]
          d = D - ((f-1)>>1)
          ud = iseven(f)
          node |= ud<<(d-1)
        end
        node += 1
        #@assert  node_to_inout[node] ∈ (inout,UNSET)
        _inout = inout
        if node_to_inout[node] ∉ (inout,UNSET)
          _inout = FACE_CUT
        end
        node_to_inout[node] = _inout
      end
    end
  end
  node_to_inout
end

function complete_facets_to_inoutcut!(facet_to_inoutcut,polys,inout,p::Polytope)
  facet_list = Int32[]
  for poly in polys
    istouch = map( i -> falses(length(i)), get_graph(poly) )
    v_to_Π = poly.data.vertex_to_planes
    for v in 1:num_vertices(poly)
      isactive(poly,v) || continue
      for i in 1:length(get_graph(poly)[v])
        !istouch[v][i] || continue
        istouch[v][i] = true
        vcurrent = v
        vnext = get_graph(poly)[v][i]
        vnext ∉ (OPEN,UNSET) || continue
        any( i -> i < 0, v_to_Π[v] ) || continue
        empty!(facet_list)
        append!( facet_list, v_to_Π[v] )
        filter!( i -> i < 0, facet_list)
        filter!( i -> i ∈ v_to_Π[vnext], facet_list )
        while vnext != v
          inext = findfirst( isequal(vcurrent), get_graph(poly)[vnext] )
          inext = ( inext % length( get_graph(poly)[vnext] ) ) + 1
          istouch[vnext][inext] = true
          vcurrent = vnext
          vnext = get_graph(poly)[vnext][inext]
          vnext ∉ (OPEN,UNSET) || break
          filter!( i -> i ∈ v_to_Π[vnext], facet_list )
          !isempty(facet_list) || break
          if vnext == v
            f = abs( only( facet_list) )
            if facet_to_inoutcut[f] == UNSET
              facet_to_inoutcut[f] = inout
            elseif facet_to_inoutcut[f] ≠ inout
              facet_to_inoutcut[f] = FACE_CUT
            end
          end
        end
      end
    end
  end
  facet_to_inoutcut
end

function Base.eps(T::Type{<:AbstractFloat},grid::Grid)
  pmin,pmax = get_bounding_box(grid)
  vmax = max(abs.(Tuple(pmin))...,abs.(Tuple(pmax))...)
  eps(T(vmax))
end

Base.eps(grid::Grid) = eps(Float64,grid)

Base.eps(model::DiscreteModel) = eps(get_grid(model))

function compute_submesh(
  bgmodel::CartesianDiscreteModel,
  stlmodel::DiscreteModel;
  threading=:spawn,
  kdtree=false,
  tolfactor=1e3)

  grid = get_grid(bgmodel)
  grid_topology = get_grid_topology(bgmodel)
  p = get_polytope(only(get_reffes(bgmodel)))
  node_to_coords = get_node_coordinates(grid)
  cell_to_nodes = get_cell_node_ids(grid)
  stl = STL(stlmodel)
  D = num_dims(grid)
  atol = eps(grid)*tolfactor

  f_to_isempty = get_facet_to_isempty(stl;atol)
  Πf = get_facet_planes(stl)
  Πf = correct_small_facets_planes!(stl,Πf,f_to_isempty;atol)
  Πr = get_reflex_planes(stl,Πf)

  c_to_stlf = compute_cell_to_facets(grid,stl)

  Γ0 = Polyhedron(stl)

  submesh = _get_threaded_empty_arrays(stl)
  io_arrays = _init_io_arrays(bgmodel,p)
  caches = _get_threaded_caches(cell_to_nodes)

  cut_cells = filter(i->!isempty(c_to_stlf[i]),1:num_cells(grid))
  if threading == :threads
    Threads.@threads for cell in cut_cells
      save_cell_submesh!(submesh,io_arrays,stl,p,cell,
        compute_polyhedra!(caches,Γ0,stl,p,f_to_isempty,Πf,Πr,
          c_to_stlf,node_to_coords,cell_to_nodes,cell;atol,kdtree)... )
    end
  elseif threading == :spawn
    @sync for cell in cut_cells
      Threads.@spawn save_cell_submesh!(submesh,io_arrays,stl,p,cell,
        compute_polyhedra!(caches,Γ0,stl,p,f_to_isempty,Πf,Πr,
          c_to_stlf,node_to_coords,cell_to_nodes,cell;atol,kdtree)... )
    end
  else
    @unreachable
  end

  submesh = _append_threaded_submesh!(submesh)
  io_arrays = _reduce_io_arrays(bgmodel,io_arrays)
  bgcell_to_ioc, bgnode_to_io, bgfacet_to_ioc = io_arrays
  T,F,X,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf = submesh

  propagate_inout!(bgmodel,bgcell_to_ioc,bgnode_to_io) 
  set_facets_as_inout!(bgmodel,bgcell_to_ioc,bgfacet_to_ioc)
  
  delete_small_subcells!(bgmodel,T,X,k_to_io,k_to_bgcell)
  delete_small_subfacets!(bgmodel,F,Xf,f_to_bgcell,f_to_stlf)
  T,X,F,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf,bgcell_to_ioc,bgfacet_to_ioc
end

function _get_empty_arrays(stl::STL)
  P = eltype( get_vertex_coordinates(stl) )
  T = Vector{Int32}[]
  F = Vector{Int32}[]
  X = P[]
  Xf = P[] 
  k_to_io = Int8[]
  k_to_bgcell = Int32[]
  f_to_bgcell = Int32[]
  f_to_stlf = Int32[]
  T,F,X,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf
end

function _get_threaded_empty_arrays(stl::STL)
  n = Threads.nthreads()
  [ _get_empty_arrays(stl::STL) for _ in 1:n ]
end

function _threaded_empty_array(::Type{T}) where T
  n = Threads.nthreads()
  [ T[] for i in 1:n ]
end

function _init_io_arrays(model::DiscreteModel,p::Polytope)
  D = num_dims(model)
  ncells = num_cells(model)
  nnodes_x_cell = num_vertices(p)
  nfacets_x_cell = num_facets(p)

  bgcell_to_ioc = fill(Int8(UNSET),ncells)
  bgcell_node_to_io = fill(Int8(UNSET),nnodes_x_cell,ncells)
  bgcell_facet_to_ioc = fill(Int8(UNSET),nfacets_x_cell,ncells)
  bgcell_to_ioc, bgcell_node_to_io, bgcell_facet_to_ioc
end

function _get_caches(c_to_n)
  nc = array_cache(c_to_n)
  facets = Int32[]
  expanded_facets = Int32[]
  empty_facets = Int32[]
  reflex_faces = Int32[]
  nc,facets,expanded_facets,empty_facets,reflex_faces
end

function _get_threaded_caches(c_to_n)
  n = Threads.nthreads()
  [ _get_caches(c_to_n) for _ in 1:n ]
end

function compute_polyhedra!(caches,Γ0,stl,p,f_to_isempty,Πf,Πr,
  c_to_stlf,node_to_coords,cell_to_nodes,cell;atol,kdtree)

  cell_coords = _get_cell_coordinates!(caches,node_to_coords,cell_to_nodes,cell)
  pmin = cell_coords[1] - atol
  pmax = cell_coords[end] + atol
   _faces = _get_cell_faces!(caches,stl,c_to_stlf,cell,f_to_isempty)
  facets,expanded_facets,empty_facets,reflex_faces = _faces

  Πk,Πk_ids,Πk_io = get_cell_planes(p,pmin,pmax)
  Πkf,Πkf_ids = filter_face_planes(stl,Πr,reflex_faces,Πf,facets) 

  K = Polyhedron(p,cell_coords)
  Γk0 = restrict(Γ0,stl,expanded_facets)

  compute_distances!(Γk0,lazy_append(Πk,Πkf),lazy_append(Πk_ids,Πkf_ids))
  compute_distances!(K,Πkf,Πkf_ids)
  Π_to_refΠ,Πs = link_planes!(Γk0,stl,Πr,Πf;atol)
  set_linked_planes!(Γk0,Π_to_refΠ,Πs)
  set_linked_planes!(K,Π_to_refΠ,Πs)

  Γk = clip(Γk0,Πk_ids,inout=Πk_io)
  if isnothing(Γk) || isempty(get_original_facets(Γk,stl))
    nothing,nothing
  else
    if kdtree
      Γv,Kv = refine_by_vertices(Γk,K,atol)
    else
      Γv,Kv = [Γk],[K]
    end
    Kn_in = typeof(K)[]
    Kn_out = typeof(K)[]
    for (γk,k) in zip(Γv,Kv)
      k_in = refine(k,γk,stl,reflex_faces,empty_facets,inside=true)
      k_out = refine(k,γk,stl,reflex_faces,empty_facets,inside=false)
      append!(Kn_in,k_in)
      append!(Kn_out,k_out)
    end
    Kn_in,Kn_out
  end
end

function _get_cell_coordinates!(caches,node_to_coords,cell_to_nodes,cell)
  i = Threads.threadid()
  nc, = caches[i]
  cell_nodes = getindex!(nc,cell_to_nodes,cell)
  cell_coords = view(node_to_coords,cell_nodes)
  map(i->Point(float.(Tuple(i))),cell_coords)
end

function _get_cell_faces!(caches,stl,c_to_stlf,cell,f_to_isempty)
  i = Threads.threadid()
  _,facets,expanded_facets,empty_facets,reflex_faces = caches[i]
  D = num_dims(stl)
  get_expanded_facets!(expanded_facets,stl,c_to_stlf,cell)
  get_reflex_faces!(reflex_faces,stl,c_to_stlf,cell)
  copy!(facets,c_to_stlf[cell])
  copy!(empty_facets,c_to_stlf[cell])
  filter!(f->f_to_isempty[f],empty_facets)
  map!(i->i+get_offset(stl,D),facets,facets)
  map!(i->i+get_offset(stl,D),empty_facets,empty_facets)
  map!(i->i+get_offset(stl,D-1),reflex_faces,reflex_faces)
  facets,expanded_facets,empty_facets,reflex_faces
end

function _get_cell_io(T_Γ,Kn_in,Kn_out)
  if length(T_Γ) > 0
    FACE_CUT
  else
    if length(Kn_in) > 0 && length(Kn_out) == 0
      FACE_IN
    elseif length(Kn_in) == 0 && length(Kn_out) > 0
      FACE_OUT
    else
      UNSET
    end
  end
end

function _append_submesh!(submesh,Xin,Tin,Xout,Tout,Xfn,Tfn,fn_to_f,bgcell)
  i = Threads.threadid()
  _submesh = submesh[i]
  _append_subcells!(_submesh,Xin,Tin,FACE_IN,bgcell)
  _append_subcells!(_submesh,Xout,Tout,FACE_OUT,bgcell)
  _append_subfacets!(_submesh,Xfn,Tfn,fn_to_f,bgcell)
end

function _append_subcells!(submesh_arrays,Xn,Tn,io,bgcell)
  T,F,X,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf = submesh_arrays
  append!(T, map(i->i.+length(X),Tn) )
  append!(X,Xn)
  append!(k_to_io,fill(io,length(Tn)))
  append!(k_to_bgcell,fill(bgcell,length(Tn)))
  submesh_arrays
end

function _append_subfacets!(submesh_arrays,Xfn,Tfn,fn_to_f,bgcell)
  T,F,X,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf = submesh_arrays
  append!(F, map(i->i.+length(Xf),Tfn) ) 
  append!(Xf,Xfn)
  append!(f_to_bgcell,fill(bgcell,length(Tfn)))
  append!(f_to_stlf,fn_to_f)
  submesh_arrays
end

function _append_threaded_submesh!(submesh)
  n = Threads.nthreads()
  for i in 2:n
    _append!(submesh[1],submesh[i])
  end
  submesh[1]
end

function _append!(submesh_a,submesh_b)
  T,F,X,Xf, = submesh_a
  T,F, = submesh_b
  nnodes_cells_a = length(X)
  nnodes_facets_a = length(Xf)
  map!(i->map!(j->j+nnodes_cells_a,i,i),T,T)
  map!(i->map!(j->j+nnodes_facets_a,i,i),F,F)
  for (a,b) in zip(submesh_a,submesh_b)
    append!(a,b)
  end
  submesh_a
end

function _reduce_io_arrays(model::DiscreteModel,io_arrays)
  bgcell_to_io, bgcell_node_to_io, bgcell_facet_to_io = io_arrays
  D = num_dims(model)
  topo = get_grid_topology(model)
  bgnode_to_io = fill(Int8(UNSET),num_vertices(topo) )
  bgfacet_to_io = fill(Int8(UNSET),num_facets(topo) )
  c_to_n = get_cell_vertices(topo)
  c_to_f = get_faces(topo,D,D-1)
  nc = array_cache(c_to_n)
  nf = array_cache(c_to_f)
  for cell in 1:num_cells(model)
    for (i,bg_n) in enumerate( getindex!(nc,c_to_n,cell) )
      io = bgcell_node_to_io[i,cell]
      io ≠ UNSET || continue
      bgnode_to_io[bg_n] ≠ FACE_CUT || continue
      if bgnode_to_io[bg_n] ∉ (UNSET,io)
        io = FACE_CUT
      end
      bgnode_to_io[bg_n] = io
    end
    for (i,bg_f) in enumerate( getindex!(nf,c_to_f,cell) )
      io = bgcell_facet_to_io[i,cell]
      if bgfacet_to_io[bg_f] ∉ (UNSET,io)
        io = FACE_CUT
      end
      bgfacet_to_io[bg_f] = io
    end
  end
  bgcell_to_io, bgnode_to_io, bgfacet_to_io
end

function save_cell_submesh!(submesh,io_arrays,stl,p,cell,Kn_in,Kn_out)
  !isnothing(Kn_in) || return
  bgcell_to_ioc, bgcell_node_to_io, bgcell_facet_to_ioc = io_arrays
  Tin,Xin = simplexify(Kn_in)
  Tout,Xout = simplexify(Kn_out)
  T_Γ,X_Γ,f_to_f = simplexify_boundary(Kn_in,stl)
  bgcell_to_ioc[cell] = _get_cell_io(T_Γ,Kn_in,Kn_out)
  bgcell_to_ioc[cell] == FACE_CUT || return
  D = num_dims(stl)
  f_to_f .-= get_offset(stl,D)
  n_to_io = get_cell_nodes_to_inout(Kn_in,Kn_out,p)
  f_to_ioc = get_cell_facets_to_inoutcut(Kn_in,Kn_out,p)
  _append_submesh!(submesh,Xin,Tin,Xout,Tout,X_Γ,T_Γ,f_to_f,cell)
  for i in 1:num_vertices(p)
    bgcell_node_to_io[i,cell] = n_to_io[i]
  end
  for i in 1:num_facets(p)
    bgcell_facet_to_ioc[i,cell] = f_to_ioc[i]
  end
end

function propagate_inout!(bgmodel,bgcell_to_ioc,bgnode_to_io)
  D = num_dims(bgmodel)
  grid_topology = get_grid_topology(bgmodel)
  stack = Int32[]
  n_to_c = get_faces(grid_topology,0,D)
  c_to_n = get_faces(grid_topology,D,0)
  node_cache = array_cache(c_to_n)
  neig_node_cache = array_cache(c_to_n)
  neig_cell_cache = array_cache(n_to_c)
  for cell in 1:num_cells(grid_topology)
    if bgcell_to_ioc[cell] ∈ (FACE_CUT,FACE_IN)
      resize!(stack,0)
      push!(stack,cell)
      while !isempty(stack)
        current_cell = pop!(stack)
        for node in getindex!(node_cache,c_to_n,cell)
          if bgnode_to_io[node] == FACE_IN || bgcell_to_ioc[cell] == FACE_IN
            for neig_cell in getindex!(neig_cell_cache,n_to_c,node)
              if bgcell_to_ioc[neig_cell] == UNSET
                bgcell_to_ioc[neig_cell] = FACE_IN
                for neig_node in getindex!(neig_node_cache,c_to_n,neig_cell)
                  if bgnode_to_io[neig_node] == UNSET
                    bgnode_to_io[neig_node] = FACE_IN
                  end
                  push!(stack,neig_cell)
                end
              end
            end
          end
        end
      end
    end
  end
  replace!(bgcell_to_ioc, UNSET => FACE_OUT )
  bgcell_to_ioc
end

function set_facets_as_inout!(bgmodel,bgcell_to_ioc,bgfacet_to_ioc)
  D = num_dims(bgmodel)
  grid_topology = get_grid_topology(bgmodel)
  c_to_f = get_faces(grid_topology,D,D-1)
  cache = array_cache(c_to_f)
  for cell in 1:num_cells(grid_topology)
    bgcell_to_ioc[cell] ≠ FACE_CUT || continue
    for facet in getindex!(cache,c_to_f,cell)
      bgfacet_to_ioc[facet] = bgcell_to_ioc[cell]
    end
  end
  bgfacet_to_ioc
end



function get_expanded_facets!(facets,stl,c_to_stlf,cell)
  D = num_dims(stl)
  f_to_v = get_faces(stl,D,0)
  v_to_f = get_faces(stl,0,D)
  nv = array_cache(f_to_v)
  nf = array_cache(v_to_f)
  empty!(facets)
  for f in c_to_stlf[cell] 
    for v in getindex!(nv,f_to_v,f)
      for _f in getindex!(nf,v_to_f,v)
        if _f ∉ facets
          push!(facets,_f)
        end
      end
    end
  end
  facets
end

function get_reflex_faces!(rfaces,stl,c_to_stlf,cell)
  D = num_dims(stl)
  f_to_r = get_faces(stl,D,D-1)
  r_to_f = get_faces(stl,D-1,D)
  rc = array_cache(f_to_r)
  fc = array_cache(r_to_f)
  empty!(rfaces)
  for f in c_to_stlf[cell]
    for r in getindex!(rc,f_to_r,f)
      facets = getindex!(fc,r_to_f,r)
      r += get_offset(stl,D-2)
      if r ∉ rfaces && all( i -> i ∈ c_to_stlf[cell], facets )
        push!(rfaces,r)
      end
    end
  end
  rfaces
end

## Kd Tree stuff

function delete_small_subcells!(bgmodel,T,X,arrays...)
  delete_small_subfaces!(bgmodel,T,X,TET,arrays...)
end

function delete_small_subfacets!(bgmodel,T,X,arrays...)
  delete_small_subfaces!(bgmodel,T,X,TRI,arrays...)
end

function delete_small_subfaces!(bgmodel,T,X,p::Polytope{D},arrays...) where D
  h = float(get_cartesian_descriptor(bgmodel).sizes[1])
  c = array_cache(T)
  ids = findall( i -> measure(get_cell!(c,T,X,p,i)) < eps(h^D), 1:length(T) )
  deleteat!(T,ids)
  for array in arrays
    deleteat!(array,ids)
  end
end

function refine_by_vertices(Γ::Polyhedron,K::Polyhedron,atol=0)
  refine_by_vertices(Γ,K,get_vertex_coordinates(Γ),atol)
end

function refine_by_vertices(Γ,K,vertices,atol)
  i = findfirst( v -> has_intersection(K,v), vertices )
  if isnothing(i)
    R = [Γ],[K]
  else
    v = vertices[i]
    d = get_direction(v,K)
    Π,Π⁻,Π⁺ = get_planes(v,d,atol)
    Πid,Πid⁻,Πid⁺ = get_new_plane_ids(Γ)
    compute_distances!(Γ,[Π⁻,Π⁺],[Πid⁻,Πid⁺])
    compute_distances!(K,[Π],[Πid])
    _,γ⁻ = split(Γ,Πid⁺)
    _,γ⁺ = split(Γ,Πid⁻)
    k⁻,k⁺ = split(K,Πid)
    Γn = [γ⁻,γ⁺]
    Kn = [k⁻,k⁺]
    if length(vertices) == i
      R = Γn,Kn
    else
      Γr,Kr = empty(Γn),empty(Kn)
      vertices = view(vertices,i+1:length(vertices))
      for (γ,k) in zip(Γn,Kn)
        Γi,Ki = refine_by_vertices(γ,k,vertices,atol)
        append!(Γr,Γi)
        append!(Kr,Ki)
      end
      R = Γr,Kr
    end
  end
  R
end

function get_bounding_box(cell::Polyhedron)
  i = findfirst(i->isactive(cell,i),1:num_vertices(cell))
  pmin,pmax = cell[i],cell[i]
  for i in 1:num_vertices(cell)
    if isactive(cell,i)
      pmin = Point( min.(Tuple(cell[i]),Tuple(pmin)) )
      pmax = Point( max.(Tuple(cell[i]),Tuple(pmax)) )
    end
  end
  pmin,pmax
end

function has_intersection(cell::Polyhedron,v::Point)
  pmin,pmax = get_bounding_box(cell)
  all( Tuple(pmin) .< Tuple(v) ) || return false
  all( Tuple(pmax) .> Tuple(v) ) || return false
  true
end

function get_direction(v,cell)
  pmin,pmax = get_bounding_box(cell)
  d0 = Tuple(v) .- Tuple(pmin)
  d1 = Tuple(pmax) .- Tuple(v)
  d = min.(d0,d1)
  d_max = maximum(d)
  findfirst(isequal(d_max),d)
end

function get_planes(v,d,atol)
  δ = VectorValue(Base.setindex(Tuple(zero(v)),atol,d))
  Π = CartesianPlane(v,d,+1)
  Π⁻ = CartesianPlane(v-δ,d,+1)
  Π⁺ = CartesianPlane(v+δ,d,-1)
  Π,Π⁻,Π⁺
end

function get_new_plane_ids(poly)
  id = minimum( get_plane_ids(poly.data) )
  id-1,id-2,id-3
end

## Plane sharing stuff

function link_planes!(surf::Polyhedron,stl::STL,Πr,Πf;atol)
  planes = surf.data.plane_to_ids
  v_to_f = surf.data.vertex_to_original_faces
  f_to_v = get_face_vertices(stl)
  c = array_cache(f_to_v)

  vertices = Int32[]
  Π_to_faces = [ Int32[] for _ in 1:length(planes) ]
  for (i,Π) in enumerate(planes)
    Π > 0 || continue
    faces = Π_to_faces[i]
    dists = get_plane_distances(surf.data,Π)
    empty!(vertices)
    for v in 1:num_vertices(surf)
      @assert isactive(surf,v)
      if abs(dists[v]) < atol/10
        push!(vertices,v)
        dists[v] = 0
      end
    end
    for v in vertices
      for f in v_to_f[v]
        Π ≠ f || continue
        f ∉ faces || continue
        f ∈ planes || continue
        face_on_plane = true
        for _v in getindex!(c,f_to_v,f)
          if !any( i-> first(v_to_f[i]) == _v, vertices )
            face_on_plane = false
            break
          end
        end
        if face_on_plane
          push!(faces,f)
        end
      end
    end
  end

  if maximum(length,Π_to_faces) == 0
    return fill(UNSET,length(planes)), planes
  end

  ## Link coplanar planes
  Π_to_coplanar_Π = [ Int32[] for _ in 1:length(planes) ]
  D = num_dims(stl)
  facedims = get_facedims(stl)
  rf_offset = get_offset(stl,D-1)
  f_offset = get_offset(stl,D)
  for (i,Π) in enumerate(planes)
    Π > 0 || continue
  #  get_facedims(stl_topo)[Π] == D-1 || continue
    faces = Π_to_faces[i]
    for f in faces
      d = facedims[f]
      if d == D-1
        j = findfirst(isequal(f),planes)
        if Π ∈ Π_to_faces[j] && distance_between_planes(surf,Π,f,abs,abs) < atol
          push!(Π_to_coplanar_Π[i],f)
        end
      elseif d == D
        facet = get_cell(stl,f-get_offset(stl,D))
        if distance_between_planes(surf,Π,f,abs,abs) < atol
          push!(Π_to_coplanar_Π[i],f)
        end
      end
    end
  end

  ## Cross link coplanar planes
  for (i,Πi) in enumerate(planes)
    Πi > 0 || continue
    for Πj in Π_to_coplanar_Π[i]
      j = findfirst(isequal(Πj),planes)
      if Πi ∉ Π_to_coplanar_Π[j]
        push!(Π_to_coplanar_Π[j],Πi)
      end
    end
  end


  ## Join planes recursively
  Π_to_ref_Π = collect(1:length(planes))
  stack = Int32[]
  for (i,Πi) in enumerate(planes)
    Πi > 0 || continue
    Π_to_ref_Π[i] == i || continue
    empty!(stack)
    push!(stack,i)
    while !isempty(stack)
      j = pop!(stack)
      for Πk in Π_to_coplanar_Π[j]
        k = findfirst(isequal(Πk),planes)
        i ≠ k || continue
        Π_to_ref_Π[k] == k || continue
        Π_to_ref_Π[k] = i
        push!(stack,k)
      end
    end
  end

  ## Check coplanarity
  for (i,Πi) in enumerate(planes)
    Πi > 0 || continue
    for Πj in Π_to_coplanar_Π[i]
      Πi < Πj || continue
    #  @assert distance_between_planes(surf,Πi,Πj,abs,abs) < atol*1e3
    end
  end

  ## Merge faces
  for (i,Πi) in enumerate(planes)
    Πi > 0 || continue
    Π_to_ref_Π[i] ≠ i || continue
    for f in Π_to_faces[i]
      if f ∉ Π_to_faces[Π_to_ref_Π[i]]
        push!(Π_to_faces[Π_to_ref_Π[i]],f)
      end
    end
  end

  ## Correct dists
  for (i,Πi) in enumerate(planes)
    Πi > 0 || continue
    Π_to_ref_Π[i] == i || continue
    !isempty(Π_to_coplanar_Π[i]) || continue
    dists = get_plane_distances(surf.data,Πi)
    for v in 1:num_vertices(surf)
      dists[v] ≠ 0 || continue
      if any( f-> f ∈ v_to_f[v], Π_to_faces[i] )
        dists[v] = 0
      end
    end
  end

  ## Do mapping
  for i in reverse(1:length(Π_to_ref_Π))
    if Π_to_ref_Π[i] > 0 && Π_to_ref_Π[i] ≠ i
      Π_to_ref_Π[ Π_to_ref_Π[i] ] = -Π_to_ref_Π[i]  
      Π_to_ref_Π[i] = -Π_to_ref_Π[i]
    end
  end
  for i in 1:length(Π_to_ref_Π)
    if Π_to_ref_Π[i] > 0
      Π_to_ref_Π[i] = UNSET
    else
      Π_to_ref_Π[i] = -Π_to_ref_Π[i]
    end
  end
  
  ## Mark inverted planes
  for i in 1:length(Π_to_ref_Π)
    Π_to_ref_Π[i] ∉ (UNSET,i) || continue
    j = Π_to_ref_Π[i]
    fi = planes[i]
    fj = planes[j]
    Πi = facedims[fi] == D ? Πf[fi-f_offset] : Πr[fi-rf_offset]
    Πj = facedims[fj] == D ? Πf[fj-f_offset] : Πr[fj-rf_offset]
    if relative_orientation(Πi,Πj;atol) < 0
      Π_to_ref_Π[i] = -Π_to_ref_Π[i]
    end
  end

  # TODO: 
  # reduce to cell plane (cartesian dist)
  # keep the facet planes as far as possible (just the distances)
  # All reflex coplanar planes must have the same relative orientation (modify)
  Π_to_ref_Π, planes
end

function relative_orientation(Π1::Plane,Π2::Plane;atol)
  d = normal(Π1) ⋅ normal(Π2)
  @assert abs(d) > atol
  sign(d)
end

function distance_between_planes(poly::Polyhedron,Π1,Π2,f1::Function=identity,f2::Function=identity)
  dist1 = get_plane_distances(poly.data,Π1)
  dist2 = get_plane_distances(poly.data,Π2)
  max_dist = 0.0
  for v in 1:num_vertices(poly)
    _d = abs(f1(dist1[v]) - f2(dist2[v]))
    max_dist = max(max_dist,_d)
  end
  max_dist 
end

function set_linked_planes!(poly::Polyhedron,Π_to_ref_Π,planes)
  _Π_to_ref_Π = poly.data.plane_to_ref_plane
  _planes = poly.data.plane_to_ids
  Π_to_v_to_dist = poly.data.plane_to_vertex_to_distances
  for (i,Π) in enumerate(planes)
    Π_to_ref_Π[i] ≠ UNSET || continue
    ref_Π = planes[abs(Π_to_ref_Π[i])]
    j = findfirst(isequal(Π),_planes)
    ref_j = findfirst(isequal(ref_Π),_planes)
    _Π_to_ref_Π[j] = ref_j * sign(Π_to_ref_Π[i])

    for v in 1:length(Π_to_v_to_dist[j])
      Π_to_v_to_dist[j][v] =  Π_to_v_to_dist[ref_j][v] * sign(Π_to_ref_Π[i])
    end
  end
end

function has_coplanars(data::PolyhedronData,Π)
  i = findfirst(isequal(Π),get_plane_ids(data))
  data.plane_to_ref_plane[i] ≠ UNSET
end

function are_coplanar(data::PolyhedronData,Πi,Πj)
  i = findfirst(isequal(Πi),get_plane_ids(data))
  j = findfirst(isequal(Πj),get_plane_ids(data))
  data.plane_to_ref_plane[i] == data.plane_to_ref_plane[j]
end

function add_plane!(data::PolyhedronData,Π)
  i = findfirst(isequal(Π),get_plane_ids(data))
  @assert data.plane_to_ref_plane[i] ≠ UNSET
  Πref = get_plane_ids(data)[abs(data.plane_to_ref_plane[i])]
  Πlast = UNSET
  for j in reverse(1:length(get_plane_ids(data)))
    jref = abs(data.plane_to_ref_plane[j])
    jref ≠ UNSET || continue
    if Πref == get_plane_ids(data)[jref]
      Πj = get_plane_ids(data)[j]
      for v in 1:length(data.vertex_to_planes)
        if Πj in data.vertex_to_planes[v]
          Πlast = Πj
          break
        end
      end
      Πlast == UNSET || break
    end
  end
  Πlast ≠ UNSET || return
  for v in 1:length(data.vertex_to_planes)
    if Πlast in data.vertex_to_planes[v]
      push!(data.vertex_to_planes[v],Π)
    end
  end
end

function contains_coplanars(data::PolyhedronData,Π)
  i = findfirst(isequal(Π),get_plane_ids(data))

  for (j,Πj) in enumerate(get_plane_ids(data))
    i ≠ j || continue
    if abs(data.plane_to_ref_plane[i]) == abs(data.plane_to_ref_plane[j])
      for (v,planes) in enumerate(data.vertex_to_planes)
        if Πj ∈ planes
          return true
        end
      end
    end
  end
  false
end

function add_missing_facets(
  surf::Polyhedron,
  stl::STL,
  facets,
  reflex_faces,
  empty_facets)

  D = num_dims(stl)
  rf_offset = get_offset(stl,D-1)
  f_offset = get_offset(stl,D)
  for f in facets
    has_coplanars(surf.data,f) || continue
    contains_coplanars(surf.data,f) || continue
    for rf in get_faces(stl,D,D-1)[f-f_offset]
      rf += rf_offset
      rf ∉ reflex_faces || continue
      has_original_reflex_face(surf,rf,empty=false) || continue
      for neig_f in get_faces(stl,D-1,D)[rf-rf_offset]
        neig_f += f_offset
        neig_f ∉ facets || continue
        if has_original_facet(surf,neig_f,empty=true)
          if neig_f ∈ empty_facets
            for e in get_faces(stl,D,D-1)[neig_f-f_offset]
              e += rf_offset
              e ≠ rf || continue
              for neig_neig_f in get_faces(stl,D-1,D)[e-rf_offset]
                neig_neig_f += f_offset
                neig_neig_f ≠ neig_f || continue
                neig_neig_f ∉ facets || continue
                if has_original_facet(surf,neig_neig_f,empty=true)
                  @notimplementedif neig_neig_f ∈ empty_facets
                  push!(facets,neig_neig_f)
                end
              end
            end
          else
            push!(facets,neig_f)
          end
        end
      end
    end
  end
  facets
end

function split_reflex_face(
  S,K,
  surf::Polyhedron,
  cell::Polyhedron,
  stl::STL,
  reflex_face::Integer,
  empty_facets::AbstractVector)

  # TODO: consider reflex faces sharing planes
  D = num_dims(stl)
  rface_offset = get_offset(stl,D-1)
  facet_offset = get_offset(stl,D)
  neig_facets = get_faces(stl,D-1,D)[reflex_face-rface_offset]
  if !any(i->has_original_facet(surf,i+facet_offset),neig_facets) ||
     !has_original_reflex_face(surf,reflex_face,empty=false)

    Sr,Kr = [surf],[cell]
  else
    @notimplementedif count(in(empty_facets),neig_facets) > 1
    cond = 
      f->!has_original_facet(surf,f+facet_offset) || 
        f+facet_offset∈empty_facets
    j = findfirst(cond,neig_facets)
    !isnothing(j) || 
      throw(
        ErrorException("One of these facets may be degenerate: $neig_facets"))

    missing_facet = neig_facets[j] + facet_offset
    _surf = one_face_polyhedron(surf,missing_facet)
    j = findfirst(i->isnothing(i)||!has_facets(i),S)
    Sr,Kr = [surf,surf],K
    Sr[j] = _surf
  end
  Sr,Kr
end

function one_face_polyhedron(poly::Polyhedron,face::Integer)
  v_to_f = get_data(poly).vertex_to_original_faces
  nodes = Int32[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    if face ∈ v_to_f[v]
      push!(nodes,v)
    end
  end
  sort!(nodes)
  r = restrict(poly,nodes)
  for v in 1:num_vertices(r) 
    r.data.vertex_to_original_faces[v] = [face]
  end
  r
end

function get_facet_to_isempty(stl::STL;atol)
  f_to_isempty = falses(num_cells(stl))
  c = get_cell_cache(stl)
  for f in 1:num_cells(stl)
    facet = get_cell!(c,stl,f)
    if min_height(facet) < atol
      f_to_isempty[f] = true
    end
  end
  f_to_isempty
end

function correct_small_facets_planes!(stl::STL,Πf,f_to_isempty;atol)
  D = num_dims(stl)
  e_to_f = get_faces(stl,D-1,D)
  f_to_e = get_faces(stl,D,D-1)
  full_facets = Int32[]
  queue = Int32[]
  num_new_planes = 0
  Πnew = empty(Πf)
  f_to_new_plane = fill(UNSET,num_cells(stl))
  ec = array_cache(f_to_e)
  fc = array_cache(e_to_f)
  for f in 1:num_cells(stl)
    f_to_isempty[f] || continue
    f_to_new_plane[f] == UNSET || continue 
    num_new_planes += 1
    f_to_new_plane[f] = num_new_planes
    empty!(full_facets)
    empty!(queue)
    push!(queue,f)
    head = 1
    while length(queue) ≥ head 
      f_curr = queue[head]
      head += 1
      for e in getindex!(ec,f_to_e,f_curr)
        for f_neig in getindex!(fc,e_to_f,e)
          f_neig ≠ f_curr || continue
          if f_to_isempty[f_neig]
            f_to_new_plane[f_neig] == UNSET || continue
            f_to_new_plane[f_neig] = num_new_planes
            push!(queue,f_neig)
          elseif f_to_new_plane[f_neig] ≠ num_new_planes 
            f_to_new_plane[f_neig] == num_new_planes 
            push!(full_facets,f_neig)
          end
        end
      end
    end
    empty_facets = queue
    !isempty(full_facets) || @unreachable
    n = zero( normal(first(Πf)) )
    c = zero( center(first(Πf)) )
    for full_f in full_facets
      n += normal(Πf[full_f]) / length(full_facets)
    end
    for empty_f in empty_facets
      c += center(Πf[empty_f]) / length(empty_facets)
    end
    Π = Plane(c,n)
    push!(Πnew,Π)
  end
  v_coords = get_vertex_coordinates(stl)
  f_to_v = get_face_vertices(stl,D)
  c = array_cache(f_to_v)
  for f in 1:num_cells(stl)
    f_to_isempty[f] || continue
    Πf[f] = Πnew[ f_to_new_plane[f] ]
    for i in getindex!(c,f_to_v,f)
      v = v_coords[i]
      @assert abs(signed_distance(v,Πf[f])) < atol
    end
  end
  Πf
end
