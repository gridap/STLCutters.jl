
## Constants

const OPEN = -1

const FACE_CUT = -1
const FACE_IN = 1
const FACE_OUT = 2

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
  vertex_to_parent_edge::Vector{Tuple{Int,Int}}
  plane_to_vertex_to_distances::Vector{Vector{Float64}}
  plane_to_ref_plane::Vector{Int}
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
    p⁻,p⁺ = split(p,Πi)
    p = inside == inout[i] ? p⁻ : p⁺
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
  if smin > 0
    return nothing,p
  end
  if smax ≤ 0
    if has_coplanars(get_data(p),Π)
      add_plane!(get_data(p),Π)
    end
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
  in_vertices = [p.vertices;new_vertices]
  out_vertices = [p.vertices;new_vertices]
  in_data = data
  out_data = deepcopy(data)
  update_data!(in_data,in_graph,num_vertices(p))
  update_data!(out_data,out_graph,num_vertices(p))
  compact!(Polyhedron(in_vertices,in_graph,isopen(p),in_data)), 
  compact!(Polyhedron(out_vertices,out_graph,isopen(p),out_data))
end

function decompose(surf::Polyhedron,cell::Polyhedron,rfaces,empty_facets,stl::DiscreteModel)
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
  vertex_coordinates = get_vertex_coordinates(poly)
  T = Vector{Int}[]
  X = empty(vertex_coordinates)
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
          kv = v_to_pv[k]
          #length(unique(kv)) == length(k) || break
          ki = length(X) .+ (1:length(k))
          push!(X,vertex_coordinates[k]...)
          push!(T,ki)
        end
      end
    end
  end
  T,X
end

function simplexify_surface(poly::Polyhedron{3})
  stack = Int[]
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

function simplexify_boundary(poly::Polyhedron{3},stl::DiscreteModel)
  stack = Int[]
  v_to_pv = get_data(poly).vertex_to_parent_vertex
  v_to_Π = get_data(poly).vertex_to_planes
  facedims = get_facedims(get_grid_topology(stl))
  facet_list = Int[]
  istouch = map( i -> zeros(Bool,length(i)), get_graph(poly) )
  vertex_coordinates = get_vertex_coordinates(poly)
  T = Vector{Int}[]
  X = empty(vertex_coordinates)
  f_to_stlf = Int[]
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
          # @assert length(facet_list) == 1 # Can have more than one (coplanars)
          k = [v,vcurrent,vnext]
          kv = v_to_pv[k]
          ki = length(X) .+ (1:length(k))
          push!(X,vertex_coordinates[k]...)
          push!(T,ki)
          push!(f_to_stlf,first(facet_list))
        end
      end
    end
  end
  T,X,f_to_stlf
end

function simplexify(polys::AbstractVector{<:Polyhedron{Dp,Tp}}) where {Dp,Tp}
  T = Vector{Int}[]
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
  stl::DiscreteModel) where {Dp,Tp}

  T = Vector{Int}[]
  X = Point{Dp,Tp}[]
  f_to_stlf = Int[]
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

get_graph(a::Polyhedron) = a.edge_vertex_graph

get_data(a::Polyhedron) = a.data

Base.isopen(a::Polyhedron) = a.isopen

get_vertex_to_planes(a::PolyhedronData) = a.vertex_to_planes

function isactive(p::Polyhedron,vertex::Integer)
  !isempty( get_graph(p)[vertex] )
end

function plot(p::Polyhedron,filename=nothing)
  vertices = Int[]
  for i in 1:num_vertices(p)
    if isactive(p,i)
      push!(vertices,i)
    end
  end
  vertex_to_node = Dict( zip( vertices, 1:length(vertices) ) )
  g = empty( get_graph(p) )
  for (i,vneigs) in enumerate( get_graph(p) )
    isactive(p,i) || continue
    nodes = Int[]
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
  graphplot(g,names=names,curves=false;kwargs...)
  isnothing(filename) || savefig(filename)
  plot!()
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

function get_original_reflex_faces(p::Polyhedron{D},stl::DiscreteModel;empty=true) where D
  get_original_faces(p,stl,Val{D-2}();empty)
end

function get_original_facets(p::Polyhedron{D},stl::DiscreteModel;empty=false) where D
  get_original_faces(p,stl,Val{D-1}();empty)
end

function get_original_faces(p::Polyhedron,stl::DiscreteModel,::Val{d};empty) where d
  faces = collect_original_faces(p,stl,d)
  filter!(f->has_original_face(p,f,Val{d}();empty),faces)
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


function compute_cell_to_facets(grid::CartesianGrid,stl::DiscreteModel)
  desc = get_cartesian_descriptor(grid)
  @assert length(get_reffes(grid)) == 1
  p = get_polytope(get_cell_reffe(grid)[1])
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
      if fast_intersection(f,_pmin,_pmax,p)
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
  v_to_e = fill((UNSET,UNSET),num_vertices)
  Π_to_v_to_d = Vector{Int}[]
  Π_to_rΠ = Int[]
  Π_to_id = Int[]
  PolyhedronData( v_to_Π, v_to_of, v_to_v, v_to_e, Π_to_v_to_d, Π_to_rΠ, Π_to_id  )
end

function polyhedron_data(p::Polytope)
  v_to_Π = -get_faces(p,0,num_dims(p)-1)
  v_to_of = [ Int[] for _ in 1:num_vertices(p) ]
  v_to_v = collect(1:num_vertices(p))
  v_to_e = fill((UNSET,UNSET),num_vertices(p))
  Π_to_v_to_d = Vector{Int}[]
  Π_to_rΠ = Int[]
  Π_to_id = Int[]
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

function disconnect_graph!(edge_graph,num_vertices,distances,mask::Bool)
  for i in 1:num_vertices
    if mask == ( distances[i] > 0 )
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

function add_vertex!(data::PolyhedronData,v1::Integer,v2::Integer,Πid::Integer)
  v_to_Π = data.vertex_to_planes
  v_to_f = data.vertex_to_original_faces
  v_to_pv = data.vertex_to_parent_vertex
  v_to_pe = data.vertex_to_parent_edge
  Πs = intersect( v_to_Π[v1], v_to_Π[v2] )
  fs = intersect( v_to_f[v1], v_to_f[v2] )
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
    for n in get_cell_node_ids(stl)[f]
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

function bisector_plane(
   stl::DiscreteModel{Dc,Dp},
   d::Integer,
   dface::Integer,
   Πf::AbstractArray) where {Dc,Dp}

  @notimplementedif Dc ≠ Dp-1
  @notimplementedif d ≠ Dc-1
  facets = get_faces(get_grid_topology(stl),d,Dc)[dface]
  length(facets) == 2 || return Πf[ only(facets) ]
  edge = get_dface(stl,dface,Val{Dc-1}())
  Π1 = Πf[ facets[1] ]
  Π2 = Πf[ facets[2] ]
  bisector_plane(edge,Π1,Π2)
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

function get_reflex_planes(stl::DiscreteModel) 
  Dc = num_dims(stl)
  Dp = num_point_dims(stl)
  T_Π = typeof( bisector_plane(stl,Dc-1,1) )
  Π_r = Vector{T_Π}(undef,num_faces(stl,Dc-1))
  for f in 1:num_faces(stl,Dc-1)
    Π_r[f] = bisector_plane(stl,Dc-1,f)
  end
  Π_r
end

function get_reflex_planes(stl::DiscreteModel,Πf) 
  Dc = num_dims(stl)
  Dp = num_point_dims(stl)
  T_Π = typeof( bisector_plane(stl,Dc-1,1,Πf) )
  Π_r = Vector{T_Π}(undef,num_faces(stl,Dc-1))
  for f in 1:num_faces(stl,Dc-1)
    Π_r[f] = bisector_plane(stl,Dc-1,f,Πf)
  end
  Π_r
end

function get_facet_plane(stl::DiscreteModel,i)
  f = get_cell(stl,i)
  Plane(f)
end

function get_facet_planes(stl::DiscreteModel)
  [ get_facet_plane(stl,i) for i in 1:num_cells(stl) ]
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
  group = Int[]
  ids = Int[]
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
  group_to_facets = [ Int[] for _ in 1:num_groups ]
  for (i,f) in enumerate(facets)
    g = p_to_group[f_to_part[i]]
    push!(group_to_facets[g],f)
  end
  group_to_facets
end


function get_disconnected_facets(poly::Polyhedron,stl::DiscreteModel)
  v_to_pv = poly.data.vertex_to_parent_vertex
  v_to_f = poly.data.vertex_to_original_faces
  facets = get_original_facets(poly,stl,empty=true)
  facet_to_part = fill(UNSET,length(facets))
  stack = Int[]
  num_parts = 0
  for (i,facet) in enumerate(facets)
    facet_to_part[i] == UNSET || continue
    num_parts += 1
    facet_to_part[i] = num_parts
    empty!(stack)
    push!(stack,facet)
    while !isempty(stack)
      current_facet = pop!(stack)
      v,vnext = get_facet_onset(poly,stl,current_facet) # reduce calls
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
  part_to_facets = [ Int[] for _ in 1:num_parts ]
  for (i,facet) in enumerate(facets)
    push!(part_to_facets[facet_to_part[i]],facet)
  end
  part_to_facets
end

function get_facet_onset(poly::Polyhedron,stl::DiscreteModel,facet::Integer)
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

function is_reflex(poly::Polyhedron,stl::DiscreteModel,reflex_face;inside)
  stl_topo = get_grid_topology(stl)
  Dc = num_dims(stl)
  rf = reflex_face - get_offset(stl_topo,Dc-1)
  length(get_faces(stl_topo,Dc-1,Dc)[rf]) == 2 || return false
  f1 = get_faces(stl_topo,Dc-1,Dc)[rf][1] + get_offset(stl_topo,Dc)
  f2 = get_faces(stl_topo,Dc-1,Dc)[rf][2] + get_offset(stl_topo,Dc)
  has_plane(poly.data,f1) || return false
  has_plane(poly.data,f2) || return false
  is_facet_in_facet(poly,f1,f2;inside) || return true
  is_facet_in_facet(poly,f2,f1;inside) || return true
  false
end

function refine(
  K::Polyhedron,
  Γ::Polyhedron,
  stl::DiscreteModel,
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
  stl::DiscreteModel,
  reflex_planes::AbstractVector,
  reflex_faces::AbstractVector{<:Integer},
  facet_planes::AbstractVector,
  facets::AbstractVector{<:Integer})
  
  Dc = num_dims(stl)
  Πr = view(reflex_planes,reflex_faces.-get_offset(get_grid_topology(stl),Dc-1))
  Πf = view(facet_planes,facets.-get_offset(get_grid_topology(stl),Dc))
  (Πr,Πf),(reflex_faces,facets)
end

function compute_distances!(poly::Polyhedron,Πn::Tuple,Πn_ids::Tuple)
  for (Π,Π_ids) in zip(Πn,Πn_ids)
    compute_distances!(poly,Π,Π_ids)
  end
end

function get_cell_nodes_to_inout(polys_in,polys_out,p::Polytope)
  node_to_inout = fill(UNSET,num_vertices(p))
  
  complete_nodes_to_inout!(node_to_inout,polys_in,FACE_IN,p)
  complete_nodes_to_inout!(node_to_inout,polys_out,FACE_OUT,p)
  #@assert UNSET ∉ node_to_inout
  node_to_inout
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

function Base.eps(T::Type{<:AbstractFloat},grid::Grid)
  pmin,pmax = get_bounding_box(grid)
  vmax = max(abs.(Tuple(pmin))...,abs.(Tuple(pmax))...)
  eps(T(vmax))
end

Base.eps(grid::Grid) = eps(Float64,grid)


function compute_submesh(grid::CartesianGrid,stl::DiscreteModel;kdtree=false)
  D = num_dims(grid)
  atol = eps(grid)*1e3

  f_to_isempty = get_facet_to_isempty(stl;atol)
  Πf = get_facet_planes(stl)
  Πf = correct_small_facets_planes!(stl,Πf,f_to_isempty;atol)
  Πr = get_reflex_planes(stl,Πf)

  c_to_stlf = compute_cell_to_facets(grid,stl)

  Γ0 = Polyhedron(stl)

  T = Vector{Int}[]
  F = Vector{Int}[]
  X = empty(get_node_coordinates(stl))
  Xf = empty(get_node_coordinates(stl))
  k_to_io = Int8[]
  k_to_bgcell = Int[]
  f_to_bgcell = Int[]
  f_to_stlf = Int[]

  cell_facets = Int[]
  bgcell_to_ioc = fill(UNSET,num_cells(grid))
  bgnode_to_io = fill(UNSET,num_nodes(grid))

  for cell in 1:num_cells(grid)
    !isempty(c_to_stlf[cell]) || continue
    
    pmin = get_cell_coordinates(grid)[cell][1]-atol
    pmax = get_cell_coordinates(grid)[cell][end]+atol

    empty!(cell_facets)
    for f in c_to_stlf[cell] 
      for v in get_faces(get_grid_topology(stl),D-1,0)[f]
        for _f in get_faces(get_grid_topology(stl),0,D-1)[v]
          if _f ∉ cell_facets
            push!(cell_facets,_f)
          end
        end
      end
    end

    p = get_polytope(get_cell_reffe(grid)[cell])
    v = map(i->Point(float.(Tuple(i))),get_cell_coordinates(grid)[cell])
    K = Polyhedron(p,v)
    Γk0 = restrict(Γ0,stl,cell_facets)

    stl_reflex_faces_k = Int[]
    f_to_r = get_faces(stl.grid_topology,D-1,D-2)
    r_to_f = get_faces(stl.grid_topology,D-2,D-1)
    for f in c_to_stlf[cell]
      for r in f_to_r[f]
        _r = r+get_offset(get_grid_topology(stl),D-2)
        if _r ∉ stl_reflex_faces_k &&
           all(i->i∈c_to_stlf[cell],r_to_f[r])

          push!(stl_reflex_faces_k,_r)
        end
      end
    end

    stl_facets_k = c_to_stlf[cell].+get_offset(get_grid_topology(stl),num_dims(stl))

    empty_facets = filter(f->f_to_isempty[f-get_offset(get_grid_topology(stl),num_dims(stl))],stl_facets_k)

    Πk,Πk_ids,Πk_io = get_cell_planes(p,pmin,pmax)
    Πkf,Πkf_ids = filter_face_planes(stl,Πr,stl_reflex_faces_k,Πf,stl_facets_k)

    compute_distances!(Γk0,(Πk,Πkf),(Πk_ids,Πkf_ids))
    compute_distances!(K,Πkf,Πkf_ids)

    Π_to_refΠ,Πs = link_planes!(Γk0,stl,Πr,Πf;atol)
    set_linked_planes!(Γk0,Π_to_refΠ,Πs)
    set_linked_planes!(K,Π_to_refΠ,Πs)

    Γk = clip(Γk0,Πk_ids,inout=Πk_io)
    !isnothing(Γk) || continue

    if isempty(get_original_facets(Γk,stl))
      continue
    end

    if kdtree
      Γv,Kv = refine_by_vertices(Γk,K,atol)
    else
      Γv,Kv = [Γk],[K]
    end

    Kn_in = typeof(K)[]
    Kn_out = typeof(K)[]
    for (γk,k) in zip(Γv,Kv)
      k_in = refine(k,γk,stl,stl_reflex_faces_k,empty_facets,inside=true)
      k_out = refine(k,γk,stl,stl_reflex_faces_k,empty_facets,inside=false)
      append!(Kn_in,k_in)
      append!(Kn_out,k_out)
    end

    Tin,Xin = simplexify(Kn_in)
    Tout,Xout = simplexify(Kn_out)
    T_Γ,X_Γ,f_to_f = simplexify_boundary(Kn_in,stl)
    T_Γk,X_Γk = simplexify(Γk)

    bgcell_to_ioc[cell] = FACE_CUT

    n_to_io = get_cell_nodes_to_inout(Kn_in,Kn_out,p)

    append!(T, map(i->i.+length(X),Tin) )
    append!(X,Xin)
    append!(k_to_io,fill(FACE_IN,length(Tin)))
    append!(k_to_bgcell,fill(cell,length(Tin)))
    append!(T, map(i->i.+length(X),Tout) )
    append!(X,Xout)
    append!(k_to_io,fill(FACE_OUT,length(Tout)))
    append!(k_to_bgcell,fill(cell,length(Tout)))
    append!(F, map(i->i.+length(Xf),T_Γ) ) 
    append!(Xf,X_Γ)
    append!(f_to_bgcell,fill(cell,length(T_Γ)))
    append!(f_to_stlf,f_to_f.-get_offset(get_grid_topology(stl),num_dims(stl)))
    
    for (i,bg_v) in enumerate(get_cell_node_ids(grid)[cell])
      n_to_io[i] ≠ UNSET || continue
      n_to_io[i] ≠ FACE_CUT || continue
      bgnode_to_io[bg_v] ≠ FACE_CUT || continue
      # @assert bgnode_to_io[bg_v] ∈ (UNSET,n_to_io[i])
      bgnode_to_io[bg_v] = n_to_io[i]
    end
  end
  grid_topology = GridTopology(grid)
  stack = Int[]
  for cell in 1:num_cells(grid)
    if bgcell_to_ioc[cell] == FACE_CUT
      resize!(stack,0)
      push!(stack,cell)
      while length(stack) > 0
        current_cell = pop!(stack)
        for node in get_cell_node_ids(grid)[current_cell]
          if bgnode_to_io[node] == FACE_IN
            for neig_cell in get_faces(grid_topology,0,num_dims(grid))[node]
              if bgcell_to_ioc[neig_cell] == UNSET
                bgcell_to_ioc[neig_cell] = FACE_IN
                for neig_node in get_cell_node_ids(grid)[neig_cell]
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
  for cell in 1:num_cells(grid)
    if bgcell_to_ioc[cell] == UNSET
      bgcell_to_ioc[cell] = FACE_OUT
    end
  end
  T,X,F,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf,bgcell_to_ioc
end

## Kd Tree stuff

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

function link_planes!(surf::Polyhedron,stl::DiscreteModel,Πr,Πf;atol)
  planes = surf.data.plane_to_ids
  v_to_f = surf.data.vertex_to_original_faces
  stl_topo = get_grid_topology(stl)
  f_to_v = get_face_vertices(stl_topo)

  vertices = Int[]
  Π_to_faces = [ Int[] for _ in 1:length(planes) ]
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
        for _v in f_to_v[f]
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
  Π_to_coplanar_Π = [ Int[] for _ in 1:length(planes) ]
  D = num_dims(stl)
  facedims = get_facedims(stl_topo)
  rf_offset = get_offset(stl_topo,D-1)
  f_offset = get_offset(stl_topo,D)
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
        facet = get_cell(stl,f-get_offset(stl_topo,D))
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
  stack = Int[]
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
  stl::DiscreteModel,
  facets,
  reflex_faces,
  empty_facets)

  D = num_dims(stl)
  stl_topo = get_grid_topology(stl)
  rf_offset = get_offset(stl_topo,D-1)
  f_offset = get_offset(stl_topo,D)
  for f in facets
    has_coplanars(surf.data,f) || continue
    contains_coplanars(surf.data,f) || continue
    for rf in get_faces(stl_topo,D,D-1)[f-f_offset]
      rf += rf_offset
      rf ∉ reflex_faces || continue
      has_original_reflex_face(surf,rf,empty=false) || continue
      for neig_f in get_faces(stl_topo,D-1,D)[rf-rf_offset]
        neig_f += f_offset
        neig_f ∉ facets || continue
        if has_original_facet(surf,neig_f,empty=true)
          if neig_f ∈ empty_facets
            for e in get_faces(stl_topo,D,D-1)[neig_f-f_offset]
              e += rf_offset
              e ≠ rf || continue
              for neig_neig_f in get_faces(stl_topo,D-1,D)[e-rf_offset]
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
  stl::DiscreteModel,
  reflex_face::Integer,
  empty_facets::AbstractVector)

  # TODO: consider reflex faces sharing planes
  D = num_dims(stl)
  stl_topo = get_grid_topology(stl)
  rface_offset = get_offset(stl_topo,D-1)
  facet_offset = get_offset(stl_topo,D)
  neig_facets = get_faces(stl_topo,D-1,D)[reflex_face-rface_offset]
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
  nodes = Int[]
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

function get_facet_to_isempty(stl::DiscreteModel;atol)
  f_to_isempty = falses(num_cells(stl))
  for f in 1:num_cells(stl)
    facet = get_cell(stl,f)
    if min_height(facet) < atol
      f_to_isempty[f] = true
    end
  end
  f_to_isempty
end

function correct_small_facets_planes!(stl::DiscreteModel,Πf,f_to_isempty;atol)
  D = num_dims(stl)
  stl_topo = get_grid_topology(stl)
  e_to_f = get_faces(stl_topo,D-1,D)
  f_to_e = get_faces(stl_topo,D,D-1)
  full_facets = Int[]
  queue = Int[]
  num_new_planes = 0
  Πnew = empty(Πf)
  f_to_new_plane = fill(UNSET,num_cells(stl))
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
      for e in f_to_e[f_curr]
        for f_neig in e_to_f[e]
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
  for f in 1:num_cells(stl)
    f_to_isempty[f] || continue
    Πf[f] = Πnew[ f_to_new_plane[f] ]
    for i in get_face_nodes(stl,D)[f]
      v = get_node_coordinates(stl)[i]
      @assert abs(signed_distance(v,Πf[f])) < atol
    end
  end
  Πf
end
