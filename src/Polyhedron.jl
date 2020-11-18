module PolyhedraTests

using Test

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using Gridap.Helpers

using STLCutters

using STLCutters: Plane,Face
import STLCutters: signed_distance 
import STLCutters: compute_stl_model 
import STLCutters: compute_grid 
import STLCutters: get_cell 
import STLCutters: get_dface 
import STLCutters: normal 
import STLCutters: center 
import STLCutters: volume 

import Gridap: num_dims
import Gridap.ReferenceFEs: num_vertices
import Gridap.ReferenceFEs: get_vertex_coordinates

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
end

struct CartesianPlane{D,T}
  d::Int8
  value::T
  positive::Bool
end

## Constructors

function Polyhedron(stl::DiscreteModel{Dc,Dp},cell_facets) where {Dc,Dp}
  𝓖,v_to_n = compute_graph( stl, cell_facets )
  X = get_node_coordinates(stl)[v_to_n]
  p = Polyhedron(X,𝓖,isopen=true)
  set_faces!(p,stl,v_to_n,Dc-1)
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

## Polyhedra operations

function clip(poly::Polyhedron,Π::AbstractVector)
  p = poly
  for (i,Πi) ∈ enumerate(Π)
    p,_ = split(p,Πi,i)
    p !== nothing || break
  end
  p
end

function split(p::Polyhedron,Π,Πid=UNSET)
  smin = Inf
  smax = -Inf
  distances = zeros(num_vertices(p))
  new_vertices = empty( get_vertex_coordinates(p) )
  for v in 1:num_vertices(p)
    isactive(p,v) || continue
    distances[v] = signed_distance(p[v],Π)
    smin = min(smin,distances[v])
    smax = max(smax,distances[v])
  end
  if smin ≥ 0
    return nothing,p
  end
  if smax ≤ 0
    return p,nothing
  end
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
      add_vertex!(data,v,vneig,Πid)

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
  disconnect_graph!(in_graph,distances,true)
  disconnect_graph!(out_graph,distances,false)
  add_open_vertices!(in_graph,data,num_vertices(p))
  add_open_vertices!(out_graph,data,num_vertices(p))
  vertices = [p.vertices;new_vertices]
  in_data = data
  out_data = deepcopy(data)
  Polyhedron(vertices,in_graph,in_data), Polyhedron(vertices,out_graph,out_data)
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
    if is_edge_oriented(p,cell,cut_edge,poly,first(e_to_v[cut_edge]))
      vstart = get_face_vertices(p,1)[cut_edge][1]
    else
      vstart = get_face_vertices(p,1)[cut_edge][2]
    end
    isinside = zeros(Bool,num_vertices(cell)) # These can be stored on Int
    istouch = zeros(Bool,num_vertices(cell))
    istouch[vstart] = true
    isinside[vstart] = true
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
        istarts = isinside[edge[1]] ? istats : reverse(istarts)
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
  delete_open_vertices!(graph_in,p)
  delete_open_vertices!(graph_out,p)
  update_open_data!(graph_in,data_in)
  update_open_data!(graph_out,data_out)
  Polyhedron(vertices,graph_in,data_in), Polyhedron(vertices,graph_out,data_out)
end

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

function simplexify(poly::Polyhedron{3})
  istouch = map( i -> zeros(Bool,length(i)), get_graph(poly) )
  vstart = UNSET
  T = Vector{Int}[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    vstart = vstart == UNSET ? v : vstart
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
        if v ∉ (vstart,vcurrent,vnext)
          k = [vstart,v,vcurrent,vnext]
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
      link!(graph,vs[1],vs[end],vs[2],UNSET)
      link!(graph,vs[end],vs[1],vs[end-1],UNSET)
      for i in 3:2:length(vs)-3
        link!(graph,vs[i],UNSET,vs[i+1],UNSET)
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

function get_stl_edges(p::Polyhedron)
  stack = Int[]
  istouch = zeros(Bool,num_vertices(p))
  v_to_f = get_data(p).vertex_to_original_faces
  stl_f = Int[]
  for v in 1:num_vertices(p)
    if !istouch[v] && isactive(p,v)
      empty!(stack)
      push!(stack,v)
      while !isempty(stack)
        vcurrent = pop!(stack)
        !istouch[vcurrent] || continue
        istouch[vcurrent] = true
        for vneig in get_graph(p)[vcurrent]
          vneig ≠ OPEN || continue
          if !istouch[vneig]
            push!(stack,vneig)
            for f in v_to_f[vcurrent]
              if f ∈ v_to_f[vneig] && f ∉ stl_f
                push!(stl_f,f)
              end
            end
          end
        end
      end
    end
  end
  stl_f
end

function has_stl_edge(p::Polyhedron,edge::Integer)
  v_to_f = get_data(p).vertex_to_original_faces
  for v in 1:num_vertices(p)
    if isactive(p,v) && edge ∈ v_to_f[v]
      for vneig in get_graph(p)[v]
        vneig ≠ OPEN || continue
        if edge ∈ v_to_f[vneig]
          # check length > 0
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

## Helpers

function polyhedron_data(num_vertices::Integer;isopen=false)
  v_to_Π = [ Int[] for _ in 1:num_vertices ]
  v_to_of = [ Int[] for _ in 1:num_vertices ]
  v_to_o = fill(isopen,num_vertices)
  PolyhedronData( v_to_Π, v_to_of, v_to_o )
end

function polyhedron_data(p::Polytope)
  v_to_Π = deepcopy( get_faces(p,0,num_dims(p)-1) )
  v_to_of = [ Int[] for _ in 1:num_vertices(p) ]
  v_to_o = fill(false,num_vertices(p))
  PolyhedronData( v_to_Π, v_to_of, v_to_o )
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

function disconnect_graph!(edge_graph,distances,mask::Bool)
  for i in 1:length(distances)
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
  data
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

function set_faces!(data::PolyhedronData,stl::DiscreteModel,nodes,d::Integer)
  topo = get_grid_topology(stl)
  v_to_f = data.vertex_to_original_faces
  for v in 1:length(v_to_f)
    v_to_f[v] = get_faces(topo,0,d)[nodes[v]]
  end
  data
end

function set_faces!(p::Polyhedron,a...)
  set_faces!(get_data(p),a...)
  p
end

function next_vertex(p::Polyhedron,vprevious::Integer,vcurrent::Integer)
  i = findfirst( isequal(vprevious), get_graph(p)[vcurrent] )
  vnext = OPEN
  while vnext == OPEN
    i = ( i % length( get_graph(p)[vcurrent] ) ) + 1
    vnext = get_graph(p)[vcurrent][ i ]
  end
  vnext
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
  [ CartesianPlane(pmin+iseven(i)*pmax,D-((i-1)÷2),(-1)^i) for i in 1:2*D ]
end

function delete_open_vertices!(graph,p::Polytope)
  delete_open_vertices!(graph,1:num_vertices(p))
end

function delete_open_vertices!(graph,vrange)
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

function bisector(stl::DiscreteModel{Dc,Dp},d::Integer,dface::Integer) where {Dc,Dp}
  @notimplementedif Dc ≠ Dp-1
  @notimplementedif d ≠ Dc-1
  facets = get_faces(get_grid_topology(stl),d,Dc)[dface]
  @assert length(facets) == 2
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
  
function get_reflex_planes(stl::DiscreteModel)
  Dc = num_dims(stl)
  Dp = num_point_dims(stl)
  T_Π = typeof( bisector(stl,Dc-1,1) )
  Π_r = Vector{T_Π}(undef,num_faces(stl,Dc-1))
  for f in 1:num_faces(stl,Dc-1)
    Π_r[f] = vertical_plane(stl,Dc-1,f,Dp)
  end
  Π_r
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

## Driver

# STL-like
vertices = [
  Point(0.1,-0.2,0.5),
  Point(1.0,-0.3,0.5),
  Point(-0.2,1.2,0.5),
  Point(-0.3,0.3,0.5),
  Point(0.5,1.3,0.5),
  Point(0.5,0.5,0.7),
  Point(1.2,1.2,0.5)]

facet_to_vertices = 
[[4,6,1],
 [6,5,7],
 [4,3,6],
 [1,6,2],
 [3,5,6],
 [6,7,2] ]

stl = compute_stl_model( Table(facet_to_vertices), vertices )
Π_r = get_reflex_planes(stl)
e_to_isconvex = get_convex_faces(stl)

cell_facets = 1:num_cells(stl)
pmin,pmax = Point(0,0,0),Point(1,1,1)

Γ0 = Polyhedron(stl,cell_facets)

Πc = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))
cell = Polyhedron(HEX)

Γ = clip(Γ0,Πc)

p⁻,p⁺ = merge(HEX,cell,Γ)

in_polys = decompose(p⁻,Π_r,e_to_isconvex,isinside=true)
out_polys = decompose(p⁺,Π_r,e_to_isconvex,isinside=false)

Tin,Xin = simplexify(in_polys)
Tout,Xout = simplexify(out_polys)

mesh_in = compute_grid(Table(Tin),Xin,TET)
mesh_out = compute_grid(Table(Tout),Xout,TET)

writevtk(edge_mesh(Γ0),"Gamma_0")
writevtk(edge_mesh(Γ),"Gamma")

writevtk(edge_mesh(p⁻),"poly_in")
writevtk(edge_mesh(p⁺),"poly_out")


for (i,poly) in enumerate(out_polys)
  writevtk(edge_mesh(poly),"poly_out_$i")
end

for (i,poly) in enumerate(in_polys)
  writevtk(edge_mesh(poly),"poly_in_$i")
end

writevtk(mesh_in,"simplices_in")
writevtk(mesh_out,"simplices_out")

@test volume(mesh_in) + volume(mesh_out) ≈ 1

# Sharp corner

vertices = [
  Point(0.1,-0.2,-1.0),
  Point(1.0,-0.3,-1.0),
  Point(-0.2,1.2,-1.0),
  Point(-0.3,0.3,-1.0),
  Point(0.5,1.3,-1.0),
  Point(0.5,0.5,0.5),
  Point(1.2,1.2,-1.0)]

facet_to_vertices = 
[[4,6,1],
 [6,5,7],
 [4,3,6],
 [1,6,2],
 [3,5,6],
 [6,7,2] ]

stl = compute_stl_model( Table(facet_to_vertices), vertices )
Π_r = get_reflex_planes(stl)
e_to_isconvex = get_convex_faces(stl)

cell_facets = 1:num_cells(stl)
pmin,pmax = Point(0,0,0),Point(1,1,1)

Γ0 = Polyhedron(stl,cell_facets)

Πc = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))
cell = Polyhedron(HEX)

Γ = clip(Γ0,Πc)

p⁻,p⁺ = merge(HEX,cell,Γ)

in_polys = decompose(p⁻,Π_r,e_to_isconvex,isinside=true)
out_polys = decompose(p⁺,Π_r,e_to_isconvex,isinside=false)

Tin,Xin = simplexify(in_polys)
Tout,Xout = simplexify(out_polys)

mesh_in = compute_grid(Table(Tin),Xin,TET)
mesh_out = compute_grid(Table(Tout),Xout,TET)

writevtk(edge_mesh(Γ0),"Gamma_0")
writevtk(edge_mesh(Γ),"Gamma")

writevtk(edge_mesh(p⁻),"poly_in")
writevtk(edge_mesh(p⁺),"poly_out")

for (i,poly) in enumerate(out_polys)
  writevtk(edge_mesh(poly),"poly_out_$i")
end

for (i,poly) in enumerate(in_polys)
  writevtk(edge_mesh(poly),"poly_in_$i")
end

writevtk(mesh_in,"simplices_in")
writevtk(mesh_out,"simplices_out")

@test volume(mesh_in) + volume(mesh_out) ≈ 1


#TODO: 
#  Workflow:
#    - filter faces with voxels
#    - embedd with global workflow
#    - gridap FE problem
#  
#  Improvement:
#    - force dist to be zero if it shoul be, do not reflex for zero length edge
#

end # module




