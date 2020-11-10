module PolyhedraTests

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using Gridap.Helpers

using STLCutters

using STLCutters: Plane,Face
import STLCutters: signed_distance 

import Gridap: num_dims
import Gridap.ReferenceFEs: num_vertices
import Gridap.ReferenceFEs: get_vertex_coordinates

struct Polyhedron{D,T}
  vertices::Vector{Point{D,T}}
  edge_vertex_graph::Vector{Vector{Int}}
  vertex_tags::Vector{Vector{Int8}}
end

struct CartesianPlane{D,T}
  d::Int8
  value::T
  orientation::Bool
end

function Polyhedron(vertices::AbstractVector{<:Point},graph::Vector{<:Vector})
  tags = [ Int8[] for _ in 1:length(vertices) ]
  Polyhedron(vertices,graph,tags)
end

function Polyhedron(p::Polytope{2},vertices::AbstractVector{<:Point})
  if p == TRI
    e_v_graph = [[2,3],[3,1],[1,2]]
  elseif p == QUAD
    e_v_graph = [[2, 3],[4, 1],[1, 4],[3, 2]]
  else
    @unreachable
  end
  Polyhedron(vertices,e_v_graph)
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
  Polyhedron(vertices,e_v_graph)
end

function Polyhedron(p::Polytope)
  Polyhedron(p,get_vertex_coordinates(p))
end

function CartesianPlane(p::Point{D,T},d::Integer,orientation::Bool) where {D,T}
  CartesianPlane{D,T}(Int8(d),p[d],orientation)
end

num_dims(::Polyhedron{D}) where D = D

num_vertices(a::Polyhedron) = length(a.vertices)

get_vertex_coordinates(a::Polyhedron) = a.vertices

Base.getindex(a::Polyhedron,i::Integer) = a.vertices[i]

get_graph(a::Polyhedron) = a.edge_vertex_graph

get_tags(a::Polyhedron) = a.vertex_tags

function signed_distance(point::Point{D},Π::CartesianPlane{D}) where D
  Π.orientation ? point[Π.d] - Π.value : Π.value - point[Π.d]
end

function split(p::Polyhedron,Π,plane_id=UNSET)
  smin = Inf
  smax = -Inf
  distances = zeros(num_vertices(p))
  new_vertices = empty( get_vertex_coordinates(p) )
  for v in 1:num_vertices(p)
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
  new_tags = empty( get_tags(p) )
  D = num_dims(p)
  for v in 1:num_vertices(p)
    distances[v] > 0 && continue
    for (i,vneig) in enumerate( get_graph(p)[v] )
      vneig != UNSET || continue
      distances[vneig] > 0 || continue
      vertex = compute_intersection(p[v],distances[v],p[vneig],distances[vneig])
      push!( new_vertices, vertex )
      tags = intersect( get_tags(p)[v], get_tags(p)[vneig] )
      push!(tags,plane_id)
      push!(new_tags,tags)

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
  compact_graph!(in_graph,distances,true)
  compact_graph!(out_graph,distances,false)
  vertices = [p.vertices;new_vertices]
  tags = [ get_tags(p); new_tags ] 
  Polyhedron(vertices,in_graph,tags), Polyhedron(vertices,out_graph,tags)
end

function complete_graph!(edge_graph,num_vertices::Integer)
  for v in num_vertices+1:length(edge_graph)
    vnext = next_vertex(edge_graph,num_vertices,v)
    vnext ≠ UNSET || continue
    edge_graph[v][end] = vnext
    edge_graph[vnext][2] = v
  end
end

function compact_graph!(edge_graph,distances,mask::Bool)
  for i in 1:length(distances)
    if mask == ( distances[i] > 0 )
      edge_graph[i] = empty!( edge_graph[i] )
    end
  end
end

function next_vertex(edge_graph::Vector,num_vertices::Integer,vstart::Integer)
  vcurrent = vstart
  vnext = first( edge_graph[vcurrent] )
  while vnext ≤ num_vertices && vnext ≠ UNSET
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

function simplexify(poly::Polyhedron)
  @notimplemented
  #TODO: use touched paths arrays for find facets this is in r3d: reduce()
  # the polyhedron must be convex
end

function refine(poly::Polyhedron,tree_id::Integer,Π::AbstractVector{<:Plane})
  Πi = first(Π)
  p⁻,p⁺ = split(poly,Πi)
  T = [p⁺,p⁻]
  tree_ids = [tree_id << 1, tree_id << 1 | 1 ]
  deleteat!(tree_ids,T.===nothing)
  deleteat!(T,T.===nothing)
  if length(Π) == 1
    R,rids = T,tree_ids
  else
    R,rids = empty(T),empty(tree_ids)
    Πn = view(Π,2:length(Π))
    for (k,id) in zip(T,tree_ids)
      T,ids = refine(k,id,Πn)
      append!(R,T)
      append!(rids,ids)
    end
  end
  return R,rids
end

function clip(poly::Polyhedron,Π::AbstractVector)
  p = poly
  for (i,Πi) ∈ enumerate(Π)
    p,_ = split(p,Πi,i)
    p !== nothing || break
  end
  p
end

function compute_graph(face_to_vertices)
  num_vertices = maximum( maximum(v) for v in face_to_vertices )
  graph = [ Int[] for _ in 1:num_vertices ]
  for face_vertices in face_to_vertices
    for (i,v) in enumerate(face_vertices)
      v⁺ = face_vertices[ i < length(face_vertices) ? i+1 : 1 ]
      v⁻ = face_vertices[ i > 1 ? i-1 : length(face_vertices) ]
      if v⁺ ∉ graph[v]
        push!(graph[v],v⁺)
      end
      if v⁻ ∉ graph[v]
        push!(graph[v],v⁻)
      end
    end
  end
  for v in 1:num_vertices
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

function next_vertex(p::Polyhedron,vprevious::Integer,vcurrent::Integer)
  i = findfirst( isequal(vprevious), get_graph(p)[vcurrent] )
  inext = ( i % length( get_graph(p)[vcurrent] ) ) + 1
  get_graph(p)[vcurrent][ inext ]
end

function orient_edge_endpoints(
  p::Polytope,
  cell::Polyhedron,
  edge::Integer,
  poly::Polyhedron,
  vertex::Integer)

  vnext = next_vertex(poly,UNSET,vertex)
  v_to_f = get_tags(poly)
  f = v_to_f[vertex][ findfirst( in(v_to_f[vnext]), v_to_f[vertex] ) ] # f[next] ∩ f[v]
  e = get_faces(p,1,0)[edge]
  vnext = next_vertex(cell,e[1],e[2])

  vnext ∈ get_face_vertices(p,num_dims(p)-1)[f] ? (e[1],e[2]) : (e[2],e[1])
end

function merge(p::Polytope,cell::Polyhedron,poly::Polyhedron)
  D = num_dims(p)
  v_to_df = fill(UNSET,num_vertices(poly))
  for v in 1:num_vertices(poly)
    d = D - length( get_tags(poly)[v] )
    for f in get_tags(poly)[v], df in get_faces(HEX,D-1,d)[f]
      if get_faces(HEX,d,D-1)[df] ⊆ get_tags(poly)[v]
        v_to_df[v] = get_dimrange(HEX,d)[df]
        break
      end
    end
  end
  graph = [ deepcopy(get_graph(cell)); map( i -> i.+num_vertices(HEX), get_graph(poly)) ]
  vertices = [ cell.vertices; poly.vertices ]
  e_to_v = get_faces(HEX,1,0)
  for v in 1:num_vertices(poly)
    v_to_df[v] != UNSET || continue
    if get_facedims(HEX)[v_to_df[v]] == 1
      e = v_to_df[v] - get_offset(HEX,1)
      
      v_in,v_out = orient_edge_endpoints(HEX,cell,e,poly,v)

      # TODO: If more than one vertex intersects the cell-edge:
      #       Join to the previous/next in the list.

      i = findfirst( isequal(UNSET), get_graph(poly)[v] )
      i_in = findfirst( isequal(v_out), get_graph(cell)[v_in] )

      graph[v_in][i_in] = v + num_vertices(HEX)
      graph[v+num_vertices(HEX)][i] = v_in
    end
  end
  vtouch = fill(false,length(graph))
  stack = Int[]
  for v in num_vertices(HEX)+1:length(graph)
    vtouch[v] = true
    empty!(stack)
    push!(stack,v)
    while length(stack) > 0
      vcurrent = pop!(stack)
      for vnext in graph[vcurrent]
        if !vtouch[vnext]
          push!(stack,vnext)
          vtouch[vnext] = true
        end
      end
    end
  end
  ## Unlink unused (OUT) vertices, we should compact
  for v in 1:num_vertices(HEX)
    if !vtouch[v]
      empty!(graph[v])
    end
  end
  Polyhedron(vertices,graph)
end

function get_cell_planes(p::Polytope,pmin::Point,pmax::Point)
  @notimplementedif !is_n_cube(p)
  D = num_dims(p)
  [ CartesianPlane(pmin+iseven(i)*pmax,D-((i-1)÷2),iseven(i)) for i in 1:2*D ]
end

## Driver

poly = Polyhedron(HEX)

n1 = VectorValue(1.0,0.7,0.5)
n2 = VectorValue(0.0,0.0,-1.0)
o = Point(0.5,0.5,0.5)

Π1 = Plane(o,n1)
Π2 = Plane(o,n2)


p⁻,p⁺ = split(poly,Π1)

p⁻⁻,p⁻⁺ = split(p⁻,Π2)

G⁻ = edge_mesh(p⁻⁻)
G⁺ = edge_mesh(p⁺)
writevtk(G⁻,"poly-")
writevtk(G⁺,"poly+")

T,ids = refine(poly,1,[Π1,Π2])

poly = clip(poly,[Π1,Π2])

G = edge_mesh(poly)
writevtk(G,"poly")


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

𝓖1 = compute_graph( facet_to_vertices )

facet_to_vertices = 
[[4,6,1],
 [4,3,6],
 [1,6,2],
 [3,5,6],
 [6,7,2] ]

𝓖2 = compute_graph( facet_to_vertices )
@show 𝓖1

@show 𝓖1[6] == 𝓖2[6]


## Start of the Algorithm

𝓖0 = [
  [4, 6, 2, 0], 
  [1, 6, 7, 0], 
  [6, 4, 0, 5], 
  [6, 1, 0, 3], 
  [7, 6, 3, 0], 
  [1, 4, 3, 5, 7, 2], 
  [6, 5, 0, 2]]

p0 = Polyhedron(vertices,𝓖0)

cell = Polyhedron(HEX)

Πc = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))

p = clip(p0,Πc)

p = merge(HEX,cell,p)

writevtk(edge_mesh(p0),"original_skin")

writevtk(edge_mesh(p),"in_graph")

writevtk( edge_mesh(cell),"cell_graph" )

end # module




