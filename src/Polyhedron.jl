module PolyhedraTests

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using Gridap.Helpers

using STLCutters

using STLCutters: Plane
using STLCutters: signed_distance 

import Gridap: num_dims
import Gridap: num_vertices
import Gridap.ReferenceFEs: get_vertex_coordinates

struct Polyhedron{D,T}
  vertices::Vector{Point{D,T}}
  edge_vertex_graph::Vector{Vector{Int}}
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

num_dims(::Polyhedron{D}) where D = D

num_vertices(a::Polyhedron) = length(a.vertices)

get_vertex_coordinates(a::Polyhedron) = a.vertices

Base.getindex(a::Polyhedron,i::Integer) = a.vertices[i]

get_graph(a::Polyhedron) = a.edge_vertex_graph

function split(p::Polyhedron,Π::Plane)
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
  D = num_dims(p)
  for v in 1:num_vertices(p)
    distances[v] > 0 && continue
    for (i,vneig) in enumerate( get_graph(p)[v] )
      distances[vneig] > 0 || continue
      vertex = compute_intersection(p[v],distances[v],p[vneig],distances[vneig])
      push!( new_vertices, vertex )

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
  Polyhedron(vertices,in_graph), Polyhedron(vertices,out_graph)
end

function complete_graph!(edge_graph,num_vertices::Integer)
  for v in num_vertices+1:length(edge_graph)
    vnext = next_vertex(edge_graph,num_vertices,v)
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

function next_vertex(edge_graph,num_vertices::Integer,vstart::Integer)
  vcurrent = vstart
  vnext = first( edge_graph[vcurrent] )
  while vnext ≤ num_vertices
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
  #TODO: use touched paths arrays for find facets
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

function clip(poly::Polyhedron,Π::AbstractVector{<:Plane})
  p = poly
  for Πi ∈ Π
    p,_ = split(p,Πi)
    p !== nothing || break
  end
  p
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


end # module




