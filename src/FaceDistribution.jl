
function distribute_vertices(
  cell_to_nodes,
  node_to_coordinates::Vector{<:Point},
  vertices::AbstractVector,
  STL_vertices)

  cell_to_vertices = Vector{Int}[]
  for (i,cell) in enumerate(cell_to_nodes)
    push!(cell_to_vertices,[])
    for vertex in vertices
      point = STL_vertices[vertex]
      if have_intersection(cell,node_to_coordinates,point)
        push!(cell_to_vertices[i],vertex)
      end
    end
  end
  cell_to_vertices
end


## Helpers

function Base.setindex(p::Point,v,idx::Integer)
  data = Base.setindex(p.data,v,idx)
  Point(data)
end

distance(a::Point,b::Point) =  norm( a - b )

function get_bounding_box(
  cell_nodes::Vector{<:Integer},
  node_to_coordinates::Vector{<:Point})

  pmin = node_to_coordinates[ first(cell_nodes) ]
  pmax = node_to_coordinates[ last(cell_nodes) ]
  pmin,pmax
end

function have_intersection(cell_nodes,node_to_coordinates,p::Point)
  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates)
  all( pmin.data .< p.data ) || return false
  all( pmax.data .> p.data ) || return false
  true
end

function distance_to_boundary(cell_nodes,node_to_coordinates,p::Point)
  @assert have_intersection(cell_nodes,node_to_coordinates,p)
  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates)
  min( minimum(p-pmin), minimum(pmax-p) )
end

