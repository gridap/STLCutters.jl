
function distribute_faces(
  cell_to_nodes,
  node_to_coordinates::Vector{<:Point{D}},
  faces::AbstractVector,
  STL_faces) where D

  p = Polytope(tfill(HEX_AXIS,Val{D}()))
  cell_to_faces = Vector{Int}[]
  for (i,cell) in enumerate(cell_to_nodes)
    push!(cell_to_faces,[])
    for face in faces
      object = STL_faces[face]
      if have_intersection(cell,node_to_coordinates,p,object)
        push!(cell_to_faces[i],face)
      end
    end
  end
  cell_to_faces
end

## Helpers

function distance_to_boundary(cell_nodes,node_to_coordinates,p::Point)
  @assert have_intersection(cell_nodes,node_to_coordinates,p)
  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates)
  min( minimum(p-pmin), minimum(pmax-p) )
end

function farthest_axis_from_boundary(cell_nodes,node_to_coordinates,p::Point)
  @assert have_intersection(cell_nodes,node_to_coordinates,p)
  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates)
  max_dists = max( p-pmin, pmax-p )
  _,d = findmax(max_dists.data)
  d
end

