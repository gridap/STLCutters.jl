
function distribute_faces(
  cell_to_nodes,
  node_to_coordinates::Vector{<:Point},
  p::Polytope,
  faces::AbstractVector,
  STL_faces)

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

function distance_to_boundary(cell_nodes,node_to_coordinates,p::Polytope,point::Point)
  @assert have_intersection(cell_nodes,node_to_coordinates,p,point)
  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates,p)
  min( minimum(point-pmin), minimum(pmax-point) )
end

function farthest_axis_from_boundary(cell_nodes,node_to_coordinates,p::Polytope,point::Point)
  @assert have_intersection(cell_nodes,node_to_coordinates,p,point)
  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates,p)
  max_dists = max( point-pmin, pmax-point )
  _,d = findmax(max_dists.data)
  d
end

