
function distribute_faces(
  cell_to_nodes,
  node_to_coordinates::Vector{<:Point},
  p::Polytope,
  faces::AbstractVector,
  STL_faces,
  query=have_intersection::Function)

  cell_to_faces = Vector{Int}[]
  for (i,cell) in enumerate(cell_to_nodes)
    push!(cell_to_faces,[])
    for face in faces
      object = STL_faces[face]
      if query(cell,node_to_coordinates,p,object)
        push!(cell_to_faces[i],face)
      end
    end
  end
  cell_to_faces
end

function distribute_face_skeleton(
  cell_to_nodes,
  node_to_coordinates::Vector{<:Point},
  p::Polytope,
  faces::AbstractVector,
  f::Face)

  cell_to_faces = Vector{Int}[]
  for (i,cell) in enumerate(cell_to_nodes)
    push!(cell_to_faces,[])
    if is_on_cell_facet(cell,node_to_coordinates,p,f)
      for face in faces
        if is_on_cell_facet(cell,node_to_coordinates,p,f,face)
          push!(cell_to_faces[i],face)
        end
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

