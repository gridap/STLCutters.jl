
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

