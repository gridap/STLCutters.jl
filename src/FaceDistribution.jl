
function distribute_faces(
  cell_to_nodes,
  node_to_coordinates::Vector{<:Point},
  p::Polytope,
  stl::DiscreteModel,
  dfaces::AbstractVector,
  ::Val{d}) where d

  cell_to_dfaces = Vector{Int}[]
  for (i,cell) in enumerate(cell_to_nodes)
    push!(cell_to_dfaces,[])
    for dface in dfaces
      object = get_dface(stl,dface,Val{d}())
      if have_intersection(cell,node_to_coordinates,p,object)
        push!(cell_to_dfaces[i],dface)
      end
    end
  end
  cell_to_dfaces
end

function distribute_vertices(
  cell_to_nodes,
  node_to_coordinates::Vector{<:Point},
  p::Polytope,
  stl::DiscreteModel,
  vertices::AbstractVector)

  distribute_faces(cell_to_nodes,node_to_coordinates,p,stl,vertices,Val{0}())
end

function distribute_edges(
  cell_to_nodes,
  node_to_coordinates::Vector{<:Point},
  p::Polytope,
  stl::DiscreteModel,
  edges::AbstractVector)

  distribute_faces(cell_to_nodes,node_to_coordinates,p,stl,edges,Val{1}())
end

function distribute_facets(
  cell_to_nodes,
  node_to_coordinates::Vector{<:Point},
  p::Polytope,
  stl::DiscreteModel{Dc,Dp},
  facets::AbstractVector) where {Dc,Dp}

  distribute_faces(cell_to_nodes,node_to_coordinates,p,stl,facets,Val{Dp-1}())
end

function distribute_faces(grid::Grid,stl::DiscreteModel,faces::AbstractVector)
  cell_to_faces = Vector{Int}[]
  for i in 1:num_cells(grid)
    push!(cell_to_faces,[])
    for face in faces
      if have_intersection(grid,i,stl,face)
        push!(cell_to_faces[i],face)
      end
    end
  end
  cell_to_faces
end
