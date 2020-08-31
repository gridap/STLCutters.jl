
function vertex_refinement(
  cell_nodes,
  node_to_coordinates::Vector{<:Point},
  point::Point{D}) where D

  p = Polytope(tfill(HEX_AXIS,Val{D}()))
  reffe = LagrangianRefFE(Float64,p,2)
  new_cells = compute_new_cells(cell_nodes,node_to_coordinates,reffe)
  new_vertices = compute_new_vertices(cell_nodes,node_to_coordinates,reffe,point)
  new_cells, new_vertices
end

## Helpers

function compute_new_cells(
  cell_nodes::Vector{<:Integer},
  node_to_coordinates::Vector{<:Point},
  reffe::ReferenceFE)

  grid, gface_to_rface = compute_linear_grid_and_facemap(reffe)
  num_nodes_per_cell = length(cell_nodes)
  num_nodes = length(node_to_coordinates)
  new_cells = Vector{Int}[]
  for lcell in 1:num_cells(grid)
    new_cell = fill(UNSET,num_nodes_per_cell)
    for lnode in 1:num_nodes_per_cell
      node = get_cell_nodes(grid)[lcell][lnode] 
      node = gface_to_rface[node]
      if node ≤ num_nodes_per_cell
        n = cell_nodes[node]
      else
        n = node - num_nodes_per_cell + num_nodes
      end
      new_cell[lnode] = n
    end
    push!(new_cells,new_cell)
  end
  new_cells
end

function compute_new_vertices(
  cell_nodes::Vector{<:Integer},
  node_to_coordinates::Vector{<:Point},
  reffe::ReferenceFE,
  point::Point{D}) where D

  p = get_polytope(reffe)
  num_nodes_per_cell = length(cell_nodes)
  new_node_to_coordinates = eltype(point)[]
  for node in num_nodes_per_cell+1:num_nodes(reffe)
    vertex = compute_vertex_coordinates(cell_nodes,node_to_coordinates,p,node,point)
    push!(node_to_coordinates,vertex)
  end
  new_node_to_coordinates
end

function compute_vertex_coordinates(
  cell_nodes,
  node_to_coordinates,
  p::Polytope{D},
  iface::Integer,
  point::Point{D}) where D

  nface = p.dface.nfaces[iface]
  dim = p.dface.dims[iface]
  node = cell_nodes[ get_faces(p)[iface][1] ]
  anchor = node_to_coordinates[ node ]
  extrusion = nface.extrusion
  vertex = anchor
  for d in 1:D
    if extrusion[d] == HEX_AXIS
      v = point[d]
      vertex = Base.setindex(vertex,v,d)
    end
  end
  vertex
end

function compute_linear_grid_and_facemap(reffe::LagrangianRefFE)
  grid = compute_linear_grid(reffe)
  desc = get_cartesian_descriptor(grid)
  model = CartesianDiscreteModel(desc)
  labels = get_face_labeling(model)
  grid_face_to_reffe_face = get_face_entity(labels)
  grid,grid_face_to_reffe_face
end

