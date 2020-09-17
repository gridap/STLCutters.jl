
function compute_stl_topology(
  cell_to_vertices::Table,
  vertex_to_coordinates::Vector{<:Point})

  grid = compute_stl_grid(cell_to_vertices,vertex_to_coordinates)
  GridTopology(grid)
end

function compute_stl_grid(
  cell_to_vertices::Table,
  vertex_to_coordinates::Vector{<:Point{D}}) where D

  p = Polytope( tfill(TET_AXIS,Val{D-1}()) )
  reffe = LagrangianRefFE(p)
  reffes = [reffe]
  cell_types = fill(1,length(cell_to_vertices))
  UnstructuredGrid(vertex_to_coordinates,cell_to_vertices,reffes,cell_types)
end

function get_edge_coordinates(stl::GridTopology)
  Tp = eltype(get_vertex_coordinates(stl))
  T = eltype(Tp)
  Dp = length(Tp)
  edge_coordinates = Vector{Segment{Dp,T}}(undef,num_edges(stl))
  for edge in 1:num_edges(stl)
    edge_coordinates[edge] = get_edge_coordinates(stl,edge)
  end
  edge_coordinates
end

function get_facet_coordinates(stl::GridTopology{Dc,2}) where Dc
  get_edge_coordinates(stl)
end

function get_facet_coordinates(stl::GridTopology{Dc,3}) where Dc
  Tp = eltype(get_vertex_coordinates(stl))
  T = eltype(Tp)
  Dp = length(Tp)
  facet_coordinates = Vector{Triangle{Dp,T}}(undef,num_cells(stl))
  for facet in 1:num_cells(stl)
    facet_coordinates[facet] = get_facet_coordinates(stl,facet)
  end
  facet_coordinates
end

function get_edge_coordinates(stl::GridTopology,edge::Integer)
  edge_vertices = get_faces(stl,1,0)[edge]
  X = get_vertex_coordinates(stl)
  Segment( X[edge_vertices[1]], X[edge_vertices[2]] )
end

function get_facet_coordinates(stl::GridTopology{Dc,3},facet::Integer) where Dc
  facet_vertices = get_faces(stl,2,0)[facet]
  X = get_vertex_coordinates(stl)
  Triangle( X[facet_vertices[1]], X[facet_vertices[2]], X[facet_vertices[3]] )
end
