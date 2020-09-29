
function read_stl(filename::String)
  mesh = load(filename)
  vertex_to_coordinates = MeshIO.decompose(MeshIO.Point{3,Float64}, mesh )
  facet_to_vertices = MeshIO.decompose(MeshIO.Face{3,Int},mesh)
  vertex_to_normals = MeshIO.decompose(MeshIO.Normal{3,Float64},mesh)
  vertex_to_coordinates = Vector{Point{3,Float64}}( vertex_to_coordinates )
  facet_to_vertices = Table( Vector{Vector{Int}}( facet_to_vertices ) )
  vertex_to_normals = Vector{VectorValue{3,Float64}}( vertex_to_normals )
  facet_to_normals = vertex_to_normals[1:3:end]
  vertex_to_coordinates,facet_to_vertices,facet_to_normals
end

function Base.convert(::Type{Point{D,T}},x::MeshIO.Point{D}) where {D,T} 
  Point{D,T}(x.data)
end

function Base.convert(::Type{VectorValue{D,T}},x::MeshIO.Normal{D})  where {D,T}
  VectorValue{D,T}(x.data)
end

Base.convert(::Type{Vector},x::MeshIO.Face) = collect(x.data)

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

function merge_nodes(stl::Grid)
  X,T = delete_repeated_vertices(stl)
  compute_stl_grid(T,X)
end

function delete_repeated_vertices(stl::Grid)
  group_to_vertices =  _group_vertices(stl)
  vertices_map = collect(1:num_nodes(stl))
  for _vertices in group_to_vertices
    for i in 1:length(_vertices), j in i+1:length(_vertices)
      vi = get_node_coordinates(stl)[_vertices[i]]
      vj = get_node_coordinates(stl)[_vertices[j]]
      if vertices_map[_vertices[i]] == _vertices[i] && 
         vertices_map[_vertices[j]] == _vertices[j]

        if distance(vi,vj) < TOL
          vertices_map[_vertices[j]] = _vertices[i]
        end
      end
    end
  end
  u = unique(vertices_map)
  m = vertices_map .== 1:length(vertices_map)
  m = cumsum(m) .* m
  vertices_map = m[vertices_map]
  X = get_node_coordinates(stl)[u]
  T = get_cell_nodes(stl)
  T = Table( map( i -> vertices_map[i], T ) )
  X,T
end

function _group_vertices(stl::Grid{Dc,D}) where {Dc,D}
  min_length = _compute_min_length(stl)
  num_digits = -Int(floor(log10(min_length)))
  pmin, pmax = get_bouding_box(stl)
  cells = Int[]
  vertices = Int[]
  ranks = ceil.(Tuple(pmax-pmin),digits=num_digits)
  ranks = Int.(exp10(num_digits).*ranks).+1
  lids = LinearIndices( ranks )
  m = CartesianIndices( tfill(2,Val{D}()) )
  for (iv,v) in enumerate(get_node_coordinates(stl))
    p = v-pmin
    for i in m
      f = d -> i.I[d] == 1 ? floor(p[d],digits=num_digits) : ceil(p[d],digits=num_digits)
      r = ntuple( f, Val{D}() )
      r = Int.(exp10(num_digits).*r).+1
      cell = lids[r...]
      push!(cells,cell)
      push!(vertices,iv)
    end
  end
  u = unique(cells)
  dict = Dict( zip(u,1:length(u)) )
  cells = map( i -> dict[i], cells )
  cell_to_vertices = [ Int[] for _ in 1:maximum(cells) ]
  for (cell,vertex) in zip(cells,vertices)
    push!(cell_to_vertices[cell],vertex)
  end
  cell_to_vertices
end

function get_bouding_box(mesh::Grid)
  pmin = get_node_coordinates(mesh)[1]
  pmax = get_node_coordinates(mesh)[1]
  for vertex in get_node_coordinates(mesh)
    pmin = Point(min.(Tuple(pmin),Tuple(vertex)))
    pmax = Point(max.(Tuple(pmax),Tuple(vertex)))
  end
  pmin,pmax
end

function _compute_min_length(mesh::Grid)
  min_length = Inf
  for nodes in get_cell_nodes(mesh)
    for i in 1:length(nodes), j in i+1:length(nodes)
      vi = get_node_coordinates(mesh)[nodes[i]]
      vj = get_node_coordinates(mesh)[nodes[j]]
      dist = distance(vi,vj)
      if dist < min_length
        min_length = dist
      end
    end
  end
  min_length
end

function is_water_tight(top::GridTopology{Dc}) where Dc
  l = length.( get_faces(top,Dc-1,Dc) )
  maximum(l) == minimum(l) == 2
end
