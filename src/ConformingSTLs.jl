
struct ConformingSTL{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  d_face_to_vertices::Vector{TableOfVectors{Int}}
  facet_normals::Vector{Point{D,T}}
  d_face_to_facets::Vector{TableOfVectors{Int}}
end

function ConformingSTL(filename::String)
  stl = RawSTL(filename)
  vertices_map = map_repeated_vertices(stl)
  vertex_coordinates = extract_unique_vertices(stl, vertices_map)
  vertices_map = compact_map(vertices_map)
  facet_to_vertices = apply_map(stl.facet_to_vertices, vertices_map)
  vertex_to_facets = compute_vertex_to_facets(facet_to_vertices,length(vertex_coordinates))
  facet_to_edge_neighbors = compute_edge_neighbors(facet_to_vertices,vertex_to_facets)
  facet_to_edges = compute_facet_to_edges(facet_to_edge_neighbors)
  edge_to_vertices = compute_edge_to_vertices(facet_to_vertices,facet_to_edges)
  edge_to_facets = compute_edge_to_facets(facet_to_edges,length(edge_to_vertices))

  vertex_to_vertex = TableOfVectors{Int}([])
  facet_to_facet = TableOfVectors{Int}([])
  for i = 1:length(vertex_coordinates)
    pushlist!(vertex_to_vertex,[i])
  end
  for i = 1:length(facet_to_vertices)
    pushlist!(facet_to_facet,[i])
  end

  d_face_to_vertices = [ vertex_to_vertex, edge_to_vertices,facet_to_vertices ]
  d_face_to_facets   = [ vertex_to_facets, edge_to_facets, facet_to_facet ]

  ConformingSTL( vertex_coordinates, d_face_to_vertices, stl.facet_normals, d_face_to_facets )
end

function get_dface_to_vertices(a::ConformingSTL,d::Integer)
  a.d_face_to_vertices[d+1]
end

function get_dface_to_facets(a::ConformingSTL,d::Integer)
  a.d_face_to_facets[d+1]
end

function compute_vertex_to_facets(facet_to_vertices::TableOfVectors{Int},num_vertices::Int)
  vertex_to_facets = TableOfVectors{Int}([])
  for i in 1:num_vertices
    pushlist!(vertex_to_facets, Vector{Int}([]))
  end
  for i in 1:length(facet_to_vertices)
    list = getlist(facet_to_vertices,i)
    for j in 1:length(list)
      push_to_list!(vertex_to_facets, list[j], i )
    end
  end
  vertex_to_facets
end

const edge_to_vertex = [[1,2],[1,3],[2,3]]
function get_local_edges( v::Vector{Int}, ledge )
  v[edge_to_vertex[ledge]]
end

const num_local_edges = 3
function compute_edge_neighbors(facet_to_vertices::TableOfVectors{Int},vertex_to_facets::TableOfVectors{Int})
  facet_to_edge_neighbors = TableOfVectors{Int}([])
  for i in 1:length(facet_to_vertices)
    pushlist!(facet_to_edge_neighbors,Vector{Int}(zeros(num_local_edges)))
  end
  for ifacet in 1:length(facet_to_vertices)
    facet = getlist(facet_to_vertices,ifacet)
    for ledge in 1:num_local_edges
      edge = get_local_edges(facet,ledge)
      @assert length(edge) == 2
      f1 = getlist(vertex_to_facets,edge[1])
      f2 = getlist(vertex_to_facets,edge[2])
      neighbor = f1[findall( in(f2),f1 )]
      neighbor = neighbor[findall(x->x != ifacet, neighbor)]
      @assert length(neighbor) == 1
      set_to_list!(facet_to_edge_neighbors,ifacet,ledge,neighbor[1])
    end
  end
  facet_to_edge_neighbors
end

function compute_facet_to_edges(facet_to_edge_neighbors::TableOfVectors{Int})
  facet_to_edges = TableOfVectors{Int}([])
  for i in 1:length(facet_to_edge_neighbors)
    pushlist!(facet_to_edges,Vector{Int}(zeros(3)))
  end
  num_edges = 0
  for ifacet in 1:length(facet_to_edge_neighbors)
    edges = getlist(facet_to_edges,ifacet)
    neighbors = getlist(facet_to_edge_neighbors,ifacet)
    for ledge in 1:num_local_edges
      if ( edges[ledge] == 0 )
        num_edges += 1
        set_to_list!(facet_to_edges, ifacet, ledge, num_edges )
        ledge_at_neighbor = findall( x->x == ifacet, getlist(facet_to_edge_neighbors,neighbors[ledge]) )
        @assert length(ledge_at_neighbor) == 1
        set_to_list!(facet_to_edges, neighbors[ledge], ledge_at_neighbor[1], num_edges )
      end
    end
  end
  facet_to_edges
end

function compute_edge_to_vertices(facet_to_vertices::TableOfVectors{Int},facet_to_edges::TableOfVectors{Int})
  edge_to_vertices = TableOfVectors{Int}([])
  num_edges = 0
  for ifacet in 1:length(facet_to_edges)
    facet = getlist(facet_to_vertices,ifacet)
    edges = getlist(facet_to_edges,ifacet)
    for ledge in 1:num_local_edges
      if ( edges[ledge] == num_edges + 1 )
        num_edges += 1
        edge = get_local_edges(facet,ledge)
        pushlist!(edge_to_vertices,edge)
      end
    end
  end
  edge_to_vertices
end

function compute_edge_to_facets(facet_to_edges::TableOfVectors,num_edges::Int)
  edge_to_facets = TableOfVectors{Int}([])
  for i in 1:num_edges
    pushlist!(edge_to_facets, Vector{Int}([]))
  end

  for i in 1:length(facet_to_edges)#, d = 1:num_dims(STL)
    list = getlist(facet_to_edges,i)
    for j in 1:length(list)
      push_to_list!(edge_to_facets, list[j], i )
    end
  end
  edge_to_facets
end
