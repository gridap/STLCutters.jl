
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

@inline function num_vertices(stl::ConformingSTL)
  length(stl.vertex_coordinates)
end

@inline function num_dfaces(stl::ConformingSTL,d::Int)
  length(stl.d_face_to_vertices[d+1])
end

function num_dfaces(stl::ConformingSTL{D}) where D
  n = 0
  for d in 0:D-1
    n += num_dfaces(stl,d)
  end
  n
end

function global_dface(stl::ConformingSTL,d::Int,lid::Int)
  gid = 0
  for i in 0:d-1
    gid += num_dfaces(stl,i)
  end
  gid += lid
end

function local_dface(stl::ConformingSTL{D},gid::Int) where D
  lid = gid
  for d in 0:D-1
    if lid > num_dfaces(stl,d)
      lid -= num_dfaces(stl,d)
    else
      return (d,lid)
    end
  end
  (D-1,lid)
end

function get_vertex(stl::ConformingSTL{D,T},i::Int) where {D,T}
  stl.vertex_coordinates[i]
end

function get_edge(stl::ConformingSTL{D},i::Int) where D
  dface_to_vertices = get_dface_to_vertices(stl,1)
  l = getlist(dface_to_vertices,i)
  v = stl.vertex_coordinates
  Segment(v[l[1]],v[l[2]])
end


function get_facet(stl::ConformingSTL{D},i::Int) where D
  dface_to_vertices = get_dface_to_vertices(stl,2)
  l = getlist(dface_to_vertices,i)
  v = stl.vertex_coordinates
  Triangle(v[l[1]],v[l[2]],v[l[3]])
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

function writevtk(stl::ConformingSTL{D,T},file_base_name) where {D,T}
  d_to_vtk_type_id = Dict(0=>1,1=>3,2=>5)
  num_points = num_vertices(stl)
  points = zeros(T,D,num_points)
  for (i ,p ) ∈ enumerate(stl.vertex_coordinates), d ∈ 1:D
    points[d,i] = p[d]
  end
  cells = MeshCell{Vector{Int64}}[]
  for d ∈ 0:D-1
    dface_to_vertices = get_dface_to_vertices(stl,d)
    num_dfaces = length(dface_to_vertices)
    vtk_type = VTKCellType(d_to_vtk_type_id[d])
    for i in 1:num_dfaces
      vertices = getlist(dface_to_vertices,i)
      push!( cells, MeshCell(vtk_type,vertices) )
    end
  end
  vtkfile = vtk_grid("out",points,cells)
  vtk_save(vtkfile)
end

function have_intersection(hex::HexaCell{D},stl::ConformingSTL{D},d::Int,i::Int) where D
  if d == 0
    p = get_vertex(stl,i)
    have_intersection(p,hex)
  elseif d == 1
    e = get_edge(stl,i)
    have_intersection(e,hex)
  elseif d == 2
    f = get_facet(stl,i)
    have_intersection(f,hex)
  else
    throw(ArgumentError("$d-face does not exist"))
  end
end

@inline have_intersection(hex::HexaCell,stl::ConformingSTL,gid::Int) = have_intersection(hex,stl,local_dface(stl,gid)...)

function BoundingBox(stl::ConformingSTL{D},d::Int,i::Int) where D
  if d == 0
    p = get_vertex(stl,i)
    BoundingBox(p)
  elseif d == 1
    e = get_edge(stl,i)
    BoundingBox(e)
  elseif d == 2
    f = get_facet(stl,i)
    BoundingBox(f)
  else
    throw(ArgumentError("$d-face does not exist"))
  end
end

@inline BoundingBox(stl::ConformingSTL,gid::Int) = BoundingBox(stl,local_dface(stl,gid)...)

function BoundingBox(stl::ConformingSTL)
  pmin = stl.vertex_coordinates[1]
  pmax = stl.vertex_coordinates[1]
  for v ∈ stl.vertex_coordinates
    pmin = min.(pmin,v)
    pmax = max.(pmax,v)
  end
  BoundingBox(pmin,pmax)
end
