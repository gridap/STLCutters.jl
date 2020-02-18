
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
  vertex_to_facets = TableOfVectors(Int,num_vertices,0)
  for i in 1:length(facet_to_vertices)
    list = getlist(facet_to_vertices,i)
    for j in 1:length(list)
      push_to_list!(vertex_to_facets, list[j], i )
    end
  end
  vertex_to_facets
end

const ledge_to_facet_vertices = ((1,2),(1,3),(2,3))
function get_local_edges( v::Vector{Int}, ledge )
  lvertices = ledge_to_facet_vertices[ledge]
  (v[lvertices[1]],v[lvertices[2]])
end

const num_local_edges = 3
function compute_edge_neighbors(facet_to_vertices::TableOfVectors{Int},vertex_to_facets::TableOfVectors{Int})
  # methods
  function set_intersection!(a::Vector{T},b::Vector{T}) where T
    i = 1
    while i <= length(a)
      ai = a[i]
      found = false
      for bi in b
        if ai == bi 
          found = true
          break
        end
      end
      if !found 
        deleteat!(a,i)
      else
        i += 1
      end
    end
    a
  end
  function set_difference!(a::Vector{T},b::T) where T
    i = 1
    while i <= length(a)
      if a[i] == b
        deleteat!(a,i)
      else
        i += 1
      end
    end
    a
  end
  # body 
  num_facets = length(facet_to_vertices)
  facet_to_edge_neighbors = TableOfVectors(Int,num_facets,num_local_edges)
  neighbor = Int[]
  for ifacet in 1:num_facets
    facet = getlist(facet_to_vertices,ifacet)
    for ledge in 1:num_local_edges
      edge = get_local_edges(facet,ledge)
      f1 = getlist(vertex_to_facets,edge[1])
      f2 = getlist(vertex_to_facets,edge[2])
      copy!(neighbor,f1)
      set_intersection!(neighbor,f2)
      set_difference!(neighbor,ifacet)
      @assert length(neighbor) == 1
      set_to_list!(facet_to_edge_neighbors,ifacet,ledge,neighbor[1])
    end
  end
  facet_to_edge_neighbors

end

const num_edges_per_facet = 3
function compute_facet_to_edges(facet_to_edge_neighbors::TableOfVectors{Int})
  num_facets = length(facet_to_edge_neighbors)
  facet_to_edges = TableOfVectors(Int,num_facets,num_edges_per_facet)
  num_edges = 0
  for ifacet in 1:length(facet_to_edge_neighbors)
    edges = getlist(facet_to_edges,ifacet)
    neighbors = getlist(facet_to_edge_neighbors,ifacet)
    for ledge in 1:num_edges_per_facet
      if ( edges[ledge] == 0 )
        num_edges += 1
        set_to_list!(facet_to_edges, ifacet, ledge, num_edges )
        ledge_at_neighbor = findfirst( x->x == ifacet, getlist(facet_to_edge_neighbors,neighbors[ledge]) )
        set_to_list!(facet_to_edges, neighbors[ledge], ledge_at_neighbor, num_edges )
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
        pushlist!(edge_to_vertices,collect(edge))
      end
    end
  end
  edge_to_vertices
end

function compute_edge_to_facets(facet_to_edges::TableOfVectors,num_edges::Int)
  edge_to_facets = TableOfVectors(Int,num_edges,0)
  for i in 1:length(facet_to_edges)
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
  for (i ,p ) in enumerate(stl.vertex_coordinates), d in 1:D
    points[d,i] = p[d]
  end
  cells = MeshCell{Vector{Int64}}[]
  for d in 0:D-1
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

function have_intersection(bb::BoundingBox{D},stl::ConformingSTL{D},d::Int,i::Int) where D
  if d == 0
    p = get_vertex(stl,i)
    have_intersection(p,bb)
  elseif d == 1
    e = get_edge(stl,i)
    have_intersection(e,bb)
  elseif d == 2
    f = get_facet(stl,i)
    have_intersection(f,bb)
  else
    throw(ArgumentError("$d-face not implemented"))
  end
end

function have_intersection(bb::BoundingBox,stl::ConformingSTL,gid::Int)
  have_intersection(bb,stl,local_dface(stl,gid)...)
end

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
