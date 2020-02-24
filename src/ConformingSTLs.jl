
struct SurfaceMesh{D,T,M}
  vertex_coordinates::Vector{Point{D,T}}
  facet_normals::Vector{Point{D,T}}
  mfaces_to_nfaces::Matrix{M}
  d_to_offset::Vector{Int}
  
  function SurfaceMesh(
    vertex_coordinates::Vector{Point{D,T}},
    facet_normals::Vector{Point{D,T}},
    facet_to_vercies::AbstractTable) where {D,T}


    mfaces_to_nfaces = _compute_mfaces_to_nfaces(facet_to_vertices) 

    d_to_offset = _compute_d_to_offset(mfaces_to_nfaces)
    new{D,T,M}(vertex_coordinates, facet_normals, mfaces_to_nfaces, d_to_offset)
  end
end

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

  vertex_to_vertex = TableOfVectors( [ [i] for i in 1:length(vertex_coordinates) ] )
  facet_to_facet = TableOfVectors( [ [i] for i in 1:length(facet_to_vertices) ] )

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

@inline function global_dface(stl::ConformingSTL,d::Int,lid::Int)
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

function get_vertex(stl::ConformingSTL,i::Int)
  get_dface(stl,Val{0}(),i)
end

function get_edge(stl::ConformingSTL,i::Int)
  get_dface(stl,Val{1}(),i)
end

function get_facet(stl::ConformingSTL{D},i::Int) where D
  get_dface(stl,Val{D-1}(),i)
end

function get_dface(smesh::ConformingSTL,::Val{d}, i::Int) where d
  error("Not implemented for dimension = $d")
end

function get_dface(smesh::ConformingSTL,::Val{0}, i::Int)
  stl.vertex_coordinates[i]
end

function get_dface(smesh::ConformingSTL,::Val{1}, i::Int)
  dface_to_vertices = get_dface_to_vertices(stl,1)
  l = getlist(dface_to_vertices,i)
  v = stl.vertex_coordinates
  Segment(v[l[1]],v[l[2]])
end

function get_dface(smesh::ConformingSTL,::Val{2}, i::Int)
  dface_to_vertices = get_dface_to_vertices(stl,2)
  l = getlist(dface_to_vertices,i)
  v = stl.vertex_coordinates
  Triangle(v[l[1]],v[l[2]],v[l[3]])
end

function compute_vertex_to_facets(facet_to_vertices::TableOfVectors{Int},num_vertices::Int)
  vertex_to_facets = [ Int[] for i in 1:num_vertices ]
  for i in 1:length(facet_to_vertices), j in 1:length(facet_to_vertices,i)
    push!(vertex_to_facets[ facet_to_vertices[i,j] ], i )
  end
  TableOfVectors(vertex_to_facets)
end

function compute_nface_to_dfaces_dual(nface_to_dfaces::Table,n_dfaces::Int)

  ptrs = zeros(Int32,n_dfaces+1)
  for nface in 1:length(nface_to_dfaces)
    for ldface in 1:length(nface_to_dfaces,face)
      dface = nface_to_dfaces[face,ldface]
      ptrs[dface+1] += 1
    end
  end

  length_to_ptrs!(ptrs)

  T = eltype(nface_to_dfaces)
  ndata = ptrs[end]-1
  data = zeros(T,ndata)

  for nface in 1:length(nface_to_dfaces)
    for ldface in 1:length(nface_to_dfaces,face)
      dface = nface_to_dfaces[face,ldface]
      p = ptrs[dface]
      data[p] = nface
      ptrs[dface] += 1
    end
  end

  rewind_ptrs!(ptrs)

  Table(data,ptrs)

end

const UNSET = 0

function compute_nface_to_dfaces_primal(
  nface_to_vertices::Table,
  vertex_to_nfaces::Table,
  ldface_to_lvertex)

  n_nfaces = length(nface_to_vertices)
  n_vertices = length(vertex_to_nfaces)
  n_ldfaces = length(ldface_to_lvertex)

  ptrs = fill(Int32(n_ldfaces),n_nfaces+1)
  length_to_ptrs!(ptrs)

  T = eltype(nface_to_vertices)
  ndata = ptrs[end]-1
  data = fill(T(UNSET),ndata)

  nface_to_dfaces  = Table(data,ptrs)

  dface = 0

  max_nfaces_around = 0
  for vertex in 1:n_nfaces
    nfaces_around = length(vertex_to_nfaces,vertex)
    max_nfaces_around = max(max_nfaces_around,nfaces_around)
  end
  nfaces_around = fill(T(UNSET),max_nfaces_around)
  nfaces_around_bis = fill(T(UNSET),max_nfaces_around)

  for nface in 1:n_nfaces
    for ldface in 1:n_ldfaces

      if nface_to_dfaces[nface,ldface] == UNSET
        dface += 1
        nface_to_dfaces[nface,ldface] = dface

        lvertices = ldface_to_lvertex[ldface]
        for (j,lvertex) in enumerate(lvertices)
          vertex = nface_to_vertices[nface,lvertex]
          if j == 1
            n_nfaces_around = length(vertex_to_nfaces,vertex)
            for i in 1:n_nfaces_around
              nface_around = vertex_to_nfaces[vertex,i]
              nfaces_around[i] = nface_around
            end
          else
            n_nfaces_around_bis = length(vertex_to_nfaces,vertex)
            for i in 1:n_nfaces_around_bis
              nface_around = vertex_to_nfaces[vertex,i]
              nfaces_around_bis[i] = nface_around
            end
            _set_intersection!(
              nfaces_around,
              nfaces_around_bis,
              n_nfaces_around,
              n_nfaces_around_bis)
          end
        end
        for i in 1:n_nfaces_around
          nface_around = nfaces_around[i]
          if (nface_around != UNSET) && (nface_around != nface)
            ldface_around = _find_nface_around_ldface(
              nface, ldface, nface_around, nface_to_vertices, ldface_to_lvertex, n_ldfaces)
            nface_to_dfaces[nface_around,ldface_around] = dface
          end
        end

      end

    end
  end

  (nface_to_dfaces, dface)

end


@inline function _find_nface_around_ldface(
  nface,
  ldface,
  nface_around,
  nface_to_vertices,
  ldface_to_lvertex,
  n_ldfaces)

  lvertices = ldface_to_lvertex[ldface]

  for ldface_around in 1:n_ldfaces

      all_found = true
      lvertices_around = ldface_to_lvertex[ldface_around]
      for lvertex_around in lvertices_around
        vertex_around = nface_to_vertices[nface_around,lvertex_around]
        found = false
        for lvertex in lvertices
          vertex = nface_to_vertices[nface,lvertex]
          if vertex == vertex_around
            found = true
            break
          end
        end
        if !found
          all_found = false
          break
        end
      end

      if all_found
        return ldface_around
      end
  end

  @check false

  return UNSET

end


@inline function _set_intersection!(
  nfaces_around,
  nfaces_around_bis,
  n_nfaces_around,
  n_nfaces_around_bis)

  for i in 1:n_nfaces_around
    if nfaces_around[i] == UNSET
      continue
    end
    _find_eq!(i,nfaces_around,nfaces_around_bis,n_nfaces_around_bis)
  end
end

@inline function _find_eq!(i,nfaces_around,nfaces_around_bis,n_nfaces_around_bis)
  for j in 1:n_nfaces_around_bis
    if nfaces_around[i] == nfaces_around_bis[j]
      return
    end
  end
  nfaces_around[i] = UNSET
  return
end


function compute_dface_to_vertices(
  nface_to_vertices::Table,
  nface_to_dfaces::Table,
  ldface_to_lvertex,
  n_dfaces)

  n_nfaces = length(nface_to_vertices)
  n_ldfaces = length(ldface_to_lvertex)

  ptrs = fill(Int32(UNSET),n_dfaces+1)
  for nface in 1:n_nfaces
    for ldface in 1:n_ldfaces
      dface = nface_to_dfaces[nface,ldface]
      if ptrs[dface+1] != UNSET
        continue
      end
      lvertices = ldface_to_lvertex[ldface]
      ptrs[dface+1] = length(lvertices)
    end
  end

  length_to_ptrs!(ptrs)
  ndata = ptrs[end]-1
  T = eltype(nface_to_vertices)
  data = fill(T(UNSET),ndata)
  dface_to_vertices = Table(data,ptrs)

  for nface in 1:n_nfaces
    for ldface in 1:n_ldfaces
      dface = nface_to_dfaces[nface,ldface]
      if dface_to_vertices[dface,1] != UNSET
        continue
      end
      lvertices = ldface_to_lvertex[ldface]
      for (i,lvertex) in enumerate(lvertices)
        vertex = nface_to_vertices[nface,lvertex]
        dface_to_vertices[dface,i] = vertex
      end
    end
  end

  dface_to_vertices

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

