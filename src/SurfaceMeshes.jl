
struct SurfaceMesh{D,T<:AbstractFloat,M}
  vertex_coordinates::Vector{Point{D,T}}
  facet_normals::Vector{VectorValue{D,T}}
  mfaces_to_nfaces::Matrix{M}
  d_to_offset::Vector{Int}
  
  function SurfaceMesh(
    vertex_coordinates::Vector{Point{D,T}},
    facet_normals::Vector{VectorValue{D,T}},
    facet_to_vertices::M) where {D,T,M}

    mfaces_to_nfaces = _compute_mfaces_to_nfaces(facet_to_vertices,D) 

    d_to_offset = _compute_d_to_offset(mfaces_to_nfaces)
    new{D,T,M}(vertex_coordinates, facet_normals, mfaces_to_nfaces, d_to_offset)
  end

  function SurfaceMesh(
    vertex_coordinates::Vector{Point{D,T}},
    facet_normals::Vector{VectorValue{D,T}},
    d_to_dface_to_vertices::Vector{M}) where {D,T,M}

    mfaces_to_nfaces = _compute_mfaces_to_nfaces(d_to_dface_to_vertices,D) 

    d_to_offset = _compute_d_to_offset(mfaces_to_nfaces)
    new{D,T,M}(vertex_coordinates, facet_normals, mfaces_to_nfaces, d_to_offset)
  end

end

function SurfaceMesh(stl::STL)
  stl = delete_repeated_vertices(stl)
  SurfaceMesh( 
    get_vertex_coordinates(stl), 
    get_facet_normals(stl), 
    get_facet_to_vertices(stl) )
end

function SurfaceMesh(vertex_coordinates::Vector{<:Point},facet_to_vertices)
  facet_normals = _compute_facet_normals(vertex_coordinates,facet_to_vertices)
  SurfaceMesh(vertex_coordinates,facet_normals,facet_to_vertices)
end

function SurfaceMesh(vertex_coordinates::Vector{<:Point},d_to_dface_to_vertices::Vector)
  facet_to_vertices = d_to_dface_to_vertices[end]
  facet_normals = _compute_facet_normals(vertex_coordinates,facet_to_vertices)
  SurfaceMesh(vertex_coordinates,facet_normals,d_to_dface_to_vertices)
end

function _compute_facet_normals(vertex_coordinates::Vector{Point{D,T}},facet_to_vertices::M) where {D,T,M}
  facet_normals = Vector{VectorValue{D,Float64}}(undef,length(facet_to_vertices))
  for i in 1:length(facet_to_vertices)
    facet = _get_facet(vertex_coordinates,facet_to_vertices,i)
    facet_normals[i] = normal(facet)
  end
  facet_normals
end


function _get_facet(vertex_coordinates::Vector{<:Point{2}},facet_to_vertices,i::Integer) 
  v = vertex_coordinates
  f = facet_to_vertices
  Segment( v[f[i,1]], v[f[i,2]] )
end

function _get_facet(vertex_coordinates::Vector{<:Point{3}},facet_to_vertices,i::Integer) 
  v = vertex_coordinates
  f = facet_to_vertices
  Triangle( v[f[i,1]], v[f[i,2]], v[f[i,3]] )
end

function _get_facet(vertex_coordinates::Vector{<:Point{4}},facet_to_vertices,i::Integer) 
  v = vertex_coordinates
  f = facet_to_vertices
  Tetrahedron( v[f[i,1]], v[f[i,2]], v[f[i,3]], v[f[i,4]] )
end

function _compute_mfaces_to_nfaces(facet_to_vertices::Table,D::Int)
  T = eltype(facet_to_vertices)
  mf_to_nf = Matrix{Table{T}}(undef,D,D)
  n_dfaces = fill(UNSET,D)

  mf_to_nf[D,0+1] = facet_to_vertices
  n_dfaces[D] = length( mf_to_nf[D,0+1] )
  n_dfaces[0+1] = maximum( mf_to_nf[D,0+1] )

  mf_to_nf[0+1,D] = compute_nface_to_dfaces_dual( mf_to_nf[D,0+1],n_dfaces[0+1])
    
  for d in reverse(2:D-1)
    ldface_to_lvertex = _get_simplex_dface_to_lvertices(d,d-1)
    mf_to_nf[d+1,d], n_dfaces[d] = compute_nface_to_dfaces_primal(mf_to_nf[d+1,0+1],mf_to_nf[0+1,d+1],ldface_to_lvertex)
    mf_to_nf[d,0+1] = compute_dface_to_vertices( mf_to_nf[d+1,0+1], mf_to_nf[d+1,d],ldface_to_lvertex,n_dfaces[d])
    mf_to_nf[0+1,d] = compute_nface_to_dfaces_dual( mf_to_nf[d,0+1],n_dfaces[0+1])
    mf_to_nf[d,d+1] = compute_nface_to_dfaces_dual( mf_to_nf[d+1,d],n_dfaces[d])
  end

  for d in 0:D-1
    mf_to_nf[d+1,d+1] = compute_dface_to_dface(n_dfaces[d+1],T)
  end
  mf_to_nf
end

function _compute_mfaces_to_nfaces(d_to_dface_to_vertices::Vector{<:Table},D) 
  T = eltype( eltype( d_to_dface_to_vertices ) )
  mf_to_nf = Matrix{Table{T}}(undef,D,D)
  n_dfaces = fill(UNSET,D)

  for d in 0:D-1
    dface_to_vertices = d_to_dface_to_vertices[d+1]
    n_dfaces[d+1] = length( dface_to_vertices )
    mf_to_nf[d+1,0+1] = dface_to_vertices 
    mf_to_nf[0+1,d+1] = compute_nface_to_dfaces_dual( mf_to_nf[d+1,0+1],n_dfaces[0+1])
    mf_to_nf[d+1,d+1] = compute_dface_to_dface(n_dfaces[d+1],T)
  end

  for d in 2:D-1, n in 1:d-1
    ldface_to_lvertex = _get_simplex_dface_to_lvertices(d,n)
    mf_to_nf[d+1,n+1] = compute_nface_to_dfaces(mf_to_nf[d+1,0+1],mf_to_nf[0+1,n+1],ldface_to_lvertex)
    mf_to_nf[n+1,d+1] = compute_nface_to_dfaces_dual( mf_to_nf[d+1,n+1], length(mf_to_nf[n+1,0+1]) )
  end 

  mf_to_nf

end

function _compute_d_to_offset(mfaces_to_nfaces::Matrix)
  D = size(mfaces_to_nfaces,1)
  d_to_offset = zeros(Int,D+1)
  for d in 0:D-1
    d_to_offset[d+2] = d_to_offset[d+1] + length(mfaces_to_nfaces[d+1,d+1])
  end
  d_to_offset
end

function _get_simplex_dface_to_lvertices(D::Integer,d::Integer)
  D_to_d_to_dface_to_vertices_for_tet[D][d]
end
  
num_dims(s::SurfaceMesh{D}) where D = D

function get_faces(s::SurfaceMesh,d::Integer,n::Integer)
  s.mfaces_to_nfaces[d+1,n+1]
end

function get_dface_to_vertices(s::SurfaceMesh,d::Integer)
  get_faces(s,d,0)
end

get_vertex_coordinates(s::SurfaceMesh) = s.vertex_coordinates

function num_vertices(s::SurfaceMesh)
  length(s.vertex_coordinates)
end

function num_faces(s::SurfaceMesh,d::Integer)
  length( get_dface_to_vertices(s,d) )
end

num_dfaces(s::SurfaceMesh,d::Integer) = num_faces(s,d)

num_edges(s::SurfaceMesh) = num_faces(s,1)

num_facets(s::SurfaceMesh{D}) where D = num_faces(s,D-1)

function num_faces(s::SurfaceMesh)
  s.d_to_offset[end]
end

@inline function global_dface(s::SurfaceMesh,d::Integer,lid::Integer)
  s.d_to_offset[d+1]+lid
end

function face_dimension(s::SurfaceMesh{D},gid::Integer) where D
  for d in 0:D-1
    if s.d_to_offset[d+2] ≥ gid
      return d
    end
  end
  @check false
  return UNSET
end


function is_face_dimension(s::SurfaceMesh,gid::Integer,d::Integer)
  s.d_to_offset[d+1] < gid ≤ s.d_to_offset[d+2]
end

is_vertex(s::SurfaceMesh,gid) = is_face_dimension(s,gid,0)

is_edge(s::SurfaceMesh,gid) = is_face_dimension(s,gid,1)

is_facet(s::SurfaceMesh{D},gid) where D = is_face_dimension(s,gid,D-1)

function local_dface(s::SurfaceMesh,gid::Integer,d::Integer)
  @check is_face_dimension(s,gid,d)
  gid - s.d_to_offset[d+1]
end

local_vertex(s::SurfaceMesh,gid::Integer) = local_dface(s,gid,0) 

local_edge(s::SurfaceMesh,gid::Integer) = local_dface(s,gid,1) 

local_facet(s::SurfaceMesh{D},gid::Integer) where D = local_dface(s,gid,D-1) 

function get_vertex_coordinates(s::SurfaceMesh,i::Integer)
  get_face_coordinates(s,Val{0}(),i)
end

function get_edge_coordinates(s::SurfaceMesh,i::Integer)
  get_face_coordinates(s,Val{1}(),i)
end

function get_facet_coordinates(s::SurfaceMesh{D},i::Integer) where D
  get_face_coordinates(s,Val{D-1}(),i)
end

function get_face_coordinates(s::SurfaceMesh,::Val{d}, i::Integer) where d
  error("Not implemented for dimension = $d")
end

function get_face_coordinates(s::SurfaceMesh,::Val{0}, i::Integer)
  s.vertex_coordinates[i]
end

function get_face_coordinates(s::SurfaceMesh,::Val{1}, i::Integer)
  df_to_v = get_dface_to_vertices(s,1)
  v = get_vertex_coordinates(s)
  Segment(v[df_to_v[i,1]],v[df_to_v[i,2]])
end

function get_face_coordinates(s::SurfaceMesh,::Val{2}, i::Integer)
  df_to_v = get_dface_to_vertices(s,2)
  v = get_vertex_coordinates(s)
  Triangle( v[df_to_v[i,1]], v[df_to_v[i,2]], v[df_to_v[i,3]] )
end

get_facet_to_vertices(s::SurfaceMesh{D}) where D = get_dface_to_vertices(s,D-1)

get_facet_normals(s::SurfaceMesh) = s.facet_normals

get_facet_normal(s::SurfaceMesh,i::Integer) = s.facet_normals[i]

function flip_normals(s::SurfaceMesh)
  SurfaceMesh(
    get_vertex_coordinates(s),
    - get_facet_normals(s),
    get_facet_to_vertices(s) )
end

function surface(s::SurfaceMesh)
  surface = 0.0
  for i in 1:num_facets(s)
    facet = get_facet_coordinates(s,i)
    surface += measure(facet)
  end
  surface
end

function move!(s::SurfaceMesh,offsets::VectorValue)
  v = get_vertex_coordinates(s)
  for i in 1:num_vertices(s)
    v[i] = v[i] + offsets
  end
  s
end

function compute_nface_to_dfaces_dual(nface_to_dfaces::Table,n_dfaces::Integer)

  ptrs = zeros(Int32,n_dfaces+1)
  for nface in 1:length(nface_to_dfaces)
    for ldface in 1:length(nface_to_dfaces,nface)
      dface = nface_to_dfaces[nface,ldface]
      ptrs[dface+1] += 1
    end
  end

  length_to_ptrs!(ptrs)

  T = eltype(nface_to_dfaces)
  ndata = ptrs[end]-1
  data = zeros(T,ndata)

  for nface in 1:length(nface_to_dfaces)
    for ldface in 1:length(nface_to_dfaces,nface)
      dface = nface_to_dfaces[nface,ldface]
      p = ptrs[dface]
      data[p] = nface
      ptrs[dface] += 1
    end
  end

  rewind_ptrs!(ptrs)

  Table(data,ptrs)

end

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
  for vertex in 1:n_vertices
    nfaces_around = length(vertex_to_nfaces,vertex)
    max_nfaces_around = max(max_nfaces_around,nfaces_around)
  end
  nfaces_around = fill(T(UNSET),max_nfaces_around)
  nfaces_around_bis = fill(T(UNSET),max_nfaces_around)
  n_nfaces_around = UNSET
  n_nfaces_around_bis = UNSET


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

function compute_dface_to_dface(n_dfaces,::Type{T}) where T
  ptrs = fill(Int32(1),n_dfaces+1)
  length_to_ptrs!(ptrs)
  data = [ T(i) for i in 1:n_dfaces ]
  Table(data,ptrs)
end

function writevtk(sm::SurfaceMesh{D,T},file_base_name) where {D,T}
  d_to_vtk_type_id = Dict(0=>1,1=>3,2=>5,3=>10)
  num_points = num_vertices(sm)
  points = zeros(T,D,num_points)
  for (i ,p ) in enumerate( get_vertex_coordinates(sm) ), d in 1:D
    points[d,i] = p[d]
  end
  cells = MeshCell{Vector{Int64}}[]
  for d in 0:D-1
    dface_to_vertices = get_dface_to_vertices(sm,d)
    vtk_type = VTKCellType(d_to_vtk_type_id[d])
    for i in 1:num_faces(sm,d)
      vertices = [ dface_to_vertices[i,j] for j in 1:vtk_type.nodes ]
      push!( cells, MeshCell(vtk_type,vertices) )
    end
  end
  normals = zeros(3,num_faces(sm))
  offset = num_faces(sm) - num_facets(sm)
  for i in 1:num_facets(sm)
    facet_normal = get_facet_normal(sm,i)
    for d in 1:D
       normals[d,i+offset] = facet_normal[d]
    end
  end
  vtkfile = vtk_grid(file_base_name,points,cells)
  vtkfile["facet_normals",VTKCellData()] = normals
  vtk_save(vtkfile)
end

function have_intersection(bb::BoundingBox,sm::SurfaceMesh,gid::Integer)
  d = face_dimension(sm,gid)
  lid = local_dface(sm,gid,d)
  have_intersection(bb,sm,d,lid)
end

function BoundingBox(sm::SurfaceMesh,gid::Integer)
  d = face_dimension(sm,gid)
  lid = local_dface(sm,gid,d)
  BoundingBox(sm,d,lid)
end

BoundingBox(s::SurfaceMesh) = BoundingBox(s.vertex_coordinates)

@generated function BoundingBox(s::SurfaceMesh{D},d::Integer,i::Integer) where D
  body = ""
  for d in 0:D-1
    body *= "if d == $d \n"
    body *= "  f$d = get_face_coordinates(s,Val{$d}(),i) \n"
    body *= "  BoundingBox(f$d) \n"
    body *= "else"
  end
  error = "throw(ArgumentError(\"\$d-face does not exist\"))"
  str = "$body \n  $error \nend"
  Meta.parse(str)
end

@generated function have_intersection(b::BoundingBox{D},s::SurfaceMesh{D},d::Integer,i::Integer) where D
  body = ""
  for d in 0:D-1
    body *= "if d == $d \n"
    body *= "  f$d = get_face_coordinates(s,Val{$d}(),i) \n"
    body *= "  have_intersection(f$d,b) \n"
    body *= "else"
  end
  error = "throw(ArgumentError(\" \$d-face does not exist\"))"
  str = "$body \n  $error \nend"
  Meta.parse(str)
end

function is_nface_in_dface(sm::SurfaceMesh,n::Integer,nface::Integer,d::Integer,dface::Integer)
  @check n ≤ d

  dface_to_nfaces = get_faces(sm,d,n)
  for lnface in 1:length(dface_to_nfaces,dface)
    _nface = dface_to_nfaces[dface,lnface]
    if _nface == nface
      return true
    end
  end
  false
end

function is_watter_tight(s::SurfaceMesh{D}) where D
  lface_to_facets = get_faces(s,D-2,D-1)
  for i in 1:length(lface_to_facets)
    if length(lface_to_facets,i) != 2
      return false
    end
  end
  true
end
