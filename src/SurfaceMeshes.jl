
struct SurfaceMesh{D,T,M}
  vertex_coordinates::Vector{Point{D,T}}
  facet_normals::Vector{VectorValue{D,T}}
  mfaces_to_nfaces::Matrix{M}
  d_to_offset::Vector{Int}
  
  function SurfaceMesh(
    vertex_coordinates::Vector{Point{D,T}},
    facet_normals::Vector{VectorValue{D,T}},
    facet_to_vertices::M) where {D,T,M}


    mfaces_to_nfaces = _compute_mfaces_to_nfaces(facet_to_vertices,Val{D}()) 

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

function _compute_mfaces_to_nfaces(facet_to_vertices::Table,::Val{D}) where D
  T = eltype(facet_to_vertices)
  mf_to_nf = Matrix{Table{T}}(undef,D,D)
  n_dfaces = fill(UNSET,D)

  mf_to_nf[D,0+1] = facet_to_vertices
  n_dfaces[D] = length( mf_to_nf[D,0+1] )
  n_dfaces[0+1] = maximum( mf_to_nf[D,0+1] )

  mf_to_nf[0+1,D] = compute_nface_to_dfaces_dual( mf_to_nf[D,0+1],n_dfaces[0+1])
    
  for d in reverse(2:D-1)
    ldface_to_lvertex = tetn_ldface_to_lvertex[d][d-1]
    mf_to_nf[d+1,d], n_dfaces[d] = compute_nface_to_dfaces_primal(mf_to_nf[d+1,0+1],mf_to_nf[0+1,d+1],ldface_to_lvertex)
    mf_to_nf[d,0+1] = compute_dface_to_vertices( mf_to_nf[d+1,0+1], mf_to_nf[d+1,d],ldface_to_lvertex,n_dfaces[d])
    mf_to_nf[0+1,d] = compute_nface_to_dfaces_dual( mf_to_nf[d,0+1],n_dfaces[0+1])
    mf_to_nf[d,d+1] = compute_nface_to_dfaces_dual( mf_to_nf[d+1,d],n_dfaces[d])
  end

  for d in 0:D-1
    mf_to_nf[d+1,d+1] = compute_dface_to_dface(n_dfaces[d+1],Val{T}())
  end

  # Build ldface_to_vertex in lookup tables/ first hard coded
  # Complete matrix
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

const tetn_ldface_to_lvertex =
[Vector{Vector{Int}}[],
 [[[1,2],[1,3],[2,3]]],
 [[[1,2],[1,3],[2,3],[1,4],[2,4],[3,4]],[[1,2,3],[1,2,4],[1,3,4],[1,2,4]]]]

num_dims(s::SurfaceMesh{D}) where D = D

function get_dface_to_nfaces(s::SurfaceMesh,d::Integer,n::Integer)
  s.mfaces_to_nfaces[d+1,n+1]
end

function get_dface_to_vertices(s::SurfaceMesh,d::Integer)
  get_dface_to_nfaces(s,d,0)
end

get_vertex_coordinates(s::SurfaceMesh) = s.vertex_coordinates

@inline function num_vertices(s::SurfaceMesh)
  length(s.vertex_coordinates)
end

@inline function num_faces(s::SurfaceMesh,d::Integer)
  length( get_dface_to_vertices(s,d) )
end

function num_faces(s::SurfaceMesh)
  s.d_to_offset[end]
end

@inline function global_dface(s::SurfaceMesh,d::Integer,lid::Integer)
  s.d_to_offset[d+1]+lid
end

function dface_dimension(s::SurfaceMesh{D},gid::Integer) where D
  for d in 0:D-1
    if s.d_to_offset[d+2] ≥ gid
      return d
    end
  end
  @check false
  return UNSET
end

function local_dface(s::SurfaceMesh,gid::Integer,d::Integer)
  gid - s.d_to_offset[d+1]
end

function get_vertex(s::SurfaceMesh,i::Int)
  get_dface(s,Val{0}(),i)
end

function get_edge(s::SurfaceMesh,i::Int)
  get_dface(s,Val{1}(),i)
end

function get_facet(s::SurfaceMesh{D},i::Int) where D
  get_dface(s,Val{D-1}(),i)
end

function get_dface(s::SurfaceMesh,::Val{d}, i::Int) where d
  error("Not implemented for dimension = $d")
end

function get_dface(s::SurfaceMesh,::Val{0}, i::Int)
  s.vertex_coordinates[i]
end

function get_dface(s::SurfaceMesh,::Val{1}, i::Int)
  df_to_v = get_dface_to_vertices(s,1)
  v = get_vertex_coordinates(s)
  Segment(v[df_to_v[i,1]],v[df_to_v[i,2]])
end

function get_dface(s::SurfaceMesh,::Val{2}, i::Int)
  df_to_v = get_dface_to_vertices(s,2)
  v = get_vertex_coordinates(s)
  Triangle( v[df_to_v[i,1]], v[df_to_v[i,2]], v[df_to_v[i,3]] )
end

function compute_nface_to_dfaces_dual(nface_to_dfaces::Table,n_dfaces::Int)

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
            n_nfaces_around
            for i in 1:n_nfaces_around
              nface_around = vertex_to_nfaces[vertex,i]
              nfaces_around[i] = nface_around
            end
          else
            n_nfaces_around
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

function compute_dface_to_dface(n_dfaces,::Val{T}) where T
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
    num_dfaces = length(dface_to_vertices)
    vtk_type = VTKCellType(d_to_vtk_type_id[d])
    for i in 1:num_dfaces
      vertices = [ dface_to_vertices[i,j] for j in 1:vtk_type.nodes ]
      push!( cells, MeshCell(vtk_type,vertices) )
    end
  end
  vtkfile = vtk_grid(file_base_name,points,cells)
  vtk_save(vtkfile)
end

function have_intersection(bb::BoundingBox,sm::SurfaceMesh,gid::Integer)
  d = dface_dimension(sm,gid)
  lid = local_dface(sm,gid,d)
  have_intersection(bb,sm,d,lid)
end

function BoundingBox(sm::SurfaceMesh,gid::Integer)
  d = dface_dimension(sm,gid)
  lid = local_dface(sm,gid,d)
  BoundingBox(sm,d,lid)
end

BoundingBox(s::SurfaceMesh) = BoundingBox(s.vertex_coordinates)

@generated function BoundingBox(s::SurfaceMesh{D},d::Integer,i::Integer) where D
  body = ""
  for d in 0:D-1
    body *= "if d == $d \n"
    body *= "  f$d = get_dface(s,Val{$d}(),i) \n"
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
    body *= "  f$d = get_dface(s,Val{$d}(),i) \n"
    body *= "  have_intersection(f$d,b) \n"
    body *= "else"
  end
  error = "throw(ArgumentError(\" \$d-face does not exist\"))"
  str = "$body \n  $error \nend"
  Meta.parse(str)
end