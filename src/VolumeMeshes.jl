struct VolumeMesh{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  mfaces_to_nfaces::Matrix{Table{Int}}
end

function VolumeMesh(m::CartesianMesh)
  VolumeMesh( 
    get_vertex_coordinates(m),
    _compute_mfaces_to_nfaces(m) )
end

num_dims(vm::VolumeMesh{D}) where D = D

function get_faces(vm::VolumeMesh,d::Integer,n::Integer)
  vm.mfaces_to_nfaces[d+1,n+1]
end

function get_dface_to_vertices(vm::VolumeMesh,d::Integer)
  get_faces(vm,d,0)
end

function num_faces(vm::VolumeMesh,d::Integer)
  length( get_dface_to_vertices(vm,d) )
end


function get_vertex_coordinates(vm::VolumeMesh)
  vm.vertex_coordinates
end

function get_vertex_coordinates(vm::VolumeMesh,i::Integer)
  get_face_coordinates(vm,Val{0}(),i)
end

function get_edge_coordinates(vm::VolumeMesh,i::Integer)
  get_face_coordinates(vm,Val{1}(),i)
end

function get_facet_coordinates(vm::VolumeMesh,i::Integer) where D
  get_face_coordinates(vm,Val{D-1}(),i)
end

function get_cell_coordinates(vm::VolumeMesh{D},i::Integer) where D
  get_face_coordinates(vm,Val{D}(),i)
end

function get_face_coordinates(vm::VolumeMesh,::Val{0},i::Integer)
  v = get_vertex_coordinates(vm)
  v[i]
end

function get_face_coordinates(vm::VolumeMesh,::Val{1},i::Integer)
  v = get_vertex_coordinates(vm)
  f_to_v = get_dface_to_vertices(vm,1)
  Segment( v[f_to_v[i,1]], v[f_to_v[i,2]] )
end

function get_face_coordinates(vm::VolumeMesh,::Val{2},i::Integer)
  v = get_vertex_coordinates(vm)
  f_to_v = get_dface_to_vertices(vm,2)
  Quadrilater( v[f_to_v[i,1]], v[f_to_v[i,2]], v[f_to_v[i,3]], v[f_to_v[i,4]] )
end

function get_face_coordinates(vm::VolumeMesh,::Val{3},i::Integer)
  v = get_vertex_coordinates(vm)
  f_to_v = get_dface_to_vertices(vm,3)
  Hexahedron( 
    v[f_to_v[i,1]], v[f_to_v[i,2]], v[f_to_v[i,3]], v[f_to_v[i,4]],
    v[f_to_v[i,5]], v[f_to_v[i,6]], v[f_to_v[i,7]], v[f_to_v[i,8]] )
end

function _compute_mfaces_to_nfaces(m::CartesianMesh{D}) where D
  mf_to_nf = Matrix{Table{Int}}(undef,D+1,D+1)
  
  mf_to_nf[D+1,0+1] = get_cell_to_vertices(m)
  mf_to_nf[0+1,D+1] = compute_nface_to_dfaces_dual( mf_to_nf[D+1,0+1], num_vertices(m) )

  for d in 1:D-1
    mf_to_nf[d+1,0+1] = face_to_vertices(m,d)
    mf_to_nf[0+1,d+1] = compute_nface_to_dfaces_dual( mf_to_nf[d+1,0+1], num_vertices(m) )
  end

  for d in 0:D
    mf_to_nf[d+1,d+1] = compute_dface_to_dface( _num_faces(m,d), Int )
  end

  for d in 2:D, n in 1:d-1
    ldface_to_lvertex = _get_hexahedral_dface_to_lvertices(d,n)
    mf_to_nf[d+1,n+1] = compute_nface_to_dfaces(mf_to_nf[d+1,0+1],mf_to_nf[0+1,n+1],ldface_to_lvertex)
    mf_to_nf[n+1,d+1] = compute_nface_to_dfaces_dual( mf_to_nf[d+1,n+1], length(mf_to_nf[n+1,0+1]) )
  end

  mf_to_nf
end

function _get_hexahedral_dface_to_lvertices(D::Integer,d::Integer)
  D_to_d_to_dface_to_vertices_for_hex[D][d]
end

function compute_nface_to_dfaces(nface_to_vertices::Table,vertex_to_dfaces::Table,ldface_to_lvertices)
  T = eltype(nface_to_vertices)
  n_vertices = length(vertex_to_dfaces)
  max_dfaces_around = 0
  for vertex in 1:n_vertices
    dfaces_around = length(vertex_to_dfaces,vertex)
    max_dfaces_around = max(max_dfaces_around,dfaces_around)
  end
  dfaces_around = fill(T(UNSET),max_dfaces_around)
  dfaces_around_bis = fill(T(UNSET),max_dfaces_around)
  n_dfaces_around = UNSET
  n_dfaces_around_bis = UNSET
  n_nfaces = length( nface_to_vertices )
  n_ldfaces = length( ldface_to_lvertices )
  nf_to_df = Table{Int}( undef, n_nfaces, n_ldfaces )
  for nface in 1:n_nfaces
    for ldface in 1:n_ldfaces
      lvertices = ldface_to_lvertices[ldface]
      for (i,lvertex) in enumerate(lvertices)
        vertex = nface_to_vertices[nface,lvertex]
        if i == 1
          n_dfaces_around = length(vertex_to_dfaces,vertex)
          for j in 1:n_dfaces_around
            dfaces_around[j] = vertex_to_dfaces[vertex,j]
          end
        else
          n_dfaces_around_bis = length(vertex_to_dfaces,vertex)
          for j in 1:n_dfaces_around_bis
            dfaces_around_bis[j] = vertex_to_dfaces[vertex,j]
          end
          _set_intersection!(
              dfaces_around,
              dfaces_around_bis,
              n_dfaces_around,
              n_dfaces_around_bis)
        end
      end
      dface = UNSET
      for i in 1:n_dfaces_around
        if dfaces_around[i] != UNSET
          dface = dfaces_around[i]
          break
        end
      end
      @check dface != UNSET
      nf_to_df[nface,ldface] = dface    
    end
  end
  nf_to_df
end

@generated function distance(m::VolumeMesh{D},d::Int,i::Int,p::Point{D}) where D
  body = ""
  for d in 0:D
    body *= "if d == $d \n"
    body *= "  f$d = get_face_coordinates(m,Val{$d}(),i) \n"
    body *= "  distance(p,f$d) \n"
    body *= "else"
  end
  error = "  throw(ArgumentError(\"\$d-face does not exist\"))\n"
  str = body * '\n' * error * "end"
  Meta.parse(str)
end

