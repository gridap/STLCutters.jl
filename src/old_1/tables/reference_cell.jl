
struct RefCell
  cell_type::String
  num_dims::Int
  coordinates::Array{Float64,2}
  dface_to_vertices::Vector{Array{Int,2}}
  d_face_to_dplus1_face::Vector{Vector{Vector{Int}}}
end

const HexCell = "hex"
const TetCell = "tet"

function RefCell(ctype::String,ndims::Int)
  f_to_dplus1f = Vector{Vector{Vector{Int}}}(undef,max(1,ndims-1))
  if ctype == TetCell
    x = zeros(Float64,ndims,ndims+1)
    for d in 1:ndims
      x[d,d+1] = 1
    end
    num_v = size(x,2)
    c2v = reshape([1:num_v;],:,1)
    nf_to_mf = compute_connectivities(c2v)
    df_to_v = nf_to_mf[:,0]
    for d in 1:ndims-1
      f_to_dplus1f[d] = dual_map( nf_to_mf[d+1,d] )
    end
  elseif ctype == HexCell

    x = zeros(Float64,ndims,2^ndims)
    for i in 1:2^ndims
      x[:,i] = digits( (i-1), base=2, pad=ndims )
    end
    x = x.*2 .- 1
    num_v = size(x,2)
    c2v = reshape([1:num_v;],:,1)
    nf_to_mf = compute_hex_connectivities(c2v,ndims)
    df_to_v = nf_to_mf[:,0]
    for d in 1:ndims-1
      f_to_dplus1f[d] = dual_map( nf_to_mf[d+1,d] )
    end
  else
    throw(ArgumentError("RefCell not implemented for $ctype"))
  end

  RefCell(ctype,ndims,x,df_to_v,f_to_dplus1f)
end


Base.ndims(c::RefCell) = c.num_dims

function num_dfaces(c::RefCell,d::Int)
  size( dface_to_vertices(c,d), 2 )
end

function dface_to_vertices(c::RefCell,d::Int)
  c.dface_to_vertices[d+1]
end

function average(x::Array{T,2}) where T
  sum(x,dims=2)./size(x,2)
end

function center(c::RefCell)
  average(c.coordinates)
end

function dface_center(c::RefCell,d::Int,i::Int)
  average( c.coordinates[ :,dface_to_vertices(c,d)[:,i] ] )
end

function dfaces_around_nface(c::RefCell,d::Int,n::Int,i::Int)
  if d < n || d == ndims(c)
    Int[] 
  elseif d == 0 || n == 0
    throw(ArgumentError("0-face (point) not refined"))
  else
    dfaces = [ i ]
    for k in n:d-1
      dfaces = union(c.d_face_to_dplus1_face[k][dfaces]...)
    end
    dfaces
  end
end

function hex_facet_to_vertices(num_dims::Int)
  num_faces = 2*num_dims
  f2v = zeros(Int,2^(num_dims-1),2*num_dims)
  num_vxf = 2^(num_dims-1)

  bin = zeros(Bool,num_dims)
  for i in 1:2*num_dims
    lv = num_dims - (i-1) ÷ 2
    b = (i-1) % 2
    for j in 1:num_vxf
      ids = setdiff(1:num_dims,lv)
      bin[ids] = digits(j-1,base=2,pad=num_dims-1)
      bin[lv] = b
      f2v[j,i] = sum( [  bin[n]*2^(n-1) for n in 1:length(bin) ] ) + 1
    end
  end
  f2v
end

function compute_face_to_nfaces(f2v::Array{Int,2},nf_to_lv::Array{Int,2})

  function faces_around(i::Int,j::Int)
    
    function nface()
      f2v[nf_to_lv[:,j],i]
    end
    fs = intersect(v2f[nface()]...)
    setdiff(fs,i)
  end
  
  v2f = dual_map(f2v)
  f2nf = zeros(Int,size(nf_to_lv,2),size(f2v,2))
  count_nf = 0
  for i in 1:size(f2v,2), j in 1:size(nf_to_lv,2)
    if f2nf[j,i] == 0
      count_nf += 1
      f2nf[j,i] = count_nf
      for ineig in faces_around(i,j), jneig in 1:size(nf_to_lv,2)
        if i in faces_around(ineig,jneig)
          f2nf[jneig,ineig] = count_nf
        end
      end
    end
  end
  f2nf
end

function compute_nface_to_vertices(f2nf::Array{Int,2},f2v::Array{Int,2},nf_to_lv::Array{Int,2})

  function nface(i::Int,j::Int)
    f2v[nf_to_lv[:,j],i]
  end

  num_nf = maximum(f2nf)
  num_vxnf = size(nf_to_lv,1)
  nf2v = zeros(Int,num_vxnf,num_nf)
  count_nf = 1
  for i in 1:size(f2nf,2), j in 1:size(f2nf,1)
    if f2nf[j,i] == count_nf
      nf2v[:,count_nf] = nface(i,j)
      count_nf += 1
    end
  end
  nf2v
end

function compute_hex_connectivities(c2v::Array{Int,2},num_dims::Int)
  
  nF_to_mF = NFaceToMFace(num_dims)
  nF_to_mF[num_dims,0] = c2v
  num_dfaces = zeros(Int,num_dims+1)
  num_dfaces[0+1] = maximum(c2v)
  num_dfaces[num_dims+1] = size(c2v,2)

  for d in reverse(1:num_dims-1)
    df_to_lv = hex_facet_to_vertices(d+1)
    nF_to_mF[d+1,d] = compute_face_to_nfaces( nF_to_mF[d+1,0], df_to_lv )
    nF_to_mF[d,0] = compute_nface_to_vertices( nF_to_mF[d+1,d], nF_to_mF[d+1,0], df_to_lv )
    num_dfaces[d+1] = size( nF_to_mF[d,0], 2 )
  end
  for d in reverse(1:num_dims-2)
    nF_to_mF[num_dims,d] = chain_maps( nF_to_mF[num_dims,d+1], nF_to_mF[d+1,d] )
  end
  for d in 0:num_dims
    nF_to_mF[d,d] = reshape([1:num_dfaces[d+1];],1,:)
  end
  nF_to_mF
end

function compute_face_to_initial_face(cell::RefCell,nf_to_mf::NFaceToMFace,d::Int,i::Int)
  nf_to_v = nf_to_mf[:,0]
  v_to_nF = [ dual_map( dface_to_vertices(cell,d) ) for d in 0:ndims(cell) ]
  for n in 0:ndims(cell) 
    push!(v_to_nF[n+1], dfaces_around_nface(cell,n,d,i) )
  end
  nf_to_nF = compute_face_to_initial_face(v_to_nF,nf_to_v)
end
