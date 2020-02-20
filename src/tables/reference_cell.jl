
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
  f_to_dplus1f = Vector{Vector{Vector{Int}}}(undef,max(1,ndims-2))
  if ctype == TetCell
    x = zeros(Float64,ndims,ndims+1)
    for d in 1:ndims
      x[d,d+1] = 1
    end
    c2v = reshape(collect(1:ndims+1),:,1)
    nf_to_mf = compute_connectivities(c2v)
    df_to_v = nf_to_mf[1:end,0]
    for d in 1:ndims-2
      f_to_dplus1f[d] = dual_map( nf_to_mf[d+1,d] )
    end
  elseif ctype == HexCell

    x = zeros(Float64,ndims,2^ndims)
    for i in 1:2^ndims
      x[:,i] = digits( (i-1), base=2, pad=ndims )
    end
    x = x.*2 .- 1 
    c2v = reshape(collect(1:2^ndims),:,1)
    f2v = hex_facet_to_vertices(ndims)
    if ndims == 2
      df_to_v = [ f2v, c2v ]
    elseif ndims == 3
      e2v = reshape( [ 1 2  3 4  5 6  7 8  1 3  2 4  5 7  6 8  1 5  2 6  3 7  4 8 ],2,: )
      df_to_v = [ e2v, f2v, c2v ]
    else
      throw(ArgumentError(""))
    end
  else
    throw(ArgumentError("RefCell not implemented for $ctype"))
  end

  RefCell(ctype,ndims,x,df_to_v,f_to_dplus1f)
end


Base.ndims(c::RefCell) = c.num_dims

function num_dfaces(c::RefCell,d::Int)
  if d == 0
    size(c.coordinates,2)
  else
    size(c.dface_to_vertices[d],2)
  end
end

function average(x::Array{T,2}) where T
  sum(x,dims=2)./size(x,2)
end

function center(c::RefCell)
  average(c.coordinates)
end

function facet_center(c::RefCell,i::Int)
  average( c.coordinates[ :,c.facet_to_vertices[:,i] ] )
end

function dface_center(c::RefCell,d::Int,i::Int)
  average( c.coordinates[ :,c.dface_to_vertices[d][:,i] ] )
end

function dfaces_around_nface(c::RefCell,d::Int,n::Int,i::Int)
  if n == ndims(c) || d < n
    Int[]
  elseif d == ndims(c) 
    throw(ArgumentError("$(d)-face (cell) not in boundary"))
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

function compute_face_to_initial_face(cell::RefCell,nf_to_mf::NFaceToMFace,d::Int,i::Int)
  nf_to_v = nf_to_mf[1:end-1,0]
  v_to_nF = [ dual_map( cell.dface_to_vertices[d] ) for d in 1:ndims(cell)-1 ]
  for n in 1:ndims(cell)-1 
    push!(v_to_nF[n], dfaces_around_nface(cell,n,d,i) )
  end
  nf_to_nF = compute_face_to_initial_face(v_to_nF,nf_to_v)
end
