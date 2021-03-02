
struct NFaceToMFace
  data::Array{Array{Int,2},2}
end

function NFaceToMFace(num_dims::Int)
  data = reshape([ zeros(Int,0,0) for i in 0:num_dims, j in 0:num_dims ],num_dims+1,num_dims+1) 
  NFaceToMFace(data)
end

Base.getindex(a::NFaceToMFace,i::Int,j::Int) = a.data[i+1,j+1]
Base.getindex(a::NFaceToMFace,I::Colon,j::Int) = a.data[I,j+1]
Base.getindex(a::NFaceToMFace,i::Int,J::Colon) = a.data[i+1,J]
Base.getindex(a::NFaceToMFace,I::UnitRange{Int},j::Int) = a.data[I.+1,j+1]
Base.getindex(a::NFaceToMFace,i::Int,J::UnitRange{Int}) = a.data[i+1,J.+1]
Base.getindex(a::NFaceToMFace,I::UnitRange{Int},J::UnitRange{Int}) = a.data[I.+1,J.+1]

Base.setindex!(a::NFaceToMFace,v::Array{Int,2},i::Int,j::Int) = a.data[i+1,j+1] = v

Base.lastindex(a::NFaceToMFace,i::Int) = lastindex(a.data,i) - 1

Base.ndims(f::NFaceToMFace) = size(f.data,1) - 1

get_data(a::NFaceToMFace) = a.data

function get_m_of_v_of_v(a::NFaceToMFace)
  vectors.(a.data)
end

function vectors(a::Matrix)
  [ a[:,i] for i in 1:size(a,2) ]
end

function vectors(a::NFaceToMFace) 
  [a.data[i,:] for i in 1:size(a.data,1)]
end

function Base.print(io::IO,a::Vector{NFaceToMFace})
  v = get_m_of_v_of_v.(a)
  print(io,v)
end

function Base.print(io::IO,a::NFaceToMFace) 
  m = get_m_of_v_of_v(a)
  T = typeof(m)
  data = string( m )
  str = "$(string(T))(\n  $data )"
  print(io,str)
end

function Base.string(a::Matrix{<:Array})
  data = ""
  for i in a 
    data *= "   $(string(i)),\n"
  end
  "reshape([\n$data  ], $(size(a)) )"
end

function Base.string(v::Vector{<:Vector})
  data = join([ "$(string(i))," for i in v ])
  "[$data]"
end

function Base.string(::Type{Vector{T}}) where T
  "Vector{$(string(T))}"
end

function Base.string(::Type{Matrix{T}}) where T
  "Matrix{$(string(T))}"
end


function Base.print(io::IO,a::Vector{T}) where T<:Array#{<:Array}
  data =join([ "$(string(a[i])) , \n \n  " for i in 1:length(a) ])
  data = replace(data, "Array{Int64}(undef,0,0)" => "_empty_matrix")
  str = "$(string(T))[\n  $(rstrip(data)) \n]"
  print(io::IO,str)
end

function delaunay(x::Array{T,2}) where T
  x = copy(transpose(x))
  d = scipy_spatial.Delaunay(x)
  p2j(a::Array{T,2}) where T = copy(transpose(a)) .+ 1
  p2j(d.simplices), p2j(d.neighbors)
end

function compute_cell_to_facets(c2n::Array{Int,2})
  c2f = zeros(Int,size(c2n))
  count_f = 0
  for i in 1:size(c2n,2), j in 1:size(c2n,1)
    if c2f[j,i] == 0
      count_f += 1
      c2f[j,i] = count_f
      if c2n[j,i] != 0
        ineig = c2n[j,i]
        jneig = findfirst( x->x == i, c2n[:,ineig])
        c2f[jneig,ineig] = count_f
      end
    end
  end
  c2f
end

function dual_map(map::Array{Int,2})
  n = maximum(map)
  d_map = [ Int[] for i in 1:n ]
  for i in 1:size(map,2), j in 1:size(map,1)
    push!( d_map[ map[j,i] ],i )
  end
  d_map
end

function chain_maps(a2b::Array{Int,2},b2c::Array{Int,2})
  
  cs(i::Int) = unique(b2c[:,a2b[:,i]])

  num_cxa = length( cs(1) )
  num_a = size(a2b,2)
  a2c = zeros(Int,num_cxa,num_a)
  for i in 1:num_a
    a2c[:,i] = cs(i)
  end
  a2c
end

function compute_face_to_nfaces(f2v::Array{Int,2})

  function faces_around(i::Int,j::Int)
    
    function nface()
      k = size(f2v,1) - j + 1
      setdiff(f2v[:,i],f2v[k,i])
    end

    fs = intersect(v2f[nface()]...)
    setdiff(fs,i)
  end
  
  v2f = dual_map(f2v)
  f2nf = zeros(Int,size(f2v))
  count_nf = 0
  for i in 1:size(f2v,2), j in 1:size(f2v,1)
    if f2nf[j,i] == 0
      count_nf += 1
      f2nf[j,i] = count_nf
      for ineig in faces_around(i,j), jneig in 1:size(f2v,1)
        if i in faces_around(ineig,jneig)
          f2nf[jneig,ineig] = count_nf
        end
      end
    end
  end
  f2nf
end

function compute_nface_to_vertices(f2nf::Array{Int,2},f2v::Array{Int,2})

  function nface(i::Int,j::Int)
    k = size(f2v,1) - j + 1
    setdiff(f2v[:,i],f2v[k,i])
  end

  num_nf = maximum(f2nf)
  num_vxnf = size(f2nf,1) - 1
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


function signed_volume(x::Array{Float64,2})
  size(x,1) + 1 == size(x,2) || throw(DimensionMismatch("volume only implemented for symplices"))
  A = zeros(size(x,1),size(x,1))
  for d in 1:size(x,1)
    A[:,d] = x[:,d+1] - x[:,1]
  end
  det(A) / factorial(size(x,1))
end


function fix_cell_volume(x::Array{Float64,2},c2v::Array{Int,2})
  num_cells = size(c2v,2)
  count = 0
  for i in 1:num_cells
    if signed_volume(x[:,c2v[:,i]]) < 0 
      c2v[[end-1,end],i] = reverse(c2v[[end-1,end],i])
    end

    @assert (signed_volume(x[:,c2v[:,i]]) ≥ 0) "Volume still negative after permutation of `x = $((x[:,c2v[:,i]]))`"

    if signed_volume(x[:,c2v[:,i]]) > 1e-5
      count += 1
      c2v[:,count] = c2v[:,i]
    end
  end
  c2v[:,1:count]
end

function compute_connectivities(c2v::Array{Int,2})
  
  num_dims = size(c2v,1) - 1
  nF_to_mF = NFaceToMFace(num_dims)
  nF_to_mF[num_dims,0] = c2v
  num_dfaces = zeros(Int,num_dims+1)
  num_dfaces[0+1] = maximum(c2v)
  num_dfaces[num_dims+1] = size(c2v,2)

  for d in reverse(1:num_dims-1)
    nF_to_mF[d+1,d] = compute_face_to_nfaces( nF_to_mF[d+1,0] )
    nF_to_mF[d,0] = compute_nface_to_vertices( nF_to_mF[d+1,d], nF_to_mF[d+1,0] )
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

function compute_mesh(x::Array{Float64,2})
  if size(x,1) == 1
    @assert size(x,2) == 3
    nf_to_mf = NFaceToMFace(1)
    nf_to_mf[0,0] = [ 1 2 3 ]
    nf_to_mf[1,0] = [ 1 3; 3 2 ]
    nf_to_mf[1,1] = [ 1 2 ]
  else
    c2v, = delaunay( x )
    c2v = fix_cell_volume(x,c2v)
    nf_to_mf = compute_connectivities(c2v)
    D = size(x,1)
    #fix_boundary_orientation!( nf_to_mf[D,D-1], nf_to_mf[D,0], nf_to_mf[D-1,0], x)
  end
  nf_to_mf
end

function compute_cell_to_facet_orientation(x::Array{Float64,2},c2v::Array{Int,2},f2v::Array{Int,2},c2f::Array{Int,2})
  num_c = size(c2f,2)
  num_fxc = size(c2f,1)
  c2f_orientation = zeros(Int,num_fxc,num_c)

  for i in 1:num_c, j in 1:num_fxc
    c = c2v[:,i]
    f = f2v[:,c2f[j,i]]
    cp = vcat( f, setdiff(c,f) )
    c2f_orientation[j,i] = -sign( signed_volume(x[:,cp]) )
  end
  c2f_orientation
end

function compute_cell_to_facet_orientation(x::Array{Float64,2},nf_to_mf::NFaceToMFace)
  D = ndims(nf_to_mf)
  c2v=nf_to_mf[D,0]
  f2v=nf_to_mf[D-1,0]
  c2f=nf_to_mf[D,D-1]
  compute_cell_to_facet_orientation(x,c2v,f2v,c2f) 
end

function compute_face_to_initial_face(v_to_nF::Vector{Vector{Vector{Int}}},nf_to_v::Vector{Array{Int,2}})
  nf_to_nF = [ zeros(Int,size( nf_to_v[d],2)) for d in 1:length(v_to_nF) ]
  for d in 1:length(v_to_nF)
    for i in 1:size( nf_to_v[d],2)
      F = intersect(v_to_nF[d][ nf_to_v[d][:,i]]...)
      if length(F) == 1
        nf_to_nF[d][i] = F[1]
      else
        length(F) == 0 || throw(ErrorException("$d-Face of $i not idenfitied: intersection = $F"))
      end
    end
  end
  nf_to_nF
end

function fix_boundary_orientation!(c_to_f::Matrix,c_to_v::Matrix,f_to_v::Matrix,x::Matrix)
  for i in 1:size(c_to_f,2), j in 1:size(c_to_f,1)
    c = c_to_v[:,i]
    f = f_to_v[:,c_to_f[j,i]]
    cp = vcat( f, setdiff(c,f) )
    if sign( signed_volume(x[:,cp]) ) > 0
      f_to_v[ [end-1,end], c_to_f[j,i] ] = reverse( f_to_v[ [end-1,end], c_to_f[j,i] ] ) 
    end
  end
  f_to_v
end

