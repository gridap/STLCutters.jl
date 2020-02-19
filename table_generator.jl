using PyCall
using LinearAlgebra

scipy_spatial = pyimport("scipy.spatial")

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

function orientation(x::Array{Float64,2})
  if size(x,1) == 2
    size(x,2) == 3 || throw(ErrorException("Not a 2D simplex"))
    v1 = x[:,2] - x[:,1]
    v2 = x[:,3] - x[:,1]
    n1 = [ -v1[2] ; v1[1] ]
    Int(sign( dot(n1,v2) ))
  elseif size(x,1) == 3
    size(x,2) == 4 || throw(ErrorException("Not a 3D simplex"))
    v1 = x[:,2] - x[:,1]
    v2 = x[:,3] - x[:,1]
    v3 = x[:,4] - x[:,1]
    Int(sign( (v1×v2) ⋅ v3 ))
  else
    throw(DimensionMismatch("Orientation of a $(size(x,1))D simplex not implemented"))
  end
end

function correct_cell_orientation(x::Array{Float64,2},c2v::Array{Int,2})
  num_cells = size(c2v,2)
  for i in 1:num_cells
    if orientation(x[:,c2v[:,i]]) < 0 
      c2v[[1,end],i] = reverse(c2v[[1,end],i])
    end
  end
  c2v
end

function compute_connectivities(c2v::Array{Int,2})
  
  num_dims = size(c2v,1) - 1
  nF_to_mF = NFaceToMFace(num_dims)
  nF_to_mF[num_dims,0] = c2v

  for d in reverse(1:num_dims-1)
    nF_to_mF[d+1,d] = compute_face_to_nfaces( nF_to_mF[d+1,0] )
    nF_to_mF[d,0] = compute_nface_to_vertices( nF_to_mF[d+1,d], nF_to_mF[d+1,0] )
  end

  for d in reverse(1:num_dims-2)
    nF_to_mF[num_dims,d] = chain_maps( nF_to_mF[num_dims,d+1], nF_to_mF[d+1,d] )
  end
  nF_to_mF
end

function compute_mesh(x::Array{Float64,2})
  c2v, = delaunay( x )
  c2v = correct_cell_orientation(x,c2v)
  compute_connectivities(c2v)
end

function compute_cell_to_facet_orientation(x::Array{Float64,2},c2v::Array{Int,2},f2v::Array{Int,2},c2f::Array{Int,2})
  num_c = size(c2f,2)
  num_fxc = size(c2f,1)
  c2f_orientation = zeros(Int,num_fxc,num_c)

  for i in 1:num_c, j in 1:num_fxc
    c = c2v[:,i]
    f = f2v[:,c2f[j,i]]
    cp = vcat( f, setdiff(c,f) )
    c2f_orientation[j,i] = -orientation(x[:,cp])
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

# RefCell

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
      dfaces = union(c.d_face_to_dplus1_face[n][dfaces]...)
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


function Base.string(a::Vector{Array{T,N}}) where {T,N}
  data =join([ "$(string(a[i])), \n  " for i in 1:length(a) ])
  length(data) >= 3 ? data = data[1:end-3] : data = ""
  "[$data]"
end

function vectors(a::NFaceToMFace) 
  [a.data[i,:] for i in 1:size(a.data,1)]
end

function Base.print(io::IO,a::Vector{T}) where T<:Array
  data =join([ "$(string(a[i])) , \n \n  " for i in 1:length(a) ])
 # data = replace(data, ";" => "; \n  ")
  data = replace(data, "Array{Int64}(undef,0,0)" => "empty_matrix")
  str = "$(string(T))[\n  $(rstrip(data)) \n]"
  print(io::IO,str)
end

function Base.print(io::IO,a::Vector{NFaceToMFace})
  v = vectors.(a)
  print(io,v)
end




cell = RefCell(TetCell,3)
prefix = "cut_$(cell.cell_type)$(ndims(cell))"
filename = "$(prefix)_tables.jl"
f = open(filename,"w")

header = 
"#
# Cut TET tables 
#  Please do not modify it by hand:
#    code automatically generated by $(basename(@__FILE__))
#
"

println(f, header )


println(f, "const empty_matrix = Array{Int64}(undef,0,0)" ); println(f)


for num_dims in 2:3
  
  cell = RefCell(TetCell,num_dims)
  prefix = "cut_$(cell.cell_type)$(ndims(cell))"
  println(f,"## Combinations for $prefix"); println(f)

  global_x = Array{Float64,2}[]
  global_nf_to_mf = NFaceToMFace[]
  global_nf_to_nF = Vector{Vector{Int}}[]
  global_c2f_orientation = Array{Int,2}[]
  for d in 1:ndims(cell)
    for i in 1:num_dfaces(cell,d)
    x = [ cell.coordinates dface_center(cell,d,i) ]
    nf_to_mf = compute_mesh( x )
    nf_to_nF = compute_face_to_initial_face(cell,nf_to_mf,d,i)
    c2f_orientation = compute_cell_to_facet_orientation(x,nf_to_mf)

    push!( global_x, x )
    push!( global_nf_to_mf, nf_to_mf )
    push!( global_nf_to_nF, nf_to_nF )
    push!( global_c2f_orientation, c2f_orientation )
    end
  end
  print(f, "const $(prefix)_nface_to_mface = " ); println(f, global_nf_to_mf ); println(f)
  print(f, "const $(prefix)_nsubface_to_nface = " ); println(f, global_nf_to_nF ); println(f)
  print(f, "const $(prefix)_orientation_cell_to_facet = " ); println(f, global_c2f_orientation ); println(f)
end
close(f)

Base.print(io::IO,a::NFaceToMFace) = print(io,vectors(a))



filename = "hex_to_tet_tables.jl"
f = open(filename,"w")

header = 
"#
# Split HEX into TETs tables 
#  Please do not modify it by hand:
#    code automatically generated by $(basename(@__FILE__))
#
"
println(f, header )

for num_dims in 2:3
  prefix = "hex_to_tet$(num_dims)"
  cell = RefCell(HexCell,num_dims)
  println(f,"## Connectivities for $prefix"); println(f)

  x = cell.coordinates
  nf_to_mf = compute_mesh( x )
  nf_to_v = nf_to_mf[1:end-1,0]
  v_to_nF = [ dual_map( cell.dface_to_vertices[d] ) for d in 1:ndims(cell)-1 ]
  nf_to_nF = compute_face_to_initial_face(v_to_nF,nf_to_v)
  c2f_orientation = compute_cell_to_facet_orientation(x,nf_to_mf)

  print(f, "const $(prefix)_nface_to_mface = " ); println(f, nf_to_mf ); println(f);
  print(f, "const $(prefix)_nsubface_to_mface = " ); println(f, nf_to_nF ); println(f);
  print(f, "const $(prefix)_orientation_cell_to_facet = " ); println(f, c2f_orientation ); println(f);
end
close(f)




