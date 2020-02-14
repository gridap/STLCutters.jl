using PyCall
using LinearAlgebra

scipy_spatial = pyimport("scipy.spatial")


function delaunay(x::Array{T,2}) where T
  x = copy(transpose(x))
  d = scipy_spatial.Delaunay(x)
  p2j(a::Array{T,2}) where T = copy(transpose(a)) .+ 1
  p2j(d.simplices), p2j(d.neighbors)
end

x=reshape([0 0  0 1  1 0  1 1],2,:)
c=delaunay(x)

x=reshape([0 0 0  0 1 0  1 0 0  0 0 1  1 1 1],3,:)
c=delaunay(x)


struct RefCell
  cell_type::String
  num_dims::Int
  coordinates::Array{Float64,2}
  facet_to_vertices::Array{Int,2}
end


const HexCell = "hex"
const TetCell = "tet"
function RefCell(ctype::String,ndims::Int)
  if ndims == 2
    if ctype == TetCell
      x = reshape([ 0. 0.  1. 0.  0. 1. ],ndims,:)
      f2v = reshape([ 1 2  1 3  2 3 ],ndims,:)
    elseif ctype == HexCell
      x = reshape([ -1. -1.  1. -1.  -1. 1.  1. 1. ],ndims,:)
      f2v = reshape([ 1 2  3 4  1 3  2 4 ],2^(ndims-1),:)
    else
      throw(ArgumentError("Cell type $ctype not defined"))
    end 
  elseif ndims == 3
    if ctype == TetCell
      x = reshape([ 0. 0. 0.  1. 0. 0.  0. 1. 0.  0. 0. 1. ],ndims,:)
      f2v = reshape([ 1 2 3  1 2 4  1 3 4  2 3 4 ],ndims,:)
    elseif ctype == HexCell
      x = reshape(
        [ -1. -1. -1.  1. -1. -1.  -1. 1. -1.  1. 1. -1.  -1. -1. 1.  1. -1. 1.  -1. 1. 1.  1. 1. 1.],
        ndims,:)
      f2v = reshape( [ 1 2 3 4  5 6 7 8  1 2 5 6  3 4 7 8  1 3 5 7  2 4 6 8],2^(ndims-1),:)
    else
      throw(ArgumentError("Cell type $ctype not defined"))
    end
  else
    throw(ArgumentError("RefCell not defined for $ndims dimensions"))
  end
  RefCell(ctype,ndims,x,f2v)
end

@show RefCell(TetCell,2)
@show RefCell(TetCell,3)
@show RefCell(HexCell,2)
@show RefCell(HexCell,3)

rh = RefCell(HexCell,3)

x = rh.coordinates
c2v,c2n = delaunay(rh.coordinates)


function average(x::Array{T,2}) where T
  sum(x,dims=2)./size(x,2)
end

function center(c::RefCell)
  average(c.coordinates)
end

function facet_center(c::RefCell,i::Int)
  average( c.coordinates[ :,c.facet_to_vertices[:,i] ] )
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
  d_map = [ [] for i in 1:n ]
  for i in 1:size(map,2), j in 1:size(map,1)
    push!( d_map[ map[j,i] ],i )
  end
  d_map
end

"""
  `compute_face_to_nfaces(f2v::Array{Int,2})`

 - Given n-faces to vertices (or 0-faces) map:

    * computes the map from n-faces to (n-1)-faces 

    * n-faces are called `faces` or just `f`

    * (n-1)-faces are called `nfaces` or just `nf`

    * All n-faces are simplicies
"""  
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
    sign( dot(n1,v2) )
  elseif size(x,1) == 3
    size(x,2) == 4 || throw(ErrorException("Not a 3D simplex"))
    v1 = x[:,2] - x[:,1]
    v2 = x[:,3] - x[:,1]
    v3 = x[:,4] - x[:,1]
    sign( (v1×v2) ⋅ v3 )
  else
    throw(DimensionMismatch("Orientation of a $(size(x,1))D simplex not implemented"))
  end
end

function correct_cell_orientation(x::Array{Float64,2},c2v::Array{Int,2})
  num_cells = size(c2v,2)
  for i in 1:num_cells
    if orientation(x[:,c2v[:,i]]) < 0 
      println("Cell $i has negative orientation")
      c2v[[1,end],i] = reverse(c2v[[1,end],i])
    end
  end
  c2v
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

struct NFaceToMFace
  data::Vector{Vector{Array{Int,2}}}
end

Base.getindex(a::NFaceToMFace,i::Int,j::Int) = a.data[i+1][j+1]
Base.setindex!(a::NFaceToMFace,v::Array{Int,2},i::Int,j::Int) = a.data[i+1][j+1] = v

function NFaceToMFace(ndims::Int)
  data = [[ zeros(Int,0,0) for i in 0:num_dims] for j in 0:num_dims ]
  NFaceToMFace(data)
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

t2 = RefCell(TetCell,2)

x = [ t2.coordinates center(t2) ]

c2v,c2n = delaunay( x )

c2f = compute_cell_to_facets(c2n)
f2v = compute_nface_to_vertices(c2f,c2v)
@show f2v


t3 = RefCell(TetCell,3)

x = [ t3.coordinates center(t3) ]

c2v,c2n = delaunay( x )

c2n = reverse(c2n,dims=1)
c2f = compute_cell_to_facets(c2n)
f2v = compute_nface_to_vertices(c2f,c2v)
v2c = dual_map(c2v)
v2f = dual_map(f2v)


@show f2e = compute_face_to_nfaces(f2v)
@show e2v = compute_nface_to_vertices(f2e,f2v)


@show c2f == compute_face_to_nfaces(c2v)






#c2v = correct_cell_orientation(x,c2v)
@show c2e=chain_maps(c2f,f2e)





t3 = RefCell(TetCell,3)

x = [ t3.coordinates center(t3) ]

f2f = compute_mesh(x)


# TODO: 
# Considering renaming when refering to nface and (n-1)face by face(f) and nface(nf) respectvely
# as cell(c) and facet(f)
# Facet to subcell orientation
# subNfacet to Nfacet
# [x] flip nface creation (opposite to N-i-1)
# [x] reorient cells
# [x] create container of submesh
# build refCell connectivities with the same functions
# test
# rename comp_c2nf as chain map or something similar








