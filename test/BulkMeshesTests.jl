module BulkMeshesTests

using Test
using STLCutter
import STLCutter: num_dims

struct StructuredBulkMesh{D,T}
  origin::Point{D,T}
  sizes::VectorValue{D,T}
  partition::NTuple{D,Int}
end

num_dims(::StructuredBulkMesh{D}) where D = D

function num_cells(m::StructuredBulkMesh{D}) where D
  n = 1
  for d in 1:D
    n *= m.partition[d]
  end
  n
end

function get_cell(m::StructuredBulkMesh{D},i::Integer) where D
  x_min = mutable(Point{3,Float64})
  x_max = mutable(Point{3,Float64})
  p_d = 1
  for d in 1:D
    n_d = ( (i-1) ÷ p_d ) % p[d]
    p_d *= p[d]
    x_min[d] = m.origin[d] + m.sizes[d] * n_d / m.partition[d]
    x_max[d] = m.origin[d] + m.sizes[d] * (n_d+1) / m.partition[d]
  end
  HexaCell(Point(x_min),Point(x_max))
end

o = Point(0.0,0.0,0.0)
s = VectorValue(1.0,1.0,1.0)
p = (2,2,2)
m = StructuredBulkMesh(o,s,p)

@test num_dims(m) == 3

@test num_cells(m) == 8

@test get_cell(m,2).bb.pmin == Point(0.5,0.0,0.0)
@test get_cell(m,2).bb.pmax == Point(1.0,0.5,0.5)
import STLCutter: num_dfaces, have_intersection, local_dface,pushlist!, push_to_list!

function num_dfaces(stl::ConformingSTL{D}) where D
  n = 0
  for d in 0:D-1
    n += num_dfaces(stl,d)
  end
  n
end

@inline have_intersection(hex::HexaCell,stl::ConformingSTL,gid::Int) = have_intersection(hex,stl,local_dface(stl,gid)...)

function compute_cell_to_stl_nfaces(m::StructuredBulkMesh{D},stl::ConformingSTL{D}) where D
  cell_to_stl_nfaces = TableOfVectors{Int}([])
  for i in 1:num_cells(m)
    pushlist!(cell_to_stl_nfaces, Vector{Int}([]))
  end
  for k in 1:num_cells(m)
    hex = get_cell(m,k)
    for stl_nface in 1:num_dfaces(stl)
      if have_intersection(hex,stl,stl_nface)
        push_to_list!(cell_to_stl_nfaces, k, stl_nface )
      end
    end
  end
  cell_to_stl_nfaces
end


o = Point(-0.5,-0.5,0.0)
s = VectorValue(1.0,1.0,1.0)
p = (1,1,1)
m = StructuredBulkMesh(o,s,p)
stl = ConformingSTL(joinpath(@__DIR__,"data/cube.stl"))

cell_to_stl_nfaces = compute_cell_to_stl_nfaces(m,stl)
@test length(getlist(cell_to_stl_nfaces,1)) == num_dfaces(stl)



function find_container(m::StructuredBulkMesh{D},p::Point{D}) where D
  p_d = 1
  k = 0
  for d in 1:D
    n_d = Int(floor( (p[d] - m.origin[d]) * m.partition[d] / m.sizes[d] ))
    k += n_d*p_d
    p_d *= m.partition[d]
  end
  k + 1
end

o = Point(0.0,0.0,0.0)
s = VectorValue(1.0,1.0,1.0)
p = (2,2,2)
m = StructuredBulkMesh(o,s,p)

p = Point(0.4,0.0,0.0)
@test find_container(m,p) == 1

p = Point(0.6,0.0,0.0)
@test find_container(m,p) == 2

p = Point(0.4,0.6,0.0)
@test find_container(m,p) == 3

p = Point(0.6,0.6,0.0)
@test find_container(m,p) == 4

p = Point(0.4,0.0,0.6)
@test find_container(m,p) == 5

p = Point(0.6,0.0,0.6)
@test find_container(m,p) == 6

p = Point(0.4,0.6,0.6)
@test find_container(m,p) == 7

p = Point(0.6,0.6,0.6)
@test find_container(m,p) == 8



end # module
