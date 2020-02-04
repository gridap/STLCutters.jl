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

end # module
