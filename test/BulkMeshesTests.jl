module BulkMeshesTests

using Test
using STLCutter

using STLCutter: BoundingBox,find_container, get_cell_id,compute_cell_to_stl_nfaces,optimized_compute_cell_to_stl_nfaces

o = Point(0.0,0.0,0.0)
s = VectorValue(1.0,1.0,1.0)
p = (2,2,2)
m = StructuredBulkMesh(o,s,p)

@test num_dims(m) == 3

@test num_cells(m) == 8

@test get_cell(m,2).pmin == Point(0.5,0.0,0.0)
@test get_cell(m,2).pmax == Point(1.0,0.5,0.5)

o = Point(-0.5,-0.5,0.0)
s = VectorValue(1.0,1.0,1.0)
p = (1,1,1)
m = StructuredBulkMesh(o,s,p)
stl = ConformingSTL(joinpath(@__DIR__,"data/cube.stl"))

cell_to_stl_nfaces = compute_cell_to_stl_nfaces(m,stl)
@test length(getlist(cell_to_stl_nfaces,1)) == num_dfaces(stl)

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

o = Point(0.0,0.0,0.0)
s = VectorValue(1.0,1.0,1.0)
p = (10,1,1)
m = StructuredBulkMesh(o,s,p)

p = Point(0.5,0.5,0.5)
@test find_container(m,p) == 6

o = Point(0.0,0.0,0.0)
s = VectorValue(1.0,1.0,1.0)
p = (2,3,4)
m = StructuredBulkMesh(o,s,p)

opt_c2s = optimized_compute_cell_to_stl_nfaces(m,stl)
non_opt_c2s = compute_cell_to_stl_nfaces(m,stl)

@test opt_c2s == non_opt_c2s

o = Point(-1.0,0.0,0.0)
s = VectorValue(2.0,1.0,1.0)
p = (10,1,1)
m = StructuredBulkMesh(o,s,p)

opt_c2s = optimized_compute_cell_to_stl_nfaces(m,stl)
non_opt_c2s = compute_cell_to_stl_nfaces(m,stl)

@test opt_c2s == non_opt_c2s


struct StructuredBulkMeshCells{D,T}
  mesh::StructuredBulkMesh{D,T}
  num_cells::Int
end

struct StructuredBulkCell{D,T}
  gid::Int
  bb::BoundingBox{D}
  mesh::StructuredBulkMesh{D,T}
end

function cells(m::StructuredBulkMesh)
  StructuredBulkMeshCells(m,num_cells(m))
end


Base.iterate(c::StructuredBulkMeshCells,state=1) = state > c.num_cells ? nothing : (get_cell(c.mesh,state),state+1)

for k in cells(m)
  #  @show k.bb
end

import STLCutter: int_coordinates

function vertex_int_coordinates(m::StructuredBulkMesh{D},n::Integer) where D
  n_d = mutable(VectorValue{D,Int})
  p_d = 1
  p = m.partition .+ 1
  for d in 1:D
    n_d[d] = ( (n-1) ÷ p_d ) % p[d] + 1
    p_d *= p[d]
  end
  n_d.data
end

function get_vertex_id(m::StructuredBulkMesh{D},n::NTuple{D,Int}) where D
  gid = 0
  p_d = 1
  for d in 1:D
    gid += (n[d]-1)*p_d
    p_d *= (m.partition[d]+1)
  end
  gid + 1
end

o = Point(-1.0,0.0,0.0)
s = VectorValue(2.0,1.0,1.0)
p = (2,2,2)
m = StructuredBulkMesh(o,s,p)


vertex_int_coordinates(m,6)
get_vertex_id(m,(3,2,1))

end # module
