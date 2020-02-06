module BulkMeshesTests

using Test
using STLCutter

import STLCutter: find_container, compute_cell_to_stl_nfaces,optimized_compute_cell_to_stl_nfaces

o = Point(0.0,0.0,0.0)
s = VectorValue(1.0,1.0,1.0)
p = (2,2,2)
m = StructuredBulkMesh(o,s,p)

@test num_dims(m) == 3

@test num_cells(m) == 8

@test get_cell(m,2).bb.pmin == Point(0.5,0.0,0.0)
@test get_cell(m,2).bb.pmax == Point(1.0,0.5,0.5)

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

o_c2s = optimized_compute_cell_to_stl_nfaces(m,stl)
c2s = compute_cell_to_stl_nfaces(m,stl)
o_c2s == o_c2s


end # modul
