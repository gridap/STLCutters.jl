module ConformingSTLsTests

using Test
using STLCutter
using STLCutter: global_dface, local_dface

vertex_coordinates = Point{2,Float64}[(0, 0), (1, 0), (1, 1), (0, 1)]
d_face_to_vertices = [TableOfVectors([[3, 4], [1, 2], [2, 3]])]
facet_normals = Point{2,Float64}[(0, 1), (0, -1), (1, 0)]
d_face_to_facets = [TableOfVectors([[1, 2, 5], [2, 6, 7]])]

stl = ConformingSTL(vertex_coordinates, d_face_to_vertices, facet_normals, d_face_to_facets)

@test stl.vertex_coordinates == vertex_coordinates
@test stl.d_face_to_vertices == d_face_to_vertices
@test stl.facet_normals == facet_normals
@test stl.d_face_to_facets == d_face_to_facets

stl = ConformingSTL(joinpath(@__DIR__,"data/cube.stl"))

outfiles = writevtk(stl,"cube")
rm(outfiles...)

h = HexaCell(Point(0.1,0.1,0.1),Point(1.0,1.0,1.0))
@test have_intersection(h,stl,2,1)

@test global_dface(stl,0,8) == 8
@test global_dface(stl,1,2) == 10
@test global_dface(stl,2,1) == 27

@test local_dface(stl,8) == (0,8)
@test local_dface(stl,10) == (1,2)
@test local_dface(stl,27) == (2,1)

@test local_dface(stl,global_dface(stl,2,5)) == (2,5)
@test global_dface(stl,local_dface(stl,15)...) == 15

using STLCutter: BoundingBox, optimized_compute_cell_to_stl_nfaces


#stl = ConformingSTL("sbunny.stl")
#bb = BoundingBox(stl)
#o = bb.pmin
#s = bb.pmax - bb.pmin
#p = (100,100,100)
#m = StructuredBulkMesh(o,s,p)
#x=optimized_compute_cell_to_stl_nfaces(m,stl);
#@time x=optimized_compute_cell_to_stl_nfaces(m,stl);
#@show @allocated optimized_compute_cell_to_stl_nfaces(m,stl)


end # module
