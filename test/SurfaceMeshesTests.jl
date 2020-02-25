module SurfaceMeshesTests

using Test
using STLCutter
using STLCutter: global_dface, local_dface, dface_dimension

#vertex_coordinates = Point{2,Float64}[(0, 0), (1, 0), (1, 1), (0, 1)]
#d_face_to_vertices = [TableOfVectors([[3, 4], [1, 2], [2, 3]])]
#facet_normals = Point{2,Float64}[(0, 1), (0, -1), (1, 0)]
#d_face_to_facets = [TableOfVectors([[1, 2, 5], [2, 6, 7]])]
#
#s_mesh = SurfaceMesh(vertex_coordinates, d_face_to_vertices, facet_normals, d_face_to_facets)
#
#@test s_mesh.vertex_coordinates == vertex_coordinates
#@test s_mesh.d_face_to_vertices == d_face_to_vertices
#@test s_mesh.facet_normals == facet_normals
#@test s_mesh.d_face_to_facets == d_face_to_facets
#

stl = STL(joinpath(@__DIR__,"data/cube.stl"))

s_mesh = SurfaceMesh(stl)

outfiles = writevtk(s_mesh,"cube")
rm(outfiles...)

box = BoundingBox(Point(0.1,0.1,0.1),Point(1.0,1.0,1.0))
@test have_intersection(box,s_mesh,2,1)

@test global_dface(s_mesh,0,8) == 8
@test global_dface(s_mesh,1,2) == 10
@test global_dface(s_mesh,2,1) == 27

@test dface_dimension(s_mesh,8) == 0
@test dface_dimension(s_mesh,10) == 1
@test dface_dimension(s_mesh,27) == 2


@test local_dface(s_mesh,8,0) == 8
@test local_dface(s_mesh,10,1) == 2
@test local_dface(s_mesh,27,2) == 1


#using STLCutter: BoundingBox, optimized_compute_cell_to_s_mesh_nfaces


#s_mesh = SurfaceMesh("sbunny.stl")
#bb = BoundingBox(s_mesh)
#o = bb.pmin
#s = bb.pmax - bb.pmin
#p = (100,100,100)
#m = StructuredBulkMesh(o,s,p)
#x=optimized_compute_cell_to_s_mesh_nfaces(m,stl);
#@time x=optimized_compute_cell_to_s_mesh_nfaces(m,stl);
#@show @allocated optimized_compute_cell_to_s_mesh_nfaces(m,stl)


end # module
