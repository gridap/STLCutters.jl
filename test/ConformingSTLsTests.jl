module ConformingSTLsTests

using Test
using STLCutter

vertex_coordinates = Point{2,Float64}[ (0,0), (1,0), (1,1), (0,1) ]
facet_vertices = [ (3,4), (1,2), (2,3) ]
facet_normals = Point{2,Float64}[ (0,1), (0,-1), (1,0)]
edge_vertices = facet_vertices

stl = ConformingSTL(vertex_coordinates,edge_vertices,facet_vertices,facet_normals)

@test stl.vertex_coordinates == vertex_coordinates
@test stl.facet_vertices == facet_vertices
@test stl.facet_normals == facet_normals

end # module
