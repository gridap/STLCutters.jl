module ConformingSTLsTests

using Test
using STLCutter

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

@inline function num_d_faces(stl::ConformingSTL,d::Int)
  length(stl.d_face_to_vertices[d])
end

function get_d_face(stl::ConformingSTL,d::Int,i::Int)
  list = getlist(stl.d_face_to_vertices[d],i)
  tuple(stl.vertex_coordinates[list]...)
end

num_d_faces(stl,1)
get_d_face(stl,2,2)

end # moduleÇ
