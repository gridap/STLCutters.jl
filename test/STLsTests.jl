module STLsTests

using Test
using STLCutters
import STLCutters: delete_repeated_vertices!, delete_repeated_vertices, flip_normals
import STLCutters: get_facet_to_vertices, get_facet_normals

stl = STL(joinpath(@__DIR__,"data/cube.stl"))
@test num_dims(stl) == 3

stl_1 = delete_repeated_vertices(stl)
@test stl_1 != stl
@test num_vertices(stl_1) == 8
@test num_vertices(stl) == 36
@test maximum(stl.facet_to_vertices) == 36
@test maximum(stl_1.facet_to_vertices) == 8


stl_2 = delete_repeated_vertices!(stl)
@test stl_2 === stl
@test num_vertices(stl) == 8
@test stl_1 == stl_2
@test stl_1 !== stl_2


vertices = [
  Point(0,0),
  Point(0,1),
  Point(1,1),
  Point(1,0) ]
f_to_v = Table(
[ 1 2 3 4;
  2 3 4 1 ] )

stl = STL(vertices,f_to_v)

@test num_dims(stl) == 2
@test num_facets(stl) == 4
@test num_vertices(stl) == 4


circ = closed_polyline(vertices)

@test get_vertex_coordinates(stl) == get_vertex_coordinates(circ)

@test get_facet_to_vertices(stl) == get_facet_to_vertices(circ)

@test get_facet_normals(stl) == get_facet_normals(circ)


naca = closed_polyline(joinpath(@__DIR__,"data/naca.dat"))

@test num_dims(naca) == 2

@test num_vertices(naca) == num_facets(naca)

@test get_facet_normals(naca) == -get_facet_normals( flip_normals(naca) )




end # module
