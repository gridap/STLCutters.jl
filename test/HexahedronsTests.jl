module HexahedronsTests

using STLCutters
using Test

using STLCutters: relative_orientation

b = BoundingBox( Point(0,0,0), Point(1,1,1) )

h1 = Hexahedron( b )
h2 = Hexahedron(
  Point(0,0,0),
  Point(1,0,0),
  Point(0,1,0),
  Point(1,1,0),
  Point(0,0,1),
  Point(1,0,1),
  Point(0,1,1),
  Point(1,1,1) )

@test num_dims(h1) == 3
@test num_facets(h1) == 6
@test num_vertices(h1) == 8

@test h1 == h2
@test BoundingBox(h1) == b
@test volume(h1) == measure(h1) == 1
p = Point(0.5,0.5,0.5)

@test have_intersection(h1,p)
@test have_intersection_point(h1,p)
@test contains_projection(h1,p)
@test projection(h1,p) == projection(p,h1) == p
@test closest_point(h1,p) == p

p = Point(0.5,0.5,1.5)

@test !have_intersection(h1,p)
@test distance(p,h1) == Inf


t = Triangle( Point(0,0,0), Point(0,1,0), Point(1,1,0) )
@test relative_orientation(t,h1) == 1

end # module
