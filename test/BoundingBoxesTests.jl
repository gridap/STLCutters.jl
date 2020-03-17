module HexaCellsTests

using Test
using STLCutter
import STLCutter: _fix_triangle

p1 = Point(1,1,1)
p2 = Point(3,3,3)
bb = BoundingBox(Point(1,1,1),Point(3,3,3))

@test have_intersection(Point(1,2,3),BoundingBox(Point(1,1,1),Point(3,3,3)))
@test !have_intersection(Point(1,2,3),BoundingBox(Point(2,2,2),Point(3,3,3)))

bb = BoundingBox(Point(0,0,0),Point(1,1,1))
@test have_intersection(Segment(Point(0.5,0.5,0.5),Point(1.0, 1.0, 1.0)),bb)
@test have_intersection(Segment(Point(0.0,0.0,0.5),Point(0.0, 0.0, 1.0)),bb)
@test !have_intersection(Segment(Point(0,2,2),Point(1, 2, 2)),bb)
@test have_intersection(Segment(Point(0.5,0.5,0.5),Point(1.0,1.0,1.0)),bb)

t = Triangle(Point(0.5,0.5,0.5), Point(0.0,0.0,1.0), Point(0.0,1.0,0.0))
bb = BoundingBox(Point(1,1,1),Point(2,2,2))
t_pos = _fix_triangle(bb,t)
@test normal(t_pos) == abs(normal(t))

t = Triangle(Point(-1.0,-1.0,-1.0), Point(5.0,-1.0,3.0), Point(-1.0,5.0,3.0))
bb = BoundingBox(Point(0,0,0),Point(2,2,2))

@test have_intersection(t,bb)

t = Triangle(Point(0.5,0.5,0.5), Point(0.0,0.0,1.0), Point(0.0,1.0,0.0))
bb = BoundingBox(Point(0,0,0),Point(2,2,2))

@test have_intersection(t,bb)

p = Point(-1e-8,-1e-8,-1e-8)
bb = BoundingBox(Point(0,0,0),Point(1,1,1))

@test have_intersection(p,bb)


t = Triangle(
  Point(2.1,16.7,32.4),
  Point(-7.6,1.5,25.0),
  Point(-8.2,14.5,44.9))
b = BoundingBox(
  Point((-13.0775749206543, 10.547343444824222, 26.705834579467776)), 
  Point((-2.2652954101562557, 19.209839248657232, 37.43181042671203)))

@test have_intersection(t,b)

#bb1 = BoundingBox(Point(0,0,0),Point(1,1,1))
#bb2 = BoundingBox(Point(0.5,0.5,0.5),Point(2.0,2.0,2.0))
#
#@test have_intersection(bb1,bb2)

#bb1 = BoundingBox(Point(0,0,0),Point(1,1,1))
#bb2 = BoundingBox(Point(1.5,1.5,1.5),Point(2.0,2.0,2.0))
#
#@test !have_intersection(bb1,bb2)

end # module
