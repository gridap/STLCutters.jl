module HexaCellsTests

using Test
using STLCutter
import STLCutter: BoundingBox, positivize_normal

p1 = Point(1,1,1)
p2 = Point(3,3,3)
bb = BoundingBox(Point(1,1,1),Point(3,3,3))

@test have_intersection(Point(1,2,3),BoundingBox(Point(1,1,1),Point(3,3,3)))
@test !have_intersection(Point(1,2,3),BoundingBox(Point(2,2,2),Point(3,3,3)))
@test have_intersection(Point(1,2,3),HexaCell(Point(1,1,1),Point(3,3,3)))
@test !have_intersection(Point(1,2,3),HexaCell(Point(2,2,2),Point(3,3,3)))

@test have_intersection(Segment(Point(0.5,0.5,0.5),Point(1, 1, 1)),BoundingBox(Point(0,0,0),Point(1,1,1)))
@test have_intersection(Segment(Point(0.0,0.0,0.5),Point(0, 0, 1)),BoundingBox(Point(0,0,0),Point(1,1,1)))
@test !have_intersection(Segment(Point(0,2,2),Point(1, 2, 2)),BoundingBox(Point(0,0,0),Point(1,1,1)))
@test have_intersection(Segment(Point(0.5,0.5,0.5),Point(1, 1, 1)),HexaCell(Point(0,0,0),Point(1,1,1)))

t=Triangle(Point(0.5,0.5,0.5), Point(0,0,1), Point(0,1,0))
bb=BoundingBox(Point(1,1,1),Point(2,2,2))
t_pos=positivize_normal(bb,t)
@test normal(t_pos) == abs(normal(t))

t=Triangle(Point(-1,-1,-1), Point(5.0,-1.0,3.0), Point(-1,5,3))
bb=BoundingBox(Point(0,0,0),Point(2,2,2))
hex=HexaCell(Point(0,0,0),Point(2,2,2))

@test have_intersection(t,bb)
@test have_intersection(t,hex)

t=Triangle(Point(0.5,0.5,0.5), Point(0,0,1), Point(0,1,0))
bb=BoundingBox(Point(0,0,0),Point(2,2,2))

@test have_intersection(t,bb)

p=Point(-1e-8,-1e-8,-1e-8)
bb=BoundingBox(Point(0,0,0),Point(1,1,1))

@test have_intersection(p,bb)

end # module
