module TrianglesTests

using Test
using STLCutter
import STLCutter: num_dims, get_edge

p1 = Point(0,0,0)
p2 = Point(1,0,0)
p3 = Point(0,1,0)

t = Triangle(p1,p2,p3)
@test isa(t,Triangle{3})

t = Triangle((p1,p2,p3))
@test isa(t,Triangle{3})
@test t[1] == p1
@test num_dims(t) == 3

n = normal(t)
@test n.data == (0,0,1)

s=t[(2,3)]
@test isa(s,Segment{3})
@test s.p == Segment(p2,p3).p

e = get_edge(t,1)
@test isa(e,Segment{3})
@test e.p == Segment(p1,p2).p


p = Point(0.5,0.5,1.0)
@test distance(p,t) ==  distance(t,p) == 1

p1 = Point(0,0,0)
p2 = Point(3,0,0)
p3 = Point(0,3,0)

t = Triangle(p1,p2,p3)
c = center(t)
@test c.data == (1.0,1.0,0.0)

have_intersection(p1,t)
@test have_intersection(t,p1)
@test have_intersection(t,p2)
@test have_intersection(t,p3)
@test have_intersection(t,p)

end # module
