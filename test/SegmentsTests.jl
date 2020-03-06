module SegmentsTests

using Test
using STLCutter

p1 = Point(1,1,1)
p2 = Point(2,2,2)

s = Segment(p1,p2)
@test isa(s,Segment{3})
@test num_dims(s) == 3

s = Segment((p1,p2))
@test isa(s,Segment{3})
@test num_dims(s) == 3

@test length(s) == measure(s) == sqrt(3)

c = center(s)
@test get_data(c) == (1.5,1.5,1.5)

@test distance(p1,s) == 0

p = Point(0,0,0)
@test distance(s,p) == distance(s,p)
@test distance(p,s) == distance(p,s[1])

p = Point(1.5,1.5,1.5)
@test distance(s,p) < 1e-15

@test contains_projection(s[1],s)

p = Point(0,0,0)
@test !contains_projection(s,p)
@test !contains_projection(s,p)

p = Point(1.5,1.5,1.5)
@test contains_projection(s,p)

p1 = Point(0.0,0.0,0.0)
p2 = Point(1.0,1.0,1.0)
p3 = Point(0.5,0.5,0.5)
p4 = Point(2.0,2.0,2.0)

s1 = Segment(p1,p2)
s2 = Segment(p3,p4)

@test distance(s1,s2) < 1e-13

p1 = Point(0,0,0)
p2 = Point(1,1,0)
p3 = Point(0,1,1)
p4 = Point(1,0,1)

s1 = Segment(p1,p2)
s2 = Segment(p3,p4)

@test distance(s1,s2) == 1.0

p1 = Point(0,0)
p2 = Point(1,1)
p3 = Point(0,1)
p4 = Point(1,0)

s1 = Segment(p1,p2)
s2 = Segment(p3,p4)
@test have_intersection(s1,s2)

x = Point(0.5,0.5)
@test norm( intersection(s1,s2) - x ) < 1e-13
@test norm( intersection(s2,s1) - x ) < 1e-13

p1 = Point(0.0,0.0)
p2 = Point(0.5,0.0)
s1 = Segment(p1,p2)
@test !have_intersection(s1,s2)


p1 = Point(0.0,0.0)
p2 = Point(1.0,0.0)
p = Point(0.5,0.5)
s = Segment(p1,p2)

@test projection(p,s) == Point(0.5,0.0)

p1 = Point(0,0,0)
p2 = Point(1,1,0)
p3 = Point(0,1,1)
p4 = Point(1,0,1)

s1 = Segment(p1,p2)
s2 = Segment(p3,p4)

x = closest_point(s1,s2)
@test x == Point(0.5,0.5,0.0)
@test distance(x,s1) < 1e-13
@test distance(x,s2) == distance(s1,s2)

p1 = Point(0,0)
p2 = Point(1,1)
p3 = Point(0,1)
p4 = Point(1,0)

p = Point(1,2)

s1 = Segment(p1,p2)
s2 = Segment(p3,p4)

@test closest_point(s1,s2) == intersection(s1,s2)
@test closest_point(p,s1) == p
@test closest_point(s1,p) == projection(p,s1)
@test have_intersection(s1,s2)

end # module
