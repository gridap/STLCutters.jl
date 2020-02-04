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

@test length(s) == sqrt(3)

c = center(s)
@test c.data == (1.5,1.5,1.5)

@test distance(p1,s) == 0

p = Point(0,0,0)
@test distance(s,p) == distance(s,p)
@test distance(p,s) == distance(p,s[1])

p = Point(1.5,1.5,1.5)
@test distance(s,p) < 1e-15

@test have_intersection(s[1],s)

p = Point(0,0,0)
@test !have_intersection(s,p)
@test !have_intersection(s,p)

p = Point(1.5,1.5,1.5)
@test have_intersection(s,p)

p1 = Point(0,0,0)
p2 = Point(1,1,1)
p3 = Point(0.5,0.5,0.5)
p4 = Point(2,2,2)

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

end # module
