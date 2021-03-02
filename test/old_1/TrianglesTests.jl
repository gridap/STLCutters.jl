module TrianglesTests

using Test
using LinearAlgebra
using STLCutters
using STLCutters: num_dims, get_edge, get_vertices, relative_orientation


p1 = Point(0,0,0)
p2 = Point(3,0,0)
p3 = Point(0,4,5)

t = Triangle(p1,p2,p3)
n = normal(t)
v1 = p3-p1
v2 = p2-p1

@test norm(n) ≈ 1
@test v1⋅n ≈ 0
@test v2⋅n ≈ 0

p1 = Point(0,0,0)
p2 = Point(1,0,0)
p3 = Point(0,1,0)

t = Triangle(p1,p2,p3)
@test isa(t,Triangle{3})

t = Triangle((p1,p2,p3))
@test num_dims(t) == 3
@test num_edges(t) == 3
@test num_vertices(t) == 3
@test isa(t,Triangle{3})
@test get_vertices(t)[1] == p1
@test num_dims(t) == 3

n = normal(t)
@test get_data(n) == (0,0,1)

e = get_edge(t,1)
@test isa(e,Segment{3})
@test e.vertices == get_vertices(Segment(p1,p2))

p = Point(0.5,0.5,1.0)
@test distance(p,t) ==  distance(t,p) == 1

p1 = Point(0,0,0)
p2 = Point(3,0,0)
p3 = Point(0,3,0)

t = Triangle(p1,p2,p3)
c = center(t)
@test get_data(c) == (1.0,1.0,0.0)

@test area(t) == measure(t) == 4.5


p1 = Point(0,0)
p2 = Point(3,0)
p3 = Point(0,3)
t = Triangle(p1,p2,p3)
@test area(t) == measure(t) == 4.5

p1 = Point(0,0,0)
p2 = Point(3,0,0)
p3 = Point(0,3,0)
p = Point(0.5,0.5,1.0)
t = Triangle(p1,p2,p3)

@test contains_projection(p1,t)
@test contains_projection(t,p1)
@test contains_projection(t,p2)
@test contains_projection(t,p3)
@test contains_projection(t,p)

s1 = Point(0.5,0.5,1.0)
s2 = Point(0.5,0.5,-1.0)

s = Segment(s1,s2)

p1 = Point(0,0,0)
p2 = Point(3,0,0)
p3 = Point(0,3,0)

t = Triangle(p1,p2,p3)

@test have_intersection(s,t)
@test have_intersection(t,s)

@test intersection(s,t) == Point(0.5,0.5,0.0)

p = Point(0.5,0.5,0.5)
@test projection(p,t) == Point(0.5,0.5,0.0)


p = Point(1,1)
p1 = Point(0,0)
p2 = Point(3,0)
p3 = Point(0,3)
t = Triangle(p1,p2,p3)

@test projection(p,t) == p
@test closest_point(p,t) == closest_point(t,p) == p

p = Point(1,1,1)
p1 = Point(0,0,0)
p2 = Point(3,0,0)
p3 = Point(0,3,0)

t = Triangle(p1,p2,p3)

@test closest_point(p,t) == p
@test closest_point(t,p) == projection(p,t)


p1 = Point(0,0)
p2 = Point(1,0)
p3 = Point(0,1)

t = Triangle(p1,p2,p3)
e = Segment(p1,p2)

@test relative_orientation(e,t) < 0

e = Segment(p2,p1)

@test relative_orientation(e,t) > 0


end # module
