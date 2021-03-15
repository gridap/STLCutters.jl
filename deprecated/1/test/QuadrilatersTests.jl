module QuadrilaterTests

using Test
using LinearAlgebra
using STLCutters

using STLCutters: relative_orientation

p1 = Point(0,0,0)
p2 = Point(3,0,0)
p3 = Point(0,4,0)
p4 = Point(3,4,0)

q = Quadrilater(p1,p2,p3,p4)
n = normal(q)
v1 = p4-p1
v2 = p3-p1
v3 = p2-p1

@test norm(n) ≈ 1
@test v1⋅n ≈ 0
@test v2⋅n ≈ 0
@test v3⋅n ≈ 0


b = BoundingBox( Point(0,0), Point(1,1) )
q1 = Quadrilater( b )
q2 = Quadrilater( Point(0,0), Point(1,0), Point(0,1), Point(1,1) )

@test num_dims(q1) == 2
@test num_vertices(q1) == 4
@test num_edges(q1) == 4

@test q1 == q2
@test area(q1) == measure(q1) == 1
@test center(q1) == Point(0.5,0.5)

@test BoundingBox(q1) == b

p = Point(0.5,0.5)

@test have_intersection(p,q1)
@test have_intersection(q1,p)

p = Point(0.5,1.5)

@test !have_intersection(p,q1)
@test !have_intersection(q1,p)

q = Quadrilater( Point(0,0,0), Point(1,0,0), Point(0,1,0), Point(1,1,0) )
p = Point(0.5,0.5,0.5)

@test contains_projection(p,q) 
@test contains_projection(q,p) 
@test distance(p,q) == 0.5
@test projection(p,q) == projection(p,q) == closest_point(q,p) == Point(0.5,0.5,0.0)
@test center(q) == Point(0.5,0.5,0.0)

q = Quadrilater( Point(0,0,0), Point(1,0,0), Point(0,1,0), Point(1,1,0) )
p = Point(1.5,0.5,0.5)
@test !contains_projection(p,q) 
@test distance(p,q) ≈ 1/√2 

s = Segment( Point(0,0), Point(0,1) )
q = Quadrilater( Point(0,0), Point(1,0), Point(0,1), Point(1,1) )
@test relative_orientation(s,q) == 1

end # module
