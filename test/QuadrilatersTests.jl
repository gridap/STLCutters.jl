module QuadrilaterTests

using STLCutter.Cutter
using Test

b = BoundingBox( Point(0,0), Point(1,1) )
q1 = Quadrilater( b )
q2 = Quadrilater( Point(0,0), Point(1,0), Point(0,1), Point(1,1) )

@test q1 == q2
@test area(q1) == measure(q1) == 1
@test center(q1) == Point(0.5,0.5)

@test BoundingBox(q1) == b

p = Point(0.5,0.5)

@test have_intersection(p,q1)

p = Point(0.5,1.5)

@test !have_intersection(p,q1)

q = Quadrilater( Point(0,0,0), Point(1,0,0), Point(0,1,0), Point(1,1,0) )
p = Point(0.5,0.5,0.5)

@test contains_projection(p,q) 
@test distance(p,q) == 0.5
@test projection(p,q) == Point(0.5,0.5,0.0)
@test center(q) == Point(0.5,0.5,0.0)

q = Quadrilater( Point(0,0,0), Point(1,0,0), Point(0,1,0), Point(1,1,0) )
p = Point(1.5,0.5,0.5)
@test !contains_projection(p,q) 
@test distance(p,q) ≈ 1/√2 

end # module
