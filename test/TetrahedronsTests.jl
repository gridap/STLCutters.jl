module TetrahedronsTest

using Test
using LinearAlgebra
using STLCutters
using STLCutters: relative_orientation

t1 = Point(0,0,0)
t2 = Point(1,0,0)
t3 = Point(0,1,0)
t4 = Point(0,0,1)
tet = Tetrahedron(t1,t2,t3,t4)

@test num_dims(tet) == 3
@test num_vertices(tet) == 4
@test num_facets(tet) == 4

@test volume(tet) == 1/6

p = Point(0.2,0.2,0.2)
@test have_intersection(p,tet)
@test have_intersection_point(p,tet)
@test contains_projection(p,tet)

p = Point(2,2,2)
@test !have_intersection(p,tet)

p = Point(0,0,0)
@test have_intersection(p,tet)

p = Point(0.9,0.0,0.0)
@test have_intersection(p,tet)

p = Point(0.0,0.9,0.0)
@test have_intersection(p,tet)

p = Point(0.0,0.0,0.9)
@test have_intersection(p,tet)

@test closest_point(p,tet) == closest_point(tet,p) == p

tri = Triangle(t1,t2,t3)
@test relative_orientation(tri,tet) < 0

tri = Triangle(t2,t3,t4)
@test relative_orientation(tri,tet) > 0

t1 = Point(0,0,0,0)
t2 = Point(1,0,0,0)
t3 = Point(0,1,0,0)
t4 = Point(0,0,1,0)
tet = Tetrahedron(t1,t2,t3,t4)

@test measure(tet) == 1/6

t1 = Point(0,0,0,0)
t2 = Point(2,0,0,0)
t3 = Point(1,4,0,0)
t4 = Point(0,0,5,6)
tet = Tetrahedron(t1,t2,t3,t4)

n = normal(tet)
v1 = t2-t1
v2 = t3-t1
v3 = t4-t1

@test norm(n) ≈ 1
@test abs(v1⋅n) < 1e-14  
@test abs(v2⋅n) < 1e-14 
@test abs(v3⋅n) < 1e-14 

end # module
