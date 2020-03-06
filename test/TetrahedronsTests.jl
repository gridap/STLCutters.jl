
using Test
using STLCutter

t1 = Point(0,0,0)
t2 = Point(1,0,0)
t3 = Point(0,1,0)
t4 = Point(0,0,1)
tet = Tetrahedron(t1,t2,t3,t4)

@test num_dims(tet) == 3

@test volume(tet) == 1/6

p = Point(0.2,0.2,0.2)
@test have_intersection(p,tet)

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
