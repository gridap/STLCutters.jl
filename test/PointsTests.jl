module PointsTests

using Test
using STLCutters

p = Point{0,Int}()

@test num_components(p) == 0
@test component_type(p) == Int

p = Point{0,Int}(())

@test num_components(p) == 0
@test component_type(p) == Int

p = Point{1,Int}((1,))

@test num_components(p) == 1
@test component_type(p) == Int

p = Point((1,))

@test num_components(p) == 1
@test component_type(p) == Int

a = Point(1,1,1)
b = Point(2,3,3)

@test distance(a,b) == 3

p1 = Point(1,2,3)
p2 = Point(0,0,0)
p3 = Point(2,1,0)

@test average(p1) == p1

a = average(p1,p2)
@test isa(a,Point{3})
@test get_data(a) == (0.5,1.0,1.5)
@test a == average((p1,p2))

a = average(p1,p2,p3)
@test isa(a,Point{3})
@test get_data(a) == (1.0,1.0,1.0)
@test a == average((p1,p2,p3))

a = average([p1,p2,p3])
@test isa(a,Point{3})
@test get_data(a) == (1.0,1.0,1.0)
@test a == average((p1,p2,p3))

a = Point(1,1,1)
b = Point(2,3,3)

@test closest_point(a,b) == a
@test closest_point(b,a) == b

a = Point(1,1)
m = mutable(a)

m[1] = 2
b = Point(m)

@test b == Point(2,1)

a = Point(1,2,3)

@test one(a) == Point(1,1,1)
@test zero(a) == Point(0,0,0)


end # module
