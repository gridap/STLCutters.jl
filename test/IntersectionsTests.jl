module IntersectionsTests

using Test
using STLCutters
using Gridap
using Gridap.ReferenceFEs

using STLCutters: orthogonal, distance, projection

p = QUAD
K = collect(1:num_vertices(p))
X = collect(get_vertex_coordinates(p))

p1 = Point(0.2,0.1)

@test distance(K,X,p,1,1,p1) == 0.1
@test projection(K,X,p,1,1,p1) == Point(0.2,0.0) 

p1 = Point(0.2,0.0)

@test distance(K,X,p,0,2,p1) == 0.8
@test projection(K,X,p,0,2,p1) == Point(1.0,0.0) 


v = VectorValue(1,2)
n = orthogonal(v)
t = TensorValue(v.data...,n.data...)

@test n == VectorValue(-2,1)
@test v⋅n == 0
@test det(t) ≠ 0

v1 = VectorValue(5,6,3)
v2 = VectorValue(1,2,3)
n = orthogonal(v1,v2)
t = TensorValue(v1.data...,v2.data...,n.data...)

@test n == v1×v2
@test v1⋅n == v2⋅n == 0
@test det(t) ≠ 0

v1 = VectorValue(5,6,3,8)
v2 = VectorValue(1,2,3,4)
v3 = VectorValue(6,9,3,7)
n = orthogonal(v1,v2,v3)
t = TensorValue(v1.data...,v2.data...,v3.data...,n.data...)

@test v1⋅n == v2⋅n == v3⋅n == 0
@test det(t) ≠ 0

p = QUAD
K = collect(1:num_vertices(p))
X = collect(get_vertex_coordinates(p))

p1 = Point(0.2,0.5)
@test have_intersection(K,X,p,p1)

p1 = Point(1.2,0.5)
@test !have_intersection(K,X,p,p1)

p1 = Point(0.5,1.5)
@test !have_intersection(K,X,p,p1)

p = HEX
K = collect(1:num_vertices(p))
X = collect(get_vertex_coordinates(p))

p1 = Point(0.2,0.3,0.5)
@test have_intersection(K,X,p,p1)

p1 = Point(1.2,0.5,0.1)
@test !have_intersection(K,X,p,p1)

p1 = Point(0.5,0.2,-0.1)
@test !have_intersection(K,X,p,p1)

p = QUAD
K = collect(1:num_vertices(p))
X = collect(get_vertex_coordinates(p))

p1 = Point(-0.1,0.5)
p2 = Point(0.5,1.1)
edge = Segment(p1,p2) 

@test have_intersection(K,X,p,edge)

@test intersection_point(K,X,p,2,edge) ≈ Point(0.4,1.0)
@test intersection_point(K,X,p,3,edge) ≈ Point(0.0,0.6)

p1 = Point(-2.1,0.5)
p2 = Point(0.5,2.1)
edge = Segment((p1,p2)) 

@test !have_intersection(K,X,p,edge)

p = HEX
K = collect(1:num_vertices(p))
X = collect(get_vertex_coordinates(p))

p1 = Point(-0.1,0.5,0.8)
p2 = Point(0.5,1.1,0.2)
edge = Segment((p1,p2)) 

@test have_intersection(K,X,p,edge)

@test intersection_point(K,X,p,4,edge) ≈ Point(0.4,1.0,0.3)
@test intersection_point(K,X,p,5,edge) ≈ Point(0.0,0.6,0.7)

end # module
