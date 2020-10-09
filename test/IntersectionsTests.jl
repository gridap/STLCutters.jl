module IntersectionsTests

using Test
using STLCutters
using Gridap
using Gridap.ReferenceFEs

using STLCutters: orthogonal
using STLCutters: distance
using STLCutters: projection
using STLCutters: is_on_boundary 
using STLCutters: contains_projection 
using STLCutters: Cell, Triangle 

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

@test have_intersection(K,X,p,edge,atol=1e-9)

face1 = get_dimrange(p,1)[2]
face2 = get_dimrange(p,1)[3]
@test intersection_point(K,X,p,face1,edge) ≈ Point(0.4,1.0)
@test intersection_point(K,X,p,face2,edge) ≈ Point(0.0,0.6)

p1 = Point(-2.1,0.5)
p2 = Point(0.5,2.1)
edge = Segment((p1,p2)) 

@test !have_intersection(K,X,p,edge,atol=1e-9)

p = HEX
K = collect(1:num_vertices(p))
X = collect(get_vertex_coordinates(p))

p1 = Point(-0.1,0.5,0.8)
p2 = Point(0.5,1.1,0.2)
edge = Segment((p1,p2)) 

@test have_intersection(K,X,p,edge,atol=1e-9)

face1 = get_dimrange(p,2)[4]
face2 = get_dimrange(p,2)[5]
@test intersection_point(K,X,p,face1,edge) ≈ Point(0.4,1.0,0.3)
@test intersection_point(K,X,p,face2,edge) ≈ Point(0.0,0.6,0.7)


p1 = Point(0.0,0.0,0.5)
p2 = Point(1.0,0.0,0.5)
p3 = Point(0.0,1.0,0.5)

t = Triangle(p1,p2,p3)

p1 = Point(0.2,0.2,0.0)
p2 = Point(0.3,0.3,1.0)
s = Segment(p1,p2)
@test have_intersection(t,s)
@test intersection_point(t,s) ≈ Point(0.25,0.25,0.5)

p1 = Point(1.0,0.2,0.0)
p2 = Point(1.0,0.3,1.0)
s = Segment(p1,p2)
@test !have_intersection(t,s)


p1 = Point(0.5,0.5,0.0)
p2 = Point(0.5,0.5,1.0)
s = Segment(p1,p2)
@test have_intersection(t,s)

p1 = Point(1.0,0.0,0.0)
p2 = Point(1.0,0.0,1.0)
s = Segment(p1,p2)
@test have_intersection(t,s)

p1 = Point(0.0,0.0,0.5)
p2 = Point(3.0,0.0,0.5)
p3 = Point(0.0,3.0,0.5)

t = Triangle(p1,p2,p3)

p = HEX
K = collect(1:num_vertices(p))
X = collect(get_vertex_coordinates(p))

@test have_intersection(K,X,p,t,atol=1e-9)

p = HEX
K = collect(1:num_vertices(p))
X = collect(get_vertex_coordinates(p))

p1 = Point(0.0,0.0,0.0)
p2 = Point(1.0,0.0,0.0)
p3 = Point(0.0,1.0,0.0)

t = Triangle(p1,p2,p3)

@test is_on_boundary(K,X,p,t,atol=1e-9)

p = QUAD
K = collect(1:num_vertices(p))
X = collect(get_vertex_coordinates(p))

p1 = Point(0.0,0.0)
p2 = Point(0.5,0.0)

s = Segment(p1,p2)

@test is_on_boundary(K,X,p,s,atol=1e-9)

p = QUAD
X = [
 Point(0.0, 1.0, 0.0)
 Point(1.0, 1.0, 0.0)
 Point(0.0, 1.0, 0.0)
 Point(0.0, 0.0, 0.0) ]
K = collect(1:length(X))
point = Point(0.6, 0.9, 0.0)
c = Cell(K,X,p)

@test contains_projection(c,point)


p = HEX

X = Point{3,Float64}[
(7.465756149291991, -28.433887672424323, 5.310378551483154),
(7.465756149291991, -26.387540672677574, 5.310378551483154),
(7.465756149291991, -28.433887672424323, 6.326480469703668),
(7.465756149291991, -26.387540672677574, 6.326480469703668),
(7.465756149291991, -28.433887672424323, 5.310378551483154),
(13.433560621622705, -28.433887672424323, 5.310378551483154),
(7.465756149291991, -28.433887672424323, 6.326480469703668),
(13.433560621622705, -28.433887672424323, 6.326480469703668) ]
K = collect(1:num_vertices(p))
s = Segment{3,Float64}((
  (1.3716708421707153, -29.638381958007812, 8.035669326782227),
  (12.697721481323242, -28.181570053100586, 5.310378551483154)))
c = Cell(K,X,p)
@test have_intersection(c,s,atol=1e-9)


X = Point{3,Float64}[
  (7.465756149291991, -28.433887672424323, -1.1817026233673147),
  (10.736066634884406, -28.433887672424323, -1.1817026233673147),
  (7.465756149291991, -28.433887672424323, 6.326480469703668),
  (10.736066634884406, -28.433887672424323, 6.326480469703668),
  (7.465756149291991, -28.433887672424323, -1.1817026233673147),
  (7.9181595660308295, -28.796340631149157, -1.1817026233673147),
  (7.465756149291991, -28.433887672424323, 6.326480469703668),
  (7.9181595660308295, -28.796340631149157, 6.326480469703668) ]
K = collect(1:num_vertices(p))
t = Triangle{3,Float64}((
  (7.337845325469971, -24.842544555664062, 15.425774574279785), 
  (1.3716708421707153, -29.638381958007812, 8.035669326782227), 
  (12.697721481323242, -28.181570053100586, 5.310378551483154)))

c = Cell(K,X,p)
@test have_intersection(c,t,atol=1e-9)






end # module

