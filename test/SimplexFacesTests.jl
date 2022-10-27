module SimplexFacesTests

using Test
using STLCutters
using Gridap
using Gridap.ReferenceFEs

using STLCutters: center
using STLCutters: distance
using STLCutters: projection
using STLCutters: measure
using STLCutters: normal
using STLCutters: origin
using STLCutters: signed_distance
using STLCutters: contains_projection
using STLCutters: expand_face
using STLCutters: has_intersection
using STLCutters: voxel_intersection
using STLCutters: simplex_face
using STLCutters: min_height
using STLCutters: displace
using STLCutters: Plane
using STLCutters: CartesianPlane

# Face{1,2}

p1 = Point(0.0,1.0)
p2 = Point(1.0,0.0)
f = simplex_face(p1,p2)

@test num_dims(f) == 1
@test num_point_dims(f) == 2
@test center(f) == Point(0.5,0.5)
@test measure(f) ≈ √2
@test normal(f) ≈ (√2)/2 * Point(1,1)
@test normal(f,1) ≈ (√2)/2 * Point(-1,1)
@test normal(f,2) ≈ (√2)/2 * Point(1,-1)

pmin,pmax = get_bounding_box(f)

@test measure(pmin,pmax) == 1

# Face{1,3}

p1 = Point(0.0,1.0,0.0)
p2 = Point(1.0,0.0,1.0)
f = simplex_face(p1,p2)

@test num_dims(f) == 1
@test num_point_dims(f) == 3
@test center(f) == Point(0.5,0.5,0.5)
@test measure(f) ≈ √3
@test normal(f,1) ≈ 1/(√3) * Point(-1,1,-1)
@test normal(f,2) ≈ 1/(√3) * Point(1,-1,1)

pmin,pmax = get_bounding_box(f)

@test measure(pmin,pmax) == 1

# Face{2,2}

p1 = Point(0.0,1.0)
p2 = Point(1.0,0.0)
p3 = Point(0.0,0.0)
f = simplex_face(p1,p2,p3)

@test num_dims(f) == 2
@test num_point_dims(f) == 2
@test center(f) ≈ Point(1/3,1/3)
@test measure(f) ≈ 1/2
@test min_height(f) ≈ 1/√2 

pmin,pmax = get_bounding_box(f)

@test measure(pmin,pmax) == 1

# Face{2,3}

p1 = Point(0.0,1.0,0.0)
p2 = Point(1.0,0.0,0.0)
p3 = Point(0.0,0.0,1.0)
f = simplex_face(p1,p2,p3)

@test center(f) ≈ Point(1/3,1/3,1/3)
@test measure(f) ≈ (√3)/2

@test num_dims(f) == 2
@test num_point_dims(f) == 3
@test normal(f) ≈ 1/√3 * Point(-1,-1,-1)
@test min_height(f) ≈ √(3/2) 

pmin,pmax = get_bounding_box(f)

@test measure(pmin,pmax) == 1

# Face{3,3}

p1 = Point(0.0,0.0,0.0)
p2 = Point(1.0,0.0,0.0)
p3 = Point(0.0,1.0,0.0)
p4 = Point(0.0,0.0,1.0)
f = simplex_face(p1,p2,p3,p4)

@test center(f) ≈ Point(1/4,1/4,1/4)
@test measure(f) ≈ 1/6
@test min_height(f) ≈ 1/√3 

pmin,pmax = get_bounding_box(f)

@test measure(pmin,pmax) == 1

# Plane{2}

p1 = Point(0.0,1.0)
p2 = Point(1.0,0.0)
f = simplex_face(p1,p2)

Π = Plane(f)
p = Point(0.0,0.0)

@test signed_distance(p,Π) ≈ -1/√2

x0 = Point(0.0,0.0)
d = 1
or = 1
Π = CartesianPlane(x0,d,or)
p = Point(1.0,1.0)

@test normal(Π) == VectorValue(1,0)
@test origin(Π) == Point(0,0)
@test signed_distance(p,Π) == 1


# Plane{3}

p1 = Point(0.0,1.0,0.0)
p2 = Point(1.0,0.0,0.0)
p3 = Point(0.0,0.0,1.0)
f = simplex_face(p1,p2,p3)

Π = Plane(f)
p = Point(0.0,0.0,0.0)

@test signed_distance(p,Π) ≈ 1/√3

x0 = Point(0.0,0.0,0.0)
d = 1
or = 1
Π = CartesianPlane(x0,d,or)
p = Point(1.0,1.0,1.0)

@test normal(Π) == VectorValue(1,0,0)
@test origin(Π) == Point(0,0,0)
@test signed_distance(p,Π) == 1

# Distance Point - Point

p1 = Point(0.0,0.0)
p2 = Point(1.0,1.0)

distance(p1,p2) ≈ √2

p1 = Point(0.0,0.0,0.0)
p2 = Point(1.0,1.0,1.0)

distance(p1,p2) ≈ √3

# Distance Point - Face{1}

p1 = Point(0.0,1.0)
p2 = Point(1.0,0.0)
f = simplex_face(p1,p2)

a = Point(0.0,0.0)
b = Point(2.0,0.0)
c = Point(0.5,0.5)

@test distance(a,f) ≈ 1/√(2)
@test distance(b,f) ≈ 1.0
@test abs( distance(c,f) ) < 1e-9

# Distance Point - Face{2}

p1 = Point(0.0,1.0,0.0)
p2 = Point(1.0,0.0,0.0)
p3 = Point(0.0,0.0,1.0)
f = simplex_face(p1,p2,p3)

a = Point(0.0,0.0,0.0)
b = Point(2.0,0.0,0.0)
c = Point(1/3,1/3,1/3)

@test distance(a,f) ≈ 1/√(3)
@test distance(b,f) ≈ 1.0
@test abs( distance(c,f) ) < 1e-9

# Distance Point - Face{3}

p1 = Point(0.0,0.0,0.0)
p2 = Point(1.0,0.0,0.0)
p3 = Point(0.0,1.0,0.0)
p4 = Point(0.0,0.0,1.0)
f = simplex_face(p1,p2,p3,p4)

a = Point(0.0,0.0,0.0)
b = Point(2.0,0.0,0.0)
c = Point(1/4,1/4,1/4)
d = Point(1.0,1.0,1.0)

@test abs(distance(a,f)) < 1e-9
@test distance(b,f) ≈ 1.0
@test distance(c,f) == 0
@test distance(d,f) ≈ 2/√3

# Projection Point - Point

a = Point(0.0,1.0)
b = Point(1.0,0.0)

@test projection(a,b) == b

# Projection Point - Face{1}

p1 = Point(0.0,1.0)
p2 = Point(1.0,0.0)
f = simplex_face(p1,p2)

a = Point(0.0,0.0)
b = Point(0.5,0.0)

@test projection(a,f) ≈ Point(0.5,0.5)
@test projection(b,f) ≈ Point(0.75,0.25)

# Projection Point - Face{2}

p1 = Point(0.0,1.0,0.0)
p2 = Point(1.0,0.0,0.0)
p3 = Point(0.0,0.0,1.0)
f = simplex_face(p1,p2,p3)

a = Point(0.0,0.0,0.0)
@test projection(a,f) ≈ Point(1/3,1/3,1/3)

# Contains projection Point - Face{1}

p1 = Point(0.0,1.0)
p2 = Point(1.0,0.0)
f = simplex_face(p1,p2)

a = Point(0.0,0.0)
b = Point(2.0,0.0)

@test contains_projection(f,a)
@test !contains_projection(f,b)

# Contains projection Point - Face{2}

p1 = Point(0.0,1.0,0.0)
p2 = Point(1.0,0.0,0.0)
p3 = Point(0.0,0.0,1.0)
f = simplex_face(p1,p2,p3)

a = Point(0.0,0.0,0.0)
b = Point(2.0,0.0,0.0)

@test contains_projection(f,a)
@test !contains_projection(f,b)

# Contains projection Point - Face{3}

p1 = Point(0.0,0.0,0.0)
p2 = Point(1.0,0.0,0.0)
p3 = Point(0.0,1.0,0.0)
p4 = Point(0.0,0.0,1.0)
f = simplex_face(p1,p2,p3,p4)

a = Point(1/4,1/4,1/4)
b = Point(1.0,1.0,1.0)

@test contains_projection(f,a)
@test !contains_projection(f,b)

# Displace planes

x = Point(1.0,2.0,3.0)
n = VectorValue(0,0,1.0)

plane0 = Plane(x,n)
dist = 0.5
plane = displace(plane0,dist)
@test normal(plane) == normal(plane0) == n
@test origin(plane) == Point(1,2,3.5)

plane = displace(plane0,dist,false)
@test normal(plane) == normal(plane0) == n
@test origin(plane) == Point(1,2,2.5)

x = Point(1,2,3)
d = 3
plane0 = CartesianPlane(x,d,1)
plane = displace(plane0,dist)
dist = 0.5
@test origin(plane)⋅normal(plane) == origin(plane0)⋅normal(plane0)+dist == 3.5

x = Point(1,2,3)
d = 3
plane0 = CartesianPlane(x,d,-1)
plane = displace(plane0,dist)
dist = 0.5
@test origin(plane)⋅normal(plane) == origin(plane0)⋅normal(plane0)+dist == -2.5

x = Point(1,2,3)
d = 3
plane0 = CartesianPlane(x,d,1)
plane = displace(plane0,dist,false)
dist = 0.5
@test origin(plane)⋅normal(plane) == origin(plane0)⋅normal(plane0)-dist == 2.5

x = Point(1,2,3)
d = 3
plane0 = CartesianPlane(x,d,-1)
plane = displace(plane0,dist,false)
dist = 0.5
@test origin(plane)⋅normal(plane) == origin(plane0)⋅normal(plane0)-dist == -3.5

# General intersection

p1 = Point(0.0,0.0,0.0)
p2 = Point(1.0,0.0,0.0)
p3 = Point(0.0,1.0,0.0)
p4 = Point(0.0,0.0,1.0)
f1 = simplex_face(p1,p2,p3,p4)

offset = 0.5
p1 = Point(0.0,0.0,0.0) + offset
p2 = Point(1.0,0.0,0.0) + offset
p3 = Point(0.0,1.0,0.0) + offset
p4 = Point(0.0,0.0,1.0) + offset
f2 = simplex_face(p1,p2,p3,p4)

f3 = expand_face(f1,0.5)
f4 = expand_face(f2,0.5)

@test !has_intersection(f1,f2)
@test has_intersection(f1,f3)
@test has_intersection(f3,f1)
@test !has_intersection(f2,f3)
@test has_intersection(f1,f4)
@test has_intersection(f2,f4)
@test has_intersection(f3,f4)
@test has_intersection(f4,f3)

#writevtk(f1,"f1")
#writevtk(f2,"f2")
#writevtk(f3,"f3")
#writevtk(f4,"f4")

# Voxel intersection Face{2}

p = QUAD

pmin,pmax = Point(0.0,0.0),Point(1.0,1.0)
pmin -= 0.1
pmax += 0.1

p1 = Point(0.5,0.5)
p2 = Point(2.0,5.0)
a = simplex_face(p1,p2)

p1 = Point(0.0,1.0)
p2 = Point(1.0,0.0)
b = simplex_face(p1,p2)

p1 = Point(0.0,2.0)
p2 = Point(2.0,0.0)
c = simplex_face(p1,p2)

p1 = Point(1.0,2.0)
p2 = Point(2.0,1.0)
d = simplex_face(p1,p2)

p1 = Point(1.5,0.5)
p2 = Point(2.0,0.0)
e = simplex_face(p1,p2)

@test voxel_intersection(a,pmin,pmax,p)
@test voxel_intersection(b,pmin,pmax,p)
@test voxel_intersection(c,pmin,pmax,p)
@test !voxel_intersection(d,pmin,pmax,p)
@test !voxel_intersection(e,pmin,pmax,p)


# Voxel intersection Face{3}

p = HEX

pmin,pmax = Point(0.0,0.0,0.0),Point(1.0,1.0,1.0)
pmin -= 0.1
pmax += 0.1

p1 = Point(1.0,0.0,0.0)
p2 = Point(0.0,1.0,0.0)
p3 = Point(0.0,0.0,1.0)
a = simplex_face(p1,p2,p3)

p1 = Point(9.0,-11.0,-11.0)
p2 = Point(-11.0,9.0,9.0)
p3 = Point(10.0,10.0,10.0)
b = simplex_face(p1,p2,p3)

p1 = Point(0.5,0.5,0.5)
p2 = Point(0.0,2.0,0.0)
p3 = Point(0.0,0.0,2.0)
c = simplex_face(p1,p2,p3)

p1 = Point(4.0,0.0,0.0)
p2 = Point(0.0,4.0,0.0)
p3 = Point(0.0,0.0,4.0)
d = simplex_face(p1,p2,p3)

p1 = Point(2.0,-0.2,-0.2)
p2 = Point(1.5,0.5,0.0)
p3 = Point(1.5,0.0,0.5)
e = simplex_face(p1,p2,p3)

@test voxel_intersection(a,pmin,pmax,p)
@test voxel_intersection(b,pmin,pmax,p)
@test voxel_intersection(c,pmin,pmax,p)
@test !voxel_intersection(d,pmin,pmax,p)
@test !voxel_intersection(e,pmin,pmax,p)

# Voxel intersection Cell

t = simplex_face(get_vertex_coordinates(TET)...)
pmin,pmax = get_bounding_box(t)
p0 = center(t)

@test voxel_intersection(t,pmin,pmax,HEX)
@test voxel_intersection(t,pmin,pmin+0.1,HEX)
@test !voxel_intersection(t,pmax,pmax+0.1,HEX)
@test voxel_intersection(t,p0,p0+0.1,HEX)

f, = writevtk(t,"tet")
rm.(f)

end # module

