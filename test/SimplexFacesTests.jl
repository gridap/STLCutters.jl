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
using STLCutters: signed_distance
using STLCutters: contains_projection
using STLCutters: voxel_intersection
using STLCutters: simplex_face
using STLCutters: min_height
using STLCutters: Plane

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

# Plane{3}

p1 = Point(0.0,1.0,0.0)
p2 = Point(1.0,0.0,0.0)
p3 = Point(0.0,0.0,1.0)
f = simplex_face(p1,p2,p3)

Π = Plane(f)
p = Point(0.0,0.0,0.0)

@test signed_distance(p,Π) ≈ 1/√3

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

end # module

