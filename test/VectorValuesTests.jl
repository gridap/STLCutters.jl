module VectorValuesTests

using STLCutter
using Test
using LinearAlgebra
using STLCutter: max_dimension, canonical_vector

data = (1,2,3)
v = VectorValue(data)
@test isa(v,VectorValue{3,Int})
@test v.data == data

data = ()
v = VectorValue{0,Int}(data)
@test isa(v,VectorValue{0,Int})
@test v.data == data

data = (1,2,3)
v = VectorValue(data...)
@test isa(v,VectorValue{3,Int})
@test v.data == data

data = ()
v = VectorValue{0,Int}()
@test isa(v,VectorValue{0,Int})
@test v.data == data

data = (1,2,3)
v = VectorValue{3,Float64}(data...)
@test isa(v,VectorValue{3,Float64})
@test v.data == data

v = zero(VectorValue{2,Int})
@test isa(v,VectorValue{2,Int})
@test v.data == (0,0)

v = zero(v)
@test isa(v,VectorValue{2,Int})
@test v.data == (0,0)

v = zero(VectorValue{0,Int})
@test isa(v,VectorValue{0,Int})
@test v.data == ()

v = zero(v)
@test isa(v,VectorValue{0,Int})
@test v.data == ()

data = (1,3)
v = VectorValue( data )
@test string(v) == "(1, 3)"
@test v[2] == 3

data = (1,2,3)
v = VectorValue(data)
@test v[1] == 1
@test num_components(v) == 3
@test component_type(v) == Int

data = ()
v = VectorValue{0,Int}()
@test num_components(v) == 0
@test component_type(v) == Int

data = ()
v = VectorValue{0,Float64}()
@test num_components(v) == 0
@test component_type(v) == Float64

# Operations

u = VectorValue(1,1,1)
v = VectorValue(1,2,3)

@test u == u
@test u != v

v = VectorValue(1,2)
u = VectorValue(1.0,3.0)

w = v + u
@test w.data == (2.0,5.0)

w = v - u
@test w.data == (0.0,-1.0)

w = v*u
@test w == 7.0

w = dot(v,u)
@test w == 7.0

w = v ⋅ u
@test w == 7.0

v = VectorValue(1,2,3)

u = v * 1
@test u.data == v.data

u = v * 2
@test u.data == (2,4,6)

u = 2 * v
@test u.data == (2,4,6)

u = v / 2
@test u.data == (0.5,1.0,1.5)

v = VectorValue(1,2,2)
@test norm(v) == 3

v = VectorValue(1,2,3)
u = VectorValue(1,1,1)

w = cross(v,v)
@test w.data == (0,0,0)

w = v × v
@test w.data == (0,0,0)

w = u × v
@test w.data == (1,-2,1)

a = VectorValue(1,2,3)
b = VectorValue(4,5,6)

c = a .* b
@test c.data == (4,10,18)

v = VectorValue(0,1,1)
@test max_dimension(v) == 2

x = canonical_vector(v,3)
@test x.data == (0,0,1)

u = VectorValue(1,1,2)
v = VectorValue(0,2,3)
w = VectorValue(2,1,4)

@test min.(u) == max.(u) == u
@test min.(u,v).data == (0,1,2)
@test max.(u,v).data == (1,2,3)
@test min.(u,v,w).data == (0,1,2)
@test max.(u,v,w).data == (2,2,4)

v1 = VectorValue(1,2)
@test orthogonal(v1) ⋅ v1 == 0

v1 = VectorValue(1,2)
v2 = VectorValue(3,4)
A = [ 1 2; 3 4 ]
@test det(v1,v2) == det(A)


v1 = VectorValue(1,2,3)
v2 = VectorValue(4,5,6)
v3 = VectorValue(1,1,2)
A = [ 1 2 3; 4 5 6; 1 1 2 ]
@test orthogonal(v1,v2) == v1 × v2 
@test det(v1,v2,v3) == (v1×v2)⋅v3 == det(A)


v1 = VectorValue(1,0,0,0)
v2 = VectorValue(0,1,0,0)
v3 = VectorValue(0,0,1,0)
v4 = orthogonal(v1,v2,v3)
@test v4 == VectorValue(0,0,0,1)

v1 = VectorValue(1,2,3,4)
v2 = VectorValue(5,6,7,8)
v3 = VectorValue(6,2,9,9)
v4 = orthogonal(v1,v2,v3)
@test v1 ⋅ v4 == 0 
@test v2 ⋅ v4 == 0
@test v3 ⋅ v4 == 0

v1 = VectorValue(1,2,3,4)
v2 = VectorValue(5,6,7,8)
v3 = VectorValue(6,2,9,9)
v4 = VectorValue(5,7,8,5)
A = [
1 2 3 4; 
5 6 7 8;
6 2 9 9;
5 7 8 5 ]
@test det(v1,v2,v3,v4) ≈ det(A)

end # module
