module VectorValuesTests

using STLCutter
using Test
using LinearAlgebra

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

# TODO test num_components, component_type, getindex


# Operations

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



function Base.show(io::IO,p::VectorValue)
  print(io,p.data)
end

function Base.show(io::IO,::MIME"text/plain",p::VectorValue)
  print(io,typeof(p))
  print(io,p.data)
end

data = (1,3)
p = VectorValue( data )
@test string(p) == "(1, 3)"

data = (1,2,3)
v = VectorValue(data)
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

v = VectorValue(1,2,2)
@test norm(v) == 3

end # module
