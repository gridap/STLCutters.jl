module MutableVectorValuesTests

using Test
using STLCutter
using STLCutter: MutableVectorValue

data = (1,2,3)
m = MutableVectorValue(data)

@test num_components(m) == 3
@test num_components(typeof(m)) == 3
@test component_type(m) == Int
@test component_type(m) == Int

m[2] = -1.0
@test m.data == (1,-1,3)
@test m[2] == -1

m = mutable(VectorValue{2,Int})
@test isa(m,MutableVectorValue{2,Int})
@test m.data == (0,0)

v = VectorValue(1,2,3)
m = mutable(v)
@test isa(m,MutableVectorValue{3,Int})
@test m.data == v.data

v = VectorValue(1,2,3)
m = convert(MutableVectorValue{3,Int},v)
@test isa(m,MutableVectorValue{3,Int})
@test m.data == v.data

m = convert(MutableVectorValue{3,Float64},v)
@test isa(m,MutableVectorValue{3,Float64})
@test m.data == v.data

m = MutableVectorValue((1,2,3))
v = convert(VectorValue{3,Int},m)
@test isa(v,VectorValue{3,Int})
@test m.data == v.data

m = MutableVectorValue((1,2,3))
v = convert(VectorValue{3,Float64},m)
@test isa(v,VectorValue{3,Float64})
@test m.data == v.data

m = MutableVectorValue((1,2,3))
v = VectorValue(m)
@test isa(v,VectorValue{3,Int})
@test m.data == v.data

m = MutableVectorValue((1,2,3))
v = VectorValue{3,Float64}(m)
@test isa(v,VectorValue{3,Float64})
@test m.data == v.data

end # module
