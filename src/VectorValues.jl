
struct VectorValue{D,T} <: Number
  data::NTuple{D,T}
end

@inline function VectorValue(v::T...) where T
  data = v
  VectorValue(data)
end

@inline function VectorValue{D,T}(v::Vararg{S,D} where S) where {D,T}
  data = convert(NTuple{D,T},v)
  VectorValue(data)
end

function VectorValue{0,T}() where T
  data = ()
  VectorValue{0,T}(data)
end

function VectorValue(m::MutableVectorValue)
  D = num_components(m)
  T = component_type(m)
  convert(VectorValue{D,T},m)
end

function VectorValue{D,T}(m::MutableVectorValue{D}) where {D,T}
  convert(VectorValue{D,T},m)
end

function Base.show(io::IO,p::VectorValue)
  print(io,p.data)
end

function Base.show(io::IO,::MIME"text/plain",p::VectorValue)
  print(io,typeof(p))
  print(io,p.data)
end

@generated function Base.zero(::Type{VectorValue{D,T}}) where {D,T}
  data = join(["zero(T), " for i in 1:D])
  str = "VectorValue{D,T}($data)"
  Meta.parse(str)
end

Base.zero(v::T) where T<:VectorValue = zero(T)

function mutable(::Type{VectorValue{D,T}}) where {D,T}
  v = zero(VectorValue{D,T})
  MutableVectorValue(v.data)
end

mutable(v::VectorValue) = MutableVectorValue(v.data)

function Base.convert(::Type{VectorValue{D,T}},v::MutableVectorValue{D}) where {D,T}
  data = convert(NTuple{D,T},v.data)
  VectorValue(data)
end

function Base.convert(::Type{MutableVectorValue{D,T}},v::VectorValue{D}) where {D,T}
  data = convert(NTuple{D,T},v.data)
  MutableVectorValue(data)
end

function Base.convert(::Type{VectorValue{D,T}},v::NTuple{D}) where {D,T}
  data = convert(NTuple{D,T},v)
  VectorValue(data)
end

function Base.convert(::Type{NTuple{D,T}},v::VectorValue{D}) where {D,T}
  convert(NTuple{D,T},v.data)
end

num_components(::VectorValue{D}) where D = D

num_components(::Type{<:VectorValue{D}}) where D = D

component_type(::VectorValue{D,T}) where {D,T} = T

component_type(::Type{VectorValue{D,T}}) where {D,T} = T

@inline Base.getindex(a::VectorValue,i::Integer) = a.data[i]

@inline Base.isequal(a::VectorValue{D},b::VectorValue{D}) where D = ( a.data == b.data )

@generated function Base.:+(a::VectorValue{D},b::VectorValue{D}) where D
  data = join(["a.data[$i] + b.data[$i]," for i in 1:D])
  str = "VectorValue(($data))"
  Meta.parse(str)
end

@generated function Base.:-(a::VectorValue{D},b::VectorValue{D}) where D
  data = join(["a.data[$i] - b.data[$i]," for i in 1:D])
  str = "VectorValue(($data))"
  Meta.parse(str)
end

@generated function Base.:*(a::VectorValue{D},b::VectorValue{D}) where D
  data = join(["+ a.data[$i] * b.data[$i]" for i in 1:D])
  Meta.parse(data)
end

@generated function Base.:*(v::VectorValue{D},a::Number) where D
  data = join(["v.data[$i] * a," for i in 1:D])
  str = "VectorValue(($data))"
  Meta.parse(str)
end

@inline Base.:*(a::Number,v::VectorValue) = v * a

@generated function Base.:/(v::VectorValue{D},a::Number) where D
  data = join(["v.data[$i] / a ," for i in 1:D])
  str = "VectorValue($data)"
  Meta.parse(str)
end

@inline function LinearAlgebra.dot(a::VectorValue{D},b::VectorValue{D}) where D
  a*b
end

@inline function LinearAlgebra.norm(v::VectorValue)
  √( v ⋅ v )
end

@inline function LinearAlgebra.cross(a::VectorValue,b::VectorValue)
  if !( num_components(a) == num_components(b) == 3)
        throw(DimensionMismatch("cross product only supports VectorValues of 3 components"))
    end
  data = (  (a[2]*b[3] - a[3]*b[2]),
            (a[3]*b[1] - a[1]*b[3]),
            (a[1]*b[2] - a[2]*b[1]))
  VectorValue(data)
end
