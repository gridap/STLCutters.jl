
"""
  const Point{D,T} = VectorValue{D,T}

  Redefine Point s.t.:

    Points{D,T} only support the follwing operations:
       Point - Point = VectorValue
       Point + VectorValue = Point
       Point - VectorValue = Point
       VectorValue + Point = Point
"""

struct Point{D,T}
  vector::VectorValue{D,T}
end

function get_data(p::Point)
  get_data(p.vector)
end

@inline function Point(data::NTuple{D,T}) where {D,T}
  Point{D,T}(VectorValue(data))
end

@inline function Point(v::T...) where T
  data = v
  Point(VectorValue(data))
end

@inline function Point{D,T}(v::Vararg{S,D} where S) where {D,T}
  data = convert(NTuple{D,T},v)
  Point(VectorValue(data))
end

@inline function Point{D}(v::Vararg{T,D}) where {D,T}
  data = v
  Point(VectorValue(data))
end

function Point{0,T}() where T
  data = ()
  Point(VectorValue{0,T}(data))
end

function Point{D}(v::VectorValue{D,T}) where {D,T}
  Point{D,T}(v)
end

function Point(v::MutableVectorValue)
  Point(VectorValue(v))
end

function Point{D}(v::MutableVectorValue{D}) where D
  Point{D}(VectorValue(v))
end

function Point{D,T}(v::MutableVectorValue{D,T}) where {D,T}
  Point{D,T}(VectorValue(v))
end

VectorValue(p::Point) = p.vector

VectorValue{D}(p::Point{D}) where D = p.vector

VectorValue{D,T}(p::Point{D}) where {D,T} = convert(Vectorvalue{D,T},p.vector)

function Base.convert(::Type{Point{D,T}},v::VectorValue{D}) where {D,T}
  w = convert(VectorValue{D,T},v)
  Point{D,T}(w)
end

function Base.convert(::Type{Point{D,T}},v::NTuple{D}) where {D,T}
  w = convert(NTuple{D,T},v)
  Point{D,T}(w)
end

function Base.convert(::Type{Point{D,T}},v::MutableVectorValue{D}) where {D,T}
  w = convert(VectorValue{D,T},v)
  Point{D,T}(w)
end

function Base.convert(::Type{VectorValue{D,T}},v::Point{D}) where {D,T}
  w = convert(Point{D,T},v)
  VectorValue{D,T}(w)
end

num_components(p::Point{D}) where D = D
num_components(::Type{<:Point{D}}) where D = D

component_type(p::Point{D,T}) where {D,T} = T
component_type(::Type{Point{D,T}}) where {D,T} = T

function Base.getindex(p::Point,i::Integer)
  p.vector[i]
end

function Base.:+(a::Point{D},b::Point{D}) where D
  a.vector + b.vector
end

function Base.:-(a::Point{D},b::Point{D}) where D
  a.vector - b.vector
end

function Base.:+(a::Point{D},b::VectorValue{D}) where D
  Point(a.vector + b)
end

function Base.:-(a::Point{D},b::VectorValue{D}) where D
  Point(a.vector - b)
end

function Base.:+(a::VectorValue{D},b::Point{D}) where D
  Point(a + b.vector)
end

function Base.:-(a::VectorValue{D},b::Point{D}) where D
  Point(a - b.vector)
end

@inline function distance(a::Point{D},b::Point{D}) where D
  norm( b - a )
end

function average(p::Point)
  p
end

@generated function average(p1::Point{D},p2::Point{D}) where D
  data = join(["(p1[$i] + p2[$i])/2," for i in 1:D])
  str = "Point(($data))"
  Meta.parse(str)
end

@generated function average(p1::Point{D},p2::Point{D},p3::Point{D}) where D
  data = join(["(p1[$i] + p2[$i] + p3[$i])/3," for i in 1:D])
  str = "Point(($data))"
  Meta.parse(str)
end

@inline average(p::NTuple{N,Point{D}}) where {N,D} = average(p...)

function average(ps::AbstractArray{Point{D,T}}) where {D,T}
  r = zero(Point{D,T})
  m = (1/length(ps))
  for p in ps
    r += m*p.vector
  end
  r
end

function Base.broadcasted(op,a::Point{D}...) where D 
  datas = get_datas(a...)
  Point(broadcast(op,datas...))
end

function mutable(p::Point{D,T}) where {D,T}
  mutable(p.vector)
end

function mutable(::Type{Point{D,T}}) where {D,T}
  mutable(VectorValue{D,T})
end

function Base.zero(p::Point{D,T}) where {D,T}
  Point(zero(p.vector))
end

function Base.zero(::Type{Point{D,T}}) where {D,T}
  Point(zero(VectorValue{D,T}))
end

