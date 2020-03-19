
struct VectorValue{D,T} <: Number
  data::NTuple{D,T}
end

function get_data(v::VectorValue)
  v.data
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

@inline function VectorValue{D}(v::Vararg{T,D}) where {D,T}
  data = v
  VectorValue(data)
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

@generated function Base.one(::Type{VectorValue{D,T}}) where {D,T}
  data = join(["one(T), " for i in 1:D])
  str = "VectorValue{D,T}($data)"
  Meta.parse(str)
end

Base.one(v::T) where T<:VectorValue = zero(T)

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

@inline function dot(a::VectorValue{D},b::VectorValue{D}) where D
  a*b
end

@inline function norm(v::VectorValue)
  √( v ⋅ v )
end

function cross(a::VectorValue,b::VectorValue)
  throw(DimensionMismatch("cross product only supports VectorValues of 3 components"))
end

@inline function cross(a::VectorValue{3},b::VectorValue{3})
  data = (
    (a[2]*b[3] - a[3]*b[2]),
    (a[3]*b[1] - a[1]*b[3]),
    (a[1]*b[2] - a[2]*b[1]) )
  VectorValue(data)
end

@generated function Base.abs(v::VectorValue{D}) where D
  data = join(["abs(v.data[$i])," for i in 1:D])
  str = "VectorValue(($data))"
  Meta.parse(str)
end

function max_dimension(v::VectorValue{D,T}) where {D,T}
  max_value = typemin(T)
  max_d = 0
  for d = 1:D
    if v[d] > max_value
      max_value = v[d]
      max_d = d
    end
  end
  max_d
end

@generated function canonical_vector(::VectorValue{D,T},d::Integer) where {D,T}
  data = join(["convert(T,1.0) * ( $i == d )," for i in 1:D])
  str = "VectorValue(($data))"
  Meta.parse(str)
end

Base.:-(v::VectorValue{D,T}) where {D,T} = VectorValue{D,T}(.-v.data)

function Base.broadcasted(op,a::VectorValue{D}...) where D 
  datas = get_datas(a...)
  VectorValue(broadcast(op,datas...))
end

function get_datas(a)
  (get_data(a),)
end

function get_datas(a,b...)
  a_data = get_data(a)
  b_data = get_datas(b...)
  (a_data, b_data...)
end

function det(a::VectorValue{1})
  a[1]
end

function det(a::VectorValue{2},b::VectorValue{2})
  a[1]*b[2] - a[2]*b[1]
end

@generated function det(a::NTuple{D,VectorValue{D}}) where D
  str = ""
  for i in 1:D
    vectors = " "
    for j in 2:D
      data = ""
      for k in 1:D
        if k != i
          data *= "a[$j][$k],"
        end
      end
      vectors *= "VectorValue($data), "
    end
    if isodd(i)
      str *= " + "
    else
      str *= " - "
    end
    str *= "a[1][$i] * det($vectors)"
    if i != D
      str *= " + \n"
    end
  end
  Meta.parse(str)
end

function det(::NTuple{N,VectorValue}) where N
  throw(ArgumentError("det(::VectorValue{D}...) only valid for D VectorValue{D}'s, i.e., square matrix"))
end

det(a::VectorValue...) = det(a,)

function orthogonal(a::VectorValue{2})
  VectorValue( -a[2], a[1] )
end

@generated function orthogonal(a::NTuple{N,VectorValue{D}}) where {N,D}
  entries = ""
  for i in 1:D
    vectors = " "
    for j in 1:D-1
      data = ""
      for k in 1:D
        if k != i
          data *= "a[$j][$k],"
        end
      end
      vectors *= "VectorValue($data), "
    end
    if iseven(D+i)
      entries *= " + "
    else
      entries *= " - "
    end
    entries *= "det($vectors),\n"
  end
  str = "VectorValue(\n$entries)"
  Meta.parse(str)
end

function orthogonal(a::VectorValue{D}...) where D
  if length(a) != D-1
    throw(ArgumentError("orthogonal(::VectorValue{D}...) only well-defined for D-1 VectorValues{D}'s"))
  end
  orthogonal(a,)
end


