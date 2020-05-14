
mutable struct MutableVectorValue{D,T} <: Number
  data::NTuple{D,T}
end

function get_data(v::MutableVectorValue)
  v.data
end

num_components(::MutableVectorValue{D}) where D = D

num_components(::Type{<:MutableVectorValue{D}}) where D = D

component_type(::MutableVectorValue{D,T}) where {D,T} = T

component_type(::Type{MutableVectorValue{D,T}}) where {D,T} = T

Base.getindex(a::MutableVectorValue,i::Integer) = a.data[i]

Base.@propagate_inbounds function Base.setindex!(v::MutableVectorValue, val, i::Integer)
  D = num_components(v)
  T = component_type(v)
  Base.@boundscheck  ((0 < i) && (i <= D)) || throw(DomainError("$i is not in range 1:$D"))
  if isbitstype(T)
    ptr = Base.unsafe_convert(Ptr{T}, pointer_from_objref(v))
    GC.@preserve v unsafe_store!(ptr, convert(T, val), convert(Int,i))
  else
    error("setindex!() with non-isbitstype eltype is not supported.")
  end
  return val
end

@generated function Base.zero(::Type{MutableVectorValue{D,T}}) where {D,T}
  data = join(["zero(T), " for i in 1:D])
  str = "MutableVectorValue{D,T}(($data))"
  Meta.parse(str)
end

Base.zero(v::T) where T<:MutableVectorValue = zero(T)
