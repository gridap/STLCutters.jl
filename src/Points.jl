
struct Point{D,T} <: Number
  data::NTuple{D,T}
end

function Point(data::Vararg{D,T}) where {D,T}
  Point(data)
end

function Base.show(io::IO,p::Point)
  print(io,p.data)
end

function Base.show(io::IO,::MIME"text/plain",p::Point)
  print(io,typeof(p))
  print(io,p.data)
end

Base.length(::Type{<:Point{D}}) where D = D

Base.length(::Point{D}) where D = D

@inline Base.getindex(p::Point,i::Integer) = p.data[i]

Base.convert(::Type{Point{D,T}},p::NTuple{D}) where {D,T} = Point{D,T}(p)


