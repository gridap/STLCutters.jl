using LinearAlgebra

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

function distance(p1::Point{D,T},p2::Point{D,T}) where {D,T}
  norm( collect(p2.data) - collect(p1.data) )
end

# VectorField to replace Point in some instances
struct VectorField{D,T} <: Number
  data::NTuple{D,T}
end

function VectorField(data::Vararg{D,T}) where {D,T}
  VectorField(data)
end

function Base.show(io::IO,v::VectorField)
  print(io,v.data)
end

function Base.show(io::IO,::MIME"text/plain",v::VectorField)
  print(io,typeof(v))
  print(io,v.data)
end

Base.length(::Type{<:VectorField{D}}) where D = D

Base.length(::VectorField{D}) where D = D

@inline Base.getindex(v::VectorField,i::Integer) = v.data[i]

Base.convert(::Type{VectorField{D,T}},v::NTuple{D}) where {D,T} = VectorField{D,T}(v)
