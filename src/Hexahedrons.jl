
struct Hexahedron{N,D,T}
  vertices::NTuple{N,Point{D,T}}
end

@generated function Hexahedron(b::BoundingBox{D,T}) where {D,T}
  N = 2^D
  d = Dict( 0 => "pmin", 1 => "pmax" )
  v_str = [ "" for i in 1:N ]
  for i in 1:N
    bin = digits(i-1,base=2,pad=D)
    data = join([ "b.$(d[b])[$i]," for (i,b) in enumerate(bin) ])
    v_str[i] = "Point{$D,$T}($data),"
  end
  vertices = join(v_str)
  str = "Hexahedron{$N,$D,$T}(($vertices))"
  Meta.parse(str)
end

function num_vertices(::Type{<:Hexahedron{N,D}}) where {N,D}
  @check N == 2^D
  N
end

num_vertices(::T) where T<:Hexahedron = num_vertices(T)

num_dims(::Type{<:Hexahedron{N,D}}) where {N,D} = D

num_dims(::T) where T<:Hexahedron = num_dims(T)

Base.getindex(c::Hexahedron,i::Int) = c.vertices[i]

get_vertices(c::Hexahedron) = c.vertices

function BoundingBox(c::Hexahedron{N,D,T}) where {N,D,T}
  @check all( get_data(c[1]) .<= get_data( c[N] ) )
  BoundingBox{D,T}( c[1], c[N] )
end

