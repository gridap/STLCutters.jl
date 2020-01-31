
"""
  const Point{D,T} = VectorValue{D,T}

  Redefine Point s.t.:

    Points{D,T} only support the follwing operations:
       Point - Point = VectorValue
       Point + VectorValue = Point
       Point - VectorValue = Point
       VectorValue + Point = Point
"""

const Point = VectorValue

@inline function distance(a::Point{D},b::Point{D}) where D
  norm( b - a )
end

function average(p::Point)
  p
end

@generated function average(p1::Point{D},p2::Point{D}) where D
  data = join(["(p1.data[$i] + p2.data[$i])/2.0," for i in 1:D])
  str = "Point(($data))"
  Meta.parse(str)
end

@generated function average(p1::Point{D},p2::Point{D},p3::Point{D}) where D
  data = join(["(p1.data[$i] + p2.data[$i] + p3.data[$i])/3.0," for i in 1:D])
  str = "Point(($data))"
  Meta.parse(str)
end

@inline average(p::NTuple{N,Point{D}}) where {N,D} = average(p...)
