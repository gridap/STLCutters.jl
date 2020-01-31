
"""
  const Point{D,T} = VectorValue{D,T}

  Redefine Point s.t.:

    Points{D,T} only support the follwing operations:
       Point - Point = VectorValue
       Point + VectorValue = Point
       VectorValue + Point = Point
"""

const Point = VectorValue

@inline function distance(a::Point{D},b::Point{D}) where D
  norm( b - a )
end
