
struct Triangle{D,T}
  points::Tuple{Point{D,T},Point{D,T},Point{D,T}}
end

@inline function Triangle(p1::Point,p2::Point,p3::Point)
  Triangle((p1,p2,p3))
end

num_dims(::Triangle{D}) where D = D

num_edges(::Type{<:Triangle}) = 3

num_edges(::T) where T<:Triangle = num_edges(T)

num_vertices(::Type{<:Triangle}) = 2

num_vertices(::T) where T<:Triangle = num_vertices(T)

@inline vertices(t::Triangle) = t.points

get_edge_to_vertices(::Type{<:Triangle}) = ((1,2),(2,3),(3,1))

get_edge_to_vertices(::T) where T<:Triangle = get_edge_to_vertices(T)

function get_edge(t::Triangle,i::Integer)
  edge_to_vertices = get_edge_to_vertices(t)
  lpoints = edge_to_vertices[i]
  Segment(t.points[lpoints[1]],t.points[lpoints[2]])
end

function normal(t::Triangle{D}) where D
  throw(ArgumentError("normal(::Triangle{$D}) not defined, only defined in 3D"))
end

function normal(t::Triangle{3})
  v1 = t.points[2] - t.points[1]
  v2 = t.points[3] - t.points[1]
  v1 × v2
end

function center(t::Triangle)
  average(t.points)
end

function volume(t::Triangle{D}) where D
  factor = 1/2
  v1 = t.points[2] - t.points[1]
  v2 = t.points[3] - t.points[1]
  if D == 2
    n1 = VectorValue(v1[2],-v1[1])
    abs( v2 ⋅ n1 ) * factor
  elseif D == 3
    norm( v1 × v2 ) * factor
  else
    throw(DimensionMismatch("volume of Triangle{D} (area) only defined in 2 and 3 dimensions"))
  end
end

function orientation(t::Triangle{2})
  v1 = t.points[2] - t.points[1]
  v2 = t.points[3] - t.points[1]
  sign( v1[1]*v2[2] - v1[2]*v2[1] )
end

function orientation(t::Triangle{D}) where D
  throw(DimensionMismatch("orientation of Triangle{D} only defined in 2"))
end

function distance(p::Point{3},t::Triangle{3})
  if !have_intersection(p,t)
    d = typemax(Float64)
    for i in 1:num_edges(t)
      e = get_edge(t,i)
      d_e = distance(p,e)
      if d_e < d
        d = d_e
      end
    end
  else
    o = center(t)
    n = normal(t)
    n = n / norm(n)
    o_p = p - o
    p_projection = o_p ⋅ n
    d = abs( p_projection )
  end
  d
end

function distance(p::Point{D},t::Triangle{D}) where D
  throw(DimensionMismatch("distance Point{D} to Triangle{D} only valid for 3 dimensions"))
end

@inline distance(t::Triangle,p::Point) = distance(p,t)

function have_intersection(p::Point{2},t::Triangle{2})
  o = orientation(t)
  o != 0 || throw(ErrorException("Triangle has no orientation"))
  for i ∈ 1:num_edges(t)
    e = get_edge(t,i)
    n = normal(e) * o
    c = center(e)
    if ( p - c ) ⋅ n < 0
      return false
    end
  end
  true
end

function have_intersection(p::Point{3},t::Triangle{3})
  n = normal(t)
  n = n / norm(n)
  for i ∈ 1:num_edges(t)
    e = get_edge(t,i)
    v_e = e[2] - e[1]
    c_e = center(e)
    v_e = v_e / norm(v_e)
    n_e = n × v_e
    if ( p - c_e ) ⋅ n_e < 0
      return false
    end
  end
  true
end

function have_intersection(p::Point{D},t::Triangle{D}) where D
  throw(DimensionMismatch("have_intersection Point{D} and Triangle{D} only defined for 2 and 3 dimensions"))
end

@inline have_intersection(t::Triangle,p::Point) = have_intersection(p,t)

function have_intersection(s::Segment{D},t::Triangle{D}) where D
  n = normal(t)
  c = center(t)
  s1_s2 = s[2] - s[1]
  s1_c = c - s[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  if α < 0 || α > 1
    false
  else
    x = s[1] + s1_s2 * α
    have_intersection(x,t)
  end
end

@inline have_intersection(t::Triangle,s::Segment) = have_intersection(s,t)

function intersection(s::Segment{D},t::Triangle{D}) where D
  @check have_intersection(s,t)
  n = normal(t)
  c = center(t)
  s1_s2 = s[2] - s[1]
  s1_c = c - s[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  s[1] + s1_s2 * α
end

function projection(p::Point{D},t::Triangle{D}) where D
  c = center(t)
  n = normal(t)
  n = n / norm(n)
  p + ( ( c - p ) ⋅ n ) * n
end
