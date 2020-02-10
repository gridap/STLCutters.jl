
struct Triangle{D}
  points::Tuple{Point{D},Point{D},Point{D}}
end

@inline function Triangle(p1::Point,p2::Point,p3::Point)
  Triangle((p1,p2,p3))
end

@inline Base.getindex(t::Triangle,i::Integer) = t.points[i]

@inline function Base.getindex(t::Triangle,i::NTuple{2,Integer})
 Segment(t[i[1]],t[i[2]])
end

num_dims(::Triangle{D}) where D = D

@inline vertices(t::Triangle) = t.points

const num_edges_per_triangle = 3
const ledge_to_triangle_point = ((1,2),(2,3),(3,1))

function get_edge(t::Triangle,i::Integer)
  t[ ledge_to_triangle_point[i] ]
end

function normal(t::Triangle{D}) where D
  v1 = t[2] - t[1]
  v2 = t[3] - t[1]
  v1 × v2
end

function center(t::Triangle)
  average(t.points)
end

const volume_factor_triangle = 1/2
function volume(t::Triangle{D}) where D
  v1 = t[2] - t[1]
  v2 = t[3] - t[1]
  if D == 2
    n1 = VectorValue(v1[2],-v1[1])
    abs( v2 ⋅ n1 ) * volume_factor_triangle
  elseif D == 3
    norm( v1 × v2 ) * volume_factor_triangle
  end
end

function distance(p::Point{D},t::Triangle{D}) where D
  if D != 3
    throw(DimensionMismatch("distance Point{D} to Triangle{D} only valid for 3 dimensions"))
  end
  if have_intersection(p,t)
    d = typemax(Float64)
    for i in 1:num_edges_per_triangle
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

@inline distance(t::Triangle,p::Point) = distance(p,t)

function have_intersection(p::Point{D},t::Triangle{D}) where D
  if D != 3
    throw(DimensionMismatch("have_intersection Point{D} abd Triangle{D} only defined for 3 dimensions"))
  end
  n = normal(t)
  n = n / norm(n)
  for i ∈ 1:num_edges_per_triangle
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
