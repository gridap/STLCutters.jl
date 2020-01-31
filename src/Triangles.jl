
const num_points_per_triangle = 3
struct Triangle{D}
  p::NTuple{num_points_per_triangle,Point{D}}
end

function Triangle(p1::Point,p2::Point,p3::Point)
  Triangle((p1,p2,p3))
end

@inline Base.getindex(t::Triangle,i::Integer) = t.p[i]

@inline function Base.getindex(t::Triangle,i::NTuple{num_points_per_segment,Integer})
 Segment(t[i[1]],t[i[2]])
end

num_dims(::Triangle{D}) where D = D

@inline vertices(t::Triangle) = t.p

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
  average(t.p)
end

function distance(p::Point{D},t::Triangle{D}) where D
  if D != 3
    throw(DimensionMismatch("distance Point{D} to Triangle{D} only valid for 3 dimensions"))
  end
  o = center(t)
  n = normal(t)
  n = n / norm(n)
  o_p = p - o
  p_projection = o_p ⋅ n
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
