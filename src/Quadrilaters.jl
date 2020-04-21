
struct Quadrilater{D,T} 
  vertices::NTuple{4,Point{D,T}}
end

Quadrilater(v::Vararg{Point,4}) = Quadrilater(v,)

Quadrilater(b::BoundingBox{2}) = Quadrilater(get_vertices(b))

num_dims(::Quadrilater{D}) where D = D

num_edges(::Type{<:Quadrilater}) = 4

num_edges(::T) where T<:Quadrilater = num_edges(T)

get_vertices(q::Quadrilater) = q.vertices

get_edge_to_vertices(::Type{<:Quadrilater}) = ((1,2),(4,3),(3,1),(2,4))

get_edge_to_vertices(::T) where T<:Quadrilater = get_edge_to_vertices(T)

function get_edge(q::Quadrilater,i::Integer)
  edge_to_vertices = get_edge_to_vertices(q)
  lpoints = edge_to_vertices[i]
  vertices = get_vertices(q)
  Segment(vertices[lpoints[1]],vertices[lpoints[2]])
end

function center(q::Quadrilater)
  average(get_vertices(q))
end

function normal(q::Quadrilater{3})
  v1 = q.vertices[2] - q.vertices[1]
  v2 = q.vertices[3] - q.vertices[1]
  v1 × v2
end

function area(q::Quadrilater{D}) where D
  v1 = q.vertices[2] - q.vertices[1]
  v2 = q.vertices[3] - q.vertices[1]
  if D == 2
    n1 = VectorValue(v1[2],-v1[1])
    abs( v2 ⋅ n1 )
  elseif D == 3
    norm( v1 × v2 )
  else
    throw(ArgumentError("Quadrilater{D} area only implemented in 2 and 3 dimensions"))
  end
end

measure(q::Quadrilater) = area(q)

function have_intersection(p::Point{2},q::Quadrilater{2})
  for i in 1:num_edges(q)
    e = get_edge(q,i)
    n = normal(e) 
    n = n / norm(n)
    c = center(e)
    if ( p - c ) ⋅ n < 0
      return false
    end
  end
  true
end

have_intersection(q::Quadrilater,p::Point) = have_intersection(p,q)

function contains_projection(p::Point{3},q::Quadrilater{3})
  n = normal(q)
  n = n / norm(n)
  for i ∈ 1:num_edges(q)
    e = get_edge(q,i)
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

contains_projection(q::Quadrilater,p::Point) = contains_projection(p,q)

function distance(p::Point{3},q::Quadrilater{3})
  if !contains_projection(p,q)
    d = typemax(Float64)
    for i in 1:num_edges(q)
      e = get_edge(q,i)
      d_e = distance(p,e)
      if d_e < d
        d = d_e
      end
    end
  else
    o = center(q)
    n = normal(q)
    n = n / norm(n)
    o_p = p - o
    p_projection = o_p ⋅ n
    d = abs( p_projection )
  end
  d
end

function distance(p::Point{2},q::Quadrilater{2})
  if have_intersection(p,q)
    distance = 0.0
  else
    distance = typemax(0.0)
  end
end

distance(q::Quadrilater,p::Point) = distance(p,q)

function have_intersection_point(s::Segment{3},q::Quadrilater{3})
  n = normal(q)
  c = center(q)
  s1_s2 = s[2] - s[1]
  s1_c = c - s[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  if α < 0 || α > 1 || isnan(α)
    false
  else
    x = s[1] + s1_s2 * α
    contains_projection(x,q)
  end
end

have_intersection(s::Segment{3},q::Quadrilater{3}) = have_intersection_point(s,q)

function intersection(s::Segment{3},q::Quadrilater{3})
  @check have_intersection_point(s,q)
  n = normal(q)
  c = center(q)
  s1_s2 = s[2] - s[1]
  s1_c = c - s[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  s[1] + s1_s2 * α
end

function projection(p::Point{3},q::Quadrilater{3})
  @check contains_projection(p,q)
  c = center(q)
  n = normal(q)
  n = n / norm(n)
  p + ( ( c - p ) ⋅ n ) * n
end

function distance(s::Segment{3},q::Quadrilater{3})
  if have_intersection(s,q)
    0.0
  else
    min_dist = typemax(0.0)
    for i in 1:num_edges(q)
      e = get_edge(q,i)
      dist = distance(s,e)
      if dist < min_dist
        min_dist = dist
      end
    end
    min_dist
  end
end

distance(q::Quadrilater,s::Segment) = distance(s,q)

closest_point(q::Quadrilater{3},p::Point{3}) = projection(p,q)

closest_point(s::Segment{3},q::Quadrilater{3}) = intersection(s,q)

closest_point(q::Quadrilater{3},s::Segment{3}) = closest_point(s,q)

