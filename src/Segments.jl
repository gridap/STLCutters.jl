
const num_points_per_segment = 2
struct Segment{D}
  p::NTuple{num_points_per_segment,Point{D}}
end

@inline function Segment(p1::Point,p2::Point)
  Segment((p1,p2))
end

@inline Base.getindex(s::Segment,i::Integer) = s.p[i]

@inline get_vertices(s::Segment) = s.p

@inline Base.length(s::Segment) = distance(s[1],s[2])

num_dims(::Segment{D}) where D = D

function distance(p::Point{D},s::Segment{D}) where D
  v = s[2] - s[1]
  l = norm(v)
  v = v / l
  s1_p = p - s[1]
  s1_projection = ( s1_p ⋅ v ) * v
  s2_projection = s1_projection
  p_projection = s1_projection - s1_p
  if norm( s1_projection ) > l || norm( s2_projection) > l
    min( distance(s[1],p), distance(s[2],p) )
  else
    norm(p_projection)
  end
end

@inline distance(s::Segment{D},p::Point{D}) where D = distance(p,s)

function have_intersection(p::Point{D},s::Segment{D}) where D
  s1_s2 = s[2] - s[1]
  s1_s2 = s1_s2 / norm(s1_s2)
  s1_p = p - s[1]
  p_s2 = s[2] - p
  if s1_p ⋅ s1_s2 < 0
    false
  elseif p_s2 ⋅ s1_s2 < 0
    false
  else
    true
  end
end

@inline have_intersection(s::Segment{D},p::Point{D}) where D = have_intersection(p,s)
