
struct Segment{D,T}
  vertices::Tuple{Point{D,T},Point{D,T}}
end

@inline function Segment(p1::Point,p2::Point)
  Segment((p1,p2))
end

@inline Base.getindex(s::Segment,i::Integer) = s.vertices[i]

@inline get_vertices(s::Segment) = s.vertices

@inline num_vertices(s::Segment) = 2

@inline Base.length(s::Segment) = distance(s[1],s[2])

@inline measure(s::Segment) = length(s)

num_dims(::Segment{D}) where D = D

function center(s::Segment)
  average( get_vertices(s) )
end

function distance_to_line(p::Point{2},s::Segment{2})
  o = center(s)
  n = normal(s)
  o_p = p - o
  p_projection = o_p ⋅ n
  abs( p_projection )

end

function distance(p::Point{D},s::Segment{D}) where D
  s1_s2 = s[2] - s[1]
  l = norm(s1_s2)
  v = s1_s2 / l
  s1_p = p - s[1]
  s1_projection = ( s1_p ⋅ v ) * v
  s2_projection = s1_projection - s1_s2
  p_projection = s1_projection - s1_p
  if norm( s1_projection ) > l || norm( s2_projection) > l
    min( distance(s[1],p), distance(s[2],p) )
  else
    norm(p_projection)
  end
end

@inline distance(s::Segment{D},p::Point{D}) where D = distance(p,s)

function contains_projection(p::Point{D},s::Segment{D}) where D
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

@inline contains_projection(s::Segment{D},p::Point{D}) where D = contains_projection(p,s)

function distance(s1::Segment{D},s2::Segment{D}) where D
  throw(ArgumentError("distance(::Segment{$D},::Segment{$D}) not implemented"))
end

function distance(s1::Segment{2},s2::Segment{2})
  if have_intersection(s1,s2)
    0.0
  else
    typemax(0.0)
  end
end


function distance(s1::Segment{3},s2::Segment{3})

  v1 = s1[2] - s1[1]
  v2 = s2[2] - s2[1]
  v1 = v1 / norm(v1)
  v2 = v2 / norm(v2)
  n = v1 × v2

  if norm(n) < 1e-14
    return min( distance(s1[1],s2), distance(s1[2],s2), distance(s2[1],s1), distance(s2[2],s1) )
  else
    n = n / norm(n)
  end

  n1 = n × v1
  n1 = n1 / norm(n1)
  c1 = center(s1)
  s1_s2 = s2[2] - s2[1]
  s1_c = c1 - s2[1]
  α = ( n1 ⋅ s1_c ) / ( n1 ⋅ s1_s2 )
  if α < 0 || α > 1
    return min( distance(s2[1],s1), distance(s2[2],s1) )
  end

  n2 = n × v2
  n2 = n2 / norm(n2)
  c2 = center(s2)
  s1_s2 = s1[2] - s1[1]
  s1_c = c2 - s1[1]
  α = ( n2 ⋅ s1_c ) / ( n2 ⋅ s1_s2 )
  if α < 0 || α > 1
    return min( distance(s1[1],s2), distance(s1[2],s2) )
  end

  abs( ( c2 - c1 ) ⋅ n )
end

function normal(s::Segment{2}) 
  v = s[2] - s[1]
  n = VectorValue( -v[2], v[1])
  n / norm(n)
end

function normal(s::Segment{D}) where D
  throw(ArgumentError("normal(::Segment{D}) only defined in 2D"))
end

function have_intersection_point(s1::Segment{2},s2::Segment{2})
  n1 = normal(s1)
  c1 = center(s1)
  if sign( ( s2[1] - c1 ) ⋅ n1 ) == sign( ( s2[2] - c1 ) ⋅ n1 )
    return false
  end
  n2 = normal(s2)
  c2 = center(s2)
  if sign( ( s1[1] - c2 ) ⋅ n2 ) == sign( ( s1[2] - c2 ) ⋅ n2 )
    return false
  end
  true
end

function have_intersection_point(s1::Segment{D},s2::Segment{D}) where D
  throw(ArgumentError("intersection(::Segment{$D},::Segment{$D} not defined, only in 2D"))
end

have_intersection(s1::Segment,s2::Segment) = have_intersection_point(s1,s2)

function intersection(s1::Segment{2},s2::Segment{2})
  @check have_intersection_point(s1,s2) "The provided segments have no intersection"
  v1 = s1[2] - s1[1]
  n2 = normal(s2)
  c2 = center(s2)
  α = abs( ( ( c2 - s1[1] ) ⋅ n2 ) / ( v1 ⋅ n2 ) )
  s1[1] + α * v1
end

function intersection(s1::Segment{D},s2::Segment{D}) where D
  throw(ArgumentError("intersection(::Segment{$D},::Segment{$D} not defined, only in 2D"))
end

function projection(p::Point{D},s::Segment{D}) where D
  c = center(s)
  v = s[2] - s[1]
  v = v / norm(v)
  c + ( ( p - c ) ⋅ v ) * v
end

function closest_point(s1::Segment{3},s2::Segment{3})
  v1 = s1[2] - s1[1]
  v2 = s2[2] - s2[1]
  v1 = v1 / norm(v1)
  v2 = v2 / norm(v2)

  n = v1 × v2
  n = n / norm(n)
  n2 = n × v2
  n2 = n2 / norm(n2)
  c2 = center(s2)
  s1_s2 = s1[2] - s1[1]
  s1_c = c2 - s1[1]
  α = ( n2 ⋅ s1_c ) / ( n2 ⋅ s1_s2 )
  if α ≤ 0
    s1[1]
  elseif α ≥ 1
    s1[2]
  else
    s1[1] + α * s1_s2
  end
end

function closest_point(s1::Segment{D},s2::Segment{D}) where D
  throw(ArgumentError("distance(:Segment{$D},s2::Segment{$D}) not implemented"))
end

function closest_point(s::Segment,p::Point)
  if contains_projection(p,s)
    p = projection(p,s)
  end
  if !contains_projection(p,s)
    min_dist = Inf
    closest_vertex = 0
    for (i,v) in enumerate(get_vertices(s))
      dist = distance(v,p)
      if dist < min_dist
        closest_vertex = i
      end
    end
    p = get_vertices(s)[closest_vertex]
  end
  p
end

function closest_point(s1::Segment{2},s2::Segment{2})
  intersection(s1,s2)
end

function relative_orientation(p::Point{1},s::Segment{1})
  max_distance = 0.0
  _v = get_vertices(s)
  for v in get_vertices(s)
    dist = distance(v,p)
    if dist ≥ max_distance
      max_distance = dist
      _v = v
    end
  end
  _s = Segment(p,_v)
  - measure_sign(_s)
end

relative_orientation(s::Segment,p::Point) = relative_orientation(p,s)

function signed_measure(s::Segment{1})
  v = s[2]-s[1]
  det(v)
end

function measure_sign(s::Segment{1})
  sign(signed_measure(s))
end

function writevtk(b::Segment{D,T},file_base_name) where {D,T}
  vtk_type_id = 3

  points = zeros(T,D,num_vertices(b))
  for (i,v) in enumerate(get_vertices(b)), d in 1:D
      points[d,i] = v[d]
  end

  vtk_type = VTKCellType(vtk_type_id)
  vertices = [1:num_vertices(b);]
  cells = [ MeshCell(vtk_type,vertices) ]

  vtkfile = vtk_grid(file_base_name,points,cells)
  vtk_save(vtkfile)
end

