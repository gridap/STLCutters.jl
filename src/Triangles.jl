
struct Triangle{D,T}
  vertices::Tuple{Point{D,T},Point{D,T},Point{D,T}}
end

@inline function Triangle(p1::Point,p2::Point,p3::Point)
  Triangle((p1,p2,p3))
end

function Triangle(p::Point,s::Segment)
  Triangle(p,get_vertices(s)...)
end

function Triangle(s::Segment,p::Point)
  Triangle(get_vertices(s)...,p)
end

num_dims(::Triangle{D}) where D = D

num_edges(::Type{<:Triangle}) = 3

num_edges(::T) where T<:Triangle = num_edges(T)

num_vertices(::Type{<:Triangle}) = 3

num_vertices(::T) where T<:Triangle = num_vertices(T)

@inline get_vertices(t::Triangle) = t.vertices

get_edge_to_vertices(::Type{<:Triangle}) = ((1,2),(2,3),(3,1))

get_edge_to_vertices(::T) where T<:Triangle = get_edge_to_vertices(T)

function get_edge(t::Triangle,i::Integer)
  edge_to_vertices = get_edge_to_vertices(t)
  lpoints = edge_to_vertices[i]
  vertices = get_vertices(t)
  Segment(vertices[lpoints[1]],vertices[lpoints[2]])
end

function normal(t::Triangle{D}) where D
  throw(ArgumentError("normal(::Triangle{$D}) not defined, only defined in 3D"))
end

function normal(t::Triangle{3})
  v1 = t.vertices[2] - t.vertices[1]
  v2 = t.vertices[3] - t.vertices[1]
  v1 × v2
end

function center(t::Triangle)
  average(get_vertices(t))
end

function area(t::Triangle{D}) where D
  factor = 1/2
  v1 = t.vertices[2] - t.vertices[1]
  v2 = t.vertices[3] - t.vertices[1]
  if D == 2
    n1 = VectorValue(v1[2],-v1[1])
    abs( v2 ⋅ n1 ) * factor
  elseif D == 3
    norm( v1 × v2 ) * factor
  else
    throw(ArgumentError("Triangle{D} area only implemented in 2 and 3 dimensions"))
  end
end

measure(t::Triangle) = area(t)

function measure_sign(t::Triangle{2})
  sign( signed_measure(t) )
end

function signed_measure(t::Triangle{2})
  factor = 1/2
  v1 = t.vertices[2] - t.vertices[1]
  v2 = t.vertices[3] - t.vertices[1]
  det(v1,v2) * factor
end

function measure_sign(t::Triangle{D}) where D
  throw(ArgumentError("measure_sign(::Triangle{$D}) not defined, only in 2D"))
end

function distance_to_plane(p::Point{3},t::Triangle{3})
  o = center(t)
  n = normal(t)
  n = n / norm(n)
  o_p = p - o
  p_projection = o_p ⋅ n
  abs( p_projection )
end

function distance(p::Point{3},t::Triangle{3})
  if !contains_projection(p,t)
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

function distance(p::Point{2},t::Triangle{2})
  if have_intersection(p,t)
    distance = 0.0
  else
    distance = typemax(0.0)
  end
end

function distance(p::Point{D},t::Triangle{D}) where D
  throw(ArgumentError("distance(::Point{$D},::Triangle{$D}) not imlemented"))
end

@inline distance(t::Triangle,p::Point) = distance(p,t)

function distance(s::Segment{3},t::Triangle{3})
  if have_intersection(s,t)
    0.0
  else
    min_dist = typemax(0.0)
    for i in 1:num_edges(t)
      e = get_edge(t,i)
      dist = distance(s,e)
      if dist < min_dist
        min_dist = dist
      end
    end
    min_dist
  end
end

distance(t::Triangle,s::Segment) = distance(s,t)

function contains_projection(p::Point{2},t::Triangle{2})
  s = measure_sign(t)
  s != 0 || throw(ErrorException("Triangle area is 0"))
  for i ∈ 1:num_edges(t)
    e = get_edge(t,i)
    n = normal(e) * s
    c = center(e)
    if ( p - c ) ⋅ n < 0
      return false
    end
  end
  true
end

have_intersection(p::Point{2},t::Triangle{2}) = contains_projection(p,t)

have_intersection_point(p::Point{2},t::Triangle{2}) = contains_projection(p,t)

function contains_projection(p::Point{3},t::Triangle{3})
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

function contains_projection(p::Point{D},t::Triangle{D}) where D
  throw(ArgumentError("contains_projection(::Point{$D},::Triangle{$D}) not implemented, only for 2D and 3D"))
end

function have_intersection(p::Point{D},t::Triangle{D}) where D
  throw(ArgumentError("have_intersection(::Point{$D},::Triangle{$D}) not implemented"))
end

function have_intersection_point(p::Point{D},t::Triangle{D}) where D
  throw(ArgumentError("have_intersection_point(::Point{D},::Triangle{D}) not defined, only in 2D"))
end

contains_projection(t::Triangle,p::Point) = contains_projection(p,t)

have_intersection(t::Triangle,p::Point) = have_intersection(p,t)

have_intersection_point(t::Triangle,p::Point) = have_intersection_point(p,t)

function have_intersection_point(s::Segment{3},t::Triangle{3})
  n = normal(t)
  c = center(t)
  s1_s2 = s[2] - s[1]
  s1_c = c - s[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  if α < 0 || α > 1 || isnan(α)
    false
  else
    x = s[1] + s1_s2 * α
    contains_projection(x,t)
  end
end

have_intersection(s::Segment{3},t::Triangle{3}) = have_intersection_point(s,t)

have_intersection_point(t::Triangle,s::Segment) = have_intersection_point(s,t)

have_intersection(t::Triangle,s::Segment) = have_intersection(s,t)

function have_intersection_point(s::Segment{D},t::Triangle{D}) where D
  throw(ArgumentError("have_intersection_point(::Segment{D},::Triangle{D}) only defined in 2D"))
end

function have_intersection(s::Segment{D},t::Triangle{D}) where D
  throw(ArgumentError("have_intersection(::Segment{$D},::Triangle{$D}) not implemented"))
end

function intersection(s::Segment{3},t::Triangle{3})
  @check have_intersection_point(s,t)
  n = normal(t)
  c = center(t)
  s1_s2 = s[2] - s[1]
  s1_c = c - s[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  s[1] + s1_s2 * α
end

function projection(p::Point{3},t::Triangle{3})
  @check contains_projection(p,t)
  c = center(t)
  n = normal(t)
  n = n / norm(n)
  p + ( ( c - p ) ⋅ n ) * n
end

function projection(p::Point{2},t::Triangle{2})
  @check contains_projection(p,t)
  p
end

function closest_point(t::Triangle,p::Point)
  projection(p,t)
end

function closest_point(s::Segment{3},t::Triangle{3})
  intersection(s,t)
end

function closest_point(t::Triangle{3},s::Segment{3})
  intersection(s,t)
end

function relative_orientation(s::Segment{2},t::Triangle{2})
  max_distance = 0.0
  _v = get_vertices(t)[1]
  for v in get_vertices(t)
    dist = distance(v,s)
    if dist ≥ max_distance
      max_distance = dist
      _v = v
    end
  end
  _t = Triangle(s,_v)
  - measure_sign(_t)
end

function writevtk(t::Triangle{D,T},file_base_name) where {D,T}
  vtk_type_id = 5

  points = zeros(T,D,num_vertices(t))
  for (i,v) in enumerate(get_vertices(t)), d in 1:D
      points[d,i] = v[d]
  end

  vtk_type = VTKCellType(vtk_type_id)
  vertices = [1:num_vertices(t);]
  cells = [ MeshCell(vtk_type,vertices) ]

  vtkfile = vtk_grid(file_base_name,points,cells)
  vtk_save(vtkfile)
end
