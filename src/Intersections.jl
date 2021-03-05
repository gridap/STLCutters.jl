
abstract type Face{Df,Dp} end

struct Segment{D,T}<:Face{1,D}
  vertices::Tuple{Point{D,T},Point{D,T}}
end

struct Triangle{D,T}<:Face{2,D}
  vertices::Tuple{Point{D,T},Point{D,T},Point{D,T}}
end

struct Tetrahedron{D,T}<:Face{3,D}
  vertices::Tuple{Point{D,T},Point{D,T},Point{D,T},Point{D,T}}
end

struct Cell{Df,Dp,Tp,P<:Polytope{Df}}<:Face{Df,Dp}
  nodes::Vector{Int}
  vertex_coordinates::Vector{Point{Dp,Tp}}
  polytope::P
end

struct CellFace{Df,Dp,C}<:Face{Df,Dp}
  dface::Int
  cell::C
  function CellFace{Df}(df::Integer,c::Face{Dc,Dp}) where {Df,Dc,Dp}
    new{Df,Dp,typeof(c)}(df,c)
  end
end

struct GridCell{Df,Dp,G}<:Face{Df,Dp}
  cell::Integer
  grid::G
  function GridCell(grid::Grid{Df,Dp},cell::Integer) where {Df,Dp}
    new{Df,Dp,typeof(grid)}(cell,grid)
  end
end

struct Plane{D,T}
  origin::Point{D,T}
  normal::VectorValue{D,T}
end
  
Segment(a::Point...) = Segment(a)

Triangle(a::Point...) = Triangle(a)

Tetrahedron(a::Point...) = Tetrahedron(a)

function Plane(a::Face{Df,Dp}) where {Df,Dp}
  o = center(a)
  n = normal(a)
  Plane(o,n)
end

origin(a::Plane) = a.origin

function Base.getindex(::Face)
  @abstractmethod
end

Base.getindex(a::Segment,i::Integer) = a.vertices[i]

Base.getindex(a::Triangle,i::Integer) = a.vertices[i]

Base.getindex(a::Tetrahedron,i::Integer) = a.vertices[i]

function Base.getindex(a::Cell,i::Integer)
  K = a.nodes
  X = a.vertex_coordinates
  X[K[i]]
end

function Base.getindex(a::CellFace,i::Integer)
  p = get_polytope(a.cell)
  d = face_dim(a)
  dface_nodes = get_face_vertices(p,d)[a.dface]
  a.cell[dface_nodes[i]]
end

function Base.getindex(a::GridCell,i::Integer)
  get_cell_coordinates(a.grid)[a.cell][i]
end

Base.lastindex(a::Face) = num_vertices(a)

function get_polytope(::Face)
  @abstractmethod
end

get_polytope(::Segment) = SEGMENT

function get_polytope(::Triangle)
  TRI.face_orientations[3] = -1
  TRI
end

get_polytope(::Tetrahedron) = TET

function get_polytope(a::Cell)
  HEX.face_orientations[3] = 1
  HEX.face_orientations[4] = -1
  TRI.face_orientations[3] = -1
  a.polytope
end

get_polytope(a::CellFace{1}) = SEGMENT

function get_polytope(a::CellFace{2})
  if is_simplex(a.cell)
    TRI.face_orientations[3] = -1
    TRI
  elseif is_n_cube(a.cell)
    QUAD
  else
    @notimplemented
  end
end

function get_polytope(a::CellFace{3})
  if is_simplex(a.cell)
    TET
  elseif is_n_cube(a.cell)
    HEX
  else
    @notimplemented
  end
end

function get_polytope(a::GridCell)
  HEX.face_orientations[3] = 1
  HEX.face_orientations[4] = -1
  TRI.face_orientations[3] = -1
  get_polytope(get_cell_reffe(a.grid)[a.cell])
end

num_dims(::Face{Df,Dp}) where {Df,Dp} = Dp

face_dim(::Face{D}) where D = D

function num_vertices(a::Face)
  p = get_polytope(a)
  num_vertices(p)
end

function num_edges(a::Face)
  p = get_polytope(a)
  num_edges(p)
end

function num_facets(a::Face)
  p = get_polytope(a)
  num_facets(p)
end

function num_faces(a::Face,d::Integer)
  p = get_polytope(a)
  num_faces(p,d)
end

function num_faces(a::Face)
  p = get_polytope(a)
  num_faces(p)
end

get_vertex(a::Face,vertex::Integer) = get_dface(a,vertex,Val{0}())

get_edge(a::Face,edge::Integer) = get_dface(a,edge,Val{1}())

get_facet(a::Face{D},facet::Integer) where D = get_dface(a,facet,Val{D-1}())

get_cell(a::Face) = a

get_cell(g::Grid,i::Integer) = GridCell(g,i)

function get_dface(a::Face,dface::Integer,::Val{d}) where d
  p = get_polytope(a)
  dface_to_nodes = get_face_vertices(p,d)[dface]
  if is_simplex(a)
    dface_vertices = ntuple( i -> a[dface_to_nodes[i]], Val{d+1}() )
    simplex_face( dface_vertices )
  elseif is_n_cube(a)
    dface_vertices = ntuple( i -> a[dface_to_nodes[i]], Val{2^d}() )
    n_cube_face( dface_vertices )
  else
    @unreachable
  end
end

function get_dface(a::Cell,dface::Integer,::Val{d}) where d
  CellFace{d}(dface,a)
end

function get_dface(a::CellFace,dface::Integer,::Val{d}) where d
  CellFace{d}(dface,a)
end

function get_dface(a::GridCell,dface::Integer,::Val{d}) where d
  CellFace{d}(dface,a)
end

function get_dface(a::Face,dface::Integer,::Val{0})
  a[dface]
end

function get_dface(a::Cell,dface::Integer,::Val{0})
  a[dface]
end

function get_dface(a::CellFace,dface::Integer,::Val{0})
  a[dface]
end

function get_dface(a::GridCell,dface::Integer,::Val{0})
  a[dface]
end

is_simplex(a::Face) = is_simplex(get_polytope(a))

is_n_cube(a::Face) = is_n_cube(get_polytope(a))

is_simplex(::Segment) = true

is_simplex(::Triangle) = true

is_simplex(::Tetrahedron) = true

is_n_cube(::Segment) = true

is_n_cube(::Triangle) = false

is_n_cube(::Tetrahedron) = false

function simplex_face(v::NTuple{N,<:Point}) where N
  if N == 1
    v[1]
  elseif N == 2
    Segment(v)
  elseif N == 3
    Triangle(v)
  elseif N == 4
    Tetrahedron(v)
  else
    @notimplemented
  end
end

function n_cube_face(v::NTuple{N,<:Point}) where N
  @notimplemented
end

function num_simplices(f::Face)
  if is_simplex(f)
    1
  else
    p = get_polytope(f)
    length(simplexify(p)[1])
  end
end

function simplex(f::Face{Df},i::Integer) where Df
  @boundscheck 0 < i ≤ num_simplices(f) || throw(BoundsError())
  if is_simplex(f)
    f
  else
    p = get_polytope(f)
    simplex_to_nodes = simplexify(p)[1][i]
    simplex_vertices = ntuple( i -> f[simplex_to_nodes[i]], Val{Df+1}() )
    simplex_face( simplex_vertices )
  end
end

function center(a::Face)
  c = zero( typeof(a[1]) )
  for i in 1:num_vertices(a)
    c += a[i] / num_vertices(a)
  end
  c
end

center(a::Point) = a

center(a::Plane) = orgin(a)

function normal(::Face)
  @notimplemented
end

normal(a::Plane) = a.normal

function normal(a::Face{1,2})
  v = a[2] - a[1]
  n = orthogonal(v)
  n / norm(n)
end

function normal(a::Face{2,3})
  for i in 1:num_simplices(a)
    s = simplex(a,i)
    v1 = s[2]-s[1]
    v2 = s[3]-s[1]
    n = orthogonal(v1,v2)
    n *= simplex_sign(a,i)
    _norm = norm(n)
    if _norm > 0
      return n / _norm
    end
  end
  @unreachable
end

function normal(a::Segment{2})
  v = a[2] - a[1]
  n = orthogonal(v)
  n / norm(n)
end

function normal(a::Triangle{3})
  v1 = a[2]-a[1]
  v2 = a[3]-a[1]
  n = orthogonal(v1,v2)
  n / norm(n)
end

function normal(::Face,::Integer)
  @abstractmethod
end

function normal(f::Face{D,D},facet::Integer) where D
  p = get_polytope(f)
  n = normal( get_facet(f,facet) )
  n * get_facet_orientations(p)[facet]
end

function normal(f::Face{1,D},facet::Integer) where D
  if facet == 1
    n = f[1]-f[2]
  elseif facet == 2
    n = f[2]-f[1]
  else
    @unreachable
  end
  n / norm(n)
end

function normal(f::Face{2,3},facet::Integer)
  p = get_polytope(f)
  n = normal(f)
  _f = get_facet(f,facet)
  v = _f[2]-_f[1]
  v /= norm(v)
  v *= get_facet_orientations(p)[facet]
  n × v
end

function measure(f::Face)
  if is_simplex(f)
    simplex_measure(f)
  else
    m = 0.0
    for i in 1:num_simplices(f)
      s = simplex(f,i)
      m += simplex_measure(s)
    end
    m
  end
end

function measure(f::Face{1})
  norm( f[2] -f[1] )
end

function measure(c::GridCell{D,D,<:CartesianGrid}) where D
  pmin,pmax = get_bounding_box(c)
  measure(pmin,pmax)
end

function measure(pmin::Point,pmax::Point)
  v = pmax - pmin
  prod(Tuple(v))
end

function simplex_measure(f::Face{Df,Dp}) where {Df,Dp}
  @notimplementedif Df ≠ Dp-1
  @notimplementedif !is_simplex(f)
  factor = 1/factorial(Df)
  v = ntuple( i -> f[i+1]-f[1], Val{Df}() )
  norm( orthogonal(v...) ) * factor
end

function simplex_measure(f::Face{2,2})
  @notimplementedif !is_simplex(f)
  factor = 1 / 2
  v1 = f[2] - f[1]
  v2 = f[3] - f[1]
  abs( v1 × v2 ) * factor
end

function simplex_measure(f::Face{3,3})
  @notimplementedif !is_simplex(f)
  factor = 1 / 6
  v1 = f[2] - f[1]
  v2 = f[3] - f[1]
  v3 = f[4] - f[1]
  abs( (v1×v2) ⋅ v3 ) * factor
end

simplex_measure(f::Face{1}) = measure(f)

Base.length(f::Face{1}) = measure(f)

function surface(f::Face{Df,Dp}) where {Df,Dp}
  @notimplementedif Df ≠ Dp-1
  measure(f)
end

volume(f::Face{D,D}) where D = measure(f)

function signed_volume(f::Face{D,D}) where D
  if is_simplex(f)
    simplex_signed_volume(f)
  else
    vol = 0.0
    for i in 1:num_simplices(f)
      s = simplex(f,i)
      vol += simplex_signed_volume(s) * simplex_sign(f,i)
    end
    vol
  end
end

function simplex_signed_volume(f::Face{D,D}) where D
  @notimplementedif !is_simplex(f)
  factor = 1/factorial(D)
  v = ntuple( i -> f[i+1]-f[1], Val{D}() )
  data = ntuple( i -> v[ ((i-1)÷D)+1 ][ ((i-1)%D)+1 ], Val{D^2}() )
  t = TensorValue(data...)
  vol = det(t) * factor
  vol
end

function simplex_sign(f::Face,i::Integer)
  simplex_sign(get_polytope(f),i)
end

function simplex_sign(p::Polytope,i::Integer)
  c = Cell(get_face_vertices(p,num_dims(p))[1],get_vertex_coordinates(p),p)
  sign( simplex_signed_volume( simplex(c,i) ) ) 
end

function signed_distance(p::Point,Π::Plane)
  o = origin( Π )
  n = normal( Π )
  (p-o) ⋅ n 
end

distance(a::Point,b::Point) =  norm( a - b )

distance(s::Face{1},p::Point) = distance(p,s)

function distance(p::Point{D},s::Face{1,D}) where D
  s1_s2 = s[2] - s[1]
  l = norm(s1_s2)
  v = s1_s2 / l
  s1_p = p - s[1]
  s1_projection = ( s1_p ⋅ v ) * v
  s2_projection = s1_projection - s1_s2
  p_projection = s1_projection - s1_p
  if norm( s1_projection ) ≤ l && norm( s2_projection) ≤ l
    norm(p_projection)
  else
    min( distance(s[1],p), distance(s[2],p) )
  end
end

function distance(s1::Face{1,3},s2::Face{1,3})
  v1 = s1[2] - s1[1]
  v2 = s2[2] - s2[1]
  v1 /= norm(v1)
  v2 /= norm(v2)
  n = v1 × v2  
  if norm(n) < 1e-14
    return min( 
      distance(s1[1],s2), 
      distance(s1[2],s2), 
      distance(s2[1],s1), 
      distance(s2[2],s1) )
  end
  n /= norm(n)
  n1 = n × v1
  n1 /= norm(n1)
  c1 = center(s1)
  s1_s2 = s2[2] - s2[1]
  s1_c = c1 - s2[1]
  α1 = ( n1 ⋅ s1_c ) / ( n1 ⋅ s1_s2 )
  n2 = n × v2
  n2 /= norm(n2)
  c2 = center(s2)
  s1_s2 = s1[2] - s1[1]
  s1_c = c2 - s1[1]
  α2 = ( n2 ⋅ s1_c ) / ( n2 ⋅ s1_s2 )
  if α1 < 0 || α2 < 0 || α1 > 1 || α2 > 1
    return min( 
      distance(s1[1],s2),
      distance(s1[2],s2),
      distance(s2[1],s1),
      distance(s2[2],s1) )
  end
  abs( ( c2 - c1 ) ⋅ n )
end

function distance(p::Point{Dp},f::Face{Df,Dp}) where {Df,Dp}
  if contains_projection(f,p)
    distance_to_infinite_face(f,p)
  else
    min_dist = Inf
    for i in 1:num_facets(f)
      facet = get_facet(f,i)
      dist = distance(p,facet)
      if dist < min_dist
        min_dist = dist
      end
    end
    min_dist
  end
end

distance(f::Face,p::Point) = distance(p,f)

function distance(f::Face{D,D},p::Point{D}) where D
  if have_intersection(f,p)
    0.0
  else
    distance_to_boundary(f,p)
  end
end

function distance(a::Face{D,D},b::Face{Df,D}) where {Df,D}
  @notimplementedif D == Df
  if have_intersection(a,b)
    0.0
  else
    distance_to_boundary(a,b)
  end
end

function distance(a::Face{Df1,Dp},b::Face{Df2,Dp}) where {Df1,Df2,Dp}
  @notimplementedif Df1+Df2 ≠ Dp
  if have_intersection(a,b)
    0.0
  else
    distance_to_boundary(a,b)
  end
end

function distance(a::Face,fa::Integer,b::Face,fb::Integer)
  dispatch_faces(distance,a,fa,b,fb)
end

function distance(a::Union{Point,Face},b::Face,d::Integer,df::Integer)
  dispatch_face(distance,a,b,d,df)
end

function distance_to_infinite_face(a::Point{D},b::Point{D}) where D
  distance(a,b)
end

function distance_to_infinite_face(s::Face{1,3},p::Point{3})
  p′ = projection(p,s)
  distance(p,p′)
end

function distance_to_infinite_face(f::Face{Df,Dp},p::Point{Dp}) where {Df,Dp}
  @notimplementedif Df ≠ Dp-1
  n = normal(f)
  c = center(f)
  abs( (p-c) ⋅ n )
end

function distance_to_infinite_face(
  f1::Face{Df1,Dp},
  f2::Face{Df2,Dp}) where {Df1,Df2,Dp}

  @notimplementedif Df1 < Df2 || Df1 == Dp
  max_dist = 0.0
  for i in 1:num_vertices(f2)
    dist = distance_to_infinite_face(f1,f2[i])
    if dist > max_dist
      max_dist = dist
    end
  end
  max_dist
end

function distance_to_boundary(a::Face,b::Union{Point,Face})
  min_dist = Inf
  for i in 1:num_facets(a)
    facet = get_facet(a,i)
    if isa(facet,Point) || measure(facet) ≠ 0 
      dist = distance(facet,b)
      if dist < min_dist
        min_dist = dist
      end
    end
  end
  min_dist
end

function projection(p::Point{D},q::Point{D}) where D
  q
end

function projection(p::Point{D},s::Face{1,D}) where D
  c = center(s)
  v = s[2] - s[1]
  v = v / norm(v)
  c + ( ( p - c ) ⋅ v ) * v
end

function projection(p::Point{Dp},f::Face{Df,Dp}) where {Df,Dp}
  @notimplementedif Df ≠ Dp-1
  n = normal(f)
  c = center(f)
  p + ( ( c - p ) ⋅ n ) * n
end

function have_intersection(f::Face{D,D},p::Point;atol=nothing) where D
  if is_cartesian(f)
    pmin,pmax = get_bounding_box(f)
    all( pmin.data .< p.data ) || return false
    all( pmax.data .> p.data ) || return false
    true
  else
    contains_projection(f,p)
  end
end

function have_intersection(f::Face{Df,Dp},e::Face{1,Dp};atol=nothing) where {Df,Dp}
  @notimplementedif Df + 1 ≠ Dp
  c = center(f)
  n = normal(f)
  s1_s2 = e[2] - e[1]
  s1_c = c - e[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  if α < 0 || α > 1 || isnan(α)
    false
  else
    x = e[1]*(1-α) + e[2]* α
    contains_projection(f,x)
  end
end

have_intersection(e::Face{1},f::Face{2};atol=nothing) = have_intersection(f,e)

function have_intersection(a::Face{D,D},b::Face{1,D};atol::Real) where D
  _have_intersection(a,b,atol=atol)
end

function have_intersection(a::Face{D,D},b::Face{Df,D};atol::Real) where {Df,D}
  _have_intersection(a,b;atol=atol)
end

function _have_intersection(a::Face{D,D},b::Face{Df,D};atol::Real) where {Df,D}
  height = intersection_height(a,b,atol=atol) 
  !isnothing(height) && height > atol
end

function intersection_height(f::Face{D,D},e::Face{1,D};atol::Real) where D
  points = zero(e[1]),zero(e[1])
  n = 0
  for i in 1:num_vertices(e)
    if distance(f,e[i]) < atol
      p = e[i]
      points,n = _max_area(p,points,n)
    end
  end
  for i in 1:num_facets(f)
    facet = get_facet(f,i)
    if surface(facet) > 0 && have_intersection(facet,e)
      p = intersection_point(facet,e)
      points,n = _max_area(p,points,n)
    end
  end
  n < length(points) ? nothing : min_height(simplex_face(points))
end

function intersection_height(a::Face{D,D},b::Face{2,D};atol::Real) where D
  points = zero(b[1]),zero(b[1]),zero(b[1])
  n = 0
  for i in 1:num_vertices(b)
    if distance(a,b[i]) < atol
      p = b[i]
      points,n = _max_area(p,points,n)
    end
  end
  for i in 1:num_faces(a,D-1), j in 1:num_edges(b)
    face = get_dface(a,i,Val{D-1}())
    edge = get_edge(b,j)
    if measure(face) ≠ 0 &&
      distance_to_infinite_face(face,edge) > atol &&
      distance(face,edge) < atol
      p = intersection_point(face,edge)
      p = projection(p,face)
      points,n = _max_area(p,points,n)
    end
  end
  for i in 1:num_faces(a,D-2)
    face = get_dface(a,i,Val{D-2}())
    if distance(b,face) < atol
      p = intersection_point(b,face)
      p = projection(p,face)
      points,n = _max_area(p,points,n)
    end
  end
  n < length(points) ? nothing : min_height(simplex_face(points))
end


function _max_area(p::Point,points::Tuple,n::Integer)
  if n < length(points)
    repeated = false
    for i in 1:n
      if points[i] == p
        repeated = true
        break
      end
    end
    if !repeated
      n += 1
      points = Base.setindex(points,p,n)
    end
  else
    max_m = measure( simplex_face(points) )
    imax = UNSET
    for i in 1:length(points)
      _points = Base.setindex(points,p,i)
      m = measure( simplex_face(_points) )
      if m > max_m
        imax = i
        max_m = m
      end
    end
    if imax ≠ UNSET
      points = Base.setindex(points,p,imax)
    end
  end
  points,n
end
      
function min_height(e::Face{1})
  length(e)
end

function min_height(f::Face{Df}) where Df
  @notimplementedif !is_simplex(f)
  factor = 1/factorial(Df)
  vol = measure(f)
  max_surf = 0.0
  for i in 1:num_facets(f)
    facet = get_facet(f,i)
    surf = measure(facet)
    if surf > max_surf
      max_surf = surf
    end
  end
  vol / max_surf / factor
end

function have_intersection(a::Face,fa::Integer,b::Face,fb::Integer;atol=nothing)
  dispatch_face(have_intersection,a,fa,b,fb;atol=atol)
end

function intersection_point(f::Face{Df,Dp},e::Face{1,Dp}) where {Df,Dp}
  @notimplementedif Df + 1 ≠ Dp
  c = center(f)
  n = normal(f)
  s1_s2 = e[2] - e[1]
  s1_c = c - e[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  e[1]*(1-α) + e[2]*α
end

intersection_point(e::Face{1},f::Face{2}) = intersection_point(f,e)

function intersection_point(s1::Face{1,3},s2::Face{1,3})
  v1 = s1[2] - s1[1]
  v2 = s2[2] - s2[1]
  v1 /= norm(v1)
  v2 /= norm(v2)
  n = v1 × v2
  n /= norm(n)
  n1 = n × v1
  n1 /= norm(n1)
  c1 = center(s1)
  s1_s2 = s2[2] - s2[1]
  s1_c = c1 - s2[1]
  α1 = ( n1 ⋅ s1_c ) / ( n1 ⋅ s1_s2 )
  n2 = n × v2
  n2 /= norm(n2)
  c2 = center(s2)
  s1_s2 = s1[2] - s1[1]
  s1_c = c2 - s1[1]
  α2 = ( n2 ⋅ s1_c ) / ( n2 ⋅ s1_s2 )
  if α1 ≤ 0 || α2 ≤ 0 || α1 ≥ 1 || α2 ≥ 1
    min_dist = Inf
    p = s1[1]
    for (a,b) in zip((s2,s1),(s1,s2))
      for i in 1:num_vertices(a)
        dist = distance(a[i],b)
        if dist < min_dist
          p = a[i]
        end
      end
    end
    p = projection(p,s1)
    @assert contains_projection(s1,p)
    return p
  end
  s1[1] + α2 * s1_s2
end

function intersection_point(f::Face{D,D},p::Point) where D
  p
end

function intersection_point(f::Face,p::Point)
  projection(p,f)
end

function intersection_point(a::Point,b::Point)
  a
end

function intersection_point(p::Point,f::Face)
  p
end

function contains_projection(f::Face,p::Point)
  for i in 1:num_facets(f)
    facet = get_facet(f,i)
    if isa(facet,Point) || measure(facet) ≠ 0
      c = center(facet)
      n = normal(f,i)
      if (p-c) ⋅ n > 0
        return false
      end
    end
  end
  true
end

function contains_projection(f1::Face,f2::Face)
  for i in 1:num_vertices(f2)
    if !contains_projection(f1,f2[i])
      return false
    end
  end
  true
end

function is_on_boundary(f1::Face,f2::Union{Point,Face};atol::Real)
  for i in 1:num_facets(f1)
    facet = get_facet(f1,i)
    if surface(facet) > 0
      if distance_to_infinite_face(facet,f2) < atol
        if are_overlapped(facet,f2,atol=atol)
          return true
        end
      end
    end
  end
  false
end

is_cartesian(::Face) = false

function is_cartesian(f::Face{D,D}) where D
  if is_n_cube(f)
    pmin,pmax = f[1],f[end]
    for i in 1:num_vertices(f)
      for d in 1:D
        pd = ((i-1) & (1<<(d-1))) == 0 ? pmin[d] : pmax[d]
        if f[i][d] ≠ pd
          return false
        end
      end
    end
    return true
  end
  false
end

get_bounding_box(p::Point) = p,p

function get_bounding_box(f::Face)
  if is_cartesian(f)
    f[1],f[end]
  else
    pmin = pmax = f[1]
    for i in 1:num_vertices(f)
      pmin = Point(min.(Tuple(pmin),Tuple(f[i])))
      pmax = Point(max.(Tuple(pmax),Tuple(f[i])))
    end
    pmin,pmax
  end
end

function overlapping_height(f1::Face{1,2},f2::Face{1,2})
  points = zero(f2[1]),zero(f2[1])
  n = 0
  for (a,b) in zip((f1,f2),(f2,f1)), i in 1:num_vertices(a)
    if contains_projection(b,a[i])
      p = a[i]
      points,n = _max_area(p,points,n)
    end
  end
  n < length(points) ? nothing : min_height(simplex_face(points))
end

function overlapping_height(f1::Face{2,3},f2::Face{2,3})
  points = tfill(zero(f2[1]),Val{3}())
  n = 0
  for (a,b) in zip((f1,f2),(f2,f1))
    for i in 1:num_vertices(a)
      if contains_projection(b,a[i])
        p = a[i]
        points,n = _max_area(p,points,n)
      end
    end
  end
  for i in 1:num_facets(f1), j in 1:num_facets(f2)
    facet1 = get_facet(f1,i)
    facet2 = get_facet(f2,j)
    _o = center(facet2)
    _n = normal(f1) × (facet2[2]-facet2[1])
    _n /= norm(_n)
    α = ( (_o-facet1[1]) ⋅ _n ) / ( (facet1[2]-facet1[1]) ⋅ _n )
    if 0 ≤ α ≤ 1
      p = facet1[1] + α*(facet1[2]-facet1[1])
      if contains_projection(facet2,p)
        points,n = _max_area(p,points,n)
      end
    end
  end
  n < length(points) ? nothing : min_height(simplex_face(points))
end

function overlapping_height(f1::Face{2,3},f2::Face{1,3})
  points = tfill(zero(f2[1]),Val{2}())
  n = 0
  for i in 1:num_vertices(f2)
    if contains_projection(f1,f2[i])
      p = f2[i]
      points,n = _max_area(p,points,n)
    end
  end
  for i in 1:num_facets(f1)
    facet = get_facet(f1,i)
    _o = center(f2)
    _n = normal(f1) × (f2[2]-f2[1])
    _n /= norm(_n)
    α = ( (_o-facet[1]) ⋅ _n ) / ( (facet[2]-facet[1]) ⋅ _n )
    if 0 ≤ α ≤ 1
      p = facet[1] + α*(facet[2]-facet[1])
      if contains_projection(f2,p)
        points,n = _max_area(p,points,n)
      end
    end
  end
  n < length(points) ? nothing : min_height(simplex_face(points))
end

function are_overlapped(a::Face,b::Face;atol::Real)
  height = overlapping_height(a,b)
  !isnothing(height) && height > atol
end

function are_overlapped(f::Face,p::Point;atol)
  contains_projection(f,p)
end

function fast_intersection(p::Point,pmin::Point,pmax::Point)
  all( Tuple(p) .> Tuple(pmin) ) || return false
  all( Tuple(p) .< Tuple(pmax) ) || return false
  true
end

function fast_intersection(e::Face{1,D},pmin::Point{D},pmax::Point{D}) where D
  t_min,t_max = 0.0,1.0
  for d in 1:D
    p_d = e[1][d]
    v_d = e[2][d] - e[1][d]
    if v_d < 0
      v_d = -v_d
      p_d = - p_d + pmin[d] + pmax[d]
    end
    if v_d != 0
      t = ( pmin[d] - p_d ) / v_d
      if t > t_min
        t_min = t
      end
      t = ( pmax[d] - p_d ) / v_d
      if t < t_max
        t_max = t
      end
    else
      if p_d < pmin[d]
        return false
      end
      if p_d > pmax[d]
        return false
      end
    end
  end
  t_min < t_max
end

function fast_intersection(Π::Plane{D},pmin::Point{D},pmax::Point{D}) where D
  n = normal(Π)
  p = origin(Π)
  f = d -> pmin[d] + pmax[d] - p[d]
  _p = Point(ntuple( d -> n[d] < 0 ? f(d) : p[d], Val{D}() ))
  _n = VectorValue(abs.(Tuple(n)))
  _Π = Plane(_p,_n)
  signed_distance(pmin,_Π) < 0 || return false
  signed_distance(pmax,_Π) > 0 || return false
  true
end


function _compute_distances(Π::Plane{D},pmin::Point{D},pmax::Point{D},p::Polytope{D}) where D
  @assert is_n_cube(p)
  ntuple(i->signed_distance(_get_vertex(p,pmin,pmax,i),Π),Val{2^D}())
end

function _get_vertex(p::Polytope{D},pmin::Point{D},pmax::Point{D},i::Integer) where D
  Point(ntuple(d-> (i-1) & (1<<(d-1)) == 0 ? pmin[d] : pmax[d], Val{D}()))
end

function fast_intersection(f::Face{2,3},pmin::Point{3},pmax::Point{3})
  D = 3
  n = normal(f)
  c = center(f)
  plane = Plane(c,n)
  fast_intersection(plane,pmin,pmax) || return false
  for i in 1:num_vertices(f)
    v = f[i]
    fast_intersection(v,pmin,pmax) && return true
  end
  for i in 1:num_edges(f)
    e = get_edge(f,i)
    fast_intersection(e,pmin,pmax) && return true
  end

  fd = (p,d) -> pmin[d] + pmax[d] - p[d]
  fp = p -> Point(ntuple( d -> n[d] < 0 ? fd(p,d) : p[d], Val{D}() ))
  vertices = ntuple(i-> fp(f[i]),3)
  t = Triangle(vertices)
  contains_projection(t,pmin) && return true
  contains_projection(t,pmax) && return true
  false
end

function fast_intersection(f::Face{2,3},pmin::Point{3},pmax::Point{3},p::Polytope{3})
  D = 3
  for i in 1:num_vertices(f)
    v = f[i]
    fast_intersection(v,pmin,pmax) && return true
  end  
  for i in 1:num_edges(f)
    e = get_edge(f,i)
    fast_intersection(e,pmin,pmax) && return true
  end
  n = normal(f)
  c = center(f)
  plane = Plane(c,n)
  v_to_dists = _compute_distances(plane,pmin,pmax,p)
  all( d -> d > 0, v_to_dists ) && return false
  all( d -> d < 0, v_to_dists ) && return false
  for e in 1:num_edges(p)
    if _is_edge_cut(p,v_to_dists,e)
      i = _intersection_point(p,pmin,pmax,v_to_dists,e)
      if contains_projection(f,i)
        return true
      end
    end
  end
  false
end

function _is_edge_cut(p::Polytope,v_to_dists::Tuple,e::Integer)
  e_to_v = get_face_vertices(p,1)[e]
  sign(v_to_dists[e_to_v[1]]) ≠ sign(v_to_dists[e_to_v[2]])
end

function _intersection_point(p::Polytope,pmin::Point,pmax::Point,v_to_dists::Tuple,e::Integer)
  e_to_v = get_face_vertices(p,1)[e]
  d1 = v_to_dists[e_to_v[1]]
  d2 = v_to_dists[e_to_v[2]]
  v1 = _get_vertex(p,pmin,pmax,e_to_v[1])
  v2 = _get_vertex(p,pmin,pmax,e_to_v[2])
  (v1*d2-v2*d1)/(d2-d1)
end

## Generated funcs

function dispatch_faces(fun,a::Face,fa::Integer,b::Face,fb::Integer)
  pa = get_polytope(a)
  da = get_facedims(pa)[fa]
  dfa = fa - get_offset(pa,da)
  pb = get_polytope(b)
  db = get_facedims(pb)[fb]
  dfb = fb - get_offset(pb,db)
  dispatch_face(fun,a,da,dfa,b,db,dfb)
end

@generated(
function dispatch_face(
  fun,a::Face{D},d::Integer,dface::Integer,b...;kargs...) where D

  str = ""
  for d in 0:D
    str *= "if d == $d \n"
    str *= "  f$d =  get_dface(a,dface,Val{$d}()) \n"
    str *= "  fun(f$d,b...;kargs...) \n"
    str *= "else"
  end
  str *= "\n  @notimplemented \nend"
  Meta.parse(str)
end
)

@generated(
function dispatch_face(
  fun,a,b::Face{D},d::Integer,dface::Integer;kargs...) where D

  str = ""
  for d in 0:D
    str *= "if d == $d \n"
    str *= "  f$d =  get_dface(b,dface,Val{$d}()) \n"
    str *= "  fun(a,f$d;kargs...) \n"
    str *= "else"
  end
  str *= "\n  @notimplemented \nend"
  Meta.parse(str)
end
)

## Interface

function distance(K,X::Vector{<:Point},p::Polytope,point::Point)
  c = Cell(K,X,p)
  distance(c,point)
end

function distance(
  K,
  X::Vector{<:Point},
  p::Polytope,
  face::Integer,
  f::Union{Point,Face})
  
  d = get_facedims(p)[face]
  dface = face - get_offset(p,d)
  distance(K,X,p,d,dface,f)
end

function distance(
  K,
  X::Vector{<:Point},
  p::Polytope,
  d::Integer,
  dface::Integer,
  point::Point)

  c = Cell(K,X,p)
  distance(c,d,dface,point)
end

function distance(c::Cell{D},d::Integer,dface::Integer,point::Point) where D
  dispatch_face(distance,c,d,dface,point)
end

function projection(K,X::Vector{<:Point},p::Polytope,point::Point)
  c = Cell(K,X,p)
  projection(point,c)
end

function projection(
  K,
  X::Vector{<:Point},
  p::Polytope,
  d::Integer,
  dface::Integer,
  point::Point)

  c = Cell(K,X,p)
  projection(c,d,dface,point)
end

function projection(c::Cell{D},d::Integer,dface::Integer,point::Point) where D
  dispatch_face(projection,point,c,d,dface)
end

function have_intersection(
  K,
  X::Vector{<:Point},
  p::Polytope,
  f::Union{Point,Face}
  ;atol=nothing) 

  c = Cell(K,X,p)
  have_intersection(c,f,atol=atol)
end

function have_intersection(
  K,
  X::Vector{<:Point},
  p::Polytope,
  face::Integer,
  f::Union{Point,Face}
  ;atol=nothing)
  
  d = get_facedims(p)[face]
  dface = face - get_offset(p,d)
  have_intersection(K,X,p,d,dface,f,atol=atol)
end

function have_intersection(
  K,
  X::Vector{<:Point},
  p::Polytope,
  d::Integer,
  dface::Integer,
  f::Union{Point,Face}
  ;atol=nothing)

  c = Cell(K,X,p)
  have_intersection(c,d,dface,f,atol=atol)
end

function have_intersection(
  c::Cell{D},
  d::Integer,
  dface::Integer,
  f::Union{Point,Face}
  ;atol=nothing) where D

  dispatch_face(have_intersection,c,d,dface,f;atol=atol)
end

function intersection_point(
  K,
  X::Vector{<:Point},
  p::Polytope,
  f::F) where F<:Union{Point,Face}

  c = Cell(K,X,p)
  intersection_point(c,f)
end

function intersection_point(
  K,
  X::Vector{<:Point},
  p::Polytope,
  d::Integer,
  dface::Integer,
  f::F) where F<:Union{Point,Face}

  c = Cell(K,X,p)
  intersection_point(c,d,dface,f)
end

function intersection_point(
  c::Cell{D},
  d::Integer,
  dface::Integer,
  f::F) where {D,F<:Union{Point,Face}}

  dispatch_face(intersection_point,c,d,dface,f)
end

function intersection_point(a::Face,fa::Integer,b::Face,fb::Integer)
  dispatch_faces(intersection_point,a,fa,b,fb)
end

function intersection_point(a::Union{Point,Face},b::Face,d::Integer,df::Integer)
  dispatch_face(intersection_point,a,b,d,df)
end

function intersection_point(K,X::Vector{<:Point},p::Polytope,face::Integer,f) 
  d = get_facedims(p)[face]
  dface = face - get_dimrange(p,d)[1] + 1
  intersection_point(K,X,p,d,dface,f)
end

function is_on_boundary(
  K,
  X::Vector{<:Point},
  p::Polytope,
  f::Union{Point,Face}
  ;atol::Real)

  c = Cell(K,X,p)
  is_on_boundary(c,f,atol=atol)
end

function get_bounding_box(K,X::Vector{<:Point},p::Polytope)
  c = Cell(K,X,p)
  get_bounding_box(c)
end

function get_facet(K,X::Vector{<:Point},p::Polytope,i::Integer)
  c = Cell(K,X,p)
  get_facet(c,i)
end

function get_bounding_box(c::Cell{D,D,<:CartesianGrid}) where D
  c[1],c[end]
end

## Orthogonal

function orthogonal(a::VectorValue{2})
  VectorValue( -a[2], a[1] )
end

@generated function orthogonal(a::NTuple{N,VectorValue{D}}) where {N,D}
  entries = ""
  for i in 1:D
    data = ""
    for j in 1:D-1
      for k in 1:D
        if k != i
          data *= "a[$j][$k],"
        end
      end
    end
    if iseven(D+i)
      entries *= " + "
    else
      entries *= " - "
    end
    entries *= "det(TensorValue($data)),\n"
  end
  str = "VectorValue(\n$entries)"
  Meta.parse(str)
end

function orthogonal(a::VectorValue{D}...) where D
  if length(a) != D-1
    msg = 
    "orthogonal(::VectorValue{D}...) only defined for D-1 VectorValues{D}'s"
    throw(ArgumentError(msg))
  end
  orthogonal(a,)
end
