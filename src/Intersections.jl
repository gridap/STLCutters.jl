
abstract type Face{Df,Dp} end

struct Segment{D,T}<:Face{1,D}
  vertices::Tuple{Point{D,T},Point{D,T}}
end

struct Triangle{D,T}<:Face{2,D}
  vertices::Tuple{Point{D,T},Point{D,T},Point{D,T}}
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

Segment(a::Point...) = Segment(a)

Triangle(a::Point...) = Triangle(a)

function Base.getindex(::Face)
  @abstractmethod
end

Base.getindex(a::Segment,i::Integer) = a.vertices[i]

Base.getindex(a::Triangle,i::Integer) = a.vertices[i]

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

Base.lastindex(a::Face) = num_vertices(a)

function get_polytope(::Face)
  @abstractmethod
end

get_polytope(::Segment) = SEGMENT

function get_polytope(::Triangle)
  TRI.face_orientations[3] = -1
  TRI
end

get_polytope(a::Cell) = a.polytope

get_polytope(a::CellFace{1}) = SEGMENT

function get_polytope(a::CellFace{2})
  if is_simplex(a.cell)
    TRI
  elseif is_n_cube(a.cell)
    QUAD
  else
    @notimplemented
  end
end

function get_polytope(a::CellFace{3})
  if is_simplex(a.cell)
    TRI
  elseif is_n_cube(a.cell)
    QUAD
  else
    @notimplemented
  end
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

function get_dface(a::Face,dface::Integer,::Val{0})
  a[dface]
end

function get_dface(a::Cell,dface::Integer,::Val{0})
  a[dface]
end

function get_dface(a::CellFace,dface::Integer,::Val{0})
  a[dface]
end

is_simplex(a::Face) = is_simplex(get_polytope(a))

is_n_cube(a::Face) = is_n_cube(get_polytope(a))

is_simplex(::Segment) = true

is_simplex(::Triangle) = true

is_n_cube(::Segment) = true

is_n_cube(::Triangle) = true

function simplex_face(v::NTuple{N,<:Point}) where N
  if N == 1
    v[1]
  elseif N == 2
    Segment(v)
  elseif N == 3
    Triangle(v)
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
    length(simplexify(p))
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

function normal(::Face)
  @notimplemented
end

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
  normal( get_facet(f,facet) )
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
  abs( v1 ⋅ v2 ) * factor
end

simplex_measure(f::Face{1}) = measure(f)

Base.length(f::Face{1}) = measure(f)

function surface(f::Face{Df,Dp}) where {Df,Dp}
  @notimplementedif Df ≠ Dp-1
  measure(f)
end

volume(f::Face{D,D}) where D = measure(f)

distance(a::Point,b::Point) =  norm( a - b )

distance(s::Face{1},p::Point) = distance(p,s)

function distance(p::Point{D},s::Face{1,D}) where D
  s1_s2 = s[2] - s[1]
  l = norm(s1_s2)
  @assert l > 0
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
  @notimplementedif Df ≠ Dp-1
  n = normal(f)
  c = center(f)
  abs( (p-c) ⋅ n )
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

function have_intersection(f::Face{Df,Dp},e::Face{1,Dp}) where {Df,Dp}
  @notimplementedif Df + 1 ≠ Dp
  c = center(f)
  n = normal(f)
  s1_s2 = e[2] - e[1]
  s1_c = c - e[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  if α < -TOL || α > 1+TOL || isnan(α)
    false
  else
    x = e[1] + s1_s2 * α
    contains_projection(f,x)
  end
end

function have_intersection(f1::Face{1,3},f2::Face{1,3})
  !are_colinear(f1,f2) && distance(f1,f2) < TOL
end

function have_intersection(f::Face{D,D},e::Face{1,D}) where D
  p0 = nothing
  for i in 1:num_facets(f)
    facet = get_facet(f,i)
    if surface(facet) > 0
      if have_intersection(facet,e)
        p = intersection_point(facet,e)
        if p0 === nothing
          p0 = p
        elseif distance(p,p0) > TOL
          return true
        end
      end
    end
  end
  false
end

function have_intersection(a::Face{D,D},b::Face{2,D}) where D
  p0,p1 = nothing,nothing
  for i in 1:num_faces(a,D-2)
    face = get_dface(a,i,Val{D-2}())
    if have_intersection(b,face)
      p = intersection_point(b,face)
      if p0 === nothing
        p0 = p
      elseif p1 === nothing
        if distance(p,p0) > TOL
          p1 = p
        end
      elseif distance(p,Segment(p0,p1)) > TOL
        return true
      end
    end
  end
  false
end

function have_intersection(f::Face{D,D},p::Point) where D
  if is_cartesian(f)
    pmin,pmax = get_bounding_box(f)
    all( pmin.data .< p.data ) || return false
    all( pmax.data .> p.data ) || return false
    true
  else
    @notimplemented
  end
end

function have_intersection(f::Face,p::Point)
  @notimplemented
end



function intersection_point(f::Face{Df,Dp},e::Face{1,Dp}) where {Df,Dp}
  @notimplementedif Df + 1 ≠ Dp
  c = center(f)
  n = normal(f)
  s1_s2 = e[2] - e[1]
  s1_c = c - e[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  e[1] + s1_s2 * α
end

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

function contains_projection(f::Face,p::Point)
  for i in 1:num_facets(f)
    facet = get_facet(f,i)
    c = center(facet)
    n = normal(f,i)
    if (p-c) ⋅ n > TOL
      return false
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

function is_on_cell_facet(c::Face{D,D},f::Face{Df,D}) where {Df,D}
  for i in 1:num_facets(c)
    facet = get_facet(c,i)
    if surface(facet) > 0
      if are_coplanar(facet,f)
        if are_overlapped(facet,f)
          return true
        end
      end
    end
  end
  false
end

function is_on_cell_facet(c::Face{D,D},f::Face{Df,D},face::Integer) where {Df,D}
  for i in 1:num_facets(c)
    facet = get_facet(c,i)
    if surface(facet) > 0
       are_coplanar(facet,f) &&
       are_overlapped(facet,f) &&
       are_overlapped(facet,f,face)

       return true
    end
  end
  false
end

function get_cell_facet_on_facet(
  c::Face{D,D},
  f::Face{Df,D},
  face::Integer) where {Df,D}

  for i in 1:num_facets(c)
    facet = get_facet(c,i)
    if surface(facet) > 0
       are_coplanar(facet,f) &&
       are_overlapped(facet,f) &&
       are_overlapped(facet,f,face)

       return i
    end
  end
  @unreachable
end

is_cartesian(::Face) = false

function is_cartesian(f::Face{D,D}) where D
  if is_n_cube(f)
    p = get_polytope(f)
    for i in 1:num_facets(f)
      n = get_facet_normals(p)[i]
      facet = get_facet(f,i)
      pd = n ⋅ facet[1]
      for j in 1:num_vertices(facet)
        if n ⋅ facet[j] ≠ pd
          return false
        end
      end
    end
    return true
  end
  false
end

function get_bounding_box(f::Face)
  @assert is_cartesian(f)
  f[1],f[end]
end

function are_coplanar(::Face{Df,D},::Face{D,D}) where {Df,D}
  @notimplemented
end

function are_coplanar(::Face{D,D},::Face{Df,D}) where {Df,D}
  @notimplemented
end

function are_coplanar(::Face{D,D},::Face{D,D}) where D
  @notimplemented
end

function are_coplanar(f1::Face{Df1,Dp},f2::Face{Df2,Dp}) where {Df1,Df2,Dp}
  @notimplementedif Df1 ≠ Dp - 1 || Df1 < Df2
  n = normal(f1)
  c = center(f1)
  for i in 1:num_vertices(f2)
    p = f2[i]
    if abs( (p-c)⋅n ) > TOL
      return false
    end
  end
  true
end

function are_colinear(s1::Face{1,3},s2::Face{1,3})
  θmin = 1e-13
  v1 = s1[2] - s1[1]
  v2 = s2[2] - s2[1]
  v1 /= norm(v1)
  v2 /= norm(v2)
  n = v1 × v2
  norm(n) < θmin
end

function are_overlapped(f1::Face{2,3},f2::Face{2,3})
  @notimplementedif !are_coplanar(f1,f2)
  p0,p1 = nothing,nothing
  for (a,b) in zip((f1,f2),(f2,f1))
    for i in 1:num_vertices(a)
      if contains_projection(b,a[i])
        p = a[i]
        if p0 === nothing
          p0 = p
        elseif p1 === nothing
          if distance(p,p0) > TOL
            p1 =  p
          end
        elseif distance(p,Segment(p0,p1)) > TOL
          return true
        end
      end
    end
  end
  for i in 1:num_facets(f1), j in 1:num_facets(f2)
    facet1 = get_facet(f1,i)
    facet2 = get_facet(f2,j)
    if have_intersection(facet1,facet2)
      p = intersection_point(facet1,facet2)
      if p0 === nothing
        p0 = p
      elseif p1 === nothing
        if distance(p,p0) > TOL
          p1 =  p
        end
      elseif distance(p,Segment(p0,p1)) > TOL
        return true
      end
    end
  end
  false
end

function are_overlapped(f1::Face{1,2},f2::Face{1,2})
  @notimplementedif !are_coplanar(f1,f2)
  p0 = nothing
  for (a,b) in zip((f1,f2),(f2,f1))
    for i in 1:num_vertices(a)
      if contains_projection(b,a[i])
        p = a[i]
        if p0 === nothing
          p0 = p
        elseif distance(p,p0) > TOL
          return true
        end
      end
    end
  end
  false
end

function are_overlapped(f::Face,p::Point)
  contains_projection(f,p)
end

function are_overlapped(f::Face{2,3},s::Face{1,3})
  p0 = nothing
  for i in 1:num_vertices(s)
    if contains_projection(f,s[i])
      if p0 == nothing
        p0 = s[i]
      elseif distance(p0,s[i]) > TOL
        return true
      end
    end
  end
  for i in 1:num_facets(f)
    facet = get_facet(f,i)
    if have_intersection(facet,s)
      p = intersection_point(facet,s)
      if p0 === nothing
        p0 = p
      elseif distance(p0,p) > TOL
        return true
      end
    end
  end
  false
end

function are_overlapped(c::Face,f::Face,face::Integer)
  p = get_polytope(f)
  d = get_facedims(p)[face]
  dface = face - get_dimrange(p,d)[1] + 1
  are_overlapped(c,f,d,dface)
end

function are_overlapped(c::Face,f::Face{D},d::Integer,dface::Integer) where D
  dispatch_face(are_overlapped,c,f,d,dface)
end

## Generated funcs

@generated(
function dispatch_face(fun,a::Face{D},d::Integer,dface::Integer,b...) where D
  str = ""
  for d in 0:D
    str *= "if d == $d \n"
    str *= "  f$d =  get_dface(a,dface,Val{$d}()) \n"
    str *= "  fun(f$d,b...) \n"
    str *= "else"
  end
  str *= "\n  @notimplemented \nend"
  Meta.parse(str)
end
)

@generated(
function dispatch_face(fun,a,b::Face{D},d::Integer,dface::Integer) where D
  str = ""
  for d in 0:D
    str *= "if d == $d \n"
    str *= "  f$d =  get_dface(b,dface,Val{$d}()) \n"
    str *= "  fun(a,f$d) \n"
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
  f::F) where F<:Union{Point,Face}

  c = Cell(K,X,p)
  have_intersection(c,f)
end

function have_intersection(
  K,
  X::Vector{<:Point},
  p::Polytope,
  d::Integer,
  dface::Integer,
  f::F) where F<:Union{Point,Face}

  c = Cell(K,X,p)
  have_intersection(c,d,dface,f)
end

function have_intersection(
  c::Cell{D},
  d::Integer,
  dface::Integer,
  f::F) where {D,F<:Union{Point,Face}}

  dispatch_face(have_intersection,c,d,dface,f)
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

function intersection_point(K,X,p,face,f) 
  d = get_facedims(p)[face]
  dface = face - get_dimrange(p,d)[1] + 1
  intersection_point(K,X,p,d,dface,f)
end

function is_on_cell_facet(
  K,
  X::Vector{<:Point},
  p::Polytope,
  f::F) where F<:Union{Point,Face}

  c = Cell(K,X,p)
  is_on_cell_facet(c,f)
end

function is_on_cell_facet(
  K,
  X::Vector{<:Point},
  p::Polytope,
  f::F,
  df::Integer) where F<:Union{Point,Face}

  c = Cell(K,X,p)
  is_on_cell_facet(c,f,df)
end

function get_cell_facet_on_facet(
  K,
  X::Vector{<:Point},
  p::Polytope,
  f::F) where F<:Union{Point,Face}

  c = Cell(K,X,p)
  get_cell_facet_on_facet(c,f)
end

function get_cell_facet_on_facet(
  K,
  X::Vector{<:Point},
  p::Polytope,
  f::F,
  df::Integer) where F<:Union{Point,Face}

  c = Cell(K,X,p)
  get_cell_facet_on_facet(c,f,df)
end

function get_bounding_box(K,X::Vector{<:Point},p::Polytope)
  c = Cell(K,X,p)
  get_bounding_box(c)
end

function get_facet(K,X::Vector{<:Point},p::Polytope,i::Integer)
  c = Cell(K,X,p)
  get_facet(c,i)
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
