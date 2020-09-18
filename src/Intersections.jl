
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
  function CellFace{Df}(df::Integer,c::Cell{Dc,Dp}) where {Df,Dc,Dp}
    new{Df,Dp,typeof(c)}(df,c)
  end
  function CellFace{Df}(df::Integer,c::CellFace{Dc,Dp}) where {Df,Dc,Dp}
    new{Df,Dp,typeof(c)}(df,c)
  end
end

Segment(a...) = Segment(a)

Triangle(a...) = Triangle(a)

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
  v1 = a[2]-a[1]
  v2 = a[3]-a[1]
  n = orthogonal(v1,v2)
  n / norm(n)
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

distance(a::Point,b::Point) =  norm( a - b )

function distance(p::Point{D},s::Face{1,D}) where D
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

function have_intersection(f::Face{D,D},e::Face{1,D}) where D
  p0 = nothing
  for i in 1:num_facets(f)
    facet = get_facet(f,i)
    if have_intersection(facet,e)
      p = intersection_point(facet,e)
      if p0 === nothing
        p0 = p
      elseif distance(p,p0) > TOL
        return true
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


## Interface

function distance(K,X::Vector{<:Point},p::Polytope,point::Point)
  c = Cell(K,X,p)
  distance(c,point)
end

function distance(K,X::Vector{<:Point},p::Polytope,d::Integer,dface,point::Point)
  # @generated 
  c = Cell(K,X,p)
  if d == 0
    f0 = get_dface(c,dface,Val{0}())
    distance(point,f0)
  elseif d == 1
    f1 = get_dface(c,dface,Val{1}())
    distance(point,f1)
  elseif d == 2
    f2 = get_dface(c,dface,Val{2}())
    distance(point,f2)
  elseif d == 3
    f3 = get_dface(c,dface,Val{3}())
    distance(point,f3)
  else
    @notimplemented
  end
end

function projection(K,X::Vector{<:Point},p::Polytope,point::Point)
  c = Cell(K,X,p)
  projection(c,point)
end

function projection(K,X::Vector{<:Point},p::Polytope,d::Integer,dface,point::Point)
  # @generated 
  c = Cell(K,X,p)
  if d == 0
    f0 = get_dface(c,dface,Val{0}())
    projection(point,f0)
  elseif d == 1
    f1 = get_dface(c,dface,Val{1}())
    projection(point,f1)
  elseif d == 2
    f2 = get_dface(c,dface,Val{2}())
    projection(point,f2)
  elseif d == 3
    f3 = get_dface(c,dface,Val{3}())
    projection(point,f3)
  else
    @notimplemented
  end
end

function have_intersection(K,X::Vector{<:Point},p::Polytope,f::F) where F<:Union{Point,Face}
  c = Cell(K,X,p)
  have_intersection(c,f)
end

function have_intersection(K,X::Vector{<:Point},p::Polytope,d::Integer,dface,f::F) where F<:Union{Point,Face}
  # @generated 
  c = Cell(K,X,p)
  if d == 0
    f0 = get_dface(c,dface,Val{0}())
    have_intersection(f0,f)
  elseif d == 1
    f1 = get_dface(c,dface,Val{1}())
    have_intersection(f1,f)
  elseif d == 2
    f2 = get_dface(c,dface,Val{2}())
    have_intersection(f2,f)
  elseif d == 3
    f3 = get_dface(c,dface,Val{3}())
    have_intersection(f3,f)
  else
    @notimplemented
  end
end

function intersection_point(K,X::Vector{<:Point},p::Polytope,f::F) where F<:Union{Point,Face}
  c = Cell(K,X,p)
  intersection_point(c,f)
end

function intersection_point(K,X::Vector{<:Point},p::Polytope,d::Integer,dface,f::F) where F<:Union{Point,Face}
  # @generated 
  c = Cell(K,X,p)
  if d == 0
    f0 = get_dface(c,dface,Val{0}())
    intersection_point(f0,f)
  elseif d == 1
    f1 = get_dface(c,dface,Val{1}())
    intersection_point(f1,f)
  elseif d == 2
    f2 = get_dface(c,dface,Val{2}())
    intersection_point(f2,f)
  elseif d == 3
    f3 = get_dface(c,dface,Val{3}())
    intersection_point(f3,f)
  else
    @notimplemented
  end
end

function intersection_point(K,X,p,face,f) 
  d = get_facedims(p)[face]
  dface = face - get_dimrange(p,d)[1] + 1
  intersection_point(K,X,p,d,dface,f)
end

function get_bounding_box(K,X::Vector{<:Point},p::Polytope)
  c = Cell(K,X,p)
  get_bounding_box(c)
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
    throw(ArgumentError("orthogonal(::VectorValue{D}...) only well-defined for D-1 VectorValues{D}'s"))
  end
  orthogonal(a,)
end

