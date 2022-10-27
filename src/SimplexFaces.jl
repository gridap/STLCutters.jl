
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

abstract type AbstractPlane{D,T} end

struct Plane{D,T} <: AbstractPlane{D,T}
  origin::Point{D,T}
  normal::VectorValue{D,T}
end

struct CartesianPlane{D,T} <: AbstractPlane{D,T}
  d::Int8
  value::T
  positive::Bool
end

# Constructors

Segment(a::Point...) = Segment(a)

Triangle(a::Point...) = Triangle(a)

Tetrahedron(a::Point...) = Tetrahedron(a)

function Plane(a::Face{Df,Dp}) where {Df,Dp}
  o = center(a)
  n = normal(a)
  Plane(o,n)
end

function CartesianPlane(p::Point{D},d::Integer,orientation::Integer) where D
  @notimplementedif abs( orientation ) ≠ 1
  T = eltype(p)
  CartesianPlane{D,T}(Int8(d),p[d],orientation>0)
end

# Getters

origin(a::Plane) = a.origin

function normal(a::CartesianPlane{D,T}) where {D,T}
  n = zero( VectorValue{D,T} )
  data = Tuple(n)
  v = a.positive ? 1 : -1
  data = Base.setindex(data,v,a.d)
  VectorValue( data )
end

function Base.zero(::Type{<:Plane{D,T}}) where {D,T}
  Plane(zero(Point{D,T}),zero(Point{D,T}))
end

function origin(p::CartesianPlane{D,T}) where {D,T}
  o = zero( Point{D,T} )
  Point( Base.setindex(o.data,p.value,p.d) )
end

function Base.getindex(::Face)
  @abstractmethod
end

Base.getindex(a::Segment,i::Integer) = a.vertices[i]

Base.getindex(a::Triangle,i::Integer) = a.vertices[i]

Base.getindex(a::Tetrahedron,i::Integer) = a.vertices[i]

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

num_dims(::Face{Df}) where Df = Df

num_point_dims(::Face{Df,Dp}) where {Df,Dp} = Dp

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

@inline function get_dface(a::Face,dface::Integer,::Val{d}) where d
  @notimplementedif !is_simplex(a)
  p = get_polytope(a)
  dface_to_nodes = get_face_vertices(p,d)[dface]
  dface_vertices = ntuple( i -> a[dface_to_nodes[i]], Val{d+1}() )
  simplex_face( dface_vertices )
end

function get_dface(a::Face,dface::Integer,::Val{0})
  a[dface]
end

is_simplex(a::Face) = is_simplex(get_polytope(a))

is_n_cube(a::Face) = is_n_cube(get_polytope(a))

simplex_face(v::Point...) = simplex_face(v)

function simplex_face(coords::AbstractVector{<:Point{D}}) where D
  tcoords = ntuple( i -> coords[i], Val{D+1}() )
  simplex_face(tcoords)
end

@inline function simplex_face(v::NTuple{N,<:Point}) where N
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

function is_cartesian(f::Face)
  @notimplementedif !is_simplex(f)
  false
end

function center(a::Face)
  c = zero( typeof(a[1]) )
  for i in 1:num_vertices(a)
    c += a[i] / num_vertices(a)
  end
  c
end

center(a::Point) = a

center(a::Plane) = origin(a)

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
  @notimplementedif !is_simplex(a)
  v1 = a[2]-a[1]
  v2 = a[3]-a[1]
  n = orthogonal(v1,v2)
  _norm = norm(n)
  @assert _norm > 0
  n / _norm
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
  @notimplementedif !is_simplex(f)
  simplex_measure(f)
end

function measure(f::Face{1})
  norm( f[2] -f[1] )
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

function min_height(e::Face{1})
  length(e)
end

function min_height(f::Face{Df}) where Df
  @notimplementedif !is_simplex(f)
  min_dist = Inf
  for i in 1:num_vertices(f)
    max_dist = 0.0
    vertex = f[i]
    for j in 1:num_facets(f)
      facet = get_facet(f,j)
      dist = distance(facet,vertex)
      max_dist = max(dist,max_dist)
    end
    min_dist = min(max_dist,min_dist)
  end
  min_dist
end

function max_length(f::Face{Df}) where Df
  max_len = 0.0
  for e in 1:num_edges(f)
    edge = get_edge(f,e)
    max_len = max(max_len,length(edge))
  end
  max_len
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

# Distances

function signed_distance(p::Point,Π::Plane)
  o = origin( Π )
  n = normal( Π )
  (p-o) ⋅ n
end

function signed_distance(point::Point{D},Π::CartesianPlane{D}) where D
  Π.positive ? point[Π.d] - Π.value : Π.value - point[Π.d]
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

function distance(p::Point{Dp},f::Face{Df,Dp}) where {Df,Dp}
  if contains_projection(f,p)
    distance_to_infinite_face(f,p)
  else
    distance_to_boundary(f,p)
  end
end

distance(f::Face,p::Point) = distance(p,f)

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

distance_to_infinite_face(a::Face{D,D},b::Point{D}) where D = 0.0

function distance_to_boundary(a::Face,b::Point)
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

# Projections

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

function projection(p::Point,Π::Plane)
  n = normal(Π)
  c = origin(Π)
  p + ( ( c - p ) ⋅ n ) * n
end

function projection(p::Point,Π::CartesianPlane)
  Point( Base.setindex(Tuple(p),Π.value,Π.d) )
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

function intersection_point(f::Face{1,D},p::AbstractPlane{D}) where D
  v1 = f[1]
  v2 = f[2]
  d1 = signed_distance(v1,p)
  d2 = signed_distance(v2,p)
  if d1-d2 == 0
    d1,d2 = 1,10
  end
  (v1*d2-v2*d1)/(d2-d1)
end

intersection_point(a::Face{1},b::Face) = _intersection_point(a,b)

intersection_point(a::Face{1},b::Face{1}) = _intersection_point(a,b)

intersection_point(a::Face,b::Face{1}) = _intersection_point(b,a)

function _intersection_point(e::Face{1},f::Face)
  @notimplementedif num_dims(f) + 1 != num_point_dims(f)
  n = normal(f)
  c = center(f)
  plane = Plane(c,n)
  intersection_point(e,plane)
end

function expand_face(f::Face,dist::Real)
  Dc = num_dims(f)
  if is_simplex(f)
    verts = ntuple(i->_pull_vertex(f,i,dist),Val{Dc+1}())
    simplex_face( verts )
  else
    @notimplemented
  end
end

function _pull_vertex(f::Face,v::Integer,d::Real)
  dir = _vertex_direction(f,v)
  f[v] + dir*d
end

function _vertex_direction(f::Face,v::Integer)
  p = get_polytope(f)
  v_to_e = get_faces(p,0,1)
  e_to_v = get_faces(p,1,0)
  dir = zero(f[1])
  n = length(v_to_e[v])
  for e in v_to_e[v]
    for vneig in e_to_v[e]
      if vneig != v
        dv = f[v] - f[vneig]
        dir += dv / n
      end
    end
  end
  @assert norm(dir) > 0
  dir / norm(dir)
end


# General predicates

function has_intersection(a::Face,b::Face)
  Dp = num_point_dims(a)
  pa = get_polytope(a)
  pb = get_polytope(b)
  for da in 0:num_dims(a)
    db = Dp-da
    if num_dims(b) ≥ db
      for fa in 1:num_faces(pa,da), fb in 1:num_faces(pb,db)
        if has_intersection_point(a,da,fa,b,db,fb)
          return true
        end
      end
    end
  end
  false
end

function has_intersection_point(p::Point,f::Face{D,D}) where D
  contains_projection(f,p)
end

has_intersection_point(f::Face,p::Point) = has_intersection_point(p,f)

function has_intersection_point(a::Face,b::Face)
  p = intersection_point(a,b)
  contains_projection(a,p) && contains_projection(b,p)
end

function has_intersection_point(
  f::Face,
  d::Integer,
  df::Integer,
  args::Vararg{<:Any,N}) where N

  if d == 0
    f0 = get_dface(f,df,Val{0}())
    has_intersection_point(args...,f0)
  elseif d == 1
    f1 = get_dface(f,df,Val{1}())
    has_intersection_point(args...,f1)
  elseif d == 2
    f2 = get_dface(f,df,Val{2}())
    has_intersection_point(args...,f2)
  elseif d == 3
    f3 = get_dface(f,df,Val{3}())
    has_intersection_point(args...,f3)
  else
    @notimplemented
  end
end

# Voxel predicates

function voxel_intersection(f::Face{1,2},pmin::Point,pmax::Point,p::Polytope)
  D = num_point_dims(f)
  for i in 1:num_vertices(f)
    v = f[i]
    voxel_intersection(v,pmin,pmax) && return true
  end
  voxel_intersection(f,pmin,pmax)
end

function voxel_intersection(f::Face{2,3},pmin::Point,pmax::Point,p::Polytope)
  D = num_point_dims(f)
  for i in 1:num_vertices(f)
    v = f[i]
    voxel_intersection(v,pmin,pmax) && return true
  end
  for i in 1:num_edges(f)
    e = get_edge(f,i)
    voxel_intersection(e,pmin,pmax) && return true
  end
  measure(f) == 0 && return false
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

function voxel_intersection(p::Point,pmin::Point,pmax::Point)
  @check all( pmin .≤ pmax ) "pmin  > pmax"
  all( p .> pmin ) && all( p .< pmax )
end

function voxel_intersection(e::Face{1,D},pmin::Point{D},pmax::Point{D}) where D
  @check all( pmin .≤ pmax ) "pmin  > pmax"
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

function voxel_intersection(
  c::Face{D,D},
  pmin::Point{D},
  pmax::Point{D},
  p::Polytope) where D

  @check all( pmin .≤ pmax ) "pmin  > pmax"
  for i in 1:num_facets(c)
    f = get_facet(c,i)
    voxel_intersection(f,pmin,pmax,p) && return true
  end
  for i in 1:num_facets(c)
    f = get_facet(c,i)
    plane = Plane( center(f), normal(c,i) )
    v_to_dists = _compute_distances(plane,pmin,pmax,p)
    all( d -> d < 0, v_to_dists) || return false
  end
  true
end

function _compute_distances(
  Π::Plane,
  pmin::Point,
  pmax::Point,
  p::Polytope{D}) where D

  @assert is_n_cube(p)
  ntuple(i->signed_distance(_get_vertex(p,pmin,pmax,i),Π),Val{2^D}())
end

function _get_vertex(p::Polytope{D},pmin::Point,pmax::Point,i::Integer) where D
  Point(ntuple(d-> (i-1) & (1<<(d-1)) == 0 ? pmin[d] : pmax[d], Val{D}()))
end

function _is_edge_cut(p::Polytope,v_to_dists::Tuple,e::Integer)
  e_to_v = get_face_vertices(p,1)[e]
  sign(v_to_dists[e_to_v[1]]) ≠ sign(v_to_dists[e_to_v[2]])
end

function _intersection_point(
  p::Polytope,
  pmin::Point,
  pmax::Point,
  v_to_dists::Tuple,
  e::Integer)

  e_to_v = get_face_vertices(p,1)[e]
  d1 = v_to_dists[e_to_v[1]]
  d2 = v_to_dists[e_to_v[2]]
  v1 = _get_vertex(p,pmin,pmax,e_to_v[1])
  v2 = _get_vertex(p,pmin,pmax,e_to_v[2])
  (v1*d2-v2*d1)/(d2-d1)
end

# Compute planes

function bisector_plane(edge::Face{1,3},Π1::AbstractPlane,Π2::AbstractPlane)
  n1 = normal(Π1)
  n2 = normal(Π2)
  n1 ⋅ n2 ≉ -1 || error("Edge too sharp")
  n = n1-n2
  if norm(n) < 1
    v = edge[2]-edge[1]
    v /= norm(v)
    _n = n1+n2
    @assert norm(_n) > 1
    _n /= norm(_n)
    n = _n × v
    @assert norm(n) ≈ 1
  end
  n /= norm(n)
  @assert norm(n) ≈ 1
  c = center(edge)
  Plane(c,n)
end

function get_cell_planes(p::Polytope,pmin::Point,pmax::Point)
  @notimplementedif !is_n_cube(p)
  D = num_dims(p)
  N = 2*D
  @assert num_facets(p) == N
  Π_cell = lazy_map(
    i -> CartesianPlane(isodd(i)*pmin+iseven(i)*pmax,D-((i-1)÷2),1),
    1:N )
  Π_ids = - ( 1:N )
  Π_inout = lazy_map( iseven, 1:N )
  Π_cell, Π_ids, Π_inout
end

function get_cell_planes(p::Polytope,coords)
  if is_simplex(p)
    f = simplex_face(coords)
    get_cell_planes(f)
  elseif is_n_cube(p)
    pmin,pmax = coords[1],coords[end]
    get_cell_planes(p,pmin,pmax)
  end
end

function get_cell_planes(face::Face)
  n_facets = num_facets(face)
  facets = lazy_map(get_facet,Fill(face,n_facets),1:n_facets)
  centers = lazy_map(center,facets)
  normals = lazy_map(normal,Fill(face,n_facets),1:n_facets)
  mask = lazy_map(n-> n > -n, normals)
  facs = lazy_map(i->(-1)^(i+1),mask)
  _normals = lazy_map(*,normals,facs)
  planes = lazy_map(Plane,centers,_normals)
  ids = -(1:n_facets)
  planes,ids,mask
end

function displace(plane::Plane,dist,oriented=true)
  c0 = center(plane)
  n0 = normal(plane)
  n = oriented ? n0 : -n0
  c = c0 + n*dist
  Plane(c,n0)
end

function displace(plane::CartesianPlane,dist,oriented=true)
  c0 = origin(plane)
  d = plane.d
  or = plane.positive
  diff = or == oriented ? dist : -dist
  c = c0 + diff
  dir = or ? 1 : -1
  CartesianPlane(c,d,dir)
end

function expand_planes(
  planes::AbstractVector{<:AbstractPlane},
  inout::AbstractVector,
  dist)

  n_planes = length(planes)
  lazy_map(displace,planes,Fill(dist,n_planes),inout)
end

function _visualization_grid(face::Face)
  X = collect(face.vertices)
  T = Table( [1:length(X)] )
  ctype = [1]
  reffes = [LagrangianRefFE(get_polytope(face))]
  UnstructuredGrid(X,T,reffes,ctype)
end

function writevtk(face::Face,filename;kwargs...)
  grid = _visualization_grid(face)
  writevtk(grid,filename;kwargs...)
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
  length(a) == D-1 || error(
      "orthogonal(::VectorValue{D}...) only defined for D-1 VectorValues{D}'s")
  orthogonal(a,)
end

