
struct Segment{D,T}
  vertices::Tuple{Point{D,T},Point{D,T}}
end

Segment(a...) = Segment(a)

Base.getindex(a::Segment,i::Integer) = a.vertices[i]

center(a::Segment) = a[1]/2+a[2]/2

get_polytope(a::Segment) = SEGMENT

struct Triangle{D,T}
  vertices::Tuple{Point{D,T},Point{D,T},Point{D,T}}
end

Triangle(a...) = Triangle(a)

Base.getindex(a::Triangle,i::Integer) = a.vertices[i]

center(a::Triangle) = a[1]/3+a[2]/3+a[3]/3

get_polytope(a::Triangle) = TRI
  
get_vertex_coordinates(a::Triangle) = a.vertices

num_vertices(a::Triangle) = 3

function normal(a::Triangle)
  v1 = a[2]-a[1]
  v2 = a[3]-a[1]
  n = orthogonal(v1,v2)
  n / norm(n)
end

function get_cell(a::Triangle)
  K = 1:num_vertices(a)
  X = get_vertex_coordinates(a)
  p = get_polytope(a)
  K,X,p
end

function Base.setindex(p::Point,v,idx::Integer)
  data = Base.setindex(p.data,v,idx)
  Point(data)
end

distance(a::Point,b::Point) =  norm( a - b )

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

distance(s::Segment,p::Point) = distance(p,s) 

function projection(p::Point{D},q::Point{D}) where D
  q
end

function projection(p::Point{D},s::Segment{D}) where D
  c = center(s)
  v = s[2] - s[1]
  v = v / norm(v)
  c + ( ( p - c ) ⋅ v ) * v
end

function distance(cell_nodes,node_to_coordinates,p::Polytope,d::Integer,dface::Integer,point::Point)
  face = get_dimrange(p,d)[dface]
  distance(cell_nodes,node_to_coordinates,p,face,point)
end

function distance(cell_nodes,node_to_coordinates,p::Polytope,face::Integer,point::Point)
  D = length(eltype(node_to_coordinates))
  d = get_facedims(p)[face]
  dface = face - get_dimrange(p,d)[1] + 1
  if d == 0
    vertex = node_to_coordinates[ cell_nodes[ dface ] ]
    distance(vertex,point)
  elseif d == 1
    dface_nodes = get_face_vertices(p,d)[dface]
    p1 = node_to_coordinates[ cell_nodes[ dface_nodes[1] ] ] 
    p2 = node_to_coordinates[ cell_nodes[ dface_nodes[2] ] ] 
    seg = Segment(p1,p2)
    distance(seg,point)
  elseif d == D-1
    c = center(cell_nodes,node_to_coordinates,p,face)
    n = normal(cell_nodes,node_to_coordinates,p,face)
    abs( (point-c) ⋅ n )
  else
    @assert false
  end
end

function projection(cell_nodes,node_to_coordinates,p::Polytope,d::Integer,dface::Integer,point::Point)
  face = get_dimrange(p,d)[dface]
  projection(cell_nodes,node_to_coordinates,p,face,point)
end

function projection(cell_nodes,node_to_coordinates,p::Polytope{D},face::Integer,point::Point) where D
  d = get_facedims(p)[face]
  dface = face - get_dimrange(p,d)[1] + 1
  if d == 0
    vertex = node_to_coordinates[ cell_nodes[ dface ] ]
    projection(point,vertex)
  elseif d == 1
    dface_nodes = get_face_vertices(p,d)[dface]
    p1 = node_to_coordinates[ cell_nodes[ dface_nodes[1] ] ] 
    p2 = node_to_coordinates[ cell_nodes[ dface_nodes[2] ] ] 
    seg = Segment(p1,p2)
    projection(point,seg)
  elseif d == D-1
    facet = dface
    c = center(cell_nodes,node_to_coordinates,p,face)
    n = normal(cell_nodes,node_to_coordinates,p,face)
    point + ( ( c - point ) ⋅ n ) * n
  else
    @assert false
  end
end

function get_bounding_box(
  cell_nodes::Vector{<:Integer},
  node_to_coordinates::Vector{<:Point})

  # check is a cartesian cell
  pmin = node_to_coordinates[ first(cell_nodes) ]
  pmax = node_to_coordinates[ last(cell_nodes) ]
  pmin,pmax
end

function have_intersection(cell_nodes,node_to_coordinates,::Polytope,point::Point)
  have_intersection(cell_nodes,node_to_coordinates,point)
end

function have_intersection(cell_nodes,node_to_coordinates,point::Point)
  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates)
  all( pmin.data .< point.data ) || return false
  all( pmax.data .> point.data ) || return false
  true
end

function have_intersection(cell_nodes,node_to_coordinates,p::Polytope,e::Segment)
  D = num_dims(p)
  tmin = 1.0
  tmax = 0.0
  for facet in 1:num_facets(p)
    face = get_dimrange(p,D-1)[facet]
    if have_intersection(cell_nodes,node_to_coordinates,p,face,e)
      t = compute_intersection_paremeter(cell_nodes,node_to_coordinates,p,face,e)
      tmin = min(t,tmin)
      tmax = max(t,tmax)
    end
  end
  tmin + TOL < tmax
end

function have_intersection(cell_nodes,node_to_coordinates,p::Polytope,face::Integer,e::Segment)
  D = length(eltype(node_to_coordinates))
  @assert get_facedims(p)[face] == D-1
  facet = face - get_dimrange(p,D-1)[1] + 1
  c = center(cell_nodes,node_to_coordinates,p,face)
  n = normal(cell_nodes,node_to_coordinates,p,face)
  s1_s2 = e[2] - e[1]
  s1_c = c - e[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  if α < -TOL || α > 1+TOL || isnan(α)
    false
  else
    x = e[1] + s1_s2 * α
    contains_projection(cell_nodes,node_to_coordinates,p,face,x)
  end
end

function contains_projection(cell_nodes,node_to_coordinates,p::Polytope,face::Integer,point::Point)
  d = get_facedims(p)[face]
  dface = face - get_dimrange(p,d)[1] + 1
  nfaces = get_faces(p,d,d-1)[dface]
  for (nlface,nface) ∈ enumerate(nfaces)
    c = center(cell_nodes,node_to_coordinates,p,d-1,nface)
    n = normal(cell_nodes,node_to_coordinates,p,d,dface,d-1,nlface)
    if ( point - c ) ⋅ n > TOL
      return false
    end
  end
  true
end

function compute_intersection_paremeter(cell_nodes,node_to_coordinates,p::Polytope,face::Integer,e::Segment)
  D = length(eltype(node_to_coordinates))
  @assert get_facedims(p)[face] == D-1
  facet = face - get_dimrange(p,D-1)[1] + 1
  c = center(cell_nodes,node_to_coordinates,p,face)
  n = normal(cell_nodes,node_to_coordinates,p,face)
  s1_s2 = e[2] - e[1]
  s1_c = c - e[1]
  ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
end

function intersection_point(cell_nodes,node_to_coordinates,p::Polytope,face::Integer,e::Segment)
  D = length(eltype(node_to_coordinates))
  @assert get_facedims(p)[face] == D-1
  facet = face - get_dimrange(p,D-1)[1] + 1
  @assert have_intersection(cell_nodes,node_to_coordinates,p,face,e)
  c = center(cell_nodes,node_to_coordinates,p,face)
  n = normal(cell_nodes,node_to_coordinates,p,face)
  s1_s2 = e[2] - e[1]
  s1_c = c - e[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  e[1] + s1_s2 * α
end


function center(cell_nodes,node_to_coordinates,p::Polytope,d::Integer,dface::Integer)
  face = get_dimrange(p,d)[dface]
  center(cell_nodes,node_to_coordinates,p,face)
end

function center(cell_nodes,node_to_coordinates,p::Polytope,face::Integer)
  d = get_facedims(p)[face]
  dface = face - get_dimrange(p,d)[1] + 1
  local_dface_nodes = get_faces(p,d,0)[dface]
  c = zero(eltype(node_to_coordinates))
  num_local_nodes = length(local_dface_nodes)
  for ln in local_dface_nodes
    n = cell_nodes[ln]
    vertex = node_to_coordinates[n]
    c += vertex/num_local_nodes
  end
  c
end

function normal(cell_nodes,node_to_coordinates,p::Polytope,d,dface::Integer,n,nlface::Integer)
  D = length(eltype(node_to_coordinates))
  @assert d == n+1
  if d == D
    facet = nlface
    normal(cell_nodes,node_to_coordinates,p,facet)
  elseif n == 0
    nface = get_faces(p,d,d-1)[dface][nlface]
    local_dface_nodes = get_face_vertices(p,d)[dface]
    p1 = node_to_coordinates[ cell_nodes[ local_dface_nodes[1] ] ]
    p2 = node_to_coordinates[ cell_nodes[ local_dface_nodes[2] ] ]
    v = p2 - p1
    v /= norm(v)
    if nface == local_dface_nodes[1] 
      -v
    else # nface == 2
      v
    end
  elseif d == D-1
    @assert D == 3 && n == 1
    if is_simplex(p)
      pf = TRI
      pf.face_orientations[3] = -1
    elseif is_n_cube(p)
      pf = QUAD
    else
      @assert false
    end
    facet = dface
    nface = get_faces(p,d,d-1)[dface][nlface]
    n_f = normal(cell_nodes,node_to_coordinates,p,D-1,facet)
    local_nface_nodes = get_face_vertices(pf,n)[nlface]
    local_facet_nodes = get_face_vertices(p,d)[facet]
    p1 = node_to_coordinates[ cell_nodes[ local_facet_nodes[ local_nface_nodes[1] ] ] ]
    p2 = node_to_coordinates[ cell_nodes[ local_facet_nodes[ local_nface_nodes[2] ] ] ]
    v = p2 - p1
    v = v / norm(v)
    v = v * pf.face_orientations[nlface] * p.face_orientations[facet]
    n_f × v
  else
    @assert false
  end

end

function normal(cell_nodes,node_to_coordinates,p::Polytope,d::Integer,dface::Integer)
  face = get_dimrange(p,d)[dface]
  normal(cell_nodes,node_to_coordinates,p,face)
end

function normal(cell_nodes,node_to_coordinates,p::Polytope,face::Integer)
  D = length(eltype(node_to_coordinates))
  @assert get_facedims(p)[face] == D-1
  facet = face - get_dimrange(p,D-1)[1] + 1
  function get_vertex(i)
    ln = local_facet_nodes[i]
    node = cell_nodes[ln]
    node_to_coordinates[node]
  end
  function get_vector(i)
    p_0 = get_vertex(1)
    p_i = get_vertex(i+1)
    p_i - p_0
  end

  local_facet_nodes = get_faces(p,D-1,0)[facet]
  @assert length(local_facet_nodes) ≥ D
  vectors = ntuple(get_vector,Val{D-1}())
  n = orthogonal(vectors...)
  n = n / norm(n)
  n * p.face_orientations[facet]
end


function have_intersection(cell_nodes,node_to_coordinates,p::Polytope,tri::Triangle)
  p1 = nothing
  p2 = nothing
  s = nothing
  D = num_dims(p)
  for dface in 1:num_faces(p,D-2)
    face = get_dimrange(p,D-2)[dface]
    if have_intersection(cell_nodes,node_to_coordinates,p,face,tri)
      point = intersection_point(cell_nodes,node_to_coordinates,p,face,tri)
      if p1 == nothing
        p1 = point
      elseif p2 === nothing
        if distance(p1,point) > TOL
          p2 = point
          s = Segment(p1,p2)
        end
      elseif distance(s,point) > TOL
        return true
      end

    end
  end
  false
end


function have_intersection(cell_nodes,node_to_coordinates,p::Polytope{3},face::Integer,t::Triangle{3})
  d = get_facedims(p)[face]
  @assert d == 1
  edge = face - get_dimrange(p,d)[1] + 1
  edge_nodes = get_face_vertices(p,d)[edge]
  p1 = node_to_coordinates[ cell_nodes[ edge_nodes[1] ] ]
  p2 = node_to_coordinates[ cell_nodes[ edge_nodes[2] ] ]
  s = Segment(p1,p2)
  have_intersection(t,s)
end

function have_intersection(t::Triangle,s::Segment)
  Kt,Xt,pt = get_cell(t)
  tface = num_faces(pt)
  have_intersection(Kt,Xt,pt,tface,s)
end

function intersection_point(cell_nodes,node_to_coordinates,p::Polytope{3},face::Integer,t::Triangle{3})
  d = get_facedims(p)[face]
  @assert d == 1
  edge = face - get_dimrange(p,d)[1] + 1
  edge_nodes = get_face_vertices(p,d)[edge]
  p1 = node_to_coordinates[ cell_nodes[ edge_nodes[1] ] ]
  p2 = node_to_coordinates[ cell_nodes[ edge_nodes[2] ] ]
  s = Segment(p1,p2)
  intersection_point(t,s)
end

function intersection_point(t::Triangle,s::Segment)
  Kt,Xt,pt = get_cell(t)
  tface = num_faces(pt)
  intersection_point(Kt,Xt,pt,tface,s)
end

function compute_intersection_paremeters(
  cell_nodes,
  node_to_coordinates,
  p::Polytope{3},
  face::Integer,
  t::Triangle{3})

  point = intersection_point(cell_nodes,node_to_coordinates,p,face,t)
  v1 = t[2]-t[1]
  v2 = t[3]-t[1]
  v1 /= norm(v1)
  v2 /= norm(v1)
  s = (point-t[1])⋅v1
  t = (point-t[1])⋅v2
  s,t
end

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

