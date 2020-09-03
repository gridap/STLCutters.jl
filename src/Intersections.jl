
struct Segment{D,T}
  vertices::Tuple{Point{D,T},Point{D,T}}
end

struct Cell{D,T}
  cell_nodes::Vector{Int}
  node_to_coordinates::Vector{Point{D,T}}
  polytope::ExtrusionPolytope{D}
end

struct Face{D,T}
  cell::Cell{D,T}
  face::Int
end

Segment(a...) = Segment(a)

Base.getindex(a::Segment,i::Integer) = a.vertices[i]

center(a::Segment) = a[1]/2+a[2]/2

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

function distance(cell_nodes,node_to_coordinates,p::Polytope,d,dface,point::Point)
  if d == 0
    vertex = node_to_coordinates[ cell_nodes[ dface ] ]
    distance(vertex,point)
  elseif d == 1
    dface_nodes = get_face_vertices(p,d)[dface]
    p1 = node_to_coordinates[ cell_nodes[ dface_nodes[1] ] ] 
    p2 = node_to_coordinates[ cell_nodes[ dface_nodes[2] ] ] 
    face = Segment(p1,p2)
    distance(face,point)
  else
    @assert false
  end
end

function projection(cell_nodes,node_to_coordinates,p::Polytope,d,dface,point::Point)
  if d == 0
    vertex = node_to_coordinates[ cell_nodes[ dface ] ]
    projection(point,vertex)
  elseif d == 1
    dface_nodes = get_face_vertices(p,d)[dface]
    p1 = node_to_coordinates[ cell_nodes[ dface_nodes[1] ] ] 
    p2 = node_to_coordinates[ cell_nodes[ dface_nodes[2] ] ] 
    face = Segment(p1,p2)
    projection(point,face)
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
  for facet in 1:num_facets(p)
    if have_intersection(cell_nodes,node_to_coordinates,p,facet,e)
      return true
    end
  end
  false
end

function have_intersection(cell_nodes,node_to_coordinates,p::Polytope,facet::Integer,e::Segment)
  c = center(cell_nodes,node_to_coordinates,p,facet)
  n = normal(cell_nodes,node_to_coordinates,p,facet)
  s1_s2 = e[2] - e[1]
  s1_c = c - e[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  if α < 0 || α > 1 || isnan(α)
    false
  else
    x = e[1] + s1_s2 * α
    contains_projection(cell_nodes,node_to_coordinates,p,facet,x)
  end
end

function intersection_point(cell_nodes,node_to_coordinates,p::Polytope,facet::Integer,e::Segment)
  @assert have_intersection(cell_nodes,node_to_coordinates,p,facet,e)
  c = center(cell_nodes,node_to_coordinates,p,facet)
  n = normal(cell_nodes,node_to_coordinates,p,facet)
  s1_s2 = e[2] - e[1]
  s1_c = c - e[1]
  α = ( n ⋅ s1_c ) / ( n ⋅ s1_s2 )
  e[1] + s1_s2 * α
end

function contains_projection(cell_nodes,node_to_coordinates,p::Polytope,facet::Integer,point::Point)
  d = num_dims(p)-1
  dface = facet
  nfaces = get_faces(p,d,d-1)[dface]
  for (nlface,nface) ∈ enumerate(nfaces)
    c = center(cell_nodes,node_to_coordinates,p,d-1,nface)
    n = normal(cell_nodes,node_to_coordinates,p,d,dface,d-1,nlface)
    if ( point - c ) ⋅ n > 0
      return false
    end
  end
  true
end


function center(cell_nodes,node_to_coordinates,p::Polytope,facet::Integer)
  D = num_dims(p)
  center(cell_nodes,node_to_coordinates,p,D-1,facet)
end

function center(cell_nodes,node_to_coordinates,p::Polytope,d,dface::Integer)
  D = num_dims(p)
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
  D = num_dims(p)
  @assert d == n+1
  if d == D
    facet = nface
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
    facet = dface
    nface = get_faces(p,d,d-1)[dface][nlface]
    n_f = normal(cell_nodes,node_to_coordinates,p,facet)
    local_nface_nodes = get_face_vertices(QUAD,n)[nlface]
    local_facet_nodes = get_face_vertices(p,d)[facet]
    p1 = node_to_coordinates[ cell_nodes[ local_facet_nodes[ local_nface_nodes[1] ] ] ]
    p2 = node_to_coordinates[ cell_nodes[ local_facet_nodes[ local_nface_nodes[2] ] ] ]
    v = p2 - p1
    v = v / norm(v)
    v = v * QUAD.face_orientations[nlface] * p.face_orientations[facet]
    n_f × v
  else
    @assert false
  end

end

function normal(cell_nodes,node_to_coordinates,p::Polytope{D},facet::Integer) where D
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
  vectors = ntuple(get_vector,D-1)
  n = orthogonal(vectors...)
  n = n / norm(n)
  n * p.face_orientations[facet]
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

