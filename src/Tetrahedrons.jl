
struct Tetrahedron{D,T}
  vertices::Tuple{Point{D,T},Point{D,T},Point{D,T},Point{D,T}}
end

function Tetrahedron(p1::Point,p2::Point,p3::Point,p4::Point)
  Tetrahedron((p1,p2,p3,p4))
end

num_dims(::Tetrahedron{D}) where D = D

num_facets(::Type{<:Tetrahedron}) = 4

num_facets(::T) where T<:Tetrahedron = num_facets(T)

get_facet_to_vertices(::Type{<:Tetrahedron}) = ((1,2,3),(1,4,2),(1,3,4),(2,4,3))

get_facet_to_vertices(::T) where T<:Tetrahedron = get_facet_to_vertices(T)

@inline get_vertices(t::Tetrahedron) = t.vertices

function get_facet(t::Tetrahedron,i::Integer)
  facet_to_vertices = get_facet_to_vertices(t)
  lpoints = facet_to_vertices[i]
  Triangle( t.vertices[lpoints[1]],  t.vertices[lpoints[2]],  t.vertices[lpoints[3]] )
end

function volume(t::Tetrahedron{3})
  factor = 1/6
  v1 = t.vertices[2] - t.vertices[1]
  v2 = t.vertices[3] - t.vertices[1]
  v3 = t.vertices[4] - t.vertices[1]
  ( v1 × v2 ) ⋅ v3 * factor 
end

measure(t::Tetrahedron) = volume(t)

function measure(t::Tetrahedron{4})
  v1 = t.vertices[2] - t.vertices[1]
  v2 = t.vertices[3] - t.vertices[1]
  v3 = t.vertices[4] - t.vertices[1]
  norm( orthogonal(v1,v2,v3) )
end

function normal(t::Tetrahedron{4})
  v1 = t.vertices[2] - t.vertices[1]
  v2 = t.vertices[3] - t.vertices[1]
  v3 = t.vertices[4] - t.vertices[1]
  orthogonal(v1,v2,v3)
end

function have_intersection(p::Point{3},t::Tetrahedron{3})
  for i in 1:num_facets(t)
    facet = get_facet(t,i)
    n = normal(facet)
    c = center(facet)
    n = n / norm(n)
    c_p = p - c
    if c_p ⋅ n < 0
      return false
    end
  end
  true
end

have_intersection_point(p::Point{3},t::Tetrahedron{3}) = have_intersection(p,t)

contains_projection(p::Point{3},t::Tetrahedron{3}) = have_intersection(p,t)

function distance(p::Point{3},t::Tetrahedron{3})
  if have_intersection(p,t)
    distance = 0.0
  else
    distance = typemax(0.0)
  end
end

function closest_point(t::Tetrahedron{3},p::Point{3})
  p
end

