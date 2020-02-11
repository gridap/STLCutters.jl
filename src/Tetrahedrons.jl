
struct Tetrahedron{D,T}
  points::Tuple{Point{D,T},Point{D,T},Point{D,T},Point{D,T}}
end

function Tetrahedron(p1::Point,p2::Point,p3::Point,p4::Point)
  Tetrahedron((p1,p2,p3,p4))
end

const num_facets_per_tetrahedron = 4
const lfacet_to_tetrahedra_points = ((1,2,3),(1,4,2),(1,3,4),(2,4,3))

num_dims(::Tetrahedron{D}) where D = D

@inline vertices(t::Tetrahedron) = t.points

function get_facet(t::Tetrahedron,i::Integer)
  lpoints = lfacet_to_tetrahedra_points[i]
  Triangle( t.points[lpoints[1]],  t.points[lpoints[2]],  t.points[lpoints[3]] )
end

function volume(t::Tetrahedron{D}) where D
  factor = 1/6
  v1 = t.points[2] - t.points[1]
  v2 = t.points[3] - t.points[1]
  v3 = t.points[4] - t.points[1]
  ( v1 × v2 ) ⋅ v3 * factor 
end

function have_intersection(p::Point{D},t::Tetrahedron{D}) where D
  for i in 1:num_facets_per_tetrahedron
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
