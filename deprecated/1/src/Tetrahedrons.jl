
struct Tetrahedron{D,T}
  vertices::Tuple{Point{D,T},Point{D,T},Point{D,T},Point{D,T}}
end

function Tetrahedron(p1::Point,p2::Point,p3::Point,p4::Point)
  Tetrahedron((p1,p2,p3,p4))
end

function Tetrahedron(p::Point,t::Triangle)
  Tetrahedron(p,get_vertices(t)...)
end

function Tetrahedron(t::Triangle,p::Point)
  Tetrahedron(get_vertices(t)...,p)
end

num_dims(::Tetrahedron{D}) where D = D

num_vertices(::Type{<:Tetrahedron}) = 4

num_vertices(::T) where T<:Tetrahedron = num_vertices(T)

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

function volume(t::Tetrahedron{3}) #TODO: abs()
  factor = 1/6
  v1 = t.vertices[2] - t.vertices[1]
  v2 = t.vertices[3] - t.vertices[1]
  v3 = t.vertices[4] - t.vertices[1]
  ( v1 × v2 ) ⋅ v3 * factor
end

measure(t::Tetrahedron) = volume(t)

function measure(t::Tetrahedron{4})
  factor = 1/6
  v1 = t.vertices[2] - t.vertices[1]
  v2 = t.vertices[3] - t.vertices[1]
  v3 = t.vertices[4] - t.vertices[1]
  norm( orthogonal(v1,v2,v3) ) * factor
end

function signed_measure(t::Tetrahedron{3})
  factor = 1/6
  v1 = t.vertices[2] - t.vertices[1]
  v2 = t.vertices[3] - t.vertices[1]
  v3 = t.vertices[4] - t.vertices[1]
  det(v1,v2,v3) * factor 
end

function measure_sign(t::Tetrahedron{3})
  sign(signed_measure(t))
end

function min_height(t::Tetrahedron)
  factor = 3
  max_surf = 0.0
  for i in 1:num_facets(t)
    f = get_facet(t,i)
    surf = measure(f)
    if surf > max_surf
      max_surf = surf
    end
  end
  abs(measure(t))/max_surf*factor
end

function normal(t::Tetrahedron{4})
  v1 = t.vertices[2] - t.vertices[1]
  v2 = t.vertices[3] - t.vertices[1]
  v3 = t.vertices[4] - t.vertices[1]
  v1 = v1 / norm(v1)
  v2 = v2 / norm(v2)
  v3 = v3 / norm(v3)
  n = orthogonal(v1,v2,v3)
  n / norm(n)
end

function have_intersection(p::Point{3},t::Tetrahedron{3})
  for i in 1:num_facets(t)
    facet = get_facet(t,i)
    n = normal(facet)
    c = center(facet)
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

function relative_orientation(tri::Triangle{3},tet::Tetrahedron{3})
  max_distance = 0.0
  _v = get_vertices(tet)[1]
  for v in get_vertices(tet)
    dist = distance_to_plane(v,tri)
    if dist ≥ max_distance
      max_distance = dist
      _v = v
    end
  end
  _tet = Tetrahedron(tri,_v)
  - measure_sign(_tet)
end

function writevtk(t::Tetrahedron{D,T},file_base_name) where {D,T}
  vtk_type_id = 10

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
