
struct Hexahedron{D,T}
  vertices::NTuple{8,Point{D,T}}
end

Hexahedron(v::Vararg{Point,8}) = Hexahedron(v,)

Hexahedron(b::BoundingBox{3}) = Hexahedron(get_vertices(b))

num_dims(::Hexahedron{D}) where D = D

num_facets(::Type{<:Hexahedron}) = 6

num_facets(::T) where T<:Hexahedron = num_facets(T)

get_facet_to_vertices(::Type{<:Hexahedron}) = 
  ( (1,2,3,4),
    (6,5,8,7),
    (2,1,6,5),
    (3,4,7,8),
    (1,3,5,7),
    (4,2,8,6) )

get_facet_to_vertices(::T) where T<:Hexahedron = get_facet_to_vertices(T)

@inline get_vertices(t::Hexahedron) = t.vertices

num_vertices(::Type{<:Hexahedron}) = 8

num_vertices(::T) where T<:Hexahedron = num_vertices(T)

function get_facet(h::Hexahedron,i::Integer)
  facet_to_vertices = get_facet_to_vertices(h)
  v = get_vertices(h)
  lpoints = facet_to_vertices[i]
  Quadrilater( v[lpoints[1]], v[lpoints[2]], v[lpoints[3]], v[lpoints[4]] )
end

function volume(h::Hexahedron{3})
  v = get_vertices(h)
  v1 = v[2] - v[1]
  v2 = v[3] - v[1]
  v3 = v[5] - v[1]
  ( v1 × v2 ) ⋅ v3 
end

measure(h::Hexahedron) = volume(h)

function have_intersection(p::Point{3},h::Hexahedron{3})
  for i in 1:num_facets(h)
    facet = get_facet(h,i)
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

have_intersection_point(p::Point{3},h::Hexahedron{3}) = have_intersection(p,h)

contains_projection(p::Point{3},h::Hexahedron{3}) = have_intersection(p,h)

function distance(p::Point{3},h::Hexahedron{3})
  if have_intersection(p,h)
    distance = 0.0
  else
    distance = typemax(0.0)
  end
end

distance(h::Hexahedron,p::Point) = distance(p,h)

projection(p::Point{3},h::Hexahedron{3}) = p

projection(h::Hexahedron,p::Point) = projection(p,h)

closest_point(h::Hexahedron{3},p::Point{3}) = projection(p,h)

function relative_orientation(tri::Triangle{3},hex::Hexahedron{3})
  max_distance = 0.0
  _v = get_vertices(hex)[1]
  for v in get_vertices(hex)
    dist = distance(v,tri)
    if dist ≥ max_distance
      max_distance = dist
      _v = v
    end
  end
  tet = Tetrahedron(tri,_v)
  - measure_sign(tet)
end

function writevtk(h::Hexahedron{3,T},file_base_name) where T
  vtk_type_id = 12
  vtk_hex_lvertices = [1,2,4,3,5,6,8,7]

  points = zeros(T,3,num_vertices(h))
  for (i,v) in enumerate(get_vertices(h)), d in 1:3
      points[d,i] = v[d]
  end

  vtk_type = VTKCellType(vtk_type_id)
  vertices = vtk_hex_lvertices[1:num_vertices(h)]
  cells = [ MeshCell(vtk_type,vertices) ]

  vtkfile = vtk_grid(file_base_name,points,cells)
  vtk_save(vtkfile)
end

function BoundingBox(h::Hexahedron{D,T}) where {D,T}
  BoundingBox{D,T}(min.(get_vertices(h)...),max.(get_vertices(h)...))
end

