
struct BoundingBox{D,T}
  pmin::Point{D,T}
  pmax::Point{D,T}
end

function BoundingBox(p::Point{D,T}) where {D,T}
  BoundingBox{D,T}(p,p)
end

function BoundingBox(s::Segment{D,T}) where {D,T}
  BoundingBox{D,T}(min.(get_vertices(s)...),max.(get_vertices(s)...))
end

function BoundingBox(t::Triangle{D,T}) where {D,T}
  BoundingBox{D,T}(min.(get_vertices(t)...),max.(get_vertices(t)...))
end

function BoundingBox(t::Tetrahedron{D,T}) where {D,T}
  BoundingBox{D,T}(min.(get_vertices(t)...),max.(get_vertices(t)...))
end

@generated function get_vertices(b::BoundingBox{D,T}) where {D,T}
  N = 2^D
  d = Dict( 0 => "pmin", 1 => "pmax" )
  v_str = [ "" for i in 1:N ]
  for i in 1:N
    bin = digits(i-1,base=2,pad=D)
    data = join([ "b.$(d[b])[$i]," for (i,b) in enumerate(bin) ])
    v_str[i] = "Point{$D,$T}($data),"
  end
  vertices = join(v_str)
  str = "($vertices)"
  Meta.parse(str)
end

num_vertices(::T) where T<:BoundingBox = num_vertices(T)

num_vertices(::Type{<:BoundingBox{D}}) where D = 2^D

num_dims(::Type{<:BoundingBox{D}}) where D = D

num_dims(::T) where T<:BoundingBox = num_dims(T)

function measure(b::BoundingBox)
  v = b.pmax - b.pmin
  prod(v)
end

const BB_tolerance = 1e-5
function expand(bb::BoundingBox,ε::Number)
  d = bb.pmax - bb.pmin
  δ = d * ε
  BoundingBox( bb.pmin - δ, bb.pmax + δ )
end

function have_intersection(p::Point{D},bb::BoundingBox{D}) where D
  bb = expand(bb,BB_tolerance)
  all( get_data( p .>= bb.pmin ) ) || return false
  all( get_data( p .<= bb.pmax ) ) || return false
  true
end

function have_intersection(s::Segment{D},bb::BoundingBox{D}) where D
  bb = expand(bb,BB_tolerance)
  t = mutable(VectorValue{D,Float64})
  t_min = 0.0
  t_max = 1.0
  for d in 1:D
    p_d = s[1][d]
    v_d = s[2][d] - s[1][d]
    if v_d < 0
      v_d = - v_d
      p_d = - p_d + bb.pmin[d] + bb.pmax[d]
    end
    if v_d != 0
      t = (  bb.pmin[d] - p_d ) / v_d
      if t > t_min
        t_min = t
      end
      t = (  bb.pmax[d] - p_d ) / v_d
      if t < t_max
        t_max = t
      end
    else
      if p_d < bb.pmin[d]
        return false
      end
      if p_d > bb.pmax[d]
        return false
      end
    end
  end
  if t_min > t_max
    false
  else
    true
  end
end

function have_intersection(t::Triangle{3},bb::BoundingBox{3})
  bb = expand(bb,BB_tolerance)
  for p ∈ get_vertices(t)
    if have_intersection(p,bb)
      return true
    end
  end
  for i ∈ 1:num_edges(t)
    e = get_edge(t,i)
    if have_intersection(e,bb)
      return true
    end
  end
  t_pos = _fix_triangle(bb,t)
  n = abs(normal(t_pos))
  n = n / norm(n)
  p0 = center(t_pos)
  main_d = max_dimension( n .* ( bb.pmax - bb.pmin ) )
  main_length = bb.pmax[main_d] - bb.pmin[main_d]

  i_int = 0
  sq = _square_vertices(bb,main_d)
  d = mutable(VectorValue{length(sq),Float64})
  for i in 1:length(sq)
    d[i] = ( n ⋅ (p0 -sq[i]) ) / n[main_d]
    if d[i] ≥ 0 && d[i] ≤ main_length
      i_int = i
    end
  end

  if d[1] < 0 || d[length(sq)] > main_length
    return false
  else
    main_x = canonical_vector(n,main_d)
    if i_int != 0
      intersection = sq[i_int] + main_x *d[i_int]
    else
      throw(ErrorException(""))
    end
    contains_projection(intersection,t_pos)
  end
end

function _fix_triangle(bb::BoundingBox{D},t::Triangle{D}) where {D}
  n = normal(t)

  function mirror(point::Point{D,T}) where {D,T}
    m = mutable(Point{D,T})
    for d in 1:D
      if n[d] < 0
        m[d] = point[d] + ( bb.pmin[d] - point[d] ) + ( bb.pmax[d] - point[d] )
      else
        m[d] = point[d]
      end
    end
    Point(m)
  end

  points = mirror.( get_vertices(t) )
  Triangle(points)
end

function _square_vertices(bb::BoundingBox{D},axis::Int) where D
  n_vertices = 2^(D-1)
  p_i = mutable(Point{D,Float64})
  v = mutable(VectorValue{n_vertices,Point{D,Float64}})
  for i in 1:n_vertices
    c = 2
    b = 0
    for j in 1:D
      if j ≠ axis
        b = (i-1)÷c - b*2
        p_i[j] = (1-b)*bb.pmin[j] + b*bb.pmax[j]
        c -= 1
      else
        p_i[j] = bb.pmin[j]
      end
    end
    v[i] = p_i
  end
  v.data
end

function have_intersection(t::Triangle{D},bb::BoundingBox{D}) where D
  throw(ArgumentError("have_intersection(::Triangle{$D},::BoundingBox{$D}) not implemented"))
end

function have_intersection(bb1::BoundingBox{D},bb2::BoundingBox{D}) where {D}
  throw(ArgumentError("have_intersection(::BoundingBox,::BoundingBox) not implemented"))
end

function transformation(dest_box::BoundingBox{D},origin_box::BoundingBox{D},p::Point{D}) where D
  dest_sizes = dest_box.pmax - dest_box.pmin
  origin_sizes = origin_box.pmax - origin_box.pmin
  ( p-origin_box.pmin ) .* ( dest_sizes ./ origin_sizes ) + dest_box.pmin
end

function transformation(
  dest_box::BoundingBox{D},
  origin_box::BoundingBox{D},
  origin_points::Vector{Point{D,T}}) where {D,T}
  
  dest_sizes = dest_box.pmax - dest_box.pmin
  origin_sizes = origin_box.pmax - origin_box.pmin
  frac_sizes = dest_sizes ./ origin_sizes
  dest_points = zeros(Point{D,T},length(origin_points))
  for (i,p) in enumerate( origin_points )
    dest_points[i] = ( p - origin_box.pmin ) .* frac_sizes + dest_box.pmin
  end
  dest_points
end

function writevtk(b::BoundingBox{D,T},file_base_name) where {D,T}
  d_to_vtk_hex_type_id = Dict(0=>1,1=>3,2=>9,3=>12)
  vtk_hex_lvertices = [1,2,4,3,5,6,8,7]

  points = zeros(T,D,num_vertices(b))
  for (i,v) in enumerate(get_vertices(b)), d in 1:D
      points[d,i] = v[d]
  end

  vtk_type = VTKCellType(d_to_vtk_hex_type_id[D])
  vertices = vtk_hex_lvertices[1:num_vertices(b)]
  cells = [ MeshCell(vtk_type,vertices) ]

  vtkfile = vtk_grid(file_base_name,points,cells)
  vtk_save(vtkfile)
end
  

