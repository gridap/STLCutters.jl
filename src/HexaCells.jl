
struct BoundingBox{D,T}
  pmin::Point{D,T}
  pmax::Point{D,T}
end

struct HexaCell{D,T}
  bb::BoundingBox{D,T}
end

function BoundingBox(p::Point{D,T}) where {D,T}
  BoundingBox{D,T}(p,p)
end

function BoundingBox(s::Segment{D,T}) where {D,T}
  BoundingBox{D,T}(min.(s.points...),max.(s.points...))
end

function BoundingBox(t::Triangle{D,T}) where {D,T}
  BoundingBox{D,T}(min.(t.points...),max.(t.points...))
end

function BoundingBox(t::Tetrahedron{D,T}) where {D,T}
  BoundingBox{D,T}(min.(t.p...),max.(t.p...))
end


function HexaCell(pmin::Point{D},pmax::Point{D}) where {D}
  for d in 1:D
    @assert pmin[d] <= pmax[d]
  end
  HexaCell(BoundingBox(pmin,pmax))
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

@inline have_intersection(p::Point,h::HexaCell) = have_intersection(p,h.bb)

function have_intersection(s::Segment{D},bb::BoundingBox{D}) where D
  bb = expand(bb,BB_tolerance)
  t = mutable(VectorValue{D,Float64})
  t_min = 0.0
  t_max = 1.0
  for d in 1:D
    p_d = s.points[1][d]
    v_d = s.points[2][d] - s.points[1][d]
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

@inline have_intersection(s::Segment,h::HexaCell) where{D} = have_intersection(s,h.bb)

function have_intersection(t::Triangle{D},bb::BoundingBox{D}) where {D}
  bb = expand(bb,BB_tolerance)
  have_intersection(bb,BoundingBox(t)) || return false
  for p ∈ vertices(t)
    if have_intersection(p,bb)
      return true
    end
  end
  for i ∈ 1:num_edges_per_triangle
    e = get_edge(t,i)
    if have_intersection(e,bb)
      return true
    end
  end
  t_pos = _fix_triangle(bb,t)
  n = abs(normal(t_pos))
  p0 = center(t_pos)
  main_d = max_dimension( n .* ( bb.pmax - bb.pmin ) )
  main_length = bb.pmax[main_d] - bb.pmin[main_d]

  i_int = 0
  sq = _square_vertices(bb,main_d)
  d = mutable(VectorValue{length(sq),Float64})
  for i in 1:length(sq)
    d[i] = ( n ⋅ (p0 -sq[i]) ) / n[main_d]
    if d[i] ≥ 0 || d[i] ≤ main_length
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
    have_intersection(intersection,t_pos)
  end
end

@inline have_intersection(t::Triangle{D},h::HexaCell{D}) where{D} = have_intersection(t,h.bb)

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

  points = mirror.(t.points)
  Triangle(points)
end

function _square_vertices(bb::BoundingBox{D},axis::Int) where D
  n_vertices = 2^(D-1)
  p_i = mutable(Point{D,Float64})
  v = mutable(VectorValue{n_vertices,Point{D,Float64}})
  for i in 1:4
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

function have_intersection(bb1::BoundingBox{D},bb2::BoundingBox{D}) where {D}
  !have_intersection(bb1.pmin,bb2) || return true
  !have_intersection(bb1.pmax,bb2) || return true
  !have_intersection(bb2.pmin,bb1) || return true
  !have_intersection(bb2.pmax,bb1) || return true
  false
end
