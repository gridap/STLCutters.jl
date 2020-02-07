
struct BoundingBox{D}
  pmin::Point{D}
  pmax::Point{D}
end

struct HexaCell{D}
  bb::BoundingBox{D}
end

function BoundingBox(p::Point{D}) where D
  BoundingBox(p,p)
end

function BoundingBox(s::Segment{D}) where D
  BoundingBox(min_bound(s.p),max_bound(s.p))
end

function BoundingBox(t::Triangle{D}) where D
  BoundingBox(min_bound(t.p),max_bound(t.p))
end

function BoundingBox(t::Tetrahedron{D}) where D
  BoundingBox(min_bound(t.p),max_bound(t.p))
end

function HexaCell(pmin::Point{D},pmax::Point{D}) where {D}
  for d in 1:D
    @assert pmin[d] <= pmax[d]
  end
  HexaCell(BoundingBox(pmin,pmax))
end

const BB_tolerance = 1e-5
function expand(bb::BoundingBox,ε::Float64)
  d = bb.pmax - bb.pmin
  δ = d * ε
  BoundingBox( bb.pmin - δ, bb.pmax + δ )
end

function have_intersection(p::Point{D},bb::BoundingBox{D}) where D
  bb = expand(bb,BB_tolerance)
  for d in 1:D
    if p[d] < bb.pmin[d]
      return false
    end
    if p[d] > bb.pmax[d]
      return false
    end
  end
  true
end

@inline have_intersection(p::Point,h::HexaCell) = have_intersection(p,h.bb)

function have_intersection(s::Segment{D},bb::BoundingBox{D}) where D
  bb = expand(bb,BB_tolerance)
  t = mutable(VectorValue{D,Float64})
  t_min = 0.0
  t_max = 1.0
  for d in 1:D
    p_d = s.p[1][d]
    v_d = s.p[2][d] - s.p[1][d]
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
  t_pos = positivize_normal(bb,t)
  n = normal(t_pos)
  p0 = center(t_pos)
  n_scal = scal(n,( bb.pmax - bb.pmin ))
  n_scal = n_scal / norm(n_scal)
  main_d = max_dimension(n_scal)
  main_length = bb.pmax[main_d] - bb.pmin[main_d]

  d_pmin = ( n ⋅ ( p0 - bb.pmin ) ) / n[main_d]
  d_pmax = ( n ⋅ ( p0 - bb.pmax ) ) / n[main_d]
  if d_pmin < 0 || d_pmax > 0
    return false
  else
    main_x = cartesian_axis(n,main_d)
    if d_pmin ≤ main_length
      intersection = bb.pmin + main_x * d_pmin
    elseif d_pmax ≥ -main_length
      intersection = bb.pmax + main_x * d_pmax
    else
      @assert false
    end
    have_intersection(intersection,t_pos)
  end
end

@inline have_intersection(t::Triangle{D},h::HexaCell{D}) where{D} = have_intersection(t,h.bb)

function positivize_normal(bb::BoundingBox{D},t::Triangle{D}) where {D}
   x = mutable(VectorValue{num_points_per_triangle,typeof(mutable(VectorValue{D,Float64}))})
   n = normal(t)
   for i ∈ 1:num_points_per_triangle
     for d in 1:D
       if n[d] < 0
         x[i][d] = t[i][d] + ( bb.pmin[d] - t[i][d] ) + ( bb.pmax[d] - t[i][d] )
       else
         x[i][d] = t[i][d]
       end
     end
   end
   data = NTuple{num_points_per_triangle,VectorValue{D,Float64}}(x.data)
   Triangle(data)
end

function have_intersection(bb1::BoundingBox{D},bb2::BoundingBox{D}) where {D}
  !have_intersection(bb1.pmin,bb2) || return true
  !have_intersection(bb1.pmax,bb2) || return true
  !have_intersection(bb2.pmin,bb1) || return true
  !have_intersection(bb2.pmax,bb1) || return true
  false
end
