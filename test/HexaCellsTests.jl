module HexaCellsTests

using Test
using STLCutter

struct BoundingBox{D}
  pmin::Point{D}
  pmax::Point{D}
end

struct HexaCell{D}
  bb::BoundingBox{D}
end

function HexaCell(pmin::Point{D},pmax::Point{D}) where {D}
  for d in 1:D
    @assert pmin[d] <= pmax[d]
  end
  HexaCell(BoundingBox(pmin,pmax))
end

function do_intersect(p::Point{D},bb::BoundingBox{D}) where D
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

@inline do_intersect(p::Point,h::HexaCell) = do_intersect(p,h.bb)

function do_intersect(s::Segment{D},bb::BoundingBox{D}) where D
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
        t_min = t
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


@inline do_intersect(s::Segment,h::HexaCell) where{D} = do_intersect(s,h.bb)


@test do_intersect(Point(1,2,3),BoundingBox(Point(1,1,1),Point(3,3,3)))
@test !do_intersect(Point(1,2,3),BoundingBox(Point(2,2,2),Point(3,3,3)))
@test do_intersect(Point(1,2,3),HexaCell(Point(1,1,1),Point(3,3,3)))
@test !do_intersect(Point(1,2,3),HexaCell(Point(2,2,2),Point(3,3,3)))

@test do_intersect(Segment(Point(0.5,0.5,0.5),Point(1, 1, 1)),BoundingBox(Point(0,0,0),Point(1,1,1)))
@test do_intersect(Segment(Point(0.0,0.0,0.5),Point(0, 0, 1)),BoundingBox(Point(0,0,0),Point(1,1,1)))
@test !do_intersect(Segment(Point(0,2,2),Point(1, 2, 2)),BoundingBox(Point(0,0,0),Point(1,1,1)))
@test do_intersect(Segment(Point(0.5,0.5,0.5),Point(1, 1, 1)),HexaCell(Point(0,0,0),Point(1,1,1)))


# function do_intersect(t::Triangle{D},bb::BoundingBox{D}) where {D}
#   for p ∈ get_vertices(t)
#     if do_intersect(t.p[1],bb)
#       return true
#     end
#   end
#   for e ∈ get_edges(t)
#     if do_intersect(Segment(t.p[1],t.p[2]),bb)
#       return true
#     end
#   end
#   n = get_normal(t)
#   p0 = collect(t[1].data)
#   n_scal = n .* ( bb.pmax - bb.pmin )
#   n_scal = n_scal / norm(n_scal)
#   main_dir = findfirst( x->x == maximum(n_scal), n_scal )
#   main_length = bb.pmax[main_dir] - bb.pmin[main_dir]
#   main_x = zeros(D)
#   main_x[main_dir] = 1
#
#   d = zeros(2)
#   p = [ collect(bb.pmin.data), collect(bb.pmax.data) ]
#
#   for i in 1:length(d)
#    d = ( ( n ⋅ p0 ) - ( n ⋅ p[i] ) ) / n[main_dir]
#   end
#
#   if d[1] < 0 || d[1] > 0
#     return false
#   else
#     @assert any( abs.(d) .≤ main_length )
#     i = findfirst( x->x ≤ main_length, abs.(d) )
#     intersection = p[i] + main_x * d[i]
#   end
#
# end
#
# function positivize_normal(bb::BoundingBox{D},t::Triangle{D}) where {D}
#    n = get_normal(t)
#    α = n .< 0
#    x = Vector{Point{D}}([])
#    for p ∈ get_vertices(t)
#      x_pos = p + α.*( -2*p + bb.pmin + bb.pmax )
#      push!(x,Point(x_pos...))
#    end
#    Triangle(x)
# end




#
#

#
# t=Triangle(Point(0.5,0.5,0.5), Point(0,0,1), Point(0,1,0))
# bb=BoundingBox(Point(1,1,1),Point(2,2,2))
# t_pos=positivize_normal(bb,t)
# @test any( get_normal(t) .< 0 )
# @test get_normal(t_pos) == abs.(get_normal(t))



end # module
