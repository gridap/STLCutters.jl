module GeometriesTests

using Test
using STLCutter

struct Segment{D}
  p1::Point{D}
  p2::Point{D}
end

struct Triangle{D}
  p1::Point{D}
  p2::Point{D}
  p3::Point{D}
end

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

do_intersect(p::Point,h::HexaCell) = do_intersect(p,h.bb)


function do_intersect(s::Segment{D},bb::BoundingBox{D}) where D
  t = zeros(D)
  t_min = 0.0
  t_max = 1.0
  for d in 1:D
    p_d = s.p1[d]
    v_d = s.p2[d] - s.p1[d]
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

do_intersect(s::Segment,h::HexaCell) where{D} = do_intersect(s,h.bb)

function do_intersect(t::Triangle{D},bb::BoundingBox{D}) where {D}
  # for p ∈ get_vertices(t)
  if do_intersect(t.p1,bb)
    return true
  end
  # for e ∈ get_edges(t)
  if do_intersect(Segment(t.p1,t.p2),bb)
    return true
  end
  false
end

t=Triangle(Point(0,0,0), Point(0,0,1), Point(0,1,0))


@test do_intersect(Point(1,2,3),BoundingBox(Point(1,1,1),Point(3,3,3)))
@test !do_intersect(Point(1,2,3),BoundingBox(Point(2,2,2),Point(3,3,3)))
@test do_intersect(Point(1,2,3),HexaCell(Point(1,1,1),Point(3,3,3)))
@test !do_intersect(Point(1,2,3),HexaCell(Point(2,2,2),Point(3,3,3)))


@test do_intersect(Segment(Point(0.5,0.5,0.5),Point(1, 1, 1)),BoundingBox(Point(0,0,0),Point(1,1,1)))
@test do_intersect(Segment(Point(0.0,0.0,0.5),Point(0, 0, 1)),BoundingBox(Point(0,0,0),Point(1,1,1)))
@test !do_intersect(Segment(Point(0,2,2),Point(1, 2, 2)),BoundingBox(Point(0,0,0),Point(1,1,1)))
@test do_intersect(Segment(Point(0.5,0.5,0.5),Point(1, 1, 1)),HexaCell(Point(0,0,0),Point(1,1,1)))

end # module
