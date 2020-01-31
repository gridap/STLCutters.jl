module GeometriesTests

using Test
using STLCutter

struct Segment{D}
  p::NTuple{2,Point{D}}
end

struct Triangle{D}
  p::NTuple{3,Point{D}}
end

struct BoundingBox{D}
  pmin::Point{D}
  pmax::Point{D}
end

struct HexaCell{D}
  bb::BoundingBox{D}
end

function Segment(p1::Point,p2::Point)
  Segment((p1,p2))
end

@inline Base.getindex(s::Segment,i::Integer) = s.p[i]

@inline get_vertices(s::Segment) = s.p

function Triangle(p1::Point,p2::Point,p3::Point)
  Triangle((p1,p2,p3))
end

function Triangle(x::Vector{Point{D}}) where D
  Triangle(x...)
end

@inline Base.getindex(t::Triangle,i::Integer) = t.p[i]

@inline get_vertices(t::Triangle) = t.p

const edge_to_vertex = [[1,2],[1,3],[2,3]]

function get_edges(t::Triangle)
  Segment.(t.p[edge_to_vertex])
end

import Base.-
import Base.+
import Base.*

@inline -(p1::Point{D},p2::Point{D}) where {D} = collect(p1.data) .- collect(p2.data)
@inline -(p::Point,v::Vector) = collect(p.data) .- v
@inline -(v::Vector,p::Point) = v .- collect(p.data)

@inline +(p1::Point{D},p2::Point{D}) where {D} = collect(p1.data) .+ collect(p2.data)
@inline +(p::Point,v::Vector) = collect(p.data) .+ v
@inline +(v::Vector,p::Point) = p + v

@inline *(p::Point,α::Number) = collect(p.data) * α
@inline *(β::Number,p::Point) = p * β

using LinearAlgebra

function get_normal(t::Triangle{D}) where {D}
  v1 = t[2] - t[1]
  v2 = t[3] - t[1]
  v1 × v2
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

do_intersect(s::Segment,h::HexaCell) where{D} = do_intersect(s,h.bb)

function do_intersect(t::Triangle{D},bb::BoundingBox{D}) where {D}
  for p ∈ get_vertices(t)
    if do_intersect(t.p[1],bb)
      return true
    end
  end
  for e ∈ get_edges(t)
    if do_intersect(Segment(t.p[1],t.p[2]),bb)
      return true
    end
  end
  n = get_normal(t)
  p0 = collect(t[1].data)
  n_scal = n .* ( bb.pmax - bb.pmin )
  n_scal = n_scal / norm(n_scal)
  main_dir = findfirst( x->x == maximum(n_scal), n_scal )
  main_length = bb.pmax[main_dir] - bb.pmin[main_dir]
  main_x = zeros(D)
  main_x[main_dir] = 1

  d = zeros(2)
  p = [ collect(bb.pmin.data), collect(bb.pmax.data) ]

  for i in 1:length(d)
   d = ( ( n ⋅ p0 ) - ( n ⋅ p[i] ) ) / n[main_dir]
  end

  if d[1] < 0 || d[1] > 0
    return false
  else
    @assert any( abs.(d) .≤ main_length )
    i = findfirst( x->x ≤ main_length, abs.(d) )
    intersection = p[i] + main_x * d[i]
  end

end

function positivize_normal(bb::BoundingBox{D},t::Triangle{D}) where {D}
   n = get_normal(t)
   α = n .< 0
   x = Vector{Point{D}}([])
   for p ∈ get_vertices(t)
     x_pos = p + α.*( -2*p + bb.pmin + bb.pmax )
     push!(x,Point(x_pos...))
   end
   Triangle(x)
end













@test do_intersect(Point(1,2,3),BoundingBox(Point(1,1,1),Point(3,3,3)))
@test !do_intersect(Point(1,2,3),BoundingBox(Point(2,2,2),Point(3,3,3)))
@test do_intersect(Point(1,2,3),HexaCell(Point(1,1,1),Point(3,3,3)))
@test !do_intersect(Point(1,2,3),HexaCell(Point(2,2,2),Point(3,3,3)))


@test do_intersect(Segment(Point(0.5,0.5,0.5),Point(1, 1, 1)),BoundingBox(Point(0,0,0),Point(1,1,1)))
@test do_intersect(Segment(Point(0.0,0.0,0.5),Point(0, 0, 1)),BoundingBox(Point(0,0,0),Point(1,1,1)))
@test !do_intersect(Segment(Point(0,2,2),Point(1, 2, 2)),BoundingBox(Point(0,0,0),Point(1,1,1)))
@test do_intersect(Segment(Point(0.5,0.5,0.5),Point(1, 1, 1)),HexaCell(Point(0,0,0),Point(1,1,1)))


t=Triangle(Point(0.5,0.5,0.5), Point(0,0,1), Point(0,1,0))
bb=BoundingBox(Point(1,1,1),Point(2,2,2))
t_pos=positivize_normal(bb,t)
@test any( get_normal(t) .< 0 )
@test get_normal(t_pos) == abs.(get_normal(t))


#do_intersect(p::Point{D},t::Triangle{D}) where {D}




end # module
