module TrianglesTests

using Test
using STLCutter
import STLCutter: num_dims, ×

const num_points_per_triangle = 3
struct Triangle{D}
  p::NTuple{num_points_per_triangle,Point{D}}
end

function Triangle(p1::Point,p2::Point,p3::Point)
  Triangle((p1,p2,p3))
end

@inline Base.getindex(t::Triangle,i::Integer) = t.p[i]

@inline get_vertices(t::Triangle) = t.p

num_dims(::Triangle{D}) where D = D

@inline function normal(t::Triangle{D}) where D
  v1 = t[2] - t[1]
  v2 = t[3] - t[1]
  v1 × v2
end

function distance(p::Point{D},t::Triangle) 

p1 = Point(0,0,0)
p2 = Point(1,0,0)
p3 = Point(0,1,0)

t = Triangle(p1,p2,p3)
@test isa(t,Triangle{3})

t = Triangle((p1,p2,p3))
@test isa(t,Triangle{3})
@test t[1] == p1
@test num_dims(t) == 3

n = normal(t)
@test n.data == (0,0,1)


# const ledge_to_triangle_point = ((1,2),(1,3),(2,3))

# t.p[ledge_to_triangle_point]
# function get_edges(t::Triangle)
#   Segment.(t.p[ledge_to_triangle_point])
# end


end # module
