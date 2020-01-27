module PointsTests

using Test
using STLCutter

data = (1,3)
p = Point( data )
@test string(p) == "(1, 3)"
@test p.data[2] == 3
@test p[2] == 3
@test length(p) == 2

x = Point(5,3)
@test x.data == (5,3)

@test length(Point{3,Int}) == 3

ps = Point{2,Float64}[(1,3),]
@test ps[1] == Point(1.0,3.0)

end # module
