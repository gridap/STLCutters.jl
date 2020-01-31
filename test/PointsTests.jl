module PointsTests

using Test
using STLCutter

a = Point(1,1,1)
b = Point(2,3,3)

@test distance(a,b) == 3

end # module
