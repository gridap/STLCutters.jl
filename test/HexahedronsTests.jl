module HexahedronsTests

using STLCutter
using Test

p1 = Point(1,2,3)
p2 = Point(4,5,6)

b = BoundingBox(p1,p2)

h = Hexahedron(b)

b2 = BoundingBox(h)

@test b2 == b

end # module
