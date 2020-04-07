module VolumeMeshesTests

using STLCutter


box = BoundingBox(
  Point( 0.2, 0.2, 0.2 ),
  Point( 0.7, 0.7, 0.7 ) )

bg_mesh = CartesianMesh( box, 2 )

VolumeMesh(bg_mesh)


# TODO: Complete test



struct Square{D,T} 
  vertices::NTuple{4,Point{D,T}}
end

struct Cube{D,T}
  vertices::NTuple{8,Point{D,T}}
end

Square(b::BoundingBox{2}) = Square(get_vertices(b))

Cube(b::BoundingBox{3}) = Cube(get_vertices(b))

function have_intersection(s::Square{3},Point{3})
end

function contains_projection(s::Square{3},Point{3})
end

function distance(s::Square{D},Point{D}) where D
end

function distance(s::Cube{D},Point{D}) where D
end


end # module
