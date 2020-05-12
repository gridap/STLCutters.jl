
struct STLGeometry <: Geometry
 stl::STL
end

get_tree(geo::STLGeometry) = Leaf( ( geo.stl, "stl", nothing )  )

function compatible_geometries(a::STLGeometry,b::STLGeometry)
  a,b
end

function similar_geometry(a::STLGeometry,tree::Leaf)
  stl, = tree.data
  STLGeometry(stl)
end
