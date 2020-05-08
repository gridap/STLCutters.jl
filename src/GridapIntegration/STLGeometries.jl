
struct STLGeometry <: GridapEmbedded.CSG.Geometry
 stl::STL
end

get_tree(geo::STLGeometry) = GridapEmbedded.CSG.Leaf( ( geo.stl, "stl", nothing )  )

function compatible_geometries(a::STLGeometry,b::STLGeometry)
  a,b
end

function similar_geometry(a::STLGeometry,tree::Leaf)
  stl, = get_data(tree)
  STLGeometry(stl)
end
