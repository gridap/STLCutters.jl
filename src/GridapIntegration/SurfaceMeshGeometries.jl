
struct SurfaceMeshGeometry <: GridapEmbedded.CSG.Geometry
 stl::SurfaceMesh
end

get_tree(geo::SurfaceMeshGeometry) = GridapEmbedded.CSG.Leaf( ( geo.stl, "stl", nothing )  )

function compatible_geometries(a::SurfaceMeshGeometry,b::SurfaceMeshGeometry)
  a,b
end

function similar_geometry(a::SurfaceMeshGeometry,tree::Leaf)
  sm, = get_data(tree)
  SurfaceMeshGeometry(sm)
end
