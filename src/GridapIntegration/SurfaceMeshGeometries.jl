
struct SurfaceMeshGeometry <: Geometry
 surface_mesh::SurfaceMesh
end

get_tree(geo::SurfaceMeshGeometry) = Leaf( ( geo.surface_mesh, "surface_mesh", nothing )  )

function compatible_geometries(a::SurfaceMeshGeometry,b::SurfaceMeshGeometry)
  a,b
end

function similar_geometry(a::SurfaceMeshGeometry,tree::Leaf)
  sm, = tree.data
  SurfaceMeshGeometry(sm)
end

function SurfaceMeshGeometry(a::STLGeometry)
  tree = get_tree(a)
  stl, = tree.data
  sm = SurfaceMesh(a.stl)
  SurfaceMeshGeometry(sm)
end
  

