
struct ConformingSTL{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  d_face_to_vertices::Vector{TableOfVectors{Int}}
  facet_normals::Vector{Point{D,T}}
  d_face_to_facets::Vector{TableOfVectors{Int}}
end

function get_dface_to_vertices(a::ConformingSTL,d::Integer)
  a.d_face_to_vertices[d+1]
end

function get_dface_to_facets(a::ConformingSTL,d::Integer)
  a.d_face_to_facets[d+1]
end

using FileIO
using MeshIO: decompose, Face, Point3, Normal, GeometryTypes

struct STLMesh{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  facet_to_vertices::TableOfVectors{Int}
  facet_normals::Vector{Point{D,T}}
end

Base.convert( ::Type{Point{D,T}}, x::GeometryTypes.Point{D} )  where {D,T} = Point{D,T}(x.data)
Base.convert( ::Type{Point{D,T}}, x::Normal{D} )  where {D,T} = Point{D,T}(x.data)
Base.convert( ::Type{Vector}, x::Face ) =  collect(x.data)

function STLMesh(filename::String)
  stl = load(filename)
  vertex_coordinates = Vector{Point{3,Float64}}( decompose(Point3{ Float32}, stl ) )
  facet_to_vertices  = TableOfVectors( Vector{Vector{Int}}( decompose(Face{3,Int32},stl.faces) ) )
  facet_normals      = Vector{Point{3,Float64}}( decompose(Normal{3, Float64}, stl) )
  STLMesh( vertex_coordinates, facet_to_vertices, facet_normals )
end
