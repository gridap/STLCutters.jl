
import MeshIO: decompose, Face, Point3, Normal, GeometryTypes

struct RawSTL{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  facet_to_vertices::TableOfVectors{Int}
  facet_normals::Vector{Point{D,T}}
end

function RawSTL(filename::String)
  stl = load(filename)
  vertex_coordinates = Vector{Point{3,Float64}}( decompose(Point3{ Float32}, stl ) )
  facet_to_vertices  = TableOfVectors{Int}(decompose(Face{3,Int32},stl.faces))
  facet_normals      = Vector{Point{3,Float64}}( decompose(Normal{3, Float64}, stl) )
  RawSTL( vertex_coordinates, facet_to_vertices, facet_normals )
end

function num_vertices( stl::RawSTL )
  length(stl.vertex_coordinates)
end

function num_facets( stl::RawSTL )
  length(stl.facet_to_vertices)
end

function num_dims( stl::RawSTL )
  length(stl.vertex_coordinates[1])
end

Base.convert( ::Type{Point{D,T}}, x::GeometryTypes.Point{D} )  where {D,T} = Point{D,T}(x.data)
Base.convert( ::Type{Point{D,T}}, x::Normal{D} )  where {D,T} = Point{D,T}(x.data)
Base.convert( ::Type{Vector}, x::Face ) = collect(x.data)

const STL_tolerance = 1e-8
function map_repeated_vertices(stl::RawSTL)
  vertices_map = Vector{Int}(1:num_vertices(stl))
  for (i, vi) in enumerate(stl.vertex_coordinates),
      (j, vj) in enumerate(stl.vertex_coordinates[i+1:end])

    if (vertices_map[i] == i)
      if (distance(vi, vj) < STL_tolerance)
        vertices_map[j+i] = i
      end
    end
  end
  vertices_map
end

function extract_unique_vertices(stl::RawSTL, vertices_map::Vector{Int})
  vertex_coordinates = typeof(stl.vertex_coordinates)([])
  for (i, v) in enumerate(stl.vertex_coordinates)
    if vertices_map[i] == i
      push!(vertex_coordinates, v)
    end
  end
  vertex_coordinates
end

function compact_map(map::Vector{Int})
  counter = 1
  for i = 1:length(map)
    if (i == map[i])
      map[i] = counter
      counter += 1
    else
      map[i] = map[map[i]]
    end
  end
  map
end

function apply_map(v::TableOfVectors{Int}, map::Vector{Int})
  mapped_v = typeof(v)([])
  for i = 1:length(v)
    pushlist!(mapped_v, map[getlist(v, i)])
  end
  mapped_v
end
