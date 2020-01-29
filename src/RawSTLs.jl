
import MeshIO

struct RawSTL{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  facet_to_vertices::TableOfVectors{Int}
  facet_normals::Vector{Point{D,T}}
end

function RawSTL(filename::String)
  stl = load(filename)
  ndims = num_dims(stl.vertices)
  nvxf  = num_vertices_per_facet(stl.faces)
  @assert num_dims(stl.vertices) == num_dims(stl.normals)
  vertex_coordinates = Vector{Point{ndims,Float64}}( MeshIO.decompose(MeshIO.Point{ndims,Float64}, stl ) )
  facet_to_vertices  = TableOfVectors{Int}(MeshIO.decompose(MeshIO.Face{nvxf,Int},stl.faces))
  facet_normals      = Vector{Point{ndims,Float64}}( MeshIO.decompose(MeshIO.Normal{ndims, Float64}, stl) )
  RawSTL( vertex_coordinates, facet_to_vertices, facet_normals )
end

num_dims(::Array{MeshIO.Point{D,T}}) where {D,T} = D

num_dims(::Array{MeshIO.Normal{D,T}}) where {D,T} = D

num_vertices_per_facet(::Array{MeshIO.Face{D,T}}) where{D,T} = D

function num_vertices( stl::RawSTL )
  length(stl.vertex_coordinates)
end

function num_facets( stl::RawSTL )
  length(stl.facet_to_vertices)
end

function num_dims( stl::RawSTL )
  length(stl.vertex_coordinates[1])
end

Base.convert( ::Type{Point{D,T}}, x::MeshIO.Point{D} )  where {D,T} = Point{D,T}(x.data)
Base.convert( ::Type{Point{D,T}}, x::MeshIO.Normal{D} )  where {D,T} = Point{D,T}(x.data)
Base.convert( ::Type{Vector}, x::MeshIO.Face ) = collect(x.data)

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
