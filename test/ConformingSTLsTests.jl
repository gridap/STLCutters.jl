module ConformingSTLsTests

using Test
using STLCutter

vertex_coordinates = Point{2,Float64}[ (0,0), (1,0), (1,1), (0,1) ]
d_face_to_vertices = [ TableOfVectors([[3,4], [1,2], [2,3]] )]
facet_normals = Point{2,Float64}[ (0,1), (0,-1), (1,0)]
d_face_to_facets = [TableOfVectors( [[ 1, 2, 5 ], [2, 6, 7]])]

stl = ConformingSTL(vertex_coordinates,d_face_to_vertices,facet_normals,d_face_to_facets)

@test stl.vertex_coordinates == vertex_coordinates
@test stl.d_face_to_vertices == d_face_to_vertices
@test stl.facet_normals == facet_normals
@test stl.d_face_to_facets == d_face_to_facets

stl = STLMesh("test/data/cube.stl")

struct VectorField{D,T} <: Number
  data::NTuple{D,T}
end

function VectorField(data::Vararg{D,T}) where {D,T}
  VectorField(data)
end

function Base.show(io::IO,v::VectorField)
  print(io,v.data)
end

function Base.show(io::IO,::MIME"text/plain",v::VectorField)
  print(io,typeof(v))
  print(io,v.data)
end

Base.length(::Type{<:VectorField{D}}) where D = D

Base.length(::VectorField{D}) where D = D

@inline Base.getindex(v::VectorField,i::Integer) = v.data[i]

Base.convert(::Type{VectorField{D,T}},v::NTuple{D}) where {D,T} = VectorField{D,T}(v)


using LinearAlgebra

function distance(p1::Point{D,T},p2::Point{D,T}) where {D,T}
  norm( collect(p2.data) - collect(p1.data) )
end



distance( stl.vertex_coordinates[1], stl.vertex_coordinates[2] )

const tolerance = 1e-8

function num_vertices( stl::STLMesh )
  length(stl.vertex_coordinates)
end

function num_facets( stl::STLMesh )
  length(stl.facet_to_vertices)
end

function num_dims( stl::STLMesh )
  length(stl.vertex_coordinates[1])
end



vertices_map = Vector{Int}( 1:num_vertices(stl) )

for ( i, vi ) in enumerate( stl.vertex_coordinates ), ( j, vj ) in enumerate( stl.vertex_coordinates[i+1:end] )
  @test i ≤ j+i
  if ( vertices_map[i] == i )
    if ( distance( vi, vj ) < tolerance )
      vertices_map[j+i] = i
    end
  end
end


function clean_map( vertices_map )
  counter = 1
  for  i in 1:length(vertices_map)
    # global counter
    if ( i == vertices_map[i] )
       vertices_map[i] = counter
       counter  += 1
    else
      vertices_map[i] = vertices_map[vertices_map[i]]
    end
  end
 vertices_map
end

vertices_map = clean_map(vertices_map)



vertex_coordinates =  []
for ( i, v ) in enumerate( stl.vertex_coordinates )
  if vertices_map[i] == i
    push!( vertex_coordinates, v )
  end
end

facet_to_vertices = []

for i in 1:num_facets(stl)
   vertices = getlist(stl.facet_to_vertices,i)
   push!( facet_to_vertices, vertices_map[ vertices ])
end
#
@show facet_to_vertices

end # module
