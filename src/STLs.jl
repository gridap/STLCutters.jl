
struct STL{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  facet_to_vertices::Table{<:Integer}
  facet_normals::Vector{VectorValue{D,T}}
end

function STL(filename::String)
  stl = load(filename)
  ndims = num_dims(stl.vertices)
  nvxf  = num_vertices_per_facet(stl.faces)
  @assert num_dims(stl.vertices) == num_dims(stl.normals) "Corrupted STL file"
  vertex_coordinates = Vector{Point{ndims,Float64}}( MeshIO.decompose(MeshIO.Point{ndims,Float64}, stl ) )
  facet_to_vertices  = Table( Vector{Vector{Int}}(MeshIO.decompose(MeshIO.Face{nvxf,Int},stl.faces)))
  vertex_normals     = Vector{VectorValue{ndims,Float64}}( MeshIO.decompose(MeshIO.Normal{ndims, Float64}, stl) )
  facet_normals      = vertex_normals[1:nvxf:length(vertex_coordinates)]
  STL{ndims,Float64}( vertex_coordinates, facet_to_vertices, facet_normals )
end

function STL(vertex_coordinates::Vector{<:Point},facet_to_vertices::Table)
  facet_normals = _compute_facet_normals(vertex_coordinates,facet_to_vertices)
  STL(vertex_coordinates,facet_to_vertices,facet_normals)
end

function circular_points(vertex_coordinates::Vector{<:Point{2}})
  n = length(vertex_coordinates)
  data = zeros(Int,n*2)
  ptrs = fill(Int32(2),n+1)
  length_to_ptrs!(ptrs)
  f_to_v = Table(data,ptrs)
  display(f_to_v)
  for f in 1:length(f_to_v)
    for lv in 1:length(f_to_v,f)
      v = f+lv-1
      if v > n
        v = 1
      end
      f_to_v[f,lv] = v
    end
  end
  display(f_to_v)
  STL(vertex_coordinates,f_to_v)
end

num_dims(::Array{MeshIO.Point{D,T}}) where {D,T} = D

num_dims(::Array{MeshIO.Normal{D,T}}) where {D,T} = D

num_vertices_per_facet(::Array{MeshIO.Face{D,T}}) where{D,T} = D

function num_vertices( stl::STL )
  length(stl.vertex_coordinates)
end

function num_facets( stl::STL )
  length(stl.facet_to_vertices)
end

num_dims( stl::STL{D} ) where D = D

function Base.:(==)(a::STL,b::STL) 
 equal = (a.vertex_coordinates == b.vertex_coordinates) 
 equal = equal && (a.facet_to_vertices == b.facet_to_vertices) 
 equal && (a.facet_normals == b.facet_normals)
end

Base.convert( ::Type{Point{D,T}}, x::MeshIO.Point{D} )  where {D,T} = Point{D,T}(x.data)
Base.convert( ::Type{VectorValue{D,T}}, x::MeshIO.Normal{D} )  where {D,T} = VectorValue{D,T}(x.data)
Base.convert( ::Type{Vector}, x::MeshIO.Face ) = collect(x.data)


function delete_repeated_vertices!(stl::STL)
  old_to_new = identify_repeated_vertices(stl.vertex_coordinates)
  delete_repeated_vertex_coordinates!( stl.vertex_coordinates, old_to_new )
  apply_map!( stl.facet_to_vertices, old_to_new )
  stl
end

function delete_repeated_vertices(stl::STL)
  old_to_new = identify_repeated_vertices(stl.vertex_coordinates)
  vertex_coordinates = extract_unique_vertex_coordinates( stl.vertex_coordinates, old_to_new )
  facet_to_vertices = apply_map( stl.facet_to_vertices, old_to_new )
  facet_normals = deepcopy( stl.facet_normals )
  STL(vertex_coordinates,facet_to_vertices,facet_normals)
end

get_vertex_coordinates(stl::STL) = stl.vertex_coordinates

get_facet_normals(stl::STL) = stl.facet_normals

get_facet_to_vertices(stl::STL) = stl.facet_to_vertices

const STL_tolerance = 1e-8

function identify_repeated_vertices(vertex_coordinates::Vector{Point{D,T}}) where {D,T}
  bb = BoundingBox(vertex_coordinates)
  n_vertices = length(vertex_coordinates)
  n = Int(round(n_vertices ^ (1/D)))
  mesh = CartesianMesh(bb,n)
  cell_to_stl_vertices = [ Int[] for i in 1:num_cells(mesh) ]
  cell_cache = Int[]
  for (i,v) ∈ enumerate(vertex_coordinates)
    for k ∈ cells_touching_bounding_box!(cell_cache,mesh,v)
      h = get_cell(mesh,k)
      if have_intersection(v,h)
        push!(cell_to_stl_vertices[k],i)
      end
    end
  end
  vertices_map = Vector{Int}(1:n_vertices)
  for k in 1:num_cells(mesh)
    list = cell_to_stl_vertices[k]
    for i in 1:length(list), j in i+1:length(list) 
      vi = vertex_coordinates[list[i]]
      vj = vertex_coordinates[list[j]]
      if vertices_map[list[i]] == list[i] && vertices_map[list[j]] == list[j]
        if distance(vi,vj) < STL_tolerance
          vertices_map[list[j]] = list[i]
        end
      end
    end
  end
  _enumerate_map( vertices_map )
end

function _enumerate_map(map::Vector{<:Integer})
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

function delete_repeated_vertex_coordinates!(vertex_coordinates::Vector, old_to_new_vertices::Vector)
  count = 1
  for (i, v) in enumerate(vertex_coordinates)
    if old_to_new_vertices[i] == count
      vertex_coordinates[count] = v
      count += 1
    end
  end
  resize!(vertex_coordinates,count-1)
end

function extract_unique_vertex_coordinates(vertex_coordinates::Vector, old_to_new_vertices::Vector)
  T = eltype(vertex_coordinates)
  unique_vertex_coordinates = T[]
  count = 1
  for (i, v) in enumerate(vertex_coordinates)
    if old_to_new_vertices[i] == count
      push!(unique_vertex_coordinates,v)
      count += 1
    end
  end
  unique_vertex_coordinates
end

function apply_map(v::Table{<:Integer},old_to_new::Vector)
  v_copy = deepcopy(v)
  apply_map!(v_copy,old_to_new)
end

function apply_map!(v::Table{<:Integer},old_to_new::Vector)
  for i = 1:length(v), j = 1:length(v,i)
    v[i,j] = old_to_new[ v[i,j] ]
  end
  v
end

BoundingBox(stl::STL) = BoundingBox(stl.vertex_coordinates)

function BoundingBox(vertex_coordinates::Vector)
  pmin = vertex_coordinates[1]
  pmax = vertex_coordinates[1]
  for v ∈ vertex_coordinates
    pmin = min.(pmin,v)
    pmax = max.(pmax,v)
  end
  BoundingBox(pmin,pmax)
end
