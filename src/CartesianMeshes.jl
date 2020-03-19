
struct CartesianMesh{D,T}
  origin::Point{D,T}
  sizes::VectorValue{D,T}
  partition::NTuple{D,Int}
end

function CartesianMesh(b::BoundingBox,partition::Tuple)
  origin = b.pmin
  sizes = b.pmax - b.pmin
  CartesianMesh(origin,sizes,partition)
end

function CartesianMesh(b::BoundingBox{D},n::Integer) where D
  p = tfill(n,Val{D}())
  CartesianMesh(b,p)
end

function get_reference_cell(::CartesianMesh{D,T}) where {D,T}
  BoundingBox( -one(Point{D,T}), one(Point{D,T}) )
end

num_dims(::CartesianMesh{D}) where D = D

function num_cells(m::CartesianMesh{D}) where D
  n = 1
  for d in 1:D
    n *= m.partition[d]
  end
  n
end

function num_vertices(m::CartesianMesh{D}) where D
  n = 1
  for d in 1:D
    n *= m.partition[d]+1
  end
  n
end

function get_cell(m::CartesianMesh,i::Integer)
  c_coords = cartesian_cell_coordinates(m,i)
  x_min = get_data(m.origin) .+ get_data(m.sizes) .* (c_coords.-1) ./ m.partition
  x_max = get_data(m.origin) .+ get_data(m.sizes) .* c_coords ./ m.partition
  BoundingBox(Point(x_min),Point(x_max))
end

function get_vertex_coordinates(m::CartesianMesh,i::Integer)
  v_coords = cartesian_vertex_coordinates(m,i)
  data = get_data(m.origin) .+ get_data(m.sizes) .* (v_coords.-1) ./ m.partition
  Point(data)
end

function cells_around_vertex!(cache::Vector,m::CartesianMesh,vertex::Integer)

  int_coord = cartesian_vertex_coordinates(m,vertex)

  min_int_coord = int_coord .- 1
  min_int_coord = max.(min_int_coord,1)

  max_int_coord = int_coord
  max_int_coord = min.(max_int_coord,m.partition)

  A=UnitRange.(min_int_coord,max_int_coord)
  resize!(cache,0)
  for i in CartesianIndices(A)
    push!(cache,get_cell_id(m,i.I))
  end
  cache
end


function cartesian_lvertex_coordinates(::CartesianMesh{D},lvertex::Integer) where D
  CartesianIndices( tfill(2,Val{D}()) )[lvertex].I
end

function get_vertex_id(mesh::CartesianMesh,cell::Integer,lvertex) 
  c_coords = cartesian_cell_coordinates(mesh,cell)
  lv_coords = cartesian_lvertex_coordinates(mesh,lvertex)
  v_coords = c_coords .+ lv_coords .- 1
  get_vertex_id(mesh,v_coords)
end

num_vertices_per_cell(::CartesianMesh{D}) where D = 2^D

function find_container(m::CartesianMesh{D},p::Point{D}) where D
  c_coords = Int.(floor.( (get_data(p) .- get_data(m.origin)) .* m.partition ./ get_data(m.sizes) ))
  c_coords = c_coords .+ 1
  c_coords = max.(c_coords,1)
  c_coords = min.(c_coords,m.partition)
  get_cell_id(m,c_coords)
end

function cartesian_vertex_coordinates(m::CartesianMesh{D},n::Integer) where D
  CartesianIndices(m.partition.+1)[n].I
end

function get_vertex_id(m::CartesianMesh{D},n::NTuple{D,Integer}) where D
  LinearIndices(m.partition.+1)[n...]
end

function cartesian_cell_coordinates(m::CartesianMesh{D},n::Integer) where D
  CartesianIndices(m.partition)[n].I
end

function get_cell_id(m::CartesianMesh{D},n::NTuple{D,Integer}) where D
  LinearIndices(m.partition)[n...]
end

function cells_touching_bounding_box!(cache::Vector,m::CartesianMesh{D},p::Point{D}) where {D}
  id = find_container(m,p)
  int_coord = cartesian_cell_coordinates(m,id)

  min_int_coord = int_coord .- 1
  min_int_coord = max.(min_int_coord,1)

  max_int_coord = int_coord .+ 1
  max_int_coord = min.(max_int_coord,m.partition)

  A=UnitRange.(min_int_coord,max_int_coord)
  resize!(cache,0)
  for i in CartesianIndices(A)
    push!(cache,get_cell_id(m,i.I))
  end
  cache
end

function cells_touching_bounding_box!(cache::Vector,m::CartesianMesh{D},bb::BoundingBox{D}) where {D}
  min_id = find_container(m,bb.pmin)
  min_int_coord = cartesian_cell_coordinates(m,min_id)
  min_int_coord = min_int_coord .- 1
  min_int_coord = max.(min_int_coord,1)

  max_id = find_container(m,bb.pmax)
  max_int_coord = cartesian_cell_coordinates(m,max_id)
  max_int_coord = max_int_coord .+ 1
  max_int_coord = min.(max_int_coord,m.partition)

  A=UnitRange.(min_int_coord,max_int_coord)
  resize!(cache,0)
  for i in CartesianIndices(A)
    push!(cache,get_cell_id(m,i.I))
  end
  cache
end

function all_to_all_compute_cell_to_surface_mesh_faces(m::CartesianMesh{D},sm::SurfaceMesh{D}) where D
  cells = Int[]
  faces = Int[]
  for k in 1:num_cells(m)
    cell = get_cell(m,k)
    for sm_face in 1:num_faces(sm)
      if have_intersection(cell,sm,sm_face)
        push!(cells,k)
        push!(faces,sm_face)
      end
    end
  end
  Table(faces,cells,num_cells(m))
end

function compute_cell_to_surface_mesh_faces(m::CartesianMesh{D},sm::SurfaceMesh{D}) where D
  cells = Int[]
  faces = Int[]
  cell_cache = Int[]
  for sm_face in 1:num_faces(sm)
    bb = BoundingBox(sm,sm_face)
    for k in cells_touching_bounding_box!(cell_cache,m,bb)
      cell = get_cell(m,k)
      if have_intersection(cell,sm,sm_face)
        push!(cells,k)
        push!(faces,sm_face)
      end
    end
  end
  Table(faces,cells,num_cells(m))
end

function writevtk(m::CartesianMesh{D,T},file_base_name) where {D,T}
  d_to_vtk_type_id = Dict(0=>1,1=>3,2=>9,3=>12)
  vtk_sorting_lvertices = [1,2,4,3,5,6,8,7]

  points = zeros(Float64,D,num_vertices(m))
  for i in 1:num_vertices(m)
    p = get_vertex_coordinates(m,i)
    for d in 1:D
      points[d,i] = p[d]
    end
  end
  
  cells = MeshCell{Vector{Int64}}[]
  vtk_type = VTKCellType(d_to_vtk_type_id[D])
  for k in 1:num_cells(m)
    vertices = zeros(Int,num_vertices_per_cell(m))
    for j in 1:num_vertices_per_cell(m)
      lv = vtk_sorting_lvertices[j]
      vertices[j] = get_vertex_id(m,k,lv)
    end
    push!( cells, MeshCell(vtk_type,vertices) )
  end

  vtkfile = vtk_grid(file_base_name,points,cells)

  vtk_save(vtkfile)
end
