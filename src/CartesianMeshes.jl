
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

num_dims(::CartesianMesh{D}) where D = D

function num_cells(m::CartesianMesh{D}) where D
  n = 1
  for d in 1:D
    n *= m.partition[d]
  end
  n
end

function get_cell(m::CartesianMesh{D},i::Integer) where D
  c_coords = cartesian_coordinates(m,i)
  x_min = get_data(m.origin) .+ get_data(m.sizes) .* (c_coords.-1) ./ m.partition
  x_max = get_data(m.origin) .+ get_data(m.sizes) .* c_coords ./ m.partition
  BoundingBox(Point(x_min),Point(x_max))
end

function find_container(m::CartesianMesh{D},p::Point{D}) where D
  c_coords = Int.(floor.( (get_data(p) .- get_data(m.origin)) .* m.partition ./ get_data(m.sizes) ))
  c_coords = c_coords .+ 1
  c_coords = max.(c_coords,1)
  c_coords = min.(c_coords,m.partition)
  get_cell_id(m,c_coords)
end

function cartesian_coordinates(m::CartesianMesh{D},n::Integer) where D
  CartesianIndices(m.partition)[n].I
end

function get_cell_id(m::CartesianMesh{D},n::NTuple{D,Integer}) where D
  LinearIndices(m.partition)[n...]
end

function cells_touching_bounding_box!(cache::Vector,m::CartesianMesh{D},p::Point{D}) where {D}
  id = find_container(m,p)
  int_coord = cartesian_coordinates(m,id)

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
  min_int_coord = cartesian_coordinates(m,min_id)
  min_int_coord = min_int_coord .- 1
  min_int_coord = max.(min_int_coord,1)

  max_id = find_container(m,bb.pmax)
  max_int_coord = cartesian_coordinates(m,max_id)
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
