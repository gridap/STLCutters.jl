
struct CartesianMesh{D,T}
  origin::Point{D,T}
  sizes::VectorValue{D,T}
  partition::NTuple{D,Int}
end

function CartesianMesh(b::BoundingBox,p::Tuple)
  origin = b.pmin
  sizes = b.pmax - b.pmin
  CartesianMesh(origin,sizes,p)
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
  n_coords = int_coordinates(m,i)
  x_min = get_data(m.origin) .+ get_data(m.sizes) .* (n_coords.-1) ./ m.partition
  x_max = get_data(m.origin) .+ get_data(m.sizes) .* n_coords ./ m.partition
  BoundingBox(Point(x_min),Point(x_max))
end

function find_container(m::CartesianMesh{D},p::Point{D}) where D
  i_coords = Int.(floor.( (get_data(p) .- get_data(m.origin)) .* m.partition ./ get_data(m.sizes) ))
  i_coords = i_coords .+ 1
  i_coords = max.(i_coords,1)
  i_coords = min.(i_coords,m.partition)
  get_cell_id(m,i_coords)
end

function int_coordinates(m::CartesianMesh{D},n::Integer) where D
  n_d = mutable(VectorValue{D,Int})
  p_d = 1
  p = m.partition
  for d in 1:D
    n_d[d] = ( (n-1) ÷ p_d ) % p[d] + 1
    p_d *= p[d]
  end
  get_data(n_d)
end

function get_cell_id(m::CartesianMesh{D},n::NTuple{D,Integer}) where D
  gid = 0
  p_d = 1
  for d in 1:D
    gid += (n[d]-1)*p_d
    p_d *= m.partition[d]
  end
  gid + 1
end

function cells_around!(cache::Vector,m::CartesianMesh{D},p::Point{D}) where {D}
  id = find_container(m,p)
  int_coord = int_coordinates(m,id)

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

function cells_around!(cache::Vector,m::CartesianMesh{D},bb::BoundingBox{D}) where {D}
  min_id = find_container(m,bb.pmin)
  min_int_coord = int_coordinates(m,min_id)
  min_int_coord = min_int_coord .- 1
  min_int_coord = max.(min_int_coord,1)

  max_id = find_container(m,bb.pmax)
  max_int_coord = int_coordinates(m,max_id)
  max_int_coord = max_int_coord .+ 1
  max_int_coord = min.(max_int_coord,m.partition)

  A=UnitRange.(min_int_coord,max_int_coord)
  resize!(cache,0)
  for i in CartesianIndices(A)
    push!(cache,get_cell_id(m,i.I))
  end
  cache
end

function all_to_all_compute_cell_to_stl_faces(m::CartesianMesh{D},stl::SurfaceMesh{D}) where D
  cells = Int[]
  faces = Int[]
  for k in 1:num_cells(m)
    cell = get_cell(m,k)
    for stl_face in 1:num_faces(stl)
      if have_intersection(cell,stl,stl_face)
        push!(cells,k)
        push!(faces,stl_face)
      end
    end
  end
  Table(faces,cells,num_cells(m))
end

function compute_cell_to_stl_faces(m::CartesianMesh{D},stl::SurfaceMesh{D}) where D
  cells = Int[]
  faces = Int[]
  cell_cache = Int[]
  for stl_face in 1:num_faces(stl)
    bb = BoundingBox(stl,stl_face)
    for k in cells_around!(cell_cache,m,bb)
      cell = get_cell(m,k)
      if have_intersection(cell,stl,stl_face)
        push!(cells,k)
        push!(faces,stl_face)
      end
    end
  end
  Table(faces,cells,num_cells(m))
end
