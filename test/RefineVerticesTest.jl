module RefineVerticesTest

using Test
using WriteVTK
using LinearAlgebra
using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Arrays


using Gridap.Helpers: tfill

## Missing from Gridap

function Base.setindex(p::Point,v,idx::Integer)
  data = Base.setindex(p.data,v,idx)
  Point(data)
end

Base.prod(p::Point) = prod(p.data)

# Helpers

distance(a::Point,b::Point) =  norm( a - b )

function get_bounding_box(
  cell_nodes::Vector{<:Integer},
  node_to_coordinates::Vector{<:Point})

  pmin = node_to_coordinates[ first(cell_nodes) ]
  pmax = node_to_coordinates[ last(cell_nodes) ]
  pmin,pmax
end

function have_intersection(cell_nodes,node_to_coordinates,p::Point)
  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates)
  all( pmin.data .< p.data ) || return false
  all( pmax.data .> p.data ) || return false
  true
end

function distance_to_boundary(cell_nodes,node_to_coordinates,p::Point)
  @assert have_intersection(cell_nodes,node_to_coordinates,p)
  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates)
  min( minimum(p-pmin), minimum(pmax-p) )
end

function compute_vertex_coordinates(
  cell_nodes,
  node_to_coordinates,
  p::Polytope{D},
  iface::Integer,
  point::Point{D}) where D

  nface = p.dface.nfaces[iface]
  dim = p.dface.dims[iface]
  node = cell_nodes[ get_faces(p)[iface][1] ]
  anchor = node_to_coordinates[ node ]
  extrusion = nface.extrusion
  vertex = anchor
  for d in 1:D
    if extrusion[d] == HEX_AXIS
      v = point[d]
      vertex = Base.setindex(vertex,v,d)
    end
  end
  vertex
end

function compute_linear_grid_and_facemap(reffe::LagrangianRefFE)
  grid = compute_linear_grid(reffe)
  desc = get_cartesian_descriptor(grid)
  model = CartesianDiscreteModel(desc)
  labels = get_face_labeling(model)
  grid_face_to_reffe_face = get_face_entity(labels)
  grid,grid_face_to_reffe_face
end


function compute_new_cells(
  cell_nodes::Vector{<:Integer},
  node_to_coordinates::Vector{<:Point},
  reffe::ReferenceFE)

  grid, gface_to_rface = compute_linear_grid_and_facemap(reffe)
  num_nodes_per_cell = length(cell_nodes)
  num_nodes = length(node_to_coordinates)
  new_cells = Vector{Int}[]
  for lcell in 1:num_cells(grid)
    new_cell = fill(UNSET,num_nodes_per_cell)
    for lnode in 1:num_nodes_per_cell
      node = get_cell_nodes(grid)[lcell][lnode] 
      node = gface_to_rface[node]
      if node ≤ num_nodes_per_cell
        n = cell_nodes[node]
      else
        n = node - num_nodes_per_cell + num_nodes
      end
      new_cell[lnode] = n
    end
    push!(new_cells,new_cell)
  end
  new_cells
end

function compute_new_vertices(
  cell_nodes::Vector{<:Integer},
  node_to_coordinates::Vector{<:Point},
  reffe::ReferenceFE,
  point::Point{D}) where D

  p = get_polytope(reffe)
  num_nodes_per_cell = length(cell_nodes)
  new_node_to_coordinates = eltype(point)[]
  for node in num_nodes_per_cell+1:num_nodes(reffe)
    vertex = compute_vertex_coordinates(cell_nodes,node_to_coordinates,p,node,point)
    push!(node_to_coordinates,vertex)
  end
  new_node_to_coordinates
end

function vertex_refinement(
  cell_nodes,
  node_to_coordinates::Vector{<:Point},
  point::Point{D}) where D

  p = Polytope(tfill(HEX_AXIS,Val{D}()))
  reffe = LagrangianRefFE(Float64,p,2)
  new_cells = compute_new_cells(cell_nodes,node_to_coordinates,reffe)
  new_vertices = compute_new_vertices(cell_nodes,node_to_coordinates,reffe,point)
  new_cells, new_vertices
end

function move_vertex_to_cell_boundary(
  cell_nodes::Vector{<:Integer},
  node_to_coordinates::Vector{<:Point},
  point::Point{D}) where D

  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates)
  for d in 1:D
    if point[d] - pmin[d] < TOL 
      point = Base.setindex(point,pmin[d],d)
    elseif pmax[d] - point[d] < TOL pmax[d]
      point = Base.setindex(point,pmax[d],d)
    end
  end
  point
end

function farthest_vertex_from_boundary(
  cell_nodes::Vector{<:Integer},
  node_to_coordinates::Vector{<:Point},
  vertices::Vector,
  STL_vertices)

  iv = UNSET
  max_dist = 0.0
  for (i,v) in enumerate(vertices)
    p = STL_vertices[v]
    dist = distance_to_boundary(cell_nodes,node_to_coordinates,p)
    if dist > max_dist
      max_dist = dist
      iv = i
    end
  end
  iv
end

function distribute_vertices(
  cell_to_nodes,
  node_to_coordinates::Vector{<:Point},
  vertices::AbstractVector,
  STL_vertices)

  cell_to_vertices = Vector{Int}[]
  for (i,cell) in enumerate(cell_to_nodes)
    push!(cell_to_vertices,[])
    for vertex in vertices
      point = STL_vertices[vertex]
      if have_intersection(cell,node_to_coordinates,point)
        push!(cell_to_vertices[i],vertex)
      end
    end
  end
  cell_to_vertices
end

const TOL = 1e6*eps()

function insert_vertices!(T,X,V,Tnew,STL_vertices,Tnew_to_v,v_in)
  for (k,Vk) in zip(T,V)

    iv = UNSET
    while length(Vk) > 0
      iv = farthest_vertex_from_boundary(k,X,Vk,STL_vertices)
      p = STL_vertices[Vk[iv]]
      dist = distance_to_boundary(k,X,p)
      if dist < TOL
        p = move_vertex_to_cell_boundary(k,X,p)
        STL_vertices[Vk[iv]] = p
        deleteat!(Vk,iv)
      else
        break
      end
    end

    if length(Vk) > 0
      v = Vk[iv]
      deleteat!(Vk,iv)
      v_in_k = [ v_in ; [v] ]
      Tk,Xnew = vertex_refinement(k,X,STL_vertices[v])
      append!(X,Xnew)
      VTk = distribute_vertices(Tk,X,Vk,STL_vertices)
      insert_vertices!(Tk,X,VTk,Tnew,STL_vertices,Tnew_to_v,v_in_k)
    else
      push!(Tnew,k)
      push!(Tnew_to_v,v_in)
    end
  end
end

# Driver

## 2D

STL_vertices = [ 
  Point(0.1,0.2),
  Point(0.5,0.5),
  Point(0.4,0.1),
  Point(0.3,0.3) ]

K = [ 1, 2, 3, 4 ]

X = [
  Point(0.0,0.0),
  Point(1.0,0.0),
  Point(0.0,1.0),
  Point(1.0,1.0) ]

T = [K]

V = distribute_vertices(T,X,1:length(STL_vertices),STL_vertices)

Tnew = eltype(T)[]

Tnew_to_v = Vector{Int}[]

v_in = Int[]

insert_vertices!(T,X,V,Tnew,STL_vertices,Tnew_to_v,v_in)

T = Tnew
T_to_v = Tnew_to_v

display(T_to_v)

D = 2
@test length(T) == length(T_to_v) == (2^D-1)*length(STL_vertices)+1
@test T_to_v[1] == T_to_v[4] == [2,4,1]
@test T_to_v[5] == T_to_v[8] == [2,4,3]
@test length(X) == (3^D-2^D)*length(STL_vertices)+2^D


T = Table(T)
reffes = [QUAD4]
cell_types = fill( 1, length(T) )
grid = UnstructuredGrid(X,T,reffes,cell_types)

writevtk(grid,"Tree")

## 2D: Move vertices

STL_vertices = [ 
  Point(0.1,0.2),
  Point(0.5,0.5),
  Point(0.4,0.1),
  Point(0.5-1e-10,0.3) ]

K = [ 1, 2, 3, 4 ]

X = [ 
  Point(0.0,0.0),
  Point(1.0,0.0),
  Point(0.0,1.0),
  Point(1.0,1.0) ]

T = [K]

V = distribute_vertices(T,X,1:length(STL_vertices),STL_vertices)

Tnew = eltype(T)[]

Tnew_to_v = Vector{Int}[]

v_in = Int[]

insert_vertices!(T,X,V,Tnew,STL_vertices,Tnew_to_v,v_in)

T = Tnew
T_to_v = Tnew_to_v

display(T_to_v)

D = 2
@test length(T) == length(T_to_v) == (2^D-1)*(length(STL_vertices)-1)+1
@test T_to_v[2] == T_to_v[5] == [2,1,3]
@test T_to_v[6] == T_to_v[7] == [2,1]
@test length(X) == (3^D-2^D)*(length(STL_vertices)-1)+2^D

T = Table(T)
reffes = [QUAD4]
cell_types = fill( 1, length(T) )
grid = UnstructuredGrid(X,T,reffes,cell_types)

writevtk(grid,"Tree")

## 3D

STL_vertices = [ 
  Point(0.1,0.2,0.3),
  Point(0.5,0.5,0.5),
  Point(0.4,0.1,0.2),
  Point(0.3,0.7,0.4) ]

K = [1,2,3,4,5,6,7,8]


X = [ 
  Point(0.0,0.0,0.0),
  Point(1.0,0.0,0.0),
  Point(0.0,1.0,0.0),
  Point(1.0,1.0,0.0),
  Point(0.0,0.0,1.0),
  Point(1.0,0.0,1.0),
  Point(0.0,1.0,1.0),
  Point(1.0,1.0,1.0) ]

T = [K]

V = distribute_vertices(T,X,1:length(STL_vertices),STL_vertices)

Tnew = eltype(T)[]

Tnew_to_v = Vector{Int}[]

v_in = Int[]

insert_vertices!(T,X,V,Tnew,STL_vertices,Tnew_to_v,v_in)

T = Tnew
T_to_v = Tnew_to_v

display(T_to_v)

D = 3
@test length(T) == length(T_to_v) == (2^D-1)*length(STL_vertices)+1
@test T_to_v[2] == T_to_v[9] == [2,1,3]
@test T_to_v[10] == T_to_v[15] == [2,1]
@test T_to_v[17] == T_to_v[24] == [2,4]
@test length(X) == (3^D-2^D)*length(STL_vertices)+2^D


T = Table(T)
reffes = [HEX8]
cell_types = fill( 1, length(T) )
grid = UnstructuredGrid(X,T,reffes,cell_types)


writevtk(grid,"3DTree")

end # module
