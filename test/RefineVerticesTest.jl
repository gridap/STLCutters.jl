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

Gridap.num_dims(::Type{<:Point{D}}) where D = D

Gridap.num_dims(::T) where T<:Point = num_dims(T)

Base.one(::Type{Point{D,T}}) where {D,T} = Point( tfill(one(T),Val{D}()) )

Base.one(::T) where T<:Point = one(T)

Base.prod(p::Point) = prod(p.data)

## Constants


# Get the info from the Polytope
const d_to_node_to_vef = [
# 1D
  [ 1 3 2],
# 2D
  [ 1 5 2 7 9 8 3 6 4 ],
# 3D
  [ 1 9 2 13 21 14 3 10 4 17 23 18 25 27 26 19 24 20 5 11 6 15 22 16 7 12 8 ] ]

# Reuse writeVTK from Gridap

function _vtkpoints(cells)
  D = num_dims(eltype(eltype(cells)))
  vertices = eltype(eltype(cells))[]
  for cell in cells, vertex in cell
    push!(vertices,vertex)
  end
  reshape(reinterpret(Float64,vertices),(D,length(vertices)))
end

function _vtkcells(cells)
  d_to_vtk_type_id = Dict(0=>1,1=>3,2=>9,3=>12)
  vtk_ncube_lvertices = [1,2,4,3,5,6,8,7]
  D = num_dims(eltype(eltype(cells)))

  vtk_cells = MeshCell{VTKCellType,Vector{Int}}[]
  vtk_type = VTKCellType(d_to_vtk_type_id[D])
  offset = 0
  for cell in cells
    vertices = zeros(Int,length(cell))
    vertices = [ vtk_ncube_lvertices[i]+offset for i in 1:length(cell) ]
    offset = maximum(vertices)
    push!( vtk_cells, MeshCell(vtk_type,vertices) )
  end
  vtk_cells
end

function writevtk(cells,base_name)
  vtkpoints = _vtkpoints(cells)
  vtkcells = _vtkcells(cells)
  vtkfile = vtk_grid(base_name,vtkpoints,vtkcells)
  vtk_save(vtkfile)
end

# Helpers

distance(a::Point,b::Point) =  norm( a - b )

function measure(cell)
  pmin = first(cell)
  pmax = last(cell)
  diag = pmax - pmin
  Base.prod(diag)
end

function have_intersection(cell,p::Point)
  pmin = first(cell)
  pmax = last(cell)
  all( pmin.data .< p.data ) || return false
  all( pmax.data .> p.data ) || return false
  true
end

function distance_to_boundary(cell,p::Point)
  @assert have_intersection(cell,p)
  pmin = first(cell)
  pmax = last(cell)
  min( minimum(p-pmin), minimum(pmax-p) )
end

function compute_vertex_coordinates(cell,p::Polytope{D},iface::Integer,point::Point{D}) where D
  nface = p.dface.nfaces[iface]
  dim = p.dface.dims[iface]
  anchor = cell[ get_faces(p)[iface][1] ]
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

function compute_vertex_coordinates(cell,grid::Grid{D},subcell::Integer,lvertex::Integer,point::Point{D}) where D
  reffe = get_reffes( grid )[1]
  p = get_polytope(reffe)
  node = grid.cell_nodes[subcell][lvertex]
  # compute_vertex_coordinates(cell,grid,node,point)
  vef = d_to_node_to_vef[D][node]
  compute_vertex_coordinates(cell,p,vef,point)
end

# Main algorithm's functions

# Function not used
function compute_vertex_coordinates(cell,grid::Grid,v::Integer,point::Point{D}) where D
  function anchor(a)
    Int.(floor.(a.data))
  end
  function extrusion(a)
    Int.( a.data .== 0.5 )
  end
  function vertex_id(a)
    id = 0
    for d in 1:D
      id |= a[d] << (d-1)
    end
    id+1
  end

  rvertex = grid.node_coords[v]
  anc = anchor(rvertex)
  vertex = cell[ vertex_id(anc) ]
  ext = extrusion(rvertex)
  for d in 1:D
    if ext[d] == HEX_AXIS
      v = point[d]
      vertex = Base.setindex(vertex,v,d)
    end
  end
  vertex
end
## 

function split_cartesian_grid(::Val{D}) where D
  pmin = zero( Point{D,Int} )
  pmax = one( Point{D,Int} )
  partition = tfill(2,Val{D}())
  desc = CartesianDescriptor( pmin, pmax, partition )
  CartesianGrid( desc ) 
end

function compute_new_vertices(cell,point::Point{D}) where D
  grid = split_cartesian_grid(Val{D}())
  cells = Vector{typeof(point)}[]
  for subcell in 1:num_cells(grid)
    vertices = typeof(point)[]
    for lvertex in 1:length(grid.cell_nodes[subcell])
      vertex = compute_vertex_coordinates(cell,grid,subcell,lvertex,point)
      push!(vertices,vertex)
    end
    push!(cells,vertices)
  end
  cells
end

function move_vertex_to_cell_boundary(cell,point::Point{D}) where D
  pmin = first(cell)
  pmax = last(cell)
  for d in 1:D
    if point[d] - pmin[d] < TOL 
      point = Base.setindex(point,pmin[d],d)
    elseif pmax[d] - point[d] < TOL pmax[d]
      point = Base.setindex(point,pmax[d],d)
    end
  end
  point
end

function farthest_vertex_from_boundary(cell,vertices,STL_vertices)
  iv = UNSET
  max_dist = TOL
  for (i,v) in enumerate(vertices)
    p = STL_vertices[v]
    dist = distance_to_boundary(cell,p)
    if dist > max_dist
      max_dist = dist
      iv = i
    end
  end
  iv
end

function vertex_refinement(cell,point::Point{D}) where D
  compute_new_vertices(cell,point)
end

function distribute_vertices(cells,vertices,STL_vertices)
  cell_to_vertices = Vector{Int}[]
  for (i,cell) in enumerate(cells)
    push!(cell_to_vertices,[])
    for vertex in vertices
      point = STL_vertices[vertex]
      if have_intersection(cell,point)
        push!(cell_to_vertices[i],vertex)
      end
    end
  end
  cell_to_vertices
end

const TOL = 1e6*eps()

function insert_vertices!(T,V,Tnew,STL_vertices,Tnew_to_v,v_in)
  for (k,Vk) in zip(T,V)

    iv = farthest_vertex_from_boundary(k,Vk,STL_vertices)

    if iv != UNSET
      v = Vk[iv]
      deleteat!(Vk,iv)
      v_in_k = [ v_in ; [v] ]
      Tk = vertex_refinement(k,STL_vertices[v])
      VTk = distribute_vertices(Tk,Vk,STL_vertices)
      insert_vertices!(Tk,VTk,Tnew,STL_vertices,Tnew_to_v,v_in_k)
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

K = [ 
  Point(0.0,0.0),
  Point(1.0,0.0),
  Point(0.0,1.0),
  Point(1.0,1.0) ]

T = [K]

V = distribute_vertices(T,1:length(STL_vertices),STL_vertices)

Tnew = eltype(T)[]

Tnew_to_v = Vector{Int}[]

v_in = Int[]

insert_vertices!(T,V,Tnew,STL_vertices,Tnew_to_v,v_in)

T = Tnew
T_to_v = Tnew_to_v

display(T_to_v)

@test length(T) == length(T_to_v) == 13
@test T_to_v[1] == T_to_v[4] == [2,4,1]
@test T_to_v[5] == T_to_v[8] == [2,4,3]

writevtk(T,"Tree")

## 2D: Flat cells

STL_vertices = [ 
  Point(0.1,0.2),
  Point(0.5,0.5),
  Point(0.4,0.1),
  Point(0.5-1e-10,0.3) ]


K = [ 
  Point(0.0,0.0),
  Point(1.0,0.0),
  Point(0.0,1.0),
  Point(1.0,1.0) ]

T = [K]

V = distribute_vertices(T,1:length(STL_vertices),STL_vertices)

Tnew = eltype(T)[]

Tnew_to_v = Vector{Int}[]

v_in = Int[]

insert_vertices!(T,V,Tnew,STL_vertices,Tnew_to_v,v_in)

T = Tnew
T_to_v = Tnew_to_v

display(T_to_v)

#writevtk(T,"Tree")

## 3D

STL_vertices = [ 
  Point(0.1,0.2,0.3),
  Point(0.5,0.5,0.5),
  Point(0.4,0.1,0.2),
  Point(0.3,0.7,0.4) ]


K = [ 
  Point(0.0,0.0,0.0),
  Point(1.0,0.0,0.0),
  Point(0.0,1.0,0.0),
  Point(1.0,1.0,0.0),
  Point(0.0,0.0,1.0),
  Point(1.0,0.0,1.0),
  Point(0.0,1.0,1.0),
  Point(1.0,1.0,1.0) ]

T = [K]

V = distribute_vertices(T,1:length(STL_vertices),STL_vertices)

Tnew = eltype(T)[]

Tnew_to_v = Vector{Int}[]

v_in = Int[]

insert_vertices!(T,V,Tnew,STL_vertices,Tnew_to_v,v_in)

T = Tnew
T_to_v = Tnew_to_v

display(T_to_v)


writevtk(T,"3DTree")

end # module
