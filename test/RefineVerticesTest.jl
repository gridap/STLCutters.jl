module RefineVerticesTest

using WriteVTK
using LinearAlgebra
using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Arrays


using Gridap.Helpers: tfill

## Include in Gridap
function Base.setindex(p::Point,v,idx::Integer)
  data = Base.setindex(p.data,v,idx)
  Point(data)
end

Gridap.num_dims(::Type{<:Point{D}}) where D = D

Gridap.num_dims(::T) where T<:Point = num_dims(T)

Base.prod(p::Point) = prod(p.data)

## Other functions

const d_to_node_to_vef = [
# 1D
  [ 1 3 2],
# 2D
  [ 1 5 2 7 9 8 3 6 4 ] ]

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
  vef = d_to_node_to_vef[D][node]
  compute_vertex_coordinates(cell,p,vef,point)
end

function compute_new_vertices(cell,point::Point{D}) where D
  pmin = first(cell)
  pmax = last(cell)
  partition = tfill(2,Val{D}())
  desc = CartesianDescriptor( pmin, pmax, partition )
  grid = CartesianGrid( desc ) 
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

function move_vertex_to_cell_boundary(cell,p::Polytope,point::Point)
  tol = 1e-6
  for iface in 1:num_faces(p)
    vertex = compute_vertex_coordinates(cell,p,iface,point)
    dist = distance(point,vertex)
    if dist < tol
      point = vertex
    end
  end
  point
end

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

function delete_zero_cells(cells)
  _cells = eltype(cells)[]
  for cell in cells
    if measure(cell) > 0
      push!(_cells,cell)
    end
  end
  _cells
end

function vertex_refinement(cell,point::Point{D}) where D
  p = Polytope(tfill(HEX_AXIS,Val{D}()))
  point = move_vertex_to_cell_boundary(cell,p,point)
  cells = compute_new_vertices(cell,point)
  cells = delete_zero_cells(cells)
  cells
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

function insert_vertices!(T,V,Tnew,STL_vertices)
  for (k,Vk) in zip(T,V)
    if length(Vk) > 0
      v = first(Vk)
      Tk = vertex_refinement(k,STL_vertices[v])
      VTk = distribute_vertices(Tk,Vk[2:end],STL_vertices)
      insert_vertices!(Tk,VTk,Tnew,STL_vertices)
    else
      push!(Tnew,k)
    end
  end
end


# Driver

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

insert_vertices!(T,V,Tnew,STL_vertices)

T = Tnew


writevtk(T,"Tree")





end # module
