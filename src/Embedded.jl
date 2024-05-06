
const UNDEF = -10

"""
    struct STLGeometry

Object that stores the SLT geometry.

# Example
```@repl
julia> filename = download_thingi10k(293137)
julia> geo = STLGeometry(filename)
julia> writevtk(geo,"geo")
```
"""
struct STLGeometry <: CSG.Geometry
  tree::Leaf{Tuple{T,String,Nothing}} where T<:DiscreteModel
end

function STLGeometry(stl::DiscreteModel;name="stl")
  tree = Leaf( ( stl, name, nothing ) )
  STLGeometry( tree )
end

function STLGeometry(filename::String;name="stl")
  stlmodel = _to_stl_model(filename)
  STLGeometry(stlmodel,name=name)
end

get_tree(geo::STLGeometry) = geo.tree

get_stl(geo::STLGeometry) = geo.tree.data[1]

function compatible_geometries(a::STLGeometry,b::STLGeometry)
  a,b
end

function similar_geometry(a::STLGeometry,tree::Leaf)
  STLGeometry(tree)
end

function closest_point(x::AbstractVector,geo::STLGeometry,args...)
  closest_point(x,get_stl(geo),args...)
end

get_bounding_box(a::STLGeometry) = get_bounding_box( get_stl(a) )

struct STLEmbeddedDiscretization <: Interfaces.AbstractEmbeddedDiscretization
  cut::EmbeddedDiscretization
  cutfacets::EmbeddedFacetDiscretization
end

function get_background_model(cut::STLEmbeddedDiscretization)
  get_background_model(cut.cut)
end

function get_geometry(cut::STLEmbeddedDiscretization)
  get_geometry(cut.cut)
end

function compute_bgcell_to_inoutcut(
  cut::STLEmbeddedDiscretization,
  geo::STLGeometry)

  compute_bgcell_to_inoutcut(cut.cut,geo)
end

function compute_bgcell_to_inoutcut(cut::STLEmbeddedDiscretization)

  geo = get_geometry(cut)
  compute_bgcell_to_inoutcut(cut,geo)
end

function compute_bgfacet_to_inoutcut(
  cut::STLEmbeddedDiscretization,
  geo::STLGeometry)

  compute_bgfacet_to_inoutcut(cut.cutfacets,geo)
end

function compute_bgfacet_to_inoutcut(cut::STLEmbeddedDiscretization)
  geo = get_geometry(cut)
  compute_bgfacet_to_inoutcut(cut,geo)
end

function Triangulation(cut::STLEmbeddedDiscretization,args...)
  Triangulation(cut.cut,args...)
end

function compute_subcell_to_inout(
  cut::STLEmbeddedDiscretization,
  geo::STLGeometry)

  compute_subcell_to_inout(cut.cut,geo)
end

function EmbeddedBoundary(cut::STLEmbeddedDiscretization,args...)
  EmbeddedBoundary(cut.cut,args...)
end

function GhostSkeleton(cut::STLEmbeddedDiscretization,args...)
  GhostSkeleton(cut.cut,args...)
end

struct STLCutter <: Interfaces.Cutter
  options::Dict{Symbol,Any}
  STLCutter(;options...) = new(options)
end

"""
    cut([cutter::STLCutter,]model::DiscreteModel,geo::STLGeometry)


Main interface of GridapEmbedded. It intersect each cell in `model` with the
geometry `geo`. It returns an embedded discretization.

# Usage
The basic usage is to call `cut` with a model and a geometry.
```@example
cutgeo = cut(model,geo)
```
We can also set custom keyword arguments, see [`subtriangulate`](@ref).
```@example
cutter = STLCutter(;tolfactor=1e3)
cutgeo = cut(cutter,model,geo)
```
Then we can extract the domain and the boundary of the embedded discretization.
```@example
Ω = Triangulation(cutgeo)
Γ = EmbeddedBoundary(cutgeo)
```
"""
function cut(cutter::STLCutter,background::DiscreteModel,geom::STLGeometry)
  data,facet_data = _cut_stl(background,geom;cutter.options...)
  STLEmbeddedDiscretization(
    EmbeddedDiscretization(background, data..., geom),
    EmbeddedFacetDiscretization(background, facet_data..., geom) )
end

function cut(background::DiscreteModel,geom::STLGeometry)
  cutter = STLCutter()
  cut(cutter,background,geom)
end

function cut_facets(cutter::STLCutter,background::DiscreteModel,geom::STLGeometry)
  cutgeo = cut(background,geom;cutter.options...)
  cut_facets(cutgeo)
end

function cut_facets(background::DiscreteModel,geom::STLGeometry)
  cutter = STLCutter()
  cut_facets(cutter,background,geom)
end

function cut_facets(cut::STLEmbeddedDiscretization,args...)
  cut.cutfacets
end

function aggregate(strategy,cut::STLEmbeddedDiscretization)
  aggregate(strategy,cut,cut.cut.geo)
end

function aggregate(strategy,cut::STLEmbeddedDiscretization,geo)
  aggregate(strategy,cut,geo,IN)
end

function aggregate(strategy,cut::STLEmbeddedDiscretization,name::String,in_or_out)
  geo = get_geometry(cut.geo,name)
  aggregate(strategy,cut,geo,in_or_out)
end

function aggregate(
  strategy,
  cut::STLEmbeddedDiscretization,
  geo::STLGeometry,
  in_or_out)

  facet_to_inoutcut = compute_bgfacet_to_inoutcut(cut,geo)
  aggregate(strategy,cut.cut,geo,in_or_out,facet_to_inoutcut)
end

function _cut_stl(model::DiscreteModel,geom::STLGeometry;kwargs...)
  subcell_grid, subface_grid, bsubface_grid, labels = subtriangulate(model,geom;kwargs...)

  D = num_dims(model)
  inout_dict = inout_dictionary(;kwargs...)
  cell_to_io = [ replace( labels.cell_to_io, inout_dict... ) ]
  cell_to_bgcell = labels.cell_to_bgcell
  cell_to_points = get_cell_node_ids( subcell_grid )
  point_to_coords = get_node_coordinates( subcell_grid )
  point_to_rcoords = send_to_ref_space(model,cell_to_bgcell,subcell_grid)
  data = cell_to_points,cell_to_bgcell,point_to_coords,point_to_rcoords
  subcells = SubCellData(data...)

  face_to_bgcell = labels.face_to_bgcell
  face_to_points = get_cell_node_ids( subface_grid )
  point_to_coords = get_node_coordinates( subface_grid )
  point_to_rcoords = send_to_ref_space(model,face_to_bgcell,subface_grid)
  face_to_normal = _normals(geom,labels.face_to_stlface)
  data = face_to_points,face_to_normal,face_to_bgcell,point_to_coords,point_to_rcoords
  subfacets = SubFacetData(data...)

  bface_to_bgcell = labels.bface_to_bgcell
  bface_to_lbgface = labels.bface_to_lbgface
  bface_to_bgface = boundary_facet_to_background_facet(model,bface_to_bgcell,bface_to_lbgface)
  newbfaces = unique_boundary_facets(bface_to_bgcell,bface_to_bgface)
  bsubface_grid = isolate_cell_coordinates(bsubface_grid,newbfaces)
  bface_to_bgface = bface_to_bgface[newbfaces]
  bface_to_io = labels.bface_to_io[newbfaces]
  bface_to_points = get_cell_node_ids( bsubface_grid )
  point_to_coords = get_node_coordinates( bsubface_grid )
  point_to_rcoords = send_to_ref_space(Val{D-1}(),model,bface_to_bgface,bsubface_grid)
  bface_to_io = [ replace( bface_to_io, inout_dict... ) ]
  data = bface_to_points,bface_to_bgface,point_to_coords,point_to_rcoords
  bsubfacets = SubCellData(data...)

  face_to_io = [ fill(Int8(INTERFACE),num_cells(subface_grid)) ]

  bgface_to_ioc = [ replace( labels.bgface_to_ioc, inout_dict... ) ]
  bgcell_to_ioc = [ replace( labels.bgcell_to_ioc, inout_dict... ) ]

  oid_to_ls = Dict{UInt,Int}( objectid( get_stl(geom) ) => 1  )
  (bgcell_to_ioc,subcells,cell_to_io,subfacets,face_to_io,oid_to_ls),
  (bgface_to_ioc,bsubfacets,bface_to_io,oid_to_ls)
end

function inout_dictionary(;all_defined=true,kwargs...)
  undefined = all_defined ? OUT : UNDEF
  Dict{Int8,Int8}(
    FACE_IN => IN,
    FACE_OUT => OUT,
    FACE_CUT => CUT,
    UNSET => undefined )
end

function isolate_cell_coordinates(grid::Grid,cells=1:num_cells(grid))
  cell_coords = Table(get_cell_coordinates(grid)[cells])
  coords = cell_coords.data
  data = 1:length(cell_coords.data)
  ptrs = cell_coords.ptrs
  cell_nodes = Table(data,ptrs)
  reffes = get_reffes(grid)
  cell_types = get_cell_type(grid)[cells]
  UnstructuredGrid(coords,cell_nodes,reffes,cell_types)
end

function boundary_facet_to_background_facet(model::DiscreteModel,f_to_bgc,f_to_lbgf)
  topo = get_grid_topology(model)
  D = num_dims(model)
  c_to_f = get_faces(topo,D,D-1)
  cache = array_cache(c_to_f)
  f_to_bgf = fill(Int32(UNSET),length(f_to_bgc))
  for (i,(c,lf)) in enumerate(zip(f_to_bgc,f_to_lbgf))
    f_to_bgf[i] = getindex!(cache,c_to_f,c)[lf]
  end
  f_to_bgf
end

function unique_boundary_facets(bface_to_bgcell,bface_to_bgface)
  num_facets = maximum(bface_to_bgface,init=0)
  bgface_to_bgcell_owner = fill(UNSET,num_facets)
  for (bgcell,bgface) in zip(bface_to_bgcell,bface_to_bgface)
    if bgface_to_bgcell_owner[bgface] == UNSET
      bgface_to_bgcell_owner[bgface] = bgcell
    end
  end
  bface_to_bgcell_owner = lazy_map(Reindex(bgface_to_bgcell_owner),bface_to_bgface)
  bface_mask = lazy_map(==,bface_to_bgcell,bface_to_bgcell_owner)
  findall(bface_mask)
end

function send_to_ref_space(
  model::DiscreteModel,
  cell_to_bgcell::Vector,
  subgrid::Grid)

  send_to_ref_space(get_grid(model),cell_to_bgcell,subgrid)
end

function send_to_ref_space(
  ::Val{D},
  model::DiscreteModel,
  cell_to_bgcell::Vector,
  subgrid::Grid) where D

  grid = Grid(ReferenceFE{D},model)
  send_to_ref_space(grid,cell_to_bgcell,subgrid)
end

function send_to_ref_space(grid::Grid,cell_to_bgcell::Vector,subgrid::Grid)
  bgcell_map = _get_cell_affine_map(grid)
  bgcell_invmap = lazy_map(inverse_map,bgcell_map)
  cell_invmap = lazy_map(Reindex(bgcell_invmap),cell_to_bgcell)
  cell_nodes = get_cell_node_ids(subgrid)
  node_coords = get_node_coordinates(subgrid)
  node_cell = _first_inverse_index_map(cell_nodes,num_nodes(subgrid))
  node_invmap = lazy_map(Reindex(cell_invmap),node_cell)
  node_rcoords = lazy_map(evaluate,node_invmap,node_coords)
  _collect(node_rcoords)
end

function _first_inverse_index_map(a_to_b,nb)
  na = length(a_to_b)
  T = eltype(eltype(a_to_b))
  b_to_a = ones(T,nb)
  c = array_cache(a_to_b)
  for a in 1:na
    bs = getindex!(c,a_to_b,a)
    for b in bs
      b_to_a[b] = a
    end
  end
  b_to_a
end

function _normals(geom::STLGeometry,face_to_stlface)
  stl = get_stl(geom)
  cache = get_cell_cache(stl)
  N = typeof(normal(get_cell!(cache,stl,1)))
  normals = zeros(N,length(face_to_stlface))
  for (i,stlf) in enumerate( face_to_stlface )
    normals[i] = normal(get_cell!(cache,stl,stlf))
  end
  normals
end

_to_stl_model(geo::STLGeometry) = get_stl(geo)

surface(geo::STLGeometry) = surface(get_grid(get_stl(geo)))

function writevtk(geo::STLGeometry,args...;kwargs...)
  grid = get_grid(get_stl(geo))
  writevtk(grid,args...;kwargs...)
end

function check_requisites(geo::STLGeometry,bgmodel::DiscreteModel;kwargs...)
  stl = get_stl(geo)
  check_requisites(stl,bgmodel)
end

function _collect(a::LazyArray)
  c = array_cache(a)
  T = eltype(a)
  b = zeros(T,size(a))
  for i in eachindex(a)
    b[i] = getindex!(c,a,i)
  end
  b
end

function _get_cell_affine_map(grid::Grid)
  @assert all(==(1),map(get_order,get_reffes(grid)))
  D = num_dims(grid)
  x0 = Point(tfill(0.0,Val{D}()))
  ncells = num_cells(grid)
  cell_map = get_cell_map(grid)
  cell_map_gradient = lazy_map(∇,cell_map)
  origins = lazy_map(evaluate,cell_map,Fill(x0,ncells))
  gradiens = lazy_map(evaluate,cell_map_gradient,Fill(x0,ncells))
  lazy_map(affine_map,gradiens,origins)
end
