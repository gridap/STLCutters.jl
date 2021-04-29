
struct STLGeometry <: CSG.Geometry
  tree::Leaf{Tuple{T,String,Nothing}} where T<:DiscreteModel 
end

function STLGeometry(stl::DiscreteModel;name="stl")
  tree = Leaf( ( stl, name, nothing ) ) 
  STLGeometry( tree )
end

function STLGeometry(filename::String;name="stl")
  @notimplemented
  STLGeometry(stl,name=name)
end

get_tree(geo::STLGeometry) = geo.tree

get_stl(geo::STLGeometry) = geo.tree.data[1]

function compatible_geometries(a::STLGeometry,b::STLGeometry)
  a,b
end

function similar_geometry(a::STLGeometry,tree::Leaf)
  STLGeometry(tree)
end

struct STLCutter <: Interfaces.Cutter end

function cut(cutter::STLCutter,background::DiscreteModel,geom::STLGeometry)
  data,bgf_to_ioc = _cut_sm(background,geom)
  EmbeddedDiscretization(background, data..., geom), bgf_to_ioc
end

function cut(background::DiscreteModel,geom::STLGeometry)
  cutter = STLCutter()
  cut(cutter,background,geom)
end


function _cut_sm(model::DiscreteModel,geom::STLGeometry)
  grid = get_grid(model)
  stl = get_stl(geom)
  out = compute_submesh(model,get_stl(geom))
  T,X,F,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf,bgcell_to_ioc,bgfacet_to_ioc = out
  bgcell_to_ioc = [ map( _convert_io, bgcell_to_ioc ) ]
  cell_to_points = Table( map(i->Int32.(i),T) )
  point_to_coords = X
  point_to_rcoords = send_to_ref_space(grid,k_to_bgcell,T,X)
  cell_to_io = [ map( _convert_io, k_to_io ) ]
  cell_to_bgcell = Int32.(k_to_bgcell)
  data = cell_to_points,cell_to_bgcell,point_to_coords,point_to_rcoords
  subcells = SubCellData(data...)

  f_to_points = Table( map(i->Int32.(i),F) )
  point_to_coords = Xf
  point_to_rcoords = send_to_ref_space(grid,f_to_bgcell,F,Xf)
  f_to_normal = [ normal(get_cell(get_stl(geom),facet)) for facet in f_to_stlf ]
  f_to_bgcell = Int32.(f_to_bgcell)
  data = f_to_points,f_to_normal,f_to_bgcell,point_to_coords,point_to_rcoords
  subfacets = SubFacetData(data...) 
  f_to_io = fill(Int8(INTERFACE),length(F))
  f_to_io = [ f_to_io ]

  replace!( bgfacet_to_ioc, FACE_IN => IN, FACE_OUT => OUT, FACE_CUT => CUT )

  oid_to_ls = Dict{UInt,Int}( objectid( get_stl(geom) ) => 1  )
  (bgcell_to_ioc,subcells,cell_to_io,subfacets,f_to_io,oid_to_ls),bgfacet_to_ioc
end

#function _cut_sm(model::DiscreteModel,geom::STLGeometry)
#  out = refine_grid(get_grid(model),get_stl(geom))
#  T,X,reffes,cell_types,cell_to_io,cell_to_bgcell,bgcell_to_ioc = out
#  bgcell_to_ioc = [ map( _convert_io, bgcell_to_ioc ) ]
#  subgrid = UnstructuredGrid(X,Table(T),reffes,cell_types)
#  out = get_subcell_coords(subgrid,get_grid(model),cell_to_bgcell,cell_to_io)
#  cell_to_points,point_to_coords,point_to_rcoords,cell_to_bgcell,cell_to_io = out
#  data = Table(cell_to_points),Int32.(cell_to_bgcell),point_to_coords,point_to_rcoords
#  subcells = SubTriangulation(data...)
#  cell_to_io = [ map( _convert_io, cell_to_io ) ]
#
#  _stl,f_to_bgcell,f_to_f = refine_surface(get_grid(model),get_stl(geom))
#  f_to_io = fill(Int8(INTERFACE),num_cells(_stl))
#  out = get_subcell_coords(get_grid(_stl),get_grid(model),f_to_bgcell,f_to_io)
#  f_to_points,point_to_coords,point_to_rcoords,f_to_bgcell,f_to_io = out 
#  f_to_normal = [ normal(get_cell(get_stl(geom),facet)) for facet in f_to_f ]
#  data = Table(f_to_points),f_to_normal,Int32.(f_to_bgcell),point_to_coords,point_to_rcoords
#  subfacets = FacetSubTriangulation(data...) 
#  f_to_io = [ f_to_io ]
#
#  oid_to_ls = Dict{UInt,Int}( objectid( get_stl(geom) ) => 1  )
#  bgcell_to_ioc,subcells,cell_to_io,subfacets,f_to_io,oid_to_ls
#end

function send_to_ref_space(grid::Grid,cell::Integer,point::Point)
  cell_coordinates = get_cell_coordinates(grid)[cell]
  pmin,pmax = cell_coordinates[1],cell_coordinates[end]
  Point(Tuple(point-pmin)./Tuple(pmax-pmin)) 
end

function send_to_ref_space(grid::Grid,cell::Integer,points::Vector)
  cell_coordinates = get_cell_coordinates(grid)[cell]
  pmin,pmax = cell_coordinates[1],cell_coordinates[end]
  [ Point(Tuple(point-pmin)./Tuple(pmax-pmin)) for point in points ]
end

function send_to_ref_space(grid::Grid,cell_to_bgcell::Vector,cell_to_points::Vector,points::Vector)
  rpoints = zero(points)
  for (i,k) in enumerate(cell_to_points)
    for p in k
      rpoints[p] = send_to_ref_space(grid,cell_to_bgcell[i],points[p])
    end
  end
  rpoints
end

function get_subcell_coords(subgrid,grid,cell_to_bgcell,cell_to_io)
  cell_to_points = Vector{Int}[]
  point_to_coords = empty(get_node_coordinates(subgrid))
  point_to_rcoords = empty(get_node_coordinates(subgrid))
  new_cell_to_bgcell = Int32[]
  new_cell_to_io = Int8[]
  for cell in 1:num_cells(subgrid)
    p = get_polytope( get_cell_reffes(subgrid)[cell] )
    bgcell = cell_to_bgcell[cell]
    coords = get_cell_coordinates(subgrid)[cell]
    rcoords = send_to_ref_space(grid,bgcell,coords)
    points = length(point_to_coords) .+ (1:length(coords))
    _cell_to_points = map( i -> points[i], simplexify(p)[1] )
    append!(cell_to_points,_cell_to_points)
    append!(point_to_coords,coords)
    append!(point_to_rcoords,rcoords)
    append!(new_cell_to_bgcell,fill(bgcell,length(_cell_to_points)))
    append!(new_cell_to_io,fill(cell_to_io[cell],length(_cell_to_points)))
  end
  cell_to_points,
  point_to_coords,
  point_to_rcoords,
  new_cell_to_bgcell,
  new_cell_to_io
end

function _convert_io(ioc::Integer)
  if ioc == FACE_IN
    Int8(IN)
  elseif ioc == FACE_OUT
    Int8(OUT)
  elseif ioc == FACE_CUT
    Int8(CUT)
  else
    Int8(OUT)
    #@unreachable
  end
end

## Remove the next lines when upgrading GridapEmbedded to v0.7

using GridapEmbedded.CSG
using GridapEmbedded.AgFEM
using GridapEmbedded.AgFEM: _aggregate_by_threshold
import GridapEmbedded.AgFEM: aggregate

function aggregate(
  strategy,
  cut::EmbeddedDiscretization,
  facet_to_inoutcut::AbstractVector)

  aggregate(strategy,cut,cut.geo,IN,facet_to_inoutcut)
end

function aggregate(
  strategy::AggregateCutCellsByThreshold,
  cut::EmbeddedDiscretization,
  geo::CSG.Geometry,
  in_or_out,
  facet_to_inoutcut::AbstractVector)

  _aggregate_by_threshold(strategy.threshold,cut,geo,in_or_out,facet_to_inoutcut)
end

