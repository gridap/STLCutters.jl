
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

get_bounding_box(a::STLGeometry) = get_bounding_box( get_stl(a) )

struct STLCutter <: Interfaces.Cutter
  options::Dict{Symbol,Any}
  STLCutter(;options...) = new(options)
end

function cut(cutter::STLCutter,background::DiscreteModel,geom::STLGeometry)
  data,bgf_to_ioc = _cut_stl(background,geom;cutter.options...)
  EmbeddedDiscretization(background, data..., geom), bgf_to_ioc
end

function cut(background::DiscreteModel,geom::STLGeometry)
  cutter = STLCutter()
  cut(cutter,background,geom)
end

function _cut_stl(model::DiscreteModel,geom::STLGeometry;kwargs...)
  subcell_grid, subface_grid, labels = subtriangulate(model,geom;kwargs...)

  inout_dict = Dict{Int8,Int8}(
    FACE_IN => IN, FACE_OUT => OUT, FACE_CUT => CUT, UNSET => OUT )

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

  face_to_io = [ fill(Int8(INTERFACE),num_cells(subface_grid)) ]

  bgface_to_ioc = replace( labels.bgface_to_ioc, inout_dict... )
  bgcell_to_ioc = [ replace( labels.bgcell_to_ioc, inout_dict... ) ]

  oid_to_ls = Dict{UInt,Int}( objectid( get_stl(geom) ) => 1  )
  (bgcell_to_ioc,subcells,cell_to_io,subfacets,face_to_io,oid_to_ls),bgface_to_ioc
end

function send_to_ref_space(grid::Grid,cell_to_bgcell::Vector,subgrid::Grid)
  send_to_ref_space(
    grid,
    cell_to_bgcell,
    get_cell_node_ids(subgrid),
    get_node_coordinates(subgrid))
end

function send_to_ref_space(
  grid::Grid,
  cell_to_bgcell::Vector,
  cell_to_points::AbstractVector,
  points::Vector)

  rpoints = zero(points)
  cache = array_cache(cell_to_points)
  bgcache = array_cache(get_cell_node_ids(grid))
  for i in 1:length(cell_to_points)
    for p in getindex!(cache,cell_to_points,i)
      rpoints[p] = send_to_ref_space!(bgcache,grid,cell_to_bgcell[i],points[p])
    end
  end
  rpoints
end

function send_to_ref_space!(cache,grid::Grid,cell::Integer,point::Point)
  cell_nodes = getindex!(cache,get_cell_node_ids(grid),cell)
  node_coordinates = get_node_coordinates(grid)
  pmin = node_coordinates[cell_nodes[1]]
  pmax = node_coordinates[cell_nodes[end]]
  Point(Tuple(point-pmin)./Tuple(pmax-pmin))
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

