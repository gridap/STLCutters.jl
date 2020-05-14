
struct STLCutter <: Interfaces.Cutter end

function cut(cutter::STLCutter,background::Gridap.DiscreteModel,geom::STLGeometry)
  data = _cut_sm(background,geom)
  EmbeddedDiscretization(background, data..., geom)
end

function cut(background::Gridap.DiscreteModel,geom::STLGeometry)
  cutter = STLCutter()
  cut(cutter,background,geom)
end

function _cut_sm(model::Gridap.DiscreteModel,geom::STLGeometry)
  mesh = CartesianMesh(model)
  stl = get_stl(geom)
  sm = SurfaceMesh(stl)

  #
  bulk = BulkMesh(mesh,sm)
  sm_to_bgcell_to_inoutcut = zeros(Int8,num_cells(mesh))
  #TODO: implement in a separated funct
  for cell in 1:num_cells(mesh)
    if is_cell_interior(bulk,cell)
      sm_to_bgcell_to_inoutcut[cell] = IN
    elseif is_cell_exterior(bulk,cell)
      sm_to_bgcell_to_inoutcut[cell] = OUT
    elseif is_cell_cut(bulk,cell)
      sm_to_bgcell_to_inoutcut[cell] = CUT
    else
      @assert false
    end
  end
  sm_to_bgcell_to_inoutcut = [ sm_to_bgcell_to_inoutcut ]
  subcells = convert(Interfaces.SubTriangulation,bulk.subtriangulation)
  sm_to_subcell_to_inout = zeros(Int8,num_cells(bulk.subtriangulation))
  #TODO: implement in a separated funct
  for subcell in 1:num_cells(bulk.subtriangulation)
    if is_cell_interior(bulk.subtriangulation,subcell)
      sm_to_subcell_to_inout[subcell] = IN
    elseif is_cell_exterior(bulk.subtriangulation,subcell)
      sm_to_subcell_to_inout[subcell] = OUT
    else
      @assert false
    end
  end
  sm_to_subcell_to_inout = [ sm_to_subcell_to_inout ]
  subfacets = convert(Interfaces.FacetSubTriangulation,bulk.facet_subtriangulations[1])
  sm_to_subfacet_to_inout = [ fill(Int8(INTERFACE),num_facets(bulk.facet_subtriangulations[1])) ]
  oid_to_ls = Dict{UInt,Int}( objectid( get_stl(geom) ) => 1  )
  data = sm_to_bgcell_to_inoutcut,subcells,sm_to_subcell_to_inout,subfacets,sm_to_subfacet_to_inout,oid_to_ls
  #

  data
end


function CartesianMesh(a::Gridap.CartesianDiscreteModel)
  desc = get_cartesian_descriptor(a)
  origin = Point(Tuple(desc.origin))
  desc_sizes = Tuple(desc.partition).*Tuple(desc.sizes)
  sizes = VectorValue(Tuple(desc_sizes))
  partition = Tuple(desc.partition)
  CartesianMesh( origin, sizes, partition )
end

Base.convert(::Type{T},a::Gridap.Point) where T<:Point = T(Tuple(a))

Base.convert(::Type{T},a::Gridap.Point) where T<:VectorValue = T(Tuple(a))

function Base.convert(::Type{T},a::Gridap.Arrays.Table) where T<:Table
  data = a.data
  ptrs = a.ptrs
  T(data,ptrs)
end

Base.convert(::Type{T},a::Point) where T<:Gridap.Point = T(get_data(a))

Base.convert(::Type{T},a::VectorValue) where T<:Gridap.Point = T(get_data(a))

function Base.convert(::Type{Gridap.Arrays.Table{D,P}},a::Table) where {D,P}
  data = convert(Vector{D},a.data)
  ptrs = convert(Vector{P},a.ptrs)
  Gridap.Arrays.Table(data,ptrs)
end

function Base.convert(::Type{S},a::STLCutters.SubTriangulation{D,T}) where {D,T,S<:SubTriangulation}
  cell_to_points = convert(Gridap.Arrays.Table{Int,Int32},a.cell_to_points)
  cell_to_bgcell = convert(Vector{Int32},a.cell_to_bgcell)
  point_to_coords = convert(Vector{Gridap.Point{D,T}},a.point_to_coords)
  point_to_rcoords = convert(Vector{Gridap.Point{D,T}},a.point_to_rcoords)
  S(cell_to_points,cell_to_bgcell,point_to_coords,point_to_rcoords)
end

function Base.convert(::Type{F},a::STLCutters.FacetSubTriangulation{D,T}) where {D,T,F<:FacetSubTriangulation}
  facet_to_points = convert(Gridap.Arrays.Table{Int,Int32},a.facet_to_points)
  facet_to_normal = convert(Vector{Gridap.Point{D,T}},a.facet_to_normal)
  facet_to_bgcell = convert(Vector{Int32},a.facet_to_bgcell)
  point_to_coords = convert(Vector{Gridap.Point{D,T}},a.point_to_coords)
  point_to_rcoords = convert(Vector{Gridap.Point{D,T}},a.point_to_rcoords)
  F(facet_to_points,facet_to_normal,facet_to_bgcell,point_to_coords,point_to_rcoords)
end

