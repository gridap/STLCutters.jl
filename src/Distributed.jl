
function cut(bgmodel::DistributedDiscreteModel,geo::STLGeometry,args...)
  cut(STLCutter(),bgmodel,geo,args...)
end

function cut(cutter::STLCutter,bgmodel::DistributedDiscreteModel,args...)
  D = map(num_dims,local_views(bgmodel)) |> PartitionedArrays.getany
  timers,cutter = setup_distributed_cutter(cutter,bgmodel)
  cell_gids = get_cell_gids(bgmodel)
  facet_gids = get_face_gids(bgmodel,D-1)

  tic!(timers["global"],barrier=true)
  bcells = map(compute_boundary_cells,local_views(bgmodel),local_views(cell_gids))
  icells = map(compute_interior_cells,local_views(bgmodel),local_views(cell_gids))

  # Setup cutters
  tic!(timers["fine"],barrier=true)
  cell_to_facets = map(local_views(bgmodel)) do bgmodel
    compute_cell_to_facets(bgmodel,args...)
  end
  toc!(timers["fine"],"setup")

  # Cut touched parts (only boundary cells)
  tic!(timers["fine"],barrier=true)
  cuts = map(
    local_views(bgmodel),
    local_views(facet_gids),
    bcells,
    cell_to_facets) do bgmodel,facet_gids,bcells,cell_to_facets

    ownmodel = DiscreteModelPortion(bgmodel,bcells)
    cell_to_pcell = get_cell_to_parent_cell(ownmodel)
    facet_to_pfacet = get_face_to_parent_face(ownmodel,D-1)
    cell_to_facets = lazy_map(Reindex(cell_to_facets),bcells)
    _cutter = STLCutter(;cell_to_facets,cutter.options...)
    cutgeo = cut(_cutter,ownmodel,args...)
    cutgeo = change_bgmodel(cutgeo,bgmodel,cell_to_pcell,facet_to_pfacet)
    remove_ghost_subfacets(cutgeo,facet_gids)
  end
  toc!(timers["fine"],"interface")

  # Setup coarse inout and graph
  facet_to_inoutcut = map(compute_bgfacet_to_inoutcut,cuts)
  facet_neighbors = map(compute_facet_neighbors,
    local_views(bgmodel),local_views(cell_gids))
  facet_neighbor_to_ioc = map(compute_facet_neighbor_to_inoutcut,
    local_views(bgmodel),local_views(cell_gids),facet_to_inoutcut,facet_neighbors)

  # Gathers
  root = find_root_part(cuts,bcells)
  part_to_parts = gather(facet_neighbors,destination=root)
  part_to_lpart_to_ioc = gather(facet_neighbor_to_ioc,destination=root)

  # Propagate at coarse level (and complete intersections)
  tic!(timers["coarse"],barrier=true)
  tic!(timers["fine"],barrier=true)
  part_to_ioc = map(
      part_to_parts,
      part_to_lpart_to_ioc,
      local_views(cell_gids)) do part_to_parts,part_to_lpart_to_ioc,ids

    part = part_id(ids)
    if part == root
      propagate_inout(part_to_parts,part_to_lpart_to_ioc)
    else
      Int8[]
    end
  end
  toc!(timers["coarse"],"coarse")

  # Complete cut (interior cells)
  icuts = map(
    local_views(bgmodel),
    local_views(facet_gids),
    icells,
    cell_to_facets) do bgmodel,facet_gids,icells,cell_to_facets

    ownmodel = DiscreteModelPortion(bgmodel,icells)
    cell_to_pcell = get_cell_to_parent_cell(ownmodel)
    facet_to_pfacet = get_face_to_parent_face(ownmodel,D-1)
    cell_to_facets = lazy_map(Reindex(cell_to_facets),icells)
    _cutter = STLCutter(;cell_to_facets,cutter.options...)
    cutgeo = cut(_cutter,ownmodel,args...)
    cutgeo = change_bgmodel(cutgeo,bgmodel,cell_to_pcell,facet_to_pfacet)
    remove_ghost_subfacets(cutgeo,facet_gids)
  end
  toc!(timers["fine"],"interior")

  # Merge discretizations
  cuts = map(cuts,icuts,bcells,icells) do bcut,icut,bcells,icells
    complete_in_or_out!(bcut,icut,bcells,icells)
    merge(bcut,icut,bcells,icells)
  end

  # Scatter
  part_ioc = scatter(part_to_ioc,source=root)

  # Set undefined parts
  map(cuts,part_ioc,cell_gids.partition) do cut,ioc,ids
    own_cells = own_to_local(ids)
    if !istouched(cut,own_cells)
      set_in_or_out!(cut,ioc)
    end
  end

  # Nearest neighbor communication
  consistent_bgcell_to_inoutcut!(cuts,cell_gids)
  consistent_bgfacet_to_inoutcut!(cuts,facet_gids)

  toc!(timers["global"],"global")
  DistributedEmbeddedDiscretization(cuts,bgmodel)
end


function setup_distributed_cutter(a::STLCutter,bgmodel::DistributedDiscreteModel)
  parts = map(part_id,get_cell_gids(bgmodel).partition)
  timers,options = _setup_distributed_cutter(parts;a.options...)
  cutter = STLCutter(;all_defined=false,options...)
  timers,cutter
end

function _setup_distributed_cutter(
  parts;
  verbose=false,
  timers=default_timers(parts;verbose),
  showprogress=false,
  kwargs...)

  kwargs = (;showprogress,kwargs...)
  timers, kwargs
end

function default_timers(parts;kwargs...)
  Dict(
    "global" => PTimer(parts;kwargs...),
    "coarse" => PTimer(parts;kwargs...),
    "fine" => PTimer(parts;kwargs...))
end

"""
    find_root_part(cuts::AbstractArray{<:AbstractEmbeddedDiscretization},cells)

  This functions sets the root processor to a potentially idling processor.
  It returns the first non-cut part.
"""
function find_root_part(
  cuts::AbstractArray{<:AbstractEmbeddedDiscretization},
  cells)

  # TODO: USE STL DATA, e.g. num local STL facets
  touched = map(cuts,cells) do cut,cells
    istouched(cut,cells)
  end
  part_to_touched = gather(touched)
  root_part = map(part_to_touched) do part_to_touched
    if !isempty(part_to_touched)
      p = findfirst(!,part_to_touched)
      isnothing(p) ? 1 : p
    else
      0
    end
  end
  emit!(root_part,root_part)
  PartitionedArrays.getany(root_part)
end

function cut_facets(cut::DistributedEmbeddedDiscretization)
  bgmodel = get_background_model(cut)
  cutfacets = map(cut_facets,local_views(cut))
  DistributedEmbeddedDiscretization(cutfacets,bgmodel)
end

function change_bgmodel(
  cutgeo::STLEmbeddedDiscretization,
  model::DiscreteModel,
  cell_to_newcell=1:num_cells(get_background_model(cutgeo)),
  facet_to_newfacet=1:num_facets(get_background_model(cutgeo)))

  _cut = change_bgmodel(cutgeo.cut,model,cell_to_newcell)
  _cutfacets = change_bgmodel(cutgeo.cutfacets,model,facet_to_newfacet)
  STLEmbeddedDiscretization(_cut,_cutfacets)
end

function change_bgmodel(
  cut::DistributedEmbeddedDiscretization{<:AbstractArray{<:STLEmbeddedDiscretization}},
  model::DistributedDiscreteModel)

  cuts = map(c->c.cut,local_views(cut))
  cutsfacets = map(c->c.cutfacets,local_views(cut))
  _cut = DistributedEmbeddedDiscretization(cuts,get_background_model(cut))
  _cutfacets = DistributedEmbeddedDiscretization(cutsfacets,get_background_model(cut))
  _cut = change_bgmodel(_cut,model)
  _cutfacets = change_bgmodel(_cutfacets,model)

  cuts = map(local_views(_cut),local_views(_cutfacets)) do cut,cutfacets
    STLEmbeddedDiscretization(cut,cutfacets)
  end
  DistributedEmbeddedDiscretization(cuts,model)
end

function remove_ghost_subfacets(cut::STLEmbeddedDiscretization,facet_gids)
  cutfacets = remove_ghost_subfacets(cut.cutfacets,facet_gids)
  STLEmbeddedDiscretization(cut.cut,cutfacets)
end

function get_ls_to_bgcell_to_inoutcut(cut::STLEmbeddedDiscretization)
  get_ls_to_bgcell_to_inoutcut(cut.cut)
end

function get_ls_to_bgfacet_to_inoutcut(cut::STLEmbeddedDiscretization)
  get_ls_to_bgfacet_to_inoutcut(cut.cutfacets)
end

"""
    compute_interior_cells(model::DiscreteModel,indices[,d=0])

  It returns the list of of _own cells_ not touching the subdomain interface.
  See [`compute_cell_to_isboundary`](@ref) for more details.
"""
function compute_interior_cells(
  model::DiscreteModel,
  ids::AbstractLocalIndices,
  args...)

  part = part_id(ids)
  cell_islocal = map(==(part),local_to_owner(ids))
  cell_isboundary = compute_cell_to_isboundary(model,ids,args...)
  cell_isnotboundary = map(!,cell_isboundary)
  cell_isinterior = map(&,cell_isnotboundary,cell_islocal)
  findall(cell_isinterior)
end


"""
    compute_boundary_cells(model::DiscreteModel,indices[,d=0])

  It returns the list of of _own cells_ touching the subdomain interface.
  See [`compute_cell_to_isboundary`](@ref) for more details.
"""
function compute_boundary_cells(
  model::DiscreteModel,
  ids::AbstractLocalIndices,
  args...)

  cell_isboundary = compute_cell_to_isboundary(model,ids,args...)
  findall(cell_isboundary)
end


"""
    compute_cell_to_isboundary(model::DiscreteModel,indices[,d=0])

  It returns a mask whether a cell touches the subdomain interface.

  # Arguments
  - `model::DiscreteModel`: Model of the subdomain.
  - `indices::AbstractLocalIndices`: Partition indices
  - `d::Integer=0`: Dimemension of the d-faces touching the subdomain interface.
"""
function compute_cell_to_isboundary(
  model::DiscreteModel,
  ids::AbstractLocalIndices,
  d::Integer=0)

  part = part_id(ids)
  D = num_dims(model)
  topology = get_grid_topology(model)
  l_to_o = local_to_owner(ids)
  f_to_c = get_faces(topology,d,D)
  c_to_f = get_faces(topology,D,d)
  f_to_p = lazy_map(Broadcasting(Reindex(l_to_o)),f_to_c)
  c1 = array_cache(c_to_f)
  c2 = array_cache(f_to_p)
  c_isboundary = fill(false,num_cells(model))
  for cell in 1:num_cells(model)
    if l_to_o[cell] == part
      for facet in getindex!(c1,c_to_f,cell)
        parts = getindex!(c2,f_to_p,facet)
        if any(==(part),parts) && any(!=(part),parts)
          c_isboundary[cell] = true
          break
        end
      end
    end
  end
  c_isboundary
end

function compute_facet_neighbors(
  model::DiscreteModel,
  ids::AbstractLocalIndices,
  args...)

  D = num_dims(model)
  compute_face_neighbors(model,ids,D-1,args...)
end

"""
    compute_face_neighbors(model::DiscreteModel,indices,d)

  It returns a neighboring graph of the subdomain's neighbors through the interfaces of dimension d.
"""
function compute_face_neighbors(
  model::DiscreteModel,
  ids::AbstractLocalIndices,
  d::Integer)

  part = part_id(ids)
  D = num_dims(model)
  topology = get_grid_topology(model)
  l_to_o = local_to_owner(ids)
  f_to_c = get_faces(topology,d,D)
  f_to_p = lazy_map(Broadcasting(Reindex(l_to_o)),f_to_c)
  c = array_cache(f_to_p)
  neig_parts = Int[]
  for face in 1:length(f_to_c)
    parts = getindex!(c,f_to_p,face)
    if any(==(part),parts) && any(!=(part),parts)
      for neig_part in parts
        if neig_part != part
          union!(neig_parts,neig_part)
        end
      end
    end
  end
  sort(neig_parts)
end

function compute_facet_neighbor_to_inoutcut(
  model::DiscreteModel,
  ids::AbstractLocalIndices,
  args...)

  D = num_dims(model)
  compute_face_neighbor_to_inoutcut(model,ids,D-1,args...)
end

"""
    compute_face_neighbor_to_inoutcut(model::DiscreteModel,indices,d,face_to_inoutcut[,face_neighbors])

  It returns a whether interfaces of dimension ``d`` are in, out, cut or undefined.

  The number of interfaces coincides with the number of neigbors given by [`compute_face_neighbors`](@ref)

!!! note
    If the subdomain is not cut, the neighbors are considered undefined.
"""
function compute_face_neighbor_to_inoutcut(
  model::DiscreteModel,
  ids::AbstractLocalIndices,
  d::Integer,
  face_to_inoutcut::AbstractVector,
  face_neighbors=compute_face_neighbors(model,ids,d))

  part = part_id(ids)
  D = num_dims(model)
  topology = get_grid_topology(model)
  l_to_o = local_to_owner(ids)
  f_to_c = get_faces(topology,d,D)
  f_to_p = lazy_map(Broadcasting(Reindex(l_to_o)),f_to_c)
  c = array_cache(f_to_p)

  n_neigs = length(face_neighbors)
  neig_to_lneig = Dict( face_neighbors .=> 1:n_neigs )
  lneig_to_ioc = fill(Int8(UNDEF),n_neigs)

  for (face,ioc) in enumerate(face_to_inoutcut)
    parts = getindex!(c,f_to_p,face)
    if any(==(part),parts) && any(!=(part),parts)
      for neig_part in parts
        if neig_part != part
          lneig = neig_to_lneig[neig_part]
          neig_ioc = lneig_to_ioc[lneig]
          if neig_ioc != ioc && neig_ioc != UNDEF
            lneig_to_ioc[lneig] = CUT
          else
            lneig_to_ioc[lneig] = ioc
          end
        end
      end
    end
  end
  lneig_to_ioc
end

function propagate_inout(
  part_to_parts::AbstractVector,
  part_to_lpart_to_ioc::AbstractVector)

  part_to_ioc = map(part_to_lpart_to_ioc) do ioc
    any(==(UNDEF),ioc) ? Int8(UNDEF) : Int8(CUT)
  end
  propagate_inout!(part_to_ioc,part_to_parts,part_to_lpart_to_ioc)
  part_to_ioc
end

function propagate_inout!(
  part_to_ioc::AbstractVector,
  part_to_parts::AbstractVector,
  part_to_lpart_to_ioc::AbstractVector)

  propagate_inout!(part_to_ioc,part_to_parts,part_to_lpart_to_ioc,IN)
  replace!(part_to_ioc,UNDEF=>OUT)
  part_to_ioc
end

function propagate_inout!(
  part_to_ioc::AbstractVector,
  part_to_parts::AbstractVector,
  part_to_lpart_to_ioc::AbstractVector,
  in_or_out::Integer)

  c1 = array_cache(part_to_parts)
  c2 = array_cache(part_to_lpart_to_ioc)
  stack = Int32[]
  for (part,ioc) in enumerate(part_to_ioc)
    ioc == CUT || continue
    resize!(stack,0)
    push!(stack,part)
    while !isempty(stack)
      current_part = pop!(stack)
      curr_ioc = part_to_ioc[current_part]
      neig_parts = getindex!(c1,part_to_parts,current_part)
      lneig_to_ioc = getindex!(c2,part_to_lpart_to_ioc,current_part)
      for (neig_part,neig_ioc) in zip(neig_parts,lneig_to_ioc)
        if part_to_ioc[neig_part] == UNDEF &&
           (neig_ioc == in_or_out || curr_ioc != CUT)
          part_to_ioc[neig_part] = in_or_out
          push!(stack,neig_part)
        end
      end
    end
  end
  part_to_ioc
end

"""
    set_in_or_out!(cut::DistributedEmbeddedDiscretization,in_or_out[,cells,facets])

  Sets all the background cells and facets (or a set of `cells` and `facets`) as in or out.
"""
function set_in_or_out!(cut::STLEmbeddedDiscretization,args...)
  set_in_or_out!(cut.cut,args...)
  set_in_or_out!(cut.cutfacets,args...)
  cut
end

function set_in_or_out!(
  cut::EmbeddedDiscretization,
  in_or_out::Integer,
  cells=1:num_cells(cut.bgmodel))

  ls_to_bgc_to_ioc = cut.ls_to_bgcell_to_inoutcut
  @assert length(ls_to_bgc_to_ioc) == 1
  ls_to_bgc_to_ioc[1][cells] .= in_or_out
  cut
end

function set_in_or_out!(
  cut::EmbeddedFacetDiscretization,
  in_or_out::Integer,
  facets=1:num_facets(cut.bgmodel))

  ls_to_f_to_ioc = cut.ls_to_facet_to_inoutcut
  @assert length(ls_to_f_to_ioc) == 1
  ls_to_f_to_ioc[1][facets] .= in_or_out
  cut
end

"""
    complete_in_or_out!(a::AbstractEmbeddedDiscretization,b::AbstractEmbeddedDiscretization,acells,bcells)

 This function considers two discretizations of the same background model on the cells `acells` and `bcells`. These sets of cells have null intersection.
 If only one of the discretizations is not cut it sets the other one as in or out.
"""
function complete_in_or_out!(
  a::AbstractEmbeddedDiscretization,
  b::AbstractEmbeddedDiscretization,
  acells,bcells)

  atouched = istouched(a,acells)
  btouched = istouched(b,bcells)
  if ( atouched || btouched ) && ( !atouched || !btouched ) &&
     ( length(acells) > 0 && length(bcells) > 0 )

    in_or_out = define_interface(a,b,acells,bcells)
    if !atouched
      set_in_or_out!(a,in_or_out,acells)
    else
      set_in_or_out!(b,in_or_out,bcells)
    end
  end
end

function Base.merge(
  a::STLEmbeddedDiscretization,
  b::STLEmbeddedDiscretization,
  acells,bcells)

  afacets,bfacets = _classify_facets(a.cut.bgmodel,acells,bcells)
  cut = merge(a.cut,b.cut,acells,bcells)
  cutfacets = merge(a.cutfacets,b.cutfacets,afacets,bfacets)
  STLEmbeddedDiscretization(cut,cutfacets)
end

function Base.merge(
  a::EmbeddedDiscretization,
  b::EmbeddedDiscretization,
  acells,bcells)

  @check a.bgmodel === b.bgmodel
  ls_to_bgc_to_ioc = deepcopy(a.ls_to_bgcell_to_inoutcut)
  for ls in 1:length(ls_to_bgc_to_ioc)
    ls_to_bgc_to_ioc[ls][acells] .= a.ls_to_bgcell_to_inoutcut[ls][acells]
    ls_to_bgc_to_ioc[ls][bcells] .= b.ls_to_bgcell_to_inoutcut[ls][bcells]
  end
  subcells = append(a.subcells,b.subcells)
  subfacets = append(a.subfacets,b.subfacets)
  ls_to_c_to_io = map(vcat,a.ls_to_subcell_to_inout,b.ls_to_subcell_to_inout)
  ls_to_f_to_io = map(vcat,a.ls_to_subfacet_to_inout,b.ls_to_subfacet_to_inout)
  EmbeddedDiscretization(
    a.bgmodel,
    ls_to_bgc_to_ioc,
    subcells,
    ls_to_c_to_io,
    subfacets,
    ls_to_f_to_io,
    a.oid_to_ls,
    a.geo)
end

function Base.merge(
  a::EmbeddedFacetDiscretization,
  b::EmbeddedFacetDiscretization,
  afacets,bfacets)

  @check a.bgmodel === b.bgmodel
  a = _restrict(a,afacets)
  b = _restrict(b,bfacets)

  ls_to_f_to_ioc = deepcopy(a.ls_to_facet_to_inoutcut)
  for ls in 1:length(ls_to_f_to_ioc)
    ls_to_f_to_ioc[ls][afacets] .= a.ls_to_facet_to_inoutcut[ls][afacets]
    ls_to_f_to_ioc[ls][bfacets] .= b.ls_to_facet_to_inoutcut[ls][bfacets]
  end
  subfacets = append(a.subfacets,b.subfacets)
  ls_to_sf_to_io = map(vcat,a.ls_to_subfacet_to_inout,b.ls_to_subfacet_to_inout)
  EmbeddedFacetDiscretization(
    a.bgmodel,
    ls_to_f_to_ioc,
    subfacets,
    ls_to_sf_to_io,
    a.oid_to_ls,
    a.geo)
end

function append(a::SubCellData,b::SubCellData)
  na = length(a.point_to_coords)
  bc_to_p = b.cell_to_points
  T = typeof(bc_to_p.data)
  bcell_to_points = Table(T(bc_to_p.data .+ na), bc_to_p.ptrs)
  SubCellData(
    append_tables_globally(a.cell_to_points,bcell_to_points),
    vcat(a.cell_to_bgcell,b.cell_to_bgcell),
    vcat(a.point_to_coords,b.point_to_coords),
    vcat(a.point_to_rcoords,b.point_to_rcoords))
end

function append(a::SubFacetData,b::SubFacetData)
  na = length(a.point_to_coords)
  bf_to_p = b.facet_to_points
  T = typeof(bf_to_p.data)
  bfacet_to_points = Table(T(bf_to_p.data .+ na), bf_to_p.ptrs)
  SubFacetData(
    append_tables_globally(a.facet_to_points,bfacet_to_points),
    vcat(a.facet_to_normal,b.facet_to_normal),
    vcat(a.facet_to_bgcell,b.facet_to_bgcell),
    vcat(a.point_to_coords,b.point_to_coords),
    vcat(a.point_to_rcoords,b.point_to_rcoords))
end

function _restrict(cut::EmbeddedFacetDiscretization,newfacets)
  facet_isactive = falses(num_facets(cut.bgmodel))
  facet_isactive[newfacets] .= true
  subfacet_to_facet = cut.subfacets.cell_to_bgcell
  subfacet_isactive = lazy_map(Reindex(facet_isactive),subfacet_to_facet)
  newsubfacets = findall(subfacet_isactive)
  subfacets = SubCellData(cut.subfacets,newsubfacets)
  ls_to_sf_to_io = map(cut.ls_to_subfacet_to_inout) do sf_to_io
    sf_to_io[newsubfacets]
  end
  EmbeddedFacetDiscretization(
    cut.bgmodel,
    cut.ls_to_facet_to_inoutcut,
    subfacets,
    ls_to_sf_to_io,
    cut.oid_to_ls,
    cut.geo)
end

function _classify_facets(model::DiscreteModel,acells,bcells)
  D = num_dims(model)
  topo = get_grid_topology(model)
  c_to_f = get_faces(topo,D,D-1)
  f_to_ab = zeros(Int8,num_facets(model))
  A,B = (1,2)
  c = array_cache(c_to_f)
  for cell in acells
    for facet in getindex!(c,c_to_f,cell)
      f_to_ab[facet] = A
    end
  end
  for cell in bcells
    for facet in getindex!(c,c_to_f,cell)
      f_to_ab[facet] = B
    end
  end
  afacets = findall(==(A),f_to_ab)
  bfacets = findall(==(B),f_to_ab)
  afacets,bfacets
end

istouched(cut::STLEmbeddedDiscretization,args...) = istouched(cut.cut,args...)


"""
    istouched(cut::AbstractEmbeddedDiscretization[,cells])

  Check is an embedded discretization `cut` (or a set of `cells`) have information about the geometry. In other words, it checks if any cell is cut.
"""
function istouched(
  cut::EmbeddedDiscretization,
  cells=1:num_cells(get_background_model(cut)),
  ls=1)

  c_to_ioc = lazy_map(Reindex(cut.ls_to_bgcell_to_inoutcut[ls]),cells)
  !all(==(UNDEF),c_to_ioc)
end

"""
    define_interface(::STLEmbeddedDiscretization,::STLEmbeddedDiscretization,acells,bcells)

  This function considers two discretizations of the same background model on the cells `acells` and `bcells`. These sets of cells have null intersection.
  See also [`complete_in_or_out!`](@ref) for more details.

  It returns whether the interface between two discretizations is in or out.
"""

function define_interface(
  a::STLEmbeddedDiscretization,
  b::STLEmbeddedDiscretization,acells,bcells)

  @check length(acells) > 0
  @check length(bcells) > 0
  amodel = get_background_model(a)
  bmodel = get_background_model(b)
  @check amodel === bmodel
  facets = interface_facets(amodel,acells,bcells)
  c = istouched(a,acells) ? a : b
  facet_to_ioc = compute_bgfacet_to_inoutcut(c)
  f_to_ioc = lazy_map(Reindex(facet_to_ioc),facets)
  n_in = count(==(IN),f_to_ioc)
  n_out = count(==(OUT),f_to_ioc)
  if n_in > 0
    @check n_out == 0
    IN
  else
    @check n_out > 0
    OUT
  end
end

function interface_facets(model::DiscreteModel,acells,bcells)
  D = num_dims(model)
  topo = get_grid_topology(model)
  c_to_f = get_faces(topo,D,D-1)
  f_to_a = falses(num_facets(model))
  f_to_b = falses(num_facets(model))
  c = array_cache(c_to_f)
  for cell in acells
    for facet in getindex!(c,c_to_f,cell)
      f_to_a[facet] = true
    end
  end
  for cell in bcells
    for facet in getindex!(c,c_to_f,cell)
      f_to_b[facet] = true
    end
  end
  f_to_ab = lazy_map(&,f_to_a,f_to_b)
  findall(f_to_ab)
end

# PR @ GridapDistributed.jl

function Gridap.ReferenceFEs.simplexify(model::DistributedDiscreteModel{D};kwargs...) where D
  models = map(local_views(model)) do m
    Gridap.ReferenceFEs.simplexify(m;kwargs...)
  end
  gids = generate_simplex_cell_gids(model)
  DistributedDiscreteModel(models,gids)
end

function generate_simplex_cell_gids(model::DistributedDiscreteModel)
  gids = get_cell_gids(model)
  ngcells = num_cells(model)
  partition = map(local_views(model),gids.partition) do model,ids
    s, = simplexify(get_polytope(only(get_reffes(model))))
    ns_x_c = length(s)
    nc = num_cells(model)
    ngscells = ngcells*ns_x_c
    sc_to_c = reduce(append!,lazy_map(i->Fill(i,ns_x_c),1:nc),init=Int[])
    sc_to_i = reduce(append!,lazy_map(i->1:ns_x_c,1:nc),init=Int[])
    sc_to_gc = map(Reindex(local_to_global(ids)),sc_to_c)
    sc_to_gsc = map((gc,i)-> (gc-1)*ns_x_c+i,sc_to_gc,sc_to_i)
    sc_to_o = map(Reindex(local_to_owner(ids)),sc_to_c)
    LocalIndices(ngscells,part_id(ids),sc_to_gsc,sc_to_o)
  end
  PRange(partition)
end

function distributed_aggregate(
  strategy::AggregateCutCellsByThreshold,
  cut::DistributedEmbeddedDiscretization,
  geo::STLGeometry,
  in_or_out=IN)

  facet_to_inoutcut = compute_bgfacet_to_inoutcut(cut,geo)
  GridapEmbedded.Distributed._distributed_aggregate_by_threshold(strategy.threshold,cut,geo,in_or_out,facet_to_inoutcut)
end
