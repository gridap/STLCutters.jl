function cut(cutter::STLCutter,bgmodel::DistributedDiscreteModel,args...)
  D = map(num_dims,local_views(bgmodel)) |> PartitionedArrays.getany
  cutter = STLCutter(;all_defined=false,cutter.options...)
  cell_gids = get_cell_gids(bgmodel)
  facet_gids = get_face_gids(bgmodel,D-1)

  bcells = map(compute_boundary_cells,local_views(bgmodel),local_views(cell_gids))
  icells = map(compute_interior_cells,local_views(bgmodel),local_views(cell_gids))

  # Cut touched parts (only boundary cells)
  cuts = map(
    local_views(bgmodel),
    local_views(cell_gids),
    local_views(facet_gids),
    bcells) do bgmodel,cell_gids,facet_gids,bcells

    # bcells = compute_boundary_cells(bgmodel,cell_gids)
    ownmodel = DiscreteModelPortion(bgmodel,bcells)
    # ownmodel = remove_ghost_cells(bgmodel,cell_gids)
    cell_to_pcell = get_cell_to_parent_cell(ownmodel)
    facet_to_pfacet = get_face_to_parent_face(ownmodel,D-1)
    cutgeo = cut(cutter,ownmodel,args...)
    cutgeo = change_bgmodel(cutgeo,bgmodel,cell_to_pcell,facet_to_pfacet)
    remove_ghost_subfacets(cutgeo,facet_gids)
  end

  # Setup coarse inout and graph
  facet_to_inoutcut = map(compute_bgfacet_to_inoutcut,cuts)
  facet_neighbors = map(compute_facet_neighbors,
    local_views(bgmodel),local_views(cell_gids))
  facet_neighbor_to_ioc = map(compute_facet_neighbor_to_inoutcut,
    local_views(bgmodel),local_views(cell_gids),facet_to_inoutcut,facet_neighbors)

  # Gathers
  root = find_root_part(cuts)
  part_to_parts = gather(facet_neighbors,destination=root)
  part_to_lpart_to_ioc = gather(facet_neighbor_to_ioc,destination=root)

  # Propagate at coarse level (and complete intersections)
  part_to_ioc = map(
      part_to_parts,
      part_to_lpart_to_ioc,
      local_views(cell_gids)) do part_to_parts,part_to_lpart_to_ioc,ids

    part = part_id(ids)
    if part == root
      propagate_inout(part_to_parts,part_to_lpart_to_ioc)
    end
  end

  # Complete cut (interior cells)
  icuts = map(
    local_views(bgmodel),
    local_views(cell_gids),
    local_views(facet_gids),
    icells) do bgmodel,cell_gids,facet_gids,icells

    # icells = compute_interior_cells(bgmodel,cell_gids)
    ownmodel = DiscreteModelPortion(bgmodel,icells)
    cell_to_pcell = get_cell_to_parent_cell(ownmodel)
    facet_to_pfacet = get_face_to_parent_face(ownmodel,D-1)
    cutgeo = cut(cutter,ownmodel,args...)
    cutgeo = change_bgmodel(cutgeo,bgmodel,cell_to_pcell,facet_to_pfacet)
    remove_ghost_subfacets(cutgeo,facet_gids)
  end
  # Intersect interior and merge discretizations

  cuts = map(cuts,icuts,bcells,icells) do bcut,icut,bcells,icells
    complete_inout!(bcut,icut,bcells,icells)
    merge(bcut,icut,bcells,icells)
  end

  # Scatter
  part_ioc = scatter(part_to_ioc,source=root)

  # Set undefined parts
  map(cuts,part_ioc) do cut,ioc
    if ioc != CUT
      set_inout!(cut,ioc)
    end
  end

  # Nearest neighbor communication
  consistent_bgcell_to_inoutcut!(cuts,cell_gids)
  consistent_bgfacet_to_inoutcut!(cuts,facet_gids)
  DistributedEmbeddedDiscretization(cuts,bgmodel)
end

function find_root_part(cuts)
  # ( SELECT AN EMPTY PART )
  1
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

function compute_boundary_cells(
  model::DiscreteModel,
  ids::AbstractLocalIndices,
  args...)

  cell_isboundary = compute_cell_to_isboundary(model,ids,args...)
  findall(cell_isboundary)
end

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
    any(==(UNDEF),ioc) ? UNDEF : CUT
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

function set_inout!(cut::STLEmbeddedDiscretization,in_or_out::Integer)
  set_inout!(cut.cut,in_or_out)
  set_inout!(cut.cutfacets,in_or_out)
  cut
end

function set_inout!(cut::EmbeddedDiscretization,in_or_out::Integer)
  ls_to_bgc_to_ioc = cut.ls_to_bgcell_to_inoutcut
  @assert length(ls_to_bgc_to_ioc) == 1
  ls_to_bgc_to_ioc[1] .= in_or_out
  cut
end

function set_inout!(cut::EmbeddedFacetDiscretization,in_or_out::Integer)
  ls_to_f_to_ioc = cut.ls_to_facet_to_inoutcut
  @assert length(ls_to_f_to_ioc) == 1
  ls_to_f_to_ioc[1] .= in_or_out
  cut
end

function complete_inout!(
  a::AbstractEmbeddedDiscretization,
  b::AbstractEmbeddedDiscretization,
  acells,bcells)

  atouched = istouched(a,acells)
  btouched = istouched(b,bcells)
  if ( atouched || btouched ) && ( !atouched || !btouched )
    in_or_out = define_interface(a,b,acell,bcells)
    if atouched
      set_inout!(a,in_or_out)
    else
      set_inout!(b,in_or_out)
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
  f_to_ab = zeros(num_facets(model))
  a,b = (1,2)
  c = array_cache(c_to_f)
  for cell in acells
    for facet in getindex!(c,c_to_f,cell)
      f_to_ab[facet] = a
    end
  end
  for cell in bcells
    for facet in getindex!(c,c_to_f,cell)
      f_to_ab[facet] = b
    end
  end
  afacets = findall(==(a),f_to_ab)
  bfacets = findall(==(b),f_to_ab)
  afacets,bfacets
end

istouched(cut::STLEmbeddedDiscretization,args...) = istouched(cut.cut,args...)

function istouched(cut::EmbeddedDiscretization,cells,ls=1)
  c_to_ioc = lazy_map(Reindex(cut.ls_to_bgcell_to_inoutcut[ls]),cells)
  !all(==(UNDEF),c_to_ioc)
end


function define_interface(
  a::EmbeddedDiscretization,
  b::EmbeddedDiscretization,acells,bcells)

  facets = interface_facets(model,acells,bcells)
  if istouched(a,acells)
    facet_to_ioc = compute_bgfacet_to_inoutcut(a)
  else
    facet_to_ioc = compute_bgfacet_to_inoutcut(b)
  end
  f_to_ioc = lazy_map(Reindex(bgfacet_to_ioc),facets)
  n_in = count(==(IN),f_to_ioc)
  n_out = count(==(OUT),f_to_ioc)

  if n_in > 0
    @assert n_out == 0
    IN
  else
    @assert n_out > 0
    OUT
  end
end
