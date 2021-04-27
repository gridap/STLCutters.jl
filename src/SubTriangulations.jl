
function compute_submesh(
  bgmodel::CartesianDiscreteModel,
  stlmodel::DiscreteModel;
  threading=:spawn,
  kdtree=false,
  tolfactor=1e3)

  grid = get_grid(bgmodel)
  grid_topology = get_grid_topology(bgmodel)
  p = get_polytope(only(get_reffes(bgmodel)))
  node_to_coords = get_node_coordinates(grid)
  cell_to_nodes = get_cell_node_ids(grid)
  stl = STL(stlmodel)
  D = num_dims(grid)
  atol = eps(grid)*tolfactor

  f_to_isempty = get_facet_to_isempty(stl;atol)
  Πf = get_facet_planes(stl)
  Πf = correct_small_facets_planes!(stl,Πf,f_to_isempty;atol)
  Πr = get_reflex_planes(stl,Πf)

  c_to_stlf = compute_cell_to_facets(grid,stl)

  Γ0 = Polyhedron(stl)

  submesh = _get_threaded_empty_arrays(stl)
  io_arrays = _init_io_arrays(bgmodel,p)
  caches = _get_threaded_caches(cell_to_nodes)

  cut_cells = filter(i->!isempty(c_to_stlf[i]),1:num_cells(grid))
  if threading == :threads
    Threads.@threads for cell in cut_cells
      save_cell_submesh!(submesh,io_arrays,stl,p,cell,
        compute_polyhedra!(caches,Γ0,stl,p,f_to_isempty,Πf,Πr,
          c_to_stlf,node_to_coords,cell_to_nodes,cell;atol,kdtree)... )
    end
  elseif threading == :spawn
    @sync for cell in cut_cells
      Threads.@spawn save_cell_submesh!(submesh,io_arrays,stl,p,cell,
        compute_polyhedra!(caches,Γ0,stl,p,f_to_isempty,Πf,Πr,
          c_to_stlf,node_to_coords,cell_to_nodes,cell;atol,kdtree)... )
    end
  else
    @unreachable
  end

  submesh = _append_threaded_submesh!(submesh)
  io_arrays = _reduce_io_arrays(bgmodel,io_arrays)
  bgcell_to_ioc, bgnode_to_io, bgfacet_to_ioc = io_arrays
  T,F,X,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf = submesh

  propagate_inout!(bgmodel,bgcell_to_ioc,bgnode_to_io) 
  set_facets_as_inout!(bgmodel,bgcell_to_ioc,bgfacet_to_ioc)
  
  delete_small_subcells!(bgmodel,T,X,k_to_io,k_to_bgcell)
  delete_small_subfacets!(bgmodel,F,Xf,f_to_bgcell,f_to_stlf)
  T,X,F,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf,bgcell_to_ioc,bgfacet_to_ioc
end

function compute_cell_to_facets(grid::CartesianGrid,stl::STL)
  desc = get_cartesian_descriptor(grid)
  @assert length(get_reffes(grid)) == 1
  p = get_polytope(get_cell_reffe(grid)[1])
  @notimplementedif desc.map !== identity
  cell_to_stl_facets = [ Int32[] for _ in 1:num_cells(grid) ]
  n = Threads.nthreads()
  thread_to_cells = [ Int32[] for _ in 1:n ]
  thread_to_stl_facets = [ Int32[] for _ in 1:n ]
  coords = get_node_coordinates(grid)
  cell_to_nodes = get_cell_node_ids(grid)
  c = [ ( get_cell_cache(stl), array_cache(cell_to_nodes) ) for _ in 1:n ]
  δ = 0.1
  Threads.@threads for stl_facet in 1:num_cells(stl)
    i = Threads.threadid()
    cc,nc = c[i]
    f = get_cell!(cc,stl,stl_facet)
    pmin,pmax = get_bounding_box(f)
    for cid in get_cells_around(desc,pmin,pmax)
      cell = LinearIndices(desc.partition)[cid.I...]
      nodes = getindex!(nc,cell_to_nodes,cell)
      _pmin = coords[nodes[1]]
      _pmax = coords[nodes[end]]
      Δ = (_pmax - _pmin) * δ
      _pmin = _pmin - Δ
      _pmax = _pmax + Δ
      if voxel_intersection(f,_pmin,_pmax,p)
        push!(thread_to_cells[i],cell)
        push!(thread_to_stl_facets[i],stl_facet)
      end
    end
  end
  cell_to_stl_facets = [ Int32[] for _ in 1:num_cells(grid) ]
  for (cells,stl_facets) in zip(thread_to_cells,thread_to_stl_facets)
    for (cell,stl_facet) in zip(cells,stl_facets)
      push!(cell_to_stl_facets[cell],stl_facet)
    end
  end
  cell_to_stl_facets
end

function compute_polyhedra!(caches,Γ0,stl,p,f_to_isempty,Πf,Πr,
  c_to_stlf,node_to_coords,cell_to_nodes,cell;atol,kdtree)

  cell_coords = _get_cell_coordinates!(caches,node_to_coords,cell_to_nodes,cell)
  pmin = cell_coords[1] - atol
  pmax = cell_coords[end] + atol
   _faces = _get_cell_faces!(caches,stl,c_to_stlf,cell,f_to_isempty)
  facets,expanded_facets,empty_facets,reflex_faces = _faces

  Πk,Πk_ids,Πk_io = get_cell_planes(p,pmin,pmax)
  Πkf,Πkf_ids = filter_face_planes(stl,Πr,reflex_faces,Πf,facets) 

  K = Polyhedron(p,cell_coords)
  Γk0 = restrict(Γ0,stl,expanded_facets)

  compute_distances!(Γk0,lazy_append(Πk,Πkf),lazy_append(Πk_ids,Πkf_ids))
  compute_distances!(K,Πkf,Πkf_ids)
  Π_to_refΠ,Πs,inv_Π = link_planes!(Γk0,stl,Πr,Πf;atol)
  invert_plane_distances!(Γk0,Πs,inv_Π)
  invert_plane_distances!(K,Πs,inv_Π)
  set_linked_planes!(Γk0,Π_to_refΠ,Πs)
  set_linked_planes!(K,Π_to_refΠ,Πs)

  Γk = clip(Γk0,Πk_ids,inout=Πk_io)
 
  if isnothing(Γk) || isempty(get_original_facets(Γk,stl))
    nothing,nothing
  else
    if kdtree
      Γv,Kv = refine_by_vertices(Γk,K,atol)
    else
      Γv,Kv = [Γk],[K]
    end
    Kn_in = typeof(K)[]
    Kn_out = typeof(K)[]
    for (γk,k) in zip(Γv,Kv)
      k_in = refine(k,γk,stl,reflex_faces,empty_facets,inside=true)
      k_out = refine(k,γk,stl,reflex_faces,empty_facets,inside=false)
      append!(Kn_in,k_in)
      append!(Kn_out,k_out)
    end

    Kn_in,Kn_out
  end
end

function save_cell_submesh!(submesh,io_arrays,stl,p,cell,Kn_in,Kn_out)
  !isnothing(Kn_in) || return
  bgcell_to_ioc, bgcell_node_to_io, bgcell_facet_to_ioc = io_arrays
  Tin,Xin = simplexify(Kn_in)
  Tout,Xout = simplexify(Kn_out)
  T_Γ,X_Γ,f_to_f = simplexify_boundary(Kn_in,stl)
  bgcell_to_ioc[cell] = _get_cell_io(T_Γ,Kn_in,Kn_out)
  bgcell_to_ioc[cell] == FACE_CUT || return
  D = num_dims(stl)
  f_to_f .-= get_offset(stl,D)
  n_to_io = get_cell_nodes_to_inout(Kn_in,Kn_out,p)
  f_to_ioc = get_cell_facets_to_inoutcut(Kn_in,Kn_out,p)
  _append_submesh!(submesh,Xin,Tin,Xout,Tout,X_Γ,T_Γ,f_to_f,cell)
  for i in 1:num_vertices(p)
    bgcell_node_to_io[i,cell] = n_to_io[i]
  end
  for i in 1:num_facets(p)
    bgcell_facet_to_ioc[i,cell] = f_to_ioc[i]
  end
end

function propagate_inout!(bgmodel,bgcell_to_ioc,bgnode_to_io)
  D = num_dims(bgmodel)
  grid_topology = get_grid_topology(bgmodel)
  stack = Int32[]
  n_to_c = get_faces(grid_topology,0,D)
  c_to_n = get_faces(grid_topology,D,0)
  node_cache = array_cache(c_to_n)
  neig_node_cache = array_cache(c_to_n)
  neig_cell_cache = array_cache(n_to_c)
  for cell in 1:num_cells(grid_topology)
    if bgcell_to_ioc[cell] ∈ (FACE_CUT,FACE_IN)
      resize!(stack,0)
      push!(stack,cell)
      while !isempty(stack)
        current_cell = pop!(stack)
        for node in getindex!(node_cache,c_to_n,cell)
          if bgnode_to_io[node] == FACE_IN || bgcell_to_ioc[cell] == FACE_IN
            for neig_cell in getindex!(neig_cell_cache,n_to_c,node)
              if bgcell_to_ioc[neig_cell] == UNSET
                bgcell_to_ioc[neig_cell] = FACE_IN
                for neig_node in getindex!(neig_node_cache,c_to_n,neig_cell)
                  if bgnode_to_io[neig_node] == UNSET
                    bgnode_to_io[neig_node] = FACE_IN
                  end
                  push!(stack,neig_cell)
                end
              end
            end
          end
        end
      end
    end
  end
  replace!(bgcell_to_ioc, UNSET => FACE_OUT )
  bgcell_to_ioc
end

function link_planes!(surf::Polyhedron,stl::STL,Πr,Πf;atol)
  planes = surf.data.plane_to_ids
  default_out = fill(UNSET,length(planes)), planes, falses(length(planes))
  Π_to_faces = find_faces_on_planes!(surf,stl;atol)
  maximum(length,Π_to_faces) > 0 || return default_out
  Π_to_coplanar_Π = link_coplanar_planes(surf,stl,Π_to_faces;atol)
  maximum(length,Π_to_coplanar_Π) > 0 || return default_out
  Π_to_ref_Π =  group_coplanar_planes(planes,Π_to_coplanar_Π)
  merge_faces!(Π_to_faces,planes,Π_to_ref_Π)
  correct_distances!(surf,Π_to_ref_Π,Π_to_faces)
  Π_to_ref_Π = mark_inverted_planes!(Π_to_ref_Π,planes,Πr,Πf,stl;atol)
  inverted_planes = invert_reflex_planes!(Π_to_ref_Π,planes,stl)
  
  Π_to_ref_Π, planes, inverted_planes
end

function _get_threaded_empty_arrays(stl::STL)
  n = Threads.nthreads()
  [ _get_empty_arrays(stl::STL) for _ in 1:n ]
end

function _get_empty_arrays(stl::STL)
  P = eltype( get_vertex_coordinates(stl) )
  T = Vector{Int32}[]
  F = Vector{Int32}[]
  X = P[]
  Xf = P[] 
  k_to_io = Int8[]
  k_to_bgcell = Int32[]
  f_to_bgcell = Int32[]
  f_to_stlf = Int32[]
  T,F,X,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf
end

function _threaded_empty_array(::Type{T}) where T
  n = Threads.nthreads()
  [ T[] for i in 1:n ]
end

function _init_io_arrays(model::DiscreteModel,p::Polytope)
  D = num_dims(model)
  ncells = num_cells(model)
  nnodes_x_cell = num_vertices(p)
  nfacets_x_cell = num_facets(p)

  bgcell_to_ioc = fill(Int8(UNSET),ncells)
  bgcell_node_to_io = fill(Int8(UNSET),nnodes_x_cell,ncells)
  bgcell_facet_to_ioc = fill(Int8(UNSET),nfacets_x_cell,ncells)
  bgcell_to_ioc, bgcell_node_to_io, bgcell_facet_to_ioc
end

function _get_caches(c_to_n)
  nc = array_cache(c_to_n)
  facets = Int32[]
  expanded_facets = Int32[]
  empty_facets = Int32[]
  reflex_faces = Int32[]
  nc,facets,expanded_facets,empty_facets,reflex_faces
end

function _get_threaded_caches(c_to_n)
  n = Threads.nthreads()
  [ _get_caches(c_to_n) for _ in 1:n ]
end

function _get_cell_coordinates!(caches,node_to_coords,cell_to_nodes,cell)
  i = Threads.threadid()
  nc, = caches[i]
  cell_nodes = getindex!(nc,cell_to_nodes,cell)
  cell_coords = view(node_to_coords,cell_nodes)
  map(i->Point(float.(Tuple(i))),cell_coords)
end

function _get_cell_faces!(caches,stl,c_to_stlf,cell,f_to_isempty)
  i = Threads.threadid()
  _,facets,expanded_facets,empty_facets,reflex_faces = caches[i]
  D = num_dims(stl)
  get_expanded_facets!(expanded_facets,stl,c_to_stlf,cell)
  get_reflex_faces!(reflex_faces,stl,c_to_stlf,cell)
  copy!(facets,c_to_stlf[cell])
  copy!(empty_facets,c_to_stlf[cell])
  filter!(f->f_to_isempty[f],empty_facets)
  map!(i->i+get_offset(stl,D),facets,facets)
  map!(i->i+get_offset(stl,D),empty_facets,empty_facets)
  map!(i->i+get_offset(stl,D-1),reflex_faces,reflex_faces)
  facets,expanded_facets,empty_facets,reflex_faces
end

function _get_cell_io(T_Γ,Kn_in,Kn_out)
  if length(T_Γ) > 0
    FACE_CUT
  else
    if length(Kn_in) > 0 && length(Kn_out) == 0
      FACE_IN
    elseif length(Kn_in) == 0 && length(Kn_out) > 0
      FACE_OUT
    else
      UNSET
    end
  end
end

function _append_submesh!(submesh,Xin,Tin,Xout,Tout,Xfn,Tfn,fn_to_f,bgcell)
  i = Threads.threadid()
  _submesh = submesh[i]
  _append_subcells!(_submesh,Xin,Tin,FACE_IN,bgcell)
  _append_subcells!(_submesh,Xout,Tout,FACE_OUT,bgcell)
  _append_subfacets!(_submesh,Xfn,Tfn,fn_to_f,bgcell)
end

function _append_subcells!(submesh_arrays,Xn,Tn,io,bgcell)
  T,F,X,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf = submesh_arrays
  append!(T, map(i->i.+length(X),Tn) )
  append!(X,Xn)
  append!(k_to_io,fill(io,length(Tn)))
  append!(k_to_bgcell,fill(bgcell,length(Tn)))
  submesh_arrays
end

function _append_subfacets!(submesh_arrays,Xfn,Tfn,fn_to_f,bgcell)
  T,F,X,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf = submesh_arrays
  append!(F, map(i->i.+length(Xf),Tfn) ) 
  append!(Xf,Xfn)
  append!(f_to_bgcell,fill(bgcell,length(Tfn)))
  append!(f_to_stlf,fn_to_f)
  submesh_arrays
end

function _append_threaded_submesh!(submesh)
  n = Threads.nthreads()
  for i in 2:n
    _append!(submesh[1],submesh[i])
  end
  submesh[1]
end

function _append!(submesh_a,submesh_b)
  T,F,X,Xf, = submesh_a
  T,F, = submesh_b
  nnodes_cells_a = length(X)
  nnodes_facets_a = length(Xf)
  map!(i->map!(j->j+nnodes_cells_a,i,i),T,T)
  map!(i->map!(j->j+nnodes_facets_a,i,i),F,F)
  for (a,b) in zip(submesh_a,submesh_b)
    append!(a,b)
  end
  submesh_a
end

function _reduce_io_arrays(model::DiscreteModel,io_arrays)
  bgcell_to_io, bgcell_node_to_io, bgcell_facet_to_io = io_arrays
  D = num_dims(model)
  topo = get_grid_topology(model)
  bgnode_to_io = fill(Int8(UNSET),num_vertices(topo) )
  bgfacet_to_io = fill(Int8(UNSET),num_facets(topo) )
  c_to_n = get_cell_vertices(topo)
  c_to_f = get_faces(topo,D,D-1)
  nc = array_cache(c_to_n)
  nf = array_cache(c_to_f)
  for cell in 1:num_cells(model)
    for (i,bg_n) in enumerate( getindex!(nc,c_to_n,cell) )
      io = bgcell_node_to_io[i,cell]
      io ≠ UNSET || continue
      bgnode_to_io[bg_n] ≠ FACE_CUT || continue
      if bgnode_to_io[bg_n] ∉ (UNSET,io)
        io = FACE_CUT
      end
      bgnode_to_io[bg_n] = io
    end
    for (i,bg_f) in enumerate( getindex!(nf,c_to_f,cell) )
      io = bgcell_facet_to_io[i,cell]
      if bgfacet_to_io[bg_f] ∉ (UNSET,io)
        io = FACE_CUT
      end
      bgfacet_to_io[bg_f] = io
    end
  end
  bgcell_to_io, bgnode_to_io, bgfacet_to_io
end


function set_facets_as_inout!(bgmodel,bgcell_to_ioc,bgfacet_to_ioc)
  D = num_dims(bgmodel)
  grid_topology = get_grid_topology(bgmodel)
  c_to_f = get_faces(grid_topology,D,D-1)
  cache = array_cache(c_to_f)
  for cell in 1:num_cells(grid_topology)
    bgcell_to_ioc[cell] ≠ FACE_CUT || continue
    for facet in getindex!(cache,c_to_f,cell)
      bgfacet_to_ioc[facet] = bgcell_to_ioc[cell]
    end
  end
  bgfacet_to_ioc
end



function get_expanded_facets!(facets,stl,c_to_stlf,cell)
  D = num_dims(stl)
  f_to_v = get_faces(stl,D,0)
  v_to_f = get_faces(stl,0,D)
  nv = array_cache(f_to_v)
  nf = array_cache(v_to_f)
  empty!(facets)
  for f in c_to_stlf[cell] 
    for v in getindex!(nv,f_to_v,f)
      for _f in getindex!(nf,v_to_f,v)
        if _f ∉ facets
          push!(facets,_f)
        end
      end
    end
  end
  facets
end

function get_reflex_faces!(rfaces,stl,c_to_stlf,cell)
  D = num_dims(stl)
  f_to_r = get_faces(stl,D,D-1)
  r_to_f = get_faces(stl,D-1,D)
  rc = array_cache(f_to_r)
  fc = array_cache(r_to_f)
  empty!(rfaces)
  for f in c_to_stlf[cell]
    for r in getindex!(rc,f_to_r,f)
      facets = getindex!(fc,r_to_f,r)
      r += get_offset(stl,D-2)
      if r ∉ rfaces && all( i -> i ∈ c_to_stlf[cell], facets )
        push!(rfaces,r)
      end
    end
  end
  rfaces
end
