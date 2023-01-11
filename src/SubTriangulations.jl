
const DEFAULT_TOL_FACTOR = 1e3

const FACE_CUT = -1
const FACE_IN = 1
const FACE_OUT = 2

struct SubtriangulationLabels
  cell_to_bgcell::Vector{Int32}
  cell_to_io::Vector{Int8}
  face_to_stlface::Vector{Int32}
  face_to_bgcell::Vector{Int32}
  face_to_ios::Vector{Int8}
  bgcell_to_ioc::Vector{Int8}
  bgface_to_ioc::Vector{Int8}
end

function subtriangulate(
  bgmodel::DiscreteModel,
  stlmodel::DiscreteModel;
  threading=:spawn,
  kdtree=false,
  tolfactor=DEFAULT_TOL_FACTOR,
  surfacesource=:skin,
  showprogress=true)

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
  progress = progress_bar(cut_cells;showprogress)
  if threading == :threads || Threads.nthreads() == 1
    Threads.@threads for cell in cut_cells
      save_cell_submesh!(submesh,io_arrays,stl,p,cell,
        compute_polyhedra!(caches,Γ0,stl,p,f_to_isempty,Πf,Πr,
          c_to_stlf,node_to_coords,cell_to_nodes,cell;atol,kdtree)...
          ;surfacesource )
      next!(progress)
    end
  elseif threading == :spawn
    @sync for cell in cut_cells
      Threads.@spawn save_cell_submesh!(submesh,io_arrays,stl,p,cell,
        compute_polyhedra!(caches,Γ0,stl,p,f_to_isempty,Πf,Πr,
          c_to_stlf,node_to_coords,cell_to_nodes,cell;atol,kdtree)...
          ;surfacesource )
      next!(progress)
    end
  else
    @unreachable
  end

  submesh = _append_threaded_submesh!(submesh)
  io_arrays = _reduce_io_arrays(bgmodel,io_arrays)
  bgcell_to_ioc, bgnode_to_io, bgfacet_to_ioc = io_arrays
  T,F,X,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf, f_to_ios = submesh

  propagate_inout!(bgmodel,bgcell_to_ioc,bgnode_to_io)
  set_facets_as_inout!(bgmodel,bgcell_to_ioc,bgfacet_to_ioc)

  delete_small_subcells!(bgmodel,T,X,k_to_io,k_to_bgcell)
  delete_small_subfacets!(bgmodel,F,Xf,f_to_bgcell,f_to_stlf,f_to_ios)

  cell_grid = compute_grid(Table(T),X,TET)
  face_grid = compute_grid(Table(F),Xf,TRI)
  labels = SubtriangulationLabels(
    k_to_bgcell, k_to_io,
    f_to_stlf, f_to_bgcell, f_to_ios,
    bgcell_to_ioc, bgfacet_to_ioc)

  cell_grid, face_grid, labels
end

function subtriangulate(bgmodel::DiscreteModel,args...;kwargs...)
  stlmodel = _to_stl_model(args...)
  subtriangulate(bgmodel,stlmodel;kwargs...)
end

function _to_stl_model(filename::AbstractString)
  X,T,N = read_stl(filename)
  stl = compute_stl_model(T,X)
  stl = merge_and_collapse(stl)
  stl
end

function progress_bar(a::AbstractVector;kw...)
  progress_bar(length(a);kw...)
end

function progress_bar(n::Integer;showprogress)
  p = Progress(n,
    desc="STLCutter: ",
    dt = 0.5,
    barglyphs=BarGlyphs("[=> ]"),
    barlen=40,
    color=:none,
    enabled=showprogress)
end

function compute_polyhedra!(caches,Γ0,stl,p,f_to_isempty,Πf,Πr,
  c_to_stlf,node_to_coords,cell_to_nodes,cell;atol,kdtree)

  cell_coords = _get_cell_coordinates!(caches,node_to_coords,cell_to_nodes,cell)

   _faces = _get_cell_faces!(caches,stl,c_to_stlf,cell,f_to_isempty)
  facets,empty_facets,reflex_faces = _faces

  Πk,Πk_ids,Πk_io = get_cell_planes(p,cell_coords)
  Πk⁺,Πk⁺_ids,Πk⁺_io = expand_planes(Πk,Πk_io,atol), Πk_ids, Πk_io
  Πk⁺_ids = Πk⁺_ids .+ Πk_ids[end]

  Πkf,Πkf_ids = filter_face_planes(stl,Πr,reflex_faces,Πf,facets)

  K = Polyhedron(p,cell_coords)
  Γk0 = restrict(Γ0,stl,c_to_stlf[cell])

  compute_distances!(Γk0,lazy_append(Πk,Πkf),lazy_append(Πk_ids,Πkf_ids))
  compute_distances!(K,lazy_append(Πk,Πkf),lazy_append(Πk_ids,Πkf_ids))

  _merge_coplanar_planes!(Γk0,K,stl,Πk,Πr,Πf;atol)

  compute_distances!(Γk0,Πk⁺,Πk⁺_ids)
  Γk = clip(Γk0,Πk⁺_ids,inout=Πk⁺_io,boundary=true)

  if isnothing(Γk)
    return nothing,nothing,nothing
  end

  planes = get_active_planes(Γk,K,stl)

  restrict_planes!(Γk0,planes)
  restrict_planes!(K,planes)

  merge_coplanar_planes!(Γk0,K,stl,Πr,Πf;atol)

  compute_distances!(Γk0,Πk⁺,Πk⁺_ids)
  Γk = clip(Γk0,Πk⁺_ids,inout=Πk⁺_io,boundary=true)

  if isnothing(Γk) || isempty(get_original_facets(Γk,stl))
    nothing,nothing,nothing
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

    Γk = clip(Γk0,Πk_ids,inout=Πk_io)

    Kn_in,Kn_out,Γk
  end
end

function save_cell_submesh!(submesh,io_arrays,stl,p,cell,Kn_in,Kn_out,Γk;
  surfacesource)

  !isnothing(Kn_in) || return
  bgcell_to_ioc, bgcell_node_to_io, bgcell_facet_to_ioc = io_arrays
  Tin,Xin = simplexify(Kn_in)
  Tout,Xout = simplexify(Kn_out)
  B = _simplexify_boundary(Kn_in,Kn_out,Γk,stl;surfacesource)
  T_Γ,X_Γ,f_to_f,f_to_ios = B
  bgcell_to_ioc[cell] = _get_cell_io(T_Γ,Kn_in,Kn_out)
  bgcell_to_ioc[cell] == FACE_CUT || return
  D = num_dims(stl)
  f_to_f .-= get_offset(stl,D)
  n_to_io = get_cell_nodes_to_inout(Kn_in,Kn_out,p)
  f_to_ioc = get_cell_facets_to_inoutcut(Kn_in,Kn_out,p)

  f_to_iscut = _get_cell_facets_to_iscut(Γk,p;surfacesource)
  for (f,iscut) in enumerate(f_to_iscut)
    if iscut
      f_to_ioc[f] = FACE_CUT
    end
  end

  _append_submesh!(submesh,Xin,Tin,Xout,Tout,X_Γ,T_Γ,f_to_f,f_to_ios,cell)
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
        for node in getindex!(node_cache,c_to_n,current_cell)
          if bgnode_to_io[node] == FACE_IN ||
             ( bgcell_to_ioc[cell] == FACE_IN &&
               bgnode_to_io[node] == UNSET )

            for neig_cell in getindex!(neig_cell_cache,n_to_c,node)
              if bgcell_to_ioc[neig_cell] == UNSET
                bgcell_to_ioc[neig_cell] = FACE_IN
                push!(stack,neig_cell)
                for neig_node in getindex!(neig_node_cache,c_to_n,neig_cell)
                  if bgnode_to_io[neig_node] == UNSET
                    bgnode_to_io[neig_node] = FACE_IN
                  end
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


function delete_small_subcells!(bgmodel,T,X,arrays...)
  delete_small_subfaces!(bgmodel,T,X,TET,arrays...)
end

function delete_small_subfacets!(bgmodel,T,X,arrays...)
  delete_small_subfaces!(bgmodel,T,X,TRI,arrays...)
end

function max_length(model::CartesianDiscreteModel)
  float(get_cartesian_descriptor(model).sizes[1])
end

function delete_small_subfaces!(bgmodel,T,X,p::Polytope{D},arrays...) where D
  h = max_length(bgmodel)
  c = array_cache(T)
  ids = findall( i -> measure(get_cell!(c,T,X,p,i)) < eps(h^D), 1:length(T) )
  deleteat!(T,ids)
  for array in arrays
    deleteat!(array,ids)
  end
end

function refine(
  K::Polyhedron,
  Γ::Polyhedron,
  stl::STL,
  reflex_faces::AbstractVector,
  empty_facets::AbstractVector=[],
  ;inside::Bool)

  reflex_faces = filter( f -> is_reflex(Γ,stl,f;inside), reflex_faces )
  K = copy(K)
  Γn,Kn = decompose(Γ,K,reflex_faces,empty_facets,stl)
  Kn_clip = empty(Kn)
  for (i,(Γi,Ki)) in enumerate(zip(Γn,Kn))
    facets = get_original_facets(Γi,stl)
    _clear_empty_facets!(facets,empty_facets)
    facets = _add_missing_facets(Γi,stl,facets,reflex_faces,empty_facets)
    @assert !isempty(facets)
    part_to_facets = get_disconnected_facets(Γi,stl)
    colour_to_facets = colour_facing_facets(Γi,facets,part_to_facets;inside)
    for facets in colour_to_facets
      Ki_clip = clip(Ki,facets;inside)
      !isnothing(Ki_clip) || continue
      push!(Kn_clip,Ki_clip)
    end
  end
  Kn_clip
end

function decompose(surf::Polyhedron,cell::Polyhedron,rfaces,empty_facets,stl)
  i = findfirst(i->has_original_reflex_face(surf,i,empty=false), rfaces )
  if i === nothing
    R = [surf],[cell]
  else
    delete_inactive_planes!(surf,cell,stl)
    rf = rfaces[i]
    if !has_coplanars(surf.data,rf) || !contains_coplanars(surf.data,rf)
      s⁻,s⁺ = split_overlapping(surf,rf)
      k⁻,k⁺ = split(cell,rf)
      S = [s⁻,s⁺]
      K = [k⁻,k⁺]
      if any(isnothing,S) || any(!has_facets,S)
        S,K = _split_reflex_face(S,K,surf,cell,stl,rf,empty_facets)
      end
      if any(isnothing,K)
        j = findfirst(!isnothing,K)
        S,K = [S[j]],[K[j]]
      end
    else
      S,K = [surf],[cell]
    end
    if length(rfaces) == i
      R = S,K
    else
      Sr,Kr = empty(S),empty(K)
      rfaces = view(rfaces,i+1:length(rfaces))
      for (s,k) in zip(S,K)
        Si,Ki = decompose(s,k,rfaces,empty_facets,stl)
        append!(Sr,Si)
        append!(Kr,Ki)
      end
      R = Sr,Kr
    end
  end
  R
end


## Threading data

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
  f_to_ios = Int8[]
  T,F,X,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf,f_to_ios
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
  empty_facets = Int32[]
  reflex_faces = Int32[]
  nc,facets,empty_facets,reflex_faces
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
  _,facets,empty_facets,reflex_faces = caches[i]
  D = num_dims(stl)
  get_reflex_faces!(reflex_faces,stl,c_to_stlf,cell)
  copy!(facets,c_to_stlf[cell])
  copy!(empty_facets,c_to_stlf[cell])
  filter!(f->f_to_isempty[f],empty_facets)
  map!(i->i+get_offset(stl,D),facets,facets)
  map!(i->i+get_offset(stl,D),empty_facets,empty_facets)
  map!(i->i+get_offset(stl,D-1),reflex_faces,reflex_faces)
  facets,empty_facets,reflex_faces
end

function _get_cell_io(T_Γ,Kn_in,Kn_out)
  if !isnothing(T_Γ) && length(T_Γ) > 0
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

function _append_submesh!(submesh,Xin,Tin,Xout,Tout,Xfn,Tfn,fn_to_f,fn_to_ios,bgcell)
  i = Threads.threadid()
  _submesh = submesh[i]
  _append_subcells!(_submesh,Xin,Tin,FACE_IN,bgcell)
  _append_subcells!(_submesh,Xout,Tout,FACE_OUT,bgcell)
  _append_subfacets!(_submesh,Xfn,Tfn,fn_to_f,fn_to_ios,bgcell)
end

function _append_subcells!(submesh_arrays,Xn,Tn,io,bgcell)
  T,F,X,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf = submesh_arrays
  append!(T, map(i->i.+length(X),Tn) )
  append!(X,Xn)
  append!(k_to_io,fill(io,length(Tn)))
  append!(k_to_bgcell,fill(bgcell,length(Tn)))
  submesh_arrays
end

function _append_subfacets!(submesh_arrays,Xfn,Tfn,fn_to_f,fn_to_ios,bgcell)
  T,F,X,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf,f_to_ios = submesh_arrays
  append!(F, map(i->i.+length(Xf),Tfn) )
  append!(Xf,Xfn)
  append!(f_to_bgcell,fill(bgcell,length(Tfn)))
  append!(f_to_stlf,fn_to_f)
  append!(f_to_ios,fn_to_ios)
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

function _simplexify_boundary(Kn_in,Kn_out,Γk,stl;surfacesource)
  if surfacesource == :both
    S = (Kn_in,Kn_out,Γk)
    l = (FACE_IN,FACE_OUT,FACE_CUT)
    T_Γ,X_Γ,f_to_f,f_to_ios = simplexify_boundary(S,l,stl)
  else
    if surfacesource == :skin
      S = Γk
      ios = FACE_CUT
    elseif surfacesource == :inside
      S = Kn_in
      ios = FACE_IN
    elseif surfacesource == :outside
      S = Kn_out
      ios = FACE_IN
    else
      @notimplemented "surfacesource=$surfacesource is not implemented"
    end
    T_Γ,X_Γ,f_to_f = simplexify_boundary(S,stl)
    f_to_ios = isnothing(T_Γ) ? nothing : fill(Int8(ios),length(T_Γ))
  end
  T_Γ,X_Γ,f_to_f,f_to_ios
end

function _get_cell_facets_to_iscut(Γk,p;surfacesource)
  f_to_iscut = falses( num_facets(p) )
  if surfacesource == :skin
    v_to_Π = Γk.data.vertex_to_planes
    for v in 1:length(v_to_Π)
      for Π in v_to_Π[v]
        if Π < 0 && 0 < abs(Π) ≤ num_facets(p)
          f_to_iscut[ abs(Π) ] = true
        end
      end
    end
  end
  f_to_iscut
end


## Face/Reflex planes

function get_facet_to_isempty(stl::STL;atol)
  f_to_isempty = falses(num_cells(stl))
  c = get_cell_cache(stl)
  for f in 1:num_cells(stl)
    facet = get_cell!(c,stl,f)
    if min_height(facet) < atol
      f_to_isempty[f] = true
    end
  end
  f_to_isempty
end

function get_facet_planes(stl::STL)
  c = get_cell_cache(stl)
  [ get_facet_plane!(c,stl,i) for i in 1:num_cells(stl) ]
end


function get_facet_plane!(c,stl::STL,i)
  f = get_cell!(c,stl,i)
  Plane(f)
end

function correct_small_facets_planes!(stl::STL,Πf,f_to_isempty;atol)
  D = num_dims(stl)
  e_to_f = get_faces(stl,D-1,D)
  f_to_e = get_faces(stl,D,D-1)
  full_facets = Int32[]
  queue = Int32[]
  num_new_planes = 0
  Πnew = empty(Πf)
  f_to_new_plane = fill(UNSET,num_cells(stl))
  ec = array_cache(f_to_e)
  fc = array_cache(e_to_f)
  for f in 1:num_cells(stl)
    f_to_isempty[f] || continue
    f_to_new_plane[f] == UNSET || continue
    num_new_planes += 1
    f_to_new_plane[f] = num_new_planes
    empty!(full_facets)
    empty!(queue)
    push!(queue,f)
    head = 1
    while length(queue) ≥ head
      f_curr = queue[head]
      head += 1
      for e in getindex!(ec,f_to_e,f_curr)
        for f_neig in getindex!(fc,e_to_f,e)
          f_neig ≠ f_curr || continue
          if f_to_isempty[f_neig]
            f_to_new_plane[f_neig] == UNSET || continue
            f_to_new_plane[f_neig] = num_new_planes
            push!(queue,f_neig)
          elseif f_to_new_plane[f_neig] ≠ num_new_planes
            f_to_new_plane[f_neig] == num_new_planes
            push!(full_facets,f_neig)
          end
        end
      end
    end
    empty_facets = queue
    Π = _get_mean_plane(Πf,empty_facets,full_facets)
    push!(Πnew,Π)
  end
  for f in 1:num_cells(stl)
    f_to_isempty[f] || continue
    Πf[f] = Πnew[ f_to_new_plane[f] ]
  end
  _check_planes(Πf,f_to_isempty,stl;atol)
  Πf
end

function _get_mean_plane(Πf,empty_facets,full_facets)
  !isempty(full_facets) || @unreachable
  n = zero( normal(first(Πf)) )
  c = zero( center(first(Πf)) )
  for full_f in full_facets
    n += normal(Πf[full_f]) / length(full_facets)
  end
  for empty_f in empty_facets
    c += center(Πf[empty_f]) / length(empty_facets)
  end
  Plane(c,n)
end

function _check_planes(Πf,f_to_isempty,stl;atol)
  D = num_dims(stl)
  v_coords = get_vertex_coordinates(stl)
  f_to_v = get_face_vertices(stl,D)
  c = array_cache(f_to_v)
  for f in 1:num_cells(stl)
    f_to_isempty[f] || continue
    for i in getindex!(c,f_to_v,f)
      v = v_coords[i]
      @assert abs(signed_distance(v,Πf[f])) < atol
    end
  end
end

function get_reflex_planes(stl::STL,Πf)
  Dc = num_dims(stl)
  Dp = num_point_dims(stl)
  c = bisector_plane_cache(stl,Dc-1)
  T_Π = typeof( bisector_plane!(c,stl,Dc-1,1,Πf) )
  Π_r = Vector{T_Π}(undef,num_faces(stl,Dc-1))
  for f in 1:num_faces(stl,Dc-1)
    Π_r[f] = bisector_plane!(c,stl,Dc-1,f,Πf)
  end
  Π_r
end

function colour_facing_facets(poly::Polyhedron,facets,part_to_facets;inside)
  length(part_to_facets) > 1 || return [facets]
  f_to_part, parts = _get_facet_to_part(facets,part_to_facets)
  length(parts) > 1 || return [facets]
  p_to_p_to_facing = _get_part_to_part_to_facing(poly,facets,f_to_part,parts;inside)
  p_to_colour = fill(UNSET,length(parts))
  colour = Int32[]
  ids = Int32[]
  num_colours = 0
  for p in 1:length(parts)
    p_to_colour[p] == UNSET || continue
    num_colours += 1
    empty!(colour)
    push!(colour,p)
    for i in 1:length(parts)
      if p_to_p_to_facing[p][i] && p_to_p_to_facing[i][p] && i ∉ colour
        push!(colour,i)
      end
    end
    empty!(ids)
    for (i,p_i) in enumerate(colour)
      for p_j in colour
        if !p_to_p_to_facing[p_j][p_i]
          push!(ids,i)
          break
        end
      end
    end
    deleteat!(colour,ids)
    for p_i in colour
      p_to_colour[p_i] == UNSET || error("
      Ambiguity in colouring, check whether the STL geometry is a polyhedron, i.e., is solid and manifold. I.e., no self-intersections, no degenerate faces, water tight, manifold and no sharp edges.
      Run `check_requisites(geo,bgmodel)` (not checking self-intersections currently)")
      p_to_colour[p_i] = num_colours
    end
  end
  colour_to_facets = [ Int32[] for _ in 1:num_colours ]
  for (i,f) in enumerate(facets)
    g = p_to_colour[f_to_part[i]]
    push!(colour_to_facets[g],f)
  end
  colour_to_facets
end

function _get_facet_to_part(facets,part_to_facets)
  f_to_part = fill(UNSET,length(facets))
  for (i,f) in enumerate(facets)
    f_to_part[i] = findfirst( p -> f ∈ p, part_to_facets )
  end
  parts = unique(f_to_part)
  length(parts) > 1 || return f_to_part,parts
  for i in 1:length(facets)
    p = findfirst(isequal(f_to_part[i]),parts)
    f_to_part[i] = p
  end
  f_to_part,parts
end

function _get_part_to_part_to_facing(poly,facets,f_to_part,parts;inside)
  p_to_p_to_facing = [ trues(length(parts)) for _ in 1:length(parts) ]
  for (i,fi) in enumerate(facets), (j,fj) in enumerate(facets)
    fi ≠ fj || continue
    f_to_part[i] ≠ f_to_part[j] || continue
    if !is_facet_in_facet(poly,fj,fi;inside)
      p_to_p_to_facing[f_to_part[i]][f_to_part[j]] = false
    end
  end
  p_to_p_to_facing
end

function get_disconnected_facets(poly::Polyhedron,stl::STL)
  v_to_pv = poly.data.vertex_to_parent_vertex
  v_to_f = poly.data.vertex_to_original_faces
  facets = get_original_facets(poly,stl,empty=true)
  facet_to_part = fill(UNSET,length(facets))
  stack = Int32[]
  num_parts = 0
  for (i,facet) in enumerate(facets)
    facet_to_part[i] == UNSET || continue
    num_parts += 1
    facet_to_part[i] = num_parts
    empty!(stack)
    push!(stack,facet)
    while !isempty(stack)
      current_facet = pop!(stack)
      v,vnext = get_facet_onset(poly,current_facet) # reduce calls
      vcurrent = v
      while true
        if v_to_pv[vcurrent] ≠ v_to_pv[vnext]
          for fi in v_to_f[vcurrent]
            fi ≠ current_facet || continue
            for fj in v_to_f[vnext]
              fj ≠ current_facet || continue
              if fi == fj
                fi ∈ facets || continue
                _i = findfirst(isequal(fi),facets)
                facet_to_part[_i] == UNSET || continue
                facet_to_part[_i] = num_parts
                push!(stack,fi)
              end
            end
          end
        end
        vcurrent,vnext = vnext,next_vertex(poly,vcurrent,vnext)
        vcurrent ≠ v || break
      end
    end
  end
  part_to_facets = [ Int32[] for _ in 1:num_parts ]
  for (i,facet) in enumerate(facets)
    push!(part_to_facets[facet_to_part[i]],facet)
  end
  part_to_facets
end

function get_facet_onset(poly::Polyhedron,facet::Integer)
  v_to_pv = poly.data.vertex_to_parent_vertex
  v_to_f = poly.data.vertex_to_original_faces
  for v in 1:num_vertices(poly)
    if isactive(poly,v) && facet ∈ v_to_f[v]
      for vneig in get_graph(poly)[v]
        vneig ∉ (OPEN,UNSET) || continue
        if facet ∈ v_to_f[vneig]
          vcurrent = v
          vnext = vneig
          while vnext ≠ v
            vcurrent,vnext = vnext,next_vertex(poly,vcurrent,vnext)
            vnext ∉ (UNSET,OPEN) || break
            facet ∈ v_to_f[vnext] || break
          end
          if vnext == v
            return v,vneig
          end
        end
      end
    end
  end
end

function is_facet_in_facet(poly::Polyhedron,facet,plane;inside,atol=0)
  if has_coplanars(poly.data,plane)
    if are_coplanar(poly.data,facet,plane)
      return true
    end
  end
  v_to_f = get_data(poly).vertex_to_original_faces
  distances = get_plane_distances(get_data(poly),plane)
  smax = -Inf
  smin = Inf
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    if facet ∈ v_to_f[v] && plane ∉ v_to_f[v]
      smin = min(smin,distances[v])
      smax = max(smax,distances[v])
    end
  end
 # @assert !(smin == smax == 0)
  ( smin ≥ 0 && smax ≥ atol ) && return !inside
  ( smin ≤ 0 && smax ≤ -atol ) && return inside
  false
end

function is_reflex(poly::Polyhedron,stl::STL,reflex_face;inside)
  Dc = num_dims(stl)
  rf = reflex_face - get_offset(stl,Dc-1)
  length(get_faces(stl,Dc-1,Dc)[rf]) == 2 || return false
  f1 = get_faces(stl,Dc-1,Dc)[rf][1] + get_offset(stl,Dc)
  f2 = get_faces(stl,Dc-1,Dc)[rf][2] + get_offset(stl,Dc)
  has_plane(poly.data,f1) || return false
  has_plane(poly.data,f2) || return false
  is_facet_in_facet(poly,f1,f2;inside) || return true
  is_facet_in_facet(poly,f2,f1;inside) || return true
  false
end

function filter_face_planes(
  stl::STL,
  reflex_planes::AbstractVector,
  reflex_faces::AbstractVector{<:Integer},
  facet_planes::AbstractVector,
  facets::AbstractVector{<:Integer})

  Dc = num_dims(stl)
  Πr = view(reflex_planes,lazy_map(i->i-get_offset(stl,Dc-1),reflex_faces))
  Πf = view(facet_planes,lazy_map(i->i-get_offset(stl,Dc),facets))
  lazy_append(Πr,Πf),lazy_append(reflex_faces,facets)
end

## Postprocess

function get_cell_nodes_to_inout(polys_in,polys_out,p::Polytope)
  node_to_inout = fill(UNSET,num_vertices(p))

  complete_nodes_to_inout!(node_to_inout,polys_in,FACE_IN,p)
  complete_nodes_to_inout!(node_to_inout,polys_out,FACE_OUT,p)
  #@assert UNSET ∉ node_to_inout
  node_to_inout
end

function get_cell_facets_to_inoutcut(polys_in,polys_out,p::Polytope)
  facet_to_inoutcut = fill(UNSET,num_facets(p))

  complete_facets_to_inoutcut!(facet_to_inoutcut,polys_in,FACE_IN,p)
  complete_facets_to_inoutcut!(facet_to_inoutcut,polys_out,FACE_OUT,p)
  @assert UNSET ∉ facet_to_inoutcut
  facet_to_inoutcut
end

function complete_nodes_to_inout!(node_to_inout,polys,inout,p::Polytope)
  D = num_dims(p)
  for poly in polys
    v_to_Π = poly.data.vertex_to_planes
    v_to_v = poly.data.vertex_to_parent_vertex
    for v in 1:num_vertices(poly)
      isactive(poly,v) || continue
    #  v = inout == FACE_CUT ? v : v_to_v[v]
      node = find_polytope_vertex(v_to_Π[v],p)
      if node != UNSET
      #@assert  node_to_inout[node] ∈ (inout,UNSET)
        _inout = inout
        if node_to_inout[node] ∉ (inout,UNSET)
          _inout = FACE_CUT
        end
        node_to_inout[node] = _inout
      end
    end
  end
  node_to_inout
end

function find_polytope_vertex(planes,p::Polytope)
  if is_n_cube(p)
    _find_n_cube_vertex(planes,p)
  else
    _find_polytope_vertex(planes,p)
  end
end

function _find_polytope_vertex(planes,p::Polytope)
  D = num_dims(p)
  n_facets = num_faces(p,D-1)
  if count(Π -> -n_facets ≤ Π < 0 ,planes) == D
    f_to_v = get_faces(p,D-1,0)
    v_to_f = get_faces(p,0,D-1)
    i = findfirst(Π -> Π < 0,planes)
    f = -planes[i]
    for v in f_to_v[f]
      if all(f->(-f) ∈ planes,v_to_f[v])
        return v
      end
    end
  end
  UNSET
end


function _find_n_cube_vertex(planes,p::Polytope)
  @notimplementedif !is_n_cube(p)
  D = num_dims(p)
  if count(Π -> -2*D ≤ Π < 0 , planes) == D
    i = 0
    node = 0
    for _ in 1:D
      i = findnext(Π -> Π < 0, planes,i+1)
      f = -planes[i]
      d = D - ((f-1)>>1)
      ud = iseven(f)
      node |= ud<<(d-1)
    end
    node += 1
  else
    node = UNSET
  end
  node
end

function complete_facets_to_inoutcut!(facet_to_inoutcut,polys,inout,p::Polytope)
  facet_list = Int32[]
  for poly in polys
    istouch = map( i -> falses(length(i)), get_graph(poly) )
    v_to_Π = poly.data.vertex_to_planes
    for v in 1:num_vertices(poly)
      isactive(poly,v) || continue
      for i in 1:length(get_graph(poly)[v])
        !istouch[v][i] || continue
        istouch[v][i] = true
        vcurrent = v
        vnext = get_graph(poly)[v][i]
        vnext ∉ (OPEN,UNSET) || continue
        any( i -> -num_facets(p) ≤ i < 0, v_to_Π[v] ) || continue
        empty!(facet_list)
        append!( facet_list, v_to_Π[v] )
        filter!( i -> -num_facets(p) ≤ i < 0, facet_list)
        filter!( i -> i ∈ v_to_Π[vnext], facet_list )
        while vnext != v
          inext = findfirst( isequal(vcurrent), get_graph(poly)[vnext] )
          inext = ( inext % length( get_graph(poly)[vnext] ) ) + 1
          istouch[vnext][inext] = true
          vcurrent = vnext
          vnext = get_graph(poly)[vnext][inext]
          vnext ∉ (OPEN,UNSET) || break
          filter!( i -> i ∈ v_to_Π[vnext], facet_list )
          !isempty(facet_list) || break
          if vnext == v
            f = abs( only( facet_list) )
            if facet_to_inoutcut[f] == UNSET
              facet_to_inoutcut[f] = inout
            elseif facet_to_inoutcut[f] ≠ inout
              facet_to_inoutcut[f] = FACE_CUT
            end
          end
        end
      end
    end
  end
  facet_to_inoutcut
end

## (Quasi-)Coplanar stuff

function _merge_coplanar_planes!(Γk0,K,stl,Πk,Πr,Πf;atol)
  snap_distances!(Γk0,K;atol=atol/10)
  Π_to_qp_Π = merge_coplanar_with_cell_planes(Γk0,K,stl,Πk,Πr,Πf;atol)
  correct_plane_distances!(Γk0,Π_to_qp_Π,stl;atol)
  set_linked_planes!(Γk0,Π_to_qp_Π)
  set_linked_planes!(K,Π_to_qp_Π)
end

function snap_distances!(P,Q;atol)
  snap_distances!(P;atol)
  snap_distances!(Q;atol)
end

function snap_distances!(P;atol)
  Π_to_v_to_dist = get_plane_distances(P.data)
  for i in 1:length(Π_to_v_to_dist)
    v_to_dist = Π_to_v_to_dist[i]
    for v in 1:length(v_to_dist)
      if abs(v_to_dist[v]) < atol
        v_to_dist[v] = 0.0
      end
    end
  end
end

function merge_coplanar_with_cell_planes(S,K,stl,Πk,Πr,Πf;atol)
  planes = S.data.plane_to_ids
  plane_to_qp_plane = fill(Int32(UNSET),length(planes))
  for (i,Πi) in enumerate(planes)
    Πi < 0 || continue
    for (j,Πj) in enumerate(planes)
      Πj > 0 || continue
      if are_quasi_coplanar(K,S,Πi,Πj,stl;atol,restrict=true)
        s = relative_orientation(Πi,Πj,Πk,Πr,Πf,stl;atol)
        plane_to_qp_plane[j] = i * s
        plane_to_qp_plane[i] = i
      end
    end
  end
  plane_to_qp_plane
end

function relative_orientation(i,j,Πk,Πr,Πf,stl;atol)
  Πi = _get_plane(i,Πk,Πr,Πf,stl)
  Πj = _get_plane(j,Πk,Πr,Πf,stl)
  relative_orientation(Πi,Πj;atol)
end

function _get_plane(fi,Πk,Πr,Πf,stl)
  D = num_dims(stl)
  facedims = get_facedims(stl)
  rf_offset = get_offset(stl,D-1)
  f_offset = get_offset(stl,D)
  if fi < 0
    Πi = Πk[ abs(fi) ]
  else
    Πi = facedims[fi] == D ? Πf[fi-f_offset] : Πr[fi-rf_offset]
  end
  Πi
end


function correct_plane_distances!(S,plane_to_qp_plane,stl;atol)
  plane = S.data.plane_to_ids
  Π_to_v_to_dist = S.data.plane_to_vertex_to_distances
  v_to_f = S.data.vertex_to_original_faces
  for (i,qp_i) in enumerate( plane_to_qp_plane )
    qp_i ≠ 0 || continue
    v_to_dist = Π_to_v_to_dist[ abs(qp_i) ]
    for v in 1:length(v_to_dist)
      if plane[i] ∈ v_to_f[v]
        v_to_dist[v] = 0.0
      end
    end
  end
end

function set_linked_planes!(poly::Polyhedron,Π_to_qp_Π)
  planes = poly.data.plane_to_ids
  set_linked_planes!(poly,Π_to_qp_Π,planes)
end

function are_quasi_coplanar(K,S,Πi,Πj,stl;atol,restrict=false)
  distance_between_planes(S,Πi,Πj,stl;restrict) < atol && return true
  distance_between_planes(K,Πi,Πj,stl;restrict) < atol && return true
  false
end

function merge_coplanar_planes!(Γk0,K,stl,Πr,Πf;atol)
  Π_to_refΠ,Πs = link_planes!(Γk0,K,stl,Πr,Πf;atol)
  set_linked_planes!(Γk0,Π_to_refΠ,Πs)
  set_linked_planes!(K,Π_to_refΠ,Πs)
end

function link_planes!(surf::Polyhedron,cell::Polyhedron,stl::STL,Πr,Πf;atol)
  planes = surf.data.plane_to_ids
  default_out = fill(UNSET,length(planes)), planes
  Π_to_faces = find_faces_on_planes!(surf,stl;atol)
  maximum(length,Π_to_faces) > 0 || return default_out
  Π_to_coplanar_Π = link_coplanar_planes(surf,cell,stl,Π_to_faces;atol)
  maximum(length,Π_to_coplanar_Π) > 0 || return default_out
  Π_to_ref_Π =  group_coplanar_planes(planes,Π_to_coplanar_Π)
  merge_faces!(Π_to_faces,planes,Π_to_ref_Π)
  correct_distances!(surf,Π_to_ref_Π,Π_to_faces)
  Π_to_ref_Π = mark_inverted_planes!(Π_to_ref_Π,planes,Πr,Πf,stl;atol)

  Π_to_ref_Π, planes
end

function find_faces_on_planes!(surf::Polyhedron,stl::STL;atol)
  planes = surf.data.plane_to_ids
  v_to_f = surf.data.vertex_to_original_faces
  f_to_v = get_face_vertices(stl)
  c = array_cache(f_to_v)
  vertices = Int32[]
  Π_to_faces = [ Int32[] for _ in 1:length(planes) ]
  for (i,Π) in enumerate(planes)
    Π > 0 || continue
    faces = Π_to_faces[i]
    dists = get_plane_distances(surf.data,Π)
    empty!(vertices)
    for v in 1:num_vertices(surf)
      @assert isactive(surf,v)
      if abs(dists[v]) < atol/10
        push!(vertices,v)
        dists[v] = 0
      end
    end
    for v in vertices
      for f in v_to_f[v]
        Π ≠ f || continue
        f ∉ faces || continue
        f ∈ planes || continue
        face_on_plane = true
        for _v in getindex!(c,f_to_v,f)
          if !any( i-> first(v_to_f[i]) == _v, vertices )
            face_on_plane = false
            break
          end
        end
        if face_on_plane
          push!(faces,f)
        end
      end
    end
  end
  Π_to_faces
end

function link_coplanar_planes(surf,cell,stl,Π_to_faces;atol)
  planes = surf.data.plane_to_ids
  Π_to_ref_Π = surf.data.plane_to_ref_plane
  Π_to_coplanar_Π = [ Int32[] for _ in 1:length(planes) ]
  D = num_dims(stl)
  facedims = get_facedims(stl)
  rf_offset = get_offset(stl,D-1)
  f_offset = get_offset(stl,D)
  for (i,Π) in enumerate(planes)
    Π > 0 || continue
    Π_to_ref_Π[i] == UNSET || continue
    faces = Π_to_faces[i]
    for f in faces
      d = facedims[f]
      j = findfirst(isequal(f),planes)
      Π_to_ref_Π[j] == UNSET || continue
      if d == D-1
        if Π ∈ Π_to_faces[j] && are_quasi_coplanar(surf,cell,Π,f,stl;atol)
          push!(Π_to_coplanar_Π[i],f)
        end
      elseif d == D
        facet = get_cell(stl,f-get_offset(stl,D))
        if are_quasi_coplanar(surf,cell,Π,f,stl;atol)
          push!(Π_to_coplanar_Π[i],f)
        end
      end
    end
  end
  for (i,Πi) in enumerate(planes)
    Πi > 0 || continue
    for Πj in Π_to_coplanar_Π[i]
      j = findfirst(isequal(Πj),planes)
      if Πi ∉ Π_to_coplanar_Π[j]
        push!(Π_to_coplanar_Π[j],Πi)
      end
    end
  end
  Π_to_coplanar_Π
end

function group_coplanar_planes(planes,Π_to_coplanar_Π)
  Π_to_ref_Π = collect(1:length(planes))
  stack = Int32[]
  for (i,Πi) in enumerate(planes)
    Πi > 0 || continue
    Π_to_ref_Π[i] == i || continue
    empty!(stack)
    push!(stack,i)
    while !isempty(stack)
      j = pop!(stack)
      for Πk in Π_to_coplanar_Π[j]
        k = findfirst(isequal(Πk),planes)
        i ≠ k || continue
        Π_to_ref_Π[k] == k || continue
        Π_to_ref_Π[k] = i
        push!(stack,k)
      end
    end
  end
  for i in reverse(1:length(Π_to_ref_Π))
    if Π_to_ref_Π[i] > 0 && Π_to_ref_Π[i] ≠ i
      Π_to_ref_Π[ Π_to_ref_Π[i] ] = -Π_to_ref_Π[i]
      Π_to_ref_Π[i] = -Π_to_ref_Π[i]
    end
  end
  for i in 1:length(Π_to_ref_Π)
    if Π_to_ref_Π[i] > 0
      Π_to_ref_Π[i] = UNSET
    else
      Π_to_ref_Π[i] = -Π_to_ref_Π[i]
    end
  end
  Π_to_ref_Π
end

function merge_faces!(Π_to_faces,planes,Π_to_ref_Π)
  for (i,Πi) in enumerate(planes)
    Πi > 0 || continue
    Π_to_ref_Π[i] ≠ i || continue
    Π_to_ref_Π[i] ≠ UNSET || continue
    for f in Π_to_faces[i]
      if f ∉ Π_to_faces[Π_to_ref_Π[i]]
        push!(Π_to_faces[Π_to_ref_Π[i]],f)
      end
    end
  end
end

function correct_distances!(surf,Π_to_ref_Π,Π_to_faces)
  planes = surf.data.plane_to_ids
  v_to_f = surf.data.vertex_to_original_faces
  for (i,Πi) in enumerate(planes)
    Πi > 0 || continue
    Π_to_ref_Π[i] == i || continue
    !isempty(Π_to_faces[i]) || continue
    dists = get_plane_distances(surf.data,Πi)
    for v in 1:num_vertices(surf)
      dists[v] ≠ 0 || continue
      if any( f-> f ∈ v_to_f[v], Π_to_faces[i] )
        dists[v] = 0
      end
    end
  end
end

function mark_inverted_planes!(Π_to_ref_Π,planes,Πr,Πf,stl;atol)
  D = num_dims(stl)
  facedims = get_facedims(stl)
  rf_offset = get_offset(stl,D-1)
  f_offset = get_offset(stl,D)
  for i in 1:length(Π_to_ref_Π)
    Π_to_ref_Π[i] ∉ (UNSET,i) || continue
    j = Π_to_ref_Π[i]
    fi = planes[i]
    fj = planes[j]
    Πi = facedims[fi] == D ? Πf[fi-f_offset] : Πr[fi-rf_offset]
    Πj = facedims[fj] == D ? Πf[fj-f_offset] : Πr[fj-rf_offset]
    if relative_orientation(Πi,Πj;atol) < 0
      Π_to_ref_Π[i] = -Π_to_ref_Π[i]
    end
  end
  Π_to_ref_Π
end

function relative_orientation(Π1,Π2;atol)
  d = normal(Π1) ⋅ normal(Π2)
  @assert abs(d) > atol
  sign(d)
end

function distance_between_planes(poly::Polyhedron,Π1,Π2,stl;restrict=false)
  v_to_Π = poly.data.vertex_to_planes
  v_to_of = poly.data.vertex_to_original_faces
  if restrict
    @assert Π1 < 0
    @assert Π2 > 0
    if get_facedims(stl)[Π2] == 1
      e = Π2 - get_offset(stl,1)
      f1,f2 = get_faces(stl,1,num_dims(stl))[e]
      f1,f2 = (f1,f2) .+ get_offset(stl,num_dims(stl))
    else
      f1,f2 = Π2,Π2
    end
  end
  dist1 = get_plane_distances(poly.data,Π1)
  dist2 = get_plane_distances(poly.data,Π2)
  max_dist_a = 0.0
  for v in 1:num_vertices(poly)
    if !restrict || any( Π -> Π ∈ v_to_Π[v] || Π ∈ v_to_of[v], (Π1,f1,f2) )
      _d = abs(dist1[v] - dist2[v])
      max_dist_a = max(max_dist_a,_d)
    end
  end
  max_dist_b = 0.0
  for v in 1:num_vertices(poly)
    if !restrict || any( Π -> Π ∈ v_to_Π[v] || Π ∈ v_to_of[v], (Π1,f1,f2) )
      _d = abs(dist1[v] + dist2[v])
      max_dist_b = max(max_dist_b,_d)
    end
  end
  min( max_dist_a, max_dist_b)
end

function set_linked_planes!(poly::Polyhedron,Π_to_ref_Π,planes)
  _Π_to_ref_Π = poly.data.plane_to_ref_plane
  _planes = poly.data.plane_to_ids
  Π_to_v_to_dist = poly.data.plane_to_vertex_to_distances
  for (i,Π) in enumerate(planes)
    Π_to_ref_Π[i] ≠ UNSET || continue
    ref_Π = planes[abs(Π_to_ref_Π[i])]
    j = findfirst(isequal(Π),_planes)
    ref_j = findfirst(isequal(ref_Π),_planes)
    _Π_to_ref_Π[j] = ref_j * sign(Π_to_ref_Π[i])

    for v in 1:length(Π_to_v_to_dist[j])
      Π_to_v_to_dist[j][v] =  Π_to_v_to_dist[ref_j][v] * sign(Π_to_ref_Π[i])
    end
  end
end

function has_coplanars(data::PolyhedronData,Π)
  i = findfirst(isequal(Π),get_plane_ids(data))
  data.plane_to_ref_plane[i] ≠ UNSET
end

function are_coplanar(data::PolyhedronData,Πi,Πj)
  i = findfirst(isequal(Πi),get_plane_ids(data))
  j = findfirst(isequal(Πj),get_plane_ids(data))
  data.plane_to_ref_plane[i] == data.plane_to_ref_plane[j]
end

function add_plane!(data::PolyhedronData,Π)
  i = findfirst(isequal(Π),get_plane_ids(data))
  @assert data.plane_to_ref_plane[i] ≠ UNSET
  Πref = get_plane_ids(data)[abs(data.plane_to_ref_plane[i])]
  Πlast = UNSET
  for j in reverse(1:length(get_plane_ids(data)))
    jref = abs(data.plane_to_ref_plane[j])
    jref ≠ UNSET || continue
    if Πref == get_plane_ids(data)[jref]
      Πj = get_plane_ids(data)[j]
      for v in 1:length(data.vertex_to_planes)
        if Πj in data.vertex_to_planes[v]
          Πlast = Πj
          break
        end
      end
      Πlast == UNSET || break
    end
  end
  Πlast ≠ UNSET || return
  for v in 1:length(data.vertex_to_planes)
    if Πlast in data.vertex_to_planes[v]
      push!(data.vertex_to_planes[v],Π)
    end
  end
end

function contains_coplanars(data::PolyhedronData,Π)
  i = findfirst(isequal(Π),get_plane_ids(data))

  for (j,Πj) in enumerate(get_plane_ids(data))
    i ≠ j || continue
    if abs(data.plane_to_ref_plane[i]) == abs(data.plane_to_ref_plane[j])
      for (v,planes) in enumerate(data.vertex_to_planes)
        if Πj ∈ planes
          return true
        end
      end
    end
  end
  false
end


## Helpers

function _clear_empty_facets!(facets::Vector,empty_facets::Vector)
  if length(facets) > 1
    ids = findall(in(empty_facets),facets)
    if length(ids) < length(facets)
      deleteat!(facets,ids)
    end
  end
  facets
end

function _add_missing_facets(
  surf::Polyhedron,
  stl::STL,
  facets,
  reflex_faces,
  empty_facets)

  D = num_dims(stl)
  rf_offset = get_offset(stl,D-1)
  f_offset = get_offset(stl,D)
  for f in facets
    has_coplanars(surf.data,f) || continue
    contains_coplanars(surf.data,f) || continue
    for rf in get_faces(stl,D,D-1)[f-f_offset]
      rf += rf_offset
      rf ∉ reflex_faces || continue
      has_original_reflex_face(surf,rf,empty=false) || continue
      for neig_f in get_faces(stl,D-1,D)[rf-rf_offset]
        neig_f += f_offset
        neig_f ∉ facets || continue
        if has_original_facet(surf,neig_f,empty=true)
          if neig_f ∈ empty_facets
            for e in get_faces(stl,D,D-1)[neig_f-f_offset]
              e += rf_offset
              e ≠ rf || continue
              for neig_neig_f in get_faces(stl,D-1,D)[e-rf_offset]
                neig_neig_f += f_offset
                neig_neig_f ≠ neig_f || continue
                neig_neig_f ∉ facets || continue
                if has_original_facet(surf,neig_neig_f,empty=true)
                  @notimplementedif neig_neig_f ∈ empty_facets
                  push!(facets,neig_neig_f)
                end
              end
            end
          else
            push!(facets,neig_f)
          end
        end
      end
    end
  end
  facets
end

function _split_reflex_face(
  S,K,
  surf::Polyhedron,
  cell::Polyhedron,
  stl::STL,
  reflex_face::Integer,
  empty_facets::AbstractVector)

  D = num_dims(stl)
  rface_offset = get_offset(stl,D-1)
  facet_offset = get_offset(stl,D)
  neig_facets = get_faces(stl,D-1,D)[reflex_face-rface_offset]
  if !any(i->has_original_facet(surf,i+facet_offset),neig_facets) ||
     !has_original_reflex_face(surf,reflex_face,empty=false)

    Sr,Kr = [surf],[cell]
  else
    @notimplementedif count(in(empty_facets),neig_facets) > 1
    cond =
      f->!has_original_facet(surf,f+facet_offset) ||
        f+facet_offset∈empty_facets
    j = findfirst(cond,neig_facets)
    !isnothing(j) || error("Degenerated facet or sharp edge")
    missing_facet = neig_facets[j] + facet_offset
    _surf = _one_face_polyhedron(surf,missing_facet)
    j = findfirst(i->isnothing(i)||!has_facets(i),S)
    Sr,Kr = [surf,surf],K
    Sr[j] = _surf
  end
  Sr,Kr
end

function _one_face_polyhedron(poly::Polyhedron,face::Integer)
  v_to_f = get_data(poly).vertex_to_original_faces
  nodes = Int32[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    if face ∈ v_to_f[v]
      push!(nodes,v)
    end
  end
  sort!(nodes)
  r = restrict(poly,nodes)
  for v in 1:num_vertices(r)
    r.data.vertex_to_original_faces[v] = [face]
  end
  r
end
