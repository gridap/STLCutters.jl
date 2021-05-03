
Rx(ϕ) = TensorValue(
  1,0,0,
  0,cos(ϕ),-sin(ϕ),
  0,sin(ϕ),cos(ϕ))

Ry(θ) = TensorValue(
  cos(θ),0,sin(θ),
  0,1,0,
  -sin(θ),0,cos(θ))

Rz(ψ) = TensorValue(
  cos(ψ),sin(ψ),0,
  -sin(ψ),cos(ψ),0,
  0,0,1)

R(ϕ,θ,ψ) = Rx(ϕ)⋅Ry(θ)⋅Rz(ψ)

R(θ) = R(θ,θ,θ)

function download_thing(thing_id;path="",verbose::Bool=true)
  url = "https://www.thingiverse.com/download:$thing_id"
  !verbose || println("---------------------------------------")
  r = Downloads.request(url)
  if 200 ≤ r.status ≤ 399
    _,ext = splitext(r.url)
    mkpath(path)
    filename = joinpath(path,"$thing_id$ext")
    !verbose || println("Downloading Thing $thing_id ...")
    Downloads.download(url,filename)
  else
    !verbose || println("Thing $thing_id is no longer available on Thingiverse.")
    filename = nothing
  end
  filename
end

function download_things(thing_ids;path="",verbose=true)
  filenames = []
  for thing_id in thing_ids
    filename = download_thing(thing_id;path,verbose)
    push!(filenames,filename)
  end
  filenames
end

function compute_sizes(pmin::Point{D},pmax::Point{D};nmin=10,nmax=100) where D
  s = pmax - pmin
  s = Tuple(s)
  h = maximum(s)/nmax
  p = ( s./maximum(s) ) .* nmax
  p = ceil.(p)
  if minimum(p) < nmin
    p = ( s./minimum(s) ) .* nmin
    h = minimum(s) / nmin
  end
  origin = pmin
  sizes = tfill(h,Val{D}())
  partition = Int.(ceil.(p))
  origin,sizes,partition
end

function run_stl_cutter(model,stl;tolfactor=1e3,kdtree=false,vtk=false,threading=:spawn,verbose=true,title="")
  t = @timed data = compute_submesh(model,stl;tolfactor,threading,kdtree)
  T,X,F,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf,bgcell_to_ioc = data

  submesh = compute_grid(Table(T),X,TET)
  facets = compute_grid(Table(F),Xf,TRI)

  grid = get_grid(model)

  if vtk
    sf = title*"_subfacets"
    sm = title*"_submesh"
    bg = title*"_bgmesh"
    writevtk(facets,sf,cellfields=["bgcell"=>f_to_bgcell])
    writevtk(submesh,sm,cellfields=["inout"=>k_to_io,"bgcell"=>k_to_bgcell])
    writevtk(grid,bg,cellfields=["inoutcut"=>bgcell_to_ioc])
  end
  
  bgmesh_in_vol = volume(grid,bgcell_to_ioc,FACE_IN)
  bgmesh_out_vol  = volume(grid,bgcell_to_ioc,FACE_OUT)
  bgmesh_cut_vol  = volume(grid,bgcell_to_ioc,FACE_CUT)
  submesh_in_vol  = volume(submesh,k_to_io,FACE_IN)
  submesh_out_vol = volume(submesh,k_to_io,FACE_OUT)

  has_leak = bgmesh_out_vol == 0 || bgmesh_in_vol == 0

  in_volume = bgmesh_in_vol + submesh_in_vol
  out_volume = bgmesh_out_vol + submesh_out_vol
  cut_volume = bgmesh_cut_vol

  num_cut_bgcells = count(isequal(FACE_CUT),bgcell_to_ioc)

  num_k = zeros(Int32,num_cells(grid))
  for bg_cell in k_to_bgcell
    num_k[bg_cell] += 1
  end
  num_k = filter(!iszero,num_k)
  max_num_k = maximum( num_k )
  min_num_k = minimum( num_k )

  domain_surf = surface(facets)
  stl_surf = surface(get_grid(stl)) 

  volume_error = (in_volume + out_volume - volume(grid)) / volume(grid)
  surface_error = ( stl_surf - domain_surf ) / stl_surf

  out = Dict{String,Any}()

  out["time"] = t.time
  out["domain_volume"] = in_volume
  out["domain_surface"] = domain_surf
  out["box_volume"] = volume(grid)
  out["volume_error"] = volume_error
  out["surface_error"] = surface_error
  out["has_leak"] = has_leak
  out["num_cells_in"] = count(isequal(FACE_IN),bgcell_to_ioc)
  out["num_cells_out"] = count(isequal(FACE_OUT),bgcell_to_ioc)
  out["num_cells_cut"] = count(isequal(FACE_CUT),bgcell_to_ioc)
  out["num_subcells"] = length(k_to_io)
  out["num_subcells_in"] = count(isequal(FACE_IN),k_to_io)
  out["num_subcells_out"] = count(isequal(FACE_OUT),k_to_io)
  out["min_subcells_x_cell"] = min_num_k
  out["max_subcells_x_cell"] = max_num_k
  out["avg_subcells_x_cell"] = length(k_to_bgcell) / num_cut_bgcells

  if verbose
    println("---------------------------------------")
    println("TIME: \t $(out["time"]) s")
    println("Domain volume:\t$(out["domain_volume"])")
    println("Domain surface:\t$(out["domain_surface"])")
    println("Box volume: \t $(out["box_volume"])")
    println("Volume error: \t $(out["volume_error"])")
    println("Surface error: \t $(out["surface_error"])")
    !out["has_leak"] || println("Has a leak!")
    println("Num subcells: \t $(out["num_subcells"])")
    println("MIN num subcells per cut cell: \t $(out["min_subcells_x_cell"])")
    println("MAX num subcells per cut cell: \t $(out["max_subcells_x_cell"])")
    println("AVG num subcells per cut cell: \t $(out["avg_subcells_x_cell"])")
    println("---------------------------------------")
  end

  out
end

function _run_and_save(
  filename::String;
  title::String="res",
  verbose::Bool=true,
  vtk::Bool=false,
  δ::Real=0.2,
  nmin::Integer=10,
  nmax::Integer=100,
  Δx::Real=0,
  θ::Real=0,
  kdtree::Bool=false,
  tolfactor::Real=1e3,
  poisson::Bool=false,
  solution_order::Integer=1,
  agfem_threshold::Real=1,
  solver::Symbol=:direct,
  threading::Symbol=:none,
  nruns::Integer=1,
  nthreads::Integer=Threads.nthreads(),
  max_num_facets::Integer=2000)

  X,T,N = read_stl(filename)

  if length(X) == 0 || length(T) == 0
    !verbose || println("Skipping run: no facets or no vertices")
    return 
  end

  stl0 = compute_stl_model(T,X)
  stl0 = merge_and_collapse(stl0)

  X0 = get_node_coordinates(get_grid(stl0))
  T0 = get_cell_node_ids(stl0)

  pmin,pmax = get_bounding_box(stl0)
  Δ = (pmax-pmin)*δ
  origin,sizes,partition = compute_sizes(pmin-Δ,pmax+Δ;nmin,nmax)

  model = CartesianDiscreteModel(origin,sizes,partition)

  if !check_requisites(stl0,model;verbose,max_num_facets)
    !verbose || println("Not matching requisites")
    return
  end

  min_h = min_height(stl0) * (eps()/eps(model))

  Δx_scaled = minimum(pmax-pmin) * Δx

  o = (pmin+pmax)/2
  Xi = map(p-> o + R(θ)⋅(p-o),X0)
  Xi = map(p-> p + Δx_scaled,Xi)
  stl = compute_stl_model(Table(T0),Xi)

  data = Dict{String,Any}()
  
  data["name"] = first(splitext(basename(filename)))
  data["min_h"] = min_h
  data["delta"] = δ
  data["nmin"] = nmin
  data["nmax"] = nmax
  data["partition"] = partition
  data["h"] = first(sizes)
  data["num_cells"] = num_cells(model)
  data["num_stl_facets"] = num_cells(stl)
  data["displacement"] = Δx
  data["scaled_displacement"] = Δx_scaled
  data["rotation"] = θ
  data["kdtree"] = kdtree

  @pack! data = threading, nruns, nthreads
  @pack! data = poisson, solution_order, agfem_threshold, solver

  if verbose
    println("---------------------------------------")
    println("Num stl facets:\t$(data["num_stl_facets"])")
    println("Minimum stl facet height (scaled):\t$(data["min_h"])")
    println("Num background grid cells:\t$(data["num_cells"])")
    println("Background grid partition:\t$(data["partition"])")
    println("Background grid size:\t$(data["h"])")
  end

  @assert nthreads == Threads.nthreads()
  threading = threading == :none ? :spawn : threading
  out = data

  for it in 1:nruns
    _title = title
    if nruns > 1
      _title *= "_it=$it"
    end
    if !poisson
      out = run_stl_cutter(model,stl;tolfactor,kdtree,threading,vtk,verbose,title)
    else
      d = solution_order
      u = x -> x[1]^d + x[2]^d - x[3]^d
      out = run_poisson(model,stl,u;agfem_threshold,solver,vtk,verbose,title)
    end
    @pack! out = it
    merge!(out,data)
    @tagsave("$_title.bson",out;safe=true)
  end
  out
end

run_and_save(::Nothing;kwargs...) = nothing

function run_and_save(
  filename;
  datapath=datadir(),
  verbose::Bool=true,
  vtk::Bool=false,
  rerun::Bool=false,
  params...)

  name = first(splitext(basename(filename)))
  title = savename(name,params,sigdigits=1)
  !verbose || println("---------------------------------------")
  !verbose || println("Running: $title ...")
  title = joinpath(datapath,title)
  if !rerun && isfile("$title.bson")
    !verbose || println("Skipping run: $title.bson is found")
    return
  end
  _run_and_save(filename;title,verbose,vtk,params...)
end

function download_run_and_save(things;kwargs...)
  for thing in things
    filename = download_thing(thing,path=tmpdir())
    run_and_save(filename;kwargs...)
  end
end

function rotations_and_displacements(
  filename;
  displacements=exp10.(-17:-1 ),
  angles=exp10.(-17:-1 ),
  verbose::Bool=true,
  params...)

  !verbose || println("BEGIN of rotations_and_displacements()")
  run_and_save(filename;Δx=0,θ=0,params...)
  for θ in angles
    run_and_save(filename;Δx=0,θ,params...)
  end
  for Δx in displacements
    run_and_save(filename;Δx,θ=0,params...)
  end
  !verbose || println("END of rotations_and_displacements()")
end

function run_geometry_list(ids,filenames;verbose::Bool=true,params...)
  !verbose || println("BEGIN of run_geometry_list()")
  for id in ids
    !isnothing(filenames[id]) || continue
    try
      run_and_save(filenames[id];verbose,params...)
    catch e
      println(stdout,e)
      !verbose || println("Failed run of geometry $id")
    end
  end
  !verbose || println("END of run_geometry_list()")
end

function run_poisson(
  bgmodel::CartesianDiscreteModel,
  stl::DiscreteModel,
  u::Function = x -> x[1] + x[2] - x[3];
  agfem_threshold = 1,
  solver=:direct,
  vtk=false,
  verbose=true,
  title="")

  # Manufactured solution
  f(x) = - Δ(u)(x)
  ud(x) = u(x)
  
  # Cut the background model
  geo = STLGeometry(stl)
  cutter = @timed cutgeo,facet_to_inoutcut = cut(bgmodel,geo)

  strategy = AggregateCutCellsByThreshold(agfem_threshold)
  aggregates = aggregate(strategy,cutgeo,facet_to_inoutcut)

  # Setup integration meshes
  Ω_bg = Triangulation(bgmodel)
  Ω = Triangulation(cutgeo)
  Γd = EmbeddedBoundary(cutgeo)
  Γg = GhostSkeleton(cutgeo)

  # Setup normal vectors
  n_Γd = get_normal_vector(Γd)
  n_Γg = get_normal_vector(Γg)

  # Setup Lebesgue measures
  order = 1
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓd = Measure(Γd,degree)
  dΓg = Measure(Γg,degree)

  integration_volume = sum( ∫(1)*dΩ  )
  integration_surface = sum( ∫(1)*dΓd )
  
  # Setup FESpace
  model = DiscreteModel(cutgeo)
  Vstd = FESpace(model,ReferenceFE(lagrangian,Float64,order),conformity=:H1)

  V = AgFEMSpace(Vstd,aggregates)
  U = TrialFESpace(V)
  # Weak form
  γd = 10.0
  γg = 0.1
  h = get_cartesian_descriptor(bgmodel).sizes[1]

  a(u,v) =
    ∫( ∇(v)⋅∇(u) ) * dΩ +
    ∫( (γd/h)*v*u  - v*(n_Γd⋅∇(u)) - (n_Γd⋅∇(v))*u ) * dΓd +
    ∫( (γg*h)*jump(n_Γg⋅∇(v))*jump(n_Γg⋅∇(u)) ) * dΓg

  l(v) =
    ∫( v*f ) * dΩ +
    ∫( (γd/h)*v*ud - (n_Γd⋅∇(v))*ud ) * dΓd

  # FE problem
  operator = @timed op = AffineFEOperator(a,l,U,V)

  system = @timed if solver == :direct
    uh = solve(op)
  elseif solver == :amg
    A = get_matrix(op)
    b = get_vector(op)
    p = AMGPreconditioner{SmoothedAggregation}(A)
    x = cg(A,b,verbose=true,Pl=p,reltol=1e-10)
    uh = FEFunction(U,x)
  else
    @unreachable
  end

  e = u - uh

  # Postprocess
  l2(u) = sqrt(sum( ∫( u*u )*dΩ ))
  h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ ))

  el2 = l2(e)
  eh1 = h1(e)
  ul2 = l2(uh)
  uh1 = h1(uh)

  error_l2 = el2/ul2
  error_h1 = eh1/uh1

  if vtk
    r = title*"_results"
    t = title*"_trian"
    colors = color_aggregates(aggregates,bgmodel)
    writevtk(Ω_bg,t,celldata=["aggregate"=>aggregates,"color"=>colors],cellfields=["uh"=>uh])
    writevtk(Ω,r,cellfields=["uh"=>uh])
  end

  time_cutter = cutter.time
  time_operator = operator.time
  time_system = system.time

  if verbose
    println("---------------------------------------")
    println("TIME CUTTER:  \t $time_cutter")
    println("TIME OPERATOR:\t $time_operator")
    println("TIME SYSTEM:  \t $time_system")
    println("Integration volume:  \t $integration_volume")
    println("Integration surface: \t $integration_surface")
    println("Error L2: \t $error_l2")
    println("Error H1: \t $error_h1")
  end

  out = Dict{String,Any}()
  @pack! out = time_cutter, time_operator, time_system
  @pack! out = integration_volume, integration_surface
  @pack! out = error_l2, error_h1
  out
end

function check_requisites(filename)
  println("---------------------------------------")
  println("Checking requisites of : $(basename(filename))")
  X,T,N = read_stl(filename)
  if length(X) == 0 || length(T) == 0
    !verbose || println("Skipping run: no facets or no vertices")
    return 
  end
  stl = compute_stl_model(T,X)
  stl = merge_and_collapse(stl)
  
  pmin,pmax = get_bounding_box(stl)
  δ = 0.2
  nmin = 10
  nmax = 100
  Δ = (pmax-pmin)*δ
  origin,sizes,partition = compute_sizes(pmin-Δ,pmax+Δ;nmin,nmax)

  model = CartesianDiscreteModel(origin,sizes,partition)
  min_h = min_height(stl) * (eps()/eps(model))
  println("Scaled min height: $min_h")
  if !STLCutters.check_requisites(stl,model)
    println("@warn Failing requisites")
  end
  println("---------------------------------------")
end

function download_and_check_requisites(ids)
  for id in ids
    filename = download_thing(id,path=tmpdir())
    try
      check_requisites(filename)
    catch e
      println("Failing $id")
      println(e)
    end
  end
end
