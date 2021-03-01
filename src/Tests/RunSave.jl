
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
  filename = nothing
  url = "https://www.thingiverse.com/download:$thing_id"
  !verbose || println("---------------------------------------")
  try 
    link = HTTP.head(url).request.target
    _,ext = splitext(link)
    mkpath(path)
    filename = joinpath(path,"$thing_id$ext")
    !verbose || println("Downloading Thing $thing_id ...")
    wget_command = `wget -q -O $filename --tries=10 $url`
    run(wget_command)
  catch
    !verbose || println("Thing $thing_id is no longer available on Thingiverse.")
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

function run_stl_cutter(grid,stl;kdtree=false,vtk=false,verbose=true,title="")
  t = @timed data = compute_submesh(grid,stl)
  T,X,F,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf,bgcell_to_ioc = data

  submesh = compute_grid(Table(T),X,TET)
  facets = compute_grid(Table(F),Xf,TRI)

  if vtk
    sf = title*"_subfacets"
    sm = title*"_submesh"
    bg = title*"_bgmesh"
    writevtk(facets,sf,cellfields=["bgcell"=>f_to_bgcell])
    writevtk(submesh,sm,cellfields=["inout"=>k_to_io,"bgcell"=>k_to_bgcell])
    writevtk(grid,bg,cellfields=["inoutcut"=>bgcell_to_ioc])
  end

  bgmesh_vols = volumes(grid)
  submesh_vols = volumes(submesh)

  bgmesh_in_vols = bgmesh_vols[findall(isequal(FACE_IN),bgcell_to_ioc)]
  bgmesh_out_vols = bgmesh_vols[findall(isequal(FACE_OUT),bgcell_to_ioc)]
  bgmesh_cut_vols = bgmesh_vols[findall(isequal(FACE_CUT),bgcell_to_ioc)]
  submesh_in_vols = submesh_vols[findall(isequal(FACE_IN),k_to_io)]
  submesh_out_vols = submesh_vols[findall(isequal(FACE_OUT),k_to_io)]

  has_leak = sum(bgmesh_out_vols) == 0 || sum(bgmesh_in_vols) == 0

  in_volume = sum(bgmesh_in_vols) + sum(submesh_in_vols)
  out_volume = sum(bgmesh_out_vols) + sum(submesh_out_vols)
  cut_volume = sum(bgmesh_cut_vols)

  cut_bgcells = unique(k_to_bgcell)
  num_k = [ count(isequal(bgcell),k_to_bgcell) for bgcell in cut_bgcells]
  num_f = [ count(isequal(bgcell),f_to_bgcell) for bgcell in cut_bgcells]

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
  out["min_subcells_x_cell"] = minimum(num_k)
  out["max_subcells_x_cell"] = maximum(num_k)
  out["avg_subcells_x_cell"] = length(k_to_bgcell) / length(cut_bgcells)

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

function run_stl_cutter(
  filename::String;
  title::String="res",
  verbose::Bool=true,
  vtk::Bool=false,
  δ::Real=0.2,
  nmin::Integer=10,
  nmax::Integer=100,
  Δx::Real=0,
  θ::Real=0,
  kdtree::Bool=false)

  X,T,N = read_stl(filename)
  stl0 = compute_stl_model(T,X)
  stl0 = merge_nodes(stl0)

  X0 = get_node_coordinates(get_grid(stl0))
  T0 = get_cell_node_ids(stl0)

  pmin,pmax = get_bounding_box(stl0)
  Δ = (pmax-pmin)*δ
  origin,sizes,partition = compute_sizes(pmin-Δ,pmax+Δ;nmin,nmax)

  grid = CartesianGrid(origin,sizes,partition)

  min_h = min_height(stl0) * (eps()/eps(grid))
  if min_h < eps()*1e3
    !verbose || println("Skipping run: min facet height = $(min_h)")
    return 
  end

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
  data["num_cells"] = num_cells(grid)
  data["num_stl_facets"] = num_cells(stl)
  data["displacement"] = Δx
  data["scaled_displacement"] = Δx_scaled
  data["rotation"] = θ
  data["kdtree"] = kdtree

  if verbose
    println("---------------------------------------")
    println("Num stl facets:\t$(data["num_stl_facets"])")
    println("Minimum stl facet height (scaled):\t$(data["min_h"])")
    println("Num background grid cells:\t$(data["num_cells"])")
    println("Background grid partition:\t$(data["partition"])")
    println("Background grid size:\t$(data["h"])")
  end

  out = run_stl_cutter(grid,stl;kdtree,vtk,verbose,title)

  merge!(out,data)
  @tagsave("$title.bson",out;safe=true)
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
  title = savename(name,params,scientific=1)
  !verbose || println("---------------------------------------")
  !verbose || println("Running: $title ...")
  title = joinpath(datapath,title)
  if !rerun && isfile("$title.bson")
    !verbose || println("Skipping run: $title.bson is found")
    return
  end
  run_stl_cutter(filename;title,verbose,vtk,params...)
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
    run_and_save(filenames[id];verbose,params...)
  end
  !verbose || println("END of run_geometry_list()")
end

