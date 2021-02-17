module drwatson

using DrWatson

using Test

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using Gridap.Helpers
using STLCutters

using STLCutters: read_stl
using STLCutters: compute_stl_model
using STLCutters: merge_nodes 
using STLCutters: get_bounding_box 
using STLCutters: compute_submesh 
using STLCutters: compute_grid 
using STLCutters: is_water_tight
using STLCutters: surface, volume, volumes 
using STLCutters:  FACE_IN, FACE_OUT, FACE_CUT

import HTTP

include(joinpath(@__DIR__,"data/thingi10k_quality_filter.jl")) # file_ids=Int[]

## Functions

function download_thing(thing_id;path="")
  filename = nothing
  url = "https://www.thingiverse.com/download:$thing_id"
  try 
    link = HTTP.head(url).request.target
    _,ext = splitext(link)
    filename = "$thing_id$ext"
    println("Downloading Thing $thing_id ...")
    wget_command = `wget -q -O $filename --tries=10 $url`
    run(wget_command)
  catch
    println("Thing $thing_id is no longer available on Thingiverse.")
  end
  filename
end

function download_things(thing_ids,path="")
  for thing_id in thing_ids
    download_thing(thing_id;path)
  end
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

function run_stl_cutter(grid,stl;kdtree=false,vtk=false,verbose=true)
  t = @timed data = compute_submesh(grid,stl)
  T,X,F,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf,bgcell_to_ioc = data

  submesh = compute_grid(Table(T),X,TET)
  facets = compute_grid(Table(F),Xf,TRI)

  if vtk
    writevtk(facets,"subfacets",cellfields=["bgcell"=>f_to_bgcell])
    writevtk(submesh,"submesh",cellfields=["inout"=>k_to_io,"bgcell"=>k_to_bgcell])
    writevtk(grid,"bgmesh",cellfields=["inoutcut"=>bgcell_to_ioc])
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
  out["heas_leak"] = has_leak
  out["min_subcells_x_cell"] = minimum(num_k)
  out["max_subcells_x_cell"] = minimum(num_k)
  out["avg_subcells_x_cell"] = length(k_to_bgcell) / length(cut_bgcells)

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
  T0 = get_cell_nodes(stl0)

  pmin,pmax = get_bounding_box(stl0)
  Δ = (pmax-pmin)*δ
  origin,sizes,partition = compute_sizes(pmin-Δ,pmax+Δ;nmin,nmax)

  grid = CartesianGrid(origin,sizes,partition)

  o = (pmin+pmax)/2
  Xi = map(p-> o + R(θ)⋅(p-o),X0)
  Xi = map(p-> p + Δx,Xi)
  stl = compute_stl_model(Table(T0),Xi)
  
  out = run_stl_cutter(grid,stl;kdtree,vtk,verbose)

  data = Dict{String,Any}()
  
  data["name"] = first(splitext(basename(filename)))
  data["delta"] = δ
  data["nmin"] = nmin
  data["nmax"] = nmax
  data["partition"] = partition
  data["h"] = first(sizes)
  data["num_cells"] = num_cells(grid)
  data["num_stl_facets"] = num_cells(stl)
  data["disp"] = Δx
  data["rot"] = θ
  data["kdtree"] = kdtree

  merge!(out,data)
  safesave("$title.bson",out)
end

run_stl_cutter(::Nothing;kwargs...) = nothing

function download_and_run(things;kwargs...)
  for thing in things
    filename = download_thing(thing)
    title = "$thing" #use DrWatson params 
    run_stl_cutter(filename;title,kwargs...)
  end
end

download_and_run(file_ids[2:3])


end
## Depth test geometries:
#  cube.stl
#  wine_glass.stl
#  standford_bunny.stl -> thing:  441708
#  chichen_itza -> thing: 47076
#  arc_de_triomph -> 551021
#  earth -> 37266
#  heart -> 65904
#  strange thing -> 54725
