module thingi10k

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
using STLCutters: surface, volume, volumes 
using STLCutters:  FACE_IN, FACE_OUT, FACE_CUT


using Pkg; Pkg.add("HTTP")
import HTTP

include(joinpath(@__DIR__,"data/thingi10k_quality_filter.jl")) # file_ids=Int[]


for file_id in file_ids

#file_id = file_ids[i]
url = "https://www.thingiverse.com/download:$file_id"
stl_filename = nothing
try 
  link = HTTP.head(url).request.target
  _,ext = splitext(link)
  stl_filename = "$file_id$ext"
  println("Downloading Thing $file_id ...")
  wget_command = `wget -q -O $stl_filename --tries=10 $url`
  run(wget_command)
catch
  println("Thing $file_id is no longer available on Thingiverse.")
  continue
end

X,T,N = read_stl(stl_filename)
stl = compute_stl_model(T,X)
stl = merge_nodes(stl)

println("\tnum facets:\t$(length(T))")
println("\tnum vertices:\t$(length(X))")

rm(stl_filename)

end
@unreachable

# NOTE: For HPC, first download a fraction of the dataset and then run.
# Only login node will have Internet access





#X,T,N = read_stl(joinpath(@__DIR__,"data/37322_sf.obj"))
#X,T,N = read_stl(joinpath(@__DIR__,"data/441708_sf.obj"))
#X,T,N = read_stl(joinpath(@__DIR__,"data/Bunny.stl"))
#X,T,N = read_stl(joinpath(@__DIR__,"data/one_Bird.stl"))
#X,T,N = read_stl(joinpath(@__DIR__,"data/chichen_itza.stl"))
#X,T,N = read_stl(joinpath(@__DIR__,"data/47076_sf.obj"))
X,T,N = read_stl(joinpath(@__DIR__,"data/32770.stl"))
#X,T,N = read_stl(joinpath(@__DIR__,"data/Octocat-v1.stl"))

stl = compute_stl_model(T,X)
stl = merge_nodes(stl)
writevtk(stl,"stl_obj")

pmin,pmax = get_bounding_box(stl)

δ = 0.2
n = 100
D = 3

Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)


grid = CartesianGrid(pmin,pmax,partition)



t = @timed data = compute_submesh(grid,stl)
T,X,F,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf,bgcell_to_ioc = data

submesh = compute_grid(Table(T),X,TET)
facets = compute_grid(Table(F),Xf,TRI)

writevtk(facets,"subfacets",cellfields=["bgcell"=>f_to_bgcell])
writevtk(submesh,"submesh",cellfields=["inout"=>k_to_io,"bgcell"=>k_to_bgcell])
writevtk(grid,"bgmesh",cellfields=["inoutcut"=>bgcell_to_ioc])

bgmesh_vols = volumes(grid)
submesh_vols = volumes(submesh)

bgmesh_in_vols = bgmesh_vols[findall(isequal(FACE_IN),bgcell_to_ioc)]
bgmesh_out_vols = bgmesh_vols[findall(isequal(FACE_OUT),bgcell_to_ioc)]
bgmesh_cut_vols = bgmesh_vols[findall(isequal(FACE_CUT),bgcell_to_ioc)]
submesh_in_vols = submesh_vols[findall(isequal(FACE_IN),k_to_io)]
submesh_out_vols = submesh_vols[findall(isequal(FACE_OUT),k_to_io)]

in_volume = sum(bgmesh_in_vols) + sum(submesh_in_vols)
out_volume = sum(bgmesh_out_vols) + sum(submesh_out_vols)
cut_volume = sum(bgmesh_cut_vols)
 
println("Time: $(t.time)")
println("Num subcells: $(num_cells(submesh))")
println("Domain volume:  $(in_volume)") 
println("Domain surface:  $(surface(facets))") 
@show surface(get_grid(stl)) - surface(facets)
@show in_volume + out_volume - volume(grid)
@show bgmesh_in_vols > 0 && bgmesh_in_vols > 0

## Data to be stored
# - Id
# - N
# - N facets
# - Domain Volume
# - Surface error
# - Diff vol error
# - Num subcells
# - Max subcells/cell
# - Avg subcells/cell
# - Min subcells/cell
#
# Save at each step
#
#  for ...
#    ...
#    push!(table,row)
#    CSV.write("out.csv",table)
#   end
# Save .vtu's
#
# Input: 
#  - Folder of .stl + .obj
#  - List of id's
#  - Output format (↑)
#  - Loop script
#    - Test few geometries
#    - Launch on HPC


## Script:
#
# Get geos list from (download py script and crete list in text file)
#
# https://ten-thousand-models.appspot.com/results.html?q=is+closed%2C+is+oriented%2C+is+manifold%2C+is+not+degenerate%2C+without+self-intersection%2C+%23df%3D0
#
# Based on py script:
#  url = "https://www.thingiverse.com/download:$file_id"
#  Check headers for extension and availability
#  
#  wget -q -O $out_file $link







end # module
