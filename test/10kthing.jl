module thing10k

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



#X,T,N = read_stl(joinpath(@__DIR__,"data/37322_sf.obj"))
X,T,N = read_stl(joinpath(@__DIR__,"data/441708_sf.obj"))
#X,T,N = read_stl(joinpath(@__DIR__,"data/Bunny.stl"))

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



@time data = compute_submesh(grid,stl)
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
 
println("Domain volume:  $(in_volume)") 
@test surface(get_grid(stl)) ≈ surface(facets)

# Solve poisson
# get timmings

end # module
