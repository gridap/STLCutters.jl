module Rotationsests
  
using Test

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using STLCutters

using STLCutters: read_stl
using STLCutters: compute_stl_model
using STLCutters: merge_nodes 
using STLCutters: get_bounding_box 
using STLCutters: compute_submesh 
using STLCutters: compute_grid 
using STLCutters: volume, volumes 
using STLCutters:  FACE_IN, FACE_OUT, FACE_CUT



R(θ) = TensorValue(cos(θ),sin(θ),0.0,-sin(θ),cos(θ),0.0,0.0,0.0,1.0)
origin = Point(0.0,0.0,0.5)

X,T,N = read_stl(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))
X,T,N = read_stl(joinpath(@__DIR__,"data/wine_glass.stl"))
X,T,N = read_stl(joinpath(@__DIR__,"data/cube.stl"))

stl = compute_stl_model(T,X)
stl = merge_nodes(stl)



p = HEX
δ = 0.2
n = 7
D = 3

pmin,pmax = get_bounding_box(stl)
Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)


#Rational
grid = CartesianGrid(pmin,pmax,partition)

X0 = get_node_coordinates(get_grid(stl))
T0 = get_cell_nodes(stl)

θs = exp10.( -17:0 )

@show θs

for θ in θs
Xi = map(p-> p + Point(0.0,0.0,0.1),X0)
Xi = map(p-> origin + R(θ)⋅(p-origin),Xi)

stl = compute_stl_model(Table(T0),Xi)
writevtk(get_grid(stl),"stl")

data = compute_submesh(grid,stl)
T,X,F,Xf,k_to_io,k_to_bgcell,f_to_bgcell,bgcell_to_ioc = data

submesh = compute_grid(Table(T),X,TET)
facets = compute_grid(Table(F),Xf,TRI)

writevtk(facets,"subfacets")
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

@test volume(submesh) ≈ cut_volume 
@test surface(get_grid(stl)) ≈ surface(facets)
@test in_volume + out_volume ≈ volume(grid)
@test in_volume ≈ 1
println("Error for θ = $θ : $(in_volume-1)")

end


# TODO:
#  * test exact matching with rational CartesianGrid
#  * Test more displacements
#  * rotate in more axis ϕ,θ,ψ (x,y,z) ( roll, pitch, yaw)
#  * more geometries
#  * more displacements
end # module
