module Rotations
  
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

function test_stl_cut(grid,stl,vol)
  data = compute_submesh(grid,stl)
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

  @test volume(submesh) ≈ cut_volume 
  @test surface(get_grid(stl)) ≈ surface(facets)
  @test in_volume + out_volume ≈ volume(grid)
  @test in_volume ≈ vol
  
  println("\t εV = $(in_volume + out_volume - volume(grid))")
  println("\t εVin = $(in_volume-vol)")
  println("\t εΓ = $(surface(get_grid(stl)) - surface(facets))")

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


#X,T,N = read_stl(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))
#X,T,N = read_stl(joinpath(@__DIR__,"data/wine_glass.stl"))
X,T,N = read_stl(joinpath(@__DIR__,"data/cube.stl"))

stl0 = compute_stl_model(T,X)
stl0 = merge_nodes(stl0)

X0 = get_node_coordinates(get_grid(stl0))
T0 = get_cell_node_ids(stl0)

## Coindident Grid-Cube
δ = 0.2
n = 7
D = 3

pmin,pmax = get_bounding_box(stl0)
O = (pmin+pmax)/2

Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)

grid = CartesianGrid(pmin,pmax,partition)


θs = [0;exp10.( -17:0 )]
Δx = VectorValue(0.0,0.0,0.1)

for θ in θs
  println("Testing θz = $θ ...")
  Oi = O + Δx
  Xi = map(p-> p + Δx,X0)
  Xi = map(p-> Oi + Rz(θ)⋅(p-Oi),Xi)
  stl = compute_stl_model(Table(T0),Xi)
  #writevtk(stl.grid,"stl")
  test_stl_cut(grid,stl,1)
end

## Rational Grid
δ = 0.2
n = 7
D = 3

pmin,pmax = get_bounding_box(stl0)
pmin = Point(rationalize.(Tuple(pmin)))
pmax = Point(rationalize.(Tuple(pmax)))
δ = rationalize(δ)

O = (pmin+pmax)/2

Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)

grid = CartesianGrid(pmin,pmax,partition)

δx = Point(0,0,1)
Δxs = [0;exp10.(-17:-1)]

for Δx in Δxs
  println("Testing Δx = $Δx ...")
  Xi = map(p-> p + Δx,X0)
  stl = compute_stl_model(Table(T0),Xi)
  writevtk(stl.grid,"stl")
  test_stl_cut(grid,stl,1)
end

δ = 0.2
n = 7
D = 3

pmin,pmax = get_bounding_box(stl0)
pmin = Point(rationalize.(Tuple(pmin)))
pmax = Point(rationalize.(Tuple(pmax)))
δ = rationalize(δ)

origins = [pmin, pmax, (pmin+pmax)/2 ]

Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)

grid = CartesianGrid(pmin,pmax,partition)

θs = [0;exp10.( -17:-1 )]

for O in origins
  for θ in θs
    println("Testing θ = $θ over $O ...")
    Xi = map(p-> O + R(θ)⋅(p-O),X0)
    stl = compute_stl_model(Table(T0),Xi)
    writevtk(stl,"stl")
    test_stl_cut(grid,stl,1)
  end
end

end # module
