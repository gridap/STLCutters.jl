module TestMatrix

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
using STLCutters: volume, volumes 
using STLCutters:  FACE_IN, FACE_OUT, FACE_CUT

const surf =  STLCutters.surface

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
  @test surf(get_grid(stl)) ≈ surf(facets)
  @test in_volume + out_volume ≈ volume(grid)
  @test in_volume ≈ vol
  
  volume_error = (in_volume + out_volume - volume(grid)) / volume(grid)
  in_volume_error = (in_volume - vol) / vol
  surface_error = ( surf(get_grid(stl)) - surf(facets) ) / surf(get_grid(stl))

  println("\t εV = $volume_error ")
  println("\t εVin = $in_volume_error")
  println("\t εΓ = $surface_error ")
  

  volume_error, in_volume_error, surface_error
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
  cos(ψ)  ,sin(ψ),0,
  -sin(ψ),cos(ψ),0,
  0,0,1)

R(ϕ,θ,ψ) = Rx(ϕ)⋅Ry(θ)⋅Rz(ψ)

R(θ) = R(θ,θ,θ)



geometries = [ "cube", "Bunny-LowPoly"] #, "wine_glass" ]
size_factors = [1, 2] #, 4, 8 ]

base_sizes = [ 5, 10, 20 ]
ref_volumes = [ 1, 273280.0337419614, 74.12595970063474 ]

M = CartesianIndices( (length(geometries), length(size_factors) ) ) 

θs = exp10.( -17:-10 )
Δxs = [] #[exp10.(-17:-5)]

δ = 0.2

displacement_errors = []
rotation_errors = []
for I in M
  geometry = geometries[ I[1] ]
  size_factor = size_factors[ I[2] ]

  base_size = base_sizes[ I[1] ]
  ref_volume = ref_volumes[ I[1] ]

  n = size_factor * base_size

  X,T,N = read_stl(joinpath(@__DIR__,"data/$geometry.stl"))

  stl0 = compute_stl_model(T,X)
  stl0 = merge_nodes(stl0)

  X0 = get_node_coordinates(get_grid(stl0))
  T0 = get_cell_nodes(stl0)

  pmin,pmax = get_bounding_box(stl0)
  O = (pmin+pmax)/2

  Δ = (pmax-pmin)*δ
  pmin = pmin - Δ
  pmax = pmax + Δ
  partition = (n,n,n)

  grid = CartesianGrid(pmin,pmax,partition)

  println("Testing `$geometry.stl` with `mesh = ($n×$n×$n)` :")

  r_errors = []
  for θ in θs
    println("Testing θ = $θ ...")
    Xi = map(p-> O + R(θ)⋅(p-O),X0)
    stl = compute_stl_model(Table(T0),Xi)
    #writevtk(stl.grid,"stl")
    errors = test_stl_cut(grid,stl,ref_volumes[ I[1] ])
    push!(r_errors,errors)
  end

  d_errors = []
  for Δx in Δxs
    println("Testing Δx = $Δx ...")
    Xi = map(p-> p + Δx,X0)
    stl = compute_stl_model(Table(T0),Xi)
    #writevtk(stl.grid,"stl")
    errors = test_stl_cut(grid,stl,ref_volumes[ I[1] ])
    push!(d_errors,errors)
  end

  push!(rotation_errors,r_errors)
  push!(displacement_errors,d_errors)

end

using Plots

err_ids = [2,3]
err_labels = [ "Domain volume variation (εᵥ)", "Domain surface variation (εₛ)" ]
err_tags = [ "vol", "surf" ]

for (err_id,err_label,err_tag) in zip(err_ids,err_labels,err_tags)

  for (igeom,geom) in enumerate(geometries)
    plot(markershape=:auto)
    for (isf,sf) in enumerate(size_factors)
      n = base_sizes[igeom]*sf
      i = LinearIndices(M)[igeom,isf]
      errors = [ e[error_id] for e in rotation_errors[i] ]
      errors = abs.(errors .- errors[1])
      scatter!(θs,errors,label="n = $n")
    end
    plot!(xscale=:log10)
    plot!(xlabel="Rotation angle (θ) [rads]")
    plot!(ylabel=err_label)
    savefig("$(geom)_$(err_tag)_error_rotation.pdf")
  end
end


end # module

# TODO: save data and separate data generation from plotting, 
# so we can recover it easily.
