module SubTriangulation

using Test
using Gridap
using Gridap.Arrays
using Gridap.Helpers
using STLCutters

using STLCutters: FACE_IN, FACE_OUT, FACE_CUT

import Downloads

function download(id;path="")
  url = "https://www.thingiverse.com/download:$id"
  r = Downloads.request(url)
  if 200 ≤ r.status ≤ 399
    _,ext = splitext(r.url)
    mkpath(path)
    filename = joinpath(path,"$id$ext")
    Downloads.download(url,filename)
  else
    @warn "$id is no longer available on Thingiverse."
    filename = nothing
  end
  filename
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

function main(filename;
  nmin=10,nmax=100,δ=0.2,tolfactor=1000,kdtree=false,output=nothing)

  println(join(fill('-',40)))
  println("Running: $(basename(filename))")

  X,T,N = read_stl(filename)
  stl = compute_stl_model(T,X)
  stl = merge_and_collapse(stl)
  pmin,pmax = get_bounding_box(stl)
  Δ = (pmax-pmin)*δ
  origin,sizes,partition = compute_sizes(pmin-Δ,pmax+Δ;nmin,nmax)
  model = CartesianDiscreteModel(origin,sizes,partition)

  t = @timed data = compute_submesh(model,stl;tolfactor,kdtree)
  T,X,F,Xf,k_to_io,k_to_bgcell,f_to_bgcell,f_to_stlf,bgcell_to_ioc = data
  
  submesh = compute_grid(Table(T),X,TET)
  facets = compute_grid(Table(F),Xf,TRI)

  grid = get_grid(model)

  if output !== nothing
    sf = output*"_subfacets"
    sm = output*"_submesh"
    bg = output*"_bgmesh"
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

  domain_surf = surface(facets)
  stl_surf = surface(get_grid(stl)) 

  volume_error = (in_volume + out_volume - volume(grid)) / volume(grid)
  surface_error = ( stl_surf - domain_surf ) / stl_surf


  println("TIME: $(t.time)")
  println("Background partition: $partition")
  println("Num background cells: $(num_cells(model))")
  println("Num subcells: $(num_cells(submesh))")
  println("ϵ_V = $volume_error")
  println("ϵ_Γ = $surface_error")

  @test abs(volume_error) < 1e-9
  @test abs(surface_error) < 1e-9
  @test !has_leak

end

end # module
