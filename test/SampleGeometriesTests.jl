module TestSampleGeometries

using Test
using Gridap
using Gridap.Arrays
using Gridap.Helpers
using STLCutters

import Downloads

using STLCutters: check_requisites
using STLCutters: get_stl
using STLCutters: volumes
using STLCutters: volume, surface

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

  println("Running: $(basename(filename))")

  stl = STLGeometry( filename )

  pmin,pmax = get_bounding_box(stl)
  Δ = (pmax-pmin)*δ
  origin,sizes,partition = compute_sizes(pmin-Δ,pmax+Δ;nmin,nmax)
  model = CartesianDiscreteModel(origin,sizes,partition)

  @test check_requisites(get_stl(stl),model)

  t = @timed subcells,subfaces,labels = subtriangulate(model,stl;tolfactor,kdtree)

  grid = get_grid(model)

  if output !== nothing
    sf = output*"_subfacets"
    sm = output*"_submesh"
    bg = output*"_bgmesh"
    writevtk(subfaces,sf,celldata=["bgcell"=>labels.face_to_bgcell])
    writevtk(subcells,sm,celldata=["inout"=>labels.cell_to_io,"bgcell"=>labels.cell_to_bgcell])
    writevtk(grid,bg,celldata=["inoutcut"=>labels.bgcell_to_ioc])
  end

  bgmesh_in_vol, bgmesh_out_vol, bgmesh_cut_vol = volumes(grid,labels.bgcell_to_ioc)
  submesh_in_vol,submesh_out_vol, = volumes(subcells,labels.cell_to_io)

  has_leak = bgmesh_out_vol == 0 || bgmesh_in_vol == 0

  in_volume = bgmesh_in_vol + submesh_in_vol
  out_volume = bgmesh_out_vol + submesh_out_vol
  cut_volume = bgmesh_cut_vol

  domain_surf = surface(subfaces)

  stl_surf = surface(stl)

  volume_error = (in_volume + out_volume - volume(grid)) / volume(grid)
  surface_error = ( stl_surf - domain_surf ) / stl_surf


  println("TIME: $(t.time)")
  println("Background partition: $partition")
  println("Num background cells: $(num_cells(model))")
  println("Num subcells: $(num_cells(subcells))")
  println("ϵ_V = $volume_error")
  println("ϵ_Γ = $surface_error")

  @test abs(volume_error) < 1e-9
  @test abs(surface_error) < 1e-9
  @test !has_leak

end

filename = joinpath(@__DIR__,"../test/data/47076.stl")
main(filename,nmax=50)
main(filename,nmax=10,nmin=10,kdtree=true)

filename = joinpath(@__DIR__,"../test/data/47076_sf.obj")
main(filename,nmax=50)

filename = download(293137)
main(filename,nmax=50)
rm(filename)

filename = download(80084)
main(filename,nmax=20,nmin=5)
rm(filename)

filename = download(65395)
main(filename,nmax=20,nmin=5)
rm(filename)

filename = download(77343)
main(filename,nmax=20,nmin=5)
rm(filename)

filename = download(95436)
main(filename,nmax=20,nmin=5)
rm(filename)

filename = download(243015)
main(filename,nmax=20,nmin=5)
rm(filename)

filename = download(57657)
main(filename,nmax=100,nmin=5,tolfactor=10000)
rm(filename)

filename = download(1452677)
main(filename,nmax=20,nmin=5)
rm(filename)

filename = download(93703)
main(filename,nmax=20,nmin=5)
rm(filename)

end # module
