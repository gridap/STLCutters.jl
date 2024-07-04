module TestSampleGeometries

using Test
using Gridap
using Gridap.Arrays
using Gridap.Helpers
using STLCutters
using Gridap.Geometry

using STLCutters: volumes
using STLCutters: volume, surface
using STLCutters: compute_cartesian_descriptor

download = download_thingi10k

function main(filename;
  nmin=10,nmax=100,δ=0.2,tolfactor=1000,
  kdtree=false,simplex=false,unstructured=false,output=nothing)

  println("Running: $(basename(filename))")

  stl = STLGeometry( filename )

  pmin,pmax = get_bounding_box(stl)
  Δ = (pmax-pmin)*δ
  desc = compute_cartesian_descriptor(pmin-Δ,pmax+Δ;nmin,nmax)
  model = CartesianDiscreteModel(desc)
  if simplex
    model = simplexify(model,positive=true)
  elseif unstructured
    model = UnstructuredDiscreteModel(model)
  end

  @test check_requisites(stl,model)

  t = @timed subcells,subfaces,_,labels = subtriangulate(model,stl;tolfactor,kdtree)

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
  println("Background partition: $(desc.partition)")
  println("Num background cells: $(num_cells(model))")
  println("Num subcells: $(num_cells(subcells))")
  println("ϵ_V = $volume_error")
  println("ϵ_Γ = $surface_error")

  @test abs(volume_error) < 1e-9
  @test abs(surface_error) < 1e-9
  @test !has_leak

end

function download_or_local(id;download=true)
  filename_dwl = download_thingi10k(id)
  filename_stl = joinpath(@__DIR__,"../test/data/$(id).stl")
  filename_obj = joinpath(@__DIR__,"../test/data/$(id).obj")
  if !isnothing(filename_dwl) && download
    println("Downloaded: $(basename(filename_dwl))")
    filename = filename_dwl
  elseif isfile(filename_stl)
    println("Local: $(basename(filename_stl))")
    filename = filename_stl
  elseif isfile(filename_obj)
    println("Local: $(basename(filename_obj))")
    filename = filename_obj
  else
    error("File id:$(id) not found")
  end
  filename
end

function rm_dwl(filename)
  if match(r"test/data$",dirname(filename)) === nothing
    rm(filename)
  end
end

filename = joinpath(@__DIR__,"../test/data/47076.stl")
main(filename,nmax=50)
main(filename,nmax=10,nmin=10,kdtree=true)
main(filename,nmax=50,nmin=10,simplex=true)
main(filename,nmax=50,nmin=10,unstructured=true)

filename = joinpath(@__DIR__,"../test/data/47076_sf.obj")
main(filename,nmax=50)

filename = download_or_local(293137)
main(filename,nmax=50)
main(filename,nmax=50,simplex=true)
main(filename,nmax=50,unstructured=true)
rm_dwl(filename)

filename = download_or_local(80084)
main(filename,nmax=20,nmin=5)
rm_dwl(filename)

filename = download_or_local(65395)
main(filename,nmax=20,nmin=5)
rm_dwl(filename)

filename = download_or_local(77343)
main(filename,nmax=20,nmin=5)
rm_dwl(filename)

filename = download_or_local(95436)
main(filename,nmax=20,nmin=5)
rm_dwl(filename)

filename = download_or_local(243015)
main(filename,nmax=20,nmin=5)
rm_dwl(filename)

filename = download_or_local(57657)
main(filename,nmax=100,nmin=5,tolfactor=10^5)
rm_dwl(filename)

filename = download_or_local(1452677)
main(filename,nmax=20,nmin=5)
rm_dwl(filename)

filename = download_or_local(93703)
main(filename,nmax=20,nmin=5)
rm_dwl(filename)

filename = download_or_local(94492)
main(filename,nmax=100)
rm_dwl(filename)

end # module
