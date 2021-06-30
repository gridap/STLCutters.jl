module SubTriangulation

using Test
using Gridap
using STLCutters

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

function main(filename;n=20,δ=0.2,output=nothing)

  println("Running: $(basename(filename))")

  stl = STLGeometry( filename )

  pmin,pmax = get_bounding_box(stl)
  Δ = (pmax-pmin)*δ
  pmin -= Δ
  pmax += Δ
  partition = (n,n,n)
  model = CartesianDiscreteModel(pmin,pmax,partition)

  subcells,subfaces,labels = subtriangulate(model,stl)

  if output !== nothing
    sf = output*"_subfaces"
    sc = output*"_subcells"
    bg = output*"_bgmesh"
    writevtk(subfaces,sf,celldata=["bgcell"=>labels.face_to_bgcell])
    writevtk(subcells,sc,celldata=["inout"=>labels.cell_to_io,"bgcell"=>labels.cell_to_bgcell])
    writevtk(get_grid(model),bg,celldata=["inoutcut"=>labels.bgcell_to_ioc])
  end

end

end # module
