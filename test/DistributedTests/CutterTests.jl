module DistributedCutterTests

using Gridap
using STLCutters
using GridapDistributed
using PartitionedArrays
using GridapEmbedded
using Test

function main(distribute;
  np = (1,1,1),
  nc = (2,2,2),
  geoname = "cube",
  δ = 0.2,
  tolfactor=10^4,
  simplex=false,
  vtk = false,
  verbose = false)


  filename = joinpath(@__DIR__,"..","data","$geoname.stl")
  if !isfile(filename)
    filename = download_thingi10k(id;path="")
  end

  ranks = distribute(LinearIndices((prod(np),)))
  geo = STLGeometry(filename)

  pmin,pmax = get_bounding_box(geo)
  diagonal = pmax-pmin
  pmin = pmin - diagonal*δ
  pmax = pmax + diagonal*δ
  bgmodel = CartesianDiscreteModel(ranks,np,pmin,pmax,nc)
  if simplex
    bgmodel = simplexify(bgmodel,positive=true)
  end

  cutter = STLCutter(;tolfactor)
  cutgeo = cut(cutter,bgmodel,geo)

  Ωbg = Triangulation(bgmodel)
  Ωin = Triangulation(cutgeo,PHYSICAL_IN)
  Ωout = Triangulation(cutgeo,PHYSICAL_OUT)
  Γ = EmbeddedBoundary(cutgeo)
  degree = 2
  dΩbg = Measure(Ωbg,degree)
  dΩin = Measure(Ωin,degree)
  dΩout = Measure(Ωout,degree)
  dΓ = Measure(Γ,degree)

  surf = ∑( ∫(1)dΓ )
  vol_in = ∑( ∫(1)dΩin )
  vol_out = ∑( ∫(1)dΩout )
  vol_box = ∑( ∫(1)dΩbg )

  if verbose
    println("Surface : ", surf)
    println("Volume in: ", vol_in)
    println("Volume out: ", vol_out)
    println("Volume box: ", vol_box)
    println("Volume error: ", vol_in + vol_out - vol_box)
  end

  if vtk
    writevtk(Ωin,"Ω")
    writevtk(Ωout,"Ωout")
    writevtk(Ωbg,"Ωbg")
    writevtk(Γ,"Γ")
  end

  @test vol_in + vol_out - vol_box < 1e-10
  ref_vol,ref_surf = reference_volume_and_surface(geoname)
  if ref_vol > 0
    @test abs(vol_in - ref_vol) < 1e-10
    @test abs(surf - ref_surf) < 1e-10
  end
end

function reference_volume_and_surface(geoname)
  volumes = Dict("cube"=>1.0,)
  surfaces = Dict("cube"=>6.0,)
  if haskey(volumes,geoname)
    volumes[geoname],surfaces[geoname]
  else
    0.0,0.0
  end
end

end # module
