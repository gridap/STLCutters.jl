module TestMatrix

using STLCutters

using STLCutters: BulkMesh
using STLCutters: surface
using STLCutters: interior_volume
using STLCutters: exterior_volume
using STLCutters: move!

geometries = [ "cube", "Bunny-LowPoly", "wine_glass" ]
tols = [ 1e-6, 1e-9, 1e-12 ]
size_factors = [1, 2, 4, 8 ]

base_sizes = [ 5, 10, 20 ]
ref_volumes = [ 1, 273280.0337419614, 74.12595970063474 ]

max_tol = 1e-5
M = CartesianIndices( (length(geometries), length(tols), length(size_factors) ) ) 
#M = CartesianIndices( (3:3, 1:3, 2:4) ) 

for I in M
  geometry = geometries[ I[1] ]
  tol = tols[ I[2] ]
  size_factor = size_factors[ I[3] ]

  base_size = base_sizes[ I[1] ]
  ref_volume = ref_volumes[ I[1] ]

  n = size_factor * base_size

  println("Computing `$geometry.stl` with `tolerance = $tol` and `mesh = ($n×$n×$n)` :") 
    
  stl = STL(joinpath(@__DIR__,"data/$geometry.stl"))
  sm = SurfaceMesh(stl)
  box = 1.2*BoundingBox(stl)
  mesh = CartesianMesh(box,n)
  
  STLCutters.tolerance(::CellMesh) = tol
  offset = tol*one(VectorValue{3,Float64})

  for d in 1:10

    println("  Displacement $d:" )
    
    bulk = BulkMesh(mesh,sm)

    i_vol = interior_volume(bulk)
    e_vol = exterior_volume(bulk)
    b_vol = measure(box)

    b_surf = surface(bulk,1)
    s_surf = surface(sm)

    ε_Ω_in = i_vol - ref_volume
    ε_Ω = i_vol + e_vol - b_vol
    ε_Γ = b_surf - s_surf

    println("    Interior Volume  error = $ε_Ω_in ")
    println("    Domain Volume    error = $ε_Ω")
    println("    Boundary Surface error = $ε_Γ")

    if abs(ε_Ω_in) > max_tol || abs(ε_Ω) > max_tol || abs(ε_Γ) > max_tol
      printstyled("WARNING: Error above the permisive tolerance: ", bold=true,color=:red); println();
      printstyled("  ε_Ω_in = $ε_Ω_in", bold=true,color=:red); println();
      printstyled("  ε_Ω    = $ε_Ω", bold=true,color=:red); println();
      printstyled("  ε_Γ    = $ε_Γ", bold=true,color=:red); println();
    end

    sm = move!(sm,offset)
  end
end
 
end # module
