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


#M = CartesianIndices( (length(geometries), length(tols), length(size_factors) ) ) 
M = CartesianIndices( (1, 3, 3) ) # simplified test

for I in M
  geometry = geometries[ I[1] ]
  tol = tols[ I[2] ]
  size_factor = size_factors[ I[3] ]

  base_size = base_sizes[ I[1] ]
  ref_volume = ref_volumes[ I[1] ]

  n = size_factor * base_size

  println("Computing `$geometry` with `tolerance = $tol` and `mesh = ($n×$n×$n)` :") 
    
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

    println("    Interior Volume  error = $(i_vol-ref_volume)")
    println("    Domain Volume    error = $(i_vol+e_vol-b_vol)")
    println("    Boundary Surface error = $(b_surf-s_surf)")

    sm = move!(sm,offset)
  end
end
 
end # module
