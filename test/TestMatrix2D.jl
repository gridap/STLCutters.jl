module TestMatrix2D

using STLCutters

using STLCutters: BulkMesh
using STLCutters: surface
using STLCutters: interior_volume
using STLCutters: exterior_volume
using STLCutters: move

function star(;center=Point(0,0),r_in=0.5,r_out=1.0,n=5)
  points = zeros(Point{2,Float64},n*2)
  Δα = 2π / 2n
  for i in 1:2n
    r = iseven(i) ? r_in : r_out
    α = (i-1) * Δα
    R = VectorValue(r*sin(α),r*cos(α))
    points[i] = center + R
  end
  points
end

stls = STL{2,Float64}[]
base_sizes = Int[]
ref_volumes = Float64[]
geo_names = String[]

v = star(n=8,r_in=0.4)
stl = closed_polyline(v)
vol = 1.224586983568287
push!(stls,stl)
push!(base_sizes,2)
push!(ref_volumes,vol)
push!(geo_names,"star_1")

v = star(center=Point(0.5,0.5),n=8,r_in=0.2,r_out=0.5)
stl = closed_polyline(v)
vol = 0.30614674589207175
push!(stls,stl)
push!(base_sizes,2)
push!(ref_volumes,vol)
push!(geo_names,"star_2")


file = joinpath(@__DIR__,"data/naca.dat")
stl = closed_polyline(file)
vol = 0.09067041510300003
push!(stls,stl)
push!(base_sizes,5)
push!(ref_volumes,vol)
push!(geo_names,"naca")


file = joinpath(@__DIR__,"data/mallorca.csv")
v = STLCutters.read_vertices(file)
#v = v[1:2:end]
stl = closed_polyline(v)
vol = 0.38154564204977465
push!(stls,stl)
push!(base_sizes,10)
push!(ref_volumes,vol)
push!(geo_names,"mallorca")


tols = [ 1e-6, 1e-9, 1e-12 ]
size_factors = [1, 10, 100 ]
max_tol = 1e-5


M = CartesianIndices( (length(stls), length(tols), length(size_factors) ) ) 
#M = CartesianIndices( (4:4, 1:3, 3:3) ) 

for I in M
  stl = stls[ I[1] ]
  tol = tols[ I[2] ]
  size_factor = size_factors[ I[3] ]

  base_size = base_sizes[ I[1] ]
  ref_volume = ref_volumes[ I[1] ]

  n = size_factor * base_size

  geometry_name = geo_names[ I[1] ]
  println("Computing `$geometry_name` with `tolerance = $tol` and `mesh = ($n×$n)` :") 
    
  sm_0 = SurfaceMesh(stl)
  box = 1.2*BoundingBox(stl)
  mesh = CartesianMesh(box,n)
  
  STLCutters.tolerance(::CellMesh) = tol
  STLCutters.tolerance(::STLCutters.IncrementalSurfaceMesh) = tol*10
  offset = tol*one(VectorValue{2,Float64})

  for d in 1:5

    println("  Displacement $d:" )
    
    sm = move(sm_0,(d-1)*offset)

    bulk = BulkMesh(mesh,sm)

    i_vol = interior_volume(bulk)
    e_vol = exterior_volume(bulk)
    b_vol = measure(box)

    b_surf = surface(bulk,1)
    s_surf = surface(sm)

    ε_Ω_in = (i_vol - ref_volume) / ref_volume
    ε_Ω = ( i_vol + e_vol - b_vol ) / b_vol
    ε_Γ = ( b_surf - s_surf ) / s_surf

    println("    Interior Volume  error = $ε_Ω_in ")
    println("    Domain Volume    error = $ε_Ω")
    println("    Boundary Surface error = $ε_Γ")

    if abs(ε_Ω_in) > max_tol || abs(ε_Ω) > max_tol || abs(ε_Γ) > max_tol
      printstyled("WARNING: Error above the permisive tolerance: ", bold=true,color=:red); println();
      printstyled("  ε_Ω_in = $ε_Ω_in", bold=true,color=:red); println();
      printstyled("  ε_Ω    = $ε_Ω", bold=true,color=:red); println();
      printstyled("  ε_Γ    = $ε_Γ", bold=true,color=:red); println();
    end

  end
end
 
end # module
