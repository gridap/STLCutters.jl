module BulkMeshesTests

using STLCutter

using STLCutter: expand, BulkMesh, surface, interior_volume, exterior_volume, move!

using Test

geometries = [ "cube", "Bunny-LowPoly", "wine_glass" ]
volumes = [ 1, 273280.0337419614, 74.12595970063474 ]
meshes = [ 5, 40, 20 ]


offset = VectorValue(1e-7,1e-7,1e-7)

for (i,geom) in enumerate( geometries )

  stl = STL(joinpath(@__DIR__,"data/$geom.stl"))

  sm = SurfaceMesh(stl)
  
 for j in 1:10

    box = BoundingBox(sm)
    box = expand(box,0.1)

    bg_mesh = CartesianMesh(box, meshes[i] )

    bm = BulkMesh(bg_mesh,sm)

    writevtk(sm,"$(geom)_sm")
    writevtk(bm,"$(geom)_bm")

    @test interior_volume(bm) ≈ volumes[i]

    @test interior_volume(bm) + exterior_volume(bm) ≈ measure(box)

    @test surface(sm) ≈ surface(bm,1)

   sm = move!(sm,offset)

 end

end

end # module
