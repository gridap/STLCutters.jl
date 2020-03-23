module BulkMeshesTests

using STLCutter

using STLCutter: expand, BulkMesh, surface, interior_volume, exterior_volume

using Test


stl = STL(joinpath(@__DIR__,"data/cube.stl"))

sm = SurfaceMesh(stl)
box = BoundingBox(sm)
box = expand(box,0.2)
bg_mesh = CartesianMesh(box,5)

bm = BulkMesh(bg_mesh,sm)

writevtk(sm,"sm")
writevtk(bm,"bm")

@test interior_volume(bm) ≈ 1

@test interior_volume(bm) + exterior_volume(bm) ≈ measure(box)

@test surface(sm) ≈ surface(bm,1)

stl = STL(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))
stl = STL(joinpath(@__DIR__,"data/wine_glass.stl"))
#
sm = SurfaceMesh(stl)
box = BoundingBox(sm)
box = expand(box,0.1)

bg_mesh = CartesianMesh(box,20)

bm = BulkMesh(bg_mesh,sm)

writevtk(sm,"sm")
writevtk(bm,"bm")

@test interior_volume(bm) + exterior_volume(bm) ≈ measure(box)

@test surface(sm) == surface(bm,1)

end # module
