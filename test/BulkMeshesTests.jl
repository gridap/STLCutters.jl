module BulkMeshesTests

using STLCutter

using STLCutter: expand, BulkMesh, NewBulkMesh

using Test


stl = STL(joinpath(@__DIR__,"data/cube.stl"))
stl = STL(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))

sm = SurfaceMesh(stl)

box = BoundingBox(sm)

#box = expand(box,0.5)
box = expand(box,0.1)

bg_mesh = CartesianMesh(box,20)


nbm = NewBulkMesh(bg_mesh,sm)
bm = BulkMesh(bg_mesh,sm)
writevtk(bm,"bulk")
writevtk(sm,"sm")

writevtk(nbm,"new_bm")
end # module
