module BulkMeshesTests

using STLCutter

using STLCutter: expand, BulkMesh

using Test


stl = STL(joinpath(@__DIR__,"data/cube.stl"))

sm = SurfaceMesh(stl)
box = BoundingBox(sm)
box = expand(box,0.2)
bg_mesh = CartesianMesh(box,5)

bm = BulkMesh(bg_mesh,sm)

#stl = STL(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))
#
#sm = SurfaceMesh(stl)
#box = BoundingBox(sm)
#box = expand(box,0.1)
#bg_mesh = CartesianMesh(box,20)
#
#bm = BulkMesh(bg_mesh,sm)
#writevtk(sm,"sm")
#writevtk(bm,"bm")

end # module
