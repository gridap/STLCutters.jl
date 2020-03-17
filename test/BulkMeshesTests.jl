module BulkMeshesTests

using STLCutter

using STLCutter: expand, FACE_UNDEF, reset!, add_surface_mesh_face!,compute_in_out!,is_surface_mesh_captured, get_vertex_coordinates, BulkMesh, FACE_CUT, num_vertices_per_cell, get_vertex_id, get_vertex_in_out_boundary, FACE_IN, cells_around_vertex!, FACE_OUT, @check

using Test


stl = STL(joinpath(@__DIR__,"data/cube.stl"))
stl = STL(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))

sm = SurfaceMesh(stl)

box = BoundingBox(sm)

#box = expand(box,0.5)
box = expand(box,0.1)

bg_mesh = CartesianMesh(box,20)

bm = BulkMesh(bg_mesh,sm)
writevtk(bm,"bulk")

end # module
