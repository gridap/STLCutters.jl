module CellMeshesTests

using STLCutter

# To be public
using STLCutter: reset!, compact!, compute_in_out!, num_dfaces, get_faces, dface_dimension, local_dface, global_dface, get_vertex, get_face, get_dface

# Private used in test
using STLCutter: _add_vertex!, UNSET, are_all_faces_defined, add_surface_mesh_face!, find_closest_face, expand

using Test

p1 = Point(1,2,3)
p2 = Point(4,5,6)

b = BoundingBox(p1,p2)

m = CellMesh(b)


p1 = Point(1,2,4)
p2 = Point(5,3,5)

box = BoundingBox(p1,p2)

reset!(m,box)

@test @allocated(reset!(m,box)) == 0


## 2D Point Intersection

p0 = Point(0.0,0.0)
p1 = Point(1.0,1.0)
box = BoundingBox(p0,p1)

mesh = CellMesh(box)

stl_points = [ Point(0.3,0.3), Point(0.25,0.5), Point(0.5,0.5), Point(0.75,0.4), Point(0.0,0.5)  ]

for point in stl_points[1:1]

  d,face = find_closest_face(mesh,point)

  if face != UNSET 
    _add_vertex!(mesh,d,face,point)
  end

end

compact!(mesh)

writevtk(mesh,"sub_mesh")

## 3D Point Intersection

p0 = Point(0.0,0.0,0.0)
p1 = Point(1.0,1.0,1.0)
box = BoundingBox(p0,p1)

mesh = CellMesh(box)

stl_points = [ 
  Point(0.3,0.3,0.0), 
  Point(0.25,0.5,0.1), 
  Point(0.4,0.6,0.0), 
  Point(0.75,0.4,0.2), 
  Point(0.0,0.5,0.0)  ]


#stl_points = [
#  Point((0.375, 1.0, 0.125)),
#  Point((0.75, 0.5, 0.25)),
#  Point((0.8999999999999999, 0.30000000000000004, 0.3)),
#  Point((1.0, 0.16666666666666674, 0.3333333333333333)),
#  Point((1.0, 0.0, 0.3749999999999999)),
#  Point((0.8999999999999999, 0.30000000000000004, 0.3)),
#  Point((1.0, 0.16666666666666674, 0.3333333333333333)),
#  Point((0.7000000000000001, 0.29999999999999993, 0.29999999999999993)) ]

for point in stl_points

  d,face = find_closest_face(mesh,point)
  
  if face != UNSET 
    _add_vertex!(mesh,d,face,point)
  end

end

compact!(mesh)

writevtk(mesh,"sub_mesh3")

## 2D Surface Mesh

points = [ Point(-0.1,0.3), Point(0.2,0.5), Point(0.4,0.4), Point(1.,0.), Point(.8,-.2)  ]

facet_to_vertices = Table( [
  1 2 3 4 5;
  2 3 4 5 1] )

sm = SurfaceMesh(points,facet_to_vertices)

p0 = Point(0.0,0.0)
p1 = Point(1.0,1.0)
box = BoundingBox(p0,p1)

bg_mesh = CartesianMesh(box,1)
cell_mesh = CellMesh(box)

c_to_sm_f = compute_cell_to_surface_mesh_faces(bg_mesh,sm)

cell_id = 1

reset!(cell_mesh,get_cell(bg_mesh,1))

for i in 1:length(c_to_sm_f,cell_id)
  sm_face = c_to_sm_f[cell_id,i]
  add_surface_mesh_face!(cell_mesh,sm,sm_face)
end

compact!(cell_mesh)

compute_in_out!(cell_mesh,sm)

writevtk(sm,"sm2")
writevtk(cell_mesh,"cell_mesh")

@test are_all_faces_defined(cell_mesh)

## 3D Surface Mesh

points = Point{3,Float64}[
  (-0.5,-0.5,.5),
  (-0.,1.5,0.),
  (1.5,-0.5,.5),
  (1.5,1.5,.5)]

facet_to_vertices = Table( [
  1 2;
  2 4;
  3 3 ] )


sm = SurfaceMesh(points,facet_to_vertices)

p0 = Point(0.,0.,0.)
p1 = Point(1.,1.,1.)
box = BoundingBox(p0,p1)

bg_mesh = CartesianMesh(box,1)
cell_mesh = CellMesh(box)

c_to_sm_f = compute_cell_to_surface_mesh_faces(bg_mesh,sm)

cell_id = 1

reset!(cell_mesh,get_cell(bg_mesh,1))

for i in 1:length(c_to_sm_f,cell_id)
  sm_face = c_to_sm_f[cell_id,i]
  add_surface_mesh_face!(cell_mesh,sm,sm_face)
end

compact!(cell_mesh)

compute_in_out!(cell_mesh,sm)
writevtk(sm,"sm3")
writevtk(cell_mesh,"cell_mesh3")

## Big 3D facet

points = Point{3,Float64}[
  (-1.0,-1.0,0.5),
  (4.0,-1.0,0.0),
  (-1.0,4.0,0.5), 
  (-2.0,-1.0,0.5)]

facet_to_vertices = Table( [
  1 1;
  2 3;
  3 4] )

sm = SurfaceMesh(points,facet_to_vertices)

c_to_sm_f = compute_cell_to_surface_mesh_faces(bg_mesh,sm)

reset!(cell_mesh,get_cell(bg_mesh,1))

for i in 1:length(c_to_sm_f,cell_id)
  sm_face = c_to_sm_f[cell_id,i]
  add_surface_mesh_face!(cell_mesh,sm,sm_face)
end

compact!(cell_mesh)

compute_in_out!(cell_mesh,sm)
writevtk(sm,"sm3")
writevtk(cell_mesh,"cell_mesh3")

# Real STL

stl = STL(joinpath(@__DIR__,"data/cube.stl"))

sm = SurfaceMesh(stl)

box = BoundingBox(sm)

box = expand(box,0.5)

bg_mesh = CartesianMesh(box,1)

c_to_sm_f = compute_cell_to_surface_mesh_faces(bg_mesh,sm)

cell_id = 1

reset!(cell_mesh,get_cell(bg_mesh,1))


c_to_sm_f = compute_cell_to_surface_mesh_faces(bg_mesh,sm)

reset!(cell_mesh,get_cell(bg_mesh,1))

for i in 1:length(c_to_sm_f,cell_id)
  sm_face = c_to_sm_f[cell_id,i]
  add_surface_mesh_face!(cell_mesh,sm,sm_face)
end

compact!(cell_mesh)

compute_in_out!(cell_mesh,sm)
writevtk(sm,"sm3")
writevtk(cell_mesh,"cell_mesh3")


end # module

