module PolyhedraTests

using Test

using Gridap
using STLCutters

using STLCutters: Polyhedron 
using STLCutters: clip 
using STLCutters: merge 
using STLCutters: decompose 
using STLCutters: edge_mesh
using STLCutters: volume 
using STLCutters: compute_stl_model 
using STLCutters: compute_grid 
using STLCutters: get_cell_planes 
using STLCutters: get_reflex_planes 
using STLCutters: get_convex_faces 

vertices = [
  Point(0.1,-0.2,0.5),
  Point(1.0,-0.3,0.5),
  Point(-0.2,1.2,0.5),
  Point(-0.3,0.3,0.5),
  Point(0.5,1.3,0.5),
  Point(0.5,0.5,0.7),
  Point(1.2,1.2,0.5)]

facet_to_vertices = 
[[4,6,1],
 [6,5,7],
 [4,3,6],
 [1,6,2],
 [3,5,6],
 [6,7,2] ]

stl = compute_stl_model( Table(facet_to_vertices), vertices )
Π_r = get_reflex_planes(stl)
e_to_isconvex = get_convex_faces(stl)

cell_facets = 1:num_cells(stl)
pmin,pmax = Point(0,0,0),Point(1,1,1)

Γ0 = Polyhedron(stl,cell_facets)

Πc = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))
cell = Polyhedron(HEX)

Γ = clip(Γ0,Πc)

p⁻,p⁺ = merge(HEX,cell,Γ)

in_polys = decompose(p⁻,Π_r,e_to_isconvex,isinside=true)
out_polys = decompose(p⁺,Π_r,e_to_isconvex,isinside=false)

Tin,Xin = simplexify(in_polys)
Tout,Xout = simplexify(out_polys)

mesh_in = compute_grid(Table(Tin),Xin,TET)
mesh_out = compute_grid(Table(Tout),Xout,TET)

writevtk(edge_mesh(Γ0),"Gamma_0")
writevtk(edge_mesh(Γ),"Gamma")

writevtk(edge_mesh(p⁻),"poly_in")
writevtk(edge_mesh(p⁺),"poly_out")


for (i,poly) in enumerate(out_polys)
  writevtk(edge_mesh(poly),"poly_out_$i")
end

for (i,poly) in enumerate(in_polys)
  writevtk(edge_mesh(poly),"poly_in_$i")
end

writevtk(mesh_in,"simplices_in")
writevtk(mesh_out,"simplices_out")

@test volume(mesh_in) + volume(mesh_out) ≈ 1

# Sharp corner

vertices = [
  Point(0.1,-0.2,-1.0),
  Point(1.0,-0.3,-1.0),
  Point(-0.2,1.2,-1.0),
  Point(-0.3,0.3,-1.0),
  Point(0.5,1.3,-1.0),
  Point(0.5,0.5,0.5),
  Point(1.2,1.2,-1.0)]

facet_to_vertices = 
[[4,6,1],
 [6,5,7],
 [4,3,6],
 [1,6,2],
 [3,5,6],
 [6,7,2] ]

stl = compute_stl_model( Table(facet_to_vertices), vertices )
Π_r = get_reflex_planes(stl)
e_to_isconvex = get_convex_faces(stl)

cell_facets = 1:num_cells(stl)
pmin,pmax = Point(0,0,0),Point(1,1,1)

Γ0 = Polyhedron(stl,cell_facets)

Πc = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))
cell = Polyhedron(HEX)

Γ = clip(Γ0,Πc)

p⁻,p⁺ = merge(HEX,cell,Γ)

in_polys = decompose(p⁻,Π_r,e_to_isconvex,isinside=true)
out_polys = decompose(p⁺,Π_r,e_to_isconvex,isinside=false)

Tin,Xin = simplexify(in_polys)
Tout,Xout = simplexify(out_polys)

mesh_in = compute_grid(Table(Tin),Xin,TET)
mesh_out = compute_grid(Table(Tout),Xout,TET)

writevtk(edge_mesh(Γ0),"Gamma_0")
writevtk(edge_mesh(Γ),"Gamma")

writevtk(edge_mesh(p⁻),"poly_in")
writevtk(edge_mesh(p⁺),"poly_out")

for (i,poly) in enumerate(out_polys)
  writevtk(edge_mesh(poly),"poly_out_$i")
end

for (i,poly) in enumerate(in_polys)
  writevtk(edge_mesh(poly),"poly_in_$i")
end

writevtk(mesh_in,"simplices_in")
writevtk(mesh_out,"simplices_out")

@test volume(mesh_in) + volume(mesh_out) ≈ 1


#TODO: 
#
#  Workflow:
#    - filter faces with voxels
#    - embedd with global workflow
#    - gridap FE problem
#  
#  Improvement:
#    - force dist to be zero if it shoul be, do not reflex for zero length edge
#

end # module
