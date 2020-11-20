module PolyhedraTests

using Test

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using STLCutters

using STLCutters: Polyhedron 
using STLCutters: restrict 
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
using STLCutters: compute_cell_to_facets
using STLCutters: read_stl, merge_nodes 
using STLCutters: get_bounding_box 
using STLCutters: FACE_IN, FACE_OUT

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

using STLCutters: get_cell, normal
stl = compute_stl_model( Table(facet_to_vertices), vertices )
N = [ normal(get_cell(stl,i)) for i in 1:num_cells(stl) ]
Π_r = get_reflex_planes(stl)
e_to_isconvex = get_convex_faces(stl)

writevtk(stl.grid,"_stl",cellfields=["N"=>N])
cell_facets = 1:num_cells(stl)
pmin,pmax = Point(0,0,0),Point(1,1,1)

Γ0 = Polyhedron(stl)
Γ0 = restrict(Γ0,stl,cell_facets)

Πc = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))
cell = Polyhedron(HEX)

Γ = clip(Γ0,Πc)

writevtk(edge_mesh(Γ0),"Gamma_0")
writevtk(edge_mesh(Γ),"Gamma")

p⁻,p⁺ = merge(HEX,cell,Γ)

in_polys = decompose(p⁻,Π_r,e_to_isconvex,isinside=true)
out_polys = decompose(p⁺,Π_r,e_to_isconvex,isinside=false)

writevtk(edge_mesh(p⁻),"poly_in")
writevtk(edge_mesh(p⁺),"poly_out")

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

Γ0 = Polyhedron(stl)

Πc = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))
cell = Polyhedron(HEX)

Γ = clip(Γ0,Πc)

p⁻,p⁺ = merge(HEX,cell,Γ)

in_polys = decompose(p⁻,Π_r,e_to_isconvex,isinside=true)
out_polys = decompose(p⁺,Π_r,e_to_isconvex,isinside=false)


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

## Real STL

X,T,N = read_stl(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))
stl = compute_stl_model(T,X)
stl = merge_nodes(stl)
p_stl = Polyhedron(stl)


for facet in 1:num_cells(stl)
  nodes = get_cell_nodes(stl)[facet]
  @assert STLCutters.next_vertex(p_stl,nodes[1],nodes[2]) == nodes[3]
end

N = [ normal(get_cell(stl,i)) for i in 1:num_cells(stl) ]
writevtk(stl.grid,"stl",cellfields=["normals"=>N])

p = HEX
δ = 0.2
n = 10
D = 3

pmin,pmax = get_bounding_box(stl)
Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)

grid = CartesianGrid(pmin,pmax,partition)


writevtk(grid,"bg_grid")

Π_r = get_reflex_planes(stl)
e_to_isconvex = get_convex_faces(stl)

c_to_stlf = compute_cell_to_facets(grid,stl)

T = Vector{Int}[]
X = empty(X)
k_to_io = Int8[]
cell = 144
cell = 123
cell = 645
cell_facets = Int[]

for cell in 1:num_cells(grid)

!isempty(c_to_stlf[cell]) || continue

pmin = get_cell_coordinates(grid)[cell][1]
pmax = get_cell_coordinates(grid)[cell][end]

empty!(cell_facets)
for f in c_to_stlf[cell] 
  for v in get_faces(get_grid_topology(stl),D-1,0)[f]
    for _f in get_faces(get_grid_topology(stl),0,D-1)[v]
      if _f ∉ cell_facets
        push!(cell_facets,_f)
      end
    end
  end
end

Γ0 = restrict(p_stl,stl,cell_facets)

Πc = get_cell_planes(p,pmin,pmax)
cell = Polyhedron(p,get_cell_coordinates(grid)[cell])

Γ = clip(Γ0,Πc)
!isnothing(Γ) || continue

#writevtk(edge_mesh(Γ0),"Gamma_0")
#writevtk(edge_mesh(Γ),"Gamma")
#writevtk(edge_mesh(cell),"cell")


p⁻,p⁺ = merge(p,cell,Γ)

#writevtk(edge_mesh(p⁻),"poly_in")
#writevtk(edge_mesh(p⁺),"poly_out")


in_polys = decompose(p⁻,Π_r,e_to_isconvex,isinside=true)
out_polys = decompose(p⁺,Π_r,e_to_isconvex,isinside=false)

Tin,Xin = simplexify(in_polys)
Tout,Xout = simplexify(out_polys)


append!(T, map(i->i.+length(X),Tin) )
append!(X,Xin)
append!(k_to_io,fill(FACE_IN,length(Tin)))
append!(T, map(i->i.+length(X),Tout) )
append!(X,Xout)
append!(k_to_io,fill(FACE_OUT,length(Tout)))


mesh_in = compute_grid(Table(Tin),Xin,TET)
mesh_out = compute_grid(Table(Tout),Xout,TET)
#
#writevtk(edge_mesh(Γ0),"Gamma_0")
#writevtk(edge_mesh(Γ),"Gamma")
#
#writevtk(edge_mesh(p⁻),"poly_in")
#writevtk(edge_mesh(p⁺),"poly_out")
#
#for (i,poly) in enumerate(out_polys)
#  writevtk(edge_mesh(poly),"poly_out_$i")
#end
#
#for (i,poly) in enumerate(in_polys)
#  writevtk(edge_mesh(poly),"poly_in_$i")
#end
#
#writevtk(mesh_in,"simplices_in")
#writevtk(mesh_out,"simplices_out")

end

submesh = compute_grid(Table(T),X,TET)
writevtk(submesh,"submesh",cellfields=["inout"=>k_to_io])

#TODO: 
#
#  Debug:
#    - Error in simplexity, maybe due to a wrong decomposition, 
#      e.g., cell 123 bunny, check paraview
#
#  Workflow:
#    - embedd with global workflow
#    - gridap FE problem
#  
#  Improvement:
#    - force dist to be zero if it shoul be, do not reflex for zero length edge
#

end # module
