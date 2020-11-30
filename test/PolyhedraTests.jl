module PolyhedraTests

using Test

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using STLCutters

using STLCutters: Polyhedron 
using STLCutters: restrict 
using STLCutters: clip 
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
using STLCutters: get_facet_planes
using STLCutters: compute_distances! 
using STLCutters: flip
using STLCutters: get_original_facets
using STLCutters: get_original_reflex_faces
using STLCutters: get_disconnected_faces 

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
stl_topo = get_grid_topology(stl)


Πr = get_reflex_planes(stl)#,bisector=true)
Πf = get_facet_planes(stl)

rf_to_isconvex = get_convex_faces(stl)
stl_k_facets = 1:num_cells(stl)

K = Polyhedron(HEX)

Γ0 = Polyhedron(stl)
Γk0 = restrict(Γ0,stl,stl_k_facets)

Πk = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))
Πk_faces = -(1:length(Πk))

Dc = num_dims(stl)
roffset = get_offset(stl_topo,Dc-1)

stl_rfaces = get_original_reflex_faces(Γk0,stl)
Πrk = view(Πr,stl_rfaces.-roffset)
Πfk = view(Πf,stl_k_facets)
Πrk_faces = stl_rfaces 
Πfk_faces = stl_k_facets .+ get_offset(stl_topo,Dc)

compute_distances!(Γk0,Πk,Πk_faces)
compute_distances!(Γk0,Πrk,Πrk_faces)
compute_distances!(Γk0,Πfk,Πfk_faces)

compute_distances!(K,Πrk,Πrk_faces)
compute_distances!(K,Πfk,Πfk_faces)

writevtk(edge_mesh(Γk0),"Gk0")

Γk = clip(Γk0,Πk_faces)

part_to_facets = get_disconnected_faces(Γk,stl,Dc)

writevtk(edge_mesh(Γk),"Gk")

Γk_out = flip(Γk)

rfaces_in = filter(i->!rf_to_isconvex[i-roffset],Πrk_faces)
rfaces_out = filter(i->rf_to_isconvex[i-roffset],Πrk_faces)


Γks,Ks = decompose(Γk,K,rfaces_in)
Ks_in = empty(Ks)
for (i,(Γki,Ki)) in enumerate(zip(Γks,Ks))
  facets = get_original_facets(Γki,stl)
  #_part_to_facets = ...(part_to_facets,facets)
  Kin = clip(Ki,facets)
  push!(Ks_in,Kin)
end

Γks,Ks = decompose(Γk_out,K,rfaces_out)
Ks_out = empty(Ks)
for (i,(Γki,Ki)) in enumerate(zip(Γks,Ks))
  facets = get_original_facets(Γki,stl)
  Kout = clip(Ki,facets,inside=false)
  push!(Ks_out,Kout)
end

T_in,X_in = simplexify(Ks_in)
T_out,X_out = simplexify(Ks_out)

mesh_in = compute_grid(T_in,X_in,TET)
mesh_out = compute_grid(T_out,X_out,TET)

@test volume(mesh_in) + volume(mesh_out) ≈ 1

for (i,Ki_in) in enumerate(Ks_in)
  writevtk(edge_mesh(Ki_in),"Kin_$i")
end

for (i,Ki_out) in enumerate(Ks_out)
  writevtk(edge_mesh(Ki_out),"Kout_$i")
end

writevtk(mesh_in,"mesh_in")

writevtk(mesh_out,"mesh_out")


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

stl_topo = get_grid_topology(stl)


Πr = get_reflex_planes(stl,bisector=true)
Πf = get_facet_planes(stl)

rf_to_isconvex = get_convex_faces(stl)
stl_k_facets = 1:num_cells(stl)

K = Polyhedron(HEX)

Γ0 = Polyhedron(stl)
Γk0 = restrict(Γ0,stl,stl_k_facets)

Πk = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))
Πk_faces = -(1:length(Πk))

Dc = num_dims(stl)
roffset = get_offset(stl_topo,Dc-1)

stl_rfaces = get_original_reflex_faces(Γk0,stl)
Πrk = view(Πr,stl_rfaces.-roffset)
Πfk = view(Πf,stl_k_facets)
Πrk_faces = stl_rfaces 
Πfk_faces = stl_k_facets .+ get_offset(stl_topo,Dc)

compute_distances!(Γk0,Πk,Πk_faces)
compute_distances!(Γk0,Πrk,Πrk_faces)
compute_distances!(Γk0,Πfk,Πfk_faces)

compute_distances!(K,Πrk,Πrk_faces)
compute_distances!(K,Πfk,Πfk_faces)

writevtk(edge_mesh(Γk0),"Gk0")

Γk = clip(Γk0,Πk_faces)

writevtk(edge_mesh(Γk),"Gk")

Γk_out = flip(Γk)

rfaces_in = filter(i->!rf_to_isconvex[i-roffset],Πrk_faces)
rfaces_out = filter(i->rf_to_isconvex[i-roffset],Πrk_faces)


Γks,Ks = decompose(Γk,K,rfaces_in)
Ks_in = empty(Ks)
for (i,(Γki,Ki)) in enumerate(zip(Γks,Ks))
  facets = get_original_facets(Γki,stl)
  Kin = clip(Ki,facets)
  push!(Ks_in,Kin)
end

Γks,Ks = decompose(Γk_out,K,rfaces_out)
Ks_out = empty(Ks)
for (i,(Γki,Ki)) in enumerate(zip(Γks,Ks))
  facets = get_original_facets(Γki,stl)
  Kout = clip(Ki,facets,inside=false)
  push!(Ks_out,Kout)
end

T_in,X_in = simplexify(Ks_in)
T_out,X_out = simplexify(Ks_out)

mesh_in = compute_grid(T_in,X_in,TET)
mesh_out = compute_grid(T_out,X_out,TET)

@test volume(mesh_in) + volume(mesh_out) ≈ 1

for (i,Ki_in) in enumerate(Ks_in)
  writevtk(edge_mesh(Ki_in),"Kin_$i")
end

for (i,Ki_out) in enumerate(Ks_out)
  writevtk(edge_mesh(Ki_out),"Kout_$i")
end

writevtk(mesh_in,"mesh_in")

writevtk(mesh_out,"mesh_out")



## Real STL

#X,T,N = read_stl(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))
#stl = compute_stl_model(T,X)
#stl = merge_nodes(stl)
#p_stl = Polyhedron(stl)
#
#
#for facet in 1:num_cells(stl)
#  nodes = get_cell_nodes(stl)[facet]
#  @assert STLCutters.next_vertex(p_stl,nodes[1],nodes[2]) == nodes[3]
#end
#
#N = [ normal(get_cell(stl,i)) for i in 1:num_cells(stl) ]
#writevtk(stl.grid,"stl",cellfields=["normals"=>N])
#
#p = HEX
#δ = 0.2
#n = 10
#D = 3
#
#pmin,pmax = get_bounding_box(stl)
#Δ = (pmax-pmin)*δ
#pmin = pmin - Δ
#pmax = pmax + Δ
#partition = (n,n,n)
#
#grid = CartesianGrid(pmin,pmax,partition)
#
#
#writevtk(grid,"bg_grid")
#
#Π_r = get_reflex_planes(stl,bisector=true)
#e_to_isconvex = get_convex_faces(stl)
#
#c_to_stlf = compute_cell_to_facets(grid,stl)
#
#T = Vector{Int}[]
#X = empty(X)
#k_to_io = Int8[]
#cell = 144
#cell = 123
#cell = 127
#cell = 133
##cell = 645
#cell_facets = Int[]
#
##for cell in 1:num_cells(grid)
#@show cell
##!isempty(c_to_stlf[cell]) || continue
#
#pmin = get_cell_coordinates(grid)[cell][1]
#pmax = get_cell_coordinates(grid)[cell][end]
#
#empty!(cell_facets)
#for f in c_to_stlf[cell] 
#  for v in get_faces(get_grid_topology(stl),D-1,0)[f]
#    for _f in get_faces(get_grid_topology(stl),0,D-1)[v]
#      if _f ∉ cell_facets
#        push!(cell_facets,_f)
#      end
#    end
#  end
#end
#
#Γ0 = restrict(p_stl,stl,cell_facets)
#
#Πc = get_cell_planes(p,pmin,pmax)
#cell = Polyhedron(p,get_cell_coordinates(grid)[cell])
#
#Γ = clip(Γ0,Πc)
##!isnothing(Γ) || continue
#
##writevtk(edge_mesh(Γ0),"Gamma_0")
##writevtk(edge_mesh(Γ),"Gamma")
##writevtk(edge_mesh(cell),"cell")
#
#
#p⁻,p⁺ = merge(p,cell,Γ)
#
#writevtk(edge_mesh(p⁻),"poly_in")
#writevtk(edge_mesh(p⁺),"poly_out")
#
#
#in_polys = decompose(p⁻,Π_r,e_to_isconvex,isinside=true)
#out_polys = decompose(p⁺,Π_r,e_to_isconvex,isinside=false)
#
#Tin,Xin = simplexify(in_polys)
#Tout,Xout = simplexify(out_polys)
#
#
#append!(T, map(i->i.+length(X),Tin) )
#append!(X,Xin)
#append!(k_to_io,fill(FACE_IN,length(Tin)))
#append!(T, map(i->i.+length(X),Tout) )
#append!(X,Xout)
#append!(k_to_io,fill(FACE_OUT,length(Tout)))
#
#
##mesh_in = compute_grid(Table(Tin),Xin,TET)
##mesh_out = compute_grid(Table(Tout),Xout,TET)
##
##writevtk(edge_mesh(Γ0),"Gamma_0")
##writevtk(edge_mesh(Γ),"Gamma")
##
##writevtk(edge_mesh(p⁻),"poly_in")
##writevtk(edge_mesh(p⁺),"poly_out")
##
##for (i,poly) in enumerate(out_polys)
##  writevtk(edge_mesh(poly),"poly_out_$i")
##end
##
##for (i,poly) in enumerate(in_polys)
##  writevtk(edge_mesh(poly),"poly_in_$i")
##end
##
##@show length(out_polys)
##writevtk(mesh_in,"simplices_in")
##writevtk(mesh_out,"simplices_out")
#
##end
#
#submesh = compute_grid(Table(T),X,TET)
#writevtk(submesh,"submesh",cellfields=["inout"=>k_to_io])


end # module
