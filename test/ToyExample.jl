module ToyExamples

using Test

using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs
using STLCutters


using STLCutters: Polyhedron

using STLCutters: compute_stl_model 
using STLCutters: get_facet_planes 
using STLCutters: get_reflex_planes 
using STLCutters: restrict, clip, refine, decompose, split 
using STLCutters: get_original_reflex_faces 
using STLCutters: get_original_facets 
using STLCutters: get_cell_planes 
using STLCutters: filter_face_planes 
using STLCutters: compute_distances! 
using STLCutters: volume
using STLCutters: writevtk 


vertices = [
  Point(-0.3,-0.5,0.2),
  Point(0.3,-0.5,0.7),
  Point(0.7,-0.5,0.2),
  Point(1.2,-0.5,0.3),
  Point(1.5,-0.5,0.8),
  Point(-0.3,1.5,0.2),
  Point(0.7,1.5,0.3),
  Point(1.5,1.5,0.2) ]

facet_to_vertices = [
 [1,2,6],
 [2,3,7],
 [3,4,7],
 [4,5,7],
# [3,5,8],
 [2,7,6],
 [5,8,7]]

stl = compute_stl_model( Table(facet_to_vertices), vertices )

Πr = get_reflex_planes(stl,bisector=true)
Πf = get_facet_planes(stl)


K = Polyhedron(HEX)
stl_facets_k = 1:num_cells(stl)

Γ0 = Polyhedron(stl)
Γk0 = restrict(Γ0,stl,stl_facets_k)


stl_reflex_faces_k = get_original_reflex_faces(Γk0,stl)
stl_facets_k = stl_facets_k.+get_offset(get_grid_topology(stl),num_dims(stl))

Πk,Πk_ids,Πk_io = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))
Πf,Πf_ids = filter_face_planes(stl,Πr,stl_reflex_faces_k,Πf,stl_facets_k)

compute_distances!(Γk0,(Πk,Πf),(Πk_ids,Πf_ids))
compute_distances!(K,Πf,Πf_ids)

Γk = clip(Γk0,Πk_ids,inout=Πk_io)

Kn_in = refine(K,Γk,stl,stl_reflex_faces_k,[],inside=true)
Kn_out = refine(K,Γk,stl,stl_reflex_faces_k,[],inside=false)
#
T_in,X_in = simplexify(Kn_in)
T_out,X_out = simplexify(Kn_out)

mesh_in = compute_grid(T_in,X_in,TET)
mesh_out = compute_grid(T_out,X_out,TET)

reflex_faces = filter(f->STLCutters.is_reflex(Γk,stl,f,inside=true),stl_reflex_faces_k)
Γn,Kn = decompose(Γk,K,reflex_faces)

poly = Γk0
Γkn = typeof(poly)[]
push!(Γkn,poly)
for (i,Πid) in enumerate(Πk_ids)
  global poly
  poly0, poly1 = split(poly,Πid)
  poly = Πk_io[i] ? poly0 : poly1
  push!(Γkn,poly)
end

Knn = Vector{typeof(K)}[]
for (i,(Γi,Ki)) in enumerate(zip(Γn,Kn))
  F = get_original_facets(Γi,stl)
  poly = Ki
  push!(Knn,[poly])
  for f in F
    poly0, poly1 = split(poly,f)
    poly = poly0
    push!(Knn[i],poly)
  end
end

println("IN + OUT - 1 = $(volume(mesh_in)+volume(mesh_out)-1)")


testdir = "toy_example"
isdir(testdir) || mkdir(testdir)
writevtk(stl,joinpath(testdir,"stl"))
writevtk(K,joinpath(testdir,"K"))

writevtk(Γk0,joinpath(testdir,"Gk0"))

writevtk(Γk,joinpath(testdir,"Gk"))

writevtk(Γkn,joinpath(testdir,"Gk_"))


for (i,Ki_in) in enumerate(Kn_in)
  writevtk(Ki_in,joinpath(testdir,"Kin_$i"))
end

for (i,Ki_out) in enumerate(Kn_out)
  writevtk(Ki_out,joinpath(testdir,"Kout_$i"))
end

writevtk(mesh_in,joinpath(testdir,"mesh_in"))

writevtk(mesh_out,joinpath(testdir,"mesh_out"))

for (i,kni) in enumerate(Knn)
  writevtk(kni,joinpath(testdir,"Knn_$(i)_"))
end


offsets = [ Point(0.0,0.0,0.0), Point(0.3,0.0,0.15), Point(0.0,0.0,-0.3) ]
for (i,(Γi,Ki,Kin)) in enumerate(zip(Γn,Kn,Kn_in))
  map!(p->p+offsets[i], get_vertex_coordinates(Γi), get_vertex_coordinates(Γi) )
  map!(p->p+offsets[i], get_vertex_coordinates(Ki), get_vertex_coordinates(Ki) )
#  map!(p->p+offsets[i], get_vertex_coordinates(Kin), get_vertex_coordinates(Kin) )
  writevtk(Γi,joinpath(testdir,"G_$i"))
  writevtk(Ki,joinpath(testdir,"K_$i"))
  ## Illustrate clipping process
#  writevtk(Kin,joinpath(testdir,"Ki_$i"))
end


end # module

