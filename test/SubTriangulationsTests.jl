module SubTriangulationsTests

using Test

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using STLCutters

using STLCutters: Polyhedron
using STLCutters: STL
using STLCutters: restrict
using STLCutters: clip
using STLCutters: split
using STLCutters: surface, volume
using STLCutters: get_facet_planes
using STLCutters: get_reflex_planes
using STLCutters: compute_distances!
using STLCutters: get_original_reflex_faces
using STLCutters: get_cell_planes
using STLCutters: filter_face_planes
using STLCutters: refine
using STLCutters: get_cell_nodes_to_inout
using STLCutters: read_stl
using STLCutters: merge_nodes
using STLCutters: compute_stl_model
using STLCutters: compute_grid
using STLCutters: FACE_IN, FACE_OUT, FACE_CUT

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

stlmodel = compute_stl_model( Table(facet_to_vertices), vertices )
stl = STL(stlmodel)

# Global setup
Πf = get_facet_planes(stl)
Πr = get_reflex_planes(stl,Πf)

# Setup cell
stl_facets_k = 1:num_cells(stl)
K = Polyhedron(HEX)

Γ0 = Polyhedron(stl)
Γk0 = restrict(Γ0,stl,stl_facets_k)

stl_reflex_faces_k = get_original_reflex_faces(Γk0,stl)
stl_facets_k = stl_facets_k.+get_offset(stl,num_dims(stl))

Πk,Πk_ids,Πk_io = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))
Πf,Πf_ids = filter_face_planes(stl,Πr,stl_reflex_faces_k,Πf,stl_facets_k)

compute_distances!(Γk0,lazy_append(Πk,Πf),lazy_append(Πk_ids,Πf_ids))
compute_distances!(K,Πf,Πf_ids)

Γk = clip(Γk0,Πk_ids,inout=Πk_io)

Kn_in = refine(K,Γk,stl,stl_reflex_faces_k,inside=true)
Kn_out = refine(K,Γk,stl,stl_reflex_faces_k,inside=false)

T_in,X_in = simplexify(Kn_in)
T_out,X_out = simplexify(Kn_out)
T_Γ,X_Γ = simplexify(Γk)

n_to_io = get_cell_nodes_to_inout(Kn_in,Kn_out,HEX)

# Test and Write

@test n_to_io == [fill(FACE_OUT,num_vertices(HEX)÷2);fill(FACE_IN,num_vertices(HEX)÷2)]
mesh_in = compute_grid(T_in,X_in,TET)
mesh_out = compute_grid(T_out,X_out,TET)
mesh_Γ = compute_grid(T_Γ,X_Γ,TRI)

#  writevtk(Γk,"Gk")
#  writevtk(Kn_in,"Kin")
#  writevtk(Kn_out,"Kin")
#  writevtk(mesh_in,"mesh_in")
#  writevtk(mesh_out,"mesh_out")
#  writevtk(mesh_Γ,"mesh_G")

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

stlmodel = compute_stl_model( Table(facet_to_vertices), vertices )
stl = STL(stlmodel)

# Global setup
Πf = get_facet_planes(stl)
Πr = get_reflex_planes(stl,Πf)

# Setup cell
stl_facets_k = 1:num_cells(stl)
K = Polyhedron(HEX)

Γ0 = Polyhedron(stl)
Γk0 = restrict(Γ0,stl,stl_facets_k)

stl_reflex_faces_k = get_original_reflex_faces(Γk0,stl)
stl_facets_k = stl_facets_k.+get_offset(stl,num_dims(stl))

Πk,Πk_ids,Πk_io = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))
Πf,Πf_ids = filter_face_planes(stl,Πr,stl_reflex_faces_k,Πf,stl_facets_k)

compute_distances!(Γk0,lazy_append(Πk,Πf),lazy_append(Πk_ids,Πf_ids))
compute_distances!(K,Πf,Πf_ids)

Γk = clip(Γk0,Πk_ids,inout=Πk_io)

Kn_in = refine(K,Γk,stl,stl_reflex_faces_k,inside=true)
Kn_out = refine(K,Γk,stl,stl_reflex_faces_k,inside=false)

T_in,X_in = simplexify(Kn_in)
T_out,X_out = simplexify(Kn_out)
T_Γ,X_Γ = simplexify(Γk)

n_to_io = get_cell_nodes_to_inout(Kn_in,Kn_out,HEX)

# Test and Write

@test n_to_io == fill(FACE_IN,num_vertices(HEX))

mesh_in = compute_grid(T_in,X_in,TET)
mesh_out = compute_grid(T_out,X_out,TET)

#  writevtk(Γk,"Gk")
#  writevtk(Kn_in,"Kin")
#  writevtk(Kn_out,"Kin")
#  writevtk(mesh_in,"mesh_in")
#  writevtk(mesh_out,"mesh_out")
#  writevtk(mesh_Γ,"mesh_G")

@test volume(mesh_in) + volume(mesh_out) ≈ 1

## Real STL

X,T,N = read_stl(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))
stl = compute_stl_model(T,X)
stl = merge_nodes(stl)

#writevtk(get_grid(stl),"stl")

p = HEX
δ = 0.2
n = 10
D = 3

pmin,pmax = get_bounding_box(stl)
Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)

model = CartesianDiscreteModel(pmin,pmax,partition)
grid = get_grid(model)

@time subcells, subfaces, labels = subtriangulate(model,stl,kdtree=false)

#writevtk(submesh,"submesh",cellfields=["inout"=>k_to_io,"bgcell"=>k_to_bgcell])
#writevtk(grid,"bgmesh",cellfields=["inoutcut"=>bgcell_to_ioc])

bgmesh_in_vol = volume(grid,labels.bgcell_to_ioc,:in)
bgmesh_out_vol  = volume(grid,labels.bgcell_to_ioc,:out)
bgmesh_cut_vol  = volume(grid,labels.bgcell_to_ioc,:cut)
submesh_in_vol  = volume(subcells,labels.cell_to_io,:in)
submesh_out_vol = volume(subcells,labels.cell_to_io,:out)

in_volume = bgmesh_in_vol + submesh_in_vol
out_volume = bgmesh_out_vol + submesh_out_vol
cut_volume = bgmesh_cut_vol


#celldata = [ "inoutcut" => labels.bgcell_to_ioc ]
#writevtk( grid, "bgcells"; celldata )
#
#celldata = [ "inout" => labels.cell_to_io, "bgcell" => labels.cell_to_bgcell ]
#writevtk( subcells, "subcells"; celldata )
#
#celldata = [ "bgcell" => labels.face_to_bgcell ]
#writevtk( subfaces, "subfaces" )

println("Num subcells: $(num_cells(subcells))")
@test surface(get_grid(stl)) ≈ surface(subfaces)
@test submesh_in_vol + submesh_out_vol ≈ cut_volume
@test in_volume + out_volume ≈ volume(grid)
@test in_volume ≈ 273280.03374196636

X,T,N = read_stl(joinpath(@__DIR__,"data/wine_glass.stl"))
stl = compute_stl_model(T,X)
stl = merge_nodes(stl)

#writevtk(get_grid(stl),"stl")

p = HEX
δ = 0.2
n = 20
D = 3

pmin,pmax = get_bounding_box(stl)
Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)

model = CartesianDiscreteModel(pmin,pmax,partition)
grid = get_grid(model)

@time subcells, subfaces, labels = subtriangulate(model,stl,kdtree=false)

#writevtk(submesh,"submesh",cellfields=["inout"=>k_to_io,"bgcell"=>k_to_bgcell])
#writevtk(grid,"bgmesh",cellfields=["inoutcut"=>bgcell_to_ioc])

bgmesh_in_vol = volume(grid,labels.bgcell_to_ioc,:in)
bgmesh_out_vol  = volume(grid,labels.bgcell_to_ioc,:out)
bgmesh_cut_vol  = volume(grid,labels.bgcell_to_ioc,:cut)
submesh_in_vol  = volume(subcells,labels.cell_to_io,:in)
submesh_out_vol = volume(subcells,labels.cell_to_io,:out)

in_volume = bgmesh_in_vol + submesh_in_vol
out_volume = bgmesh_out_vol + submesh_out_vol
cut_volume = bgmesh_cut_vol

println("Num subcells: $(num_cells(subcells))")
@test surface(get_grid(stl)) ≈ surface(subfaces)
@test submesh_in_vol + submesh_out_vol ≈ cut_volume
@test in_volume + out_volume ≈ volume(grid)
@test in_volume ≈ 74.12595970214333

# Simplex background grid

X,T,N = read_stl(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))
stl = compute_stl_model(T,X)
stl = merge_nodes(stl)

#writevtk(get_grid(stl),"stl")

p = HEX
δ = 0.2
n = 10
D = 3

pmin,pmax = get_bounding_box(stl)
Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)

model = CartesianDiscreteModel(pmin,pmax,partition)
model = simplexify(model,positive=true)
grid = get_grid(model)

@time subcells, subfaces, labels = subtriangulate(model,stl,kdtree=false)

#writevtk(submesh,"submesh",cellfields=["inout"=>k_to_io,"bgcell"=>k_to_bgcell])
#writevtk(grid,"bgmesh",cellfields=["inoutcut"=>bgcell_to_ioc])

bgmesh_in_vol = volume(grid,labels.bgcell_to_ioc,:in)
bgmesh_out_vol  = volume(grid,labels.bgcell_to_ioc,:out)
bgmesh_cut_vol  = volume(grid,labels.bgcell_to_ioc,:cut)
submesh_in_vol  = volume(subcells,labels.cell_to_io,:in)
submesh_out_vol = volume(subcells,labels.cell_to_io,:out)

in_volume = bgmesh_in_vol + submesh_in_vol
out_volume = bgmesh_out_vol + submesh_out_vol
cut_volume = bgmesh_cut_vol


#celldata = [ "inoutcut" => labels.bgcell_to_ioc ]
#writevtk( grid, "bgcells"; celldata )
#
#celldata = [ "inout" => labels.cell_to_io, "bgcell" => labels.cell_to_bgcell ]
#writevtk( subcells, "subcells"; celldata )
#
#celldata = [ "bgcell" => labels.face_to_bgcell ]
#writevtk( subfaces, "subfaces" )

println("Num subcells: $(num_cells(subcells))")
@test surface(get_grid(stl)) ≈ surface(subfaces)
@test submesh_in_vol + submesh_out_vol ≈ cut_volume
@test in_volume + out_volume ≈ volume(grid)
@test in_volume ≈ 273280.03374196636

X,T,N = read_stl(joinpath(@__DIR__,"data/wine_glass.stl"))
stl = compute_stl_model(T,X)
stl = merge_nodes(stl)

#writevtk(get_grid(stl),"stl")

p = HEX
δ = 0.2
n = 20
D = 3

pmin,pmax = get_bounding_box(stl)
Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)

model = CartesianDiscreteModel(pmin,pmax,partition)
model = simplexify(model,positive=true)
grid = get_grid(model)

@time subcells, subfaces, labels = subtriangulate(model,stl,kdtree=false)

#writevtk(submesh,"submesh",cellfields=["inout"=>k_to_io,"bgcell"=>k_to_bgcell])
#writevtk(grid,"bgmesh",cellfields=["inoutcut"=>bgcell_to_ioc])

bgmesh_in_vol = volume(grid,labels.bgcell_to_ioc,:in)
bgmesh_out_vol  = volume(grid,labels.bgcell_to_ioc,:out)
bgmesh_cut_vol  = volume(grid,labels.bgcell_to_ioc,:cut)
submesh_in_vol  = volume(subcells,labels.cell_to_io,:in)
submesh_out_vol = volume(subcells,labels.cell_to_io,:out)

in_volume = bgmesh_in_vol + submesh_in_vol
out_volume = bgmesh_out_vol + submesh_out_vol
cut_volume = bgmesh_cut_vol

println("Num subcells: $(num_cells(subcells))")
@test surface(get_grid(stl)) ≈ surface(subfaces)
@test submesh_in_vol + submesh_out_vol ≈ cut_volume
@test in_volume + out_volume ≈ volume(grid)
@test in_volume ≈ 74.12595970214333
## Kd-Tree

X,T,N = read_stl(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))
stl = compute_stl_model(T,X)
stl = merge_nodes(stl)

#writevtk(get_grid(stl),"stl")

p = HEX
δ = 0.2
n = 10
D = 3

pmin,pmax = get_bounding_box(stl)
Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)

model = CartesianDiscreteModel(pmin,pmax,partition)
grid = get_grid(model)

@time subcells, subfaces, labels = subtriangulate(model,stl,kdtree=true)

#writevtk(submesh,"submesh",cellfields=["inout"=>k_to_io,"bgcell"=>k_to_bgcell])
#writevtk(grid,"bgmesh",cellfields=["inoutcut"=>bgcell_to_ioc])

bgmesh_in_vol = volume(grid,labels.bgcell_to_ioc,:in)
bgmesh_out_vol  = volume(grid,labels.bgcell_to_ioc,:out)
bgmesh_cut_vol  = volume(grid,labels.bgcell_to_ioc,:cut)
submesh_in_vol  = volume(subcells,labels.cell_to_io,:in)
submesh_out_vol = volume(subcells,labels.cell_to_io,:out)

in_volume = bgmesh_in_vol + submesh_in_vol
out_volume = bgmesh_out_vol + submesh_out_vol
cut_volume = bgmesh_cut_vol

println("Num subcells: $(num_cells(subcells))")
@test surface(get_grid(stl)) ≈ surface(subfaces)
@test submesh_in_vol + submesh_out_vol ≈ cut_volume
@test in_volume + out_volume ≈ volume(grid)
@test in_volume ≈ 273280.03374196636

X,T,N = read_stl(joinpath(@__DIR__,"data/wine_glass.stl"))
stl = compute_stl_model(T,X)
stl = merge_nodes(stl)

#writevtk(get_grid(stl),"stl")

p = HEX
δ = 0.2
n = 20
D = 3

pmin,pmax = get_bounding_box(stl)
Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)

model = CartesianDiscreteModel(pmin,pmax,partition)
grid = get_grid(model)

@time subcells, subfaces, labels = subtriangulate(model,stl,kdtree=true)

#writevtk(submesh,"submesh",cellfields=["inout"=>k_to_io,"bgcell"=>k_to_bgcell])
#writevtk(grid,"bgmesh",cellfields=["inoutcut"=>bgcell_to_ioc])


bgmesh_in_vol = volume(grid,labels.bgcell_to_ioc,:in)
bgmesh_out_vol  = volume(grid,labels.bgcell_to_ioc,:out)
bgmesh_cut_vol  = volume(grid,labels.bgcell_to_ioc,:cut)
submesh_in_vol  = volume(subcells,labels.cell_to_io,:in)
submesh_out_vol = volume(subcells,labels.cell_to_io,:out)

in_volume = bgmesh_in_vol + submesh_in_vol
out_volume = bgmesh_out_vol + submesh_out_vol
cut_volume = bgmesh_cut_vol

println("Num subcells: $(num_cells(subcells))")
@test surface(get_grid(stl)) ≈ surface(subfaces)
@test submesh_in_vol + submesh_out_vol ≈ cut_volume
@test in_volume + out_volume ≈ volume(grid)
@test in_volume ≈ 74.12595970214333

# Test simplexify in cube

X,T,N = read_stl(joinpath(@__DIR__,"data/cube.stl"))
stl = compute_stl_model(T,X)
stl = merge_nodes(stl)

p = HEX
δ = 1.5
n = 8
D = 3

pmin,pmax = get_bounding_box(stl)
Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)

model = CartesianDiscreteModel(pmin,pmax,partition)
model = simplexify(model,positive=true)
grid = get_grid(model)

@time subcells, subfaces, labels = subtriangulate(model,stl)

k_to_io = labels.cell_to_io
k_to_bgcell = labels.cell_to_bgcell
bgcell_to_ioc = labels.bgcell_to_ioc
writevtk(subcells,"submesh",celldata=["inout"=>k_to_io,"bgcell"=>k_to_bgcell])
writevtk(grid,"bgmesh",celldata=["inoutcut"=>bgcell_to_ioc])

bgmesh_in_vol = volume(grid,labels.bgcell_to_ioc,:in)
bgmesh_out_vol  = volume(grid,labels.bgcell_to_ioc,:out)
bgmesh_cut_vol  = volume(grid,labels.bgcell_to_ioc,:cut)
submesh_in_vol  = volume(subcells,labels.cell_to_io,:in)
submesh_out_vol = volume(subcells,labels.cell_to_io,:out)
submesh_surf = surface(subfaces)

in_volume = bgmesh_in_vol + submesh_in_vol
out_volume = bgmesh_out_vol + submesh_out_vol

@test submesh_surf ≈ 6
@test in_volume ≈ 1
@test in_volume + out_volume ≈ volume(grid)


end # module
