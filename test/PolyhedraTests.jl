module PolyhedraTests

using Test

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using STLCutters

using STLCutters: Polyhedron 
using STLCutters: restrict 
using STLCutters: clip 
using STLCutters: decompose 
using STLCutters: edge_mesh
using STLCutters: volume, volumes 
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
using STLCutters: group_facing_facets 
using STLCutters: filter_face_planes 
using STLCutters: refine 
using STLCutters: FACE_CUT 
using STLCutters: get_cell 
using STLCutters: get_cell_nodes_to_inout 

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

# Global setup
Πr = get_reflex_planes(stl)#,bisector=true)
Πf = get_facet_planes(stl)

reflex_face_to_isconvex = get_convex_faces(stl)

# Setup cell
stl_facets_k = 1:num_cells(stl)
K = Polyhedron(HEX)

Γ0 = Polyhedron(stl)
Γk0 = restrict(Γ0,stl,stl_facets_k)

stl_reflex_faces_k = get_original_reflex_faces(Γk0,stl)
stl_facets_k = stl_facets_k.+get_offset(get_grid_topology(stl),num_dims(stl))

Πk,Πk_ids,Πk_io = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))
Πf,Πf_ids = filter_face_planes(stl,Πr,stl_reflex_faces_k,Πf,stl_facets_k)

compute_distances!(Γk0,(Πk,Πf),(Πk_ids,Πf_ids))
compute_distances!(K,Πf,Πf_ids)

Γk = clip(Γk0,Πk_ids,inout=Πk_io)

Kn_in = refine(K,Γk,stl,stl_reflex_faces_k,reflex_face_to_isconvex,inside=true)
Kn_out = refine(K,Γk,stl,stl_reflex_faces_k,reflex_face_to_isconvex,inside=false)

T_in,X_in = simplexify(Kn_in)
T_out,X_out = simplexify(Kn_out)
T_Γ,X_Γ = simplexify(Γk)

n_to_io = get_cell_nodes_to_inout(Kn_in,Kn_out,HEX)

# Test and Write

@test n_to_io == [fill(FACE_OUT,num_vertices(HEX)÷2);fill(FACE_IN,num_vertices(HEX)÷2)]
mesh_in = compute_grid(T_in,X_in,TET)
mesh_out = compute_grid(T_out,X_out,TET)
mesh_Γ = compute_grid(T_Γ,X_Γ,TRI)  

writevtk(edge_mesh(Γk),"Gk")

for (i,Ki_in) in enumerate(Kn_in)
  writevtk(edge_mesh(Ki_in),"Kin_$i")
end

for (i,Ki_out) in enumerate(Kn_out)
  writevtk(edge_mesh(Ki_out),"Kout_$i")
end

writevtk(mesh_in,"mesh_in")

writevtk(mesh_out,"mesh_out")

writevtk(mesh_Γ,"mesh_G")

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
stl_topo = get_grid_topology(stl)

# Global setup
Πr = get_reflex_planes(stl)#,bisector=true)
Πf = get_facet_planes(stl)

reflex_face_to_isconvex = get_convex_faces(stl)

# Setup cell
stl_facets_k = 1:num_cells(stl)
K = Polyhedron(HEX)

Γ0 = Polyhedron(stl)
Γk0 = restrict(Γ0,stl,stl_facets_k)

stl_reflex_faces_k = get_original_reflex_faces(Γk0,stl)
stl_facets_k = stl_facets_k.+get_offset(get_grid_topology(stl),num_dims(stl))

Πk,Πk_ids,Πk_io = get_cell_planes(HEX,Point(0,0,0),Point(1,1,1))
Πf,Πf_ids = filter_face_planes(stl,Πr,stl_reflex_faces_k,Πf,stl_facets_k)

compute_distances!(Γk0,(Πk,Πf),(Πk_ids,Πf_ids))
compute_distances!(K,Πf,Πf_ids)

Γk = clip(Γk0,Πk_ids,inout=Πk_io)

Kn_in = refine(K,Γk,stl,stl_reflex_faces_k,reflex_face_to_isconvex,inside=true)
Kn_out = refine(K,Γk,stl,stl_reflex_faces_k,reflex_face_to_isconvex,inside=false)

T_in,X_in = simplexify(Kn_in)
T_out,X_out = simplexify(Kn_out)
T_Γ,X_Γ = simplexify(Γk)

n_to_io = get_cell_nodes_to_inout(Kn_in,Kn_out,HEX)

# Test and Write

@test n_to_io == fill(FACE_IN,num_vertices(HEX))

mesh_in = compute_grid(T_in,X_in,TET)
mesh_out = compute_grid(T_out,X_out,TET)

writevtk(edge_mesh(Γk),"Gk")

for (i,Ki_in) in enumerate(Kn_in)
  writevtk(edge_mesh(Ki_in),"Kin_$i")
end

for (i,Ki_out) in enumerate(Kn_out)
  writevtk(edge_mesh(Ki_out),"Kout_$i")
end

writevtk(mesh_in,"mesh_in")

writevtk(mesh_out,"mesh_out")


@test volume(mesh_in) + volume(mesh_out) ≈ 1

## Real STL

X,T,N = read_stl(joinpath(@__DIR__,"data/Bunny-LowPoly.stl"))
X,T,N = read_stl(joinpath(@__DIR__,"data/wine_glass.stl"))
stl = compute_stl_model(T,X)
stl = merge_nodes(stl)
Γ0 = Polyhedron(stl)

writevtk(get_grid(stl),"stl")


p = HEX
δ = 0.2
n = 20
D = 3

pmin,pmax = get_bounding_box(stl)
Δ = (pmax-pmin)*δ
pmin = pmin - Δ
pmax = pmax + Δ
partition = (n,n,n)

grid = CartesianGrid(pmin,pmax,partition)


writevtk(grid,"bg_grid")
Πr = get_reflex_planes(stl,bisector=true)
Πf = get_facet_planes(stl)

reflex_face_to_isconvex = get_convex_faces(stl)

c_to_stlf = compute_cell_to_facets(grid,stl)

T = Vector{Int}[]
F = Vector{Int}[]
X = empty(X)
Xf = empty(X)
k_to_io = Int8[]
k_to_bgcell = Int[]
f_to_bgcell = Int[]

cell_facets = Int[]
bgcell_to_ioc = fill(UNSET,num_cells(grid))
bgnode_to_io = fill(UNSET,num_nodes(grid))

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

  K = Polyhedron(p,get_cell_coordinates(grid)[cell])
  Γk0 = restrict(Γ0,stl,cell_facets)


  stl_reflex_faces_k = get_original_reflex_faces(Γk0,stl)
  stl_facets_k = c_to_stlf[cell].+get_offset(get_grid_topology(stl),num_dims(stl))

  Πk,Πk_ids,Πk_io = get_cell_planes(p,pmin,pmax)
  Πkf,Πkf_ids = filter_face_planes(stl,Πr,stl_reflex_faces_k,Πf,stl_facets_k)

  compute_distances!(Γk0,(Πk,Πkf),(Πk_ids,Πkf_ids))
  compute_distances!(K,Πkf,Πkf_ids)

  Γk = clip(Γk0,Πk_ids,inout=Πk_io)
  !isnothing(Γk) || continue

  Kn_in = refine(K,Γk,stl,stl_reflex_faces_k,reflex_face_to_isconvex,inside=true)
  Kn_out = refine(K,Γk,stl,stl_reflex_faces_k,reflex_face_to_isconvex,inside=false)

  Tin,Xin = simplexify(Kn_in)
  Tout,Xout = simplexify(Kn_out)
  T_Γ,X_Γ = simplexify(Γk)
  
  n_to_io = get_cell_nodes_to_inout(Kn_in,Kn_out,p)

  append!(T, map(i->i.+length(X),Tin) )
  append!(X,Xin)
  append!(k_to_io,fill(FACE_IN,length(Tin)))
  append!(k_to_bgcell,fill(cell,length(Tin)))
  append!(T, map(i->i.+length(X),Tout) )
  append!(X,Xout)
  append!(k_to_io,fill(FACE_OUT,length(Tout)))
  append!(k_to_bgcell,fill(cell,length(Tout)))
  append!(F, map(i->i.+length(Xf),T_Γ) ) 
  append!(Xf,X_Γ)
  append!(f_to_bgcell,fill(cell,length(T_Γ)))

  bgcell_to_ioc[cell] = FACE_CUT
  for (i,bg_v) in enumerate(get_cell_nodes(grid)[cell])
    @assert bgnode_to_io[bg_v] ∈ (UNSET,n_to_io[i])
    bgnode_to_io[bg_v] = n_to_io[i]
  end
end

grid_topology = GridTopology(grid)
stack = Int[]
for cell in 1:num_cells(grid)
  if bgcell_to_ioc[cell] == FACE_CUT
    resize!(stack,0)
    push!(stack,cell)
    while length(stack) > 0
      current_cell = pop!(stack)
      for node in get_cell_nodes(grid)[current_cell]
        if bgnode_to_io[node] == FACE_IN
          for neig_cell in get_faces(grid_topology,0,num_dims(grid))[node]
            if bgcell_to_ioc[neig_cell] == UNSET
              bgcell_to_ioc[neig_cell] = FACE_IN
              for neig_node in get_cell_nodes(grid)[neig_cell]
                if bgnode_to_io[neig_node] == UNSET
                  bgnode_to_io[neig_node] = FACE_IN
                end
                push!(stack,neig_cell)
              end
            end
          end
        end
      end
    end
  end
end
for cell in 1:num_cells(grid)
  if bgcell_to_ioc[cell] == UNSET
    bgcell_to_ioc[cell] = FACE_OUT
  end
end

submesh = compute_grid(Table(T),X,TET)
facets = compute_grid(Table(F),Xf,TRI)

writevtk(submesh,"submesh",cellfields=["inout"=>k_to_io,"bgcell"=>k_to_bgcell])
writevtk(grid,"bgmesh",cellfields=["inoutcut"=>bgcell_to_ioc])

bgmesh_vols = volumes(grid)
submesh_vols = volumes(submesh)

bgmesh_in_vols = bgmesh_vols[findall(isequal(FACE_IN),bgcell_to_ioc)]
bgmesh_out_vols = bgmesh_vols[findall(isequal(FACE_OUT),bgcell_to_ioc)]
bgmesh_cut_vols = bgmesh_vols[findall(isequal(FACE_CUT),bgcell_to_ioc)]
submesh_in_vols = submesh_vols[findall(isequal(FACE_IN),k_to_io)]
submesh_out_vols = submesh_vols[findall(isequal(FACE_OUT),k_to_io)]

in_volume = sum(bgmesh_in_vols) + sum(submesh_in_vols)
out_volume = sum(bgmesh_out_vols) + sum(submesh_out_vols)
cut_volume = sum(bgmesh_cut_vols)

@test volume(submesh) ≈ cut_volume 
@test surface(get_grid(stl)) ≈ surface(facets)
@test in_volume + out_volume ≈ volume(grid)
@show in_volume
@show in_volume-74.12595970214333 # Glass
@show in_volume-273280.03374196636 # Bunny

end # module
