module RefineVerticesEdgesTests

using Test
using STLCutters


using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Helpers
using Gridap.Arrays

using STLCutters: get_default_directions
using STLCutters: compute_stl_topology 
using STLCutters: compute_stl_grid 
using STLCutters: get_edge_coordinates 


p = QUAD
T,X = initial_mesh(p)


STL_v = [ 
  Point(-0.1,0.2),
  Point(0.5,0.3),
  Point(0.8,0.5),
  Point(1.1,0.5) ]
STL_faces = Table( [[1,2],[2,3],[3,4]] )

stl_grid = compute_stl_grid(STL_faces,STL_v)
stl = GridTopology(stl_grid)

STL_vertices = get_vertex_coordinates(stl)
STL_edges = get_edge_coordinates(stl)



T,X = initial_mesh(p)

V = distribute_faces(T,X,p,1:length(STL_vertices),STL_vertices)

Tnew = eltype(T)[]

Tnew_to_v = Vector{Int}[]

v_in = Int[]

insert_vertices!(T,X,p,V,Tnew,STL_vertices,Tnew_to_v,v_in)

T = Tnew
T_to_v = Tnew_to_v


E = distribute_faces(T,X,p,1:length(STL_edges),STL_edges)

vs = get_default_directions(1:length(STL_edges),STL_edges)

Tnew = eltype(T)[]

Tnew_to_e = Vector{Int}[]

e_in = Int[]

insert_edges!(T,X,p,E,Tnew,STL_edges,Tnew_to_e,e_in,vs)

T = Tnew
T_to_e = Tnew_to_e

@test length(T) == 6

grid = compute_grid(T,X,QUAD)

writevtk(grid,"Tree")
writevtk(stl_grid,"STL")


p = HEX
T,X = initial_mesh(p)


STL_vertices = [ 
  Point(0.4,0.4,0.5),
  Point(1.1,0.3,0.4),
  Point(1.1,0.9,0.7),
  Point(-0.1,0.5,0.3),
  Point(0.5,-0.1,0.5) ]

STL_T = Table( [[1,2,3],[1,2,5],[1,4,5],[1,4,3]] )

stl_grid = compute_stl_grid(STL_T,STL_vertices)
stl = GridTopology(stl_grid)

STL_vertices = get_vertex_coordinates(stl)
STL_edges = get_edge_coordinates(stl)

# 3D

T,X = initial_mesh(p)

V = distribute_faces(T,X,p,1:length(STL_vertices),STL_vertices)

Tnew = eltype(T)[]

Tnew_to_v = Vector{Int}[]

v_in = Int[]

insert_vertices!(T,X,p,V,Tnew,STL_vertices,Tnew_to_v,v_in)

T = Tnew
T_to_v = Tnew_to_v


E = distribute_faces(T,X,p,1:length(STL_edges),STL_edges)

vs = get_default_directions(1:length(STL_edges),STL_edges)

Tnew = eltype(T)[]

Tnew_to_e = Vector{Int}[]

e_in = Int[]

insert_edges!(T,X,p,E,Tnew,STL_edges,Tnew_to_e,e_in,vs)

T = Tnew
T_to_e = Tnew_to_e

grid = compute_grid(T,X,HEX)

writevtk(grid,"Tree3D")
writevtk(stl_grid,"STL3D")

end # module
