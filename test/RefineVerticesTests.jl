module RefineVerticesTest

using Test
using STLCutters


using Gridap
using Gridap.ReferenceFEs

using STLCutters: compute_stl_model
using STLCutters: compute_model
using STLCutters: distribute_vertices 

## 2D

p = QUAD
T,X = initial_mesh(p)
stl_vertices = [ 
  Point(0.1,0.2),
  Point(0.5,0.5),
  Point(0.4,0.1),
  Point(0.3,0.3) ]
stl_faces = Table( [[1],[2],[3],[4]] )
stl = compute_model(stl_faces,stl_vertices,VERTEX)
Tnew = eltype(T)[]
Tnew_to_v = Vector{Int}[]
v_in = Int[]
V = distribute_vertices(T,X,p,stl,1:num_vertices(stl))
insert_vertices!(T,X,p,stl,V,Tnew,Tnew_to_v,v_in)
T = Tnew
T_to_v = Tnew_to_v
D = 2
@test length(T) == length(T_to_v) == num_vertices(stl)+1
@test T_to_v[1] == T_to_v[2] == [2,4,1,3]
@test T_to_v[3] == [2,4,1]
@test length(X) == (2^(D-1))*num_vertices(stl)+2^D
grid = compute_grid(T,X,p)
writevtk(grid,"Tree")

## 2D: Move vertices

p = QUAD
T,X = initial_mesh(p)
stl_vertices = [ 
  Point(0.1,0.2),
  Point(0.5,0.5),
  Point(0.4,0.1),
  Point(0.5-1e-10,0.3) ]
stl_faces = Table( [[1],[2],[3],[4]] )
stl = compute_model(stl_faces,stl_vertices,VERTEX)
Tnew = eltype(T)[]
Tnew_to_v = Vector{Int}[]
v_in = Int[]
V = distribute_vertices(T,X,p,stl,1:num_vertices(stl))
insert_vertices!(T,X,p,stl,V,Tnew,Tnew_to_v,v_in)
T = Tnew
T_to_v = Tnew_to_v
D = 2
@test length(T) == length(T_to_v) == num_vertices(stl)
@test T_to_v[1] == T_to_v[2] == [2,1,3]
@test T_to_v[3] == [2,1]
@test length(X) == (2^(D-1))*(num_vertices(stl)-1)+2^D
grid = compute_grid(T,X,p)
#writevtk(grid,"Tree")

## 3D

p = HEX
T,X = initial_mesh(p)
stl_vertices = [ 
  Point(0.1,0.2,0.3),
  Point(0.5,0.5,0.5),
  Point(0.4,0.1,0.2),
  Point(0.3,0.7,0.4) ]
stl_faces = Table( [[1],[2],[3],[4]] )
stl = compute_model(stl_faces,stl_vertices,VERTEX)
Tnew = eltype(T)[]
Tnew_to_v = Vector{Int}[]
v_in = Int[]
V = distribute_vertices(T,X,p,stl,1:num_vertices(stl))
insert_vertices!(T,X,p,stl,V,Tnew,Tnew_to_v,v_in)
T = Tnew
T_to_v = Tnew_to_v
D = 3
@test length(T) == length(T_to_v) == num_vertices(stl)+1
@test T_to_v[1] == T_to_v[2] == [2,4,1,3]
@test T_to_v[3] == [2,4,1]
@test length(X) == (2^(D-1))*num_vertices(stl)+2^D
grid = compute_grid(T,X,p)
writevtk(grid,"3DTree")

end # module
