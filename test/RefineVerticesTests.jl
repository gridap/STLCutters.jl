module RefineVerticesTest

using Test
using STLCutters


using Gridap
using Gridap.ReferenceFEs

## 2D

STL_vertices = [ 
  Point(0.1,0.2),
  Point(0.5,0.5),
  Point(0.4,0.1),
  Point(0.3,0.3) ]

p = QUAD

T,X = initial_mesh(p)

V = distribute_faces(T,X,p,1:length(STL_vertices),STL_vertices)

Tnew = eltype(T)[]

Tnew_to_v = Vector{Int}[]

v_in = Int[]

insert_vertices!(T,X,p,V,Tnew,STL_vertices,Tnew_to_v,v_in)

T = Tnew
T_to_v = Tnew_to_v

#display(T_to_v)

D = 2
@test length(T) == length(T_to_v) == length(STL_vertices)+1
@test T_to_v[1] == T_to_v[2] == [2,4,1,3]
@test T_to_v[3] == [2,4,1]
@test length(X) == (2^(D-1))*length(STL_vertices)+2^D

grid = compute_grid(T,X,p)

writevtk(grid,"Tree")

## 2D: Move vertices

STL_vertices = [ 
  Point(0.1,0.2),
  Point(0.5,0.5),
  Point(0.4,0.1),
  Point(0.5-1e-10,0.3) ]

p = QUAD

T,X = initial_mesh(p)

V = distribute_faces(T,X,p,1:length(STL_vertices),STL_vertices)

Tnew = eltype(T)[]

Tnew_to_v = Vector{Int}[]

v_in = Int[]

insert_vertices!(T,X,p,V,Tnew,STL_vertices,Tnew_to_v,v_in)

T = Tnew
T_to_v = Tnew_to_v

#display(T_to_v)

D = 2
@test length(T) == length(T_to_v) == length(STL_vertices)
@test T_to_v[1] == T_to_v[2] == [2,1,3]
@test T_to_v[3] == [2,1]
@test length(X) == (2^(D-1))*(length(STL_vertices)-1)+2^D

grid = compute_grid(T,X,p)

#writevtk(grid,"Tree")

## 3D

STL_vertices = [ 
  Point(0.1,0.2,0.3),
  Point(0.5,0.5,0.5),
  Point(0.4,0.1,0.2),
  Point(0.3,0.7,0.4) ]

p = HEX

T,X = initial_mesh(p)

V = distribute_faces(T,X,p,1:length(STL_vertices),STL_vertices)

Tnew = eltype(T)[]

Tnew_to_v = Vector{Int}[]

v_in = Int[]

insert_vertices!(T,X,p,V,Tnew,STL_vertices,Tnew_to_v,v_in)

T = Tnew
T_to_v = Tnew_to_v

#display(T_to_v)

D = 3
@test length(T) == length(T_to_v) == length(STL_vertices)+1
@test T_to_v[1] == T_to_v[2] == [2,4,1,3]
@test T_to_v[3] == [2,4,1]
@test length(X) == (2^(D-1))*length(STL_vertices)+2^D

grid = compute_grid(T,X,p)

writevtk(grid,"3DTree")

end # module
