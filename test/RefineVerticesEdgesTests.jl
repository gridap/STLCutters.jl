module RefineVerticesEdgesTests

using Test
using STLCutters


using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Helpers
using Gridap.Arrays

using STLCutters: get_default_directions
using STLCutters: compute_stl_model 
using STLCutters: distribute_vertices 
using STLCutters: distribute_edges 

## 2D

p = QUAD
T,X = initial_mesh(p)
stl_vertices = [ 
  Point(-0.1,0.2),
  Point(0.5,0.3),
  Point(0.8,0.5),
  Point(1.1,0.5) ]
stl_faces = Table( [[1,2],[2,3],[3,4]] )
stl = compute_stl_model(stl_faces,stl_vertices)
# Vertices
Tnew = eltype(T)[]
fnew = Vector{Int}[]
f = fill(Int[],length(T))
V = distribute_vertices(T,X,p,stl,1:num_vertices(stl))
insert_vertices!(T,X,p,stl,V,f,Tnew,fnew)
T = Tnew
f = fnew
# Edges
Tnew = eltype(T)[]
fnew = Vector{Int}[]
E = distribute_edges(T,X,p,stl,1:num_edges(stl))
vs = get_default_directions(1:num_edges(stl),stl)
insert_edges!(T,X,p,stl,E,f,Tnew,fnew,vs)
T = Tnew
f = fnew
@test length(T) == 6
grid = compute_grid(T,X,p)
writevtk(grid,"Tree")
writevtk(get_grid(stl),"stl")

## 3D 

p = HEX
T,X = initial_mesh(p)
stl_vertices = [ 
  Point(0.4,0.4,0.5),
  Point(1.1,0.3,0.4),
  Point(1.1,0.9,0.7),
  Point(-0.1,0.5,0.3),
  Point(0.5,-0.1,0.5) ]
stl_faces = Table( [[1,2,3],[1,2,5],[1,4,5],[1,4,3]] )
stl = compute_stl_model(stl_faces,stl_vertices)
# Vertices
Tnew = eltype(T)[]
fnew = Vector{Int}[]
f = fill(Int[],length(T))
V = distribute_vertices(T,X,p,stl,1:num_vertices(stl))
insert_vertices!(T,X,p,stl,V,f,Tnew,fnew)
T = Tnew
f = fnew
# Edges
Tnew = eltype(T)[]
fnew = Vector{Int}[]
E = distribute_edges(T,X,p,stl,1:num_edges(stl))
vs = get_default_directions(1:num_edges(stl),stl)
insert_edges!(T,X,p,stl,E,f,Tnew,fnew,vs)
T = Tnew
f = fnew
grid = compute_grid(T,X,p)
writevtk(grid,"Tree3D")
writevtk(get_grid(stl),"stl3D")

end # module
