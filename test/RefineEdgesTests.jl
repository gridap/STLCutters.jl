module RefineEdgesTests

using Test
using STLCutters


using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Helpers
using Gridap.Arrays

using STLCutters: get_default_directions
using STLCutters: compute_stl_model 
using STLCutters: compute_model 
using STLCutters: distribute_edges 


p = QUAD
T,X = initial_mesh(p)
stl_v = [ 
  Point(-0.1,0.5),
  Point(0.5,1.1),
  Point(1.1,0.5) ]
stl_faces = Table( [[1,2],[1,3] ] )
stl = compute_stl_model(stl_faces,stl_v)
Tnew = eltype(T)[]
fnew = Vector{Int}[]
f = fill(Int[],length(T))
E = distribute_edges(T,X,p,stl,1:num_edges(stl))
vs = get_default_directions(E[1],stl)
insert_edges!(T,X,p,stl,E,f,Tnew,fnew,vs)
T = Tnew
f = fnew
@test length(T) == 6
grid = compute_grid(T,X,p)
writevtk(grid,"Tree")

p = QUAD
T,X = initial_mesh(p)
stl_v = [ 
  Point(0.0,0.0),
  Point(0.5,1.1),
  Point(1.1,1.1) ]
stl_faces = Table( [[1,2],[1,3] ] )
stl = compute_stl_model(stl_faces,stl_v)
Tnew = eltype(T)[]
fnew = Vector{Int}[]
f = fill(Int[],length(T))
E = distribute_edges(T,X,p,stl,1:num_edges(stl))
vs = get_default_directions(E[1],stl)
insert_edges!(T,X,p,stl,E,f,Tnew,fnew,vs)
T = Tnew
f = fnew
grid = compute_grid(T,X,p)
@test length(T) == 3
writevtk(grid,"Tree")


p = HEX
T,X = initial_mesh(p)
stl_vertices = [
  Point(-0.1,0.5,0.6),
  Point(0.5,1.1,0.8),
  Point(-0.1,0.5,0.1),
  Point(1.1,0.5,0.6) ]
stl_faces = Table( [[1,2],[3,4]] )
stl = compute_model(stl_faces,stl_vertices,SEGMENT)
Tnew = eltype(T)[]
fnew = Vector{Int}[]
f = fill(Int[],length(T))
E = distribute_edges(T,X,p,stl,1:num_edges(stl))
vs = get_default_directions(E[1],stl)
insert_edges!(T,X,p,stl,E,f,Tnew,fnew,vs)
T = Tnew
f = fnew
grid = compute_grid(T,X,p)
writevtk(grid,"Tree3D")
writevtk(get_grid(stl),"stl")
@test length(T) == 3

p = HEX
T,X = initial_mesh(p)
stl_vertices = [
  Point(-0.1,0.5,0.5),
  Point(1.1,0.55,0.55),
  Point(0.4,-0.1,0.4),
  Point(0.45,1.1,0.45),
  Point(0.3,0.3,-0.1),
  Point(0.35,0.35,1.1) ]
stl_faces = Table( [[1,2],[3,4],[5,6]] )
stl = compute_model(stl_faces,stl_vertices,SEGMENT)
Tnew = eltype(T)[]
fnew = Vector{Int}[]
f = fill(Int[],length(T))
E = distribute_edges(T,X,p,stl,1:num_edges(stl))
vs = get_default_directions(E[1],stl)
insert_edges!(T,X,p,stl,E,f,Tnew,fnew,vs)
T = Tnew
f = fnew
grid = compute_grid(T,X,p)
writevtk(grid,"Tree3D")
writevtk(get_grid(stl),"stl")

end # module
