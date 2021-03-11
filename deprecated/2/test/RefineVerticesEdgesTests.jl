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
using STLCutters: compute_model 
using STLCutters: distribute_vertices 
using STLCutters: distribute_edges 
using STLCutters: volume 

function test_vertices_edges(p,vertices,faces,pf)
  stl = compute_model(Table(faces),vertices,pf)
  T,X = initial_mesh(p)
  Tnew,fnew = eltype(T)[], Vector{Int}[]
  f = fill(Int[],length(T))
  V = distribute_vertices(T,X,p,stl,1:num_vertices(stl))
  insert_vertices!(T,X,p,stl,V,f,Tnew,fnew)
  T,f = Tnew,fnew
  Tnew,fnew = eltype(T)[], Vector{Int}[]
  E = distribute_edges(T,X,p,stl,1:num_edges(stl))
  vs = get_default_directions(1:num_edges(stl),stl)
  insert_edges!(T,X,p,stl,E,f,Tnew,fnew,vs)
  Tnew,X,fnew,stl
end

## 2D

p = QUAD

vertices = [
  Point(-0.1,0.2),
  Point(0.5,0.3),
  Point(0.8,0.5),
  Point(1.1,0.5) ]
edges = [[1,2],[2,3],[3,4]]
T,X,f,stl = test_vertices_edges(p,vertices,edges,SEGMENT)
grid = compute_grid(T,X,p)
@test length(T) == 6
@test volume(grid) ≈ 1
#writevtk(grid,"Tree")
#writevtk(get_grid(stl),"stl")

## 3D 

p = HEX
vertices = [
  Point(0.4,0.4,0.5),
  Point(1.1,0.3,0.4),
  Point(1.1,0.9,0.7),
  Point(-0.1,0.5,0.3),
  Point(0.5,-0.1,0.5) ]
faces = [[1,2,3],[1,2,5],[1,4,5],[1,4,3]]
T,X,f,stl = test_vertices_edges(p,vertices,faces,TRI)
grid = compute_grid(T,X,p)
@test length(T) == 11
@test volume(grid) ≈ 1
#writevtk(grid,"Tree")
#writevtk(get_grid(stl),"stl")

end # module
