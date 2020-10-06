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
using STLCutters: volume

function test_edges(p,vertices,edges)
  stl = compute_model(Table(edges),vertices,SEGMENT)
  T,X = initial_mesh(p)
  Tnew,fnew = eltype(T)[], Vector{Int}[]
  f = fill(Int[],length(T))
  E = distribute_edges(T,X,p,stl,1:num_edges(stl))
  vs = get_default_directions(E[1],stl)
  insert_edges!(T,X,p,stl,E,f,Tnew,fnew,vs)
  Tnew,X,fnew,stl
end

p = QUAD

vertices = [
  Point(-0.1,0.6),
  Point(1.1,0.4) ]
edges = [[1,2],[1,2]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 2
@test volume(grid) ≈ 1

vertices = [
  Point(-0.1,0.5),
  Point(0.5,1.1) ]
edges = [[1,2]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 3
@test volume(grid) ≈ 1

vertices = [
  Point(0.0,-0.1),
  Point(0.0,1.1) ]
edges = [[1,2]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1

vertices = [
  Point(-0.1,0.0),
  Point(1.1,0.0) ]
edges = [[1,2]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1

vertices = [
  Point(-0.1,0.5),
  Point(0.5,1.1),
  Point(1.1,0.5) ]
edges = [[1,2],[1,3]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 6
@test volume(grid) ≈ 1

vertices = [
  Point(0.0,0.0),
  Point(0.5,1.1),
  Point(1.1,1.1) ]
edges = [[1,2],[1,3]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 3
@test volume(grid) ≈ 1

## 3D

p = HEX

vertices = [
  Point(-0.1,0.5,0.5),
  Point(1.1,0.5,0.4) ]
edges = [[1,2]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 2
@test volume(grid) ≈ 1

vertices = [
  Point(-0.1,0.0,0.5),
  Point(1.1,0.0,0.4) ]
edges = [[1,2]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1

vertices = [
  Point(-0.1,0.5,0.0),
  Point(1.1,0.4,0.0) ]
edges = [[1,2]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1

vertices = [
  Point(-0.1,1e-10,0.5),
  Point(1.1,0.0,0.4) ]
edges = [[1,2]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1

vertices = [
  Point(-0.1,0.0,0.0),
  Point(1.1,0.0,0.0) ]
edges = [[1,2]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1

vertices = [
  Point(-0.1,1e-10,1e-10),
  Point(1.1,0.0,0.0) ]
edges = [[1,2]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1

vertices = [
  Point(-0.1,0.5,0.6),
  Point(0.5,1.1,0.8),
  Point(-0.1,0.5,0.1),
  Point(1.1,0.5,0.6) ]
edges = [[1,2],[3,4]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 3
@test volume(grid) ≈ 1

vertices = [
  Point(-0.1,0.5,0.5),
  Point(1.1,0.55,0.55),
  Point(0.4,-0.1,0.4),
  Point(0.45,1.1,0.45),
  Point(0.3,0.3,-0.1),
  Point(0.35,0.35,1.1) ]
edges = [[1,2],[3,4],[5,6]]
T,X,f,stl = test_edges(p,vertices,edges)
grid = compute_grid(T,X,p)
@test length(T) == 8
@test volume(grid) ≈ 1
#writevtk(grid,"Tree")
#writevtk(get_grid(stl),"stl")

end # module
