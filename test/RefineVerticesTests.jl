module RefineVerticesTest

using Test
using STLCutters


using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry

using STLCutters: compute_stl_model
using STLCutters: compute_model
using STLCutters: distribute_vertices 
using STLCutters: volume 

function test_vertices(p,vertices)
  T,X = initial_mesh(p)
  stl_vertices = vertices
  stl_faces = Table( [ [i] for i in 1:length(vertices) ] )
  stl = compute_model(stl_faces,stl_vertices,VERTEX)
  Tnew,fnew = eltype(T)[], Vector{Int}[]
  f = fill(Int[],length(T))
  V = distribute_vertices(T,X,p,stl,1:num_vertices(stl))
  insert_vertices!(T,X,p,stl,V,f,Tnew,fnew)
  Tnew,X,fnew
end

# 2D

p = QUAD

vertices = [ Point(0.1,0.2) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
@test T == [[1,5,3,6],[5,2,6,4]]
@test volume(grid) ≈ 1

vertices = [ Point(0.2,0.1) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
@test T == [[1,2,5,6],[5,6,3,4]]
@test volume(grid) ≈ 1

vertices = [ Point(0.0,0.2) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1

vertices = [ Point(1e-10,0.2) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1
@test vertices[1] == Point(0.0,0.2)

vertices = [ Point(-1e-10,0.2) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1
@test vertices[1] == Point(0.0,0.2)

vertices = [ Point(0.2,1e-10) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1
@test vertices[1] == Point(0.2,0.0)

vertices = [ 
  Point(0.1,0.2),
  Point(0.5,0.5),
  Point(0.4,0.1),
  Point(0.3,0.3) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
D = num_dims(p)
@test length(T) == length(f) == length(vertices)+1
@test length(X) == (2^(D-1))*length(vertices)+2^D
@test volume(grid) ≈ 1

vertices = [ 
  Point(0.1,0.2),
  Point(0.5,0.5),
  Point(0.4,0.1),
  Point(0.5-1e-10,0.3) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
D = num_dims(p)
@test length(T) == length(f) == length(vertices)
@test length(X) == (2^(D-1))*(length(vertices)-1)+2^D
@test volume(grid) ≈ 1

## 3D

p = HEX

vertices = [ Point(0.3,0.2,0.1) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
@test T == [[1,2,3,4,9,10,11,12],[9,10,11,12,5,6,7,8]]
@test volume(grid) ≈ 1

vertices = [ Point(0.2,0.1,0.3) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
@test T == [[1,2,9,10,5,6,11,12],[9,10,3,4,11,12,7,8]]
@test volume(grid) ≈ 1

vertices = [ Point(0.1,0.3,0.2) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
@test T == [[1,9,3,10,5,11,7,12],[9,2,10,4,11,6,12,8]]
@test volume(grid) ≈ 1

vertices = [ Point(0.0,0.2,0.2) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1

vertices = [ Point(1e-10,0.2,0.2) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1
@test vertices[1] == Point(0.0,0.2,0.2)

vertices = [ Point(1e-10,1e-10,0.2) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1
@test vertices[1] == Point(0.0,0.0,0.2)

vertices = [ Point(1e-10,1e-10,1e-10) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
@test length(T) == 1
@test volume(grid) ≈ 1
@test vertices[1] == Point(0.0,0.0,0.0)

vertices = [ 
  Point(0.1,0.2,0.3),
  Point(0.5,0.5,0.5),
  Point(0.4,0.1,0.2),
  Point(0.3,0.7,0.4) ]
T,X,f = test_vertices(p,vertices)
grid = compute_grid(T,X,p)
D = num_dims(p)
@test length(T) == length(f) == length(vertices)+1
@test length(X) == (2^(D-1))*length(vertices)+2^D
@test volume(grid) ≈ 1

end # module
