module RefineFacetsTests

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
using STLCutters: distribute_facets 
using STLCutters: define_cells!
using STLCutters: volume 
using STLCutters: volumes 
using STLCutters: FACE_IN

function test_facets(p,vertices,facets)
  stl = compute_stl_model(Table(facets),vertices)
  T,X = initial_mesh(p)
  f = fill(Int[],length(T))
  Tnew,fnew = empty(T), empty(f)
  V = distribute_vertices(T,X,p,stl,1:num_vertices(stl))
  insert_vertices!(T,X,p,stl,V,f,Tnew,fnew)
  T,f = Tnew,fnew
  vertex_grid = compute_grid(T,X,p)
  if num_dims(p) == 3
    Tnew,fnew = empty(T), empty(f)
    E = distribute_edges(T,X,p,stl,1:num_edges(stl))
    vs = get_default_directions(vcat(E...),stl)
    insert_edges!(T,X,p,stl,E,f,Tnew,fnew,vs)
    T,f = Tnew,fnew
    edge_grid = compute_grid(T,X,p)
  end
  Tnew,fnew = empty(T), empty(f)
  cell_types,cell_to_io = Int8[], Int8[]
  F = distribute_facets(T,X,p,stl,1:num_cells(stl))
  insert_facets!(T,X,p,stl,F,f,Tnew,fnew,cell_types,cell_to_io)
  T,f = Tnew,fnew
  if num_dims(p) == 2
    reffes = [QUAD4,TRI3]
  else
    reffes = [HEX8,TET4]
  end
  grid = UnstructuredGrid(X,Table(T),reffes,cell_types)
  define_cells!(grid,f,cell_to_io)
  if num_dims(p) == 2
    grids = vertex_grid,grid
  else
    grids = vertex_grid,edge_grid,grid
  end
  grids...,cell_to_io,stl
end

# 2D

p = QUAD

vertices = [
  Point(-0.5,1.0),
  Point(0.4,0.3),
  Point(1.5,1.0) ]
facets = [[1,2],[2,3]]
vgrid,grid,cell_to_io,stl = test_facets(p,vertices,facets)
@test count(isequal(UNSET),cell_to_io) == 0
@test volume(grid) ≈ 1
@test sum(volumes(grid) .* (cell_to_io .== FACE_IN)) ≈ 0.5232323232323232
#writevtk(grid,"mesh2D",cellfields=["IO"=>cell_to_io])
#writevtk(get_grid(stl),"stl2D")

vertices = [
  Point(-0.5,1.0),
  Point(0.4,0.3),
  Point(1.5,0.3) ]
facets = [[1,2],[2,3]]
vgrid,grid,cell_to_io,stl = test_facets(p,vertices,facets)
@test count(isequal(UNSET),cell_to_io) == 0
@test volume(grid) ≈ 1
@test sum(volumes(grid) .* (cell_to_io .== FACE_IN)) ≈ 0.637777777777777
#writevtk(grid,"mesh2D",cellfields=["IO"=>cell_to_io])
#writevtk(get_grid(stl),"stl2D")

vertices = [
  Point(-0.5,1.0),
  Point(0.4,0.3),
  Point(1.5,1.0) ]
facets = [[1,2],[2,3]]
vgrid,grid,cell_to_io,stl = test_facets(p,vertices,facets)
@test count(isequal(UNSET),cell_to_io) == 0
@test volume(grid) ≈ 1
@test sum(volumes(grid) .* (cell_to_io .== FACE_IN)) ≈ 0.5232323232323232
#writevtk(grid,"mesh2D",cellfields=["IO"=>cell_to_io])
#writevtk(get_grid(stl),"stl2D")

vertices = [
  Point(1.1,0.3),
  Point(0.6,0.3),
  Point(0.5,0.2),
  Point(0.2,0.1),
  Point(0.0,-0.1) ]
facets = [[2,1],[3,2],[4,3],[5,4]]
vgrid,grid,cell_to_io,stl = test_facets(p,vertices,facets)
@test count(isequal(UNSET),cell_to_io) == 0
@test volume(grid) ≈ 1
@test sum(volumes(grid) .* (cell_to_io .== FACE_IN)) ≈ 0.805
#writevtk(grid,"mesh2D",cellfields=["IO"=>cell_to_io])
#writevtk(get_grid(stl),"stl2D")

vertices = [
  Point(0.0,-0.1),
  Point(0.7,0.2),
  Point(0.9,0.6),
  Point(0.7,0.8),
  Point(1.0,1.1) ]
facets = [[1,2],[2,3],[3,4],[4,5]]
vgrid,grid,cell_to_io,stl = test_facets(p,vertices,facets)
@test count(isequal(UNSET),cell_to_io) == 0
@test volume(grid) ≈ 1
@test sum(volumes(grid) .* (cell_to_io .== FACE_IN)) ≈ 0.733333333333
writevtk(grid,"mesh2D",cellfields=["IO"=>cell_to_io])
writevtk(get_grid(stl),"stl2D")

# 3D

p = HEX

vertices = [
  Point(-2.0,-2.0,3.0),
  Point(0.4,-2.0,0.3),
  Point(0.4,2.0,0.3),
  Point(3.0,3.0,3.0) ]
facets = [[1,2,3],[2,4,3]]
vgrid,egrid,grid,cell_to_io,stl = test_facets(p,vertices,facets)
@test count(isequal(UNSET),cell_to_io) == 0
@test volume(grid) ≈ 1
#writevtk(egrid,"edge_mesh")
#writevtk(grid,"mesh",cellfields=["IO"=>cell_to_io])
#writevtk(get_grid(stl),"stl")

vertices = [
  Point(-0.1,0.5,0.5),
  Point(0.5,1.1,0.5),
  Point(1.1,0.5,0.5),
  Point(0.5,0.1,0.6),
  Point(1.1,1.1,0.6),
  Point(-0.1,1.1,0.5),
  Point(-0.1,-0.1,0.3),
  Point(0.5,-0.1,0.2),
  Point(1.1,-0.1,0.3) ]
facets = [[1,2,3],[1,3,4],[1,6,2],[3,2,5],[1,4,7],[4,8,7],[4,9,8],[3,9,4]]
vgrid,egrid,grid,cell_to_io,stl = test_facets(p,vertices,facets)
@test count(isequal(UNSET),cell_to_io) == 0
@test volume(grid) ≈ 1
#writevtk(egrid,"edge_mesh")
#writevtk(grid,"mesh",cellfields=["IO"=>cell_to_io])
#writevtk(get_grid(stl),"stl")

end # module
