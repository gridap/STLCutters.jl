module CellMeshesTests

using Test
using Gridap
using Gridap.Geometry

using STLCutters
using STLCutters.CellMeshes
using STLCutters.CellMeshes: BoundingBox, CellMesh, SimplexFacet
using STLCutters.CellMeshes: compute_cell_mesh!, reset! 
using STLCutters.CellMeshes: get_cell_to_vertices 
using STLCutters.CellMeshes: get_facet_to_vertices 
using STLCutters.CellMeshes: get_vertex_coordinates 
using STLCutters.CellMeshes: get_cell_to_inout 
using STLCutters.CellMeshes: get_facet_to_inout 

f1 = [ Point(0.0,0.5), Point(1.0,1.1) ]
f2 = [ Point(0.6,0.0), Point(0.0,0.3) ] 
levelsets = [ SimplexFacet(f1...), SimplexFacet(f2...) ]

box = BoundingBox(Point(0.0,0.0),Point(1.0,1.0))

mesh = CellMesh(box)

compute_cell_mesh!(mesh,levelsets)

X = get_vertex_coordinates(mesh)

T = get_cell_to_vertices(mesh)
F = get_facet_to_vertices(mesh)

cell_to_io = get_cell_to_inout(mesh)
facet_to_io = get_facet_to_inout(mesh)

grid = UnstructuredGrid(X,T,[TRI3],ones(length(T)))
facets = UnstructuredGrid(X,F,[SEG2],ones(length(F)))

#writevtk(grid,"mesh1",cellfields=["io"=>cell_to_io])
#writevtk(facets,"facets1",cellfields=["io"=>facet_to_io])


X = [
  Point(0.0,0.2),
  Point(1.0,0.0),
  Point(0.3,0.7),
  Point(0.7,1.0) ]
K = [1,2,3,4]

p = QUAD

mesh = CellMesh(X,K,p)
compute_cell_mesh!(mesh,levelsets)

X = get_vertex_coordinates(mesh)

T = get_cell_to_vertices(mesh)
F = get_facet_to_vertices(mesh)

cell_to_io = get_cell_to_inout(mesh)
facet_to_io = get_facet_to_inout(mesh)

grid = UnstructuredGrid(X,T,[TRI3],ones(length(T)))
facets = UnstructuredGrid(X,F,[SEG2],ones(length(F)))

#writevtk(grid,"mesh2",cellfields=["io"=>cell_to_io])
#writevtk(facets,"facets2",cellfields=["io"=>facet_to_io])


v = [ Point(1.1,0.6),Point(0.5,1.0),Point(0.0,0.5),Point(1.1,0.4) ]
levelsets = [ SimplexFacet(v[i],v[i+1]) for i in 1:length(v)-1 ]

X = [
  Point(0.0,0.0),
  Point(1.0,0.0),
  Point(0.0,1.0),
  Point(1.0,1.0) ]
K = [1,2,3,4]

p = QUAD

mesh = CellMesh(X,K,p)
compute_cell_mesh!(mesh,levelsets)

X = get_vertex_coordinates(mesh)

T = get_cell_to_vertices(mesh)
F = get_facet_to_vertices(mesh)

cell_to_io = get_cell_to_inout(mesh)
facet_to_io = get_facet_to_inout(mesh)

grid = UnstructuredGrid(X,T,[TRI3],ones(length(T)))
facets = UnstructuredGrid(X,F,[SEG2],ones(length(F)))

#writevtk(grid,"mesh3",cellfields=["io"=>cell_to_io])
#writevtk(facets,"facets3",cellfields=["io"=>facet_to_io])

X = Point{3,Float64}[
  (0.0, 0.0, 10.0),                         
  (0.8, 0.8, 0.0), 
  (0.5, 0.5, 0.0), 
  (0.8, 0.8, 0.6), 
  (0.5, 0.5, 0.6), 
  (0.7, 0.9, 0.0), 
  (0.0, 0.6, 0.0), 
  (0.7, 0.9, 0.6), 
  (0.0, 0.6, 0.6) ]
K = 2:length(X)

p = HEX

mesh = CellMesh(X,K,p)

p1 = Point(-0.1,0.5,0.5)
p2 = Point(0.5,1.1,0.5)
p3 = Point(1.1,0.5,0.5)

levelsets = [ SimplexFacet(p1,p2,p3) ]


mesh = CellMesh(X,K,p)
compute_cell_mesh!(mesh,levelsets)

X = get_vertex_coordinates(mesh)

T = get_cell_to_vertices(mesh)
F = get_facet_to_vertices(mesh)

cell_to_io = get_cell_to_inout(mesh)
facet_to_io = get_facet_to_inout(mesh)

grid = UnstructuredGrid(X,T,[TET4],ones(length(T)))
facets = UnstructuredGrid(X,F,[TRI3],ones(length(F)))

#writevtk(grid,"mesh",cellfields=["io"=>cell_to_io])
#writevtk(facets,"facets",cellfields=["io"=>facet_to_io])

end # module
