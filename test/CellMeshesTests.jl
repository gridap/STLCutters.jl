module CellMeshesTests

using Test
using Gridap
using Gridap.Geometry

using STLCutters.CellMeshes
using STLCutters.CellMeshes: Plane,  BoundingBox, CellMesh
using STLCutters.CellMeshes: compute_cell_mesh!, reset! 
using STLCutters.CellMeshes: get_cell_to_vertices 
using STLCutters.CellMeshes: get_facet_to_vertices 
using STLCutters.CellMeshes: get_vertex_coordinates 
using STLCutters.CellMeshes: get_cell_to_inout 
using STLCutters.CellMeshes: get_facet_to_inout 

p1 = Point(0.0,0.5)
p2 = Point(0.6,0.0)
n1 = VectorValue(-1.0,1.0)
n2 = VectorValue(-0.5,-1.0)
n1 /= norm(n1)
n2 /= norm(n2)

levelsets = [ Plane(p1,n1),Plane(p2,n2) ]

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

#writevtk(grid,"mesh",cellfields=["io"=>cell_to_io])
#writevtk(facets,"facets",cellfields=["io"=>facet_to_io])


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

#writevtk(grid,"mesh",cellfields=["io"=>cell_to_io])
#writevtk(facets,"facets",cellfields=["io"=>facet_to_io])


end # module
