module HexahedraIntersectionTests

using STLCutters

using Gridap
using Gridap.ReferenceFEs
using Gridap.Arrays


using STLCutters: intersection
using STLCutters: compute_stl_model 
using STLCutters: compute_grid 
using STLCutters: get_cell 
using STLCutters: Cell, Plane 

p = HEX
X0 = collect( get_vertex_coordinates(p) )
K = collect( 1:length(X0) )

c = Cell( K,X0,p )



vertices = [
  Point(-0.5,-0.5,0.5),
  Point(-0.5,1.5,-0.5),
  Point(0.5,0.5,1.2),
  Point(1.5,0.5,0.5) ]
faces = [ [1,2,3],[2,1,4] ]
stl = compute_stl_model(Table(faces),vertices)

Π = [ Plane( get_cell(stl,i) ) for i in 1:num_cells(stl) ]

T,X,io,p = intersection(c,Π,UNSET,1e-9)

X = [X0;X]

grid = compute_grid(T,X,p)

writevtk(get_grid(stl),"intersection_stl")
writevtk(grid,"intersection_grid",cellfields=["IO"=>io])

end # module
