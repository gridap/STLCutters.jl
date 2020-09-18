module RefineEdgesTests

using Test
using STLCutters


using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Helpers
using Gridap.Arrays

using STLCutters: get_default_directions
using STLCutters: compute_stl_topology 
using STLCutters: get_edge_coordinates 


p = QUAD
T,X = initial_mesh(p)

stl_v = [ 
  Point(-0.1,0.5),
  Point(0.5,1.1),
  Point(1.1,0.5) ]
stl_faces = Table( [[1,2],[1,3] ] )

stl = compute_stl_topology(stl_faces,stl_v)

STL_edges = get_edge_coordinates(stl)

E = distribute_faces(T,X,p,1:length(STL_edges),STL_edges)

vs = get_default_directions(E[1],STL_edges)

Tnew = eltype(T)[]

Tnew_to_e = Vector{Int}[]

e_in = Int[]

insert_edges!(T,X,p,E,Tnew,STL_edges,Tnew_to_e,e_in,vs)

T = Tnew
T_to_e = Tnew_to_e

@test length(T) == 6

grid = compute_grid(T,X,QUAD)

writevtk(grid,"Tree")


p = QUAD
T,X = initial_mesh(p)

stl_v = [ 
  Point(0.0,0.0),
  Point(0.5,1.1),
  Point(1.1,1.1) ]
stl_faces = Table( [[1,2],[1,3] ] )

stl = compute_stl_topology(stl_faces,stl_v)
STL_edges = get_edge_coordinates(stl)

E = distribute_faces(T,X,p,1:length(STL_edges),STL_edges)

vs = get_default_directions(E[1],STL_edges)

Tnew = eltype(T)[]

Tnew_to_e = Vector{Int}[]

e_in = Int[]

insert_edges!(T,X,p,E,Tnew,STL_edges,Tnew_to_e,e_in,vs)

T = Tnew
T_to_e = Tnew_to_e

grid = compute_grid(T,X,QUAD)

@test length(T) == 3

writevtk(grid,"Tree")


p = HEX
T,X = initial_mesh(p)

p1 = Point(-0.1,0.5,0.6)
p2 = Point(0.5,1.1,0.8)
edge1 = Segment(p1,p2) 

p3 = Point(-0.1,0.5,0.1)
p4 = Point(1.1,0.5,0.6)

edge2 = Segment(p3,p4) 

STL_edges = [edge1,edge2]

E = distribute_faces(T,X,p,1:length(STL_edges),STL_edges)

vs = get_default_directions(E[1],STL_edges)

Tnew = eltype(T)[]

Tnew_to_e = Vector{Int}[]

e_in = Int[]

insert_edges!(T,X,p,E,Tnew,STL_edges,Tnew_to_e,e_in,vs)

T = Tnew
T_to_e = Tnew_to_e

grid = compute_grid(T,X,HEX)


writevtk(grid,"Tree3D")


X = [p1,p2,p3,p4]
T = [[1,2],[3,4]]

stl = compute_grid(T,X,SEGMENT)


writevtk(stl,"stl")

p = HEX
T,X = initial_mesh(p)

p1 = Point(-0.1,0.5,0.5)
p2 = Point(1.1,0.55,0.55)
edge1 = Segment(p1,p2) 

p3 = Point(0.4,-0.1,0.4)
p4 = Point(0.45,1.1,0.45)

edge2 = Segment(p3,p4) 

p5 = Point(0.3,0.3,-0.1)
p6 = Point(0.35,0.35,1.1)

edge3 = Segment(p5,p6) 

STL_edges = [edge3,edge2,edge1]

E = distribute_faces(T,X,p,1:length(STL_edges),STL_edges)

vs = get_default_directions(E[1],STL_edges)

Tnew = eltype(T)[]

Tnew_to_e = Vector{Int}[]

e_in = Int[]

insert_edges!(T,X,p,E,Tnew,STL_edges,Tnew_to_e,e_in,vs)

T = Tnew
T_to_e = Tnew_to_e

grid = compute_grid(T,X,HEX)

@test length(T) == 6

writevtk(grid,"Tree3D")

X = [p1,p2,p3,p4,p5,p6]
T = [[1,2],[3,4],[5,6]]

stl = compute_grid(T,X,SEGMENT)


writevtk(stl,"stl")

end # module
