module RefineEdgesTests

using Test
using STLCutters


using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Helpers
using Gridap.Arrays

using STLCutters: get_default_directions


p = QUAD
T,X = initial_mesh(p)

p1 = Point(-0.1,0.5)
p2 = Point(0.5,1.1)
edge1 = Segment(p1,p2) 

p1 = Point(-0.1,0.5)
p2 = Point(1.1,0.5)
edge2 = Segment(p1,p2) 

STL_edges = [edge1,edge2]

E = distribute_faces(T,X,p,1:length(STL_edges),STL_edges)

vs = get_default_directions(E[1],STL_edges)

Tnew = eltype(T)[]

Tnew_to_e = Vector{Int}[]

e_in = Int[]

insert_edges!(T,X,p,E,Tnew,STL_edges,Tnew_to_e,e_in,vs)

T = Tnew
T_to_e = Tnew_to_e

grid = compute_grid(T,X,QUAD)

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

end # module
