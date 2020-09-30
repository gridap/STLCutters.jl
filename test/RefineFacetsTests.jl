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
using STLCutters: is_on_cell_facet 

p = HEX
T,X = initial_mesh(p)
stl_vertices = [ 
  Point(-0.1,0.5,0.5),
  Point(0.5,1.1,0.5),
  Point(1.1,0.5,0.5),
  Point(0.5,0.1,0.6),
  Point(1.1,1.1,0.6),
  Point(-0.1,1.1,0.5),
  Point(-0.1,-0.1,0.3),
  Point(0.5,-0.1,0.2),
  Point(1.1,-0.1,0.3) ]
stl_faces = [[1,2,3],[1,3,4],[1,6,2],[3,2,5],[1,4,7],[4,8,7],[4,9,8],[3,9,4]]
stl_faces = Table(stl_faces)
stl = compute_stl_model(stl_faces,stl_vertices)
# Vertices
Tnew = eltype(T)[]
Tnew_to_v = Vector{Int}[]
v_in = Int[]
V = distribute_vertices(T,X,p,stl,1:num_vertices(stl))
insert_vertices!(T,X,p,stl,V,Tnew,Tnew_to_v,v_in)
T = Tnew
T_to_v = Tnew_to_v
# Edges
Tnew = eltype(T)[]
Tnew_to_e = Vector{Int}[]
e_in = Int[]
E = distribute_edges(T,X,p,stl,1:num_edges(stl))
vs = get_default_directions(1:num_edges(stl),stl)
insert_edges!(T,X,p,stl,E,Tnew,Tnew_to_e,e_in,vs)
T = Tnew
T_to_e = Tnew_to_e
grid = compute_grid(T,X,p)
writevtk(grid,"edge_mesh")
# Facets
Tnew = eltype(T)[]
cell_types = Int8[]
cell_to_io = Int8[]
Fk = 1:num_cells(stl)
F = distribute_facets(T,X,p,stl,1:num_cells(stl))
@test maximum( length.(F) ) == 1
insert_facets!(T,X,p,stl,F,Tnew,cell_types,cell_to_io)
T = Tnew
grid = UnstructuredGrid(X,Table(T),[HEX8,TET4],cell_types)
writevtk(grid,"mesh",cellfields=["IO"=>cell_to_io])
writevtk(get_grid(stl),"stl")

end # module
