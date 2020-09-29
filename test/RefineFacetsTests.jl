module RefineFacetsTests

using Test
using STLCutters


using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Helpers
using Gridap.Arrays


using STLCutters: get_default_directions
using STLCutters: compute_stl_grid 
using STLCutters: get_edge_coordinates 
using STLCutters: get_facet_coordinates 
using STLCutters: is_on_cell_facet 

stl_v = [ 
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

stl_grid = compute_stl_grid(stl_faces,stl_v)
stl = GridTopology(stl_grid)


STL_vertices = get_vertex_coordinates(stl)
STL_edges = get_edge_coordinates(stl)
STL_facets = get_facet_coordinates(stl)


p = HEX
T,X = initial_mesh(p)

V = distribute_faces(T,X,p,1:length(STL_vertices),STL_vertices)

Tnew = eltype(T)[]

Tnew_to_v = Vector{Int}[]

v_in = Int[]

insert_vertices!(T,X,p,V,Tnew,STL_vertices,Tnew_to_v,v_in)

T = Tnew
T_to_v = Tnew_to_v


E = distribute_faces(T,X,p,1:length(STL_edges),STL_edges)

vs = get_default_directions(1:length(STL_edges),STL_edges)

Tnew = eltype(T)[]

Tnew_to_e = Vector{Int}[]

e_in = Int[]

insert_edges!(T,X,p,E,Tnew,STL_edges,Tnew_to_e,e_in,vs)

T = Tnew
T_to_e = Tnew_to_e


grid = compute_grid(T,X,HEX)

writevtk(grid,"edge_mesh")

writevtk(stl_grid,"stl")

F = distribute_faces(T,X,p,1:length(STL_facets),STL_facets)

@test maximum( length.(F) ) == 1

Tnew = eltype(T)[]

cell_types = Int8[]
cell_to_io = Int8[]

insert_facets!(T,X,p,F,Tnew,cell_types,cell_to_io,STL_facets)

T = Tnew

grid = UnstructuredGrid(X,Table(T),[HEX8,TET4],cell_types)

writevtk(grid,"mesh",cellfields=["IO"=>cell_to_io])

end # module
