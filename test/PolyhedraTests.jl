module PolyhedraTests

using Test

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using STLCutters


using STLCutters: Polyhedron
using STLCutters: Polygon
using STLCutters: STL
using STLCutters: restrict
using STLCutters: clip
using STLCutters: split
using STLCutters: surface, volume
using STLCutters: get_facet_planes
using STLCutters: compute_distances!
using STLCutters: compute_stl_model
using STLCutters: read_stl
using STLCutters: merge_nodes
using STLCutters: simplexify_cell_boundary
using STLCutters: compute_grid
using STLCutters: clipping
using STLCutters: check_polytope_graph


p = Polyhedron(TET,metadata=clipping)
@test check_polytope_graph(p)
@test volume(p) ≈ 1/6
@test surface(p) ≈ 3/2 + (√3)/2

p = Polyhedron(HEX,metadata=clipping)
p.metadata.vertex_to_planes

@test check_polytope_graph(p)
@test volume(p) ≈ 1
@test surface(p) ≈ 6

f, = writevtk(p,"p")
rm.(f)

vertices = [
  Point(0.1,-0.2,0.5),
  Point(1.0,-0.3,0.5),
  Point(-0.2,1.2,0.5)]

facet_to_vertices =
[[1,2,3]]

stlmodel = compute_stl_model( vertices, Table(facet_to_vertices) )
stl = STL(stlmodel)

Πf = get_facet_planes(stl)
Π = [get_offset(stl,2)+1]
compute_distances!(p,Πf,Π)

p⁻,p⁺ = split(p,Π[1])

pn = [p⁻,p⁺]

X,T,bgf = simplexify_cell_boundary(pn,HEX)

@test surface(compute_grid(X,T,TRI)) ≈ 6
@test volume(pn) ≈ 2*volume(p⁻) ≈ 1
@test surface(pn) ≈ 2*surface(p⁻) ≈ 8
@test surface(p⁻,stl) ≈ 1

X,T,N = read_stl(joinpath(@__DIR__,"data/cube.stl"))
stlmodel = compute_stl_model(X,T)
stlmodel = merge_nodes(stlmodel)
stltopo = get_grid_topology(stlmodel)
stl = STL(stlmodel)

p1 = Polyhedron(stl)
p2 = Polyhedron(stltopo)

@test check_polytope_graph(p1)
@test surface(p1) ≈ 6

@test check_polytope_graph(p2)
@test surface(p2) ≈ 6


end # module
