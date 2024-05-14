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
using STLCutters: check_graph
using STLCutters: compute_stl_model
using STLCutters: read_stl
using STLCutters: merge_nodes
using STLCutters: simplexify_cell_boundary
using STLCutters: compute_grid
using STLCutters: clipping


p = Polygon(TRI)
@test check_graph(p)
@test get_faces(p,0,0) == [[1],[2],[3]]
@test get_faces(p,1,0) == [[1,2],[2,3],[3,1]]
@test get_faces(p,2,0) == [[1,2,3]]
@test get_faces(p,0,1) == [[1,3],[1,2],[2,3]]
@test get_faces(p,1,1) == [[1],[2],[3]]
@test get_faces(p,2,1) == [[1,2,3]]
@test get_faces(p,0,2) == [[1],[1],[1]]
@test get_faces(p,1,2) == [[1],[1],[1]]
@test get_faces(p,2,2) == [[1]]
@test get_facedims(p) == [0,0,0,1,1,1,2]
@test Polytope{2}(p,1) === p
@test Polytope{1}(p,1) == SEGMENT
@test Polytope{0}(p,1) == VERTEX

p = Polygon(QUAD)
@test check_graph(p)
@test get_faces(p,0,0) == [[1],[2],[3],[4]]
@test get_faces(p,1,0) == [[1,2],[2,3],[3,4],[4,1]]
@test get_faces(p,2,0) == [[1,2,3,4]]
@test get_faces(p,0,1) == [[1,4],[1,2],[2,3],[3,4]]
@test get_faces(p,1,1) == [[1],[2],[3],[4]]
@test get_faces(p,2,1) == [[1,2,3,4]]
@test get_faces(p,0,2) == [[1],[1],[1],[1]]
@test get_faces(p,1,2) == [[1],[1],[1],[1]]
@test get_faces(p,2,2) == [[1]]
@test get_facedims(p) == [0,0,0,0,1,1,1,1,2]
@test Polytope{2}(p,1) === p
@test Polytope{1}(p,1) == SEGMENT
@test Polytope{0}(p,1) == VERTEX

p = Polyhedron(TET)
@test check_graph(p)
@test get_faces(p,0,0) == [[1],[2],[3],[4]]
@test get_faces(p,1,0) == [[1,2],[1,4],[1,3],[2,3],[2,4],[3,4]]
@test get_faces(p,2,0) == [[1,2,3],[1,4,2],[1,3,4],[2,4,3]]
@test get_faces(p,3,0) == [[1,2,3,4]]
@test get_faces(p,0,1) == [[1,2,3],[1,4,5],[3,4,6],[2,5,6]]
@test get_faces(p,1,1) == [[1],[2],[3],[4],[5],[6]]
@test get_faces(p,2,1) == [[1,4,3],[2,5,1],[3,6,2],[5,6,4]]
@test get_faces(p,3,1) == [[1,2,3,4,5,6]]
@test get_faces(p,0,2) == [[1,2,3],[1,2,4],[1,3,4],[2,3,4]]
@test get_faces(p,1,2) == [[1,2],[2,3],[1,3],[1,4],[2,4],[3,4]]
@test get_faces(p,2,2) == [[1],[2],[3],[4]]
@test get_faces(p,3,2) == [[1,2,3,4]]
@test get_faces(p,0,3) == [[1],[1],[1],[1]]
@test get_faces(p,1,3) == [[1],[1],[1],[1],[1],[1]]
@test get_faces(p,2,3) == [[1],[1],[1],[1]]
@test get_faces(p,3,3) == [[1]]
@test get_facedims(p) == [0,0,0,0,1,1,1,1,1,1,2,2,2,2,3]
@test Polytope{3}(p,1) === p
@test isa(Polytope{2}(p,1),Polygon)
@test Polytope{1}(p,1) == SEGMENT
@test Polytope{0}(p,1) == VERTEX


p = Polyhedron(HEX)
@test check_graph(p)
@test get_faces(p,0,0) == [[1],[2],[3],[4],[5],[6],[7],[8]]
@test get_faces(p,1,0) == [
  [1,5],[1,2],[1,3],[2,6],[2,4],[3,7],[3,4],[4,8],[5,7],[5,6],[6,8],[7,8]]
@test get_faces(p,2,0) == [
  [1,5,7,3],[1,2,6,5],[1,3,4,2],[2,4,8,6],[3,7,8,4],[5,6,8,7]]
@test get_faces(p,3,0) == [[1,2,3,4,5,6,7,8]]
@test get_faces(p,0,1) == [
  [1,2,3],[2,4,5],[3,6,7],[5,7,8],[1,9,10],[4,10,11],[6,9,12],[8,11,12]]
@test get_faces(p,1,1) == [[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12]]
@test get_faces(p,2,1) == [
  [1,9,6,3],[2,4,10,1],[3,7,5,2],[5,8,11,4],[6,12,8,7],[10,11,12,9]]
@test get_faces(p,3,1) == [[1,2,3,4,5,6,7,8,9,10,11,12]]
@test get_faces(p,0,2) == [
  [1,2,3],[2,3,4],[1,3,5],[3,4,5],[1,2,6],[2,4,6],[1,5,6],[4,5,6]]
@test get_faces(p,1,2) == [
  [1,2],[2,3],[1,3],[2,4],[3,4],[1,5],[3,5],[4,5],[1,6],[2,6],[4,6],[5,6]]
@test get_faces(p,2,2) == [[1],[2],[3],[4],[5],[6]]
@test get_faces(p,3,2) == [[1,2,3,4,5,6]]
@test get_faces(p,0,3) == [[1],[1],[1],[1],[1],[1],[1],[1]]
@test get_faces(p,1,3) == [[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]]
@test get_faces(p,2,3) == [[1],[1],[1],[1],[1],[1]]
@test get_faces(p,3,3) == [[1]]
@test get_facedims(p) == [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3]
@test Polytope{3}(p,1) === p
@test isa(Polytope{2}(p,1),Polygon)
@test Polytope{1}(p,1) == SEGMENT
@test Polytope{0}(p,1) == VERTEX

#

p = Polyhedron(TET,metadata=clipping)
@test check_graph(p)
@test volume(p) ≈ 1/6
@test surface(p) ≈ 3/2 + (√3)/2

p = Polyhedron(HEX,metadata=clipping)
p.metadata.vertex_to_planes

@test check_graph(p)
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

@test check_graph(p1)
@test surface(p1) ≈ 6

@test check_graph(p2)
@test surface(p2) ≈ 6


end # module
