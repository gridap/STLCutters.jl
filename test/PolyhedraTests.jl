module PolyhedraTests

using Test

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using STLCutters

using STLCutters: Polyhedron
using STLCutters: STL
using STLCutters: restrict
using STLCutters: clip
using STLCutters: split
using STLCutters: surface, volume
using STLCutters: get_facet_planes
using STLCutters: compute_distances!
using STLCutters: plot
using STLCutters: check_graph
using STLCutters: compute_stl_model

p = Polyhedron(TRI)
@test check_graph(p)
#@test volume(p) ≈ 1/2
#@test surface(p) ≈ 2+√2

p = Polyhedron(QUAD)
@test check_graph(p)
#@test volume(p) ≈ 1
#@test surface(p) ≈ 4

p = Polyhedron(TET)
@test check_graph(p)
@test volume(p) ≈ 1/6
@test surface(p) ≈ 3/2 + (√3)/2

p = Polyhedron(HEX)
@test check_graph(p)
@test volume(p) ≈ 1
@test surface(p) ≈ 6

plot(p)

f, = writevtk(p,"p")
rm.(f)

vertices = [
  Point(0.1,-0.2,0.5),
  Point(1.0,-0.3,0.5),
  Point(-0.2,1.2,0.5)]

facet_to_vertices =
[[1,2,3]]

stlmodel = compute_stl_model( Table(facet_to_vertices), vertices )
stl = STL(stlmodel)

Πf = get_facet_planes(stl)
Π = [get_offset(stl,2)+1]
compute_distances!(p,Πf,Π)

p⁻,p⁺ = split(p,Π[1])

pn = [p⁻,p⁺]

@test volume(pn) ≈ 2*volume(p⁻) ≈ 1
@test surface(pn) ≈ 2*surface(p⁻) ≈ 8
@test surface(p⁻,stl) ≈ 1

end # module
