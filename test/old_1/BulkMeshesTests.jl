module BulkMeshesTests

using STLCutters

using STLCutters: expand, BulkMesh, surface, interior_volume, exterior_volume, move!, get_reference_cell, transformation

using Test


stl = STL(joinpath(@__DIR__,"data/cube.stl"))
sm = SurfaceMesh(stl)

bb = BoundingBox(stl)
bb = expand( bb, .5 )
cm = CartesianMesh(bb,5)

bulk = BulkMesh(cm,sm)

subcells = bulk.subtriangulation

c_to_bgc = subcells.cell_to_bgcell
c_to_v = subcells.cell_to_points
points = subcells.point_to_coords
rpoints = subcells.point_to_rcoords

ref_cell = get_reference_cell(cm)
for subcell in 1:num_cells(subcells)
  cell = get_cell(cm,c_to_bgc[subcell])
  for lv in 1:length(c_to_v,subcell)
    v = c_to_v[subcell,lv]
    p = points[v]
    rp = rpoints[v]
    @test have_intersection(p,cell)
    @test have_intersection(rp,ref_cell)
    @test transformation(ref_cell,cell,p) == rp
    @test transformation(cell,ref_cell,rp) ≈ p
  end
end

subfacets = bulk.facet_subtriangulations[1]

f_to_bgc = subfacets.facet_to_bgcell
f_to_v = subfacets.facet_to_points
points = subfacets.point_to_coords
rpoints = subfacets.point_to_rcoords

ref_cell = get_reference_cell(cm)
for subfacet in 1:num_facets(subfacets)
  cell = get_cell(cm,f_to_bgc[subfacet])
  for lv in 1:length(f_to_v,subfacet)
    v = f_to_v[subfacet,lv]
    p = points[v]
    rp = rpoints[v]
    @test have_intersection(p,cell)
    @test have_intersection(rp,ref_cell)
    @test transformation(ref_cell,cell,p) == rp
    @test transformation(cell,ref_cell,rp) ≈ p
  end
end


v = [
  Point(0.0,0.0),
  Point(0.0,1.0),
  Point(1.0,1.0),
  Point(1.0,0.0) ]

f2v = Table( [
  1 2 3 4;
  2 3 4 1 ] )

sm = SurfaceMesh(v,f2v)

box = BoundingBox(sm)
box = expand(box,0.1)

cm = CartesianMesh(box,2)
bulk = BulkMesh(cm,sm)

@test surface(bulk,1) == 4
@test interior_volume(bulk) == 1



geometries = [ "cube", "Bunny-LowPoly" ]
volumes = [ 1, 273280.0337419614 ]
meshes = [ 5, 20 ]

offset = VectorValue(1e-7,1e-7,1e-7)

for (i,geom) in enumerate( geometries )
#  println( "testing $geom.stl..." )
  stl = STL(joinpath(@__DIR__,"data/$geom.stl"))
  sm = SurfaceMesh(stl)
  for j in 1:1
    box = BoundingBox(sm)
    box = expand(box,0.1)
    bg_mesh = CartesianMesh(box, meshes[i] )
    bm = BulkMesh(bg_mesh,sm)

    writevtk(sm,"$(geom)_sm")
    writevtk(bm,"$(geom)_bm")


    @test interior_volume(bm) ≈ volumes[i]
    @test interior_volume(bm) + exterior_volume(bm) ≈ measure(box)
    @test surface(sm) ≈ surface(bm,1)
 
    sm = move!(sm,offset)
  end
end




## Selective tests
#
# geo = "wine_glass"
# n = 20
# tol = 1e-6 / sm tol -> tol*10
# step = 2
# box = 1.2*BoundingBox(stl)
# cell = 599

# 

end # module
