module VolumeMeshesTests

using STLCutters
using Test
using STLCutters: get_face_coordinates, get_vertex_coordinates, get_dface_to_vertices, get_cell_coordinates
using STLCutters: get_reference_cell

box = BoundingBox(
  Point( 0.2, 0.2, 0.2 ),
  Point( 0.7, 0.7, 0.7 ) )

cm = CartesianMesh( box, 10 )

vm = VolumeMesh(cm)

@test num_vertices(vm) == num_vertices(cm)
@test num_cells(vm) == num_cells(cm)

@test get_reference_cell(vm) == get_reference_cell(cm)

for cell in 1:num_cells(vm)
  cell_coordinates = get_cell_coordinates(vm,cell)
  @test get_cell(cm,cell) == BoundingBox( cell_coordinates )
end

end # module