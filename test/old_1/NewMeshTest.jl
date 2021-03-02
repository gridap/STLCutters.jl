module _Gridap_CartesianMeshestTests

using Test

using Gridap
using Gridap.Geometry

import STLCutters
using STLCutter.GridapIntegration

using STLCutters: CartesianMesh
using STLCutters: VolumeMesh

using STLCutters: get_cell_coordinates
using STLCutters: get_vertices


pmin = Point(-1.0,-1.0,-1.0)
pmax = Point(2.0,2.2,2.5)

partition = (10,12,15)

bgmodel = CartesianDiscreteModel(pmin,pmax,partition)

cm = CartesianMesh(bgmodel)
vm = VolumeMesh(cm)

for cell in 1:num_cells(bgmodel)
  cell_points_mesh = convert(Vector{Point},collect(get_vertices(get_cell_coordinates(vm,cell))))

  cell_points_model = Geometry.get_cell_coordinates(get_triangulation(bgmodel))[cell]
  @test cell_points_mesh == cell_points_model
end


end # module
