module GridapIntegration

using STLCutter.Cutter
const _Cutter = Cutter

import Gridap

using GridapEmbedded.CSG
using GridapEmbedded.Interfaces


using STLCutter.Cutter: BulkMesh
import STLCutter.Cutter: CartesianMesh

using STLCutter.Cutter: is_cell_interior
using STLCutter.Cutter: is_cell_exterior
using STLCutter.Cutter: is_cell_cut

using Gridap.Geometry: get_cartesian_descriptor

import GridapEmbedded.Interfaces: cut
import GridapEmbedded.CSG: get_tree
import GridapEmbedded.CSG: compatible_geometries
import GridapEmbedded.CSG: similar_geometry

export STLGeometry
export SurfaceMeshGeometry
export SurfaceMeshCutter

export cut

include("STLGeometries.jl")
include("SurfaceMeshGeometries.jl")
include("SurfaceMeshCutters.jl")

end # module
