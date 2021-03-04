module GridapIntegration

using STLCutters

import Gridap

using GridapEmbedded.CSG
using GridapEmbedded.Interfaces
import GridapEmbedded.LevelSetCutters

using STLCutters: BulkMesh
import STLCutters: CartesianMesh

using STLCutters: is_cell_interior
using STLCutters: is_cell_exterior
using STLCutters: is_cell_cut

import STLCutters: surface 

using Gridap.Geometry: get_cartesian_descriptor

import GridapEmbedded.Interfaces: cut
import GridapEmbedded.CSG: get_tree
import GridapEmbedded.CSG: compatible_geometries
import GridapEmbedded.CSG: similar_geometry

export STLGeometry

export cut
export surface
export square

include("STLGeometries.jl")
include("STLCutters.jl")

end # module