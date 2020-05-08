module GridapIntegration

using STLCutter.Cutter
import GridapEmbedded



include("STLGeometries.jl")
include("SurfaceMeshGeometries.jl")
include("SurfaceMeshCutters.jl")


#TODO: 
#
# Get the same functionallities with Gridap.CartesianDiscreteModel than in VolumeMesh or CartesianMesh, addapt them to if
# num_faces, get_faces seem similar
# tables do not have acces as t[i,j]

# Convert Gridap.CartesianDiscreteModel into CartesianMesh, through get_cartesiant_descriptor
# be sure that the resulting meshes are equivalent in terms of cells

end # module
