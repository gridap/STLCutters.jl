module issue_22_mwe


using Gridap
using GridapEmbedded
using STLCutters

geo = STLGeometry(joinpath(@__DIR__,"data/cube.stl"))

pmin = Point(-1.0, -1.0, 0.5)
pmax = Point(1.0, 1.0, 1.5)
n = 5
nz = 3
partition = (n,n,nz)

model = CartesianDiscreteModel(pmin, pmax, partition)
cutgeo = cut_facets(model, geo)
cutgeo_facets = cut_facets(cutgeo)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"top",[21])


Ωin = Triangulation(cutgeo,PHYSICAL_IN)
writevtk(Ωin,"trian_in")

Ωout = Triangulation(cutgeo, PHYSICAL_OUT)
writevtk(Ωout,"trian_out")

Γ = EmbeddedBoundary(cutgeo)
writevtk(Γ,"bnd")

Γbin = BoundaryTriangulation(cutgeo_facets,CUT_IN)
Γbout = BoundaryTriangulation(cutgeo_facets,CUT_OUT)
Λin =  SkeletonTriangulation(cutgeo_facets,CUT_IN)
Λout = SkeletonTriangulation(cutgeo_facets,CUT_OUT)

writevtk(Γbin,"bnd_in")
writevtk(Γbout,"bnd_out")
writevtk(Λin,"skel_in")
writevtk(Λout,"skel_out")

Γ⁺_act = BoundaryTriangulation(cutgeo_facets,ACTIVE_OUT,tags=["top"])
Γ⁺_cut = BoundaryTriangulation(cutgeo_facets,CUT_OUT,tags=["top"])

writevtk(Γ⁺_act,"bnd_actout_top")
writevtk(Γ⁺_cut,"bnd_cutout_top")



end # module
