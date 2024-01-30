module CutFacetsTests

using Gridap
using GridapEmbedded
using STLCutters
using Test

geo = STLGeometry(joinpath(@__DIR__,"data/cube.stl"))
geoₐ = cube(x0=Point(0.0,0.0,0.5))

pmin = Point(-1.0, -1.0, 0.5)
pmax = Point(1.0, 1.0, 1.5)
n = 5
nz = 3
partition = (n,n,nz)

model = CartesianDiscreteModel(pmin, pmax, partition)
cutgeo_facets = cut_facets(model, geo)
cutgeo_facetsₐ = cut_facets(model, geoₐ)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"top",[21])

Γbin = BoundaryTriangulation(cutgeo_facets,CUT_IN)
Γbout = BoundaryTriangulation(cutgeo_facets,CUT_OUT)
Λin =  SkeletonTriangulation(cutgeo_facets,CUT_IN)
Λout = SkeletonTriangulation(cutgeo_facets,CUT_OUT)

Γbinₐ = BoundaryTriangulation(cutgeo_facetsₐ,CUT_IN)
Γboutₐ = BoundaryTriangulation(cutgeo_facetsₐ,CUT_OUT)
Λinₐ =  SkeletonTriangulation(cutgeo_facetsₐ,CUT_IN)
Λoutₐ = SkeletonTriangulation(cutgeo_facetsₐ,CUT_OUT)

Γt = BoundaryTriangulation(model,tags=["top"])

Γ⁺_act = BoundaryTriangulation(cutgeo_facets,ACTIVE_OUT,tags=["top"])
Γ⁺_cut = BoundaryTriangulation(cutgeo_facets,CUT_OUT,tags=["top"])
Γ⁺_phy = BoundaryTriangulation(cutgeo_facets,PHYSICAL_OUT,tags=["top"])

Γ⁻_act = BoundaryTriangulation(cutgeo_facets,ACTIVE_IN,tags=["top"])
Γ⁻_cut = BoundaryTriangulation(cutgeo_facets,CUT_IN,tags=["top"])
Γ⁻_phy = BoundaryTriangulation(cutgeo_facets,PHYSICAL_IN,tags=["top"])

Γ⁺_actₐ = BoundaryTriangulation(cutgeo_facetsₐ,ACTIVE_OUT,tags=["top"])
Γ⁺_cutₐ = BoundaryTriangulation(cutgeo_facetsₐ,CUT_OUT,tags=["top"])
Γ⁺_phyₐ = BoundaryTriangulation(cutgeo_facetsₐ,PHYSICAL_OUT,tags=["top"])

Γ⁻_actₐ = BoundaryTriangulation(cutgeo_facetsₐ,ACTIVE_IN,tags=["top"])
Γ⁻_cutₐ = BoundaryTriangulation(cutgeo_facetsₐ,CUT_IN,tags=["top"])
Γ⁻_phyₐ = BoundaryTriangulation(cutgeo_facetsₐ,PHYSICAL_IN,tags=["top"])

degree = 2
dΓt = Measure(Γt,degree)

dΓ⁺_act = Measure(Γ⁺_act,degree)
dΓ⁺_cut = Measure(Γ⁺_cut,degree)
dΓ⁺_phy = Measure(Γ⁺_phy,degree)
dΓ⁻_act = Measure(Γ⁻_act,degree)
dΓ⁻_cut = Measure(Γ⁻_cut,degree)
dΓ⁻_phy = Measure(Γ⁻_phy,degree)

dΓ⁺_actₐ = Measure(Γ⁺_actₐ,degree)
dΓ⁺_cutₐ = Measure(Γ⁺_cutₐ,degree)
dΓ⁺_phyₐ = Measure(Γ⁺_phyₐ,degree)
dΓ⁻_actₐ = Measure(Γ⁻_actₐ,degree)
dΓ⁻_cutₐ = Measure(Γ⁻_cutₐ,degree)
dΓ⁻_phyₐ = Measure(Γ⁻_phyₐ,degree)

dΓbin = Measure(Γbin,degree)
dΓbout = Measure(Γbout,degree)
dΛin⁺ = Measure(Λin.⁺,degree)
dΛin⁻ = Measure(Λin.⁻,degree)
dΛout⁺ = Measure(Λout.⁺,degree)
dΛout⁻ = Measure(Λout.⁻,degree)

dΓbinₐ = Measure(Γbinₐ,degree)
dΓboutₐ = Measure(Γboutₐ,degree)
dΛinₐ⁺ = Measure(Λinₐ.⁺,degree)
dΛinₐ⁻ = Measure(Λinₐ.⁻,degree)
dΛoutₐ⁺ = Measure(Λoutₐ.⁺,degree)
dΛoutₐ⁻ = Measure(Λoutₐ.⁻,degree)

# test area vs. known value
@test sum(∫(1)dΓt) ≈ 4
@test sum(∫(1)dΓ⁺_phy) + sum(∫(1)dΓ⁻_phy) ≈ 4
@test sum(∫(1)dΓ⁺_phy) ≈ 3
@test sum(∫(1)dΓ⁺_act) - sum(∫(1)dΓ⁻_cut) ≈ 3
@test sum(∫(1)dΓ⁻_phy) ≈ 1
@test sum(∫(1)dΓ⁻_act) - sum(∫(1)dΓ⁺_cut) ≈ 1

# test STLCutter vs GridapEmbedded
@test sum(∫(1)dΓ⁺_phy) ≈ sum(∫(1)dΓ⁺_phyₐ)
@test sum(∫(1)dΓ⁺_act) ≈ sum(∫(1)dΓ⁺_actₐ)
@test sum(∫(1)dΓ⁺_cut) ≈ sum(∫(1)dΓ⁺_cutₐ)
@test sum(∫(1)dΓ⁻_phy) ≈ sum(∫(1)dΓ⁻_phyₐ)
@test sum(∫(1)dΓ⁻_phy) ≈ sum(∫(1)dΓ⁻_phyₐ)
@test sum(∫(1)dΓ⁻_phy) ≈ sum(∫(1)dΓ⁻_phyₐ)
@test sum(∫(1)dΓbin) ≈ sum(∫(1)dΓbinₐ)
@test sum(∫(1)dΓbout) ≈ sum(∫(1)dΓboutₐ)
@test sum(∫(1)dΛin⁺) ≈ sum(∫(1)dΛinₐ⁺)
@test sum(∫(1)dΛin⁻) ≈ sum(∫(1)dΛinₐ⁻)
@test sum(∫(1)dΛout⁺) ≈ sum(∫(1)dΛoutₐ⁺)
@test sum(∫(1)dΛout⁻) ≈ sum(∫(1)dΛoutₐ⁻)


using Gridap.ReferenceFEs
cutgeo = cut(model, geo)
Ω = Triangulation(cutgeo,PHYSICAL)
Ω_act_in = Triangulation(cutgeo,ACTIVE_IN,geo)
Ω_act_out = Triangulation(cutgeo,ACTIVE_OUT,geo)
dΩᵐ_in = Measure(Ω_act_in,Quadrature(momentfitted,cutgeo,degree,in_or_out=IN))
dΩᵐ_out = Measure(Ω_act_out,Quadrature(momentfitted,cutgeo,degree,in_or_out=OUT))

dΩ = Measure(Ω,degree)
f = x -> x[1] + 1
f = 1
# @test ∑(∫(f)dΩᵐ_in) ≈ ∑(∫(f)dΩ)


end # module
