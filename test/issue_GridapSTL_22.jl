
# module embedded_2d_CYLINDER
using Gridap
using Gridap.Geometry
using Gridap.TensorValues
using Gridap.Fields
using GridapEmbedded.Interfaces
# using WriteVTK
using GridapEmbedded
using GridapEmbedded.LevelSetCutters
using Gridap.Arrays
using Gridap.Adaptivity
using STLCutters



# MWE

geo = STLGeometry("test/data/cube.stl")

pmin,pmax = get_bounding_box(geo)
L = pmax-pmin

pmin -= L/2
pmax += L/2
pmin = Point(pmin[1],pmin[2],0.5)


n = 5
partition = (n,n,n)
model = CartesianDiscreteModel(pmax, pmin, partition)

cutgeo = cut(model, geo)


# TODO: do tests


# TODO: reduce to a MWE to test the boundary
aaaaa




# variables

g = 9.81        # [kg/s²] gravitational constant
L = 2.0        # [m] length of domain
d = 1.0        # [m] depth
nx = 5        # [-] number of elements horizontally
ny = 5
nz = 5         # [-] number of elements vertically


# background model ~ reference domain
pmin = Point(L, L, 0.0)
pmax = Point(0.0, 0.0, -d)
partition = (nx, ny, nz)

# assign labels to mesh surfaces and corner points
# labels = get_face_labeling(model)
# add_tag_from_tags!(labels, "seabed", [13, 14, 15, 16])
# add_tag_from_tags!(labels, "surface", [1, 2, 3, 4, 5,6,7,8, 9])
# add_tag_from_tags!(labels, "inlet", [3, 7])
# add_tag_from_tags!(labels, "outlet", [4, 8])
# add_tag_from_tags!(labels, "water", [9])

# initial conditions (zero displacement & velocity potential)
uᵢ(x,t) = R/3               # [m] heave displacement Ito 1977
uᵢ(t::Real) = x -> uᵢ(x,t)
ϕᵢ(x,t) = 0.0               # [m²/s] velocity potential
ϕᵢ(t::Real) = x -> ϕᵢ(x,t)
ηᵢ(x,t) = 0.0               # [m] wave height
ηᵢ(t) = x -> ηᵢ(x,t)

# define geometry & apply embedded boundary
# geo = STLGeometry("/Users/jmodderman1/Documents/meshes/shrunk2.stl")
# pmin,pmax = get_bounding_box(geo)
# δ = 0.2
# diagonal = pmax-pmin
# pmin = pmin - diagonal*δ
# pmax = pmax + diagonal*δ
model = CartesianDiscreteModel(pmax, pmin, partition)
pmid = 0.5*(pmin+pmax)+VectorValue(0.0,0.0,d/2)
# geo = sphere(0.5;x0=pmid,name="sphere")
geo = STLGeometry("/Users/jmodderman1/Documents/meshes/spherestl_debug.stl")




labels = get_face_labeling(model)
wall_tags = union(STLCutters.get_faces(HEX)[ [23,24] ]...)
add_tag_from_tags!(labels,"wall",wall_tags)
add_tag_from_tags!(labels,"seabed",[21])
add_tag_from_tags!(labels,"surface",[22])
add_tag_from_tags!(labels,"inlet",[25])
add_tag_from_tags!(labels,"outlet",[26])





# writevtk(geo,"data/sims/issue_22/geo")

# cutgeo, facet_to_inoutcut = cut(model, geo)
# cutgeo_facets = STLCutters._cut_stl(model, geo)
cutgeo = cut(model, geo)

println("CUTGEO")
println(cutgeo.ls_to_bgcell_to_inoutcut)
# println(cutgeo.subcells.cell_to_points)
# println(cutgeo.subcells.cell_to_bgcell)
println(cutgeo.ls_to_subcell_to_inout)
println(cutgeo.subfacets)
println(cutgeo.ls_to_subfacet_to_inout)
# println(cutgeo.oid_to_ls)


Ω = Interior(model)
writevtk(model,"data/sims/issue_22/model")

Ω⁻ = Interior(cutgeo, PHYSICAL_OUT)
writevtk(Ω⁻,"data/sims/issue_22/omgmin")

Ω⁻act = Interior(cutgeo, ACTIVE_OUT)
Γ = EmbeddedBoundary(cutgeo)
n_Γ = -get_normal_vector(Γ)

facets = BoundaryTriangulation(cutgeo.bgmodel;tags=["surface"])
writevtk(facets,"data/sims/issue_22/Gammaf")
#Γf = BoundaryTriangulation(cutgeo, PHYSICAL_OUT, tags=["surface"])
cutgeo_facets = cut_facets(model, geo)
println("FACETS")
println(cutgeo_facets.ls_to_facet_to_inoutcut)
# println(cutgeo_facets.subfacets.cell_to_points)
# println(cutgeo_facets.subfacets.cell_to_bgcell)
println(cutgeo_facets.ls_to_subfacet_to_inout)
# println(cutgeo_facets.oid_to_ls)


stop
Γw = BoundaryTriangulation(model, tags=["wall"])
Γi = BoundaryTriangulation(model, tags=["inlet"])
Γo = BoundaryTriangulation(model, tags=["outlet"])
Γs = BoundaryTriangulation(model, tags=["seabed"])

writevtk(Γf,"data/sims/issue_22/Gammaf")
writevtk(Γw,"data/sims/issue_22/Gammaw")
writevtk(Γi,"data/sims/issue_22/Gammai")
writevtk(Γo,"data/sims/issue_22/Gammao")
writevtk(Γs,"data/sims/issue_22/Gammas")


writevtk(Γ,"data/sims/issue_22/Gamma",cellfields=["normal"=>n_Γ])

# create integration space (triangulation) & Gauss quadratures (measure)
degree = 2
#cutgeo_facets = cut_facets(model, geo)
#Γf⁻ = BoundaryTriangulation(cutgeo_facets, PHYSICAL, tags=["surface"])
#Γf_act = BoundaryTriangulation(cutgeo_facets, ACTIVE, tags=["surface"])
#Γf = BoundaryTriangulation(model, tags=["surface"])
#dΓf_act = Measure(Γf_act, degree)
#dΓf⁻ = Measure(Γf⁻, degree)
#n_Γf = get_normal_vector(Γf)
#writevtk(Γf,"GridapPF.jl/data/sims/linPF_gmsh_FSI_3d/Gammaf")


# copied the degree of 8 from transientstokes - set to 2
dΩ⁻ = Measure(Ω⁻, 2)
dΩ = Measure(Ω, 2)
dΓ = Measure(Γ, 2)

# definition of FE spaces
order = 1
reffeᵤ = ReferenceFE(lagrangian, Float64, order)
reffeᵩ = ReferenceFE(lagrangian, Float64, order)

stop

# AgFEM
Vstd = FESpace(Γf_act, reffeᵤ)
Ustd = TransientTrialFESpace(Vstd)
Wstd = FESpace(Ω⁻act, reffeᵩ)
Φstd = TransientTrialFESpace(Wstd)
Dstd = ConstantFESpace(model)
Rstd = TransientTrialFESpace(Dstd)

strategy = AggregateAllCutCells()
aggregates = aggregate(strategy, cutgeo, facet_to_inoutcut)

# final FE spaces
W = AgFEMSpace(Wstd, aggregates)
Φ = TransientTrialFESpace(W)
X = TransientMultiFieldFESpace([Φ, Ustd, Rstd])
Y = MultiFieldFESpace([W, Vstd, Dstd])


# checks on wetted and projected surface
wet_area = sum(∫(1.0)dΓ)
println("Wet area: ", wet_area , " m")
print("    Should approximate πR: ", π*R, " m")

projected_area = sum(∫(n_Γ⋅VectorValue(0.0,1.0))dΓ)
println("Wet area projection at Free Surface", projected_area, " m")
print("    Should approximate 2R: ", 2*R, " m")

# waterline
z(x, t) = 0.0   # [m] waterline
z(t) = x -> z(x, t)

ρ = 1000    # [kg/m³] density water
ρᵦ = 500    # [kg/m³] density cylinder (half that of water)
R=0.1

m_cyl = π*R^2 * ρᵦ/ρ /projected_area    # [kg/m] mass
a_cyl = 0.0                             # [kg/m] hydrodynamic mass coefficient
b_cyl = 0.0                             # [kg/ms] hydrodynamic damping coefficient
c_cyl = g                               # [kg/s²] restoring spring coefficient

println("Radius: ", R, " m")
println("Initial displacement: ", uᵢ(VectorValue(0.0,0.0),0), " m")
println("Mass: ", m_cyl, " kg/m")
println("Stiffness: ", c_cyl, " kg/s²")
println("Frequency: ", √(c_cyl/m_cyl / (4*π^2)), " Hz")

jac(t, (ϕ,η,u), (dϕ,dη,du), (w,v,r)) = ∫( ∇(dϕ)⋅∇(w) )dΩ⁻ +
                                        ∫( v * g * dη )dΓf⁻ +
                                        ∫( (1.0) * c_cyl * r * du * (VectorValue(0.0,1.0)⋅n_Γ))dΓ

jac_t(t, (ϕ,η,u), (dϕₜ,dηₜ,duₜ), (w,v,r)) = ∫( v * dϕₜ )dΓf⁻ -
                                        ∫( dηₜ * w )dΓf⁻ -
                                        ∫( w * duₜ * (VectorValue(0.0,1.0)⋅n_Γ) )dΓ +
                                        ∫((1.0) * b_cyl * duₜ * r * (VectorValue(0.0,1.0)⋅n_Γ))dΓ -
                                        ∫((-1.0) * dϕₜ * r * (VectorValue(0.0,1.0)⋅n_Γ))dΓ

jac_tt(t, (ϕ,η,u), (dϕₜₜ,dηₜₜ,duₜₜ), (w,v,r)) = ∫((1.0) * duₜₜ * r  * m_cyl * (VectorValue(0.0,1.0)⋅n_Γ))dΓ

res(t, (ϕ, η, u), (w, v, r)) = ∫((1.0) * ∂tt(u) * r * m_cyl * (VectorValue(0.0,1.0)⋅n_Γ))dΓ +
                                ∫( v * ∂t(ϕ) )dΓf⁻ -
                                ∫( ∂t(η) * w )dΓf⁻ -
                                ∫( w * ∂t(u) * (VectorValue(0.0,1.0)⋅n_Γ) )dΓ +
                                ∫((1.0) * b_cyl * ∂t(u) * r * (VectorValue(0.0,1.0)⋅n_Γ))dΓ -
                                ∫((-1.0) * ∂t(ϕ) * r * (VectorValue(0.0,1.0)⋅n_Γ))dΓ +
                                ∫( ∇(ϕ)⋅∇(w) )dΩ⁻ +
                                ∫( v * g * η )dΓf⁻ +
                                ∫( (1.0) * c_cyl * r * u * (VectorValue(0.0,1.0)⋅n_Γ))dΓ -
                                ∫( (1.0) * r * c_cyl * z(t) * (VectorValue(0.0,1.0)⋅n_Γ))dΓ

# time integrator variables
Δt = 1/50                          # [s]: time step
t₀ = 0.0                            # [s]: zero time
Tf = 5*Δt                           # [s]: final time
γn = 0.5                            # [-]: γ factor Newmark
βn = 0.25                           # [-]: β factor Newmark

# solver definition
ls = LUSolver()
op = TransientFEOperator(res, jac, jac_t, jac_tt, X, Y)
ode_solver = Newmark(ls, Δt, γn, βn)

# interpolate initial conditions
x₀ = interpolate_everywhere([ϕᵢ(0.0), ηᵢ(0.0), uᵢ(0.0)], X(0.0))
v₀ = interpolate_everywhere([ϕᵢ(0.0), ηᵢ(0.0), ηᵢ(0.0)], X(0.0))
a₀ = interpolate_everywhere([ϕᵢ(0.0), ηᵢ(0.0), ηᵢ(0.0)], X(0.0))
stop
# solve
ϕhₜ = Gridap.solve(ode_solver, op, (x₀, v₀, a₀), t₀, Tf)

# pvd = createpvd("GridapPF.jl/data/sims/linPF_gmsh_FSI_3d/test_case_embeddedcyl1")
# pvd1 = createpvd("GridapPF.jl/data/sims/linPF_gmsh_FSI_3d/test_case_embeddedcyl1_g")
# pvd2 = createpvd("GridapPF.jl/data/sims/linPF_gmsh_FSI_3d/test_case_embeddedcyl1_u")


#     pvd[0] = createvtk(Ω⁻, "GridapPF.jl/data/sims/linPF_gmsh_FSI_3d/test_case_embeddedcyl1_0.vtu", cellfields=[ "phih"=>ϕᵢ(0.0) ])
#     pvd1[0] = createvtk(Γf⁻, "GridapPF.jl/data/sims/linPF_gmsh_FSI_3d/test_case_embeddedcyl1_g_0.vtu", cellfields=["etah"=>ηᵢ(0.0)])
#     pvd2[0] = createvtk(Γ, "GridapPF.jl/data/sims/linPF_gmsh_FSI_3d/test_case_embeddedcyl1_u_0.vtu", cellfields=["uh"=>uᵢ(0.0)])

#     for ((ϕh,ηh,uh),t) in ϕhₜ
#       pvd[t] = createvtk(Ω⁻,"GridapPF.jl/data/sims/linPF_gmsh_FSI_3d/test_case_embeddedcyl1_$t"*".vtu",cellfields=["phih"=>ϕh])
#       pvd1[t] = createvtk(Γf⁻,"GridapPF.jl/data/sims/linPF_gmsh_FSI_3d/test_case_embeddedcyl1_g_$t"*".vtu",cellfields=["etah"=>ηh])
#       pvd2[t] = createvtk(Γ,"GridapPF.jl/data/sims/linPF_gmsh_FSI_3d/test_case_embeddedcyl1_u_$t"*".vtu",cellfields=["uh"=>uh])
#     end
# vtk_save(pvd)
# vtk_save(pvd1)
# vtk_save(pvd2)
# end
