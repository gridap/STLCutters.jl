module Stokes

import AlgebraicMultigrid

using STLCutters
using Gridap
using GridapEmbedded
using Test
using Gridap

import Gridap: ∇
using GridapEmbedded
using Test
using LinearAlgebra: tr
using Gridap.ReferenceFEs

function main(filename;n,outputfile=nothing)

  geo = STLGeometry( filename )

  δ = 0.2
  pmin,pmax = get_bounding_box(geo)
  L = pmax-pmin
  zmin = pmin[3]+0.001*norm(L)
  pmin = pmin - L*δ - (1-δ)*L.*Point(0,1,0)
  pmax = pmax + L*δ + (4-δ)*L.*Point(0,1,0)
  pmin = VectorValue(pmin[1],pmin[2],zmin)

  # Background model
  cells = (4*n,n,n)
  bgmodel = CartesianDiscreteModel(pmin,pmax,cells)

  labels = get_face_labeling(bgmodel)
#  wall_tags = [(1:20)...,(21:24)...]
  wall_tags = union(get_faces(HEX)[ [21,22,25,26] ]...)
  #wall_tags = get_faces(HEX)[21]
  add_tag_from_tags!(labels,"wall",wall_tags)
  add_tag_from_tags!(labels,"inlet",23)

  # Cut the background model
  cutgeo,facet_to_inoutcut = cut(bgmodel,geo)

  Ωb = Triangulation(bgmodel)
  Ω = Triangulation(cutgeo,PHYSICAL_OUT)
  Ωa = Triangulation(cutgeo,ACTIVE_OUT)
  Γ = EmbeddedBoundary(cutgeo)
  n_Γ = -get_normal_vector(Γ)

  k = 2
  reffe_u = ReferenceFE(lagrangian,VectorValue{3,Float64},k)
  reffe_p = ReferenceFE(lagrangian,Float64,k-1,space=:P)
  Vstd = TestFESpace(Ωa,reffe_u,dirichlet_tags=["inlet","wall"])
  Qstd = FESpace(Ωa,reffe_p,conformity=:L2)

  threshold = 0.5
  strategy = AggregateCutCellsByThreshold(threshold)
  aggregates = aggregate(strategy,cutgeo,geo,OUT,facet_to_inoutcut)
  V = AgFEMSpace(Vstd,aggregates)
  Q = AgFEMSpace(Qstd,aggregates)

  #function u_z(z,zmin,zmax)
  #  L = zmax - zmin
  #  z = (z-zmin)/L
  #  fac = 10
  #  (1 - exp(-fac*z))
  #end
  #vmax = VectorValue(0.2,0.,0.)
  #u_in(x) = vmax*u_z(x[3],pmin[3],pmax[3])

  function u_1d(x,xmin,xmax)
    L = xmax - xmin
    y = (x-xmin)/L
    4*y-4*y^2
  end
  vmax = 0.2*VectorValue(0,1,0)
 # u_in(x) = vmax*u_1d(x[2],pmin[2],pmax[2])*u_1d(x[3],pmin[3],pmax[3])
  u_in(x) = vmax*u_1d(x[1],pmin[1],pmax[1])*u_1d(x[3],pmin[3],pmax[3])
  u_wall = VectorValue(0.,0.,0.)

  U = TrialFESpace(V,[u_in,u_wall])
  P = Q

  X = MultiFieldFESpace([U,P])
  Y = MultiFieldFESpace([V,Q])

  # perhaps for just a picture 2*k is enough also in cut cells
  dΩ = Measure(Ω,2*k) #2*(k-1),3*2*(k-1)
  dΓ = Measure(Γ,2*k) #2*3*k)

  β = 10
  h =  minimum(Tuple(pmax-pmin)./cells)
  τ = β*k^2/h

  a((u,p),(v,q)) =
    ∫( ∇(u)⊙∇(v) - p*(∇⋅v) - (∇⋅u)*q )dΩ +
    ∫(
      τ*v⋅u  -
      v⋅(n_Γ⋅∇(u)) -
      (n_Γ⋅∇(v))⋅u +
      (p*n_Γ)⋅v +
      (q*n_Γ)⋅u )dΓ

  l((v,q)) = 0

  op = AffineFEOperator(a,l,X,Y)
  uh,ph = solve(op)

  if !isnothing(outputfile)
    writevtk(Ωb,outputfile*"_Ωb")
    writevtk(Ω,outputfile*"_Ω",order=2,cellfields=["uh"=>uh,"ph"=>ph])
    writevtk(Ωa,outputfile*"_Ωa",cellfields=["uh"=>uh,"ph"=>ph])
    writevtk(Γ,outputfile*"_Γ",cellfields=["n"=>n_Γ,"uh"=>uh,"ph"=>ph])
  end


#
#
#
#  threshold = 1.0
#  strategy = AggregateByThreshold()
#  aggregates = aggregate(strategy,cutgeo,geo,OUT,facet_to_inoutcut)
#
#  Ω_bg = Triangulation(bgmodel)
#
#  # Generate the "active" mesh
#  Ω_act = Triangulation(cutgeo,ACTIVE_OUT)
#
#  colors = color_aggregates(aggregates,bgmodel)
#  writevtk(Ω_bg,"trian",celldata=["aggregate"=>aggregates,"color"=>colors])
#
#  # Setup integration meshes
#  Ω = Triangulation(cutgeo,PHYSICAL_OUT)
#
#  writevtk(Ω_act,"trian_act")
#  writevtk(Ω,"trian_phys")
#
##  Γw = EmbeddedBoundary(cutgeo)
#  Γg = GhostSkeleton(cutgeo,ACTIVE_OUT)
##
##  # Setup normal vectors
###  n_Γi = get_normal_vector(Γi)
##  n_Γw = get_normal_vector(Γw)
#  n_Γg = get_normal_vector(Γg)
#
#  #writevtk(Ω,"trian_O")
#  #writevtk(Γi,"trian_Gi",cellfields=["uin"=>uin,"normal"=>n_Γi])
#  #writevtk(Γw,"trian_Gw",cellfields=["normal"=>n_Γw])
#  #writevtk(Triangulation(bgmodel),"bgtrian")
#
#  # Setup Lebesgue measures
#  order = 1
#  degree = 2*order
#  dΩ = Measure(Ω,degree)
##  dΓi = Measure(Γi,degree)
##  dΓw = Measure(Γw,degree)
#  dΓg = Measure(Γg,degree)
#
#  # Setup FESpace
#
#  reffe_u = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
#  reffe_p = ReferenceFE(lagrangian,Float64,order)
#
#  Vstd = TestFESpace(Ω_act,reffe_u,labels=labels,dirichlet_tags=["boundary"])#,conformity=:H1)
#  Qstd = TestFESpace(Ω_act,reffe_p)#,conformity=:H1)
#
#  V = AgFEMSpace(Vstd,aggregates)
#  Q = AgFEMSpace(Qstd,aggregates)
#
##  V = Vstd
##  Q = Qstd
#
#  U = TrialFESpace(V,[uin])
#  P = TrialFESpace(Q)
#
#  X = MultiFieldFESpace([U,P])
#  Y = MultiFieldFESpace([V,Q])
#
#  # Stabilization parameters
#  β0 = 0.25
#  β1 = 0.2
#  β2 = 0.1
#  β3 = 0.05
#  γ = 10.0
#  h = (pmax-pmin)[1]/partition[1]
#
#  c_Ω(p,q) = (β1*h^2)*∇(p)⋅∇(q)
#  i_Γg(u,v) = (β2*h)*jump(n_Γg⋅∇(u))⋅jump(n_Γg⋅∇(v))
#  j_Γg(p,q) = (β3*h^3)*jump(n_Γg⋅∇(p))*jump(n_Γg⋅∇(q))
#
#  # Weak form
##  a((u,p),(v,q)) = ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) + c_Ω(p,q) )dΩ +
##    ∫( i_Γg(u,v) - j_Γg(p,q) ) * dΓg
##
#  a((u,p),(v,q)) = ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )dΩ
#  l((v,q)) = 0
#
#  # FE problem
#@time  op = AffineFEOperator(a,l,X,Y)
#
#@time  uh, ph = solve(op)
#
#  # Postprocess
#  if outputfile !== nothing
#    writevtk(Ω,outputfile, cellfields=["uh"=>uh,"ph"=>ph])
#  end
#
end

end # module
