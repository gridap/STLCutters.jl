module GridRefinementTests
  
using Test
using Gridap
using STLCutters
using DelimitedFiles

using Gridap.Geometry

using STLCutters: compute_stl_model
using STLCutters: refine_grid 
using STLCutters: volume,volumes
using STLCutters: FACE_CUT, FACE_IN, FACE_OUT
using STLCutters: get_bounding_box 


vertices = [
  Point(-0.1,0.5),
  Point(0.1,0.8),
  Point(0.4,0.4),
  Point(0.7,0.7),
  Point(0.8,0.3),
  Point(1.1,0.3) ]
faces = [[1,2],[2,3],[3,4],[4,5],[5,6]]

stl = compute_stl_model(Table(faces),vertices)
bg_mesh = compute_linear_grid(QUAD4)
bg_mesh = CartesianGrid( Point(0.,0.), (0.3,0.3), (3,3) )

out = refine_grid(bg_mesh,stl)
T,X,reffes,cell_types,cell_to_io,cell_to_bgcell,bgcell_to_ioc = out 

submesh = UnstructuredGrid(X,Table(T),reffes,cell_types)
vol = volume(submesh) + sum(volumes(bg_mesh).*(bgcell_to_ioc.≠FACE_CUT))
@test volume(bg_mesh) ≈ vol

writevtk(submesh,"submesh",cellfields=["io"=>cell_to_io,"bgcell"=>cell_to_bgcell])
writevtk(bg_mesh,"bgmesh",cellfields=["io"=>bgcell_to_ioc])
writevtk(stl,"stl")

vertices = [
  Point(-2.0,-2.0,3.0),
  Point(0.4,-2.0,0.3),
  Point(0.4,2.0,0.3),
  Point(3.0,3.0,3.0) ]
facets = [[1,2,3],[2,4,3]]

stl = compute_stl_model(Table(facets),vertices)
bg_mesh = compute_linear_grid(HEX8)
bg_mesh = CartesianGrid( Point(0.,0.,0.), (0.3,0.3,0.3), (3,3,3) )

out = refine_grid(bg_mesh,stl)
T,X,reffes,cell_types,cell_to_io,cell_to_bgcell,bgcell_to_ioc = out 

submesh = UnstructuredGrid(X,Table(T),reffes,cell_types)
vol = volume(submesh) + sum(volumes(bg_mesh).*(bgcell_to_ioc.≠FACE_CUT))
@test volume(bg_mesh) ≈ vol

writevtk(submesh,"submesh",cellfields=["io"=>cell_to_io,"bgcell"=>cell_to_bgcell])
writevtk(bg_mesh,"bgmesh",cellfields=["io"=>bgcell_to_ioc])
writevtk(stl,"stl")


vertices = [
  Point(-0.1,0.5,0.5),
  Point(0.5,1.1,0.5),
  Point(1.1,0.5,0.5),
  Point(0.5,0.1,0.6),
  Point(1.1,1.1,0.6),
  Point(-0.1,1.1,0.5),
  Point(-0.1,-0.1,0.3),
  Point(0.5,-0.1,0.2),
  Point(1.1,-0.1,0.3) ]
facets = [[1,2,3],[1,3,4],[1,6,2],[3,2,5],[1,4,7],[4,8,7],[4,9,8],[3,9,4]]

stl = compute_stl_model(Table(facets),vertices)
bg_mesh = compute_linear_grid(HEX8)
bg_mesh = CartesianGrid( Point(0.,0.,0.), (0.3,0.3,0.3), (3,3,3) )

out = refine_grid(bg_mesh,stl)
T,X,reffes,cell_types,cell_to_io,cell_to_bgcell,bgcell_to_ioc = out 

submesh = UnstructuredGrid(X,Table(T),reffes,cell_types)
vol = volume(submesh) + sum(volumes(bg_mesh).*(bgcell_to_ioc.≠FACE_CUT))
@test volume(bg_mesh) ≈ vol

writevtk(submesh,"submesh",cellfields=["io"=>cell_to_io,"bgcell"=>cell_to_bgcell])
writevtk(bg_mesh,"bgmesh",cellfields=["io"=>bgcell_to_ioc])
writevtk(stl,"stl")


function read_vertices(filename::String)
  D = 2
  T = Float64
  if endswith(filename,".csv") || endswith(filename,".CSV") 
    data = readdlm(filename,',',T,skipstart=1)
  else
    data = readdlm(filename,T,skipstart=1)
    @assert size(data,2) == D "read_vertices() only supports 2D coordinates to define a polyline"
  end
  vertices = zeros(Point{D,T},size(data,1))
  for i in 1:size(data,1)
    vertices[i] = Point(data[i,1],data[i,2])
  end
  vertices
end


function degrees_to_kilometers(a::Point{2})
  φ = deg2rad(a[2])
  y_km_per_deg = 111132.92 - 559.82*cos(2φ) + 1.175*cos(4φ) - 0.0023*cos(6φ)
  x_km_per_deg = 111412.85*cos(φ) - 93.5*cos(3φ) + 0.118*cos(5φ)
  x_km_per_deg,y_km_per_deg
end

∘(a::Tuple,b::Tuple) = a.*b

file = joinpath(@__DIR__,"data/mallorca.csv")
vertices = read_vertices(file)
vertices = [ vertices[i] for i in 1:2:length(vertices) ]
deleteat!(vertices,[1659,1660])
faces = [ [ i , i==length(vertices) ? 1 : i+1 ] for i in 1:length(vertices) ]

center = sum(vertices)/length(vertices)
vertices = Point.( Tuple.(vertices.-center) .∘ degrees_to_kilometers.(vertices) )


stl = compute_stl_model(Table(faces),vertices)

writevtk(stl,"mallorca")
n = 100
pmin,pmax = get_bounding_box(stl)
diagonal = pmax-pmin
origin = pmin - diagonal*0.1
sizes = Tuple( diagonal*1.2/n )
partion = (n,n)

bg_mesh = CartesianGrid(origin,sizes,partion)

using STLCutters: compute_face_to_cells
using STLCutters: get_face_lists 
using STLCutters: init_cell_mesh 
using STLCutters: distribute_vertices, distribute_facets 
using STLCutters: define_cells! 

grid = bg_mesh
cell = 1 
_,cell_to_stlf = compute_face_to_cells(grid,stl)

reffes = [QUAD4,TRI3]
V,F = get_face_lists(stl,cell,cell_to_stlf)
if length(F) > 0
Tk,Xk,p = init_cell_mesh(grid,cell) 
fk = [Int[]]
Tknew,fknew = empty(Tk), empty(fk)
VTk = distribute_vertices(Tk,Xk,p,stl,V)
insert_vertices!(Tk,Xk,p,stl,VTk,fk,Tknew,fknew)
Tk,fk = Tknew,fknew
vgrid = compute_grid(Table(Tk),Xk,QUAD)
writevtk(vgrid,"vertex_grid")
Tknew,fknew = empty(Tk), empty(fk)
k_cell_types,k_cell_to_io = Int8[], Int8[]
FTk = distribute_facets(Tk,Xk,p,stl,F)
insert_facets!(Tk,Xk,p,stl,FTk,fk,Tknew,fknew,k_cell_types,k_cell_to_io)
Tk,fk = Tknew,fknew
grid_k = UnstructuredGrid(Xk,Table(Tk),reffes,k_cell_types)

define_cells!(grid_k,fk,k_cell_to_io)

writevtk(grid_k,"facet_grid",cellfields=["io"=>k_cell_to_io])
  
end

writevtk(bg_mesh,"mesh")
out = refine_grid(bg_mesh,stl)
T,X,reffes,cell_types,cell_to_io,cell_to_bgcell,bgcell_to_ioc = out 

submesh = UnstructuredGrid(X,Table(T),reffes,cell_types)

writevtk(submesh,"submesh",cellfields=["io"=>cell_to_io,"bgcell"=>cell_to_bgcell])
writevtk(bg_mesh,"bgmesh",cellfields=["io"=>bgcell_to_ioc])



end # module
