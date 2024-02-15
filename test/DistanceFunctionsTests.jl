module DistanceFunctionsTests

using Test
using STLCutters
using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs

using STLCutters: closest_point

using STLCutters: compute_cell_to_facets
using STLCutters: compose_index_map
using STLCutters: get_stl

using STLCutters: STL

filename = joinpath(@__DIR__,"data/cube.stl")
geo = STLGeometry(filename)

x = [
  Point(1,1,1),
  Point(0,0,0),
  Point(0,0,0.1),
  Point(1,1,0.5)]
xp = closest_point(x,geo)

@test xp[1] == Point(0.5,0.5,1)
@test xp[2] == Point(0,0,0)
@test xp[3] == Point(0,0,0)
@test xp[4] == Point(0.5,0.5,0.5)


# Setup points from CartesianGrid
n = 5
pmin = Point(-0.7,-0.7,-0.2)
pmax = Point(0.7,0.7,1.2)
bgmodel = CartesianDiscreteModel(pmin,pmax,(n,n,n))
topo = get_grid_topology( bgmodel)
D = num_dims(topo)
X = get_node_coordinates(bgmodel)

# Compute point_to_facets
point_to_cells = get_faces(topo,0,D)
cell_to_facets = compute_cell_to_facets(bgmodel,get_stl(geo))
point_to_facets = compose_index_map(point_to_cells,cell_to_facets)

# Filter only points with facets around
newpoints = findall(!isempty,point_to_facets)
X = X[newpoints]
point_to_facets = point_to_facets[newpoints]

Xc = closest_point(X,geo,point_to_facets)
dist = map(-,X,Xc)
D = map(norm,dist)
N = map(/,dist,D)

H = Tuple(pmax-pmin) ./ (n,n,n)

@test maximum(D) < norm(H)

writevtk(X,"points")
writevtk(Xc,"closest_points")
writevtk(geo,"geo")


end # module
