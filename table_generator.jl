using PyCall

scipy_spatial = pyimport("scipy.spatial")


function delaunay(x::Array{T,2}) where T
  scipy_spatial.Delaunay(x).simplices .+ 1
end

x=[0 0; 0 1; 1 0; 1 1]
c=delaunay(x)

x=[0 0 0; 0 1 0; 1 0 0; 0 0 1; 1 1 1]
c=delaunay(x)
