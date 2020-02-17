module TablesTests

using STLCutter
using Test

vectors = [[1,3],[3,4,5],[2,3,1],[1,]]

table = TableOfVectors(vectors)

faces_to_vertices = TableOfVectors(vectors)
nfaces = length(faces_to_vertices)
cache = table_cache(faces_to_vertices)
for face in 1:nfaces
  vertices = getlist!(cache,faces_to_vertices,face)
  @test vertices == vectors[face]
end
@test nfaces == length(vectors)

x=[(1,2),(3,5),(5,10),(3,4)]
n=7
t=TableOfLists(x,n)
@test getlist(t,1) == [2]
@test getlist(t,3) == [5,4]
@test getlist(t,5) == [10]
@test length(t) == n

t = TableOfVectors(Int,5,0)

@test length(t) == 5
@test length(getlist(t,1)) == 0

t = TableOfVectors(Int,5,2)

@test getlist(t,1) == [0,0]




cache = table_cache(cell_to_edges)
for cell in 1:length(cell_to_edges)
  edges = getlist!(cache,cell_to_edges,cell)
  for (ledge,edge) in enumerate(edges)
  end
end


for cell in 1:length(cell_to_edges)
  for ledge in 1:length(cell_to_edges,cell)
    edge = cell_to_edges[cell,ledge]
    cell_to_edges[cell,ledge] = edge
  end
end






end # module
