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

t = TableOfVectors(Int,5,0)

@test length(t) == 5
@test length(getlist(t,1)) == 0

t = TableOfVectors(Int,5,2)

@test getlist(t,1) == [0,0]


end # module
