module TablesTests

using STLCutters
using Test

data = [ 1, 2, 3, 4, 7 ]
rows = [ 2, 4, 1, 3, 2 ]
n = 5

t1 = Table(data,rows,n)

@test t1[1,1] == 3
@test t1[2,1] == 1
@test t1[2,2] == 7
@test length(t1) == n
@test length(t1,1) == 1
@test length(t1,2) == 2
@test length(t1,5) == 0

@test maximum(t1) == 7

@test eltype(t1) == eltype(Table{Int}) == Int

data = [ 3, 1, 7, 4, 2 ]
ptrs = [ 1, 2, 4, 5, 6, 6 ]

t2 = Table(data,ptrs)

@test t2[1,1] == 3
@test t2[2,1] == 1
@test t2[2,2] == 7
@test length(t2) == n
@test length(t2,1) == 1
@test length(t2,2) == 2
@test length(t2,5) == 0

@test t1 == t2

data = [ [3], [1,7], [4], [2], Int[] ]

t3 = Table(data)

@test t1 == t2 == t3

t = t1

push!(t,[1,2,3])
@test length(t) == n+1
@test length(t,n+1) == 3

@test isactive(t,2)

remove!(t,2)
@test !isactive(t,2)

compact!(t)
@test length(t) == n
@test length(t,2) == 1
@test t[2,1] == 4

data = [ 1, 2, 3, 4, 7 ]
rows = [ 2, 4, 1, 3, 2 ]
n = 5

t = Table(data,rows,n)
A = 
  [ 1 4 
    2 5 
    3 6 ]

append!(t,A)
@test length(t) == 5+2
@test length(t,n+1) == length(t,n+2) == 3
@test t[n+1,1] == 1
@test t[n+1,2] == 2
@test t[n+1,3] == 3

resize!(t,n) 

@test length(t) == n

v = Vector{Int}[ [ 1,2,3], [], [4,5] ]
append!(t,v)
@test length(t) == 5+3
@test length(t,n+1) == 3
@test length(t,n+2) == 0
@test length(t,n+3) == 2
@test t[n+1,1] == 1
@test t[n+1,2] == 2
@test t[n+1,3] == 3
@test t[n+3,1] == 4
@test t[n+3,2] == 5

resize!(t,n)

t2 = Table(v)
append!(t,t2)
@test length(t) == 5+3
@test length(t,n+1) == 3
@test length(t,n+2) == 0
@test length(t,n+3) == 2
@test t[n+1,1] == 1
@test t[n+1,2] == 2
@test t[n+1,3] == 3
@test t[n+3,1] == 4
@test t[n+3,2] == 5

resize!(t,0)
@test t == Table(eltype(t)) == Table{eltype(t)}()


end # module
