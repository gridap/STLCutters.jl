module CellSubMeshes

using STLCutter

using STLCutter: @check

import STLCutter: distance

using Test

include("src/tables/lookup_cut_tables.jl");

struct CellSubMesh{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  dface_to_vertices::Vector{Vector{Vector{Int}}}
end

struct HexaCell{N,D,T}
  vertices::NTuple{N,Point{D,T}}
end

@generated function HexaCell(b::BoundingBox{D,T}) where {D,T}
  N = 2^D
  d = Dict( 0 => "pmin", 1 => "pmax" )
  v_str = [ "" for i in 1:N ]
  for i in 1:N
    bin = digits(i-1,base=2,pad=D)
    data = join([ "b.$(d[b])[$i]," for (i,b) in enumerate(bin) ])
    v_str[i] = "Point{$D,$T}($data),"
  end
  vertices = join(v_str)
  str = "HexaCell{$N,$D,$T}(($vertices))"
  Meta.parse(str)
end

function num_vertices(::Type{<:HexaCell{N,D}}) where {N,D}
  @check N == 2^D
  N
end

num_vertices(::T) where T<:HexaCell = num_vertices(T)

num_dims(::Type{<:HexaCell{N,D}}) where {N,D} = D

num_dims(::T) where T<:HexaCell = num_dims(T)

Base.getindex(c::HexaCell,i::Int) = c.vertices[i]

get_vertices(c::HexaCell) = c.vertices

function BoundingBox(c::HexaCell{N,D,T}) where {N,D,T}
  @check all( get_data(c[1]) .<= get_data( c[N] ) )
  BoundingBox{D,T}( c[1], c[N] )
end

function initialize(c::HexaCell{N,D,T}) where {N,D,T}
  coordinates = [ v for v in get_vertices(c) ]
  dface_to_vertices = [ Vector{Int}[] for d in 0:D ]

  for d in 0:D
    df_to_v = hex_to_tet_nface_to_mface[D][d+1][0+1]
    for i in 1:size(df_to_v,2)
      push!( dface_to_vertices[d+1], Int[] )
      for j in 1:size(df_to_v,1)
        push!( dface_to_vertices[d+1][i], df_to_v[j,i] )
      end
    end
  end
  CellSubMesh{D,T}(coordinates,dface_to_vertices)
end

function initialize!(m::CellSubMesh{D,T},c::HexaCell{N,D,T}) where {N,D,T}
  resize!(m.vertex_coordinates,N)
  resize!(m.dface_to_vertices,D+1)

  for d in 0:D
    df_to_v = hex_to_tet_nface_to_mface[D][d+1][0+1]
    resize!(m.dface_to_vertices[d+1],size(df_to_v,2))
    for i in 1:size(df_to_v,2)
      resize!(m.dface_to_vertices[d+1][i],size(df_to_v,1))
      for j in 1:size(df_to_v,1)
        m.dface_to_vertices[d+1][i][j] = df_to_v[j,i]
      end
    end
  end

  m
end

p1 = Point(1,2,3)
p2 = Point(4,5,6)

b=BoundingBox(p1,p2)

h = HexaCell(b)

b2 = BoundingBox(h)

@test b2 == b



m = initialize(h)

p1 = Point(1,2,4)
p2 = Point(5,3,5)

box = BoundingBox(p1,p2)

cell = HexaCell(box)

initialize!(m,cell)

@test @allocated(initialize!(m,cell))== 0



p0 = Point(0.0,0.0)
p1 = Point(1.0,1.0)
box = BoundingBox(p0,p1)
cell = HexaCell(box)

mesh = initialize(cell)

stl_point = Point(0.5,0.25)

const intersection_tolerance = 1e-10

num_dfaces(m::CellSubMesh,d::Int) = length(m.dface_to_vertices[d+1])

function dface_vertices(m::CellSubMesh,d::Int) 
  m.dface_to_vertices[d+1]
end

function vertex_coordinates(m::CellSubMesh) 
  m.vertex_coordinates
end

function get_0face(m::CellSubMesh{D,T},i::Int) where {D,T}
  vertex_coordinates(m)[i]
end

function get_1face(m::CellSubMesh{D,T},i::Int) where {D,T}
  l = dface_vertices(m,1)[i]
  v = vertex_coordinates(m)
  Segment(v[l[1]],v[l[2]])
end

function get_2face(m::CellSubMesh{D,T},i::Int) where {D,T}
  l = dface_vertices(m,2)[i]
  v = vertex_coordinates(m)
  Triangle(v[l[1]],v[l[2]],v[l[3]])
end

function get_3face(m::CellSubMesh{D,T},i::Int) where {D,T}
  l = dface_vertices(m,3)[i]
  v = vertex_coordinates(m)
  Tetrahedron(v[l[1]],v[l[2]],v[l[3]],v[l[4]])
end

function distance(p::Point{2},t::Triangle{2})
  if have_intersection(p,t)
    distance = 0.0
  else
    distance = typemax(0.0)
  end
end

function distance(m::CellSubMesh{D,T},d::Int,i::Int,p::Point{D,T}) where {D,T}
  @check d <= D
  if d == 0
    distance(p,get_0face(m,i))
  elseif d == 1
    distance(p,get_1face(m,i))
  elseif d == 2
    distance(p,get_2face(m,i))
  elseif d == 3
    distance(p,get_3face(m,i))
  else
    throw(ErrorException(""))
  end
end


function refine!(m::CellSubMesh{D,T},d::Int,i::Int,p::Point{D,T}) where {D,T}
  if d == 0
    return m
  else
    return m
  end
end


point = stl_point
D = 2
min_distance = intersection_tolerance 
idface = 0
for d in 0:D
  global min_distance
  global idface
  for i in 1:num_dfaces(mesh,d)
    println(distance(mesh,d,i,point))
    if distance(mesh,d,i,point) ≤ min_distance
      min_distance = distance(mesh,d,i,point)
      idface = i
    end
  end
  if idface != 0
    println("refine!(m,$d,$idface,p)")
    break
  end
end





#@show m

end # module
