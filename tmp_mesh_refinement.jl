module CellSubMeshes

using STLCutter

using WriteVTK

using STLCutter: @check

import STLCutter: distance

using Test

include("src/tables/lookup_cut_tables.jl");

struct CellSubMesh{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  dface_to_nfaces::Vector{Vector{Vector{Vector{Int}}}}
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
  dface_to_nfaces = [ Vector{Vector{Int}}[] for d in 0:D ]

  for d in 0:D
    dface_to_nfaces[d+1] = [ Vector{Int}[] for n in 0:d-1 ]
    for n in 0:d-1
      df_to_nf = hex_to_tet_nface_to_mface[D][d+1][n+1]
      for i in 1:size(df_to_nf,2)
        push!( dface_to_nfaces[d+1][n+1], Int[] )
        for j in 1:size(df_to_nf,1)
          push!( dface_to_nfaces[d+1][n+1][i], df_to_nf[j,i] )
        end
      end
    end
  end
  CellSubMesh{D,T}(coordinates,dface_to_nfaces)
end

function initialize!(m::CellSubMesh{D,T},c::HexaCell{N,D,T}) where {N,D,T}
  resize!(m.vertex_coordinates,N)
  resize!(m.dface_to_nfaces,D+1)

  for d in 0:D
    resize!(m.dface_to_nfaces,d+1,d)
    for n in 0:d-1
      df_to_nf = hex_to_tet_nface_to_mface[D][d+1][n+1]
      resize!(m.dface_to_nfaces[d+1],n+1,size(df_to_nf,2))
      for i in 1:size(df_to_nf,2)
        resize!(m.dface_to_nfaces[d+1][n+1],i,size(df_to_nf,1))
        for j in 1:size(df_to_nf,1)
          m.dface_to_nfaces[d+1][n+1][i][j] = df_to_nf[j,i]
        end
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

function num_vertices(m::CellSubMesh) 
  length(m.vertex_coordinates)
end

function num_dfaces(m::CellSubMesh,d::Int)
  if d == 0
    num_vertices(m)
  else
    length(m.dface_to_nfaces[d+1])
  end
end

function dface_vertices(m::CellSubMesh,d::Int) 
  m.dface_to_nfaces[d+1][0+1]
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

function Base.resize!(v::Vector{Vector{T}},i::Int,n::Int) where T
  if !isassigned(v,i)
    v[i] = T[]
  end
  resize!(v[i],n)
end

function add_dface!(m::CellSubMesh,d::Int,new_df_to_v::Vector{Vector{Int}}) 
  df_to_v = dface_vertices(m,d)
  resize!(df_to_v, length(df_to_v) + length(new_df_to_v) )
  for i in 1:length(new_df_to_v)
    k = i + length(df_to_v) - length(new_df_to_v) 
    resize!(df_to_v,k,length(new_df_to_v[i]))
    for j in 1:length(new_df_to_v[i])
      df_to_v[k][j] = new_df_to_v[i][j]
    end
  end
  m
end

function add_vertex!(m::CellSubMesh{D,T},p::Point{D,T}) where {D,T}
  push!( vertex_coordinates(m), p )
end

function delete_dface!(m::CellSubMesh,d::Int,i::Int) 
  deleteat!(m.dface_to_nfaces[d+1],i)
end

function refine!(m::CellSubMesh{D,T},d::Int,id::Int,p::Point{D,T}) where {D,T}
  if d == 0
    return m
  else
    case = cut_tet_dface_to_case[D][d][1]
    dface = dface_vertices(m,d)[id]
    vertices_cache = Int[]
    resize!(vertices_cache, length(dface) + 1 )
    vertices_cache[1:length(dface)] = dface
    vertices_cache[end] = num_dfaces(m,0) + 1
    nf_to_mf = cut_tet_nface_to_mface[D][case]
    df_to_v = nf_to_mf[d+1][0+1]
    #add_dface!()
    df_to_v_cache = Vector{Int}[]
    resize!(df_to_v_cache,size(df_to_v,2))
    for i in 1:size(df_to_v,2)
      resize!(df_to_v_cache,i,size(df_to_v,1))
      for j in 1:size(df_to_v,1)
        df_to_v_cache[i][j] = vertices_cache[ df_to_v[j,i] ] 
      end
    end
    delete_dface!(m,d,id) 
    add_vertex!(m,point)
    add_dface!(m,d,df_to_v_cache,)
    return m
  end
end

## TODO: do it for df_to_nf ∀ n ≤ d (instead of df_to_v)
# Adde more cube connectivities to reach this
# manage cache data in a dface_cutter
point = stl_point
D = 2
min_distance = intersection_tolerance 
idface = 0
for d in 0:D
  global min_distance
  global idface
  for i in 1:num_dfaces(mesh,d)
    if distance(mesh,d,i,point) ≤ min_distance
      min_distance = distance(mesh,d,i,point)
      idface = i
    end
  end
  if idface != 0
    refine!(mesh,d,idface,point)
    break
  end
end



function writevtk(m::CellSubMesh{D,T},file_base_name) where {D,T}
  d_to_vtk_type_id = Dict(0=>1,1=>3,2=>5) # tets
  num_points = num_vertices(m)
  points = zeros(T,D,num_points)
  for (i ,p ) in enumerate(vertex_coordinates(m)), d in 1:D
    points[d,i] = p[d]
  end
  cells = MeshCell{Vector{Int64}}[]
  for d in 1:D
    dface_to_vertices = dface_vertices(m,d)
    num_dfaces = length(dface_to_vertices)
    vtk_type = VTKCellType(d_to_vtk_type_id[d])
    for i in 1:num_dfaces
      push!( cells, MeshCell(vtk_type,dface_to_vertices[i]) )
    end
  end
  vtkfile = vtk_grid(file_base_name,points,cells)
  vtk_save(vtkfile)
end

writevtk(mesh,"sub_mesh")
# struct DFaceCutter end

end # module
