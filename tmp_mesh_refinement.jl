module CellSubMeshes

using STLCutter

using WriteVTK

using STLCutter: @check

import STLCutter: distance

using Test

include("src/tables/lookup_cut_tables.jl");

function Base.resize!(v::Vector{Vector{T}},i::Int,n::Int) where T
  if !isassigned(v,i)
    v[i] = T[]
  end
  resize!(v[i],n)
end

struct CellSubMesh{D,T}
  vertex_coordinates::Vector{Point{D,T}}
  dface_to_nfaces::Vector{Vector{Vector{Vector{Int}}}}
end

struct Hexahedron{N,D,T}
  vertices::NTuple{N,Point{D,T}}
end

@generated function Hexahedron(b::BoundingBox{D,T}) where {D,T}
  N = 2^D
  d = Dict( 0 => "pmin", 1 => "pmax" )
  v_str = [ "" for i in 1:N ]
  for i in 1:N
    bin = digits(i-1,base=2,pad=D)
    data = join([ "b.$(d[b])[$i]," for (i,b) in enumerate(bin) ])
    v_str[i] = "Point{$D,$T}($data),"
  end
  vertices = join(v_str)
  str = "Hexahedron{$N,$D,$T}(($vertices))"
  Meta.parse(str)
end

function num_vertices(::Type{<:Hexahedron{N,D}}) where {N,D}
  @check N == 2^D
  N
end

num_vertices(::T) where T<:Hexahedron = num_vertices(T)

num_dims(::Type{<:Hexahedron{N,D}}) where {N,D} = D

num_dims(::T) where T<:Hexahedron = num_dims(T)

Base.getindex(c::Hexahedron,i::Int) = c.vertices[i]

get_vertices(c::Hexahedron) = c.vertices

function BoundingBox(c::Hexahedron{N,D,T}) where {N,D,T}
  @check all( get_data(c[1]) .<= get_data( c[N] ) )
  BoundingBox{D,T}( c[1], c[N] )
end

function initialize(c::Hexahedron{N,D,T}) where {N,D,T}
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

function initialize!(m::CellSubMesh{D,T},c::Hexahedron{N,D,T}) where {N,D,T}
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

const intersection_tolerance = 1e-10

function num_vertices(m::CellSubMesh) 
  length(m.vertex_coordinates)
end

function num_dfaces(m::CellSubMesh,d::Int)
  if d == 0
    num_vertices(m)
  else
    length(m.dface_to_nfaces[d+1][0+1])
  end
end

function dface_vertices(m::CellSubMesh,d::Int) 
  m.dface_to_nfaces[d+1][0+1]
end

function vertex_coordinates(m::CellSubMesh) 
  m.vertex_coordinates
end

function get_dface(m::CellSubMesh{D,T},::Val{0},i::Int) where {D,T}
  vertex_coordinates(m)[i]
end

function get_dface(m::CellSubMesh{D,T},::Val{1},i::Int) where {D,T}
  l = dface_vertices(m,1)[i]
  v = vertex_coordinates(m)
  Segment(v[l[1]],v[l[2]])
end

function get_dface(m::CellSubMesh{D,T},::Val{2},i::Int) where {D,T}
  l = dface_vertices(m,2)[i]
  v = vertex_coordinates(m)
  Triangle(v[l[1]],v[l[2]],v[l[3]])
end

function get_dface(m::CellSubMesh{D,T},::Val{3},i::Int) where {D,T}
  l = dface_vertices(m,3)[i]
  v = vertex_coordinates(m)
  Tetrahedron(v[l[1]],v[l[2]],v[l[3]],v[l[4]])
end

function get_vertex(m::CellSubMesh,i)
  get_dface(m,Val{0}(),i)
end

function get_edge(m::CellSubMesh,i)
  get_dface(m,Val{1}(),i)
end

function get_facet(m::CellSubMesh{D},i) where D
  get_dface(m,Val{D-1}(),i)
end

function get_cell(m::CellSubMesh{D},i) where D
  get_dface(m,Val{D}(),i)
end

function distance(p::Point{2},t::Triangle{2})
  if have_intersection(p,t)
    distance = 0.0
  else
    distance = typemax(0.0)
  end
end

@generated function distance(m::CellSubMesh{D,T},d::Int,i::Int,p::Point{D,T}) where {D,T}
  body = ""
  for d in 0:D
    body *= "if d == $d \n"
    body *= "  f$d = get_dface(m,Val{$d}(),i) \n"
    body *= "  distance(p,f$d) \n"
    body *= "else"
  end
  error = "  throw(ErrorException(\"\$d-face does not exist\"))\n"
  str = body * '\n' * error * "end"
  Meta.parse(str)
end

function dface_nfaces(m::CellSubMesh,d::Int,n::Int)
  m.dface_to_nfaces[d+1][n+1]
end

function add_dface!(m::CellSubMesh,d::Int,new_df_to_nf::Vector{Vector{Vector{Int}}})
  for n in 0:d-1
    nfaces = dface_nfaces(m,d,n)
    new_nfaces = new_df_to_nf[n+1]
    resize!( nfaces, length(nfaces) + length(new_nfaces) )
    for i in 1:length(new_nfaces)
      k = i + length(nfaces) - length(new_nfaces) 
      resize!( nfaces, k, length(new_nfaces[i]) )
      for j in 1:length(new_nfaces[i])
        nfaces[k][j] = new_nfaces[i][j]
      end
    end
  end
  m
end

function add_vertex!(m::CellSubMesh{D,T},p::Point{D,T}) where {D,T}
  push!( vertex_coordinates(m), p )
end

function delete_dface!(m::CellSubMesh,d::Int,i::Int)
  for n in 0:d-1
    deleteat!(m.dface_to_nfaces[d+1][n+1],i)
  end
end

function local_to_global_nfaces!(nfaces_cache::Vector{Vector{Int}},m::CellSubMesh,nf_to_nF::Vector{Vector{Int}},d::Int,id::Int) 

  dface = dface_vertices(m,d)[id]
  resize!(nfaces_cache,d)
  resize!(nfaces_cache,0+1, length(dface) + 1 )
  nfaces_cache[0+1][1:length(dface)] = dface
  nfaces_cache[0+1][end] = num_dfaces(m,0) + 1 
  for n in 1:d-1
    resize!(nfaces_cache,n+1,length(nf_to_nF[n]))
    n_dfaces = num_dfaces(m,n)
    for i in 1:length(nf_to_nF[n])
      if nf_to_nF[n][i] == 0
        n_dfaces += n_dfaces 
        nfaces_cache[n+1][i] = n_dfaces
      else
        nfaces_cache[n+1][i] = dface_nfaces(m,d,n)[id][ nf_to_nF[n][i] ]
      end
    end
  end
  nfaces_cache
end

function refine!(m::CellSubMesh{D,T},d::Int,id::Int,p::Point{D,T}) where {D,T}
  if d == 0
    return m
  else
    @check d == D "smaller cuts not implemented yet"
    case = cut_tet_dface_to_case[D][d][1]
    
    nf_to_mf = cut_tet_nface_to_mface[D][case]
    nf_to_nF = cut_tet_nsubface_to_nface[D][case]
    
    nfaces_cache = Vector{Int}[]

    local_to_global_nfaces!(nfaces_cache, m, nf_to_nF,d,id )

    nf_to_mf_cache = Vector{Vector{Vector{Int}}}[]
    resize!(nf_to_mf_cache,d+1)

    for dim in 1:d
      df_to_nf = nf_to_mf[dim+1]
      resize!(nf_to_mf_cache,dim+1,dim)
      df_to_nf_cache = nf_to_mf_cache[dim+1]
      for n in 0:dim-1
        resize!(df_to_nf_cache,n+1,size(df_to_nf[n+1],2)) 
        for i in 1:size(df_to_nf[n+1],2)
          resize!(df_to_nf_cache[n+1],i,size(df_to_nf[n+1],1))
          for j in 1:size(df_to_nf[n+1],1)
            df_to_nf_cache[n+1][i][j] = nfaces_cache[n+1][ df_to_nf[n+1][j,i] ] 
          end
        end
      end
    end

    delete_dface!(m,d,id)

    add_vertex!(m,point)
    for dim in 1:d
      add_dface!(m,dim,nf_to_mf_cache[dim+1])
    end

    return m
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

p1 = Point(1,2,3)
p2 = Point(4,5,6)

b=BoundingBox(p1,p2)

h = Hexahedron(b)

b2 = BoundingBox(h)

@test b2 == b

m = initialize(h)

p1 = Point(1,2,4)
p2 = Point(5,3,5)

box = BoundingBox(p1,p2)

cell = Hexahedron(box)

initialize!(m,cell)

@test @allocated(initialize!(m,cell)) == 0

p0 = Point(0.0,0.0)
p1 = Point(1.0,1.0)
box = BoundingBox(p0,p1)
cell = Hexahedron(box)

mesh = initialize(cell)

stl_points = [ Point(0.5,0.25), Point(0.25,0.5), Point(0.75,0.25), Point(0.75,0.5)  ]
point = stl_points[1]

for k in 1:4
  global stl_points, point
  point = stl_points[k]
  D = 2
  min_distance = intersection_tolerance 
  idface = 0
  for d in 0:D
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
end

writevtk(mesh,"sub_mesh")
# struct DFaceCutter end


## TODO encapsulate refinement procedures
## encapsulate cache in DFaceCutter
## allow to cut edges and dfaces with d < D
## revise why is not cut properly in the example of more than one point

end # module
