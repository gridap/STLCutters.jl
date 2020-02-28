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
  dface_to_nfaces::Matrix{Table{Int}}
  dface_to_new_dfaces::Vector{Vector{Vector{Int}}}
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
  dface_to_nfaces = Matrix{Table{Int}}(undef,D+1,D+1)

  for d in 0:D, n in 0:d
    df_to_nf = hex_to_tet_nface_to_mface[D][d+1][n+1]
    dface_to_nfaces[d+1,n+1] = Table( df_to_nf )
  end

  nf_to_new_df = [ Vector{Int}[] for d in 1:D ]

  for d in 1:D
    nf_to_new_df[d] = [ Int[] for i in 1:length( dface_to_nfaces[d+1,0+1] ) ]
  end

  CellSubMesh{D,T}(coordinates,dface_to_nfaces,nf_to_new_df)
end

function initialize!(m::CellSubMesh{D,T},c::Hexahedron{N,D,T}) where {N,D,T}
  resize!(m.vertex_coordinates,N)
  for (i,v) in enumerate( get_vertices(c) )
    m.vertex_coordinates[i] = v
  end
  for d in 0:D, n in 0:d
    df_to_nf = hex_to_tet_nface_to_mface[D][d+1][n+1]
    resize!( m.dface_to_nfaces[d+1,n+1], 0 )
    push!( m.dface_to_nfaces[d+1,n+1], df_to_nf )
  end
  for d in 1:D
    resize!( m.dface_to_new_dfaces[d], num_dfaces(m,d) )
    for i in 1:num_dfaces(m,d)
      resize!( m.dface_to_new_dfaces[d][i],0)
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
    length(m.dface_to_nfaces[d+1,0+1])
  end
end

function dface_to_nfaces(m::CellSubMesh,d::Int,n::Int)
  m.dface_to_nfaces[d+1,n+1]
end

function dface_vertices(m::CellSubMesh,d::Int) 
  m.dface_to_nfaces[d+1,0+1]
end

function get_vertex_coordinates(m::CellSubMesh) 
  m.vertex_coordinates
end

function get_dface(m::CellSubMesh,::Val{0},i::Integer)
  v = get_vertex_coordinates(m)
  v[i]
end

function get_dface(m::CellSubMesh,::Val{1},i::Integer)
  df_to_v = dface_vertices(m,1)
  v = get_vertex_coordinates(m)
  Segment( v[ df_to_v[i,1] ], v[ df_to_v[i,2] ] )
end

function get_dface(m::CellSubMesh,::Val{2},i::Integer)
  df_to_v = dface_vertices(m,2)
  v = get_vertex_coordinates(m)
  Triangle( v[ df_to_v[i,1] ], v[ df_to_v[i,2] ], v[ df_to_v[i,3] ] )
end

function get_dface(m::CellSubMesh,::Val{3},i::Integer)
  df_to_v = dface_vertices(m,3)
  v = get_vertex_coordinates(m)
  Tetrahedron( v[ df_to_v[i,1]], v[ df_to_v[i,2] ], v[ df_to_v[i,3] ], v[ df_to_v[i,4] ] )
end

function get_dface(m::CellSubMesh,::Val{d},i::Integer) where d
  throw(ArgumentError("get_dface(::CellSubMesh,::Val($d),::Integer) not implemented"))
end

function STLCutter.isactive(m::CellSubMesh,d::Integer,i::Integer)
  d != 0 || return true
  isactive(dface_vertices(m,d),i)
end

function get_vertex(m::CellSubMesh,i::Integer)
  get_dface(m,Val{0}(),i)
end

function get_edge(m::CellSubMesh,i::Integer)
  get_dface(m,Val{1}(),i)
end

function get_facet(m::CellSubMesh{D},i::Integer) where D
  get_dface(m,Val{D-1}(),i)
end

function get_cell(m::CellSubMesh{D},i::Integer) where D
  get_dface(m,Val{D}(),i)
end

function distance(p::Point{2},t::Triangle{2})
  if have_intersection(p,t)
    distance = 0.0
  else
    distance = typemax(0.0)
  end
end

function distance(p::Point{3},t::Tetrahedron{3})
  if have_intersection(p,t)
    distance = 0.0
  else
    distance = typemax(0.0)
  end
end

@generated function distance(m::CellSubMesh{D},d::Int,i::Int,p::Point{D}) where D
  body = ""
  for d in 0:D
    body *= "if d == $d \n"
    body *= "  f$d = get_dface(m,Val{$d}(),i) \n"
    body *= "  distance(p,f$d) \n"
    body *= "else"
  end
  error = "  throw(ArgumentError(\"\$d-face does not exist\"))\n"
  str = body * '\n' * error * "end"
  Meta.parse(str)
end



function add_vertex!(m::CellSubMesh{D,T},p::Point{D,T}) where {D,T}
  v = get_vertex_coordinates(m)
  push!( v, p )
end

function delete_dface!(m::CellSubMesh,d::Integer,i::Int)
  for n in 0:d-1
    df_to_nf = dface_to_nfaces(m,d,n)
    remove!(df_to_nf,i)
  end
end

function refine!(m::CellSubMesh{D,T},dim::Integer,id::Integer,p::Point{D,T}) where {D,T}
  if dim == 0
    return m
  else
    @check dim == D "smaller cuts not implemented yet"
    case = cut_tet_dface_to_case[D][dim][1]

    tables_nf_to_mf = cut_tet_nface_to_mface[D][case]
    tables_nf_to_nF = cut_tet_nsubface_to_nface[D][case]

    nfaces_cache = zero(Table{Int})
    nf_to_mf_cache = Matrix(undef,D+1,D+1)
    for i in 1:D+1, j in 1:D+1
      nf_to_mf_cache[i,j] = zero(Table{Int})
    end

    resize!(nfaces_cache,0)
    push!(nfaces_cache,tables_nf_to_nF)
    for n in 0:dim
      df_to_nf = dface_to_nfaces(m,dim,n)
      n_dfaces = num_dfaces(m,n)
      for i in 1:length(nfaces_cache,n+1)
        if nfaces_cache[n+1,i] == 0
          n_dfaces += 1
          nfaces_cache[n+1,i] = n_dfaces
        else
          nfaces_cache[n+1,i] = df_to_nf[ id, nfaces_cache[n+1,i] ]
        end
      end
    end

    for d in 0:dim, n in 0:d
      df_to_nf = nf_to_mf_cache[d+1,n+1]
      resize!( df_to_nf, 0 )
      push!( df_to_nf, tables_nf_to_mf[d+1][n+1] )
      for i in 1:length(df_to_nf), j in 1:length(df_to_nf,i)
        df_to_nf[i,j] = nfaces_cache[ n+1, df_to_nf[i,j] ]
      end
    end

    for d in 0:dim-1, n in 0:d
      df_to_nf = nf_to_mf_cache[d+1,n+1]
      for i in 1:length(nfaces_cache,d+1)
        if nfaces_cache[d+1,i] <= num_dfaces(m,d)
          nfaces_cache[d+1,i], num_dfaces(m,d) 
          remove!(df_to_nf,i)
        end
      end
    end
    
    for d in 1:dim
      df_to_new_df = m.dface_to_new_dfaces[d]
      for i in 1:length(nfaces_cache,d+1)
        if nfaces_cache[d+1,i] > length( df_to_new_df )
          resize!( df_to_new_df,nfaces_cache[d+1,i] )
        end
        if !isassigned(df_to_new_df,nfaces_cache[d+1,i])
          df_to_new_df[nfaces_cache[d+1,i]] = Int[]
        elseif nfaces_cache[d+1,i] > num_dfaces(m,d)
          resize!(df_to_new_df[nfaces_cache[d+1,i]],0)
        end
      end
    end
    for i in 1:length(nfaces_cache,dim+1)
      push!( m.dface_to_new_dfaces[dim][id], nfaces_cache[dim+1,i] )
    end

    add_vertex!(m,point)

    for d in 0:dim, n in 0:d
      push!( dface_to_nfaces(m,d,n), nf_to_mf_cache[d+1,n+1] )
    end
    delete_dface!(m,dim,id)
    return m
  end
end

function case_ref(m::CellSubMesh{D},dim::Integer,id::Integer) where D
  for d in dim:D
    df_to_nf = dface_to_nfaces(m,d,dim)
    for i in 1:length(df_to_nf)
      for j in 1:length(df_to_nf,i)
        if df_to_nf[i,j] == id
          #refine
          cut_tet_dface_to_case[D][d][j]
          break
        end
      end
    end
  end
end

function writevtk(m::CellSubMesh{D,T},file_base_name) where {D,T}
  d_to_vtk_type_id = Dict(0=>1,1=>3,2=>5) # tets
  num_points = num_vertices(m)
  points = zeros(T,D,num_points)
  for (i ,p ) in enumerate(get_vertex_coordinates(m)), d in 1:D
    points[d,i] = p[d]
  end
  cells = MeshCell{Vector{Int64}}[]
  for d in 1:D
    dface_to_vertices = dface_vertices(m,d)
    num_dfaces = length(dface_to_vertices)
    vtk_type = VTKCellType(d_to_vtk_type_id[d])
    for i in 1:num_dfaces
      vertices = [ dface_to_vertices[i,j] for j in 1:vtk_type.nodes ]
      push!( cells, MeshCell(vtk_type,vertices) )
    end
  end
  vtkfile = vtk_grid(file_base_name,points,cells)
  vtk_save(vtkfile)
end

function compact(m::CellSubMesh{D}) where D
  for d in 0:D, n in 0:d-1
    compact!( dface_to_nfaces(mesh,d,n) )
  end
  #TODO: Update connectivities with new indexes
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
      if isactive(mesh,d,i)
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
end

compact(mesh)

writevtk(mesh,"sub_mesh")

# struct DFaceCutter end


## TODO encapsulate refinement procedures
## encapsulate cache in DFaceCutter
## allow to cut edges and dfaces with d < D
## revise why is not cut properly in the example of more than one point

end # module
