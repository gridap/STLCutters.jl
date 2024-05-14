
## Constants

const OPEN = -1

"""
    struct GraphPolytope{D,Dp,Tp} <: Polytope{D}

  A graph polytope is a polytope defined by a set of vertices and a rototation
  system (a planar oriented graph). This polytopal representation can represent
  any polytope in 2 and 3 dimensions.

  In 2 dimensions ([`Polygon`](@ref)), the representation of the polygon is a closed polyline.

  In 3 dimensions ([`Polyhedron`](@ref)), the rotation system generates the connectivities, each   facet is a closed cycle of the graph.
  This construction allows complex geometrical operations, e.g., intersecting polytopes by halfspaces.
  See also,

  > K. Sugihara, "A robust and consistent algorithm for intersecting convex polyhedra", Comput. Graph. Forum 13 (3) (1994) 45–54, doi: [10.1111/1467-8659.1330045](https://doi.org/10.1111/1467-8659.1330045)

  > D. Powell, T. Abel, "An exact general remeshing scheme applied to physically conservative voxelization", J. Comput. Phys. 297 (Sept. 2015) 340–356, doi: [10.1016/j.jcp.2015.05.022](https://doi.org/10.1016/j.jcp.2015.05.022.

  > S. Badia, P. A. Martorell, F. Verdugo. "Geometrical discretisations for unfitted finite elements on explicit boundary representations", J.Comput. Phys. 460 (2022): 111162. doi: [10.1016/j.jcp.2022.111162](https://doi.org/10.1016/j.jcp.2022.111162)
"""
struct GraphPolytope{D,Dp,Tp,Td} <: Polytope{D}
  vertices::Vector{Point{Dp,Tp}}
  edge_vertex_graph::Vector{Vector{Int32}}
  n_m_to_nface_to_mfaces::Matrix{Vector{Vector{Int32}}}
  facedims::Vector{Int32}
  dimranges::Vector{UnitRange{Int}}
  dface_nfaces::Vector{Vector{Int32}}
  isopen::Bool
  metadata::Td
  function GraphPolytope{D}(
    vertices::Vector{Point{Dp,Tp}},
    edge_vertex_graph::Vector{Vector{Int32}},
    n_m_to_nface_to_mfaces::Matrix{Vector{Vector{Int32}}},
    facedims::Vector{Int32},
    dimranges::Vector{UnitRange{Int}},
    dface_nfaces::Vector{Vector{Int32}},
    isopen::Bool,
    metadata::Td) where {D,Dp,Tp,Td}
    new{D,Dp,Tp,Td}(
      vertices,
      edge_vertex_graph,
      n_m_to_nface_to_mfaces,
      facedims,
      dimranges,
      dface_nfaces,
      isopen,
      metadata)
  end
end

"""
    Polygon = GraphPolytope{2}

  A polygon is a [`GraphPolytope`](@ref) in 2 dimensions.
"""
const Polygon = GraphPolytope{2}

"""
    Polyhedron = GraphPolytope{3}

  A polyhedron is a [`GraphPolytope`](@ref) in 3 dimensions.
"""
const Polyhedron = GraphPolytope{3}

# Constructors

function GraphPolytope{D}(
  vertices::Vector{<:Point},
  graph::Vector{Vector{Int32}},
  isopen::Bool,
  data) where D

  n = D+1
  n_m_to_nface_to_mfaces = Matrix{Vector{Vector{Int32}}}(undef,n,n )
  dimranges = Vector{UnitRange{Int}}(undef,0)
  dface_nfaces = Vector{Vector{Int32}}(undef,0)
  facedims = Vector{Int32}(undef,0)
  GraphPolytope{D}(
    vertices,
    graph,
    n_m_to_nface_to_mfaces,
    facedims,
    dimranges,
    dface_nfaces,
    isopen,
    data)
end

function GraphPolytope{D}(
  vertices::AbstractVector{<:Point},
  graph::Vector{<:Vector},
  isopen::Bool,
  data) where D

  GraphPolytope{D}(collect(vertices),graph,isopen,data)
end

function GraphPolytope{D}(
  vertices::AbstractVector{<:Point},
  graph::Vector{<:Vector}
  ;isopen=false::Bool,
  metadata=nothing) where D

  GraphPolytope{D}(vertices,graph,isopen,metadata)
end

function GraphPolytope{D}(
  vertices::AbstractVector{<:Point},
  graph::Vector{<:Vector},
  data) where D

  isopen = false
  GraphPolytope{D}(vertices,graph,isopen,data)
end

function GraphPolytope{D}(t::GridTopology{Dc,D};metadata=nothing) where {Dc,D}
  graph = compute_graph(t)
  X = get_vertex_coordinates(t)
  isopen = is_open_surface(t)
  p = GraphPolytope{D}(X,graph;isopen,metadata)
  set_polytope_data!(p,t,metadata)
  p
end

function set_polytope_data!(p::GraphPolytope,t::GridTopology,::Nothing)
  p
end

function get_polytope_data(p::Polytope;metadata=nothing)
  get_polytope_data(p,metadata)
end

function get_polytope_data(p::Polytope,::Nothing)
  nothing
end

function Polygon(p::Polytope{2},vertices::AbstractVector{<:Point};kwargs...)
  if p == TRI
    e_v_graph = [[2,3],[3,1],[1,2]]
    perm = [1,2,3]
  elseif p == QUAD
    e_v_graph = [[2, 3],[4, 1],[1, 4],[3, 2]]
    perm = [1,2,4,3]
  else
    @unreachable
  end
  vertices = map(Reindex(vertices),perm)
  e_v_graph = map(Reindex(e_v_graph),perm)
  e_v_graph = map(i->replace(i, Dict(perm .=> 1:length(perm))...),e_v_graph)
  e_v_graph = map(i->Int32.(i),e_v_graph)
  data = get_polytope_data(p;kwargs...)
  Polygon(vertices,e_v_graph,data)
end

function Polygon(vertices::AbstractVector{<:Point};kwargs...)
  graph = map(1:length(vertices)) do i
    inext = i == length(vertices) ? 1 : i+1
    iprev = i == 1 ? length(vertices) : i-1
    Int32[iprev,inext]
  end
  Polygon(vertices,graph;kwargs...)
end

function Polyhedron(p::Polytope{3},vertices::AbstractVector{<:Point};kwargs...)
  if p == TET
    e_v_graph = [[2,4,3],[3,4,1],[1,4,2],[1,2,3]]
  elseif p == HEX
    e_v_graph = [
      [5, 2, 3],
      [6, 4, 1],
      [7, 1, 4],
      [8, 3, 2],
      [1, 7, 6],
      [2, 5, 8],
      [3, 8, 5],
      [4, 6, 7] ]
  else
    @unreachable
  end
  e_v_graph = map(i->Int32.(i),e_v_graph)
  data = get_polytope_data(p;kwargs...)
  Polyhedron(vertices,e_v_graph,data)
end

function GraphPolytope{D}(p::Polytope;kwargs...) where D
  GraphPolytope{D}(p,get_vertex_coordinates(p);kwargs...)
end

# Interface

num_dims(::T) where T<:GraphPolytope = num_dims(T)

num_dims(::Type{<:GraphPolytope{D}}) where D = D

num_cell_dims(a::GraphPolytope) = num_dims(a)

point_eltype(::T) where T<:GraphPolytope = point_eltype(T)

point_eltype(::Type{<:GraphPolytope{D,Dp,T}}) where {D,Dp,T} = T

num_point_dims(::Type{<:GraphPolytope{D,Dp}}) where {D,Dp} = Dp

num_vertices(a::GraphPolytope) = length(a.vertices)

get_vertex_coordinates(a::GraphPolytope) = a.vertices

Base.getindex(a::GraphPolytope,i::Integer) = a.vertices[i]

"""
    get_graph(p::GraphPolytope) -> Vector{Vector{Int32}}

  It returns the edge-vertex graph of the polytope `p`.
"""
@inline get_graph(a::GraphPolytope) = a.edge_vertex_graph


"""
    get_data(p::GraphPolytope)

  It return the metadata stored in the polytope `p`.
"""
get_data(a::GraphPolytope) = a.metadata


"""
    isopen(p::GraphPolytope) -> Bool

  In return whether the polytope is watter tight or not.
"""
Base.isopen(a::GraphPolytope) = a.isopen

"""
    isactive(p::GraphPolytope,vertex::Integer) -> Bool

  It returns whether a vertex is connected to any other vertex.
"""
function isactive(p::Polyhedron,vertex::Integer)
  !isempty( get_graph(p)[vertex] )
end

"""
    check_graph(p::GraphPolytope) -> Bool

  It checks whether the graph is well-constructed. The graph must be oriented
  and planar.
"""
function check_graph(p::GraphPolytope)
  check_graph(get_graph(p))
end

function check_graph(graph::AbstractVector{<:AbstractVector})
  for v in 1:length(graph)
    !isempty(graph[v]) || continue
    for vneig in graph[v]
      vneig ∉ (OPEN,UNSET) || continue
      v ∈ graph[vneig] || return false
    end
  end
  true
end

is_simplex(::GraphPolytope) = false

is_n_cube(::GraphPolytope) = false

function simplexify(p::GraphPolytope{D}) where D
  @assert !isopen(p)
  X,T = simplexify_interior(p)
  @check X == get_vertex_coordinates(p)
  T, simplex_polytope(Val{D}())
end

simplex_polytope(::Val{0}) = VERTEX

simplex_polytope(::Val{1}) = SEGMENT

simplex_polytope(::Val{2}) = TRI

simplex_polytope(::Val{3}) = TET

function Polytope{0}(p::GraphPolytope,faceid::Integer)
  VERTEX
end

function Polytope{1}(p::GraphPolytope,faceid::Integer)
  SEGMENT
end

function Polytope{D}(p::GraphPolytope{D},faceid::Integer) where D
  p
end

function Polytope{2}(p::Polyhedron,faceid::Integer)
  f_to_v = get_faces(p,2,0)
  coords = get_vertex_coordinates(p)
  vertices = coords[f_to_v[faceid]]
  Polygon(vertices)
end

Base.:(==)(a::GraphPolytope,b::GraphPolytope) = false

function num_faces(p::GraphPolytope{D},d::Integer) where D
  if d == 0
    num_vertices(p)
  elseif d == D
    1
  else
    length( get_faces(p,d,0) )
  end
end

function get_faces(p::GraphPolytope,n::Integer,m::Integer)
  setup_faces!(p,n,m)
  p.n_m_to_nface_to_mfaces[n+1,m+1]
end


function get_facet_orientations(p::GraphPolytope)
  Ones(num_facets(p))
end

function get_facet_normal(p::Polyhedron)
  D = 3
  f_to_v = get_faces(p,D-1,0)
  coords = get_vertex_coordinates(p)
  map(f_to_v) do v
    v1 = coords[v[2]]-coords[v[1]]
    v2 = coords[v[3]]-coords[v[1]]
    v1 /= norm(v1)
    v2 /= norm(v2)
    n = v1 × v2
    n /= norm(n)
  end
end

function get_facet_normal(p::Polygon)
  D = 2
  f_to_v = get_faces(p,D-1,0)
  coords = get_vertex_coordinates(p)
  @notimplementedif num_components(eltype(coords)) != 2
  map(f_to_v) do v
    e = coords[v[2]]-coords[v[1]]
    n = VectorValue( e[2], -e[1] )
    n /= norm(n)
  end
end

function get_edge_tangent(p::GraphPolytope)
  e_to_v = get_faces(p,1,0)
  coords = get_vertex_coordinates(p)
  map(e_to_v) do v
    e = coords[v[2]]-coords[v[1]]
    e / norm(e)
  end
end

function get_dimranges(p::GraphPolytope)
  setup_dimranges!(p)
  p.dimranges
end

function get_dimrange(p::GraphPolytope,d::Integer)
  setup_dimranges!(p)
  p.dimranges[d+1]
end

function get_faces(p::GraphPolytope)
  setup_face_to_nfaces!(p)
  p.dface_nfaces
end

function get_facedims(p::GraphPolytope)
  setup_facedims!(p)
  p.facedims
end

function setup_dimranges!(p::GraphPolytope{D}) where D
  if length(p.dimranges) < D+1
    lens = map(i->num_faces(p,i),0:D)
    sums = cumsum(lens)
    resize!(p.dimranges,D+1)
    for (i,(l,s)) in enumerate(zip(lens,sums))
      p.dimranges[i] = s-l+1 : s
    end
  end
end

function setup_face_to_nfaces!(p::GraphPolytope)
  if length(p.dface_nfaces) < num_vertices(p)
    facedims = get_facedims(p)
    dface_nfaces = Vector{Vector{Int32}}(undef,length(facedims))
    ofsets = get_offsets(p)
    for (f,d) in enumerate(facedims)
      df = f - ofsets[d+1]
      nfs = Int[]
      for n in 0:d
        union!(nfs,get_faces(p,d,n)[df] .+ ofsets[n+1])
      end
      dface_nfaces[f] = nfs
    end
    copy!(p.dface_nfaces,dface_nfaces)
  end
  nothing
end

function setup_facedims!(p::GraphPolytope)
  if length(p.facedims) < num_vertices(p)
    dimranges = get_dimranges(p)
    n_faces = dimranges[end][end]
    facedims = _nface_to_nfacedim(n_faces,dimranges)
    copy!(p.facedims,facedims)
  end
end

function setup_faces!(p::GraphPolytope{D},dimfrom,dimto) where D
  if isassigned(p.n_m_to_nface_to_mfaces,dimfrom+1,dimto+1)
    return nothing
  end
  if dimfrom == dimto
    setup_nface_to_nface!(p,dimfrom)
  elseif dimfrom == D
    setup_cell_to_faces!(p,dimto)
  elseif dimto == 0
    setup_face_to_vertices!(p,dimfrom)
  elseif dimfrom > dimto
    setup_nface_to_mface!(p,dimfrom,dimto)
  else
    setup_nface_to_mface_dual!(p,dimto,dimfrom)
  end
  nothing
end

function setup_face_to_vertices!(p::GraphPolytope,d)
  if !isassigned(p.n_m_to_nface_to_mfaces,d+1,1)
    df_to_v = generate_face_to_vertices(p,d)
    p.n_m_to_nface_to_mfaces[d+1,1] = df_to_v
  end
end

function setup_cell_to_faces!(p::GraphPolytope{D},d) where D
  if !isassigned(p.n_m_to_nface_to_mfaces,D+1,d+1)
    num_f = num_faces(p,d)
    c_to_f = [ collect(1:num_f) ]
    p.n_m_to_nface_to_mfaces[D+1,d+1] = c_to_f
  end
end

function setup_nface_to_nface!(p::GraphPolytope,n)
  if !isassigned(p.n_m_to_nface_to_mfaces,n+1,n+1)
    num_nf = num_faces(p,n)
    nf_to_nf = map(i->Int[i],1:num_nf)
    p.n_m_to_nface_to_mfaces[n+1,n+1] = nf_to_nf
  end
end

function setup_nface_to_mface!(p::Polyhedron,n,m)
  if !isassigned(p.n_m_to_nface_to_mfaces,n+1,m+1)
    @notimplementedif n != 2 && m != 1
    nf_to_v = get_faces(p,n,0)
    v_to_mf = get_faces(p,0,m)
    nf_to_ftype = map( length, nf_to_v )
    ftype_to_lmf_to_lv = map(1:maximum(nf_to_ftype)) do ftype
      map(1:ftype) do i
        inext = i == ftype ? 1 : i+1
        Int[i,inext]
      end
    end
    nf_to_mf = find_cell_to_faces(
      Table(nf_to_v),
      ftype_to_lmf_to_lv,
      nf_to_ftype,
      Table(v_to_mf))
    p.n_m_to_nface_to_mfaces[n+1,m+1] = Vector(nf_to_mf)
  end
end

function setup_nface_to_mface_dual!(p::GraphPolytope,dimto,dimfrom)
  if !isassigned(p.n_m_to_nface_to_mfaces,dimfrom+1,dimto+1)
    @assert dimfrom < dimto
    nf_to_mf = get_faces(p,dimto,dimfrom)
    n_mf = num_faces(p,dimfrom)
    mf_to_nf = generate_cells_around(Table(nf_to_mf),n_mf)
    p.n_m_to_nface_to_mfaces[dimfrom+1,dimto+1] = Vector(mf_to_nf)
  end
end

function generate_face_to_vertices(p::Polyhedron,d::Integer)
  if d == 1
    generate_edge_to_vertices(p)
  elseif d == 2
    generate_facet_to_vertices(p)
  else
    @unreachable
  end
end

function generate_face_to_vertices(p::Polygon,d::Integer)
  if d == 1
    generate_facet_to_vertices(p)
  else
    @unreachable
  end
end

function generate_facet_to_vertices(poly::Polyhedron)
  D = 3
  istouch = map( i -> falses(length(i)), get_graph(poly) )
  T = Vector{Int32}[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    for i in 1:length(get_graph(poly)[v])
      !istouch[v][i] || continue
      istouch[v][i] = true
      vcurrent = v
      vnext = get_graph(poly)[v][i]
      vnext ∉ (OPEN,UNSET) || continue
      k = [v]
      while vnext != v
        inext = findfirst( isequal(vcurrent), get_graph(poly)[vnext] )
        inext = ( inext % length( get_graph(poly)[vnext] ) ) + 1
        istouch[vnext][inext] = true
        vcurrent = vnext
        vnext = get_graph(poly)[vnext][inext]
        vnext ∉ (OPEN,UNSET) || break
        push!(k,vcurrent)
      end
      if length(k) >=D
        push!(T,k)
      end
    end
  end
  T
end

function generate_edge_to_vertices(poly::GraphPolytope)
  graph = get_graph(poly)
  T = Vector{Int32}[]
  for v in 1:length(graph)
    for vneig in graph[v]
      if vneig > v
        push!(T,[v,vneig])
      end
    end
  end
  T
end

function generate_facet_to_vertices(poly::Polygon)
  graph = get_graph(poly)
  T = Vector{Int32}[]
  for v in 1:length(graph)
    vnext = v == length(graph) ? 1 : v+1
    @assert vnext ∈ graph[v]
    push!(T,[v,vnext])
  end
  T
end

function Base.copy(poly::Polyhedron)
  vertices = get_vertex_coordinates(poly)
  graph = get_graph(poly)
  open = isopen(poly)
  data = get_data(poly)
  data = !isnothing(data) ? copy(data) : nothing
  Polyhedron(vertices,graph,open,data)
end

function simplexify_interior(p::Polygon)
  @assert !isopen(p)
  e_to_v = generate_facet_to_vertices(p)
  T = Vector{Int32}[]
  if length(e_to_v) > 0
    v0 = e_to_v[1][1]
    for verts in e_to_v
      if v0 ∉ verts
        push!(T,[v0,verts[1],verts[2]])
      end
    end
  end
  get_vertex_coordinates(p),T
end

function simplexify_interior(poly::Polyhedron)
  !isopen(poly) || return simplexify_surface(poly)
  vstart = fill(UNSET,num_vertices(poly))
  stack = Int32[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    vstart[v] == UNSET || continue
    vstart[v] = v
    empty!(stack)
    push!(stack,v)
    while !isempty(stack)
      vcurrent = pop!(stack)
      for vneig in get_graph(poly)[vcurrent]
        if vstart[vneig] == UNSET
          vstart[vneig] = v
          push!(stack,vneig)
        end
      end
    end
  end
  istouch = map( i -> falses(length(i)), get_graph(poly) )
  vertex_coordinates = get_vertex_coordinates(poly)
  T = Vector{Int32}[]
  X = vertex_coordinates
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    for i in 1:length(get_graph(poly)[v])
      !istouch[v][i] || continue
      istouch[v][i] = true
      vcurrent = v
      vnext = get_graph(poly)[v][i]
      while vnext != v
        inext = findfirst( isequal(vcurrent), get_graph(poly)[vnext] )
        !isnothing(inext) || break
        inext = ( inext % length( get_graph(poly)[vnext] ) ) + 1
        istouch[vnext][inext] = true
        vcurrent = vnext
        vnext = get_graph(poly)[vnext][inext]
        @assert vcurrent ≠ vnext
        vcurrent ≠ vnext || break
        if v ∉ (vstart[v],vcurrent,vnext)
          k = [vstart[v],v,vcurrent,vnext]
          push!(T,k)
        end
      end
    end
  end
  X,T
end

function simplexify_surface(poly::Polyhedron)
  istouch = map( i -> falses(length(i)), get_graph(poly) )
  T = Vector{Int32}[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    for i in 1:length(get_graph(poly)[v])
      !istouch[v][i] || continue
      istouch[v][i] = true
      vcurrent = v
      vnext = get_graph(poly)[v][i]
      vnext ∉ (OPEN,UNSET) || continue
      while vnext != v
        inext = findfirst( isequal(vcurrent), get_graph(poly)[vnext] )
        inext = ( inext % length( get_graph(poly)[vnext] ) ) + 1
        istouch[vnext][inext] = true
        vcurrent = vnext
        vnext = get_graph(poly)[vnext][inext]
        vnext ∉ (OPEN,UNSET) || break
        if v ∉ (vcurrent,vnext)
          k = [v,vcurrent,vnext]
          push!(T,k)
        end
      end
    end
  end
  get_vertex_coordinates(poly),T
end

"""
    PolyhedronData

  Metadata for [`GraphPolytope`](@ref) that serves for performing geometrical
  operations.

  The metadata stores the following information:
    * `vertex_to_planes`: A list of planes that intersect each vertex.
    * `vertex_to_original_faces`: A list of d-faces of the original polytope or STL
    * `vertex_to_parent_vertex`: It maps a vertex to a vertex that has the same coordinates
    * `vertex_to_parent_edge`: It maps a vertex a vertex pair that generated that vertex
    * `plane_to_vertex_to_distances`: A list of distances from each plane to each vertex
    * `plane_to_ref_plane`: It maps a plane to a plane that (almost) co-planar
    * `plane_to_ids`: It maps the plane to the plane id
"""
struct PolyhedronData
  vertex_to_planes::Vector{Vector{Int32}}
  vertex_to_original_faces::Vector{Vector{Int32}}
  vertex_to_parent_vertex::Vector{Int32}
  vertex_to_parent_edge::Vector{Tuple{Int32,Int32}}
  plane_to_vertex_to_distances::Vector{Vector{Float64}}
  plane_to_ref_plane::Vector{Int32}
  plane_to_ids::Vector{Int32}
end

struct ClipPolytopeData end

const clipping = ClipPolytopeData()

function GraphPolytope{D}(
    vertices::Vector{<:Point},
    graph::Vector{Vector{Int32}},
    isopen::Bool,
    ::ClipPolytopeData) where D

  data = polyhedron_data(length(vertices))
  GraphPolytope{D}(vertices,graph,isopen,data)
end

function set_polytope_data!(p::GraphPolytope,t::GridTopology,::ClipPolytopeData)
  set_original_faces!(p,t)
  p
end

function get_polytope_data(p::Polytope,::ClipPolytopeData)
  polyhedron_data(p)
end

# Operations

"""
    clip(p::Polyhedron,planes;kwargs...) -> Polyhedron

  It clips a polyhedron by the halfspace of a plane or a set of planes.

  # Optional keyword arguments
  * `inside::Bool=true`: It clips the polyhedron by the inside ourside of the union of halfspace.
  * `inout::Vector{Bool}=trues(length(planes))`: In reverses the halfspaces with `inside[i]=false`.
  * `boundary::nothing`:
      * If `boundary=true`, it preserves the vertices on the planes (zero distance)
      * If `boundary=false`, it removes the vertices on the planes (zero distance)
      * If `boundary=nothing`, it preserves the vertices on the planes if the normal of the plane points inwards.
"""
function clip(poly::Polyhedron,Π;inside=true,inout=trues(length(Π)),boundary=nothing)
  p = poly
  for (i,Πi) in enumerate(Π)
    side = inside == inout[i] ? :left : :right
    p = clip(p,Πi,side;boundary)
    !isnothing(p) || break
  end
  p
end

function split_overlapping(p::Polyhedron,Π)
  p⁻ = clip(p,Π,:left,boundary=true)
  p⁺ = clip(p,Π,:right,boundary=true)
  p⁻,p⁺
end

function split_gapping(p::Polyhedron,Π)
  p⁻ = clip(p,Π,:left,boundary=false)
  p⁺ = clip(p,Π,:right,boundary=false)
  p⁻,p⁺
end

function clip(p::Polyhedron,Π,side;boundary=nothing::Union{Nothing,Bool},invert=false)
  @assert side ∈ (:right,:left)
  if boundary === nothing
    p⁻,p⁺ = split(p,Π,side;invert)
    p  = side == :left ? p⁻ : p⁺
  else
    invert = false
    if side == :right && !boundary
      invert = true
      side = :left
    elseif side == :left && boundary
      invert = true
      side = :right
    end
    p = clip(p,Π,side;invert)
  end
  p
end


"""
    split(p::Polyhedron,plane;kwargs...)

  It splits a polyhedron by a plane into two polyhedra.

  It returns a tuple of `Union{Polyhedron,Nothing}`.
  If one side is empty, it returns `nothing` for that side.

  # Optional keyword arguments
  * `side::Symbol=:both`: It returns `:both` sides, the `:left` side, or the `:right` side
  * `invert::Bool=false`: It inverts the plane
"""
function split(p::Polyhedron,Π,side=:both;invert=false)
  @assert side ∈ (:both,:left,:right)
  distances = get_plane_distances(get_data(p),Π)
  _sign = invert ? (-) : (+)
  distances = lazy_map(_sign,distances)
  smin = Inf
  smax = -Inf
  for v in 1:num_vertices(p)
    isactive(p,v) || continue
    smin = min(smin,distances[v])
    smax = max(smax,distances[v])
  end
  if smin ≥ 0
    return nothing,p
  end
  if smax < 0
    if has_coplanars(get_data(p),Π)
      add_plane!(get_data(p),Π)
    end
    return p,nothing
  end
  new_vertices = empty( get_vertex_coordinates(p) )
  in_graph = _copy(get_graph(p))
  if side == :both
    out_graph = _copy(get_graph(p))
  end
  ≶ = side ∈ (:both,:left) ? (≥) : (<)
  data = copy( get_data(p) )
  D = num_dims(p)
  for v in 1:num_vertices(p)
    isactive(p,v) || continue
    distances[v] ≶ 0 && continue
    for (i,vneig) in enumerate( get_graph(p)[v] )
      vneig ∉ (UNSET,OPEN) || continue
      distances[vneig] ≶ 0 || continue
      vertex = compute_intersection(p[v],distances[v],p[vneig],distances[vneig])
      push!( new_vertices, vertex )
      add_vertex!(data,v,vneig,Π)
      push!( in_graph, fill(UNSET,D) )
      in_graph[v][i] = num_vertices(p) + length(new_vertices)
      in_graph[end][1] = v
      if side == :both
        push!( out_graph, fill(UNSET,D) )
        ineig = findfirst( isequal(v), get_graph(p)[vneig] )
        out_graph[vneig][ineig] = num_vertices(p) + length(new_vertices)
        out_graph[end][1] = vneig
      end
    end
  end
  if side == :both
    p_out = split_postprocess!(out_graph,copy(data),new_vertices,p,distances,(<))
  else
    p_out = nothing
  end
  p_in = split_postprocess!(in_graph,data,new_vertices,p,distances,(≶))
  p⁻ = side ≠ :right ? p_in : p_out
  p⁺ = side ≠ :right ? p_out : p_in
  p⁻,p⁺
end

"""
    simplexify_boundary(p::Polyhedron,t::GridTopology)

  It generates a simplex mesh of the surface of the polyhedron that touches the
  grid topology.

  It returns the coordinates of the vertices, the connectivities and a mapping
  to the topology facets.
"""
function simplexify_boundary(poly::Polyhedron{3},stl::GridTopology)
  D = 3
  stack = Int32[]
  v_to_pv = get_data(poly).vertex_to_parent_vertex
  if isopen(poly)
    v_to_Π = get_data(poly).vertex_to_original_faces
  else
    v_to_Π = get_data(poly).vertex_to_planes
  end
  facedims = get_facedims(stl)
  facet_list = Int32[]
  istouch = map( i -> falses(length(i)), get_graph(poly) )
  vertex_coordinates = get_vertex_coordinates(poly)
  T = Vector{Int32}[]
  X = empty(vertex_coordinates)
  f_to_stlf = Int32[]
  for v in 1:num_vertices(poly)
    isactive(poly,v) || continue
    for i in 1:length(get_graph(poly)[v])
      !istouch[v][i] || continue
      istouch[v][i] = true
      vcurrent = v
      vnext = get_graph(poly)[v][i]
      vnext ∉ (OPEN,UNSET) || continue
      empty!(facet_list)
      append!( facet_list, v_to_Π[v] )
      filter!( i -> i > 0, facet_list)
      filter!( i -> facedims[i] == num_dims(stl), facet_list)
      filter!( i -> i ∈ v_to_Π[vnext], facet_list )
      while vnext != v
        inext = findfirst( isequal(vcurrent), get_graph(poly)[vnext] )
        inext = ( inext % length( get_graph(poly)[vnext] ) ) + 1
        istouch[vnext][inext] = true
        vcurrent = vnext
        vnext = get_graph(poly)[vnext][inext]
        vnext ∉ (OPEN,UNSET) || break
        filter!( i -> i ∈ v_to_Π[vnext], facet_list )
        !isempty(facet_list) || break
        if v ∉ (vcurrent,vnext)
          kx = (v,vcurrent,vnext)
          k = (1:D) .+ length(X)
          for x in kx
            push!(X,vertex_coordinates[x])
          end
          push!(T,k)
          push!(f_to_stlf,first(facet_list))
        end
      end
    end
  end
  X,T,f_to_stlf
end

function simplexify_boundary(::Nothing,args...)
  nothing,nothing,nothing
end

function simplexify_interior(polys::AbstractVector{<:Polyhedron{Dp,Tp}}) where {Dp,Tp}
  T = Vector{Int32}[]
  X = Point{Dp,Tp}[]
  for poly in polys
    Xi,Ti = simplexify_interior(poly)
    append!(T, map(i->i.+length(X),Ti) )
    append!(X,Xi)
  end
  X,T
end

function simplexify_boundary(
  polys::AbstractVector{<:Polyhedron{Dp,Tp}},
  stl::GridTopology) where {Dp,Tp}

  T = Vector{Int32}[]
  X = Point{Dp,Tp}[]
  f_to_stlf = Int32[]
  for poly in polys
    Xi,Ti,f_to_f_i = simplexify_boundary(poly,stl)
    append!(T, map(i->i.+length(X),Ti) )
    append!(X,Xi)
    append!(f_to_stlf,f_to_f_i)
  end
  X,T,f_to_stlf
end

function simplexify_boundary(::NTuple{N,Nothing},args...) where N
  nothing,nothing,nothing,nothing
end

function simplexify_boundary(p::Tuple,labels::Tuple,stl::GridTopology)
  P = eltype(p[1])
  Dp = num_dims(P)
  Tp = point_eltype(P)
  T = Vector{Int32}[]
  X = Point{Dp,Tp}[]
  f_to_stlf = Int32[]
  f_to_ios = Int8[]
  for (p_i,ios) in zip( p, labels )
    Xi, Ti, f_to_f_i = simplexify_boundary(p_i,stl)
    if !isnothing(Ti)
      append!(T, map(i->i.+length(X),Ti) )
      append!(X,Xi)
      append!(f_to_stlf,f_to_f_i)
      append!(f_to_ios, fill(Int8(ios),length(Ti)) )
    end
  end
  X,T,f_to_stlf,f_to_ios
end

function simplexify_cell_boundary(P::Polyhedron,p::Polytope)
  simplexify_cell_boundary(P,num_facets(p))
end

function simplexify_cell_boundary(p::Polyhedron,nf::Integer)
  X,T = simplexify_surface(p)
  v_to_planes = p.metadata.vertex_to_planes
  T_to_planes = lazy_map(Broadcasting(Reindex(v_to_planes)),T)
  c_to_planes = map(i->reduce(intersect,i),T_to_planes)
  c_to_plane = map(first,c_to_planes)
  c_to_facet = map(abs,c_to_plane)
  mask_one = lazy_map(==(1)∘length,c_to_planes)
  mask_bg = lazy_map( i-> -nf ≤ i < 0, c_to_plane )
  mask = map(&,mask_one,mask_bg)
  X,T[mask],c_to_facet[mask]
end

function simplexify_cell_boundary(
  polys::AbstractVector{<:Polyhedron{Dp,Tp}},
  args...) where {Dp,Tp}

  T = Vector{Int32}[]
  X = Point{Dp,Tp}[]
  c_to_bgf = Int32[]
  for poly in polys
    Xi,Ti,c_to_bgf_i = simplexify_cell_boundary(poly,args...)
    append!(T, map(i->i.+length(X),Ti) )
    append!(X,Xi)
    append!(c_to_bgf,c_to_bgf_i)
  end
  X,T,c_to_bgf
end

function surface(poly::Polyhedron{3})
  X,T = simplexify_surface(poly)
  p = TRI
  c = array_cache(T)
  !isempty(T) || return 0.0
  sum( i -> measure(get_cell!(c,X,T,p,i)), 1:length(T) )
end

function surface(poly::Polyhedron{3},stl::GridTopology)
  X,T = simplexify_boundary(poly,stl)
  p = TRI
  c = array_cache(T)
  !isempty(T) || return 0.0
  sum( i -> measure(get_cell!(c,X,T,p,i)), 1:length(T) )
end

function surface(polys::AbstractVector{<:Polyhedron},args...)
  s = 0.0
  for p in polys
    s += surface(p,args...)
  end
  s
end

function volume(poly::Polyhedron{3})
  X,T = simplexify_interior(poly)
  p = TET
  c = array_cache(T)
  sum( i -> measure(get_cell!(c,X,T,p,i)), 1:length(T) )
end

function volume(polys::AbstractVector{<:Polyhedron},args...)
  v = 0.0
  for p in polys
    v += volume(p,args...)
  end
  v
end

function restrict(poly::Polyhedron,stl::GridTopology,stl_facets)
  cell_to_nodes = get_cell_vertices(stl)
  c = array_cache(cell_to_nodes)
  nodes = Int32[]
  for f in stl_facets
    for n in getindex!(c,cell_to_nodes,f)
      if n ∉ nodes
        push!(nodes,n)
      end
    end
  end
  sort!(nodes)
  restrict(poly,nodes)
end

function restrict(p::Polyhedron,nodes)
  graph = get_graph(p)[nodes]
  f = i -> Int32( i ∈ nodes ? findfirst(isequal(i),nodes) : OPEN )
  graph = map(i->map(f,i),graph)
  data = restrict(get_data(p),nodes)
  vertices = get_vertex_coordinates(p)[nodes]
  isopen = true
  Polyhedron(vertices,graph,isopen,data)
end

function restrict(data::PolyhedronData,nodes)
  v_to_Π = data.vertex_to_planes[nodes]
  v_to_of = data.vertex_to_original_faces[nodes]
  v_to_v = collect(1:length(nodes))
  v_to_e = fill((UNSET,UNSET),length(nodes))
  Π_to_rΠ = copy(data.plane_to_ref_plane)
  Π_to_id = copy(data.plane_to_ids)
  _Π_to_v_to_d = data.plane_to_vertex_to_distances
  Π_to_v_to_d = [ _Π_to_v_to_d[i][nodes] for i in 1:length(Π_to_id) ]
  args = v_to_Π, v_to_of, v_to_v, v_to_e, Π_to_v_to_d, Π_to_rΠ, Π_to_id
  PolyhedronData(args...)
end

# Printers

function writevtk(p::Polyhedron,filename;kwargs...)
  writevtk(edge_mesh(p),filename;kwargs...)
end

function writevtk(ps::Array{<:Polyhedron},filename)
  for (i,p) in enumerate(ps)
    writevtk(p,filename*"$i")
  end
end

# Getters

get_vertex_to_planes(a::PolyhedronData) = a.vertex_to_planes

function get_original_reflex_faces(p::Polyhedron{D},stl::GridTopology;empty=true) where D
  get_original_faces(p,stl,Val{D-2}();empty)
end

function get_original_facets(p::Polyhedron{D},stl::GridTopology;empty=false) where D
  get_original_faces(p,stl,Val{D-1}();empty)
end

function get_original_faces(p::Polyhedron,stl::GridTopology,::Val{d};empty) where d
  faces = collect_original_faces(p,stl,d)
  filter!(f->has_original_face(p,f,Val{d}();empty),faces)
  faces
end

function collect_original_faces(p::Polyhedron,stl::GridTopology,d::Integer)
  v_to_f = get_data(p).vertex_to_original_faces
  facedims = get_facedims(stl)
  faces = Int32[]
  for v in 1:num_vertices(p)
    if isactive(p,v)
      for f ∈ v_to_f[v]
        if facedims[f] == d && f ∉ faces
          push!(faces,f)
        end
      end
    end
  end
  faces
end

function has_original_reflex_face(p::Polyhedron{D},face::Integer;empty=true) where D
  has_original_face(p,face,Val{D-2}();empty)
end

function has_original_facet(p::Polyhedron{D},face::Integer;empty=false) where D
  has_original_face(p,face,Val{D-1}();empty)
end

function has_original_face(p::Polyhedron,face::Integer,::Val{0};empty=false)
  v_to_f = get_data(p).vertex_to_original_faces
  for v in 1:num_vertices(p)
    if isactive(p,v) && face ∈ v_to_f[v]
      return true
    end
  end
  false
end

function has_original_face(p::Polyhedron,face::Integer,::Val{1};empty=true)
  v_to_f = get_data(p).vertex_to_original_faces
  v_to_v = get_data(p).vertex_to_parent_vertex
  for v in 1:num_vertices(p)
    if isactive(p,v) && face ∈ v_to_f[v]
      for vneig in get_graph(p)[v]
        vneig ∉ (OPEN,UNSET) || continue
        if face ∈ v_to_f[vneig] && ( empty || v_to_v[v] ≠ v_to_v[vneig] )
          return true
        end
      end
    end
  end
  false
end

function has_original_face(p::Polyhedron,face::Integer,::Val{2};empty=false)
  d = 2
  v_to_f = get_data(p).vertex_to_original_faces
  v_to_v = get_data(p).vertex_to_parent_vertex
  for v in 1:num_vertices(p)
    if isactive(p,v) && face ∈ v_to_f[v]
      for vneig in get_graph(p)[v]
        vneig ∉ (OPEN,UNSET) || continue
        if face ∈ v_to_f[vneig]
          num_v = 1
          vcurrent = v
          vnext = vneig
          while vnext ≠ v
            if v_to_v[vnext] ∉ (v_to_v[v],v_to_v[vcurrent]) || empty
              num_v += 1
            end
            vcurrent,vnext = vnext,next_vertex(p,vcurrent,vnext)
            vnext ∉ (UNSET,OPEN) || break
            face ∈ v_to_f[vnext] || break
          end
          if vnext == v && num_v ≥ d+1
            return true
          end
        end
      end
    end
  end
  false
end

function has_faces(p::Polyhedron,::Val{0})
  for v in 1:num_vertices(p)
    isactive(p,v) || continue
    return true
  end
  false
end

function has_faces(p::Polyhedron,::Val{1})
  for v in 1:num_vertices(p)
    isactive(p,v) || continue
    for vneig in get_graph(p)[v]
      @assert isactive(p,vneig)
      vneig ∉ (OPEN,UNSET) || continue
      v_to_v[v] ≠ v_to_v[vneig] || continue
      return true
    end
    @unreachable
  end
  false
end

function has_faces(p::Polyhedron,::Val{2})
  d = 2
  v_to_v = get_data(p).vertex_to_parent_vertex
  for v in 1:num_vertices(p)
    isactive(p,v) || continue
    for vneig in get_graph(p)[v]
      vneig ∉ (OPEN,UNSET) || continue
      num_v = 1
      vcurrent = v
      vnext = vneig
      while vnext ≠ v
        if v_to_v[vnext] ∉ (v_to_v[v],v_to_v[vcurrent])
          num_v += 1
        end
        vcurrent,vnext = vnext,next_vertex(p,vcurrent,vnext)
        vnext ∉ (UNSET,OPEN) || break
        if vnext == v && num_v ≥ d+1
          return true
        end
      end
    end
  end
  false
end

has_vertices(p::Polyhedron) = has_faces(p,Val{0}())

has_edges(p::Polyhedron) = has_faces(p,Val{1}())

has_facets(p::Polyhedron{D}) where D = has_faces(p,Val{D-1}())

function get_plane_distances(data::PolyhedronData)
  data.plane_to_vertex_to_distances
end

function get_plane_distances(data::PolyhedronData,Π)
  i = findfirst(isequal(Π),get_plane_ids(data))
  get_plane_distances(data)[i]
end

function get_plane_ids(data::PolyhedronData)
  data.plane_to_ids
end

function has_plane(data::PolyhedronData,Π)
  Π ∈ get_plane_ids(data)
end

## Helpers

function compact!(p::Polyhedron)
  ids = findall(i->!isactive(p,i),1:num_vertices(p))
  old_to_new = fill(UNSET,num_vertices(p))
  new = 0
  for i in 1:num_vertices(p)
    isactive(p,i) || continue
    new += 1
    old = i
    old_to_new[old] = new
  end
  vertices = get_vertex_coordinates(p)
  graph = get_graph(p)
  deleteat!(vertices,ids)
  deleteat!(graph,ids)
  f = i -> i ∈ (OPEN,UNSET) ? i : old_to_new[i]
  map!(i->map!(f,i,i),graph,graph)
  compact!(p.metadata,ids,old_to_new)
  p
end

function compact!(data::PolyhedronData,ids,old_to_new)
  v_to_Π = data.vertex_to_planes
  v_to_of = data.vertex_to_original_faces
  v_to_pv = data.vertex_to_parent_vertex
  v_to_pe = data.vertex_to_parent_edge
  Π_to_v_to_dist = data.plane_to_vertex_to_distances
  Π_to_ref_Π = data.plane_to_ref_plane
  Π_to_id = data.plane_to_ids

  deleteat!(v_to_Π,ids)
  deleteat!(v_to_of,ids)
  deleteat!(v_to_pv,ids)
  for (i,pv) in enumerate(v_to_pv)
    if old_to_new[pv] == UNSET
      old_to_new[pv] = i
    end
    v_to_pv[i] = old_to_new[pv]
  end
  deleteat!(v_to_pe,ids)
  f = i -> i == UNSET ? i : old_to_new[i]
  for (i,pe) in enumerate(v_to_pe)
    v_to_pe[i] = (f(pe[1]),f(pe[2]))
  end
  for v_to_dist in Π_to_v_to_dist
    deleteat!(v_to_dist,ids)
  end
  data
end

function Base.copy(data::PolyhedronData)
  v_to_Π = copy(data.vertex_to_planes)
  v_to_of = copy(data.vertex_to_original_faces)
  v_to_pv = copy(data.vertex_to_parent_vertex)
  v_to_pe = copy(data.vertex_to_parent_edge)
  Π_to_v_to_dist = _copy(data.plane_to_vertex_to_distances)
  Π_to_ref_Π = copy(data.plane_to_ref_plane)
  Π_to_id = copy(data.plane_to_ids)
  args = v_to_Π, v_to_of, v_to_pv, v_to_pe, Π_to_v_to_dist, Π_to_ref_Π, Π_to_id
  PolyhedronData(args...)
end

"""
    compute_distances!(p::Polyhedron,planes,plane_ids[;atol=0])

  Compute the distances from the vertices of the polyhedron to each plane of
  the list `planes`. If the distance is below `atol`, it sets the distance to
  zero.
"""
function compute_distances!(p::Polyhedron,Π,faces;atol=0)
  data = get_data(p)
  v_to_f = get_data(p).vertex_to_original_faces
  Π_to_v_to_d = get_plane_distances(data)
  Π_ids = get_plane_ids(data)
  Π_to_rΠ = data.plane_to_ref_plane
  num_Π = length(Π_ids)
  append!(Π_ids,faces)
  append!(Π_to_rΠ,fill(UNSET,length(faces)))
  append!(Π_to_v_to_d, [ zeros(num_vertices(p)) for _ in faces ] )
  for (i,(Πi,fi)) in enumerate(zip(Π,faces))
    v_to_d = Π_to_v_to_d[i+num_Π]
    for v in 1:num_vertices(p)
      if fi ∈ v_to_f[v]
        dist = 0.0
      else
        dist = signed_distance(p[v],Πi)
        if abs(dist) < atol
          dist = 0.0
        end
      end
      v_to_d[v] = dist
    end
  end
end

function edge_mesh(p::Polyhedron)
  edge_to_vertices = Vector{Int32}[]
  for v in 1:num_vertices(p)
    for vneig in get_graph(p)[v]
      if vneig > v
        push!(edge_to_vertices,[v,vneig])
      end
    end
  end
  X = get_vertex_coordinates(p)
  T = Table(edge_to_vertices)
  UnstructuredGrid(X,T,[SEG2],fill(1,length(T)))
end

function compute_graph(stl::GridTopology{2})
  @notimplementedif only(get_polytopes(stl)) ≠ TRI
  D = 2
  v_to_e = get_faces(stl,0,1)
  v_to_f = get_faces(stl,0,D)
  f_to_v = get_faces(stl,D,0)
  ec = array_cache(v_to_e)
  fc = array_cache(v_to_f)
  vc = array_cache(f_to_v)
  graph = [ Int32[] for _ in 1:num_vertices(stl) ]
  for v in 1:num_vertices(stl)
    vfacets = getindex!(fc,v_to_f,v)
    !isempty(vfacets) || continue
    f0 = first(vfacets)
    fnext = f0
    while true
      i = findfirst(isequal(v), getindex!(vc,f_to_v,fnext) )
      inext = i == D+1 ? 1 : i+1
      vnext = getindex!(vc,f_to_v,fnext)[inext]
      push!(graph[v],vnext)
      faces = getindex!(fc,v_to_f,vnext)
      i = findfirst( f -> f ≠ fnext && v ∈ getindex!(vc,f_to_v,f), faces )
      fnext = isnothing(i) ? UNSET : faces[i]
      fnext ≠ UNSET || break
      fnext ≠ f0 || break
    end
    if fnext == UNSET
      n_edges = length( getindex!(ec,v_to_e,v) )
      if length(graph[v]) < n_edges
        v0 = first(graph[v])
        fnext = f0
        while fnext ≠ UNSET
          i = findfirst(isequal(v), getindex!(vc,f_to_v,fnext) )
          inext = i == 1 ? D+1 : i-1
          vnext = getindex!(vc,f_to_v,fnext)[inext]
          pushfirst!(graph[v],vnext)
          faces = getindex!(fc,v_to_f,vnext)
          i = findfirst( f -> f ≠ fnext && v ∈ getindex!(vc,f_to_v,f), faces )
          fnext = isnothing(i) ? UNSET : faces[i]
        end
        @assert length(graph[v]) == n_edges
      end
      push!(graph[v],OPEN)
    end
  end
  check_graph(graph) || error("Unable to build vertex-edge graph")
  graph
end

function set_original_faces!(p::Polyhedron,a...)
  set_original_faces!(get_data(p),a...)
  p
end

function polyhedron_data(num_vertices::Integer)
  v_to_Π = [ Int32[] for _ in 1:num_vertices ]
  v_to_of = [ Int32[] for _ in 1:num_vertices ]
  v_to_v = collect(1:num_vertices)
  v_to_e = fill((UNSET,UNSET),num_vertices)
  Π_to_v_to_d = Vector{Int32}[]
  Π_to_rΠ = Int32[]
  Π_to_id = Int32[]
  args = v_to_Π, v_to_of, v_to_v, v_to_e, Π_to_v_to_d, Π_to_rΠ, Π_to_id
  PolyhedronData(args...)
end

function polyhedron_data(p::Polytope)
  v_to_Π = map( i-> Int32.(i), -get_faces(p,0,num_dims(p)-1) )
  v_to_of = [ Int32[] for _ in 1:num_vertices(p) ]
  v_to_v = Int32.(collect(1:num_vertices(p)))
  v_to_e = fill((Int32(UNSET),Int32(UNSET)),num_vertices(p))
  Π_to_v_to_d = Vector{Int32}[]
  Π_to_rΠ = Int32[]
  Π_to_id = Int32[]
  args = v_to_Π, v_to_of, v_to_v, v_to_e, Π_to_v_to_d, Π_to_rΠ, Π_to_id
  PolyhedronData(args...)
end

@inline function split_postprocess!(graph,data,new_vertices,input_poly,distances,(≶))
  complete_graph!(graph,num_vertices(input_poly))
  disconnect_graph!(graph,num_vertices(input_poly),distances,(≶))
  add_open_vertices!(graph,input_poly)
  vertices = [input_poly.vertices;new_vertices]
  update_data!(data,graph,num_vertices(input_poly))
  poly = Polyhedron(vertices,graph,isopen(input_poly),data)
  compact!(poly)
end

function complete_graph!(edge_graph,num_vertices::Integer)
  for v in num_vertices+1:length(edge_graph)
    vnext = next_vertex(edge_graph,num_vertices,v)
    vnext ∉ (UNSET,OPEN) || continue
    edge_graph[v][end] = vnext
    edge_graph[vnext][2] = v
  end
end

function add_open_vertices!(graph,poly::Polyhedron)
  if isopen(poly)
    for v in num_vertices(poly)+1:length(graph)
      if UNSET ∉ graph[v]
        push!(graph[v],graph[v][end])
        graph[v][end-1] = OPEN
      end
    end
  end
end

function disconnect_graph!(edge_graph,num_vertices,distances,(≶)::Function)
  for i in 1:num_vertices
    if distances[i] ≶ 0
      edge_graph[i] = empty!( edge_graph[i] )
    end
  end
end

function update_data!(data,graph,num_vertices)
  v_to_pv = data.vertex_to_parent_vertex
  v_to_pe = data.vertex_to_parent_edge
  Π_to_v_to_d = get_plane_distances(data)
  for v in num_vertices+1:length(graph)
    v_to_pv[v] == v || continue
    for vnext in num_vertices+1:length(graph)
      v_to_pv[vnext] == vnext || continue
      if v_to_pe[v] == v_to_pe[vnext] && v_to_pe[v] ≠ (UNSET,UNSET)
        v_to_pv[vnext] = v
        for v_to_d in Π_to_v_to_d
          v_to_d[v] = v_to_d[v_to_pv[v]]
        end
      end
    end
  end
end

function next_vertex(edge_graph::Vector,num_vertices::Integer,vstart::Integer)
  vcurrent = vstart
  vnext = first( edge_graph[vcurrent] )
  while vnext ≤ num_vertices && vnext ∉ (UNSET,OPEN)
    i = findfirst( isequal(vcurrent), edge_graph[vnext] )
    vcurrent = vnext
    inext = ( i % length( edge_graph[vcurrent] ) ) + 1
    vnext = edge_graph[vcurrent][ inext ]
  end
  vnext
end

function compute_intersection(p1::Point,d1::Real,p2::Point,d2::Real)
  ( p1*(abs(d2))+p2*(abs(d1)) ) / ( abs(d1) + abs(d2) )
end

add_vertex!(data::Nothing,a...) = data

function _intersect(a::AbstractVector{T},b::AbstractVector{T}) where T
  n = count(in(a),b)
  r = Vector{T}(undef,n)
  n = 0
  for i in a
    if i ∈ b
      n += 1
      r[n] = i
    end
  end
  r
end

function _copy(a::AbstractVector{<:AbstractVector})
  [ copy(i) for i in a ]
end

function add_vertex!(data::PolyhedronData,v1::Integer,v2::Integer,Πid::Integer)
  v_to_Π = data.vertex_to_planes
  v_to_f = data.vertex_to_original_faces
  v_to_pv = data.vertex_to_parent_vertex
  v_to_pe = data.vertex_to_parent_edge
  Πs = _intersect( v_to_Π[v1], v_to_Π[v2] )
  fs = _intersect( v_to_f[v1], v_to_f[v2] )
  pe = (v_to_pv[v1],v_to_pv[v2])
  push!( Πs, Πid )
  push!( v_to_Π, Πs )
  push!( v_to_f, fs )
  push!( v_to_pe, pe)
  Π_to_v_to_d = get_plane_distances(data)
  Π_to_id = get_plane_ids(data)
  i = findfirst(isequal(Πid),Π_to_id)
  d1 = Π_to_v_to_d[i][v1]
  d2 = Π_to_v_to_d[i][v2]
  pv = iszero(d1) ? v_to_pv[v1] : length(v_to_pv)+1
  pv = iszero(d2) ? v_to_pv[v2] : pv
  push!(v_to_pv,pv)
  @assert d1 ≠ d2
  ref_i = data.plane_to_ref_plane[i]
  for i in 1:length(Π_to_v_to_d)
    v_to_d = Π_to_v_to_d[i]
    if pv ≠ length(v_to_pv)
      dist = v_to_d[pv]
    elseif ref_i == UNSET || data.plane_to_ref_plane[i] ∈ (UNSET,i)
      if Π_to_id[i] == Πid
        dist = 0.0
      else
        dist = (d1*v_to_d[v2]-d2*v_to_d[v1])/(d1-d2)
      end
    else
      j = data.plane_to_ref_plane[i]
      @assert abs(j) < i
      dist = sign(j) * Π_to_v_to_d[abs(j)][end]
    end
    push!(v_to_d,dist)
  end
  data
end

function set_original_faces!(data::PolyhedronData,stl::GridTopology)
  Dc = num_dims(stl)
  v_to_f = data.vertex_to_original_faces
  for d in 0:Dc
    offset = get_offset(stl,d)
    _v_to_f = get_faces(stl,0,d)
    c = array_cache(_v_to_f)
    for v in 1:length(v_to_f)
      append!( v_to_f[v], map( i -> i.+offset, getindex!(c,_v_to_f,v) ) )
    end
  end
  data
end

function next_vertex(p::Polyhedron,vprevious::Integer,vcurrent::Integer)
  i = findfirst( isequal(vprevious), get_graph(p)[vcurrent] )
  i = ( i % length( get_graph(p)[vcurrent] ) ) + 1
  get_graph(p)[vcurrent][ i ]
end

## Kd Tree stuff

function refine_by_vertices(Γ::Polyhedron,K::Polyhedron,atol=0)
  refine_by_vertices(Γ,K,get_vertex_coordinates(Γ),atol)
end

function refine_by_vertices(Γ,K,vertices,atol)
  i = findfirst( v -> has_intersection(K,v), vertices )
  if isnothing(i)
    R = [Γ],[K]
  else
    v = vertices[i]
    d = get_direction(v,K)
    Π,Π⁻,Π⁺ = get_planes(v,d,atol)
    Πid,Πid⁻,Πid⁺ = get_new_plane_ids(Γ)
    compute_distances!(Γ,[Π⁻,Π⁺],[Πid⁻,Πid⁺])
    compute_distances!(K,[Π],[Πid])
    _,γ⁻ = split(Γ,Πid⁺)
    _,γ⁺ = split(Γ,Πid⁻)
    k⁻,k⁺ = split(K,Πid)
    round_distances!(k⁻,atol/10)
    round_distances!(k⁺,atol/10)
    Γn = [γ⁻,γ⁺]
    Kn = [k⁻,k⁺]
    if length(vertices) == i
      R = Γn,Kn
    else
      Γr,Kr = empty(Γn),empty(Kn)
      vertices = view(vertices,i+1:length(vertices))
      for (γ,k) in zip(Γn,Kn)
        Γi,Ki = refine_by_vertices(γ,k,vertices,atol)
        append!(Γr,Γi)
        append!(Kr,Ki)
      end
      R = Γr,Kr
    end
  end
  R
end

function get_bounding_box(cell::Polyhedron)
  i = findfirst(i->isactive(cell,i),1:num_vertices(cell))
  pmin,pmax = cell[i],cell[i]
  for i in 1:num_vertices(cell)
    if isactive(cell,i)
      pmin = Point( min.(Tuple(cell[i]),Tuple(pmin)) )
      pmax = Point( max.(Tuple(cell[i]),Tuple(pmax)) )
    end
  end
  pmin,pmax
end

function has_intersection(cell::Polyhedron,v::Point)
  pmin,pmax = get_bounding_box(cell)
  all( Tuple(pmin) .< Tuple(v) ) || return false
  all( Tuple(pmax) .> Tuple(v) ) || return false
  true
end

function get_direction(v,cell)
  pmin,pmax = get_bounding_box(cell)
  d0 = Tuple(v) .- Tuple(pmin)
  d1 = Tuple(pmax) .- Tuple(v)
  d = min.(d0,d1)
  d_max = maximum(d)
  findfirst(isequal(d_max),d)
end

function get_planes(v,d,atol)
  δ = VectorValue(Base.setindex(Tuple(zero(v)),atol,d))
  Π = CartesianPlane(v,d,+1)
  Π⁻ = CartesianPlane(v-δ,d,+1)
  Π⁺ = CartesianPlane(v+δ,d,-1)
  Π,Π⁻,Π⁺
end

function get_new_plane_ids(poly)
  id = minimum( get_plane_ids(poly.metadata) )
  id-1,id-2,id-3
end

function round_distances!(poly,atol)
  Π_to_v_to_dist = poly.metadata.plane_to_vertex_to_distances
  for i in 1:length(Π_to_v_to_dist)
    v_to_dist = Π_to_v_to_dist[i]
    for v in 1:length(v_to_dist)
      if abs(v_to_dist[v]) < atol
        v_to_dist[v] = 0.0
      end
    end
  end
end

function delete_inactive_planes!(Γ::Polyhedron,K::Polyhedron,stl::GridTopology)
  @assert isopen(Γ)
  planes = get_active_planes(Γ,K,stl)
  restrict_planes!(Γ,planes)
  restrict_planes!(K,planes)
end

function get_active_planes(Γ::Polyhedron,K::Polyhedron,stl::GridTopology)
  used_planes_Γ = _get_used_planes(K)
  used_planes_K = _get_used_planes(Γ)
  used_planes = lazy_append( used_planes_Γ, used_planes_K )
  face_planes = _get_face_planes(Γ,stl)
  planes = lazy_append( used_planes, face_planes )
  _append_ref_planes( Γ, planes )
end

function restrict_planes!(poly::Polyhedron,planes)
  restrict_planes!(get_data(poly),planes)
end

function restrict_planes!(data::PolyhedronData,planes)
  Π_to_id = data.plane_to_ids
  Π_to_v_to_dist = data.plane_to_vertex_to_distances
  Π_to_ref_Π = data.plane_to_ref_plane
  is_inactive = trues(length(Π_to_id))
  is_touch = falses(length(Π_to_id))
  for Π in planes
    i = findfirst(isequal(Π),Π_to_id)
    !isnothing(i) || continue
    is_inactive[i] = false
  end
  id = 0
  for i in 1:length(Π_to_id)
    if !is_inactive[i]
      id += 1
      if !is_touch[i]
        ref_plane = abs( Π_to_ref_Π[i] )
        if ref_plane != UNSET
          @assert i == ref_plane
          for j in i:length(Π_to_id)
            if !is_inactive[j] &&
               !is_touch[j] &&
               abs( Π_to_ref_Π[j] ) == ref_plane

              _sign = sign( Π_to_ref_Π[j] )
              Π_to_ref_Π[j] = _sign * id
              is_touch[j] = true
            end
          end
        end
      end
    end
  end
  deleteat!( Π_to_v_to_dist, is_inactive)
  deleteat!( Π_to_ref_Π, is_inactive)
  deleteat!( Π_to_id, is_inactive)
end

function _get_used_planes(poly::Polyhedron)
  v_to_Π = get_data(poly).vertex_to_planes
  planes = Int32[]
  for v in 1:length(v_to_Π)
    for Π in v_to_Π[v]
      if Π ∉ planes
        push!(planes,Π)
      end
    end
  end
  planes
end

function _get_face_planes(poly::Polyhedron,stl::GridTopology)
  rfaces = get_original_reflex_faces(poly,stl,empty=true)
  facets = get_original_facets(poly,stl,empty=true)
  lazy_append(rfaces,facets)
end

function _append_ref_planes(poly::Polyhedron,planes)
  ref_planes = Int32[]
  Π_to_ref_Π = get_data(poly).plane_to_ref_plane
  Π_to_id = get_plane_ids( get_data(poly) )
  for Π in planes
    i = findfirst( isequal(Π), Π_to_id )
    i !== nothing || continue
    iref = abs(Π_to_ref_Π[i])
    if iref ≠ UNSET
      Π_ref = Π_to_id[iref]
      if Π_ref ∉ ref_planes && Π_ref ∉ planes
        push!(ref_planes,Π_ref)
      end
    end
  end
  lazy_append(planes,ref_planes)
end
