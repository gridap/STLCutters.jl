
function vertex_refinement(
  cell_nodes,
  node_to_coordinates::Vector{<:Point},
  p::Polytope,
  point::Point{D}) where D

  cell = cell_nodes,node_to_coordinates,p
  d = farthest_axis_from_boundary(cell...,point)
  case = compute_case(cell...,d)
  new_cells = compute_new_cells(cell...,case)
  new_vertices, new_to_old_vertices = compute_new_vertices(cell...,point,d,case)
  new_cells, new_vertices #, new_to_old_vertices
end

function edge_refinement(
  K,
  X::Vector{<:Point},
  p::Polytope,
  e::Segment,
  directions)

  cell = K,X,p
  plane = compute_plane_from_edge(cell...,e,directions)
  case = compute_case(cell...,plane,directions)
  Tnew = compute_new_cells(cell...,case)
  Xnew,Xnew_to_old = compute_new_vertices!(Tnew,cell...,plane,case)
  delete_empty_cells!(Tnew,[X;Xnew],p)
  Tnew,Xnew #,Xnew_to_old
end

function facet_refinement(
  K,
  X::Vector{<:Point},
  p::Polytope,
  stl::DiscreteModel,
  facets::Vector{<:Integer})

  levelsets = get_levelsets(K,X,p,stl,facets)
  if length(levelsets) > 0
    mesh = CellMeshes.CellMesh(X,K,p)
    CellMeshes.compute_cell_mesh!(mesh,levelsets)

    Xnew = CellMeshes.get_vertex_coordinates(mesh)
    Tnew = CellMeshes.get_cell_to_vertices(mesh)
    cell_to_io = CellMeshes.get_cell_to_inout(mesh)
    node_to_io = CellMeshes.get_vertex_to_inout(mesh)
    Xnew_to_old = CellMeshes.get_vertex_to_cell_vertices(mesh)
    
    cell_type = TET_AXIS

    update_connectivities!(Tnew,K,X,p)
    update_vertex_coordinates!(Xnew,p)

    node_to_io = node_to_io[num_vertices(p)+1:end]
    for i in 1:length(Xnew_to_old),j in 1:length(Xnew_to_old[i])
      Xnew_to_old[i,j] = Xnew_to_old[i]+length(X)
    end
    deleteat!(Xnew_to_old,1:num_vertices(p))


    #Tnew = [ Tnew[i] for i in findall(cell_to_io .≠ CellMeshes.IGNORE_FACE) ]
    #deleteat!(cell_to_io,cell_to_io .== CellMeshes.IGNORE_FACE)
  else
    Tnew = [K]
    Xnew = empty(X)
    cell_to_io = [ define_cell(K,X,p,get_cell(stl,facets[1])) ]
    cell_type = HEX_AXIS
  end

  cell_types = fill(cell_type,length(Tnew))
  Tnew,Xnew,cell_to_io,cell_types
end

##

function compute_new_cells(K,X,p,case)
  cell_to_lnodes = get_connectivities_from_case(num_dims(p),case)
  Tnew = deepcopy(cell_to_lnodes)
  for (icell,cell) in enumerate(Tnew)
    for (inode,lnode) in enumerate(cell)
      if lnode > num_vertices(p)
        node = length(X) + lnode - num_vertices(p)
      else
        node = K[lnode]
      end
      Tnew[icell][inode] = node
    end
  end
  Tnew
end

## Vertex

function farthest_axis_from_boundary(
  cell_nodes,
  node_to_coordinates::Vector{<:Point},
  p::Polytope,
  point::Point)

  pmin,pmax = get_bounding_box(cell_nodes,node_to_coordinates,p)
  max_dists = max( point-pmin, pmax-point )
  _,d = findmax(max_dists.data)
  d
end

function compute_case(K,X,p::Polytope,d::Integer)
  @assert is_n_cube(p)
  case = 0
  for node in 1:num_vertices(p)
    if (node-1) & (1<<(d-1)) ≠ 0
      case |= (1<<(node-1))
    end
  end
  case + 1
end

function compute_new_vertices(
  K,
  X::Vector{<:Point},
  p::Polytope,
  point::Point,
  d::Integer,
  case::Integer)

  v_to_cv = get_vertex_to_cell_vertices_from_case(num_dims(p),case)
  vertices = zeros(eltype(X),length(v_to_cv)-num_vertices(p))
  vertex_to_older_vertices = Vector{Vector{Int}}(undef,length(v_to_cv)-num_vertices(p)) 
  ivertex = 0
  for i in num_vertices(p)+1:length(v_to_cv)
    nodes = v_to_cv[i]
    @assert length(nodes) == 2    
    p1 = X[K[nodes[1]]]
    p2 = X[K[nodes[2]]]

    v = i - num_vertices(p)
    vertex = Base.setindex(p1,point[d],d)
    vertices[v] = vertex
    vertex_to_older_vertices[v] = K[nodes]
  end
  vertices, vertex_to_older_vertices
end

## Edges

function get_default_directions(E,stl::DiscreteModel{Dc,Dp}) where {Dc,Dp}
  sum_v = 0
  for e in E
    edge = get_edge(stl,e)
    v = edge[2] - edge[1]
    sum_v += abs(v)
  end
  if Dp == 2
    ()
  elseif Dp == 3
    _,d = findmin(Tuple(sum_v))
    v = zero( VectorValue{Dp,Int} )
    v = Base.setindex(v,1,d)
    (v,)
  end
end

function compute_plane_from_edge(K,X,p,e::Segment,default_directions)
  v = e[2] - e[1]
  n = orthogonal(v,default_directions...)
  if norm(n) < TOL
    v_ϵ = perturbation(K,X,p,e)
    v += v_ϵ
    n = orthogonal(v,default_directions...)
    @assert norm(n) > TOL
  end
  n /= norm(n)
  o = center(e)
  (o,n)
end

function compute_case(K,X,p,plane,vs)
  case = 0
  ref_ds = reference_directions(K,X,p,vs)
  for node in 1:num_vertices(p)
    n = project_node(K,X,p,node,ref_ds)
    point = X[K[n]]
    dist = signed_distance(point,plane)
    if dist < 0
      case |= (1<<(node-1))
    end
  end
  case += 1
  @assert is_case_possible(num_dims(p),case)
  case
end

function compute_new_vertices!(T,K,X,p,plane,case;atol=TOL)
  v_to_cv = get_vertex_to_cell_vertices_from_case(num_dims(p),case)
  vertices = eltype(X)[]
  vertex_to_older_vertices = Vector{Int}[]
  ivertex = 0
  for i in num_vertices(p)+1:length(v_to_cv)
    nodes = v_to_cv[i]
    @assert length(nodes) == 2    
    p1 = X[K[nodes[1]]]
    p2 = X[K[nodes[2]]]
    d1 = signed_distance(p1,plane)
    d2 = signed_distance(p2,plane)

    if abs(d1) < atol || abs(d2) < atol
      if abs(d1) < abs(d2)
        new_node = K[nodes[1]]
      else
        new_node = K[nodes[2]]
      end
      old_node = length(X)+ivertex+1
      update_connectivities!(T,old_node=>new_node)
    else
      α = abs(d1) / (abs(d1)+abs(d2))
      vertex = p1 + (p2-p1)*α
      push!(vertices,vertex)
      push!(vertex_to_older_vertices,K[nodes])
      ivertex += 1
    end
  end
  vertices, vertex_to_older_vertices
end

function delete_empty_cells!(T,X,p)
  icell = 1
  while icell ≤ length(T)
    cell = T[icell]
    if is_cell_empty(cell,p)
      deleteat!(T,icell)
    else
      icell += 1
    end
  end
end

## Facets

function get_levelsets(
  K,
  X::Vector{<:Point},
  p::Polytope,
  stl::DiscreteModel,
  facets::Vector{<:Integer})

  masks = [ is_on_boundary(K,X,p,get_cell(stl,f),atol=TOL) for f in facets] 
  levelsets = [ _get_levelset(stl,f) for f in facets] 
  deleteat!(levelsets,masks)
end

function _get_levelset(stl::DiscreteModel{Df,Dp},facet::Integer) where {Df,Dp}
  f = get_cell(stl,facet)
  vertices = ntuple( i -> f[i], Val{Df+1}() )
  CellMeshes.SimplexFacet( vertices )
end

function update_connectivities!(Tnew,K,X::Vector{<:Point},p::Polytope)
  for (i,Knew) in enumerate(Tnew), (j,node) in enumerate(Knew)
    if node ≤ num_vertices(p)
      Tnew[i,j] = K[node]
    else
      Tnew[i,j] = node - num_vertices(p) + length(X)
    end
  end
  Tnew
end

function update_vertex_coordinates!(Xnew::Vector{<:Point},p::Polytope)
  for _ in 1:num_vertices(p)
    popfirst!(Xnew)
  end
  Xnew
end

function define_cell(
  K,
  X::Vector{<:Point},
  p::Polytope,
  f::Face)

  plane = center(f),normal(f)
  max_dist = 0.0
  for i in 1:num_vertices(p)
    v = X[K[i]]
    dist = signed_distance(v,plane)
    if abs(dist) > abs(max_dist)
      max_dist = dist
    end
  end
  if max_dist > 0
    CellMeshes.FACE_OUT
  else
    CellMeshes.FACE_IN
  end
end

## Helpers

function Base.setindex(a::VectorValue,val,idx::Integer) 
  VectorValue(Base.setindex(Tuple(a),val,idx))
end

Base.abs(a::VectorValue) = VectorValue( abs.(Tuple(a)) )

function perturbation(K,X::Vector{<:Point},p::Polytope,e::Segment)
  point = nothing
  for facet in 1:num_facets(p)
    face = get_dimrange(p,num_dims(p)-1)[facet]
    if have_intersection(K,X,p,face,e)
      point = intersection_point(K,X,p,face,e)
      break
    end
  end
  @assert point !== nothing
  min_dist = Inf
  closest_facet = UNSET
  for i in 1:num_facets(p)
    facet = get_facet(K,X,p,i)
    if surface(facet) > 0 && !have_intersection(facet,e)
      dist = distance(facet,point)
      if dist < min_dist
        min_dist = dist
        closest_facet = i
      end
    end
  end
  @assert closest_facet ≠ UNSET
  normal( get_facet(K,X,p,closest_facet) )
end

function signed_distance(point::Point,plane::Tuple{VectorValue,VectorValue})
  o,n = plane
  (point-o)⋅n
end

function reference_directions(K,X,p::Polytope,vs::NTuple{N}) where N
  ntuple( i -> _get_direction(K,X,p,vs[i]), Val{N}() )
end

function _get_direction(K,X,p,v)
  facet0 = _reference_facet(K,X,p,v)
  _get_direction(get_facet_normals(p)[facet0])
end

function _reference_facet(K,X,p,v)
  @assert is_n_cube(p)
  d = _get_direction(v)
  for facet in 1:num_facets(p)
    facet_nodes = get_face_vertices(p,num_dims(p)-1)[facet]
    pd = X[K[facet_nodes[1]]][d]
    next_facet = false
    for node in facet_nodes
      if X[K[node]][d] ≠ pd
        next_facet = true
        break
      end
    end
    if !next_facet
      return facet
    end
  end
  @assert false
end

function _get_direction(v::VectorValue)
  @assert norm(v) == 1
  @assert maximum(abs(v)) == 1
  findfirst( i -> abs(i) == 1, Tuple(v) )
end

function project_node(K,X::Vector{<:Point},p::Polytope,node::Integer,ref_ds::Tuple)
  @assert is_n_cube(p)
  n = node
  for d in ref_ds
    n = ( (n-1) & ~(1<<(d-1)) ) + 1
  end
  n
end

function update_connectivities!(T,node_to_node::Pair)
  old_node,new_node = node_to_node
  for (icell,cell) in enumerate(T), (inode,node) in enumerate(cell)
    if node == old_node
      T[icell][inode] = new_node
    elseif node > old_node
      T[icell][inode] -= 1
    end
  end
end

function is_cell_empty(K,p)
  face = num_faces(p)
  is_face_empty(K,p,face)
end

function is_face_empty(K,p::Polytope,face::Integer)
  num_facets = 0
  d = get_facedims(p)[face]
  dface = face - get_dimrange(p,d)[1] + 1
  nfaces = get_faces(p,d,d-1)[dface]
  for nface in nfaces
    if d-1 == 1
      edge = nface
      if !is_edge_empty(K,p,edge)
        num_facets += 1
      end
    else
      if !is_face_empty(K,p,get_dimrange(p,d-1)[nface])
        num_facets += 1
      end
    end
  end
  num_facets < d+1
end

function is_edge_empty(K,p,edge)
  nodes = get_face_vertices(p,1)[edge]
  K[nodes[1]] == K[nodes[2]]
end

