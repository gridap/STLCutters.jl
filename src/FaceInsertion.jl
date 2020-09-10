
function vertex_refinement(
  cell_nodes,
  node_to_coordinates::Vector{<:Point},
  p::Polytope,
  point::Point{D}) where D

  d = farthest_axis_from_boundary(cell_nodes,node_to_coordinates,point)
  n = tfill(1,Val{D}())
  n = Base.setindex(n,2,d)
  reffe = LagrangianRefFE(Float64,p,n)
  new_cells = compute_new_cells(cell_nodes,node_to_coordinates,reffe)
  new_vertices = compute_new_vertices(cell_nodes,node_to_coordinates,reffe,point)
  new_cells, new_vertices
end

function edge_refinement(
  K,
  X::Vector{<:Point},
  p::Polytope,
  e::Segment,
  directions)

  plane = compute_plane_from_edge(K,X,p,e,directions)
  case = compute_case(K,X,p,plane)
  Tnew = compute_new_cells(K,X,p,case)
  Xnew = compute_new_vertices(K,X,p,plane,case)
  Tnew,Xnew
end

## Helpers

function compute_case(K,X,p,plane)
  case = 0
  for node in 1:num_vertices(p)
    point = X[K[node]]
    dist = signed_distance(point,plane)
    if dist > 0
      case |= (1<<(node-1))
    end
  end
  case += 1
  @assert is_case_possible(num_dims(p),case)
  case
end

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

function compute_new_vertices(K,X,p,plane,case)
  v_to_cv = get_vertex_to_cell_vertices_from_case(num_dims(p),case)
  vertices = eltype(X)[]
  for i in num_vertices(p)+1:length(v_to_cv)
    nodes = v_to_cv[i]
    @assert length(nodes) == 2    
    p1 = X[K[nodes[1]]]
    p2 = X[K[nodes[2]]]
    d1 = signed_distance(p1,plane)
    d2 = signed_distance(p2,plane)

    if abs(d1) < TOL || abs(d2) < TOL
      vertex = abs(d1) < abs(d2) ? p1 : p2
      #update T, delete v
    else
      α = abs(d1) / (abs(d1)+abs(d2))
      vertex = p1 + (p2-p1)*α
    end
    push!(vertices,vertex)
  end
  vertices
end

function compute_new_cells(
  cell_nodes::Vector{<:Integer},
  node_to_coordinates::Vector{<:Point},
  reffe::ReferenceFE)

  grid, gface_to_rface = compute_linear_grid_and_facemap(reffe)
  num_nodes_per_cell = length(cell_nodes)
  num_nodes = length(node_to_coordinates)
  new_cells = Vector{Int}[]
  for lcell in 1:num_cells(grid)
    new_cell = fill(UNSET,num_nodes_per_cell)
    for lnode in 1:num_nodes_per_cell
      gnode = get_cell_nodes(grid)[lcell][lnode] 
      face = gface_to_rface[gnode]
      rnode = get_face_own_nodes(reffe)[face][1]
      if rnode ≤ num_nodes_per_cell
        n = cell_nodes[rnode]
      else
        n = rnode - num_nodes_per_cell + num_nodes
      end
      new_cell[lnode] = n
    end
    push!(new_cells,new_cell)
  end
  new_cells
end

function compute_new_vertices(
  cell_nodes::Vector{<:Integer},
  node_to_coordinates::Vector{<:Point},
  reffe::ReferenceFE,
  point::Point{D}) where D

  p = get_polytope(reffe)
  new_node_to_coordinates = Vector{typeof(point)}(undef,num_nodes(reffe)-num_vertices(p))
  for face in num_vertices(p)+1:num_faces(p)
    if length(get_face_own_nodes(reffe)[face]) > 0
      node = get_face_own_nodes(reffe)[face][1]
      vertex = compute_vertex_coordinates(cell_nodes,node_to_coordinates,p,face,point)
      new_node_to_coordinates[node-num_vertices(p)] = vertex
    end
  end
  new_node_to_coordinates
end


function compute_vertex_coordinates(
  cell_nodes,
  node_to_coordinates,
  p::Polytope{D},
  iface::Integer,
  point::Point{D}) where D

  nface = p.dface.nfaces[iface]
  dim = p.dface.dims[iface]
  node = cell_nodes[ get_faces(p)[iface][1] ]
  anchor = node_to_coordinates[ node ]
  extrusion = nface.extrusion
  vertex = anchor
  for d in 1:D
    if extrusion[d] == HEX_AXIS
      v = point[d]
      vertex = Base.setindex(vertex,v,d)
    end
  end
  vertex
end

function compute_linear_grid_and_facemap(reffe::LagrangianRefFE)
  grid = compute_linear_grid(reffe)
  desc = get_cartesian_descriptor(grid)
  model = CartesianDiscreteModel(desc)
  labels = get_face_labeling(model)
  grid_face_to_reffe_face = get_face_entity(labels)
  grid,grid_face_to_reffe_face
end

function compute_facemap(grid::CartesianGrid)
  desc = get_cartesian_descriptor(grid)
  model = CartesianDiscreteModel(desc)
  labels = get_face_labeling(model)
  grid_face_to_reffe_face = get_face_entity(labels)
  grid_face_to_reffe_face
end

function get_default_directions(E,STL_edges::Vector{<:Segment{D}}) where D
  acc_v = 0
  for e in E
    edge = STL_edges[e]
    v = edge[2] - edge[1]
    acc_v += v
  end
  if D == 2
    ()
  elseif D == 3
    _,d = findmin(acc_v.data)
    v = zero( VectorValue{D,Int} )
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

function perturbation(K,T::Vector{<:Point},p::Polytope,e::Segment)
  point = nothing
  for facet in 1:num_facets(p)
    if have_intersection(K,T,p,facet,e)
      point = intersection_point(K,T,p,facet,e)
      break
    end
  end
  @assert point !== nothing
  min_dist = Inf
  closest_facet = UNSET
  for facet in 1:num_facets(p)
    if !have_intersection(K,T,p,facet,e)
      dist = distance(K,T,p,facet,p)
      if dist < min_dist
        min_dist = dist
        closest_facet = facet
      end
    end
  end
  normal(K,T,p,closest_facet)
end

function signed_distance(point::Point,plane::Tuple{VectorValue,VectorValue})
  o,n = plane
  (point-o)⋅n
end
