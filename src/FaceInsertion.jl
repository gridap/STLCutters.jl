
function vertex_refinement(
  cell_nodes,
  node_to_coordinates::Vector{<:Point},
  point::Point{D}) where D

  p = Polytope(tfill(HEX_AXIS,Val{D}()))
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
  X::Vector{<:Point{D}},
  e::Segment,
  directions) where D

  p = Polytope(tfill(HEX_AXIS,Val{D}()))
  plane = compute_plane_from_edge(K,X,p,e,directions)
  face0 = get_dface_perpendicular_to_vectors(K,X,p,2,directions)
  faces = compute_intersected_faces(K,X,p,face0,plane)
  if are_sharing_a_facet(face0,faces,p)
    return Vector{Int}[], eltype(X)[]
  end
  facets = get_farthest_facets(face0,faces,p)
  Tnew = get_connectivities(K,X,p,facets,directions)

  Xnew = compute_new_vertices(K,X,p,plane) 

  Tnew,Xnew
end

## Helpers

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

function are_facets_opposite(facet1::Integer,facet2::Integer,p::Polytope)
  D = num_dims(p)
  offset = first(get_dimrange(p,D-1))-1
  ext1 = p.dface.nfaces[offset+facet1].extrusion
  ext2 = p.dface.nfaces[offset+facet2].extrusion
  ext1 == ext2
end

function get_wedge(cell_nodes::Vector,facet1::Integer,facet2::Integer,p::Polytope)
  @assert !are_facets_opposite(facet1,facet2,p)
  D = num_dims(p)
  wedge_nodes = fill(UNSET,num_vertices(p))
  _,perm = _compute_facet_permutations(wedge_nodes,facet1,facet2,p)
  for (i,node) in enumerate( get_faces(p,D-1,0)[facet1] )
    wedge_nodes[i] = cell_nodes[ node ]
  end
  for i in 1:length( get_faces(p,D-1,0)[facet2] )
    inode = get_face_vertex_permutations(p,D-1)[facet2][perm][i] 
    node = get_faces(p,D-1,0)[facet2][inode]
    wedge_nodes[i+length(get_faces(p,D-1,0)[facet1])] = cell_nodes[ node ]
  end
  wedge_nodes
end

function _compute_facet_permutations(cache,facet1::Integer,facet2::Integer,p::Polytope)
 D = num_dims(p)
 fill!(cache,UNSET)
 for (i,node1) in enumerate( get_faces(p,D-1,0)[facet1] )
    for (j,node2) in enumerate( get_faces(p,D-1,0)[facet2] )
      if node1 == node2
        cache[i] = j
      end
    end
  end
  for (i,perm) in enumerate( get_face_vertex_permutations(p,D-1)[facet2] )
    next_perm = false
    for j in 1:length(perm)
      if cache[j] ≠ UNSET && cache[j] ≠ perm[j]
        next_perm = true
        break
      end
    end
    if !next_perm
      return (1,i)
    end
  end
  @assert false
end

function get_opposite_facets(facet1::Integer,facet2::Integer,p::Polytope)
  op1 = opposite_facet(facet1)
  op2 = opposite_facet(facet2)
  if op1 < op2
    op1,op2
  else
    op2,op1
  end
end

function opposite_facet(facet::Integer,p::Polytope)
  ((facet-1)⊻1) + 1 
end

function get_wedges(cell_nodes,facet1,facet2,p::Polytope)
  wedge1 = get_wedge(cell_nodes,facet1,facet2,p::Polytope)
  facet1 = opposite_facet(facet1,p)
  facet2 = opposite_facet(facet2,p)
  wedge2 = get_wedge(cell_nodes,facet1,facet2,p::Polytope)
  wedge1,wedge2
end

function get_connectivities(
  K,
  X::Vector{<:Point},
  p::Polytope,
  facets,vs::Tuple)

  new_cells = Vector{Int}[]
  if !are_facets_opposite(facets...,p)
    wedge1,wedge2 = get_wedges(K,facets...,p)
    facets = (1,opposite_facet(1,p))
    K = wedge1
    push!(new_cells,wedge2)
  end
  facet1,facet2 = facets
  partition = get_partition(K,X,facet1,p,vs) 
  reffe = LagrangianRefFE(Float64,p,partition)
  [ new_cells; compute_new_cells(K,X,reffe) ]
end

function compute_new_vertices(K,X::Vector{<:Point},p::Polytope,plane)
  vertices = eltype(X)[]
  edge_to_nodes = get_face_vertices(p,1)
  for edge in 1:num_faces(p,1)
    p1 = X[K[edge_to_nodes[edge][1]]]
    p2 = X[K[edge_to_nodes[edge][2]]]
    d1 = signed_distance(p1,plane)
    d2 = signed_distance(p2,plane)
    s1 = d1 ≥ 0 ? 1 : -1
    s2 = d2 ≥ 0 ? 1 : -1
    if s1 ≠ s2
      if abs(d1) < TOL || abs(d2) < TOL
        vertex = abs(d1) < abs(d2) ? p1 : p2
      else
        α = abs(d1) / (abs(d1)+abs(d2))
        vertex = p1 + (p2-p1)*α
      end
      push!(vertices,vertex)
    end
  end
  vertices
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

function get_dface_perpendicular_to_vectors(
  K,
  X::Vector{<:Point},
  p::Polytope,
  d::Integer,
  vs::Tuple)

  for dface in 1:num_faces(p,d)
    face = get_dimrange(p,d)[dface]
    if is_face_perpendicular_to_vectors(K,X,p,face,vs)
      return face
    end
  end
  @assert false
end

function compute_intersected_faces(
  K,
  X::Vector{<:Point},
  p::Polytope,
  face::Integer,
  plane)

  d = get_facedims(p)[face]
  @assert d == 2
  dface = face - get_dimrange(p,d)[1] + 1
  dface_edges = get_faces(p,d,1)[dface]
  edge_to_nodes = get_face_vertices(p,1)
  face1,face2 = UNSET,UNSET
  for edge in dface_edges
    p1 = X[K[edge_to_nodes[edge][1]]]
    p2 = X[K[edge_to_nodes[edge][2]]]
    d1 = signed_distance(p1,plane)
    d2 = signed_distance(p2,plane)
    s1 = d1 ≥ 0 ? 1 : -1
    s2 = d2 ≥ 0 ? 1 : -1
    if s1 ≠ s2
      if abs(d1) < TOL || abs(d2) < TOL
        if abs(d1) < abs(d2)
          face = edge_to_nodes[edge][1]
        else
          face = edge_to_nodes[edge][2]
        end
      else
        face = get_dimrange(p,1)[edge]
      end
      if face1 == UNSET
        face1 = face
      elseif face2 == UNSET
        face2 = face
      else
        @assert false
      end
    end
  end
  (face1,face2)
end

function get_farthest_facets(face0::Integer,faces,p::Polytope)
  @assert get_facedims(p)[face0] == 2
  D = num_dims(p)
  face1,face2 = faces
  dim1 = get_facedims(p)[face1]
  dim2 = get_facedims(p)[face2]
  nface1 = face1 - get_dimrange(p,dim1)[1] + 1
  nface2 = face2 - get_dimrange(p,dim2)[1] + 1
  for facet1 in get_faces(p,dim1,D-1)[nface1]
    if !is_facet_on_face(facet1,face0,p)
      for facet2 in get_faces(p,dim2,D-1)[nface2]
        if !is_facet_on_face(facet2,face0,p)
          @assert facet1 ≠ facet2
          if are_facets_opposite(facet1,facet2,p)
            return (facet1,facet2)
          end
        end
      end
    end
  end
  facet1 = first(get_faces(p,dim1,D-1)[nface1])
  facet2 = first(get_faces(p,dim2,D-1)[nface2])
  facet1,facet2
end

function are_sharing_a_facet(face0::Integer,faces,p::Polytope)
  @assert get_facedims(p)[face0] == 2
  D = num_dims(p)
  face1,face2 = faces
  dim1 = get_facedims(p)[face1]
  dim2 = get_facedims(p)[face2]
  nface1 = face1 - get_dimrange(p,dim1)[1] + 1
  nface2 = face2 - get_dimrange(p,dim2)[1] + 1
  for facet1 in get_faces(p,dim1,D-1)[nface1]
    if !is_facet_on_face(facet1,face0,p)
      for facet2 in get_faces(p,dim2,D-1)[nface2]
        if !is_facet_on_face(facet2,face0,p)
          if facet1 == facet2
            return true
          end
        end
      end
    end
  end
  false
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

function is_face_perpendicular_to_vectors(
  K,
  X::Vector{<:Point},
  p::Polytope,
  face::Integer,
  vs::Tuple)

  face_lnodes = get_face_vertices(p)[face]
  x1 = X[K[face_lnodes[1]]]
  for v in vs
    @assert maximum(v) == 1 && norm(v) == 1
    x0 = x1 ⋅ v
    for ln in face_lnodes
      x = X[K[ln]]
      if x ⋅ v ≠ x0
        return false
      end 
    end
  end
  true
end

function signed_distance(point::Point,plane::Tuple{VectorValue,VectorValue})
  o,n = plane
  (point-o)⋅n
end

function is_facet_on_face(facet::Integer,face::Integer,p::Polytope)
  D = num_dims(p)
  face1 = get_dimrange(p,D-1)[facet]
  face2 = face
  is_face_on_face(face1,face2,p)
end

function is_face_on_face(face1::Integer,face2::Integer,p::Polytope)
  d1 = get_facedims(p)[face1]
  d2 = get_facedims(p)[face2]
  dface1 = face1 - get_dimrange(p,d1)[1] + 1 
  dface2 = face2 - get_dimrange(p,d2)[1] + 1
  dface1 ∈ get_faces(p,d2,d1)[dface2]
end

function get_partition(K,X::Vector{<:Point},facet::Integer,p::Polytope{D},vs::Tuple) where D
  facet_edges = get_faces(p,D-1,1)[facet]
  for edge in facet_edges
    face = get_dimrange(p,1)[edge]
    if is_face_perpendicular_to_vectors(K,X,p,face,vs)
      ext = p.dface.nfaces[ face ].extrusion
      d = findfirst(isequal(1),ext.data)
      n = tfill(1,Val{D}())
      n = Base.setindex(n,2,d)
      return n
    end
  end
  @assert false
end


