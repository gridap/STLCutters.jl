
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
  cell_nodes,
  node_to_coordinates::Vector{<:Point{D}},
  edge::Segment) where D

  new_cells = typeof(cell_nodes)[]
  p = Polytope(tfill(HEX_AXIS,Val{D}()))
  face1,face2,point1,point2 = compute_intersections(cell_nodes,node_to_coordinates,p,edge)
  if are_sharing_a_facet(face1,face2,p)
    return new_cells
  end
  facet1,facet2 = get_farthest_facets(face1,face2,p)
  if !are_facets_opposite(facet1,facet2,p)
    wedge1,wedge2 = get_wedges(cell_nodes,facet1,facet2,p)
    facet1,facet2 = (1,opposite_facet(1,p))
    cell_nodes = wedge1
    push!(new_cells,wedge2)
  end
  offset = first(get_dimrange(p,num_dims(p)-1))
  d = findfirst(isequal(1),p.dface.nfaces[offset+facet1].extrusion.data)
  n = tfill(1,Val{D}())
  n = Base.setindex(n,2,d)
  reffe = LagrangianRefFE(Float64,p,n)
  new_cells = [ new_cells; compute_new_cells(cell_nodes,node_to_coordinates,reffe) ]
  new_vertices = [point1,point2]
  new_cells, new_vertices
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
  for face in num_vertices(p)+1:num_faces(reffe)
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

function are_sharing_a_facet(face1::Integer,face2::Integer,p::Polytope)
  D = num_dims(p)
  dim1 = get_facedims(p)[face1]
  dim2 = get_facedims(p)[face2]
  nface1 = face1 - first(get_dimrange(p,dim1)) + 1
  nface2 = face2 - first(get_dimrange(p,dim2)) + 1
  for facet1 in get_faces(p,dim1,D-1)[nface1]
    for facet2 in get_faces(p,dim2,D-1)[nface2]
      if facet1 == facet2
        return true
      end
    end
  end
  false
end

function get_farthest_facets(face1::Integer,face2::Integer,p::Polytope)
  D = num_dims(p)
  dim1 = get_facedims(p)[face1]
  dim2 = get_facedims(p)[face2]
  nface1 = face1 - first(get_dimrange(p,dim1)) + 1
  nface2 = face2 - first(get_dimrange(p,dim2)) + 1
  for facet1 in get_faces(p,dim1,D-1)[nface1]
    for facet2 in get_faces(p,dim2,D-1)[nface2]
      @assert facet1 ≠ facet2
      if are_facets_opposite(facet1,facet2,p)
        return (facet1,facet2)
      end
    end
  end
  facet1 = first(get_faces(p,dim1,D-1)[nface1])
  facet2 = first(get_faces(p,dim2,D-1)[nface2])
  facet1,facet2
end

function are_facets_opposite(facet1::Integer,facet2::Integer,p::Polytope)
  D = num_dims(p)
  offset = first(get_dimrange(p,D-1))-1
  ext1 = p.dface.nfaces[offset+facet1].extrusion
  ext2 = p.dface.nfaces[offset+facet2].extrusion
  ext1 == ext2
end

function get_partition(facet1::Integer,facet2::Integer,p::Polytope)
  D = num_dims(p)
  offset = first(get_dimrange(p,D-1))-1
  ext1 = p.dface.nfaces[offset+facet1].extrusion
  ext2 = p.dface.nfaces[offset+facet2].extrusion
  @assert ext1 == ext2
  (ext1+1).data
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

function compute_intersections(cell_nodes,node_to_coordinates,p::Polytope,e::Segment)
  D = num_dims(p)
  point1 = zero(eltype(node_to_coordinates))
  point2 = zero(eltype(node_to_coordinates))
  face1 = UNSET
  face2 = UNSET
  for facet in 1:num_facets(p)
    if have_intersection(cell_nodes,node_to_coordinates,p,facet,e)
      point = intersection_point(cell_nodes,node_to_coordinates,p,facet,e)
      for d in 0:D-1
        dfaces = get_faces(p,D-1,d)[facet]
        dface_to_nodes = get_face_vertices(p,d) 
        for dface in dfaces
          if distance(cell_nodes,node_to_coordinates,p,d,dface,point) < TOL
            point = projection(cell_nodes,node_to_coordinates,p,d,dface,point)
            face = get_dimrange(p,d)[dface]
            if face1 == UNSET
              face1 = face
              point1 = point
            elseif face2 == UNSET
              face2 = face
              point2 = point
            else
              @assert false
            end
          end
        end
      end
    end
  end
  @assert face1 ≠ UNSET && face2 ≠ UNSET
  (face1,face2,point1,point2)
end
