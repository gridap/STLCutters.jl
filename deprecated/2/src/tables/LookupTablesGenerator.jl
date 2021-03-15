module LookupTablesGenerator

using MiniQhull
using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using Gridap.Helpers

num_combinations(p::Polytope) = 2^num_vertices(p)


function get_node_value(p::Polytope,node::Integer,case::Integer)
  ( (case-1) & (1<<(node-1)) ) ≠ 0
end

function is_cell_cut(K,p::Polytope,case::Integer)
  face = num_faces(p)
  is_face_cut(K,p,face,case)
end

function is_facet_cut(K,p::Polytope,facet::Integer,case::Integer)
  D = num_dims(p)
  face = get_dimrange(p,D-1)[facet]
  is_face_cut(K,p,face,case)
end

function is_edge_cut(K,p::Polytope,edge::Integer,case::Integer)
  face = get_dimrange(p,1)[edge]
  is_face_cut(K,p,face,case)
end

function is_edge_cut(p::Polytope,edge::Integer,case::Integer)
  K = 1:num_vertices(p)
  is_edge_cut(K,p,edge,case)
end

function is_face_cut(K,p::Polytope,face::Integer,case::Integer)
  nodes = get_face_vertices(p)[face]
  value1 = get_node_value(p,K[nodes[1]],case)
  for node in nodes
    value = get_node_value(p,K[node],case)
    if value ≠ value1
      return true
    end
  end
  false
end

function is_cell_in(K,p::Polytope,case::Integer)
  for i in 1:num_vertices(p)
    node = K[i]
    if node ≤ num_vertices(p)
      if get_node_value(p,node,case)
        return true
      else
        return false
      end
    end
  end
  @assert false
end

function is_cell_out(K,p::Polytope,case::Integer)
  !is_cell_cut(K,p::Polytope,case)
end

function get_face_extrusion(p::Polytope,face::Integer)
  p.dface.nfaces[face].extrusion
end

function get_edge_extrusion(p::Polytope,edge::Integer)
  face = get_dimrange(p,1)[edge]
  get_face_extrusion(p,face)
end

function num_cut_edges(p::Polytope,case::Integer)
  K = 1:num_vertices(p)
  num_cut_edges(K,p,case)
end

function num_cut_edges(K,p::Polytope,case::Integer)
  c = 0
  for edge in 1:num_faces(p,1)
    if is_edge_cut(K,p,edge,case)
      c += 1
    end
  end
  c
end

function need_wedges(K,p::Polytope,case::Integer)
  is_cell_cut(K,p,case) && !are_cut_edges_parallel(K,p,case)
end


function are_cut_edges_parallel(p::Polytope,case::Integer)
  K = 1:num_vertices(p)
  are_cut_edges_parallel(K,p,case)
end

function are_cut_edges_parallel(K,p::Polytope,case::Integer)
  face = num_faces(p)
  are_cut_edges_parallel(K,p,face,case)
end


function are_cut_edges_parallel(K::AbstractVector,p::Polytope,face::Integer,case::Integer)
  @assert is_n_cube(p)
  d = get_facedims(p)[face]
  dface = face-get_dimrange(p,d)[1]+1
  ext0 = get_edge_extrusion(p,1)
  first_cut = true
  for edge in get_faces(p,d,1)[dface] 
    if is_edge_cut(K,p,edge,case)
      ext = get_edge_extrusion(p,edge)
      if first_cut
        first_cut = false
        ext0 = ext
      elseif ext ≠ ext0
        return false
      end
    end
  end
  true
end

function get_cut_edge_extrusion(K,p::Polytope,case::Integer)
  for edge in 1:num_faces(p,1)
    if is_edge_cut(K,p,edge,case)
      return get_edge_extrusion(p,edge)
    end
  end
  @assert false
end

function get_face_at_wedge_vertex(K,p::Polytope,case::Integer)
  D = num_dims(p)
  for dface in 1:num_faces(p,D-2)
    face = get_dimrange(p,D-2)[dface]
    if !is_face_cut(K,p,face,case) 
      facets = get_faces(p,D-2,D-1)[dface]
      if is_facet_cut(K,p,facets[1],case) && is_facet_cut(K,p,facets[2],case)
        return face
      end
    end
  end
  for dface in 1:num_faces(p,D-2)
    facets = get_faces(p,D-2,D-1)[dface]
    if is_facet_cut(K,p,facets[1],case) && is_facet_cut(K,p,facets[2],case)
      return get_dimrange(p,D-2)[dface]
    end
  end
  @assert false
end

function get_facets_at_wedge(K,p::Polytope,case::Integer)
  face = get_face_at_wedge_vertex(K,p::Polytope,case::Integer)
  d = get_facedims(p)[face]
  dface = face-get_dimrange(p,d)[1]+1
  D = num_dims(p)
  @assert d == D-2
  facets = get_faces(p,D-2,D-1)[dface]
  @assert length(facets) == 2
  (facets[1],facets[2])
end

function opposite_facet(p::Polytope,facet::Integer)
  @assert is_n_cube(p)
  ((facet-1)⊻1) + 1 
end

function are_facets_opposite(p::Polytope,facet1::Integer,facet2::Integer)
  D = num_dims(p)
  face1 = get_dimrange(p,D-1)[facet1]
  face2 = get_dimrange(p,D-1)[facet2]
  ext1 = get_face_extrusion(p,face1)
  ext2 = get_face_extrusion(p,face2)
  ext1 == ext2
end

function get_wedge(p::Polytope,facet1::Integer,facet2::Integer)
  @assert !are_facets_opposite(p,facet1,facet2)
  D = num_dims(p)
  wedge_nodes = fill(UNSET,num_vertices(p))
  _,perm = _compute_facet_permutations(wedge_nodes,facet1,facet2,p)
  for (i,node) in enumerate( get_faces(p,D-1,0)[facet1] )
    wedge_nodes[i] = node
  end
  for i in 1:length( get_faces(p,D-1,0)[facet2] )
    inode = get_face_vertex_permutations(p,D-1)[facet2][perm][i] 
    node = get_faces(p,D-1,0)[facet2][inode]
    wedge_nodes[i+length(get_faces(p,D-1,0)[facet1])] = node
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

function get_wedges(p::Polytope,facet1::Integer,facet2::Integer)
  wedge1 = get_wedge(p::Polytope,facet1,facet2)
  facet1 = opposite_facet(p,facet1)
  facet2 = opposite_facet(p,facet2)
  wedge2 = get_wedge(p::Polytope,facet1,facet2)
  wedge1,wedge2
end

function num_vertices_in(p::Polytope,case::Integer)
  c = 0
  for node in 1:num_vertices(p)
    if get_node_value(p,node,case) == 1
      c += 1
    end
  end
  c
end

function num_vertices_out(p::Polytope,case::Integer)
  c = 0
  for node in 1:num_vertices(p)
    if get_node_value(p,node,case) == 0
      c += 1
    end
  end
  c
end

function first_vertex_in(p::Polytope,case::Integer)
  for node in 1:num_vertices(p)
    if get_node_value(p,node,case) == 1
      return node
    end
  end
  UNSET
end

function first_vertex_out(p::Polytope,case::Integer)
  for node in 1:num_vertices(p)
    if get_node_value(p,node,case) == 0
      return node
    end
  end
  UNSET
end

function is_cut_by_one_plane(p::Polytope,case::Integer)
  cache = Int[]
  are_in_nodes_connected(cache,p,case) &&
  are_out_nodes_connected(cache,p,case)
end

function are_in_nodes_connected(cache,p::Polytope,case::Integer)
  n_nodes_in = num_vertices_in(p,case)
  if n_nodes_in > 0
    queue = cache
    resize!(queue,1)
    head = 1
    queue[1] = first_vertex_in(p,case)
    node_to_edges = get_faces(p,0,1)
    edge_to_nodes = get_faces(p,1,0)
    while length(queue) ≥ head
      node = queue[head]
      head += 1
      for edge in node_to_edges[node]
        if !is_edge_cut(p,edge,case)
          for _node in edge_to_nodes[edge]
            if _node ≠ node && _node ∉ queue
              push!(queue,_node)
            end
          end
        end
      end
    end
    @assert length(queue) ≤ n_nodes_in
    if length(queue) == n_nodes_in
      true
    else
      false
    end
  else
    true
  end
end

function are_out_nodes_connected(cache,p::Polytope,case::Integer)
  n_nodes_out = num_vertices_out(p,case)
  if n_nodes_out > 0
    queue = cache
    resize!(queue,1)
    head = 1
    queue[1] = first_vertex_out(p,case)
    node_to_edges = get_faces(p,0,1)
    edge_to_nodes = get_faces(p,1,0)
    while length(queue) ≥ head
      node = queue[head]
      head += 1
      for edge in node_to_edges[node]
        if !is_edge_cut(p,edge,case)
          for _node in edge_to_nodes[edge]
            if _node ≠ node && _node ∉ queue
              push!(queue,_node)
            end
          end
        end
      end
    end
    @assert length(queue) ≤ n_nodes_out
    if length(queue) == n_nodes_out
      true
    else
      false
    end
  else
    true
  end
end

function are_intersections_coplanar(p::Polytope,case::Integer)
  tol = 1e-3
  K = 1:num_vertices(p)
  num_intersections = num_cut_edges(p,case)
  if num_intersections == 0
    return true
  end
  @assert num_intersections ≥ num_dims(p)
  intersections = zeros(eltype(get_vertex_coordinates(p)),num_intersections)
  i = 0
  for edge in 1:num_faces(p,1)
    if is_edge_cut(K,p,edge,case)
      nodes = get_face_vertices(p,1)[edge]
      i += 1
      v1 = get_vertex_coordinates(p)[nodes[1]]
      v2 = get_vertex_coordinates(p)[nodes[2]]
      intersections[i] = (v1+v2)/2
    end
  end
  o = intersections[1]
  n = normal(intersections[1:num_dims(p)])
  for point in intersections
    if abs((point-o)⋅n) > tol
      return false
    end
  end
  true
end

function normal(points::Vector{<:Point{2}}) 
  @notimplementedif length(points) ≠ 2
  v = points[2]-points[1]
  n = VectorValue(-v[2],v[1])
  n / norm(n)
end

function normal(points::Vector{<:Point{3}}) 
  @notimplementedif length(points) ≠ 3
  v1 = points[2]-points[1]
  v2 = points[3]-points[1]
  n = v1 × v2
  n / norm(n)
end

function normal(::Vector{<:Point})
  @notimplemented
end

function compute_linear_grid_and_facemap(reffe::LagrangianRefFE)
  grid = compute_linear_grid(reffe)
  desc = get_cartesian_descriptor(grid)
  model = CartesianDiscreteModel(desc)
  labels = get_face_labeling(model)
  grid_face_to_reffe_face = get_face_entity(labels)
  grid,grid_face_to_reffe_face
end

function get_vertex_references(reffe::LagrangianRefFE)
  p = get_polytope(reffe)
  vertices = Vector{Vector{Int}}(undef,num_nodes(reffe))
  for face in 1:num_faces(reffe)
    if length(get_face_own_nodes(reffe)[face]) > 0
      @assert length(get_face_own_nodes(reffe)[face]) == 1
      own_node = get_face_own_nodes(reffe)[face][1]
      vertices[own_node] = copy(get_face_vertices(p)[face])
    end
  end
  vertices
end

function get_vertex_references(p::Polytope)
  reffe = LagrangianRefFE(p)
  get_vertex_references(reffe)
end

function get_new_vertex_references(reffe::LagrangianRefFE)
  p = get_polytope(reffe)
  vertices = Vector{Vector{Int}}(undef,num_nodes(reffe)-num_vertices(reffe))
  for face in num_vertices(reffe)+1:num_faces(reffe)
    if length(get_face_own_nodes(reffe)[face]) > 0
      @assert length(get_face_own_nodes(reffe)[face]) == 1
      own_node = get_face_own_nodes(reffe)[face][1]
      inode = own_node - num_vertices(reffe)
      vertices[inode] = copy(get_face_vertices(p)[face])
    end
  end
  vertices
end
      
function compute_new_cells(reffe::ReferenceFE)
  grid, gface_to_rface = compute_linear_grid_and_facemap(reffe)
  num_nodes_per_cell = num_vertices(reffe)
  new_cells = Vector{Vector{Int}}(undef,num_cells(grid))
  for lcell in 1:num_cells(grid)
    new_cell = fill(UNSET,num_nodes_per_cell)
    for lnode in 1:num_nodes_per_cell
      gnode = get_cell_nodes(grid)[lcell][lnode] 
      face = gface_to_rface[gnode]
      rnode = get_face_own_nodes(reffe)[face][1]
      new_cell[lnode] = rnode
    end
    new_cells[lcell] = new_cell
  end
  new_cells
end

function update_cells!(cells,K,p::Polytope,n_vertices::Integer)
  for cell in cells
    for (i,lnode) in enumerate(cell)
      if lnode > num_vertices(p)
        node = lnode - num_vertices(p) + n_vertices
      else
        node = K[lnode]
      end
      cell[i] = node
    end
  end
  cells
end

function split_cell(K,X,p,case)
  @assert are_cut_edges_parallel(K,p,case) 
  ext = get_cut_edge_extrusion(K,p,case)
  partition = Tuple(ext).+1
  reffe = LagrangianRefFE(Float64,p,partition)
  cells = compute_new_cells(reffe)
  vertices = get_new_vertex_references(reffe)
  update_cells!(cells,K,p,length(X))
  update_cells!(vertices,K,p,length(X))
  cells,vertices
end

function compute_wedges(K,X,p::Polytope,case::Integer)
  facets = get_facets_at_wedge(K,p,case)
  wedges = get_wedges(p,facets...)
  update_cells!(wedges,K,p,length(X))
  wedges
end


function initialize_mesh(p::Polytope)
  X = get_vertex_references(p)
  K = 1:num_vertices(p)
  T = [K]
  Tnew = Vector{Int}[]
  T,X,Tnew
end

function refine_cells!(T,X,p,case,Tnew)
  for K in T
    if need_wedges(K,p,case)
      wedges = compute_wedges(K,X,p,case)
      refine_cells!(wedges,X,p,case,Tnew)
    else
      if is_cell_cut(K,p,case)
        cells,vertices = split_cell(K,X,p,case) 
        append!(Tnew,cells)
        append!(X,vertices)
      else
        push!(Tnew,K)
      end
    end
  end
end

function simplex_refinement(p::Polytope,case::Integer)
  Ts,ps = simplexify(p)
  X = get_vertex_coordinates(p)
  X_to_ref = [ [i] for i in 1:num_vertices(p) ]
  T = empty(Ts)
  for Ks in Ts
    Tk,Xk = refine_cell(X,Ks,ps,case,X_to_ref)
    append!(T,Tk)
    append!(X_to_ref,Xk)
  end
  T,X_to_ref = update_vertices(T,X_to_ref)
  T,X_to_ref,ps
end

function refine_cell(X,K,p,case,X_to_ref)
  Xk = X[K]
  Xnew = empty(X_to_ref)
  for edge in 1:num_edges(p)
    if is_edge_cut(K,p,edge,case)
      nodes = get_face_vertices(p,1)[edge]
      point = sum( Xk[nodes] ) / length(nodes)
      push!(Xk,point)
      push!(Xnew,K[nodes])
    end
  end
  coordinates = zeros(num_dims(p)*length(Xk))
  i = 0
  for point in Xk, d in 1:num_dims(p)
    i += 1
    coordinates[i] = point[d]
  end
  T = delaunay(num_dims(p),length(Xk),coordinates)
  Tnew = Vector{Int}[]
  for i in 1:size(T,2)
    Knew = T[:,i]
    for (j,node) in enumerate(Knew)
      if node ≤ length(K)
        Knew[j] = K[node]
      else
        Knew[j] = length(X_to_ref) + node - length(K)
      end
    end
    push!(Tnew,Knew)
  end
  Tnew,Xnew
end

function update_vertices(T,X)
  _map = collect(1:length(X))
  for i in 1:length(X), j in i+1:length(X)
    if X[i] ⊆ X[j] && X[i] ⊇ X[j]
      _map[j] = _map[i]
    end
  end
  Xnew = X[unique(_map)]
  m = _map .== 1:length(_map)
  m = cumsum(m) .* m
  _map = m[_map]
  Tnew = map( i -> _map[i], T )
  Tnew,Xnew
end

function simplex_to_n_cubes(T,p::Polytope)
  Tnew = empty(T)
  for K in T
    push!(Tnew,K[ simplex_to_n_cube(p) ])
  end
  Tnew
end

function simplex_to_n_cube(p::Polytope{1})
  @assert is_simplex(p)
  [1,2]
end

function simplex_to_n_cube(p::Polytope{2})
  @assert is_simplex(p)
  [1,2,3,3]
end

function simplex_to_n_cube(p::Polytope{3})
  @assert is_simplex(p)
  [1,2,3,3,4,4,4,4]
end

function compute_grid(
  cell_to_nodes::AbstractArray,
  node_to_coordinates::Vector{<:Point},
  p::Polytope)

  T = Table(cell_to_nodes)
  X = node_to_coordinates
  reffes = [LagrangianRefFE(p)]
  cell_types = fill(1,length(T))
  UnstructuredGrid(X,T,reffes,cell_types)
end

function compute_vertex_coordinates(X,p::Polytope)
  vertices = zeros(eltype(get_vertex_coordinates(p)),length(X))
  for (i,nodes) in enumerate(X)
    for node in nodes
      point = get_vertex_coordinates(p)[node]
      vertices[i] += point/length(nodes)
    end
  end
  vertices
end

function compute_in_out(T,p::Polytope,case::Integer)
  in_out = fill(Int8(UNSET),length(T))
  for (i,K) in enumerate(T)
    if is_cell_in(K,p,case)
      in_out[i] = 1
    else
      in_out[i] = 0
    end
  end
  in_out
end

function signed_volume(K,X::Vector{<:Point},p::Polytope)
  simplices,ps = simplexify(p)
  Xk = X[K]
  vol = 0.0
  for simplex in simplices
    rvol = simplex_signed_volume(simplex,get_vertex_coordinates(p),ps)
    vol += simplex_signed_volume(simplex,Xk,ps) * sign(rvol)
  end
  vol
end

function correct_volumes(T,X,p::Polytope)
  @assert is_n_cube(p)
  vol = 0.0
  fv = get_face_vertices(p,num_dims(p)-1)
  Tnew = deepcopy(T)
  for (i,K) in enumerate(T)
    vol_k = signed_volume(K,X,p)
    if vol_k < 0
      Kp = [K[fv[2]];K[fv[1]]]
      vol_k = signed_volume(Kp,X,p)
      Tnew[i] = Kp
    end
    @assert vol_k > 0
    vol += vol_k
  end
  @assert vol ≈ 1
  Tnew
end

function simplex_signed_volume(K,X,p::Polytope{D}) where D
  @assert is_simplex(p)
  factor = 1/factorial(D)
  v = ntuple( i -> X[K[i+1]]-X[K[1]], Val{D}() )
  data = ntuple( i -> v[ ((i-1)÷D)+1 ][ ((i-1)%D)+1 ], Val{D^2}() )
  t = TensorValue(data...)
  vol = det(t) * factor
  vol
end


function compute_tables(p::Polytope,write=false)
  case_to_table_id = fill(UNSET,num_combinations(p))
  table_id_to_connectivities = Vector{Vector{Int}}[]
  table_id_to_vertices = Vector{Vector{Int}}[] 
  table_id_to_cell_to_inout = Vector{Int8}[]
  table_id = 0
  for case in 1:num_combinations(p)
    if is_cut_by_one_plane(p,case) 

      table_id += 1
      case_to_table_id[case] = table_id
      
      if are_intersections_coplanar(p,case)
        T,X,Tnew = initialize_mesh(p)
        refine_cells!(T,X,p,case,Tnew)
        T = Tnew
      else
        T,X,ps = simplex_refinement(p,case)
        T = simplex_to_n_cubes(T,ps)
      end

      Xp = compute_vertex_coordinates(X,p)
      T = correct_volumes(T,Xp,p)

      in_out = compute_in_out(T,p,case)

      push!( table_id_to_connectivities, T )
      push!( table_id_to_vertices, X )
      push!( table_id_to_cell_to_inout, in_out)

      if write 
        println("$p case $case: $(length(T)) cells, $(length(X)) vertices")
        grid = compute_grid(T,Xp,p)
        writevtk(grid,"$(num_dims(p))D_case_$case",cellfields=["IN_OUT"=>in_out])
      end
    end
  end
  case_to_table_id,
  table_id_to_connectivities,
  table_id_to_vertices,
  table_id_to_cell_to_inout
end

function compute_tables(dims=2:3,write=false)
  d_to_case_to_table_id = Vector{Int}[]
  d_to_table_id_to_connectivities = Vector{Vector{Vector{Int}}}[]
  d_to_table_id_to_vertices = Vector{Vector{Vector{Int}}}[] 
  d_to_table_id_to_cell_to_inout = Vector{Vector{Int8}}[] 
  for d in dims
    p = Polytope(tfill(HEX_AXIS,Val(d)))
    output = compute_tables(p,write)

    case_to_table_id,
    table_id_to_connectivities,
    table_id_to_vertices,
    table_id_to_cell_to_inout = output

    push!( d_to_case_to_table_id, case_to_table_id )
    push!( d_to_table_id_to_connectivities, table_id_to_connectivities )
    push!( d_to_table_id_to_vertices, table_id_to_vertices )
    push!( d_to_table_id_to_cell_to_inout, table_id_to_cell_to_inout )
  end

  d_to_case_to_table_id,
  d_to_table_id_to_connectivities,
  d_to_table_id_to_vertices,
  d_to_table_id_to_cell_to_inout
end


const header = 
"# Lookup tables generated automatically
# Do not modify by hand!"

const getters = 
"
function is_case_possible(dim::Integer,case::Integer)
  d_to_case_to_table_id[dim-1][case] ≠ UNSET
end

function get_connectivities_from_case(dim::Integer,case::Integer)
  table_id = d_to_case_to_table_id[dim-1][case]
  d_to_table_id_to_connectivities[dim-1][table_id]
end

function get_vertex_to_cell_vertices_from_case(dim::Integer,case::Integer)
  table_id = d_to_case_to_table_id[dim-1][case]
  d_to_table_id_to_vertices[dim-1][table_id]
end
"

function write_tables(filename,write_ouput=false)
  dims = 2:3
  tables = compute_tables(dims,write_ouput)

  d_to_case_to_table_id,
  d_to_table_id_to_connectivities,
  d_to_table_id_to_vertices,
  d_to_table_id_to_cell_to_inout = tables

  f = open(filename,"w")

  println(f,header)
  println(f,getters)


  _print_data(f,dims,d_to_case_to_table_id,"#_to_case_to_table_id")

  _print_data(f,dims,d_to_table_id_to_connectivities,"#_to_table_id_to_connectivities")
  
  _print_data(f,dims,d_to_table_id_to_vertices,"#_to_table_id_to_vertices")
  
  _print_data(f,dims,d_to_table_id_to_cell_to_inout,"#_to_table_id_to_cell_to_inout")

  close(f)
end

function _print_data(f::IO,dims,tables,name)
  for (d,table) in zip(dims,tables)
    println(f,"const "*replace(name,'#'=>"d$d")*" =")
    println(f,table)
    println(f)
  end
  println(f,"const "*replace(name,'#'=>"d")*" =")
  print(f," [ ")
  for d in dims
    print(f,replace(name,'#'=>"d$d")*", ")
  end
  println(f,"]")
  println(f)
end

end # module
