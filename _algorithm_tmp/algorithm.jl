

This algorithm is not true julia code, therefore there are typos, incomplete function arguments and syntax inconsistencies.


## Input Data

surface_mesh, bg_mesh
sm_face_to_bg_cells
bg_mface_to_bgnface # implicit

## Output Data (to be grow during the refinement)

d_to_dface_to_vertices = [ Table{Int}() for d in 0:D-1 ]
d_to_dface_to_sm_n = [ Int[] for d in 0:D-1 ]
d_to_dface_to_sm_nface = [ Int[] for d in 0:D-1 ]

vertex_to_bg_d = [ Int[] for d in 0:D ]
vertex_to_bg_dface = [ Int[] for d in 0:D ]

## Main Algorithm

SMFace_to_BgCells = compute_surface_mesh_face_to_bg_cells( bg_mesh, surface_mesh )

for vertex in vertices(sm)
  face = find_closest_face( bg_mesh, vertex, SMFace_to_BgCells )
  add_face!( face )
end

for d in 1:D-1
  for SMFace in faces(d,surface_mesh)
    
    n = d-1
    bg_cells = cells_containg_nfaces( bg_mesh, _sm, n )
    
    while length(bg_cells) > 0
      bg_cell = pop!(bg_cells)
        
      ## 1) Fetch sm_face-bg_cell intersection already computed
      ##    ( ∂Face_SM ∩ Cell_BgM ) ∩ ( 2) of neighbor cells )

      nfaces = fetch_boundary_nfaces( bg_cell, _sm, n )

      ## 2) Complete intersections in cell boundary, to define a close domain
      ##    Face_SM ∩ ∂Cell_BgM 
      
      append!(nfaces, complete_boundary(bg_cell,_sm,SMFace) )

      ## 3) Now we have a convex domain on n-dim space bounded by nfaces
      ##    Create symplex mesh locally (dface = nface ∪ vertex):

      add_faces!( sm, nfaces, n )
      
      ## *) Feed depth-first traversal

      for _bg_cell in cells_touching( bg_mesh, nfaces ) #new nfaces actually, what is computed in 2)
        push!( bg_cells, _bg_cell )
      end

    end
  end
end


## Functions

## 1) Fetch sm_face-bg_cell intersection already computed
##    ( ∂Face_SM ∩ Cell_BgM ) ∩ ( 2) of neighbor cells )

function fetch_boundary_nfaces( bg_cell, _sm, n )
  ## NOTE: for efficienty we may store, bg_face to new_sm_faces
  resize!( nfaces, 0 )
  for nface in faces(_sm, n)
    if nface ∈ bg_cell
      push!( nfaces, nface )
    end
  end
  nfaces
end


## 2) Complete intersections in cell boundary, to define a close domain
##    Face_SM ∩ ∂Cell_BgM 

function complete_boundary( bg_cell, _sm, sm_face, nfaces )
  d = face_dimension( sm_face )
  n = d-1

  vertices = get_vertices(nfaces)

  while length(vertices) > 0
    vertex = pop!(vertices)
    bg_face, bg_d = get_background_face(vertex)
    _bg_face,_d,point = _find_next_intersection(bg_cell,bg_face,bg_d,vertex,d)
    new_nface = vertex ∩ point
    add_face(new_nface)
    push!(nfaces, new_nface )
    push!(point)
  end
end


function _find_next_intersection(bg_cell,bg_face,d,vertex,sm_d)
  min_dist = tol()
  closest_bg_face = UNSET

  for _bg_face in faces(bg_cell)
    if sm_d == 1 || is_around( bg_face, _bg_face )
      if _bg_face ∉ bg_face &&
         _bg_face ∉ any(nfaces)

        dist = distance(_bg_face,sm_face) #_bg_face ∩ sm_face  ≠ ∅
        if dist < min_dist 
          min_dist = dist
          closest_bg_face = _bg_face
        end
      end
    end
  end

  if closest_bg_face != UNSET
    point = intersection(closest_bg_face,sm_face)
    return closest_bg_face, point
  end
  return UNSET,zero(Point)
end


## 3) Now we have a convex domain on n-dim space bounded by nfaces
##    Create symplex mesh locally (dface = nface ∪ vertex):

function add_faces!( _sm, nfaces, n )
  new_d_to_dfaces = compute_local_simplex_mesh( nfaces, n )
  for d in 0:n+1
    append!( _sm, d, new_d_to_dfaces[d] )
  end
  _sm
end

function compute_local_simplex_mesh( nfaces, n )
  vertex = get_vertex( nfaces[1], 1 ) # Arbitrary vertex
  d = n+1
  for nface in nfaces
    if vertex ∉ nface
      dface = nface ∩ vertex
      push!( dfaces, dface )
    end
  end
  for _d in 1:n
    _compute_connectivities(_d,dfaces) # integer based
  end
  new_d_to_dfaces
end

## Other functions
##

function find_closest_bg_face( bg_mesh, surface_mesh, vertex, sm_face_to_bg_cells )
  min_dist = tol()
  closest_cell = UNSET
  for cell in sm_face_to_bg_cells(vertex)
    dist = distance(cell, vertex)
    if dist < min_dist
      min_dist = dist
      closest_cell = cell
    end
  end
  closest_cell != UNSET || return (UNSET,UNSET)
  min_dist = tol()
  closest_dface = UNSET
  for d in 0:D
    for dface in faces(cell,d)
      dist = distance( dface, vertex )
      if dist < min_dist
        min_dist = dist
        closest_dface = dface
      end
    end
    if closest_dface != UNSET
      return d, closest_dface
    end
  end
  @assert false
end










