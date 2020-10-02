

## Quadratic cost version
function compute_face_to_cells(grid::Grid,stl::DiscreteModel)
  stl_face_to_cells = [ Int[] for _ in 1:num_faces(stl) ] 
  cell_to_stl_faces = [ Int[] for _ in 1:num_cells(grid) ]
  for cell in 1:num_cells(grid)
    for stl_facet in 1:num_cells(stl)
      stl_face = stl_facet + get_offset(get_grid_topology(stl),num_dims(stl))
      if have_intersection(grid,cell,stl,stl_face)
        push!(stl_face_to_cells[stl_face],cell)
        push!(cell_to_stl_faces[cell],stl_face)
      end
    end
  end
  stl_face_to_cells, cell_to_stl_faces
end

function refine_surface(grid::Grid,stl::DiscreteModel,stl_face_to_cells)
  @notimplementedif num_dims(grid) ≠ num_dims(stl)+1
  T = Vector{Int}[]
  X = empty( get_vertex_coordinates(get_grid_topology(stl)) )
  f_to_bgcell = Int[]
  for stl_facet in 1:num_cells(stl)
    stl_face = get_dimrange(get_grid_topology(stl),num_dims(stl))[stl_facet]
    for cell in stl_face_to_cells[stl_face]
      if is_on_boundary(grid,cell,stl,stl_face) && 
         is_cell_out(grid,cell,stl,stl_face)
        continue
      end
      Tk,Xk = refine_stl_face(grid,cell,stl,stl_face)
      Tk = [ f .+ length(X) for f in Tk ]
      append!(T,Tk)
      append!(X,Xk)
      append!(f_to_bgcell, fill(cell,length(Tk)) )
    end
  end
  T = Table(T)
  stl = compute_stl_model(T,X)
  stl,f_to_bgcell
end

function is_cell_out(
  grid::Grid,
  cell::Integer,
  m::DiscreteModel,
  face::Integer)

  @assert num_dims(m) == get_facedims(get_grid_topology(m))[face]
  facet = get_offset(get_grid_topology(m),num_dims(m))[face]
  f = get_cell(m,facet)
  plane = center(f),normal(f)
  max_dist = 0.0
  for v in get_cell_coordinates(cell)
    dist = signed_distance(v,plane)
    if abs(dist) > abs(max_dist)
      max_dist = dist
    end
  end
  if max_dist < 0
    true
  else
    false
  end
end

function refine_stl_face(grid::Grid,cell::Integer,m::DiscreteModel,f::Integer)
  @assert num_dims(m) == get_facedims(get_grid_topology(m))[f]
  ifacet = f - get_offset(get_grid_topology(m),num_dims(m))
  facet = get_cell(m,ifacet)
  p = get_polytope(get_cell_reffes(grid)[cell])
  X = empty(get_vertex_coordinates(get_grid_topology(m)))
  point_to_face = Int[]
  cell_face_to_points = [ Int[] for _ in 1:num_faces(p) ]
  D = num_dims(grid)
  for face in 1:num_faces(facet)
    d = get_facedims(get_polytope(facet))[face]
    for cell_face in 1:num_faces(p)
      if get_facedims(p)[cell_face] ≤ D-d &&
         distance(grid,cell,cell_face,facet,face) < TOL &&
         !is_face_in_cell_face(
           p,
           cell_face,
           get_polytope(facet),
           face,
           cell_face_to_points,
           point_to_face)

        point = intersection_point(grid,cell,cell_face,facet,face)
        push!(X,point)
        push!(point_to_face,face)
        push!(cell_face_to_points[cell_face], length(X) )
      end
    end
  end
  faces = [ Vector{Int}[] for _ in 0:D-1 ]
  cell_face_to_ds = [ Int8[] for _ in 1:num_faces(p) ]
  cell_face_to_dfaces = [ Int[] for _ in 1:num_faces(p) ]
  for face in 1:num_faces(facet)
    d = get_facedims(get_polytope(facet))[face]
    dface = face - get_offset(get_polytope(facet),d)
    for n in 0:d
      for cell_face in 1:num_faces(p,D-(d-n))
        points = _get_points(p,D-(d-n),cell_face,get_polytope(facet),d,dface,
                             cell_face_to_points,point_to_face)
        if length(points) == n+1
          push!(faces[n+1],points)
          push!(cell_face_to_ds[cell_face],n)
          push!(cell_face_to_dfaces[cell_face],length(faces[n+1]))
        elseif length(points) > n+1

          dfaces = _get_dfaces(n-1,p,D-(d-n),cell_face,
                               cell_face_to_ds,cell_face_to_dfaces)
          for df in dfaces
            if points[1] ∉ faces[n][df]
              push!(faces[n+1],[points[1];faces[n][df]])
              push!(cell_face_to_ds[cell_face],n)
              push!(cell_face_to_dfaces[cell_face],length(faces[n+1]))
            end
          end
        end
      end
    end
  end
  T = faces[D]
  T,X
end

function is_face_in_cell_face(
  cell_polytope::Polytope,
  cell_face::Integer,
  polytope::Polytope,
  face::Integer,
  cell_face_to_points,
  point_to_face)

  d = get_facedims(polytope)[face]
  dface = face - get_offset(polytope,d)
  for cf in get_faces(cell_polytope)[cell_face]
    for point in cell_face_to_points[cf]
      f = point_to_face[point]
      n = get_facedims(polytope)[f]
      nf = f - get_offset(polytope,n)
      for df in get_faces(polytope,n,d)[nf]
        if df == dface
          return true
        end
      end
    end
  end
  false
end

function _get_points(
  cell_polytope::Polytope,
  cell_face_d::Integer,
  cell_dface::Integer,
  polytope::Polytope,
  d::Integer,
  dface::Integer,
  cell_face_to_points,
  point_to_face)

  points = Int[]
  cell_face = get_dimrange(cell_polytope,cell_face_d)[cell_dface]
  for cf in get_faces(cell_polytope)[cell_face]
    for point in cell_face_to_points[cf]
      f = point_to_face[point]
      n = get_facedims(polytope)[f]
      nf = f - get_offset(polytope,n)
      if n ≤ d
        for df in get_faces(polytope,n,d)[nf]
          if df == dface
            push!(points,point)
          end
        end
      end
    end
  end
  points
end

function _get_dfaces(
  dim::Integer,
  p::Polytope,
  d::Integer,
  dface::Integer,
  cell_face_to_ds,
  cell_face_to_dfaces)

  dfaces = Int[]
  face = get_dimrange(p,d)[dface]
  for f in get_faces(p)[face]
    for (_d,_dface) in zip(cell_face_to_ds[f],cell_face_to_dfaces[f])
      if _d == dim
        push!(dfaces,_dface)
      end
    end
  end
  dfaces
end

function distance(
  grid::Grid,
  cell::Integer,
  cell_face::Integer,
  facet::Face,
  face::Integer) 

  c = get_cell(grid,cell)
  distance(c,cell_face,facet,face)
end

function intersection_point(
  grid::Grid,
  cell::Integer,
  cell_face::Integer,
  facet::Face,
  face::Integer) 

  c = get_cell(grid,cell)
  intersection_point(c,cell_face,facet,face)
end








