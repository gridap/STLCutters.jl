module CellMeshesTests

using STLCutter

# To be public
using STLCutter: reset!, compact!, compute_in_out!, num_dfaces, get_faces, dface_dimension, local_dface, global_dface, get_vertex, get_face, get_dface

# Private used in test
using STLCutter: _add_vertex!, UNSET, are_all_faces_defined

# Temporally included
using STLCutter: @check, get_cell_dface, get_vector, get_vector_bis, find_closest_face, intersection_tolerance, get_cache

using Test


function is_boundary_face(mesh::CellMesh,d::Integer,dface::Integer,n::Integer,nface::Integer)
  @check d > n 
  df_to_nf = get_faces(mesh,d,n)
  for lnface in 1:length(df_to_nf,dface)
    if df_to_nf[dface,lnface] == nface
      return true
    end
  end
  false
end

function get_opposite_face(mesh::CellMesh,d::Integer,dface::Integer,n::Integer,nface::Integer)
  opposite_face_dim = d-n-1
  df_to_of = get_faces(mesh,d,opposite_face_dim)
  of_to_v = get_dface_to_vertices(mesh,opposite_face_dim)
  nf_to_v = get_dface_to_vertices(mesh,n)

  for loface in 1:length(df_to_of,dface)
    opposite = true
    oface = df_to_of[dface,loface]
    for oface_lvertex in 1:length(of_to_v,oface), nface_lvertex in 1:length(nf_to_v,nface)
      if of_to_v[oface,oface_lvertex] == nf_to_v[nface,nface_lvertex]
        opposite = false
        break
      end
    end
    if opposite
      return oface
    end
  end
  @check false
  return UNSET
end

function find_first_boundary_face(mesh::CellMesh,cache,sm::SurfaceMesh,sm_face_dim::Integer,sm_face::Integer)

  D = num_dims(mesh)
  for d in 1:D-1
    for dface in 1:num_dfaces(mesh,d)
      if isactive(mesh,d,dface)
        if get_cell_dface(cache.cell,d,dface) != 0 # UNSET, is_cell_boundary_face
          for n in 1:d
            df_to_nf = get_faces(mesh,d,n)
            iface = UNSET
            min_distance = intersection_tolerance
            for lnface in 1:length(df_to_nf,dface)
              nface = df_to_nf[dface,lnface]
              dist = distance(mesh,n,nface,sm,sm_face_dim,sm_face) 
              if dist < min_distance
                min_distance = dist
                iface = nface
              end
            end
            if iface != UNSET
              return (d,iface)
            end
          end
        end
      end
    end
  end
  (UNSET,UNSET)
end

function add_surface_mesh_face!(mesh,sm,sm_face)
  cache = mesh.cache
  d = dface_dimension(sm,sm_face)
  sm_dface = local_dface(sm,sm_face,d)

  boundary_vertices = get_vector(cache)
  for n in 0:d-1 
    df_to_nf = get_faces(sm,d,n)
    for vertex in 1:num_vertices(mesh)
      sm_faces = cache.cell.vertex_to_surface_mesh_faces[vertex]
      for sm_f in sm_faces
        for sm_lnface in 1:length(df_to_nf,sm_dface)
          if sm_f == global_dface(sm,n,df_to_nf[sm_dface,sm_lnface])
            push!(boundary_vertices,vertex)
          end
        end
      end
    end
  end

  if length(boundary_vertices) == 0
    ## add one starting point
    dim,iface = find_first_boundary_face(mesh,cache,sm,d,sm_dface)
    if iface != UNSET
      point = closest_point(mesh,dim,iface, sm,d,sm_dface)
      add_vertex!(mesh,dim,iface,point,sm_face)
    end

    
  end

  D = num_dims(mesh)
  queue = get_vector_bis(cache)
  append!(queue,boundary_vertices)

  head = 1
  tail = length(queue)

  queue_dim = d -1
  sm_face_dim = d
  mesh_dim = D
  while head <= tail
    current_face = queue[head]
    head += 1

    iface = UNSET
    opposite_dim = UNSET

    for n in queue_dim+1:mesh_dim
      opposite_dim = n-queue_dim-1
      iface = UNSET
      min_distance = intersection_tolerance
      for nface in 1:num_dfaces(mesh,n)
        if isactive(mesh,n,nface)
          if is_boundary_face(mesh,n,nface,queue_dim,current_face)
            opposite_face = get_opposite_face(mesh,n,nface,queue_dim,current_face)
            dist = distance(mesh,opposite_dim,opposite_face, sm,d,sm_dface)
            if dist < min_distance
              min_distance = dist
              iface = opposite_face
            end
          end
        end
      end
      if iface != UNSET
        break
      end
    end
    
    if iface != UNSET
      point = closest_point(mesh,opposite_dim,iface, sm,d,sm_dface)
      add_vertex!(mesh,opposite_dim,iface,point,sm_face)
      push!(queue,num_vertices(mesh))
    end

  end
end

p1 = Point(1,2,3)
p2 = Point(4,5,6)

b = BoundingBox(p1,p2)

m = CellMesh(b)


p1 = Point(1,2,4)
p2 = Point(5,3,5)

box = BoundingBox(p1,p2)

reset!(m,box)

@test @allocated(reset!(m,box)) == 0


## 2D Point Intersection

p0 = Point(0.0,0.0)
p1 = Point(1.0,1.0)
box = BoundingBox(p0,p1)

mesh = CellMesh(box)

stl_points = [ Point(0.3,0.3), Point(0.25,0.5), Point(0.5,0.5), Point(0.75,0.4), Point(0.0,0.5)  ]

for point in stl_points[1:1]

  d,face = find_closest_face(mesh,point)

  if face != UNSET 
    _add_vertex!(mesh,d,face,point)
  end

end

compact!(mesh)

writevtk(mesh,"sub_mesh")

## 3D Point Intersection

p0 = Point(0.0,0.0,0.0)
p1 = Point(1.0,1.0,1.0)
box = BoundingBox(p0,p1)

mesh = CellMesh(box)

stl_points = [ 
  Point(0.3,0.3,0.0), 
  Point(0.25,0.5,0.1), 
  Point(0.4,0.6,0.0), 
  Point(0.75,0.4,0.2), 
  Point(0.0,0.5,0.0)  ]

for point in stl_points

  d,face = find_closest_face(mesh,point)
  
  if face != UNSET 
    _add_vertex!(mesh,d,face,point)
  end

end

compact!(mesh)

writevtk(mesh,"sub_mesh3")

## 2D Surface Mesh
function find_captured_surface_mesh_face_points!(mesh,sm,sm_face)
  cache = get_cache(mesh)
  d = dface_dimension(sm,sm_face)
  sm_dface = local_dface(sm,sm_face,d)
  vertices = get_vector(cache)
  for n in 0:d-1 
    df_to_nf = get_faces(sm,d,n)
    for vertex in 1:num_vertices(mesh)
      sm_faces = cache.cell.vertex_to_surface_mesh_faces[vertex]
      for sm_f in sm_faces
        for sm_lnface in 1:length(df_to_nf,sm_dface)
          if sm_f == global_dface(sm,n,df_to_nf[sm_dface,sm_lnface])
            push!(vertices,vertex)
          end
        end
      end
    end
  end
  vertices
end


function find_intersection_point_on_boundary_face(mesh::CellMesh,sm::SurfaceMesh,sm_face::Integer)
  D = num_dims(mesh)
  cache = get_cache(mesh)
  for d in 1:D-1
    for dface in 1:num_dfaces(mesh,d)
      if isactive(mesh,d,dface)
        if get_cell_dface(cache.cell,d,dface) != 0 # UNSET, is_cell_boundary_face
          for n in 1:d
            df_to_nf = get_faces(mesh,d,n)
            iface = UNSET
            min_distance = intersection_tolerance
            for lnface in 1:length(df_to_nf,dface)
              nface = df_to_nf[dface,lnface]
              dist = distance(mesh,n,nface,sm,sm_face) 
              if dist < min_distance
                min_distance = dist
                iface = nface
              end
            end
            if iface != UNSET
              point = closest_point(mesh,d,iface,sm,sm_face)
              return (d,iface,point)
            end
          end
        end
      end
    end
  end
  point = zero(get_vertex(mesh,1))
  (UNSET,UNSET,point)
end

function _is_nface_in_dface_boundary(mesh::CellMesh,d::Integer,dface::Integer,n::Integer,nface::Integer)
  @check d > n 
  df_to_nf = get_faces(mesh,d,n)
  for lnface in 1:length(df_to_nf,dface)
    if df_to_nf[dface,lnface] == nface
      return true
    end
  end
  false
end

function _get_opposite_face(mesh::CellMesh,d::Integer,dface::Integer,n::Integer,nface::Integer)
  opposite_face_dim = d-n-1
  df_to_of = get_faces(mesh,d,opposite_face_dim)
  of_to_v = get_dface_to_vertices(mesh,opposite_face_dim)
  nf_to_v = get_dface_to_vertices(mesh,n)

  for loface in 1:length(df_to_of,dface)
    opposite = true
    oface = df_to_of[dface,loface]
    for oface_lvertex in 1:length(of_to_v,oface), nface_lvertex in 1:length(nf_to_v,nface)
      if of_to_v[oface,oface_lvertex] == nf_to_v[nface,nface_lvertex]
        opposite = false
        break
      end
    end
    if opposite
      return oface
    end
  end
  @check false
  return UNSET
end

function _is_surface_mesh_face_on_vertex(mesh,vertex,sm,sm_face)
  cache = get_cache(mesh)
  sm_faces = cache.cell.vertex_to_surface_mesh_faces[vertex]
  sm_face ∈ sm_faces
end

function get_surface_mesh_faces_on_vertex(mesh::CellMesh,vertex::Integer)
  cache = get_cache(mesh)
  cache.cell.vertex_to_surface_mesh_faces[vertex]
end

function is_touching_surface_mesh_face(mesh::CellMesh,d::Integer,dface::Integer,sm_face::Integer)
  df_to_v = get_dface_to_vertices(mesh,d)
  for lvertex in 1:length(df_to_v,dface)
    vertex = df_to_v[dface,lvertex]
    if sm_face ∈ get_surface_mesh_faces_on_vertex(mesh,vertex)
      return true
    end
  end
  false
end

function find_next_point(mesh::CellMesh,dface::Integer,sm::SurfaceMesh,sm_face::Integer)
  D = num_dims(mesh)
  sm_d = dface_dimension(sm,sm_face)
  d = sm_d-1
  for n in d+1:D
    dim = n-d-1
    iface = UNSET
    min_distance = intersection_tolerance
    for nface in 1:num_dfaces(mesh,n)
      if isactive(mesh,n,nface)
        if _is_nface_in_dface_boundary(mesh,n,nface,d,dface)
          opposite_face = _get_opposite_face(mesh,n,nface,d,dface)
          if !is_touching_surface_mesh_face(mesh,dim,opposite_face,sm_face)
            dist = distance(mesh,dim,opposite_face,sm,sm_face)
            if dist < min_distance
              min_distance = dist
              iface = opposite_face
            end
          end
        end
      end
    end
    if iface != UNSET
      point = closest_point(mesh,dim,iface,sm,sm_face)
      return (dim,iface,point)
    end
  end
  point = zero(get_vertex(mesh,1))
  (UNSET,UNSET,point)
end

points = [ Point(-0.1,0.3), Point(0.2,0.5), Point(0.4,0.4), Point(1.,0.), Point(.8,-.2)  ]

facet_to_vertices = Table( [
  1 2 3 4 5;
  2 3 4 5 1] )

sm = SurfaceMesh(points,facet_to_vertices)

p0 = Point(0.0,0.0)
p1 = Point(1.0,1.0)
box = BoundingBox(p0,p1)

mesh = CellMesh(box)


for sm_face in 1:num_faces(sm)
  if dface_dimension(sm,sm_face) == 0
    
    sm_dface = local_dface(sm,sm_face,0)
    point = get_vertex(sm,sm_dface)

    _d, iface = find_closest_face(mesh,point)
    add_vertex!(mesh,_d,iface,point,sm_face)
  else

    d = dface_dimension(sm,sm_face)
    @check d == 1

    vertices = find_captured_surface_mesh_face_points!(mesh,sm,sm_face)

    if length(vertices) == 0
      _d, iface, point = find_intersection_point_on_boundary_face(mesh,sm,sm_face)
      _vertex = add_vertex!(mesh,_d,iface,point,sm_face)
      if _vertex != UNSET
        push!(vertices,_vertex)
      end
    end

    while length(vertices) > 0 
      vertex = pop!(vertices)
      _d, iface, point = find_next_point(mesh,vertex,sm,sm_face)
      _vertex = add_vertex!(mesh,_d,iface,point,sm_face)
      if _vertex != UNSET
        push!(vertices,_vertex)
      end
    end
  end
end

compact!(mesh)

compute_in_out!(mesh,sm)

writevtk(sm,"sm")
writevtk(mesh,"sub_mesh")

@test are_all_faces_defined(mesh)

end # module

