module CellMeshesTests

using STLCutter

using STLCutter: initialize!, num_dfaces, get_faces, get_cache, compact!, initialize_cache!, UNSET, dface_dimension, local_dface, global_dface, get_vertex,@check, get_face, get_dface, compute_in_out!, are_all_faces_defined, get_cell_dface

import STLCutter: compact!

using Test

const intersection_tolerance = 1e-10
function find_closest_face(mesh::CellMesh{D},point::Point{D}) where D
  min_distance = intersection_tolerance 
  iface = UNSET
  for d in 0:D
    for i in 1:num_dfaces(mesh,d)
      if isactive(mesh,d,i)
        if distance(mesh,d,i,point) ≤ min_distance
          min_distance = distance(mesh,d,i,point)
          iface = i
        end
      end
    end
    if iface != UNSET
      return (d,iface)
    end
  end
  (UNSET,UNSET)
end

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

function add_surface_mesh_face!(mesh,cache,sm,sm_face)
  d = dface_dimension(sm,sm_face)
  sm_dface = local_dface(sm,sm_face,d)

  boundary_vertices = Int[]
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
      add_vertex!(mesh,cache,dim,iface,point,sm_face)
    end

    
  end

  D = num_dims(mesh)
  queue = Int[]
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
      add_vertex!(mesh,cache,opposite_dim,iface,point,sm_face)
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

initialize!(m,box)

@test @allocated(initialize!(m,box)) == 0


## 2D Point Intersection

p0 = Point(0.0,0.0)
p1 = Point(1.0,1.0)
box = BoundingBox(p0,p1)

mesh = CellMesh(box)

cache = get_cache(mesh)

initialize_cache!(cache,mesh)

stl_points = [ Point(0.3,0.3), Point(0.25,0.5), Point(0.5,0.5), Point(0.75,0.4), Point(0.0,0.5)  ]


for point in stl_points[1:1]

  d,face = find_closest_face(mesh,point)

  if face != UNSET 
    add_vertex!(mesh,cache,d,face,point)
  end

end

compact!(mesh,cache)

writevtk(mesh,"sub_mesh")

## 3D Point Intersection

p0 = Point(0.0,0.0,0.0)
p1 = Point(1.0,1.0,1.0)
box = BoundingBox(p0,p1)

mesh = CellMesh(box)

cache = get_cache(mesh)

initialize_cache!(cache,mesh)

stl_points = [ 
  Point(0.3,0.3,0.0), 
  Point(0.25,0.5,0.1), 
  Point(0.4,0.6,0.0), 
  Point(0.75,0.4,0.2), 
  Point(0.0,0.5,0.0)  ]

for point in stl_points

  d,face = find_closest_face(mesh,point)
  
  if face != UNSET 
    add_vertex!(mesh,cache,d,face,point)
  end

end

compact!(mesh,cache)

writevtk(mesh,"sub_mesh3")

## 2D Surface Mesh

points = [ Point(-0.1,0.3), Point(0.2,0.5), Point(0.4,0.4), Point(1.,0.), Point(.8,-.2)  ]

facet_to_vertices = Table( [
  1 2 3 4 5;
  2 3 4 5 1] )

sm = SurfaceMesh(points,facet_to_vertices)

p0 = Point(0.0,0.0)
p1 = Point(1.0,1.0)
box = BoundingBox(p0,p1)

mesh = CellMesh(box)

cache = get_cache(mesh)

for sm_face in 1:num_faces(sm)
  if dface_dimension(sm,sm_face) == 0
    
    sm_dface = local_dface(sm,sm_face,0)
    sm_vertex = get_vertex(sm,sm_dface)

    d,face = find_closest_face(mesh,sm_vertex)

    if face != UNSET 
      add_vertex!(mesh,cache,d,face,sm_vertex,sm_face)
    end

  else
    add_surface_mesh_face!(mesh,cache,sm,sm_face)
  end
end


compact!(mesh,cache)

compute_in_out!(mesh,cache,sm)

writevtk(sm,"sm")
writevtk(mesh,"sub_mesh")

@test are_all_faces_defined(mesh)

end # module

