module CellMeshesTests

using STLCutter

using STLCutter: initialize!, num_dfaces, get_faces, get_cache, compact!, initialize_cache!, UNSET

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

p1 = Point(1,2,3)
p2 = Point(4,5,6)

b = BoundingBox(p1,p2)

m = CellMesh(b)


p1 = Point(1,2,4)
p2 = Point(5,3,5)

box = BoundingBox(p1,p2)

initialize!(m,box)

@test @allocated(initialize!(m,box)) == 0


p0 = Point(0.0,0.0)
p1 = Point(1.0,1.0)
box = BoundingBox(p0,p1)

mesh = CellMesh(box)

cache = get_cache(mesh)

initialize_cache!(cache,mesh)

stl_points = [ Point(0.3,0.3), Point(0.25,0.5), Point(0.5,0.5), Point(0.75,0.4), Point(0.0,0.5)  ]


for point in stl_points

  d,face = find_closest_face(mesh,point)

  if face != UNSET 
    add_vertex!(mesh,cache,d,face,point)
  end

end

compact!(mesh)

writevtk(mesh,"sub_mesh")



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

compact!(mesh)

writevtk(mesh,"sub_mesh3")



end # module

