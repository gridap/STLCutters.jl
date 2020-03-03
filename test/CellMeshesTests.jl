module CellMeshesTests

using STLCutter

using STLCutter: initialize!, num_dfaces, get_faces, get_cache

import STLCutter: compact!

using Test

function STLCutter.compact!(m::CellMesh{D}) where D
  for d in 0:D, n in 0:d
    compact!( get_faces(mesh,d,n) )
  end
  #TODO: Update connectivities with new indexes
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

stl_points = [ Point(0.3,0.3), Point(0.25,0.5), Point(0.5,0.5), Point(0.75,0.4)  ]
point = stl_points[1]

const intersection_tolerance = 1e-10
for k in 1:4
  global stl_points, point
  point = stl_points[k]
  D = 2
  min_distance = intersection_tolerance 
  iface = 0
  for d in 0:D
    for i in 1:num_dfaces(mesh,d)
      if isactive(mesh,d,i)
        if distance(mesh,d,i,point) ≤ min_distance
          min_distance = distance(mesh,d,i,point)
          iface = i
        end
      end
    end
    if iface != 0
      add_vertex!(mesh,cache,d,iface,point)
      break
    end
  end
end

compact!(mesh)


writevtk(mesh,"sub_mesh")

end # module

