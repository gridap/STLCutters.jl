module VolumeMeshesTests

using STLCutter
using STLCutter: get_face_coordinates, get_vertex_coordinates, get_faces_to_vertices


box = BoundingBox(
  Point( 0.2, 0.2, 0.2 ),
  Point( 0.7, 0.7, 0.7 ) )

bg_mesh = CartesianMesh( box, 2 )

vm = VolumeMesh(bg_mesh)

@show get_face_coordinates(vm,Val(2),1)

i = 1
v = get_vertex_coordinates(vm) 
f_to_v = get_faces_to_vertices(vm,2) 
@show ( v[f_to_v[i,1]], v[f_to_v[i,2]], v[f_to_v[i,3]], v[f_to_v[i,4]],  )
# TODO: Complete test



end # module
