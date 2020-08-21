# Vertex insertion meshes

# Octree mesh (1 case)
# Define here the resulting mesh

# Create a geometrical map from the reference partition
# to the physical cell such that:

# 1. The interior vertices is mapped to the inserted vertex
# 2. Cube vertices are mapped to the cell vertices
# 3. For the other vertices, they belong to a face f, and we extract the
# non-extruded coordinates from its anchor vertex and the extruded one from
# the interior vertex. E.g., if the vertex belongs to edge with ext = (001)
# in 3D and ach = 000, we take the coordinates x and y from the 000 vertex
# of the cell and the z component from the inserted vertex

# Line insertion meshes
# We have v1,v2 that intersect f1,f2 faces of the cell

# At least in 2D 3D works, I will think about 4D later, but I think easy generalisation
# Two cases:
# Case 1: The two faces share an edge of the cell
# Case 2: Otherwise, opposite faces

# Case 2:

# Insert the point v1 in f1 (previous algorithm)
# Extrude the resulting dim-1 mesh to create a D mesh
# Enforce that the vertices on the shared edge to be the same (merge them)
# * all this is symbolic not even introduced coordinates
# Enforce the coordinates of the vertices of the mesh in f1 as in the
# previous algorithm
# Enforce the interior vertex in f2 to be v2 as in the previous algorithm
# (when inserting v2 in f2)

# Case 1: Idem as case 2 but last step

# ...
# Enforce the interior vertex in f2 to be v2 and the still not defined
# nodes in this face as in the previous algorithm (when inserting v2 in f2)
# *but do not change the vertices in the shared edge that already have a
# coordinate

# Create a new cell as the original one but symbolically moving the shared
# edge to its opposite edge in f1, i.e., enforcing the two edges have the
# nodes of the non-shared edge. Keep the coordinates as in the original cell.
# * This way, we create a cell that covers the half of the cell not refined
# above

# Face insertion mesh:

# Create a mesh for 6 possible cases + rotations of the cube-plane intersection
# The idea is simple,
