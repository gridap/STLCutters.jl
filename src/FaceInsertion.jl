# VEFs distributions (floating point algorithm, not symbolic)

function vertices_distribution(Tk,Vk)
  # This part can be done easily with > < == (no operations, just check the box)
end

function edges_distribution(Tk,Ek)
  # Intersection of edges against Tk skeleton+boundary
  # Owner cells of the resulting points (previous algorithm)
end

function faces_distribution(Tk,Fk)
  # Intersection of faces against wire basket of Tk
  # Owner cells of the resulting points (previous algorithm)
end

function d_face_distribution(Tk,Fk)
  # Intersection of d-faces against the (D-d)-facets of Tk mesh
  # Owner cells of the resulting points (previous algorithm)
end

# Bullet-proof (?) intersection algorithms
# Some notes (probably better to add pictures)

# vertex-cell (1D,2D,3D)
# more abstract, 0-hyperplane (d=0) (D)-face intersections 

Lookup table mesh (only one case, the octree)

Aggregate nodes at <= TOL distance and move them to lowest dim face (e.g., never move the cell vertices) and eliminate zero-volume cells

# line-faces (2D, 3D)
# more abstract, 1-hyperplane (d=1) (D-1)-face intersections 

Compute line face intersections (if any).

If two different faces and two different points -> Well-posed intersection

If not, choose the two more distant points and their faces

If same face, the line belongs to the face, move the points to the face 

* We could perform line-line intersections on the face (D-1 space), but not sure it is needed, we could just say that the edge does not intersect the cell and belongs to the face, so there is nothing to do

Otherwise -> Well-posed intersection

Well posed intersection == exist f1 and f2 such that e intersect these two

Call the symbolic mesh refinement for 1-2

Aggregate nodes at <= TOL distance and move them to lowest dim face (e.g., never move the cell vertices) and eliminate zero-volume cells

# face-edges (3D)
# more abstract, 2-hyperplane (d=2) (D-2)-face intersections

Compute face edges intersections (if any)

Merge intersection vertices within <= TOL distance
* Since a plane can cut all edges in a cube, there is no ill-posed case here

Based on the number of edges being cut and the cut edges -> Lookup table mesh

Merge vertices of the mesh (always moving to vertices in lowest dim, e.g., never move cell vertices) and eliminate zero-volume cells

# 4D case
# volume-edges intersection
# 3-hyperplane (D-3)-face intersections
# it does not seem so difficult...



# Insert an edge 1-dim in D-cube

# Splits the D-1 face in the D-cube s.t. its closure contains intersects the edge

# Idem for the other endpoint of the edge

# Insert a vertex in the D-1-face 1 using previous algorithm (bisection of the D-1-face and all its faces)

# If the two faces do not share a D-2 face, create an extrusion mesh from face 1 partition to face 2 partition. We have the symbolic mesh.

# Next, generate the geometrical map such that 1) corners of the cube are mapped to the cube coordinates (in the physical space), 2) respect the vertices of face 1 using the algorithm for the D-1 vertex insertion, 3) idem for face 2.

# Else if they share a D-2 face, 

# do the same but gluing together the vertices in the shared edge (e.g., taking the coordinate of face 1)

# Define the other half of the cell as a degenerated cell (e.g., as the extrusion of the other D-2 face with the same direction, _opposite_) and the _opposite_ D-1 face


