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


# Some notes (probably better to add pictures)

# Insert a vertex 0-dim in D-cube

# Nothing to split in 0-dim

# Create a bisection partition of all 1:D objects



# Insert an edge 1-dim in D-cube

# Splits the D-1 face in the D-cube s.t. its closure contains intersects the edge

# Idem for the other endpoint of the edge

# Insert a vertex in the D-1-face 1 using previous algorithm (bisection of the D-1-face and all its faces)

# If the two faces do not share a D-2 face, create an extrusion mesh from face 1 partition to face 2 partition. We have the symbolic mesh.

# Next, generate the geometrical map such that 1) corners of the cube are mapped to the cube coordinates (in the physical space), 2) respect the vertices of face 1 using the algorithm for the D-1 vertex insertion, 3) idem for face 2.

# Else if they share a D-2 face, 

# do the same but gluing together the vertices in the shared edge (e.g., taking the average of the two coordinates)

# Define the other half of the cell as a degenerated cell (e.g., as the extrusion of the other D-2 face with the same direction, _opposite_) and the _opposite_ D-1 face


