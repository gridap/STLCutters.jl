# I propose a two-fold algorithm:

# Part 1. Vertex insertion using a binary tree. As a result of this step we have
# a hybrid tree that involves a first level, the original mesh T, which can be
# represented as a structured mesh, an octree, etc. At each cell of this mesh,
# we create a binary tree (if the cell includes a vertex in its interior). In this
# step, we must create a vertex local indexing (at each cell in T) that will
# be suitable for nonconforming meshes, using the concept of slave and master
# vertices. As a result of this step, each cell in T that contains one vertex of
# the STL returns a new tree mesh TR. TR is simplexified in order to create S0.

# Part 2. For each cell s in SR (for each cell in T) we insert edges. There are
# many different ways to insert an edge. An edge is shared by two faces. One
# approach is to consider the plane that is the bisection of these two planes.
# Another approach (with a sligth reduction in the number of cells) is to take!
# one vertex of s and cut it by the plane that contains the edge and the vertex.
# There are multiple vertices, so we could define a way to get a _well-posed_
# one. On the other hand, the newly created ones can be freely moved if close
# to cell vertices and eliminate the cells freely. Edge insertion can be performed
# exactly as above but taking edges instead of planes and defined the bisection
# plane of an edge. After this process, we have a new mesh S1.

# Part 3. Cut s against the facets (or planes, since no edges) using the algorithm
# below, to SF.

# Some discussion:

# I note that in the previous
# discussions, we were considering a step for edges (via cut planes) and another
# step for faces. However, when the cut plane for the edge contains a face, we
# had issues that we tried to solve using hacks. This is not going to happen
# taking as the plane the bisection plane.

# Another option, radically different, is to insert edges using a plane
# determined by one of the STL facets that
# contain that edge. So we know that we are always in the _ill_ situation above.
# In fact, doing that, we reduce the number of cuts, since we reuse the edge
# insertion and one facet insertion and all the other edges on that face.

# But at the end of the day, doing this is just inserting the facets. Which is the
# problem with this? The plane can also be perturbed as above, since we are not
# really touching the STL but the cell / STL intersection. But we should modify
# the plane (locally) to make this change _exact_.

# The next part is IN/OUT.

# I have to remember
# how it is being finally done. In any case, we could consider two different
# strategies. One option is a ray-shooting technique (as in Fenics article).
# Another option is the following. Since you know which are the planes
# facets of the cells in SF are on, you can
# just check (using coordinates) whether the midpoint of the facet is on the cell
# facet (potentially after modification). If they are, you can mark IN/OUT. If none is on the real STL boundary,
# we cannot mark the vertices through this cell.

# If a vertex receives in and out from different cells, we should re-run Part 2
# with more precision. Only if we observe that this is really needed. Since we
# are not going to use flat cells for IN/OUT (unless absolutely necessary) it is
# probably not that problematic. In any case, this part has to be elaborated if
# needed.

# Track

# Pseudo-code to perform tetrahedral intersections
# The following code should be called for each cell obtained after
# simplexify the kd-tree cells that have intersection with planes
# We assume that we already have a set of planes that contain all the
# planes that cut that cell and possibly more
function intersection(t::Tetrahedron,Π::Vector{Plane},tol)
  π = first(Π)
  T = cut(t,π,tol)
  if length(Π) == 1
    R = T
  else
    R = SimplicialMesh()
    Πn = view(Π,2:end)
    for s ∈ get_cells(T)
      S = intersection(s,Πn,tol)
      append!(R,S)
    end
  end
  return R
end

# Plane tet intersection with tol
function intersection(π::Plane,t::Tetrahedron,tol)
  # TBI
  # return R
end

# Not intersection, better cut

# Tet plane intersection with tol using marching tetrahedra
function cut(t::Tetrahedron,π::Plane,tol)
  R = marching_tetrahedra(t,π,tol)
  merge_vertices!(R)
  return R
end

# In the previous method, we must mark facets on π
# R is a simplicial mesh, e.g., a Vector{Tetrahedron} + metadata
function marching_tetrahedra(t::Tetrahedron,π::Plane,tol)
  # TBI
  return R
end

# In this method, we take the vertices of the mesh R and merge the
# vertices within tol distance, e.g., creating a map that for a node
# i provides the node j with min index. Next, we replace i by j in all
# structs. Probably, we mark the tets in which this happen as _flat_.
# In any case, better not to eliminate these cells, even though they
# won't be used for integration purposes, they can be needed in the
# in / out algorithm. The in/out simply requires to mark vertices as
# in /out. I would only use flat cells for this if there is no other
# choice. It is hard to believe it will ever be needed. A node that
# only belongs to flat cells? I think it is provably it cannot happen.
# Since we are just cutting by faces (no vertex or edge insertion)
# I think we just need to compute vertex-to-vertex distance. We
# cannot create cells such that all vertices are far away but
# on the same plane.
function merge_vertices(R::Mesh)
  # TBI
  return R
end

# In-out for cells with non-empty Π is immediate from the face plane ownership
# We could create a try - catch + some validity test and increase tolerance
# Validity tests can be based on contradictory in-out

# We could create a cell type that has a cell + metadata (face owner)
struct CellFromIntersection{P} <: P where P <: Polytope
  p::P
  facet_owners::Vector{Int}
  is_flat::Bool
end

# As I said above, we want to keep track to the STL that owns the tet
# faces been created. As a result, when we intersect a plane with a tet
# there will be new facets that belong to the inserted face. On the other
# hand, if an existing cell is being cut, we must do the following. We have
# to check whether the new sub-cells are still on the STL face or are out of
# the geometry. That is obvious using the signed distance.


# the bang (!) comes from the fact that f can be modified within this subroutine
function marching_tetrahedra!(t:Tetrahedron,f::Face,tol)
  π = plane(f)
  f = intersect(π,f)     # trimmed face (overwrite the face)
  f, m = merge_vertices!(f,t)
  # merge ft vertices with t vertices
  # based on this aggregation, we can also return
  # whether this face is still cutting f or not
  # e.g., if m is -1 it still cuts, if m > 0, it returns the
  # face (vertex, ..., facet) containing it (integer computation)

  if m == -1 # not cutting
    πn = plane(f)
    T = cut(t,πn)
    for fT in get_faces(T)
      if belongs_to_face(fT,f)
        # mark ft as owned by f
      else if belongs_to_boundary(fT,t)
        # keep owner (if any)
      end
    end
    # new faces of T without t vertices are boundary facets (possibly only
    # part of them on the boundary) and we mark them
    # as boundary. For the facets on the boundary of t that were marked as
    # boundary, we keep the owner face in the STL
  else
    T = mesh(t)
    in_out_vertices!(T,t,πn)
    # if m is a facet of t and not marked yet, mark it as belonging to f
  end
  return T
end

# T0 is a global mesh, T1 is a local mesh
function in_out_sub_cells(T1::Grid,T0::Grid)
  is_def = false
  for s in T
    for f in get_facets(s)
      is_def = true
      mark_in_out!(s,f)
      # simply use signed distance of an s vertex not on f wrt s
    end
  end
  if is_def
    in_out_vertices(T1,T0) # obvious, all in -> in, all out -> out, in out -> boundary
    # only for the vertices of the root mesh T0, not the local sub-mesh T1
  end
  # by construction, all cells in the Grid are marked or
  # T is one tet that is not cut, in / out
end

# The T0 info can be later propagated on undef cells of T0