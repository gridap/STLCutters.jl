# Pseudo-code to perform tetrahedral intersections
# The following code should be called for each cell obtained after
# simplexify the kd-tree cells that have intersection with planes
# We assume that we already have a set of planes that contain all the
# planes that cut that cell and possibly more
function intersection(t::Tetrahedron,Π::Vector{Plane},tol)
  T = intersection(t,first(Π),tol)
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

# Tet plane intersection with tol using marching tetrahedra
function intersection(t::Tetrahedron,π::Plane,tol)
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
# in / out algorithm
function merge_vertices(R::Mesh)
  # TBI
  return R
end

# In-out for cells with non-empty Π is immediate from the face plane ownership
# We could create a try - catch + some validity test and increase tolerance
# Validity tests can be based on contradictory in-out
