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
