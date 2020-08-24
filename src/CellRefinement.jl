
# Global face distribution algorithm

# For each _root_ cell K0 in the background mesh:

# F0 the STL faces with non-zero intersection with K0
# E0 the STL edges with non-zero intersection with K0
# V0 the STL edges that belong to K0
#* This is a global algorithm

T = {K0}    # background mesh
VK[K0] = V0 # cell vertices 
V = V0      # STL vertices to be inserted

# Local refinement algorithm

T = [K]
Tnew = []

insert_vertices!(T,V,Tnew)

# Insert vertices
V = [[V0]]
# After calling this function, we have new mesh Tnew in which there are no STL vertices
function insert_vertices!(T,V,Tnew)
  for k,Vk in (T,V)
    if length(Vk) != 0
      v = Vk[1]
      Tk = vertex_refinement(k,v)
      VTk = distribute_vertices(Tk,Vk[2:end])
      insert_vertices(Tk,VTk,Tnew)
    else 
      push(Tnew,k) # a cell is defined by its parent, the Morton index, and v
  end
end
T = Tnew
Tnew = []

# It is the same as the one for vertices. We can just parameterise wrt d !
# E1[k] the STL edges in E0 with non-zero intersection with k in T
E1 = edges_distribution(T,E0)
function insert_edges!(T,E,Tnew)
  for k,Ek in (T,E)
    if length(Ek) != 0
      e = Ek[1]
      Tk = edge_refinement(k,e)
      ETk = distribute_edges(Tk,Ek[2:end])
      insert_edges!(Tk,ETk,Tnew)
    else 
      push(Tnew,k) # a cell is defined by its parent, the Morton index, and v
  end
end
T = Tnew
Tnew = []
