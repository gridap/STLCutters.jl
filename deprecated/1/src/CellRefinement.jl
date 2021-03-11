
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
    
# @santiagobadia: I don't think you need a Morton index (tree struct) for IN-OUT. 
# I would just keep the cell-wise vertices IDs consistent among cells 
# (i.e., two cells that share a vertex have the same background cell-wise Id). With this, 
# keeping the vef to STL ownership (i.e., when inserting a vertex-edge, explictly keep the 
# info that this vertex-edge is vertex-edge X of the STL), and doing the face-cell intersection
# with care (i.e., the case the face is almost aligned with a face, which could be solved using 
# standard level set and probably extending the cell edges with some epsilon, in order to handle 
# these cases) will be enough. For the moment, I would not care about Morton indices and trees, 
# if I am wrong, we could consider this in the future. I think I can prove it works.
